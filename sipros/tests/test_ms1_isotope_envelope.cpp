#include "PSMfeatureExtractor.h"
#include "SiprosSearchRunner.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

void check(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

long long processId()
{
#ifdef _WIN32
    return static_cast<long long>(_getpid());
#else
    return static_cast<long long>(getpid());
#endif
}

struct TargetCase
{
    const char *label;
    int atomIndex;
    int isotopeIndex;
    int nominalShift;
    double expectedMassDelta;
    double naturalPct;
    std::string peptide;
};

const std::vector<TargetCase> &targetCases()
{
    static const std::vector<TargetCase> cases{
        {"C13", 0, 1, 1, 1.003355, 1.07, "CCCCCCK"},
        {"H2", 1, 1, 1, 1.006277, 0.0115, "CCCCCCK"},
        {"N15", 3, 1, 1, 0.997035, 0.368, "CCCCCCK"},
        {"O18", 2, 2, 2, 2.004245, 0.205, "CCCCCCK"},
        {"S34", 5, 2, 2, 1.995796, 4.29, "CCCCCCK"},
    };
    return cases;
}

struct NominalDistribution
{
    std::vector<double> probability{1.0};
    std::vector<double> massMoment{0.0};
};

void multiplyByAtom(NominalDistribution &distribution,
                    const std::vector<double> &probabilities,
                    const std::vector<double> &masses)
{
    check(probabilities.size() == masses.size() && !masses.empty(),
          "invalid atomic isotope distribution");
    int maximumShift = 0;
    std::vector<int> nominalShifts(masses.size(), 0);
    for (size_t isotope = 0; isotope < masses.size(); ++isotope)
    {
        nominalShifts[isotope] = static_cast<int>(
            std::lround(masses[isotope] - masses[0]));
        maximumShift =
            std::max(maximumShift, nominalShifts[isotope]);
    }

    std::vector<double> nextProbability(
        distribution.probability.size() +
            static_cast<size_t>(maximumShift),
        0.0);
    std::vector<double> nextMassMoment(
        nextProbability.size(), 0.0);
    for (size_t index = 0;
         index < distribution.probability.size();
         ++index)
    {
        for (size_t isotope = 0;
             isotope < probabilities.size();
             ++isotope)
        {
            if (!(probabilities[isotope] > 0.0))
                continue;
            const size_t nextIndex =
                index +
                static_cast<size_t>(nominalShifts[isotope]);
            const double massDelta =
                masses[isotope] - masses[0];
            nextProbability[nextIndex] +=
                distribution.probability[index] *
                probabilities[isotope];
            nextMassMoment[nextIndex] +=
                distribution.massMoment[index] *
                    probabilities[isotope] +
                distribution.probability[index] *
                    probabilities[isotope] * massDelta;
        }
    }
    distribution.probability = std::move(nextProbability);
    distribution.massMoment = std::move(nextMassMoment);
}

NominalDistribution peptideDistribution(
    const sipros::SourcedComposition &composition,
    const TargetCase &target,
    double targetFraction)
{
    NominalDistribution distribution;
    const auto &configured = ProNovoConfig::configIsotopologue
                                 .vNaturalAtomIsotopicDistribution;
    for (size_t source = 0;
         source < sipros::IsotopeSourceCount;
         ++source)
    {
        const bool biosynthetic =
            source == static_cast<size_t>(
                          sipros::IsotopeSource::Biosynthetic);
        for (size_t atom = 0; atom < sipros::ElementCount; ++atom)
        {
            std::vector<double> probabilities =
                ProNovoConfig::getNaturalAtomIsotopeProbabilities(atom);
            if (biosynthetic &&
                static_cast<int>(atom) == target.atomIndex)
            {
                check(PeptideIsotopeCalculator::changeAtomProbability(
                          probabilities,
                          target.label[0],
                          targetFraction),
                      "failed to set target isotope abundance");
            }
            for (int count = 0;
                 count < composition.atoms[source][atom];
                 ++count)
            {
                multiplyByAtom(
                    distribution,
                    probabilities,
                    configured[atom].vMass);
            }
        }
    }

    const double sum = std::accumulate(
        distribution.probability.begin(),
        distribution.probability.end(),
        0.0);
    check(sum > 0.0, "peptide distribution is empty");
    for (size_t index = 0;
         index < distribution.probability.size();
         ++index)
    {
        distribution.probability[index] /= sum;
        distribution.massMoment[index] /= sum;
    }
    return distribution;
}

double centroidMassDelta(const NominalDistribution &distribution,
                         int nominalIndex)
{
    if (nominalIndex < 0 ||
        nominalIndex >=
            static_cast<int>(distribution.probability.size()) ||
        !(distribution.probability[
              static_cast<size_t>(nominalIndex)] > 0.0))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return distribution.massMoment[
               static_cast<size_t>(nominalIndex)] /
           distribution.probability[
               static_cast<size_t>(nominalIndex)];
}

struct Fixture
{
    const TargetCase *target = nullptr;
    std::string peptide;
    int charge = 0;
    int targetAtomCount = 0;
    double expectedPct = 0.0;
    double baseMass = 0.0;
    double baseMz = 0.0;
    int modeIndex = 0;
    sipros::SourcedComposition composition;
    NominalDistribution distribution;
    sipros::RaxportMs1Data ms1;
};

Fixture makeFixture(const TargetCase &target,
                    double targetFraction,
                    int charge)
{
    Fixture fixture;
    fixture.target = &target;
    fixture.peptide = target.peptide;
    fixture.charge = charge;
    fixture.expectedPct = targetFraction * 100.0;

    PeptideIsotopeCalculator calculator;
    fixture.baseMass =
        calculator.calPrecursorBaseMass(fixture.peptide);
    fixture.composition = calculator.pepComposition;
    fixture.targetAtomCount =
        fixture.composition[sipros::IsotopeSource::Biosynthetic]
                           [static_cast<size_t>(target.atomIndex)];
    check(fixture.targetAtomCount > 0,
          std::string(target.label) +
              " test peptide has no target atoms");
    fixture.baseMz =
        fixture.baseMass / charge +
        ProNovoConfig::getProtonMass();
    fixture.distribution = peptideDistribution(
        fixture.composition, target, targetFraction);
    fixture.modeIndex = static_cast<int>(
        std::max_element(
            fixture.distribution.probability.begin(),
            fixture.distribution.probability.end()) -
        fixture.distribution.probability.begin());
    const double maximumProbability =
        *std::max_element(
            fixture.distribution.probability.begin(),
            fixture.distribution.probability.end());

    sipros::RaxportMs1Scan scan;
    scan.scanNumber = 100;
    for (size_t nominalIndex = 0;
         nominalIndex <
             fixture.distribution.probability.size();
         ++nominalIndex)
    {
        const double probability =
            fixture.distribution.probability[nominalIndex];
        if (probability < maximumProbability * 1e-10)
            continue;
        const double massDelta = centroidMassDelta(
            fixture.distribution,
            static_cast<int>(nominalIndex));
        scan.mz.push_back(
            fixture.baseMz + massDelta / charge);
        scan.intensity.push_back(probability * 1e8);
        scan.charge.push_back(charge);
    }
    fixture.ms1.scanNumberToIndex[100] = 0;
    fixture.ms1.scans.push_back(std::move(scan));
    return fixture;
}

void sortScan(sipros::RaxportMs1Scan &scan)
{
    std::vector<size_t> order(scan.mz.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(
        order.begin(), order.end(),
        [&](size_t left, size_t right)
        { return scan.mz[left] < scan.mz[right]; });

    std::vector<double> mz;
    std::vector<double> intensity;
    std::vector<int> charge;
    for (size_t index : order)
    {
        mz.push_back(scan.mz[index]);
        intensity.push_back(scan.intensity[index]);
        charge.push_back(scan.charge[index]);
    }
    scan.mz = std::move(mz);
    scan.intensity = std::move(intensity);
    scan.charge = std::move(charge);
}

double fixtureMzAt(const Fixture &fixture,
                   int nominalIndex)
{
    return fixture.baseMz +
           centroidMassDelta(
               fixture.distribution, nominalIndex) /
               fixture.charge;
}

std::pair<std::vector<double>, std::vector<double>>
fixturePrecursorEnvelope(const Fixture &fixture)
{
    std::vector<double> masses;
    std::vector<double> probabilities;
    for (size_t nominalIndex = 0;
         nominalIndex < fixture.distribution.probability.size();
         ++nominalIndex)
    {
        const double probability =
            fixture.distribution.probability[nominalIndex];
        if (!(probability > 0.0))
            continue;
        masses.push_back(
            fixture.baseMass + centroidMassDelta(
                fixture.distribution,
                static_cast<int>(nominalIndex)));
        probabilities.push_back(probability);
    }
    return {std::move(masses), std::move(probabilities)};
}

PSMfeatureExtractor::Ms1AbundanceResult fitFixture(
    Fixture &fixture,
    std::vector<isotopicPeak> *found = nullptr,
    double initializerPct = -1.0)
{
    if (initializerPct < 0.0)
        initializerPct = fixture.expectedPct;
    int scanNumber = 100;
    const double matchedMz =
        fixtureMzAt(fixture, fixture.modeIndex);
    const auto tolerance = [&](double)
    { return 0.01 / fixture.charge; };
    std::vector<isotopicPeak> peaks =
        PSMfeatureExtractor::findMs1IsotopicPeaks(
            &fixture.ms1,
            scanNumber,
            fixture.charge,
            fixture.baseMz,
            matchedMz,
            fixture.composition,
            fixture.target->label,
            initializerPct,
            tolerance);
    if (found)
        *found = peaks;
    return PSMfeatureExtractor::
        getSIPelementAbundanceFromMS1Peaks(
            peaks,
            fixture.baseMass,
            fixture.peptide,
            fixture.charge,
            fixture.target->label,
            initializerPct);
}

sipPSM extractFixtureFeatures(
    Fixture &fixture,
    double precursorRtDiffSeconds = 0.0)
{
    const double matchedMz =
        fixtureMzAt(fixture, fixture.modeIndex);

    sipPSM psm;
    psm.scanNumbers.push_back(200);
    psm.precursorScanNumbers.push_back(100);
    psm.parentCharges.push_back(fixture.charge);
    psm.measuredParentMasses.push_back(
        (matchedMz - ProNovoConfig::getProtonMass()) *
        fixture.charge);
    psm.calculatedParentMasses.push_back(
        fixture.baseMass);
    psm.precursorNeutronMasses.push_back(
        PeptideIsotopeCalculator().calPrecursorEstimate(fixture.peptide).neutronMass);
    psm.identifiedPeptides.push_back(
        "[" + fixture.peptide + "]");
    psm.originalPeptides.push_back(
        "[" + fixture.peptide + "]");
    psm.MS2IsotopicAbundances.push_back(
        fixture.expectedPct);
    psm.ranks.push_back(1);
    psm.WDPscores.push_back(1.0f);
    psm.isolationWindowCenterMZs.push_back(matchedMz);
    psm.precursorRtDiffSeconds.push_back(precursorRtDiffSeconds);

    PSMfeatureExtractor extractor;
    extractor.ms1Data = fixture.ms1;
    extractor.mSipPSM = &psm;
    extractor.initializeFeatureVectors(psm);
    extractor.extractFeaturesOfEachPSM();
    return psm;
}

void checkMissingPrecursorUsesIsolationCenter()
{
    Fixture fixture = makeFixture(targetCases()[0], 0.5, 2);
    const sipPSM psm = extractFixtureFeatures(fixture, -1.0);
    check(psm.isotopicPeakNumbers[0] == 0 &&
              psm.MS1IsotopeFitScores[0] == 0.0,
          "isolation-center fallback reported peptide-specific MS1 evidence");
    check(psm.MS1IsotopicAbundances[0] > 0.0,
          "isolation-center fallback did not retain the imputed MS1 abundance");
    check(std::fabs(psm.isotopicAbundanceDiffs[0] -
                    (psm.MS1IsotopicAbundances[0] -
                     psm.MS2IsotopicAbundances[0])) < 1e-12,
          "isolation-center fallback produced an inconsistent abundance difference");

    fixture.ms1.scans[0].mz.clear();
    fixture.ms1.scans[0].intensity.clear();
    fixture.ms1.scans[0].charge.clear();
    const sipPSM emptyPsm = extractFixtureFeatures(fixture, -1.0);
    const double expectedCenterEstimate =
        PSMfeatureExtractor::estimateSIPelementAbundanceFromIsolationCenter(
            fixture.peptide, fixture.charge,
            fixtureMzAt(fixture, fixture.modeIndex),
            fixture.target->label, fixture.expectedPct);
    check(emptyPsm.isotopicPeakNumbers[0] == 0 &&
              emptyPsm.MS1IsotopeFitScores[0] == 0.0,
          "theoretical isolation-center fallback reported MS1 evidence");
    check(expectedCenterEstimate > 0.0 &&
              std::fabs(emptyPsm.MS1IsotopicAbundances[0] -
                        expectedCenterEstimate) < 1e-12,
          "missing center anchor did not use the closest theoretical isotope");
    check(std::fabs(emptyPsm.isotopicAbundanceDiffs[0] -
                    (expectedCenterEstimate -
                     emptyPsm.MS2IsotopicAbundances[0])) < 1e-12,
          "theoretical isolation-center fallback did not report MS1 minus MS2");
}

void checkWhitelistAndConfiguredMasses()
{
    for (const TargetCase &target : targetCases())
    {
        check(
            PSMfeatureExtractor::canonicalSipIsotope(
                target.label) == target.label,
            std::string("canonical mapping failed for ") +
                target.label);
        check(
            PSMfeatureExtractor::sipAtomIndex(
                target.label) == target.atomIndex,
            std::string("atom mapping failed for ") +
                target.label);
        check(
            PSMfeatureExtractor::sipIsotopeIndex(
                target.label) == target.isotopeIndex,
            std::string("isotope mapping failed for ") +
                target.label);
        check(
            PSMfeatureExtractor::sipNominalShiftPerAtom(
                target.label) == target.nominalShift,
            std::string("nominal shift failed for ") +
                target.label);

        const auto &distribution =
            ProNovoConfig::configIsotopologue
                .vAtomIsotopicDistribution[
                    static_cast<size_t>(target.atomIndex)];
        const double massDelta =
            distribution.vMass[
                static_cast<size_t>(target.isotopeIndex)] -
            distribution.vMass[0];
        check(std::fabs(
                  massDelta - target.expectedMassDelta) <
                  1e-6,
              std::string("configured mass delta failed for ") +
                  target.label);

        char atom = '\0';
        int massNumber = 0;
        check(sipros::TextUtils::parseSipAtomSpec(
                  target.label, atom, massNumber),
              std::string("CLI parser rejected ") +
                  target.label);
        check(
            PSMfeatureExtractor::canonicalSipIsotope(
                std::string(1, atom), massNumber) ==
                target.label,
            std::string("metadata mapping failed for ") +
                target.label);
    }

    const std::vector<std::string> unsupported{
        "", "C", "H", "N", "O", "S", "P", "P31",
        "C12", "C14", "H1", "H3", "N14", "O17",
        "S33", "S36", "X99"};
    for (const std::string &label : unsupported)
    {
        char atom = '\0';
        int massNumber = 0;
        check(!sipros::TextUtils::parseSipAtomSpec(
                  label, atom, massNumber),
              "CLI parser accepted unsupported " + label);
    }
    check(PSMfeatureExtractor::sipAtomIndex("P") == -1,
          "phosphorus was accepted as a SIP target");
    check(
        PSMfeatureExtractor::canonicalSipIsotope(
            "O", 17).empty(),
        "mismatched O17 metadata mapped to O18");
    check(!ProNovoConfig::applySipAbundance('P', 0.5),
          "phosphorus abundance mutation was accepted");

    bool rejectedS33 = false;
    try
    {
        (void)ProNovoConfig::resolveSipIsotopeIndex(
            ProNovoConfig::configIsotopologue, 'S', 33);
    }
    catch (const std::exception &)
    {
        rejectedS33 = true;
    }
    check(rejectedS33,
          "configuration resolver accepted S33");

    bool rejectedPhosphorus = false;
    try
    {
        (void)ProNovoConfig::resolveSipIsotopeIndex(
            ProNovoConfig::configIsotopologue, 'P', 31);
    }
    catch (const std::exception &)
    {
        rejectedPhosphorus = true;
    }
    check(rejectedPhosphorus,
          "configuration resolver accepted phosphorus");
}

void checkModalPrecursorAndWindows()
{
    PeptideIsotopeCalculator calculator;
    const std::string peptide = "KVDVTGTSK";

    const auto checkModalShift = [&calculator, &peptide](
        const PeptideIsotopeCalculator::PrecursorEstimate &estimate,
        const std::string &label)
    {
        IsotopeDistribution envelope;
        check(ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
                  calculator.pepComposition, envelope),
              label + " precursor convolution failed");
        const auto apex = std::max_element(
            envelope.vProb.begin(), envelope.vProb.end());
        check(apex != envelope.vProb.end(),
              label + " precursor envelope is empty");
        const size_t apexIndex = static_cast<size_t>(
            std::distance(envelope.vProb.begin(), apex));
        const double baseMass = calculator.calPrecursorBaseMass(peptide);
        const int modalShift = static_cast<int>(std::lround(
            envelope.vMass[apexIndex] - baseMass));
        check(estimate.nominalShift == modalShift,
              label + " precursor estimate did not use the modal shift");
        check(std::fabs(
                  estimate.mass -
                  (baseMass + modalShift * estimate.neutronMass)) < 1e-9,
              label + " modal precursor identity failed");
    };

    check(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
          "failed to load Regular profile for precursor estimate test");
    const auto regular = calculator.calPrecursorEstimate(peptide);
    check(std::fabs(regular.neutronMass - 1.002786935) < 1e-6,
          "Regular composition-weighted neutron mass is wrong");
    checkModalShift(regular, "Regular");

    std::vector<std::pair<double, double>> windows;
    check(ProNovoConfig::getPeptideMassWindows(
              regular.mass, regular.neutronMass, windows),
          "failed to build peptide-specific precursor windows");
    check(windows.size() == 5,
          "peptide-specific precursor windows have the wrong count");
    for (int shift = -2; shift <= 2; ++shift)
    {
        const auto &window = windows[static_cast<size_t>(shift + 2)];
        const double center = (window.first + window.second) / 2.0;
        check(std::fabs(
                  center -
                  (regular.mass + shift * regular.neutronMass)) < 1e-9,
              "precursor window did not use peptide-specific spacing");
    }
    std::vector<std::pair<double, double>> invalidWindows;
    check(!ProNovoConfig::getPeptideMassWindows(
               regular.mass, 0.0, invalidWindows),
          "invalid peptide spacing fell back to the global neutron mass");

    check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
          "failed to restore SIP profile for precursor estimate test");
    check(ProNovoConfig::applySipAbundance('C', 0.5),
          "failed to set 50% C13 for precursor estimate test");
    const auto sip = calculator.calPrecursorEstimate(peptide);
    check(std::fabs(sip.neutronMass - 1.003339560) < 1e-6,
          "SIP composition-weighted neutron mass is wrong");
    checkModalShift(sip, "SIP");
}

void checkMassErrorUsesPeptideSpecificSpacing()
{
    PSMfeatureExtractor extractor;
    constexpr double calculatedMass = 1000.0;
    for (double peptideSpacing : {1.002786935, 1.003339560})
    {
        const auto result = extractor.getMassWindowShiftAndError(
            calculatedMass + peptideSpacing,
            calculatedMass,
            peptideSpacing);
        check(result.first == 1,
              "PIN isotope shift did not use peptide-specific spacing");
        check(std::fabs(result.second) < 1e-9,
              "PIN mass error did not use peptide-specific spacing");
    }
}

void checkAllTargetsProduceModalEstimates()
{
    PeptideIsotopeCalculator calculator;
    for (const TargetCase &target : targetCases())
    {
        check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
              "failed to reset SIP profile for modal target test");
        check(ProNovoConfig::applySipAbundance(target.label[0], 0.55),
              "failed to set modal target abundance");
        const auto estimate =
            calculator.calPrecursorEstimate(target.peptide);
        check(estimate.mass > 0.0 && estimate.neutronMass > 0.0 &&
                  std::isfinite(estimate.mass) &&
                  std::isfinite(estimate.neutronMass),
              std::string("invalid modal estimate for ") +
                  target.label);
        const double baseMass =
            calculator.calPrecursorBaseMass(target.peptide);
        IsotopeDistribution envelope;
        check(ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
                  calculator.pepComposition, envelope),
              std::string("precursor convolution failed for ") +
                  target.label);
        const auto apex = std::max_element(
            envelope.vProb.begin(), envelope.vProb.end());
        const size_t apexIndex = static_cast<size_t>(
            std::distance(envelope.vProb.begin(), apex));
        check(estimate.nominalShift == static_cast<int>(std::lround(
                  envelope.vMass[apexIndex] - baseMass)),
              std::string("modal shift failed for ") + target.label);
        check(std::fabs(
                  estimate.mass -
                  (baseMass + estimate.nominalShift *
                      estimate.neutronMass)) < 1e-9,
              std::string("modal identity failed for ") +
                  target.label);
    }
}

void checkParallelModalEstimates()
{
    check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
          "failed to load SIP profile for parallel modal test");
    check(ProNovoConfig::applySipAbundance('C', 0.25),
          "failed to set C13 abundance for parallel modal test");

    std::vector<std::string> peptides;
    peptides.reserve(64 * targetCases().size());
    for (int repeat = 0; repeat < 64; ++repeat)
        for (const TargetCase &target : targetCases())
            peptides.push_back(target.peptide);

    std::vector<PeptideIsotopeCalculator::PrecursorEstimate> serial(
        peptides.size());
    std::vector<PeptideIsotopeCalculator::PrecursorEstimate> parallel(
        peptides.size());
    PeptideIsotopeCalculator serialCalculator;
    for (size_t i = 0; i < peptides.size(); ++i)
        serial[i] = serialCalculator.calPrecursorEstimate(peptides[i]);

#pragma omp parallel
    {
        PeptideIsotopeCalculator calculator;
#pragma omp for schedule(guided)
        for (int i = 0; i < static_cast<int>(peptides.size()); ++i)
            parallel[static_cast<size_t>(i)] =
                calculator.calPrecursorEstimate(
                    peptides[static_cast<size_t>(i)]);
    }

    for (size_t i = 0; i < peptides.size(); ++i)
    {
        check(parallel[i].mass == serial[i].mass &&
                  parallel[i].neutronMass == serial[i].neutronMass &&
                  parallel[i].nominalShift == serial[i].nominalShift,
              "parallel modal estimate differs from serial result");
    }
}

void checkDirectSearchUsesFullPtmComposition()
{
    std::string error;
    check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
          "failed to load SIP profile for direct-search PTM mass test");
    check(ProNovoConfig::configureVariablePtms(
              {"deamidation"}, 1, error),
          error);
    check(ProNovoConfig::applySipAbundance('N', 0.55),
          "failed to set N15 abundance for direct-search PTM mass test");

    PeptideIsotopeCalculator calculator;
    const std::string plain = "NAAAAHK";
    const std::string modified = "N!AAAAHK";
    const auto plainEstimate =
        calculator.calPrecursorEstimate(plain);
    const auto modifiedEstimate =
        calculator.calPrecursorEstimate(modified);

    struct TemporaryFasta
    {
        std::string path;
        ~TemporaryFasta()
        {
            std::remove(path.c_str());
        }
    } temporary{
        (std::filesystem::temp_directory_path() /
         ("sipros5_direct_ptm_precursor_" +
          std::to_string(processId()) + ".fasta"))
            .string()};
    {
        std::ofstream output(temporary.path);
        check(output.good(), "failed to create direct-search PTM fixture");
        output << ">ptm_mass_fixture\n" << plain << "\n";
    }

    ProNovoConfig::setFASTAfilename(temporary.path);
    ProteinDatabase database;
    database.loadDatabase();
    check(database.getFirstProtein(),
          "failed to read direct-search PTM fixture");

    bool foundPlain = false;
    bool foundModified = false;
    Peptide peptide;
    for (int candidate = 0;
         candidate < 8 && database.getNextPeptide(&peptide);
         ++candidate)
    {
        if (peptide.getPeptideSeq() == "[NAAAAHK]")
        {
            foundPlain = true;
            check(peptide.getPeptideMass() == 0.0 &&
                      peptide.getPrecursorNeutronMass() == 0.0,
                  "database generation did not defer precursor estimation");
        }
        else if (peptide.getPeptideSeq() == "[N!AAAAHK]")
        {
            foundModified = true;
            check(peptide.getPeptideMass() == 0.0 &&
                      peptide.getPrecursorNeutronMass() == 0.0,
                  "PTM generation did not defer precursor estimation");
        }
    }
    check(foundPlain && foundModified,
          "direct search did not generate both PTM mass fixtures");
    check(plainEstimate.mass != modifiedEstimate.mass,
          "full modal estimator ignored the PTM composition");
    check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
          "failed to restore SIP profile after direct-search PTM mass test");
}

void checkTargetElementCategories()
{
    const auto &naturalOxygen =
        ProNovoConfig::getNaturalAtomIsotopeProbabilities(2);
    std::vector<double> oxygen = naturalOxygen;
    check(PeptideIsotopeCalculator::changeAtomProbability(
              oxygen, 'O', 0.5),
          "failed to prepare O18 categories");
    const double oxygenScale =
        0.5 / (1.0 - naturalOxygen[2]);
    check(std::fabs(
              oxygen[0] -
              naturalOxygen[0] * oxygenScale) < 1e-12 &&
              std::fabs(
                  oxygen[1] -
                  naturalOxygen[1] * oxygenScale) < 1e-12 &&
              std::fabs(oxygen[2] - 0.5) < 1e-12,
          "O16/O17 natural ratio was not retained for O18");

    const auto &naturalSulfur =
        ProNovoConfig::getNaturalAtomIsotopeProbabilities(5);
    std::vector<double> sulfur = naturalSulfur;
    check(PeptideIsotopeCalculator::changeAtomProbability(
              sulfur, 'S', 0.5),
          "failed to prepare S34 categories");
    const double sulfurScale =
        0.5 / (1.0 - naturalSulfur[2]);
    check(std::fabs(
              sulfur[0] -
              naturalSulfur[0] * sulfurScale) < 1e-12 &&
              std::fabs(
                  sulfur[1] -
                  naturalSulfur[1] * sulfurScale) < 1e-12 &&
              std::fabs(sulfur[2] - 0.5) < 1e-12 &&
              std::fabs(
                  sulfur[3] -
                  naturalSulfur[3] * sulfurScale) < 1e-12 &&
              std::fabs(
                  sulfur[4] -
                  naturalSulfur[4] * sulfurScale) < 1e-12,
          "S32/S33/S36 natural ratios were not retained for S34");

    check(PeptideIsotopeCalculator::changeAtomProbability(
              oxygen, 'O', 1.0),
          "failed to prepare O18 endpoint");
    check(oxygen[0] == 0.0 && oxygen[1] == 0.0 &&
              oxygen[2] == 1.0,
          "O18 endpoint did not suppress other isotopes");
    check(PeptideIsotopeCalculator::changeAtomProbability(
              sulfur, 'S', 1.0),
          "failed to prepare S34 endpoint");
    check(sulfur[0] == 0.0 && sulfur[1] == 0.0 &&
              sulfur[2] == 1.0 &&
              sulfur[3] == 0.0 && sulfur[4] == 0.0,
          "S34 endpoint did not suppress other isotopes");
}

void checkCleanEnrichments()
{
    for (const TargetCase &target : targetCases())
    {
        std::vector<double> percentages{
            0.0, target.naturalPct, 10.0,
            50.0, 90.0};
        if (std::string(target.label) == "O18")
            percentages.push_back(99.99);
        if (std::string(target.label) == "S34")
            percentages.push_back(99.5);
        percentages.push_back(100.0);
        double previous = -1.0;
        for (double percentage : percentages)
        {
            Fixture fixture = makeFixture(
                target, percentage / 100.0, 2);
            const auto result = fitFixture(fixture);
            check(result.valid,
                  std::string(target.label) +
                      " clean fit was invalid at " +
                      std::to_string(percentage));
            check(result.isotopicPeakCount > 0,
                  std::string(target.label) +
                      " clean fit returned no peaks");
            check(std::isfinite(result.abundancePct) &&
                      result.abundancePct >= 0.0 &&
                      result.abundancePct <= 100.0,
                  std::string(target.label) +
                      " clean fit escaped [0,100]");
            check(std::fabs(
                      result.abundancePct - percentage) <=
                      0.35,
                  std::string(target.label) +
                      " clean fit missed " +
                      std::to_string(percentage) +
                      "%: " +
                      std::to_string(result.abundancePct));
            check(result.abundancePct + 1e-9 >= previous,
                  std::string(target.label) +
                      " estimates were not monotonic");
            previous = result.abundancePct;
        }
    }
}

void checkBiasedInitializersAndCharges()
{
    for (const TargetCase &target : targetCases())
    {
        double reference = -1.0;
        std::vector<int> referenceIndices;
        for (int charge = 1; charge <= 4; ++charge)
        {
            Fixture fixture =
                makeFixture(target, 0.5, charge);
            std::vector<isotopicPeak> peaks;
            const auto result =
                fitFixture(fixture, &peaks, 45.0);
            std::vector<int> indices;
            for (const isotopicPeak &peak : peaks)
                indices.push_back(peak.isotopeIndex);

            check(result.valid,
                  std::string(target.label) +
                      " biased-initializer fit invalid");
            check(std::fabs(
                      result.abundancePct - 50.0) <= 1.0,
                  std::string(target.label) +
                      " fit echoed 45% initializer: " +
                      std::to_string(result.abundancePct));
            if (charge == 1)
            {
                reference = result.abundancePct;
                referenceIndices = indices;
            }
            else
            {
                check(indices == referenceIndices,
                      std::string(target.label) +
                          " charge changed isotope indices");
                check(std::fabs(
                          result.abundancePct - reference) <=
                          0.05,
                      std::string(target.label) +
                          " charge changed abundance");
            }
        }
    }
}

void eraseNominalPeak(Fixture &fixture, int nominalIndex)
{
    auto &scan = fixture.ms1.scans.front();
    const double targetMz =
        fixtureMzAt(fixture, nominalIndex);
    for (size_t index = scan.mz.size(); index-- > 0;)
    {
        if (std::fabs(scan.mz[index] - targetMz) < 1e-8)
        {
            scan.mz.erase(
                scan.mz.begin() +
                static_cast<std::ptrdiff_t>(index));
            scan.intensity.erase(
                scan.intensity.begin() +
                static_cast<std::ptrdiff_t>(index));
            scan.charge.erase(
                scan.charge.begin() +
                static_cast<std::ptrdiff_t>(index));
        }
    }
}

void checkMissingPeaksAndInterference()
{
    const TargetCase &carbon = targetCases()[0];
    Fixture missing = makeFixture(carbon, 0.5, 2);
    eraseNominalPeak(missing, 16);
    eraseNominalPeak(missing, 20);
    std::vector<isotopicPeak> peaks;
    const auto missingResult =
        fitFixture(missing, &peaks);
    const auto hasIndex = [&](int nominalIndex)
    {
        return std::any_of(
            peaks.begin(), peaks.end(),
            [&](const isotopicPeak &peak)
            {
                return peak.isotopeIndex ==
                       nominalIndex;
            });
    };
    check(hasIndex(15) && hasIndex(21),
          "finder stopped at an interior gap");
    check(!hasIndex(16) && !hasIndex(20),
          "finder invented missing peaks");
    check(std::fabs(
              missingResult.abundancePct - 50.0) <= 0.6,
          "missing peaks biased C13 fit");

    Fixture noisy = makeFixture(carbon, 0.5, 2);
    auto &scan = noisy.ms1.scans.front();
    const int interferenceBins[7]{
        0, 1, 2, 3, 4, 5, 6};
    const int weights[7]{
        1, 6, 15, 20, 15, 6, 1};
    for (size_t i = 0; i < 7; ++i)
    {
        const int nominalIndex = interferenceBins[i];
        const double mz =
            fixtureMzAt(noisy, nominalIndex);
        bool merged = false;
        for (size_t peakIndex = 0;
             peakIndex < scan.mz.size();
             ++peakIndex)
        {
            if (std::fabs(
                    scan.mz[peakIndex] - mz) < 1e-8)
            {
                scan.intensity[peakIndex] +=
                    weights[i] * 5e7 / 64.0;
                merged = true;
                break;
            }
        }
        if (!merged)
        {
            scan.mz.push_back(mz);
            scan.intensity.push_back(
                weights[i] * 5e7 / 64.0);
            scan.charge.push_back(noisy.charge);
        }
    }
    sortScan(scan);
    const auto noiseResult = fitFixture(noisy);
    check(std::fabs(
              noiseResult.abundancePct - 50.0) <= 0.5,
          "low-mass interference captured C13 fit");

    Fixture wrongCharge = makeFixture(carbon, 0.5, 2);
    auto &chargeScan = wrongCharge.ms1.scans.front();
    const size_t originalSize = chargeScan.mz.size();
    for (size_t i = 0; i < originalSize; ++i)
    {
        chargeScan.mz.push_back(chargeScan.mz[i]);
        chargeScan.intensity.push_back(
            chargeScan.intensity[i] * 10.0);
        chargeScan.charge.push_back(3);
    }
    sortScan(chargeScan);
    std::vector<isotopicPeak> chargePeaks;
    const auto chargeResult =
        fitFixture(wrongCharge, &chargePeaks);
    check(std::all_of(
              chargePeaks.begin(), chargePeaks.end(),
              [](const isotopicPeak &peak)
              { return peak.charge == 2; }),
          "wrong-charge peaks entered envelope");
    check(std::fabs(
              chargeResult.abundancePct - 50.0) <= 0.35,
          "wrong-charge interference biased fit");
}

void checkLargeHydrogenEnvelope()
{
    TargetCase hydrogen = targetCases()[1];
    hydrogen.peptide = std::string(20, 'W');
    Fixture fixture = makeFixture(hydrogen, 0.5, 3);
    std::vector<isotopicPeak> peaks;
    const auto result = fitFixture(fixture, &peaks);
    check(result.valid && !peaks.empty(),
          "large H2 envelope was invalid");
    check(peaks.back().isotopeIndex -
                  peaks.front().isotopeIndex >
              40,
          "large H2 envelope was truncated");
    check(std::fabs(
              result.abundancePct - 50.0) <= 0.35,
          "large H2 envelope missed 50%");
}

void checkInvalidEvidenceAndUnsupportedTarget()
{
    const TargetCase &carbon = targetCases()[0];
    Fixture fixture = makeFixture(carbon, 0.5, 2);
    const int isolatedIndex =
        static_cast<int>(
            fixture.distribution.probability.size()) +
        20;
    const isotopicPeak isolated{
        fixture.baseMz + 100.0,
        fixture.charge,
        1e6,
        isolatedIndex};
    const auto invalid =
        PSMfeatureExtractor::
            getSIPelementAbundanceFromMS1Peaks(
                {isolated},
                fixture.baseMass,
                fixture.peptide,
                fixture.charge,
                carbon.label,
                fixture.expectedPct);
    check(!invalid.valid,
          "incompatible peak was marked valid");
    check(invalid.fitScore == 0.0,
          "incompatible peak received a positive fit score");
    check(invalid.rawIsotopicPeakCount == 1 &&
              invalid.isotopicPeakCount == 1,
          "invalid fit discarded its raw peak count");
    check(std::fabs(
              invalid.abundancePct -
              fixture.expectedPct) > 1.0,
          "invalid fit echoed initializer");

    const auto unsupported =
        PSMfeatureExtractor::
            getSIPelementAbundanceFromMS1Peaks(
                {isolated},
                fixture.baseMass,
                fixture.peptide,
                fixture.charge,
                "S33",
                fixture.expectedPct);
    check(!unsupported.valid &&
              unsupported.fitScore == 0.0 &&
              unsupported.isotopicPeakCount == 0 &&
              unsupported.rawIsotopicPeakCount == 1,
          "unsupported isotope entered fitter");
}

void checkRawCountAndFitScoreFeatures()
{
    const TargetCase &carbon = targetCases()[0];

    Fixture sparse = makeFixture(carbon, 0.5, 2);
    auto &sparseScan = sparse.ms1.scans.front();
    const double anchorMz =
        fixtureMzAt(sparse, sparse.modeIndex);
    for (size_t index = sparseScan.mz.size();
         index-- > 0;)
    {
        if (std::fabs(sparseScan.mz[index] - anchorMz) <
            1e-8)
            continue;
        sparseScan.mz.erase(
            sparseScan.mz.begin() +
            static_cast<std::ptrdiff_t>(index));
        sparseScan.intensity.erase(
            sparseScan.intensity.begin() +
            static_cast<std::ptrdiff_t>(index));
        sparseScan.charge.erase(
            sparseScan.charge.begin() +
            static_cast<std::ptrdiff_t>(index));
    }
    check(sparseScan.mz.size() == 1,
          "failed to create sparse one-peak fixture");

    const sipPSM sparsePsm =
        extractFixtureFeatures(sparse);
    check(sparsePsm.isotopicPeakss[0].size() == 1 &&
              sparsePsm.isotopicPeakNumbers[0] == 1,
          "invalid MS1 fit did not preserve raw peak count");
    check(sparsePsm.MS1IsotopeFitScores[0] == 0.0,
          "one-peak MS1 fit received a positive score");

    Fixture complete = makeFixture(carbon, 0.5, 2);
    const sipPSM completePsm =
        extractFixtureFeatures(complete);
    check(completePsm.isotopicPeakNumbers[0] ==
              static_cast<int>(
                  completePsm.isotopicPeakss[0].size()),
          "valid MS1 fit did not report raw peak count");
    check(completePsm.isotopicPeakNumbers[0] > 1 &&
              completePsm.MS1IsotopeFitScores[0] >=
                  PSMfeatureExtractor::MinMs1IsotopeFitScore &&
              completePsm.MS1IsotopeFitScores[0] <= 1.0,
          "complete MS1 envelope received an invalid score");
}


void checkFollowingScanNeighborhood()
{
    const TargetCase &carbon = targetCases()[0];
    Fixture fixture = makeFixture(carbon, 0.5, 2);
    check(fixture.ms1.scans.size() == 1,
          "following-neighborhood fixture has unexpected scan count");

    sipros::RaxportMs1Scan signalScan =
        std::move(fixture.ms1.scans.front());
    signalScan.scanNumber = 101;
    sipros::RaxportMs1Scan parentScan;
    parentScan.scanNumber = 100;
    fixture.ms1.scans.clear();
    fixture.ms1.scans.push_back(std::move(parentScan));
    fixture.ms1.scans.push_back(std::move(signalScan));
    fixture.ms1.scanNumberToIndex.clear();
    fixture.ms1.scanNumberToIndex[100] = 0;
    fixture.ms1.scanNumberToIndex[101] = 1;

    int scanNumber = 100;
    const double matchedMz = fixtureMzAt(fixture, fixture.modeIndex);
    const auto tolerance = [&](double)
    { return 0.01 / fixture.charge; };
    const std::vector<isotopicPeak> peaks =
        PSMfeatureExtractor::findMs1IsotopicPeaks(
            &fixture.ms1,
            scanNumber,
            fixture.charge,
            fixture.baseMz,
            matchedMz,
            fixture.composition,
            fixture.target->label,
            fixture.expectedPct,
            tolerance);
    check(!peaks.empty() && scanNumber == 101,
          "MS1 feature extraction did not find a following-scan candidate");
}

void checkProductConvolvedPrecursorLookup()
{
    Fixture fixture = makeFixture(targetCases()[0], 0.5, 2);
    auto envelope = fixturePrecursorEnvelope(fixture);
    int scanNumber = 100;
    const double matchedMz = fixtureMzAt(fixture, fixture.modeIndex);
    const auto tolerance = [&](double)
    { return 0.01 / fixture.charge; };
    const std::vector<isotopicPeak> peaks =
        PSMfeatureExtractor::findMs1IsotopicPeaksFromEnvelope(
            &fixture.ms1,
            scanNumber,
            fixture.charge,
            fixture.baseMass,
            matchedMz,
            envelope.first,
            envelope.second,
            tolerance);
    check(!peaks.empty(),
          "product-convolved precursor envelope found no MS1 peaks");
    check(std::any_of(
              peaks.begin(), peaks.end(), [&](const isotopicPeak &peak)
              { return peak.isotopeIndex == fixture.modeIndex; }),
          "product-convolved precursor envelope lost its matched apex");
    const auto result =
        PSMfeatureExtractor::getSIPelementAbundanceFromMS1PeaksWithEnvelope(
            peaks,
            fixture.baseMass,
            fixture.peptide,
            fixture.charge,
            fixture.target->label,
            fixture.expectedPct,
            envelope.first,
            envelope.second);
    check(result.valid &&
              std::fabs(result.abundancePct - fixture.expectedPct) <= 0.35,
          "product-convolved precursor envelope produced the wrong MS1 abundance");
}


void checkEndpointStateIndependence()
{
    const TargetCase &oxygen = targetCases()[3];
    check(ProNovoConfig::applySipAbundance('O', 1.0),
          "failed to apply global O18 endpoint");
    check(ProNovoConfig::configIsotopologue
                  .vAtomIsotopicDistribution[2]
                  .vProb[2] == 1.0,
          "global O18 endpoint was not applied");
    Fixture oxygenMidpoint =
        makeFixture(oxygen, 0.5, 2);
    const auto oxygenResult = fitFixture(oxygenMidpoint);
    check(oxygenResult.valid &&
              std::fabs(
                  oxygenResult.abundancePct - 50.0) <= 0.35,
          "O18 endpoint state contaminated a later 50% fit");

    const TargetCase &sulfur = targetCases()[4];
    check(ProNovoConfig::applySipAbundance('S', 1.0),
          "failed to apply global S34 endpoint");
    check(ProNovoConfig::configIsotopologue
                  .vAtomIsotopicDistribution[5]
                  .vProb[2] == 1.0,
          "global S34 endpoint was not applied");
    Fixture sulfurMidpoint =
        makeFixture(sulfur, 0.5, 2);
    const auto sulfurResult = fitFixture(sulfurMidpoint);
    check(sulfurResult.valid &&
              std::fabs(
                  sulfurResult.abundancePct - 50.0) <= 0.35,
          "S34 endpoint state contaminated a later 50% fit");

    const double naturalCarbon =
        ProNovoConfig::getNaturalAtomIsotopeProbabilities(0)[1];
    check(ProNovoConfig::applySipAbundance(
              'C', naturalCarbon),
          "failed to restore natural isotope state");
}


int main(int argc, char **argv)
{
    try
    {
        (void)argv;
        check(argc == 1,
              "usage: ms1_isotope_envelope_test");
        check(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
              "failed to load built-in SIP profile");

        checkWhitelistAndConfiguredMasses();
        checkModalPrecursorAndWindows();
        checkMassErrorUsesPeptideSpecificSpacing();
        checkAllTargetsProduceModalEstimates();
        checkParallelModalEstimates();
        checkDirectSearchUsesFullPtmComposition();
        checkTargetElementCategories();
        checkCleanEnrichments();
        checkBiasedInitializersAndCharges();
        checkMissingPeaksAndInterference();
        checkLargeHydrogenEnvelope();
        checkInvalidEvidenceAndUnsupportedTarget();
        checkRawCountAndFitScoreFeatures();
        checkMissingPrecursorUsesIsolationCenter();
        checkFollowingScanNeighborhood();
        checkProductConvolvedPrecursorLookup();
        checkEndpointStateIndependence();

        std::cout
            << "ok: C13/H2/N15/O18/S34 exact-mass "
               "MS1 envelope fits pass 0-100%"
            << std::endl;
        return EXIT_SUCCESS;
    }
    catch (const std::exception &error)
    {
        std::cerr << "FAIL: " << error.what()
                  << std::endl;
        return EXIT_FAILURE;
    }
}
