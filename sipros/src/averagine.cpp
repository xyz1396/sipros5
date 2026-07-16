#include "averagine.h"
#include "proNovoConfig.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{

struct IsotopeShiftEstimate
{
    double exactMass = 0.0;
    int nominalShift = 0;
};

sipros::SourcedComposition scaledComposition(
    const sipros::SourcedComposition &composition,
    int factor)
{
    sipros::SourcedComposition result;
    for (size_t source = 0; source < sipros::IsotopeSourceCount; ++source)
        for (size_t element = 0; element < sipros::ElementCount; ++element)
            result.atoms[source][element] =
                composition.atoms[source][element] * factor;
    return result;
}

// Return an exact mode of the multinomial isotope-count distribution. The
// log-probability is separable and concave, so repeatedly assigning an atom
// to the largest p_i / (count_i + 1) marginal gives a global mode. If two
// allocations are tied, prefer the lower isotope index; this is the same
// lower-mass tie break used by IsotopeDistribution::getMostAbundantMass().
std::vector<int> multinomialMode(
    int atomCount,
    const IsotopeDistribution &distribution)
{
    const size_t isotopeCount = std::min(
        distribution.vMass.size(), distribution.vProb.size());
    std::vector<int> counts(isotopeCount, 0);
    for (int atom = 0; atom < atomCount && isotopeCount > 0; ++atom)
    {
        size_t best = isotopeCount;
        for (size_t isotope = 0; isotope < isotopeCount; ++isotope)
        {
            const double probability = distribution.vProb[isotope];
            if (!(probability > 0.0))
                continue;
            if (best == isotopeCount)
            {
                best = isotope;
                continue;
            }

            // Compare p_i/(count_i+1) without division. A small tolerance
            // preserves the deterministic lower-mass choice at exact ties.
            const double left = probability *
                static_cast<double>(counts[best] + 1);
            const double right = distribution.vProb[best] *
                static_cast<double>(counts[isotope] + 1);
            const double tolerance =
                16.0 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::fabs(left), std::fabs(right)});
            if (left > right + tolerance ||
                (std::fabs(left - right) <= tolerance && isotope < best))
                best = isotope;
        }
        if (best == isotopeCount)
            break;
        ++counts[best];
    }
    return counts;
}

void addIsotopeCount(IsotopeShiftEstimate &estimate,
                     int count,
                     const IsotopeDistribution &distribution,
                     size_t isotope)
{
    if (count <= 0 || isotope == 0 ||
        isotope >= distribution.vMass.size())
        return;
    const double exactDelta =
        distribution.vMass[isotope] - distribution.vMass.front();
    estimate.exactMass += static_cast<double>(count) * exactDelta;
    estimate.nominalShift += count * static_cast<int>(std::lround(exactDelta));
}

void addAtomGroup(IsotopeShiftEstimate &estimate,
                  int atomCount,
                  const IsotopeDistribution &distribution)
{
    const size_t isotopeCount = std::min(
        distribution.vMass.size(), distribution.vProb.size());
    if (atomCount <= 0 || isotopeCount <= 1)
        return;
    const std::vector<int> counts =
        multinomialMode(atomCount, distribution);
    for (size_t isotope = 1; isotope < isotopeCount; ++isotope)
        addIsotopeCount(estimate, counts[isotope], distribution, isotope);
}

bool sameProbabilities(const IsotopeDistribution &left,
                       const IsotopeDistribution &right)
{
    if (left.vProb.size() != right.vProb.size() ||
        left.vMass != right.vMass)
        return false;
    for (size_t isotope = 0; isotope < left.vProb.size(); ++isotope)
    {
        if (std::fabs(left.vProb[isotope] - right.vProb[isotope]) >
            8.0 * std::numeric_limits<double>::epsilon())
            return false;
    }
    return true;
}

// Estimate the exact modal isotope-count vector while retaining the mass
// defect of every element. Thus one O18 or S34 atom always adds its full +2
// isotope delta; it is never represented as two averaged one-neutron events.
IsotopeShiftEstimate estimateIsotopeShift(
    const sipros::SourcedComposition &composition)
{
    const auto &active = ProNovoConfig::configIsotopologue
                             .vAtomIsotopicDistribution;
    const auto &natural = ProNovoConfig::configIsotopologue
                              .vNaturalAtomIsotopicDistribution;
    IsotopeShiftEstimate estimate;

    for (size_t element = 0; element < sipros::ElementCount; ++element)
    {
        if (element >= active.size() || element >= natural.size())
            continue;
        const int biosyntheticCount =
            composition[sipros::IsotopeSource::Biosynthetic][element];
        const int naturalCount =
            composition[sipros::IsotopeSource::ReagentNatural][element] +
            composition[sipros::IsotopeSource::DigestionSolvent][element];
        // When both sources have the same distribution (the regular/natural
        // case), combine their counts before finding the mode. Source
        // provenance must not split one natural binomial into smaller modes.
        if (sameProbabilities(active[element], natural[element]))
        {
            addAtomGroup(estimate,
                         biosyntheticCount + naturalCount,
                         natural[element]);
            continue;
        }

        // Under SIP enrichment only biosynthetic atoms use the active target
        // distribution. Reagent and digestion-solvent atoms remain one
        // combined natural-abundance group.
        addAtomGroup(estimate, biosyntheticCount, active[element]);
        addAtomGroup(estimate, naturalCount, natural[element]);
    }
    return estimate;
}

} // namespace

averagine::averagine(const int minimumPeptideLength,
                     const int maximumPeptideLength)
    : minPepLen(minimumPeptideLength),
      maxPepLen(maximumPeptideLength),
      pepLenRange(maximumPeptideLength - minimumPeptideLength + 1)
{
    const auto residue = ProNovoConfig::configIsotopologue
                             .mResidueSourcedComposition.find(
                                 averagineResidue);
    if (residue != ProNovoConfig::configIsotopologue
                       .mResidueSourcedComposition.end())
        averagineComposition = residue->second;
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
        averagineComposition, averagineSIPdistribution);
}

averagine::averagine()
{
}

averagine::~averagine() = default;

double averagine::lightMass(size_t element)
{
    const auto &atoms = ProNovoConfig::configIsotopologue
                            .vNaturalAtomIsotopicDistribution;
    if (element >= atoms.size() || atoms[element].vMass.empty())
        return 0.0;
    return atoms[element].vMass.front();
}

double averagine::baseMass(
    const sipros::SourcedComposition &composition)
{
    const sipros::AtomCounts total = composition.total();
    double mass = 0.0;
    for (size_t element = 0; element < total.size(); ++element)
        mass += static_cast<double>(total[element]) * lightMass(element);
    return mass;
}

void averagine::changeAtomSIPabundance(const char sipAtom,
                                       const double fraction)
{
    const char atom = static_cast<char>(std::toupper(
        static_cast<unsigned char>(sipAtom)));
    if (SIPatoms.find(atom) == string::npos ||
        !ProNovoConfig::applySipAbundance(atom, fraction))
    {
        cerr << atom << " element is not a supported SIP target." << endl;
        return;
    }
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
        averagineComposition, averagineSIPdistribution);
}

bool averagine::changeAtomProbability(
    std::vector<double> &probabilities,
    char atom,
    const double fraction)
{
    const double targetFraction =
        std::max(0.0, std::min(1.0, fraction));
    if (targetFraction != fraction)
    {
        cerr << "Warning: SIP abundance " << fraction
             << " is outside [0,1], clamped to "
             << targetFraction << "." << endl;
    }

    const size_t targetIsotopeIndex =
        (atom == 'O' || atom == 'S') ? 2u : 1u;
    if (probabilities.empty() ||
        targetIsotopeIndex >= probabilities.size())
    {
        cout << atom
             << " target isotope is unavailable in atom distribution!"
             << endl;
        return false;
    }

    double nonTargetTotal = 0.0;
    for (size_t isotope = 0; isotope < probabilities.size(); ++isotope)
        if (isotope != targetIsotopeIndex)
            nonTargetTotal += probabilities[isotope];

    if (!(nonTargetTotal > 0.0))
    {
        std::fill(probabilities.begin(), probabilities.end(), 0.0);
        probabilities[0] = 1.0 - targetFraction;
    }
    else
    {
        const double scale =
            (1.0 - targetFraction) / nonTargetTotal;
        for (size_t isotope = 0; isotope < probabilities.size(); ++isotope)
            if (isotope != targetIsotopeIndex)
                probabilities[isotope] *= scale;
    }
    probabilities[targetIsotopeIndex] = targetFraction;
    return true;
}

void averagine::calAveraginePepAtomCounts()
{
    averaginePepCompositions.clear();
    averaginePepCompositions.reserve(
        static_cast<size_t>(std::max(0, pepLenRange)));
    for (int peptideLength = minPepLen;
         peptideLength <= maxPepLen;
         ++peptideLength)
        averaginePepCompositions.push_back(
            scaledComposition(averagineComposition, peptideLength));
}

sipros::SourcedComposition *averagine::getAveraginePepComposition(
    const int peptideLength)
{
    return &averaginePepCompositions[static_cast<size_t>(
        peptideLength - minPepLen)];
}

void averagine::calAveraginePepSIPdistributions()
{
    averaginePepSIPdistributions.clear();
    averaginePepSIPdistributions.reserve(
        static_cast<size_t>(std::max(0, pepLenRange)));
    const auto nTerm = ProNovoConfig::configIsotopologue
                           .mResidueSourcedComposition.at("Nterm");
    const auto cTerm = ProNovoConfig::configIsotopologue
                           .mResidueSourcedComposition.at("Cterm");
    for (int peptideLength = minPepLen;
         peptideLength <= maxPepLen;
         ++peptideLength)
    {
        sipros::SourcedComposition composition =
            scaledComposition(averagineComposition, peptideLength);
        composition += nTerm;
        composition += cTerm;
        IsotopeDistribution distribution;
        ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
            composition, distribution);
        averaginePepSIPdistributions.push_back(std::move(distribution));
    }
}

IsotopeDistribution *averagine::getAveraginePepSIPdistribution(
    const int peptideLength)
{
    return &averaginePepSIPdistributions[static_cast<size_t>(
        peptideLength - minPepLen)];
}

void averagine::calPepAtomCounts(const string &peptideSequence)
{
    pepComposition = {};
    if (!ProNovoConfig::configIsotopologue.computeSourcedComposition(
            peptideSequence, pepComposition))
        return;
}

void averagine::calBYionsAtomCounts(const string &peptideSequence)
{
    BionsCompositions.clear();
    YionsCompositions.clear();
    BionsCompositions.reserve(peptideSequence.size());
    YionsCompositions.reserve(peptideSequence.size());

    sipros::SourcedComposition bIon;
    sipros::SourcedComposition yIon =
        ProNovoConfig::configIsotopologue
            .mResidueSourcedComposition.at("Nterm") +
        ProNovoConfig::configIsotopologue
            .mResidueSourcedComposition.at("Cterm");

    for (size_t index = 0; index < peptideSequence.size(); ++index)
    {
        const string symbol = peptideSequence.substr(index, 1);
        const auto residue = ProNovoConfig::configIsotopologue
                                 .mResidueSourcedComposition.find(symbol);
        if (residue == ProNovoConfig::configIsotopologue
                           .mResidueSourcedComposition.end())
        {
            cerr << "ERROR: cannot find " << symbol
                 << " residue or PTM in the built-in chemistry." << endl;
            BionsCompositions.clear();
            YionsCompositions.clear();
            return;
        }
        bIon += residue->second;
        if (std::isalpha(static_cast<unsigned char>(peptideSequence[index])))
            BionsCompositions.push_back(bIon);
        else if (!BionsCompositions.empty())
            BionsCompositions.back() = bIon;
    }
    if (!BionsCompositions.empty())
        BionsCompositions.pop_back();

    for (size_t reverse = peptideSequence.size(); reverse > 0; --reverse)
    {
        const size_t index = reverse - 1;
        const string symbol = peptideSequence.substr(index, 1);
        const auto residue = ProNovoConfig::configIsotopologue
                                 .mResidueSourcedComposition.find(symbol);
        if (residue == ProNovoConfig::configIsotopologue
                           .mResidueSourcedComposition.end())
            return;
        yIon += residue->second;
        if (std::isalpha(static_cast<unsigned char>(peptideSequence[index])))
            YionsCompositions.push_back(yIon);
    }
    if (!YionsCompositions.empty())
        YionsCompositions.pop_back();
}

double averagine::calNetronMass(const string &peptideSequence)
{
    calPepAtomCounts(peptideSequence);
    const IsotopeShiftEstimate estimate =
        estimateIsotopeShift(pepComposition);
    return estimate.nominalShift > 0
               ? estimate.exactMass /
                     static_cast<double>(estimate.nominalShift)
               : ProNovoConfig::getNeutronMass();
}

double averagine::calPrecursorBaseMass(const string &peptideSequence)
{
    calPepAtomCounts(peptideSequence);
    return baseMass(pepComposition);
}

void averagine::calBYionBaseMasses(const string &peptideSequence)
{
    calBYionsAtomCounts(peptideSequence);
    BionsBaseMasses.clear();
    YionsBaseMasses.clear();
    BionsBaseMasses.reserve(BionsCompositions.size());
    YionsBaseMasses.reserve(YionsCompositions.size());
    for (const auto &composition : BionsCompositions)
        BionsBaseMasses.push_back(baseMass(composition));
    for (const auto &composition : YionsCompositions)
        YionsBaseMasses.push_back(baseMass(composition));
}

double averagine::calPrecursorMass(const string &peptideSequence)
{
    calPepAtomCounts(peptideSequence);
    const IsotopeShiftEstimate estimate =
        estimateIsotopeShift(pepComposition);
    return baseMass(pepComposition) + estimate.exactMass;
}

void averagine::calDiffAtomCounts(const string &peptideSequence)
{
    // Retained for source compatibility. Exact sourced convolution below no
    // longer needs the old flat averagine-difference approximation.
    calPepAtomCounts(peptideSequence);
}

void averagine::calPrecursorIsotopeDistribution(
    const string &peptideSequence,
    IsotopeDistribution &distribution)
{
    calPepAtomCounts(peptideSequence);
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
        pepComposition, distribution);
}
