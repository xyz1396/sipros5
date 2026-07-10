#include "PSMfeatureExtractor.h"
#include "RaxportHdf5Reader.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdexcept>

PSMfeatureExtractor::PSMfeatureExtractor() : mAveragine(averagine()), mSipPSM(nullptr)
{}

struct SupportedSipIsotope
{
    char atom = '\0';
    int atomIndex = -1;
    int isotopeIndex = -1;
    int massNumber = 0;
    int nominalShift = 0;
    const char *canonical = "";
};

static std::string normalizedSipLabel(const std::string &value)
{
    std::string normalized;
    normalized.reserve(value.size());
    for (unsigned char c : value)
    {
        if (!std::isspace(c))
            normalized.push_back(static_cast<char>(std::toupper(c)));
    }
    return normalized;
}

static bool resolveSupportedSipIsotope(const std::string &value,
                                       SupportedSipIsotope &spec,
                                       int isotopeMassNumber = -1)
{
    const std::string label = normalizedSipLabel(value);
    if (label == "C" || label == "C13")
        spec = {'C', 0, 1, 13, 1, "C13"};
    else if (label == "H" || label == "H2")
        spec = {'H', 1, 1, 2, 1, "H2"};
    else if (label == "N" || label == "N15")
        spec = {'N', 3, 1, 15, 1, "N15"};
    else if (label == "O" || label == "O18")
        spec = {'O', 2, 2, 18, 2, "O18"};
    else if (label == "S" || label == "S34")
        spec = {'S', 5, 2, 34, 2, "S34"};
    else
        return false;

    return isotopeMassNumber <= 0 || isotopeMassNumber == spec.massNumber;
}

static std::string peptideSequenceForComposition(const std::string &peptideSeq)
{
    const std::size_t start = peptideSeq.find_first_of('[');
    const std::size_t end = peptideSeq.find_last_of(']');
    if (start != std::string::npos && end != std::string::npos && end > start)
    {
        return peptideSeq.substr(start + 1, end - start - 1);
    }

    const std::size_t firstDot = peptideSeq.find_first_of('.');
    const std::size_t lastDot = peptideSeq.find_last_of('.');
    if (firstDot != std::string::npos && lastDot != std::string::npos && lastDot > firstDot)
    {
        return peptideSeq.substr(firstDot + 1, lastDot - firstDot - 1);
    }

    std::string stripped;
    stripped.reserve(peptideSeq.size());
    for (unsigned char c : peptideSeq)
    {
        if (c != '[' && c != ']' && c != '.' && c != '-' && !std::isspace(c))
        {
            stripped.push_back(static_cast<char>(c));
        }
    }
    return stripped;
}

std::pair<int, int> PSMfeatureExtractor::getSeqLengthAndMissCleavageSiteNumber(const std::string &peptideSeq)
{
    std::size_t start = peptideSeq.find_first_of('[');
    std::size_t end = peptideSeq.find_last_of(']');
    std::string seq;
    if (start == std::string::npos || end == std::string::npos || end <= start + 1)
    {
        for (char c : peptideSeq)
            if (std::isalpha(static_cast<unsigned char>(c)))
                seq.push_back(c);
    }
    else
    {
        seq = peptideSeq.substr(start + 1, end - start - 1);
    }
    std::string aminoAcids;
    aminoAcids.reserve(seq.size());
    for (char c : seq)
    {
        if (std::isalpha(static_cast<unsigned char>(c)))
            aminoAcids.push_back(c);
    }

    int count = 0;
    for (size_t i = 1; i + 1 < aminoAcids.size(); ++i)
    {
        if (std::find(cleavageSites.begin(), cleavageSites.end(), aminoAcids[i]) != cleavageSites.end())
            ++count;
    }
    return {static_cast<int>(aminoAcids.size()), count};
}

int PSMfeatureExtractor::getPTMnumber(const std::string &peptideSeq)
{
    int count = 0;
    for (unsigned char c : peptideSeq)
    {
        if (std::isalpha(c) || c == '[' || c == ']' || c == '-' || c == '.')
        {
            continue;
        }
        count++;
    }
    return count;
}

std::pair<int, double> PSMfeatureExtractor::getMassWindowShiftAndError(const double observedPrecursorMass,
                                                                       const double calculatedPrecursorMass)
{
    int massWindowShift = static_cast<int>(round(std::abs(observedPrecursorMass - calculatedPrecursorMass) /
                                                 ProNovoConfig::getNeutronMass()));
    double massError = std::fmod(std::abs(observedPrecursorMass - calculatedPrecursorMass),
                                 ProNovoConfig::getNeutronMass());
    if (massError > ProNovoConfig::getNeutronMass() / 2)
    {
        massError = ProNovoConfig::getNeutronMass() - massError;
    }
    // convert it to ppm
    massError = massError / calculatedPrecursorMass * 1000000;
    return {massWindowShift, massError};
}


std::string PSMfeatureExtractor::peptideBodyWithPtms(const std::string &decorated)
{
    std::string out;
    out.reserve(decorated.size());
    for (char c : decorated)
    {
        if (c != '[' && c != ']')
            out.push_back(c);
    }
    return out;
}

int PSMfeatureExtractor::countMissCleavage(const std::string &naked)
{
    if (naked.size() <= 1)
        return 0;
    int n = 0;
    for (size_t i = 1; i + 1 < naked.size(); ++i)
    {
        if (naked[i] == 'K' || naked[i] == 'R')
            ++n;
    }
    return n;
}

int PSMfeatureExtractor::countPTM(const std::string &decorated)
{
    int n = 0;
    for (char c : decorated)
    {
        if (!std::isalpha(static_cast<unsigned char>(c)) && c != '[' && c != ']')
            ++n;
    }
    return n;
}

std::string PSMfeatureExtractor::canonicalSipIsotope(const std::string &sipAtom,
                                                         int isotopeMassNumber)
{
    SupportedSipIsotope spec;
    return resolveSupportedSipIsotope(sipAtom, spec, isotopeMassNumber)
               ? std::string(spec.canonical)
               : std::string();
}

int PSMfeatureExtractor::sipAtomIndex(const std::string &sipAtom)
{
    SupportedSipIsotope spec;
    return resolveSupportedSipIsotope(sipAtom, spec) ? spec.atomIndex : -1;
}

int PSMfeatureExtractor::sipIsotopeIndex(const std::string &sipAtom)
{
    SupportedSipIsotope spec;
    return resolveSupportedSipIsotope(sipAtom, spec) ? spec.isotopeIndex : -1;
}

int PSMfeatureExtractor::sipNominalShiftPerAtom(const std::string &sipAtom)
{
    SupportedSipIsotope spec;
    return resolveSupportedSipIsotope(sipAtom, spec) ? spec.nominalShift : 0;
}

double PSMfeatureExtractor::expectedNaturalNominalShiftExceptTarget(
    const std::array<int, 6> &atomCounts,
    int targetAtomIndex,
    int targetIsotopeIndex,
    double targetFraction)
{
    double expectedShift = 0.0;
    const auto &atomDistributions =
        ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution;
    const size_t atomCount = std::min(atomCounts.size(), atomDistributions.size());
    static const std::string atomLetters = "CHONPS";

    for (size_t atomIdx = 0; atomIdx < atomCount; ++atomIdx)
    {
        if (atomCounts[atomIdx] <= 0)
            continue;

        std::vector<double> probabilities =
            ProNovoConfig::getNaturalAtomIsotopeProbabilities(atomIdx);
        if (static_cast<int>(atomIdx) == targetAtomIndex &&
            std::isfinite(targetFraction))
        {
            if (!averagine::changeAtomProbability(
                    probabilities, atomLetters[atomIdx], targetFraction))
                return std::numeric_limits<double>::quiet_NaN();
        }

        const auto &masses = atomDistributions[atomIdx].vMass;
        const size_t isotopeCount = std::min(probabilities.size(), masses.size());
        for (size_t isotopeIdx = 1; isotopeIdx < isotopeCount; ++isotopeIdx)
        {
            if (static_cast<int>(atomIdx) == targetAtomIndex &&
                static_cast<int>(isotopeIdx) == targetIsotopeIndex)
                continue;
            const int nominalShift = static_cast<int>(std::lround(
                masses[isotopeIdx] - masses[0]));
            expectedShift += static_cast<double>(atomCounts[atomIdx]) *
                             nominalShift * probabilities[isotopeIdx];
        }
    }
    return expectedShift;
}

int PSMfeatureExtractor::ms1PeakCharge(const sipros::RaxportMs1Scan &scan, size_t idx)
{
    if (idx < scan.charge.size())
        return scan.charge[idx];
    return 0;
}

size_t PSMfeatureExtractor::findMs1Peak(const sipros::RaxportMs1Scan &scan,
                                        double targetMz,
                                        const std::function<double(double)> &mzToleranceDaAt,
                                        int requiredCharge)
{
    if (scan.mz.empty())
        return std::numeric_limits<size_t>::max();

    size_t best = std::numeric_limits<size_t>::max();
    double bestIntensity = 0.0;
    const double mzTolerance = mzToleranceDaAt(targetMz);
    auto first = std::lower_bound(scan.mz.begin(), scan.mz.end(), targetMz - mzTolerance);
    for (auto it = first; it != scan.mz.end() && *it <= targetMz + mzTolerance; ++it)
    {
        const size_t idx = static_cast<size_t>(it - scan.mz.begin());
        if (requiredCharge >= 0 && ms1PeakCharge(scan, idx) != requiredCharge)
            continue;
        if (std::fabs(scan.mz[idx] - targetMz) <= mzTolerance &&
            scan.intensity[idx] > bestIntensity)
        {
            bestIntensity = scan.intensity[idx];
            best = idx;
        }
    }
    return best;
}

static size_t findEnvelopePeak(const sipros::RaxportMs1Scan &scan,
                               double targetMz,
                               const std::function<double(double)> &mzToleranceDaAt,
                               int precursorCharge)
{
    size_t idx = PSMfeatureExtractor::findMs1Peak(
        scan, targetMz, mzToleranceDaAt, precursorCharge);
    if (idx == std::numeric_limits<size_t>::max())
    {
        // Raxport leaves charge unknown for many low-intensity envelope peaks.
        // Accept those, but never substitute a peak assigned to another charge.
        idx = PSMfeatureExtractor::findMs1Peak(scan, targetMz, mzToleranceDaAt, 0);
    }
    return idx;
}

struct NominalEnvelopeModel
{
    std::vector<double> probability;
    std::vector<double> massDelta;
};

static NominalEnvelopeModel theoreticalSipEnvelope(
    const std::array<int, 6> &atomCounts,
    const SupportedSipIsotope &spec,
    double targetFraction);

static int isotopeHalfWindow(int targetAtomCount,
                             int targetNominalShift,
                             int matchedIsotopeIndex)
{
    const int atomCount = std::max(1, targetAtomCount);
    const int nominalShift = std::max(1, targetNominalShift);
    const double fullyLabeledShift = static_cast<double>(atomCount * nominalShift);
    const double anchorFraction = fullyLabeledShift > 0.0
                                      ? std::max(0.0, std::min(1.0,
                                            matchedIsotopeIndex / fullyLabeledShift))
                                      : 0.0;
    // The matched precursor supplies an enrichment-independent location. Use
    // the binomial variance at that location: narrow near 0/100%, widest at
    // 50%. Six sigma plus a calibration/tail margin covers the full envelope
    // without sampling a worst-case-wide field of unrelated MS1 peaks.
    const double variance = atomCount * anchorFraction * (1.0 - anchorFraction) *
                            nominalShift * nominalShift;
    const double sixSigma = 6.0 * std::sqrt(std::max(0.0, variance));
    const int minimumWindow = std::max(6, 3 * nominalShift + 2);
    return std::max(minimumWindow, static_cast<int>(std::ceil(sixSigma)) + 4);
}

std::vector<isotopicPeak> PSMfeatureExtractor::findMs1IsotopicPeaks(
    const sipros::RaxportMs1Data *ms1Data,
    int &ms1ScanNumber,
    int precursorCharge,
    double monoPrecursorMz,
    double matchedPrecursorMz,
    const std::array<int, 6> &atomCounts,
    const std::string &sipAtom,
    double expectedEnrichmentPct,
    const std::function<double(double)> &mzToleranceDaAt)
{
    std::vector<isotopicPeak> peaks;
    SupportedSipIsotope spec;
    if (!ms1Data || precursorCharge <= 0 ||
        !resolveSupportedSipIsotope(sipAtom, spec) ||
        spec.atomIndex < 0 ||
        atomCounts[static_cast<size_t>(spec.atomIndex)] <= 0 ||
        !std::isfinite(expectedEnrichmentPct))
        return peaks;

    auto scanIt = ms1Data->scanNumberToIndex.find(ms1ScanNumber);
    if (scanIt == ms1Data->scanNumberToIndex.end())
        return peaks;

    const double targetFraction =
        std::max(0.0, std::min(1.0, expectedEnrichmentPct / 100.0));
    const NominalEnvelopeModel model =
        theoreticalSipEnvelope(atomCounts, spec, targetFraction);
    if (model.probability.empty() ||
        model.massDelta.size() != model.probability.size())
        return peaks;

    const double modelMaximum =
        *std::max_element(model.probability.begin(), model.probability.end());
    const double assignmentFloor = modelMaximum * 1e-12;
    const double observedNeutralDelta =
        (matchedPrecursorMz - monoPrecursorMz) * precursorCharge;
    int assignedIndex = -1;
    double closestDelta = std::numeric_limits<double>::infinity();
    for (size_t index = 0; index < model.probability.size(); ++index)
    {
        if (model.probability[index] < assignmentFloor ||
            !std::isfinite(model.massDelta[index]))
            continue;
        const double delta =
            std::fabs(model.massDelta[index] - observedNeutralDelta);
        if (delta < closestDelta)
        {
            closestDelta = delta;
            assignedIndex = static_cast<int>(index);
        }
    }
    if (assignedIndex < 0)
        return peaks;

    size_t scanIdx = scanIt->second;
    size_t anchorIdx = std::numeric_limits<size_t>::max();
    const sipros::RaxportMs1Scan *anchorScan = nullptr;
    for (int attempt = 0; attempt < 3; ++attempt)
    {
        const sipros::RaxportMs1Scan &scan = ms1Data->scans[scanIdx];
        anchorIdx = findEnvelopePeak(
            scan, matchedPrecursorMz, mzToleranceDaAt, precursorCharge);
        if (anchorIdx != std::numeric_limits<size_t>::max())
        {
            ms1ScanNumber = scan.scanNumber;
            anchorScan = &scan;
            break;
        }
        if (scanIdx == 0)
            break;
        --scanIdx;
    }
    if (!anchorScan)
        return peaks;

    const sipros::RaxportMs1Scan &scan = *anchorScan;
    const double anchorModelMz =
        monoPrecursorMz +
        model.massDelta[static_cast<size_t>(assignedIndex)] / precursorCharge;
    const double anchorResidual = scan.mz[anchorIdx] - anchorModelMz;
    const int targetAtomCount =
        atomCounts[static_cast<size_t>(spec.atomIndex)];
    const int halfWindow = isotopeHalfWindow(
        targetAtomCount, spec.nominalShift, assignedIndex);
    const int firstIsotopeIndex = std::max(0, assignedIndex - halfWindow);
    const int lastIsotopeIndex = std::min(
        static_cast<int>(model.probability.size()) - 1,
        assignedIndex + halfWindow);
    peaks.reserve(static_cast<size_t>(
        std::max(0, lastIsotopeIndex - firstIsotopeIndex + 1)));

    const double extractionFloor = modelMaximum * 1e-12;
    for (int isotopeIndex = firstIsotopeIndex;
         isotopeIndex <= lastIsotopeIndex;
         ++isotopeIndex)
    {
        const size_t modelIndex = static_cast<size_t>(isotopeIndex);
        if (model.probability[modelIndex] < extractionFloor ||
            !std::isfinite(model.massDelta[modelIndex]))
            continue;
        const double expectedMz =
            monoPrecursorMz + model.massDelta[modelIndex] / precursorCharge +
            anchorResidual;
        const size_t idx = findEnvelopePeak(
            scan, expectedMz, mzToleranceDaAt, precursorCharge);
        if (idx == std::numeric_limits<size_t>::max())
            continue;
        peaks.push_back({scan.mz[idx],
                         ms1PeakCharge(scan, idx),
                         scan.intensity[idx],
                         isotopeIndex});
    }

    const bool hasAnchor = std::any_of(
        peaks.begin(), peaks.end(), [&](const isotopicPeak &peak)
        { return peak.isotopeIndex == assignedIndex; });
    if (!hasAnchor)
    {
        peaks.push_back({scan.mz[anchorIdx],
                         ms1PeakCharge(scan, anchorIdx),
                         scan.intensity[anchorIdx],
                         assignedIndex});
    }

    std::sort(peaks.begin(), peaks.end(),
              [](const isotopicPeak &a, const isotopicPeak &b)
              { return a.isotopeIndex < b.isotopeIndex; });
    peaks.erase(
        std::unique(peaks.begin(), peaks.end(),
                    [](const isotopicPeak &a, const isotopicPeak &b)
                    { return a.isotopeIndex == b.isotopeIndex; }),
        peaks.end());
    return peaks;
}

static std::vector<double> binomialProbabilities(int atomCount,
                                                    double targetFraction)
{
    std::vector<double> distribution(
        static_cast<size_t>(std::max(0, atomCount) + 1), 0.0);
    if (atomCount <= 0)
    {
        distribution[0] = 1.0;
        return distribution;
    }

    const double p = std::max(0.0, std::min(1.0, targetFraction));
    if (p <= 0.0)
    {
        distribution[0] = 1.0;
        return distribution;
    }
    if (p >= 1.0)
    {
        distribution[static_cast<size_t>(atomCount)] = 1.0;
        return distribution;
    }

    const int mode = std::min(
        atomCount, static_cast<int>(std::floor((atomCount + 1) * p)));
    const double logModeProbability =
        std::lgamma(atomCount + 1.0) - std::lgamma(mode + 1.0) -
        std::lgamma(atomCount - mode + 1.0) +
        mode * std::log(p) + (atomCount - mode) * std::log1p(-p);
    distribution[static_cast<size_t>(mode)] =
        std::exp(logModeProbability);

    for (int k = mode; k > 0; --k)
    {
        distribution[static_cast<size_t>(k - 1)] =
            distribution[static_cast<size_t>(k)] *
            static_cast<double>(k) /
            static_cast<double>(atomCount - k + 1) *
            (1.0 - p) / p;
    }
    for (int k = mode; k < atomCount; ++k)
    {
        distribution[static_cast<size_t>(k + 1)] =
            distribution[static_cast<size_t>(k)] *
            static_cast<double>(atomCount - k) /
            static_cast<double>(k + 1) *
            p / (1.0 - p);
    }

    const double sum =
        std::accumulate(distribution.begin(), distribution.end(), 0.0);
    if (sum > 0.0)
        for (double &probability : distribution)
            probability /= sum;
    return distribution;
}

struct EnvelopeComponent
{
    int nominalOffset = 0;
    std::vector<double> probability;
    std::vector<double> massMoment;
};

static bool normalizeAndTrim(EnvelopeComponent &component)
{
    if (component.probability.empty() ||
        component.probability.size() != component.massMoment.size())
        return false;

    const double maximum = *std::max_element(
        component.probability.begin(), component.probability.end());
    if (!(maximum > 0.0))
        return false;
    const double cutoff = maximum * 1e-14;

    size_t first = 0;
    while (first + 1 < component.probability.size() &&
           component.probability[first] < cutoff)
        ++first;
    size_t last = component.probability.size();
    while (last > first + 1 && component.probability[last - 1] < cutoff)
        --last;

    if (first > 0 || last < component.probability.size())
    {
        component.probability = std::vector<double>(
            component.probability.begin() + static_cast<std::ptrdiff_t>(first),
            component.probability.begin() + static_cast<std::ptrdiff_t>(last));
        component.massMoment = std::vector<double>(
            component.massMoment.begin() + static_cast<std::ptrdiff_t>(first),
            component.massMoment.begin() + static_cast<std::ptrdiff_t>(last));
        component.nominalOffset += static_cast<int>(first);
    }

    const double sum = std::accumulate(
        component.probability.begin(), component.probability.end(), 0.0);
    if (!(sum > 0.0))
        return false;
    for (size_t i = 0; i < component.probability.size(); ++i)
    {
        component.probability[i] /= sum;
        component.massMoment[i] /= sum;
    }
    return true;
}

static EnvelopeComponent targetElementEnvelope(
    int atomCount,
    const SupportedSipIsotope &spec,
    double targetFraction)
{
    EnvelopeComponent component;
    if (atomCount <= 0)
        return component;

    const auto &configured =
        ProNovoConfig::configIsotopologue
            .vAtomIsotopicDistribution[static_cast<size_t>(spec.atomIndex)];
    std::vector<double> probabilities =
        ProNovoConfig::getNaturalAtomIsotopeProbabilities(
            static_cast<size_t>(spec.atomIndex));
    if (!averagine::changeAtomProbability(
            probabilities, spec.atom, targetFraction))
        return component;
    if (probabilities.size() != configured.vMass.size() ||
        spec.isotopeIndex >= static_cast<int>(probabilities.size()))
        return component;

    bool binaryTarget = true;
    for (size_t isotope = 1; isotope < probabilities.size(); ++isotope)
    {
        if (static_cast<int>(isotope) != spec.isotopeIndex &&
            probabilities[isotope] > 0.0)
        {
            binaryTarget = false;
            break;
        }
    }

    const double targetMassDelta =
        configured.vMass[static_cast<size_t>(spec.isotopeIndex)] -
        configured.vMass[0];
    if (binaryTarget)
    {
        const std::vector<double> binomial =
            binomialProbabilities(atomCount, targetFraction);
        component.probability.assign(
            static_cast<size_t>(atomCount * spec.nominalShift + 1), 0.0);
        component.massMoment.assign(component.probability.size(), 0.0);
        for (int heavyCount = 0; heavyCount <= atomCount; ++heavyCount)
        {
            const size_t index =
                static_cast<size_t>(heavyCount * spec.nominalShift);
            const double probability =
                binomial[static_cast<size_t>(heavyCount)];
            component.probability[index] = probability;
            component.massMoment[index] =
                probability * heavyCount * targetMassDelta;
        }
        normalizeAndTrim(component);
        return component;
    }

    struct AtomicEntry
    {
        int nominalShift = 0;
        double probability = 0.0;
        double massDelta = 0.0;
    };
    std::vector<AtomicEntry> entries;
    int maximumAtomicShift = 0;
    for (size_t isotope = 0; isotope < probabilities.size(); ++isotope)
    {
        if (!(probabilities[isotope] > 0.0))
            continue;
        const double massDelta =
            configured.vMass[isotope] - configured.vMass[0];
        const int nominalShift =
            static_cast<int>(std::lround(massDelta));
        if (nominalShift < 0)
            continue;
        entries.push_back(
            {nominalShift, probabilities[isotope], massDelta});
        maximumAtomicShift =
            std::max(maximumAtomicShift, nominalShift);
    }
    if (entries.empty())
        return {};

    component.probability = {1.0};
    component.massMoment = {0.0};
    for (int atom = 0; atom < atomCount; ++atom)
    {
        std::vector<double> nextProbability(
            component.probability.size() +
                static_cast<size_t>(maximumAtomicShift),
            0.0);
        std::vector<double> nextMassMoment(
            nextProbability.size(), 0.0);
        for (size_t index = 0;
             index < component.probability.size();
             ++index)
        {
            for (const AtomicEntry &entry : entries)
            {
                const size_t nextIndex =
                    index + static_cast<size_t>(entry.nominalShift);
                nextProbability[nextIndex] +=
                    component.probability[index] * entry.probability;
                nextMassMoment[nextIndex] +=
                    component.massMoment[index] * entry.probability +
                    component.probability[index] * entry.probability *
                        entry.massDelta;
            }
        }
        component.probability = std::move(nextProbability);
        component.massMoment = std::move(nextMassMoment);
        if (!normalizeAndTrim(component))
            return {};
    }
    return component;
}

static EnvelopeComponent naturalBackgroundEnvelope(
    const std::array<int, 6> &atomCounts,
    int excludedAtomIndex,
    int maxNominalShift)
{
    const int maximum = std::max(0, maxNominalShift);
    std::vector<double> lambda(
        static_cast<size_t>(maximum + 1), 0.0);
    std::vector<double> massLambda(
        static_cast<size_t>(maximum + 1), 0.0);
    const auto &atomDistributions =
        ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution;
    const size_t atomLimit =
        std::min(atomCounts.size(), atomDistributions.size());

    for (size_t atom = 0; atom < atomLimit; ++atom)
    {
        if (static_cast<int>(atom) == excludedAtomIndex ||
            atomCounts[atom] <= 0)
            continue;
        const auto &probabilities =
            ProNovoConfig::getNaturalAtomIsotopeProbabilities(atom);
        const auto &masses = atomDistributions[atom].vMass;
        const size_t isotopeLimit =
            std::min(probabilities.size(), masses.size());
        for (size_t isotope = 1;
             isotope < isotopeLimit;
             ++isotope)
        {
            const double massDelta = masses[isotope] - masses[0];
            const int nominalShift =
                static_cast<int>(std::lround(massDelta));
            if (nominalShift <= 0 || nominalShift > maximum)
                continue;
            const double eventRate =
                static_cast<double>(atomCounts[atom]) *
                probabilities[isotope];
            lambda[static_cast<size_t>(nominalShift)] += eventRate;
            massLambda[static_cast<size_t>(nominalShift)] +=
                eventRate * massDelta;
        }
    }

    const double totalLambda =
        std::accumulate(lambda.begin() + 1, lambda.end(), 0.0);
    EnvelopeComponent background;
    background.probability.assign(
        static_cast<size_t>(maximum + 1), 0.0);
    background.massMoment.assign(
        background.probability.size(), 0.0);
    background.probability[0] = std::exp(-totalLambda);
    for (int shift = 1; shift <= maximum; ++shift)
    {
        double probability = 0.0;
        for (int isotopeShift = 1;
             isotopeShift <= shift;
             ++isotopeShift)
        {
            probability += static_cast<double>(isotopeShift) *
                           lambda[static_cast<size_t>(isotopeShift)] *
                           background.probability[
                               static_cast<size_t>(
                                   shift - isotopeShift)];
        }
        background.probability[static_cast<size_t>(shift)] =
            probability / static_cast<double>(shift);
    }
    for (int shift = 0; shift <= maximum; ++shift)
    {
        double moment = 0.0;
        for (int isotopeShift = 1;
             isotopeShift <= shift;
             ++isotopeShift)
        {
            moment += massLambda[
                          static_cast<size_t>(isotopeShift)] *
                      background.probability[
                          static_cast<size_t>(
                              shift - isotopeShift)];
        }
        background.massMoment[static_cast<size_t>(shift)] = moment;
    }
    normalizeAndTrim(background);
    return background;
}

static NominalEnvelopeModel theoreticalSipEnvelope(
    const std::array<int, 6> &atomCounts,
    const SupportedSipIsotope &spec,
    double targetFraction)
{
    constexpr int naturalTail = 32;
    NominalEnvelopeModel model;
    if (spec.atomIndex < 0 ||
        spec.atomIndex >= static_cast<int>(atomCounts.size()))
        return model;

    const int targetAtomCount =
        atomCounts[static_cast<size_t>(spec.atomIndex)];
    EnvelopeComponent target = targetElementEnvelope(
        targetAtomCount, spec, targetFraction);
    EnvelopeComponent background = naturalBackgroundEnvelope(
        atomCounts, spec.atomIndex, naturalTail);
    if (target.probability.empty() ||
        background.probability.empty())
        return model;

    const size_t modelSize =
        static_cast<size_t>(target.nominalOffset) +
        target.probability.size() +
        background.probability.size() - 1;
    model.probability.assign(modelSize, 0.0);
    std::vector<double> massMoment(modelSize, 0.0);
    for (size_t targetIndex = 0;
         targetIndex < target.probability.size();
         ++targetIndex)
    {
        for (size_t backgroundIndex = 0;
             backgroundIndex < background.probability.size();
             ++backgroundIndex)
        {
            const size_t modelIndex =
                static_cast<size_t>(target.nominalOffset) +
                targetIndex + backgroundIndex;
            model.probability[modelIndex] +=
                target.probability[targetIndex] *
                background.probability[backgroundIndex];
            massMoment[modelIndex] +=
                target.massMoment[targetIndex] *
                    background.probability[backgroundIndex] +
                target.probability[targetIndex] *
                    background.massMoment[backgroundIndex];
        }
    }

    const double sum = std::accumulate(
        model.probability.begin(), model.probability.end(), 0.0);
    if (!(sum > 0.0))
        return {};
    model.massDelta.assign(modelSize,
                           std::numeric_limits<double>::quiet_NaN());
    for (size_t index = 0; index < modelSize; ++index)
    {
        model.probability[index] /= sum;
        massMoment[index] /= sum;
        if (model.probability[index] > 0.0)
            model.massDelta[index] =
                massMoment[index] / model.probability[index];
    }
    return model;
}

static double weightedMedian(std::vector<std::pair<double, double>> values)
{
    if (values.empty())
        return 0.0;
    std::sort(values.begin(), values.end(),
              [](const auto &left, const auto &right)
              { return left.first < right.first; });
    double totalWeight = 0.0;
    for (const auto &value : values)
        totalWeight += value.second;
    const double middle = totalWeight * 0.5;
    double cumulative = 0.0;
    for (const auto &value : values)
    {
        cumulative += value.second;
        if (cumulative >= middle)
            return value.first;
    }
    return values.back().first;
}

PSMfeatureExtractor::Ms1AbundanceResult
PSMfeatureExtractor::getSIPelementAbundanceFromMS1Peaks(
    const std::vector<isotopicPeak> &peaks,
    double baseMass,
    const std::string &peptide,
    int precursorCharge,
    const std::string &sipAtom,
    double expectedEnrichmentPct)
{
    SupportedSipIsotope spec;
    if (peaks.empty() || precursorCharge <= 0 ||
        !resolveSupportedSipIsotope(sipAtom, spec) ||
        !std::isfinite(expectedEnrichmentPct))
        return {};

    averagine avg;
    avg.calPepAtomCounts(peptideBodyWithPtms(peptide));
    if (spec.atomIndex < 0 ||
        spec.atomIndex >= static_cast<int>(avg.pepAtomCounts.size()))
        return {};
    const double atomNumber =
        avg.pepAtomCounts[static_cast<size_t>(spec.atomIndex)];
    if (atomNumber <= 0.0)
        return {};

    const double baseMz =
        baseMass / precursorCharge + ProNovoConfig::getProtonMass();
    const double mzThreshold = baseMz - 0.5 / precursorCharge;
    std::vector<const isotopicPeak *> observed;
    observed.reserve(peaks.size());
    double rawIntensity = 0.0;
    double rawWeightedShift = 0.0;
    for (const isotopicPeak &peak : peaks)
    {
        if (peak.intensity <= 0.0 ||
            peak.mz <= mzThreshold ||
            peak.isotopeIndex < 0)
            continue;
        observed.push_back(&peak);
        rawIntensity += peak.intensity;
        rawWeightedShift +=
            peak.intensity * static_cast<double>(peak.isotopeIndex);
    }
    if (observed.empty() || rawIntensity <= 0.0)
        return {};

    const double initialPct =
        std::max(0.0, std::min(100.0, expectedEnrichmentPct));
    const double initialFraction = initialPct / 100.0;
    const auto shiftToPct = [&](double meanNominalShift,
                                double targetFraction)
    {
        const double naturalOtherShift =
            expectedNaturalNominalShiftExceptTarget(
                avg.pepAtomCounts,
                spec.atomIndex,
                spec.isotopeIndex,
                targetFraction);
        if (!std::isfinite(naturalOtherShift))
            return std::numeric_limits<double>::quiet_NaN();
        return (meanNominalShift - naturalOtherShift) /
               (atomNumber * spec.nominalShift) * 100.0;
    };

    double rawPct = shiftToPct(
        rawWeightedShift / rawIntensity, initialFraction);
    if (!std::isfinite(rawPct))
        return {};
    rawPct = std::max(0.0, std::min(100.0, rawPct));

    double fittedPct = initialPct;
    const double binomialVariance =
        atomNumber * initialFraction * (1.0 - initialFraction);
    const bool broadEnvelope = binomialVariance >= 2.0;
    const int reweightIterations = broadEnvelope ? 2 : 1;
    const double interferenceCapFactor = broadEnvelope ? 4.0 : 2.0;
    int compatiblePeakCount = 0;
    double modelCoverage = 0.0;
    bool fitted = false;

    // The selected isotope supplies the target-element distribution; all
    // remaining isotopes, including O17 and S33/S36, remain in the model.
    // Broad envelopes support two robust moment updates. Sparse envelopes get
    // one conservative update so positive noise cannot cause fixed-point drift.
    for (int iteration = 0;
         iteration < reweightIterations;
         ++iteration)
    {
        const double modelFraction = fittedPct / 100.0;
        const NominalEnvelopeModel model =
            theoreticalSipEnvelope(
                avg.pepAtomCounts, spec, modelFraction);
        if (model.probability.empty())
            break;
        const double modelMaximum =
            *std::max_element(
                model.probability.begin(), model.probability.end());
        const double probabilityFloor = modelMaximum * 1e-8;

        std::vector<std::pair<double, double>> ratios;
        ratios.reserve(observed.size());
        for (const isotopicPeak *peak : observed)
        {
            const size_t index =
                static_cast<size_t>(peak->isotopeIndex);
            if (index >= model.probability.size() ||
                model.probability[index] < probabilityFloor)
                continue;
            ratios.emplace_back(
                peak->intensity / model.probability[index],
                model.probability[index]);
        }
        const double scale = weightedMedian(std::move(ratios));
        if (!(scale > 0.0) || !std::isfinite(scale))
            break;

        double cappedIntensity = 0.0;
        double cappedWeightedShift = 0.0;
        double observedModelProbability = 0.0;
        double observedModelWeightedShift = 0.0;
        compatiblePeakCount = 0;
        for (const isotopicPeak *peak : observed)
        {
            const size_t index =
                static_cast<size_t>(peak->isotopeIndex);
            if (index >= model.probability.size() ||
                model.probability[index] < probabilityFloor)
                continue;

            const double capped = std::min(
                peak->intensity,
                interferenceCapFactor * scale *
                    model.probability[index]);
            cappedIntensity += capped;
            cappedWeightedShift +=
                capped * static_cast<double>(peak->isotopeIndex);
            observedModelProbability += model.probability[index];
            observedModelWeightedShift +=
                model.probability[index] *
                static_cast<double>(peak->isotopeIndex);
            ++compatiblePeakCount;
        }
        if (!(cappedIntensity > 0.0) ||
            !(observedModelProbability > 0.0))
            break;

        double fullModelMean = 0.0;
        for (size_t index = 0;
             index < model.probability.size();
             ++index)
        {
            fullModelMean +=
                model.probability[index] *
                static_cast<double>(index);
        }
        const double observedMean =
            cappedWeightedShift / cappedIntensity;
        const double conditionalModelMean =
            observedModelWeightedShift / observedModelProbability;
        const double correctedMean =
            observedMean + fullModelMean - conditionalModelMean;
        double update =
            shiftToPct(correctedMean, modelFraction);
        if (!std::isfinite(update))
            break;

        update = std::max(0.0, std::min(100.0, update));
        modelCoverage = observedModelProbability;
        fittedPct = update;
        fitted = true;
    }

    if (!fitted)
        return {rawPct, static_cast<int>(observed.size()), false};

    const bool valid =
        compatiblePeakCount >= 2 && modelCoverage >= 0.02;
    if (!valid)
        return {rawPct, static_cast<int>(observed.size()), false};
    return {fittedPct, compatiblePeakCount, true};
}

static std::filesystem::path resolveHdf5FeaturePath(const std::string &hdf5BasePath)
{
    std::filesystem::path base(hdf5BasePath);
    if (std::filesystem::exists(base) && std::filesystem::is_regular_file(base))
    {
        return base;
    }
    std::filesystem::path h5 = base;
    h5 += ".h5";
    if (std::filesystem::exists(h5))
    {
        return h5;
    }
    std::filesystem::path hdf5 = base;
    hdf5 += ".hdf5";
    if (std::filesystem::exists(hdf5))
    {
        return hdf5;
    }
    throw std::runtime_error("Cannot find Raxport HDF5 file for feature extraction: " + hdf5BasePath);
}

void PSMfeatureExtractor::loadHdf5Ms1(const std::string &hdf5BasePath)
{
    const std::filesystem::path hdf5Path = resolveHdf5FeaturePath(hdf5BasePath);
    ms1Data.clear();
    std::string error;
    if (!sipros::readRaxportHdf5Ms1Data(hdf5Path.string(), ms1Data, error))
    {
        throw std::runtime_error(error);
    }
}

void PSMfeatureExtractor::initializeFeatureVectors(sipPSM &psm)
{
    const size_t count = psm.scanNumbers.size();
    psm.isotopicPeakss = std::vector<std::vector<isotopicPeak>>(count);
    psm.isotopicPeakNumbers = std::vector<int>(count);
    psm.MS1IsotopicAbundances = std::vector<double>(count);
    if (psm.MS2IsotopicAbundances.size() != count)
        psm.MS2IsotopicAbundances = std::vector<double>(count, 1.07);
    psm.isotopicAbundanceDiffs = std::vector<double>(count);
    psm.peptideLengths = std::vector<int>(count);
    psm.missCleavageSiteNumbers = std::vector<int>(count);
    psm.PTMnumbers = std::vector<int>(count);
    psm.diffScores = std::vector<float>(count);
    psm.mzShiftFromisolationWindowCenters = std::vector<double>(count);
    psm.isotopicMassWindowShifts = std::vector<int>(count);
    psm.massErrors = std::vector<double>(count);
    psm.precursorIntensities = std::vector<double>(count, 0);
}

void PSMfeatureExtractor::extractFeaturesForPsm(const std::string &hdf5Path, sipPSM &psm)
{
    loadHdf5Ms1(hdf5Path);
    mSipPSM = &psm;
    initializeFeatureVectors(psm);
    extractFeaturesOfEachPSM();
}

void PSMfeatureExtractor::extractFeaturesOfEachPSM()
{
    float topScore = 0;
    for (size_t i = 0; i < mSipPSM->isotopicPeakss.size(); i++)
    {
        const std::string compositionPeptide = peptideSequenceForComposition(mSipPSM->identifiedPeptides[i]);
        const int precursorCharge = mSipPSM->parentCharges[i];
        if (precursorCharge > 0)
        {
            const double baseMass = mAveragine.calPrecursorBaseMass(compositionPeptide);
            const double monoPrecursorMz = baseMass / precursorCharge + ProNovoConfig::getProtonMass();
            const double matchedPrecursorMz =
                mSipPSM->measuredParentMasses[i] / precursorCharge + ProNovoConfig::getProtonMass();
            const std::string sipIsotope =
                canonicalSipIsotope(ProNovoConfig::getSetSIPelement());
            const auto mzToleranceDaAt = [precursorCharge](double)
            { return 0.01 / precursorCharge; };
            mSipPSM->isotopicPeakss[i] = findMs1IsotopicPeaks(
                &ms1Data,
                mSipPSM->precursorScanNumbers[i],
                precursorCharge,
                monoPrecursorMz,
                matchedPrecursorMz,
                mAveragine.pepAtomCounts,
                sipIsotope,
                mSipPSM->MS2IsotopicAbundances[i],
                mzToleranceDaAt);
            const Ms1AbundanceResult ms1Abundance = getSIPelementAbundanceFromMS1Peaks(
                mSipPSM->isotopicPeakss[i],
                baseMass,
                compositionPeptide,
                precursorCharge,
                sipIsotope,
                mSipPSM->MS2IsotopicAbundances[i]);
            mSipPSM->isotopicPeakNumbers[i] =
                ms1Abundance.valid ? ms1Abundance.isotopicPeakCount : 0;
            mSipPSM->MS1IsotopicAbundances[i] = ms1Abundance.abundancePct;
        }
        std::tie(mSipPSM->peptideLengths[i], mSipPSM->missCleavageSiteNumbers[i]) =
            getSeqLengthAndMissCleavageSiteNumber(mSipPSM->originalPeptides[i]);
        mSipPSM->PTMnumbers[i] = getPTMnumber(mSipPSM->identifiedPeptides[i]);

        if (mSipPSM->ranks[i] == 1)
        {
            // topScore = mSipPSM->MVHscores[i];
            topScore = mSipPSM->WDPscores[i];
        }
        // mSipPSM->MVHdiffScores[i] = topMVHscore - mSipPSM->MVHscores[i];
        mSipPSM->diffScores[i] = topScore - mSipPSM->WDPscores[i];

        mSipPSM->mzShiftFromisolationWindowCenters[i] = std::abs(
            mSipPSM->isolationWindowCenterMZs[i] -
            mSipPSM->measuredParentMasses[i] / mSipPSM->parentCharges[i] - ProNovoConfig::getProtonMass());
        std::tie(mSipPSM->isotopicMassWindowShifts[i], mSipPSM->massErrors[i]) = getMassWindowShiftAndError(
            mSipPSM->measuredParentMasses[i], mSipPSM->calculatedParentMasses[i]);

        mSipPSM->precursorIntensities[i] = 0;
        for (auto &peak : mSipPSM->isotopicPeakss[i])
        {
            mSipPSM->precursorIntensities[i] += peak.intensity;
        }

        mSipPSM->isotopicAbundanceDiffs[i] = mSipPSM->MS1IsotopicAbundances[i] - mSipPSM->MS2IsotopicAbundances[i];
    }
}
