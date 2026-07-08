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

static char sipAtomLetter(const std::string &sipAtom)
{
    for (unsigned char c : sipAtom)
    {
        if (std::isalpha(c))
            return static_cast<char>(std::toupper(c));
    }
    return 'C';
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
    int count = 0;
    for (char &A : cleavageSites)
    {
        for (char &S : seq)
        {
            if (S == A)
                count++;
        }
    }
    return {static_cast<int>(seq.size()), count};
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
    for (size_t i = 0; i + 1 < naked.size(); ++i)
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

int PSMfeatureExtractor::sipAtomIndex(const std::string &sipAtom)
{
    const char atom = sipAtomLetter(sipAtom);
    if (atom == 'H')
        return 1;
    if (atom == 'O')
        return 2;
    if (atom == 'N')
        return 3;
    if (atom == 'S')
        return 5;
    return 0;
}

int PSMfeatureExtractor::sipNominalShiftPerAtom(const std::string &sipAtom)
{
    const char atom = sipAtomLetter(sipAtom);
    return (atom == 'O' || atom == 'S') ? 2 : 1;
}

double PSMfeatureExtractor::expectedNaturalNominalShiftExceptTarget(const std::array<int, 6> &atomCounts,
                                                                    int targetAtomIndex,
                                                                    int targetIsotopeIndex)
{
    double expectedShift = 0.0;
    const auto &atomDistributions = ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution;
    const size_t atomCount = std::min(atomCounts.size(), atomDistributions.size());
    for (size_t atomIdx = 0; atomIdx < atomCount; ++atomIdx)
    {
        if (atomCounts[atomIdx] <= 0)
            continue;
        const auto &probs = atomDistributions[atomIdx].vProb;
        for (size_t isotopeIdx = 1; isotopeIdx < probs.size(); ++isotopeIdx)
        {
            if (static_cast<int>(atomIdx) == targetAtomIndex &&
                static_cast<int>(isotopeIdx) == targetIsotopeIndex)
            {
                continue;
            }
            expectedShift += static_cast<double>(atomCounts[atomIdx]) *
                             static_cast<double>(isotopeIdx) * probs[isotopeIdx];
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

std::vector<isotopicPeak> PSMfeatureExtractor::findMs1IsotopicPeaks(
    const sipros::RaxportMs1Data *ms1Data,
    int &ms1ScanNumber,
    int precursorCharge,
    double monoPrecursorMz,
    double matchedPrecursorMz,
    int targetNominalShift,
    const std::function<double(double)> &mzToleranceDaAt)
{
    std::vector<isotopicPeak> peaks;
    if (!ms1Data || precursorCharge <= 0)
        return peaks;
    auto scanIt = ms1Data->scanNumberToIndex.find(ms1ScanNumber);
    if (scanIt == ms1Data->scanNumberToIndex.end())
        return peaks;

    const double neutronMz = ProNovoConfig::getNeutronMass() / precursorCharge;
    const int isotopeWindow = Ms1IsotopeWindow * std::max(1, targetNominalShift);
    constexpr int kMs1AssignmentIndexTolerance = 1;
    const int assignedIndex = static_cast<int>(std::round((matchedPrecursorMz - monoPrecursorMz) / neutronMz));
    if (assignedIndex < 0)
        return peaks;
    const int firstIsotopeIndex = std::max(0, assignedIndex - isotopeWindow);
    const int lastIsotopeIndex = assignedIndex + isotopeWindow;

    size_t scanIdx = scanIt->second;
    for (int attempt = 0; attempt < 3; ++attempt)
    {
        const sipros::RaxportMs1Scan &scan = ms1Data->scans[scanIdx];
        peaks.clear();
        for (int isotopeIndex = firstIsotopeIndex; isotopeIndex <= lastIsotopeIndex; ++isotopeIndex)
        {
            const double expectedMz = monoPrecursorMz + isotopeIndex * neutronMz;
            size_t idx = findMs1Peak(scan, expectedMz, mzToleranceDaAt, precursorCharge);
            if (idx == std::numeric_limits<size_t>::max())
                idx = findMs1Peak(scan, expectedMz, mzToleranceDaAt);
            if (idx != std::numeric_limits<size_t>::max())
            {
                peaks.push_back({scan.mz[idx],
                                 ms1PeakCharge(scan, idx),
                                 scan.intensity[idx],
                                 isotopeIndex});
            }
        }
        bool hasAssignmentAnchor = false;
        bool hasExactAssignmentAnchor = false;
        for (const isotopicPeak &peak : peaks)
        {
            const int indexDelta = std::abs(peak.isotopeIndex - assignedIndex);
            if (indexDelta <= kMs1AssignmentIndexTolerance)
                hasAssignmentAnchor = true;
            if (indexDelta == 0)
                hasExactAssignmentAnchor = true;
        }
        if (!peaks.empty() && hasAssignmentAnchor &&
            (peaks.size() > 1 || hasExactAssignmentAnchor))
        {
            ms1ScanNumber = scan.scanNumber;
            break;
        }
        peaks.clear();
        if (scanIdx == 0)
            break;
        --scanIdx;
    }
    if (peaks.empty())
        return peaks;

    std::sort(peaks.begin(), peaks.end(), [](const isotopicPeak &a, const isotopicPeak &b)
              { return a.mz < b.mz; });
    peaks.erase(std::unique(peaks.begin(), peaks.end(),
                            [](const isotopicPeak &a, const isotopicPeak &b)
                            { return std::fabs(a.mz - b.mz) <= 1e-12; }),
                peaks.end());
    return peaks;
}

PSMfeatureExtractor::Ms1AbundanceResult PSMfeatureExtractor::getSIPelementAbundanceFromMS1Peaks(
    const std::vector<isotopicPeak> &peaks,
    double baseMass,
    const std::string &peptide,
    int precursorCharge,
    const std::string &sipAtom)
{
    if (peaks.empty() || precursorCharge <= 0)
        return {};

    averagine avg;
    const char atom = sipAtomLetter(sipAtom);
    const int targetIsotopeIndex = sipNominalShiftPerAtom(std::string(1, atom));
    avg.calPepAtomCounts(peptideBodyWithPtms(peptide));
    const int atomIndex = sipAtomIndex(sipAtom);
    const double atomNumber = avg.pepAtomCounts[atomIndex];
    if (atomNumber <= 0.0)
        return {};

    const double maxIsotopeIndex = atomNumber * targetIsotopeIndex;
    const double baseMz = baseMass / precursorCharge + ProNovoConfig::getProtonMass();

    double sumIntensity = 0.0;
    double weightedIsotopeIndex = 0.0;
    int validPeakCount = 0;
    for (size_t i = 0; i < peaks.size(); ++i)
    {
        const isotopicPeak &peak = peaks[i];
        if (peak.intensity <= 0.0)
            continue;
        const int isotopeIndex = peak.isotopeIndex >= 0
                                     ? peak.isotopeIndex
                                     : static_cast<int>(std::round((peak.mz - baseMz) /
                                                                  ProNovoConfig::getNeutronMass() * precursorCharge));
        if (isotopeIndex < 0)
            continue;
        sumIntensity += peak.intensity;
        weightedIsotopeIndex += peak.intensity * std::min(static_cast<double>(isotopeIndex), maxIsotopeIndex);
        ++validPeakCount;
    }
    if (validPeakCount == 0 || sumIntensity <= 0.0)
        return {};

    const double meanIsotopeIndex = weightedIsotopeIndex / sumIntensity;
    const double naturalOtherShift = expectedNaturalNominalShiftExceptTarget(
        avg.pepAtomCounts, atomIndex, targetIsotopeIndex);
    double pct = (meanIsotopeIndex - naturalOtherShift) / maxIsotopeIndex * 100.0;
    if (!std::isfinite(pct))
        return {};
    pct = std::min(100.0, std::max(0.0, pct));
    return {pct, validPeakCount};
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
            const int targetNominalShift = sipNominalShiftPerAtom(ProNovoConfig::getSetSIPelement());
            const auto mzToleranceDaAt = [precursorCharge](double)
            { return 0.01 / precursorCharge; };
            mSipPSM->isotopicPeakss[i] = findMs1IsotopicPeaks(&ms1Data,
                                                              mSipPSM->precursorScanNumbers[i],
                                                              precursorCharge,
                                                              monoPrecursorMz,
                                                              matchedPrecursorMz,
                                                              targetNominalShift,
                                                              mzToleranceDaAt);
            const Ms1AbundanceResult ms1Abundance = getSIPelementAbundanceFromMS1Peaks(
                mSipPSM->isotopicPeakss[i],
                baseMass,
                compositionPeptide,
                precursorCharge,
                ProNovoConfig::getSetSIPelement());
            mSipPSM->isotopicPeakNumbers[i] = ms1Abundance.isotopicPeakCount;
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
