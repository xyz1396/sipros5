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
    const int assignedIndex = std::max(0, static_cast<int>(std::round((matchedPrecursorMz - monoPrecursorMz) / neutronMz)));

    size_t scanIdx = scanIt->second;
    size_t anchorIdx = std::numeric_limits<size_t>::max();
    const sipros::RaxportMs1Scan *anchorScan = nullptr;
    for (int attempt = 0; attempt < 3; ++attempt)
    {
        const sipros::RaxportMs1Scan &scan = ms1Data->scans[scanIdx];
        anchorIdx = findMs1Peak(scan, matchedPrecursorMz, mzToleranceDaAt, precursorCharge);
        if (anchorIdx == std::numeric_limits<size_t>::max())
            anchorIdx = findMs1Peak(scan, matchedPrecursorMz, mzToleranceDaAt);
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
    peaks.push_back({scan.mz[anchorIdx],
                     ms1PeakCharge(scan, anchorIdx),
                     scan.intensity[anchorIdx],
                     assignedIndex});

    constexpr int kLegacyIsotopePeaksEachSide = 20;
    for (int direction : {-1, 1})
    {
        if ((direction < 0 && anchorIdx == 0) ||
            (direction > 0 && anchorIdx + 1 >= scan.mz.size()))
        {
            continue;
        }
        int currentIndex = static_cast<int>(anchorIdx) + direction;
        for (int iso = 1; iso <= kLegacyIsotopePeaksEachSide; ++iso)
        {
            bool foundPeak = false;
            size_t foundIdx = 0;
            double maxIntensity = 0.0;
            const double expectedMz = scan.mz[anchorIdx] + direction * iso * neutronMz;
            const double mzTolerance = mzToleranceDaAt(expectedMz);
            while (currentIndex >= 0 && currentIndex < static_cast<int>(scan.mz.size()))
            {
                const size_t idx = static_cast<size_t>(currentIndex);
                if (direction * (scan.mz[idx] - expectedMz) >= mzTolerance)
                    break;
                if (std::fabs(expectedMz - scan.mz[idx]) < mzTolerance &&
                    scan.intensity[idx] > maxIntensity)
                {
                    foundPeak = true;
                    foundIdx = idx;
                    maxIntensity = scan.intensity[idx];
                }
                currentIndex += direction;
            }
            if (!foundPeak)
                break;
            peaks.push_back({scan.mz[foundIdx],
                             ms1PeakCharge(scan, foundIdx),
                             scan.intensity[foundIdx],
                             assignedIndex + direction * iso});
        }
    }

    std::sort(peaks.begin(), peaks.end(), [](const isotopicPeak &a, const isotopicPeak &b)
              { return a.mz < b.mz; });
    if (peaks.size() > 2)
    {
        std::vector<int> vertexIndices;
        double lastIntensityDiff = 1.0;
        for (size_t i = 1; i < peaks.size(); ++i)
        {
            const double intensityDiff = peaks[i].intensity - peaks[i - 1].intensity;
            if (lastIntensityDiff > 0.0 && intensityDiff < 0.0)
                vertexIndices.push_back(static_cast<int>(i - 1));
            lastIntensityDiff = intensityDiff;
        }
        if (lastIntensityDiff > 0.0)
            vertexIndices.push_back(static_cast<int>(peaks.size() - 1));

        int closestVertexIndex = 0;
        double closestMzDiff = std::numeric_limits<double>::max();
        for (int vertexIndex : vertexIndices)
        {
            const double mzDiff = std::fabs(peaks[vertexIndex].mz - monoPrecursorMz);
            if (mzDiff < closestMzDiff)
            {
                closestMzDiff = mzDiff;
                closestVertexIndex = vertexIndex;
            }
        }

        int currentPeakIndex = closestVertexIndex - 1;
        while (currentPeakIndex >= 0)
        {
            const double intensityDiff = peaks[currentPeakIndex].intensity - peaks[currentPeakIndex + 1].intensity;
            if (intensityDiff > 0.0)
            {
                peaks.erase(peaks.begin(), peaks.begin() + currentPeakIndex + 1);
                closestVertexIndex -= currentPeakIndex + 1;
                break;
            }
            --currentPeakIndex;
        }

        currentPeakIndex = closestVertexIndex + 1;
        while (currentPeakIndex < static_cast<int>(peaks.size()))
        {
            const double intensityDiff = peaks[currentPeakIndex].intensity - peaks[currentPeakIndex - 1].intensity;
            if (intensityDiff > 0.0)
            {
                peaks.erase(peaks.begin() + currentPeakIndex, peaks.end());
                break;
            }
            ++currentPeakIndex;
        }
    }
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
    avg.calPepAtomCounts(peptideBodyWithPtms(peptide));
    const int atomIndex = sipAtomIndex(sipAtom);
    const double atomNumber = avg.pepAtomCounts[atomIndex];
    if (atomNumber <= 0.0)
        return {};

    const double baseMz = baseMass / precursorCharge + ProNovoConfig::getProtonMass();
    const double mzThreshold = baseMz - 0.5 / precursorCharge;
    std::vector<double> usefulPeakIntensities;
    usefulPeakIntensities.reserve(peaks.size());
    double firstUsefulPeakMz = 0.0;
    bool foundUsefulPeak = false;

    for (const isotopicPeak &peak : peaks)
    {
        if (peak.intensity <= 0.0 || peak.mz <= mzThreshold)
            continue;
        if (!foundUsefulPeak)
        {
            firstUsefulPeakMz = peak.mz;
            foundUsefulPeak = true;
        }
        usefulPeakIntensities.push_back(peak.intensity);
    }
    if (usefulPeakIntensities.empty())
        return {};

    const int firstDeltaNeutron = static_cast<int>(std::round(
        (firstUsefulPeakMz - baseMz) / ProNovoConfig::getNeutronMass() * precursorCharge));
    const double sumIntensity = std::accumulate(usefulPeakIntensities.begin(),
                                                usefulPeakIntensities.end(), 0.0);
    if (sumIntensity <= 0.0)
        return {};

    double abundance = 0.0;
    for (size_t i = 0; i < usefulPeakIntensities.size(); ++i)
    {
        abundance += usefulPeakIntensities[i] / sumIntensity *
                     static_cast<double>(static_cast<int>(i) + firstDeltaNeutron);
    }

    double nominalShiftPerAtom = 1.0;
    const char atom = sipAtomLetter(sipAtom);
    if (atom == 'O' || atom == 'S')
        nominalShiftPerAtom = 2.0;

    const double pct = abundance / (atomNumber * nominalShiftPerAtom) * 100.0;
    if (!std::isfinite(pct))
        return {};
    return {pct, static_cast<int>(usefulPeakIntensities.size())};
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
            mSipPSM->isotopicPeakNumbers[i] = static_cast<int>(mSipPSM->isotopicPeakss[i].size());
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
