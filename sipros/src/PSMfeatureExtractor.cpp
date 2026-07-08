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

size_t PSMfeatureExtractor::binarySearchPeak(const Scan *mScan, double Mz, int charge)
{
    // init peakIX to invalid value
    size_t peakIX = std::numeric_limits<size_t>::max();
    if (mScan == nullptr || mScan->mz.empty() || charge <= 0)
        return peakIX;
    size_t low = 0;
    size_t high = mScan->mz.size() - 1;
    size_t mid = 0;
    double diff = 0;
    const double mzTolerance = 0.01 / charge;
    double currentIntensity = 0;
    while (low <= high)
    {
        mid = (low + high) / 2;
        diff = std::abs(mScan->mz[mid] - Mz);
        if (diff <= mzTolerance)
        {
            // find the peak with highest intensity in tolerance range
            // if (mScan->charge[mid] == charge && mScan->intensity[mid] > currentIntensity)
            if (mScan->intensity[mid] > currentIntensity)
            {
                peakIX = mid;
                currentIntensity = mScan->intensity[mid];
            }
        }
        if (mScan->mz[mid] < Mz) // Search the right half
        {
            low = mid + 1;
        }
        else // Search the left half
        {
            // in case searched the first peak in mScan
            if (mid == 0)
                break;
            high = mid - 1;
        }
    }
    return peakIX;
}

void PSMfeatureExtractor::filterIsotopicPeaks(std::vector<isotopicPeak> &isotopicPeaks,
                                              const double calculatedPrecursorMZ)
{
    double intensityDiff = 0, lastIntensityDiff = 1.0;
    vertexIXs.clear();
    for (size_t i = 1; i < isotopicPeaks.size(); i++)
    {
        intensityDiff = isotopicPeaks[i].intensity - isotopicPeaks[i - 1].intensity;
        if (lastIntensityDiff > 0 && intensityDiff < 0)
            vertexIXs.push_back(i - 1);
        lastIntensityDiff = intensityDiff;
    }
    if (lastIntensityDiff > 0)
        vertexIXs.push_back(isotopicPeaks.size() - 1);
    // find closest vertexIX to calculatedPrecursorMZ
    int closestVertexIX = 0;
    double closestVertexMZdiff = std::numeric_limits<double>::max(), MZdiff = 0;
    for (size_t i = 0; i < vertexIXs.size(); i++)
    {
        MZdiff = std::abs(isotopicPeaks[vertexIXs[i]].mz - calculatedPrecursorMZ);
        if (MZdiff < closestVertexMZdiff)
        {
            closestVertexMZdiff = MZdiff;
            closestVertexIX = vertexIXs[i];
        }
    }
    // go left slope of isotopic envelope
    int currentPeakIX = closestVertexIX - 1;
    while (currentPeakIX >= 0)
    {
        intensityDiff = isotopicPeaks[currentPeakIX].intensity - isotopicPeaks[currentPeakIX + 1].intensity;
        if (intensityDiff > 0)
        {
            isotopicPeaks.erase(isotopicPeaks.begin(), isotopicPeaks.begin() + currentPeakIX + 1);
            // shift closestVertexIX because of erasing
            closestVertexIX = closestVertexIX - currentPeakIX - 1;
            break;
        }
        currentPeakIX--;
    }
    // go right slope of isotopic envelope
    currentPeakIX = closestVertexIX + 1;
    while (currentPeakIX < (int)isotopicPeaks.size())
    {
        intensityDiff = isotopicPeaks[currentPeakIX].intensity - isotopicPeaks[currentPeakIX - 1].intensity;
        if (intensityDiff > 0)
        {
            isotopicPeaks.erase(isotopicPeaks.begin() + currentPeakIX, isotopicPeaks.end());
            break;
        }
        currentPeakIX++;
    }
}

void PSMfeatureExtractor::filterIsotopicPeaksTopN(
    std::vector<isotopicPeak> &isotopicPeaks, const double observedPrecursorMZ,
    const size_t topN) {
    // if (isotopicPeaks.empty()) return;
    // Find the highest intensity peak within calculatedPrecursorMZ ± mzWindow
    // the best DDA isolation window is 5 m/z, so we use a ±5 m/z window.
    double mzWindow = 5;
    double maxIntensity = 0.0;
    size_t maxIdx = 0;
    for (size_t i = 0; i < isotopicPeaks.size(); ++i) {
        if (std::abs(isotopicPeaks[i].mz - observedPrecursorMZ) <= mzWindow) {
            if (isotopicPeaks[i].intensity > maxIntensity) {
                maxIntensity = isotopicPeaks[i].intensity;
                maxIdx = i;
            }
        }
    }
    // Collect peaks within highest peak mz ± mzWindow
    double centerMz = isotopicPeaks[maxIdx].mz;
    std::vector<isotopicPeak> peaksInWindow;
    peaksInWindow.reserve(isotopicPeaks.size());
    for (const auto &peak : isotopicPeaks) {
        if (std::abs(peak.mz - centerMz) <= mzWindow) {
            peaksInWindow.push_back(peak);
        }
    }
    if (peaksInWindow.size() > topN) {
        // Sort by intensity descending and keep topN
        std::sort(peaksInWindow.begin(), peaksInWindow.end(),
                  [](const isotopicPeak &a, const isotopicPeak &b) {
                      return a.intensity > b.intensity;
                  });
        peaksInWindow.resize(topN);
        std::sort(peaksInWindow.begin(), peaksInWindow.end(),
                  [](const isotopicPeak &a, const isotopicPeak &b) {
                      return a.mz < b.mz;
                  });
    }
    isotopicPeaks = std::move(peaksInWindow);
}

std::vector<isotopicPeak> PSMfeatureExtractor::
    findIsotopicPeaks(int &MS1ScanNumber,
                      const int precursorCharge,
                      const double observedPrecursorMass,
                      const double calculatedPrecursorMass)
{
    std::vector<isotopicPeak> isotopicPeaks = {};
    double observedPrecursorMZ = observedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    double calculatedPrecursorMZ = calculatedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    const auto ms1It = scanNumberMS1ScanMap.find(MS1ScanNumber);
    if (ms1It == scanNumberMS1ScanMap.end() || ms1It->second == nullptr)
        return isotopicPeaks;
    Scan *MS1Scan = ms1It->second;
    size_t peakIX = std::numeric_limits<size_t>::max();
    // if first search failed, search 2 scans before it
    for (int i = 0; i < 3; i++)
    {
        peakIX = binarySearchPeak(MS1Scan, observedPrecursorMZ, precursorCharge);
        if (peakIX != std::numeric_limits<size_t>::max())
        {
            MS1ScanNumber = MS1Scan->scanNumber;
            break;
        }
        else if (MS1Scan != MS1Scans.data())
            MS1Scan--;
        else
            break;
    }
    size_t currentIX = 0;
    size_t foundIX = 0;
    double currentMass = 0;
    double expectedMass = 0;
    double maxIntensity = 0;
    double massTolerance = 0.01;
    bool foundIsotopicPeak = false;
    if (peakIX != std::numeric_limits<size_t>::max())
    {
        // try NisotopicPeak isotopic peaks each side
        int tryISO = NisotopicPeak;
        isotopicPeaks.reserve(tryISO);
        isotopicPeaks.push_back({MS1Scan->mz[peakIX],
                                 precursorCharge, MS1Scan->intensity[peakIX]});
        // Go left and right side
        for (int direction : {-1, 1})
        {
            if ((direction < 0 && peakIX == 0) ||
                (direction > 0 && peakIX + 1 >= MS1Scan->mz.size()))
            {
                continue;
            }
            int currentIndex = static_cast<int>(peakIX) + direction;
            for (int iso = 1; iso <= tryISO; iso++)
            {
                foundIsotopicPeak = false;
                maxIntensity = 0;
                expectedMass = MS1Scan->mz[peakIX] * precursorCharge +
                               direction * iso * ProNovoConfig::getNeutronMass();
                while (currentIndex >= 0 && currentIndex < static_cast<int>(MS1Scan->mz.size()))
                {
                    currentIX = static_cast<size_t>(currentIndex);
                    currentMass = MS1Scan->mz[currentIX] * precursorCharge;
                    if (direction * (currentMass - expectedMass) >= massTolerance)
                        break;
                    // find the matched isotopic peak with max intensity
                    if (std::abs(expectedMass - currentMass) < massTolerance &&
                        // MS1Scan->charge[currentIX] == precursorCharge &&
                        MS1Scan->intensity[currentIX] > maxIntensity)
                    {
                        foundIsotopicPeak = true;
                        foundIX = currentIX;
                        maxIntensity = MS1Scan->intensity[currentIX];
                    }
                    currentIndex += direction;
                }
                if (foundIsotopicPeak)
                    isotopicPeaks.push_back({MS1Scan->mz[foundIX],
                                             //  precursorCharge,
                                             foundIX < MS1Scan->charge.size() ? MS1Scan->charge[foundIX] : precursorCharge,
                                             MS1Scan->intensity[foundIX]});
                else
                    break;
            }
        }
        std::sort(isotopicPeaks.begin(), isotopicPeaks.end(), [](const isotopicPeak &a, const isotopicPeak &b)
                  { return a.mz < b.mz; });
    }
    if (isotopicPeaks.size() > 2)
    // if (isotopicPeaks.size() > NisotopicPeak / 4) 
    {
        filterIsotopicPeaks(isotopicPeaks, calculatedPrecursorMZ);
        // filterIsotopicPeaksTopN(isotopicPeaks, observedPrecursorMZ,
        //                         NisotopicPeak);
    }
    return isotopicPeaks;
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

double PSMfeatureExtractor::getSIPelementAbundanceFromMS1(const std::string &peptideSeq,
                                                          const std::vector<isotopicPeak> &isotopicPeaks,
                                                          const int precursorCharge)
{
    // in case of no isotopic peaks
    if (isotopicPeaks.size() == 0)
        return 0.0;
    double baseMass = mAveragine.calPrecursorBaseMass(peptideSeq);
    int charge = precursorCharge;
    double baseMZ = baseMass / charge + ProNovoConfig::getProtonMass();
    double MZthreshold = baseMZ - 0.5 / charge;
    std::vector<double> usefulIsotopicPeakIntensity;
    usefulIsotopicPeakIntensity.reserve(NisotopicPeak);
    for (const auto &peak : isotopicPeaks)
    {
        if (peak.mz > MZthreshold)
        {
            usefulIsotopicPeakIntensity.push_back(peak.intensity);
        }
    }
    if (usefulIsotopicPeakIntensity.empty())
        return 0.0;
    int firstUsefulIsotopicPeakIX = isotopicPeaks.size() - usefulIsotopicPeakIntensity.size();
    if (firstUsefulIsotopicPeakIX < 0 || firstUsefulIsotopicPeakIX >= static_cast<int>(isotopicPeaks.size()))
        return 0.0;
    double firstUsefulIsotopicPeakMZ = isotopicPeaks[firstUsefulIsotopicPeakIX].mz;
    int firstDeltaNeutron = static_cast<int>(std::round((firstUsefulIsotopicPeakMZ - baseMZ) /
                                                        ProNovoConfig::getNeutronMass() * charge));
    double sumOfIntensities = std::accumulate(usefulIsotopicPeakIntensity.begin(),
                                              usefulIsotopicPeakIntensity.end(), 0.0);
    if (sumOfIntensities <= 0.0)
        return 0.0;
    for (auto &intensity : usefulIsotopicPeakIntensity)
    {
        intensity /= sumOfIntensities;
    }
    double pct = 0.0;
    for (size_t i = 0; i < usefulIsotopicPeakIntensity.size(); i++)
    {
        pct += usefulIsotopicPeakIntensity[i] * (i + firstDeltaNeutron);
    }
    double atomNumber = mAveragine.pepAtomCounts[mAveragine.SIPatomIX];
    if (atomNumber <= 0.0)
        return 0.0;
    double deltaNeutron = 1.0;
    const std::string &sipElement = ProNovoConfig::getSetSIPelement();
    // // Decrease the atom number for C element to fit the real data
    // // for the isotopic peaks filtered by filterIsotopicPeaksTopN
    // if (sipElement == "C") {
    //     atomNumber --;
    // }
    if (sipElement == "O" || sipElement == "S") {
        deltaNeutron = 2.0;
    }
    // if (sipElement != "C") {
    //     pct -= mAveragine.pepAtomCounts[0] * 0.0107; // Adjust for non-carbon elements
    // }
    pct = pct / (atomNumber * deltaNeutron) * 100.;
    // In case Non-carbon elements SIP exceed 100
    // pct = std::clamp(pct, 0.0, 100.0); // Clamp to [0, 100]
    return pct;
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

double PSMfeatureExtractor::getMS2IsotopicAbundance(const std::string &searchName)
{
    if (searchName == "SE")
        return 1.07;
    std::string pct = searchName;
    std::string delimiter = "_";
    size_t pos = 0;
    while ((pos = pct.find(delimiter)) != std::string::npos)
    {
        pct.erase(0, pos + delimiter.length());
    }
    pct = pct.substr(0, pct.size() - 3);
    double pct_num = std::stod(pct);
    pct_num = pct_num;
    return pct_num;
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
    sipros::RaxportMs1Data ms1Data;
    std::string error;
    if (!sipros::readRaxportHdf5Ms1Data(hdf5Path.string(), ms1Data, error))
    {
        throw std::runtime_error(error);
    }

    MS1Scans.clear();
    scanNumberMS1ScanMap.clear();
    MS1Scans.reserve(ms1Data.scans.size());
    for (const sipros::RaxportMs1Scan &source : ms1Data.scans)
    {
        Scan scan(static_cast<int>(source.scanNumber), static_cast<float>(source.retentionTime), source.tic);
        scan.mz = source.mz;
        scan.intensity = source.intensity;
        scan.charge = source.charge;
        scan.mass = source.mz;
        MS1Scans.push_back(std::move(scan));
    }
    scanNumberMS1ScanMap.reserve(MS1Scans.size());
    for (size_t i = 0; i < MS1Scans.size(); i++)
    {
        scanNumberMS1ScanMap.insert({MS1Scans[i].scanNumber, MS1Scans.data() + i});
    }
}

void PSMfeatureExtractor::initializeFeatureVectors(sipPSM &psm)
{
    const size_t count = psm.scanNumbers.size();
    psm.isotopicPeakss = std::vector<std::vector<isotopicPeak>>(count);
    psm.isotopicPeakNumbers = std::vector<int>(count);
    psm.MS1IsotopicAbundances = std::vector<double>(count);
    psm.MS2IsotopicAbundances = std::vector<double>(count);
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

void ompForDemo()
{
    int numThreads = omp_get_max_threads();
    int numIterations = 10000;
    double sum = 0.0;

    // #pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < numIterations; i++)
    {
        double result = std::sqrt(std::pow(i, 20) + std::pow(i + 1, 20));
        sum += result;
    }

    printf("Number of threads: %d\n", numThreads);
    printf("Number of iterations: %d\n", numIterations);
    printf("Sum: %f\n", sum);
}

void PSMfeatureExtractor::extractFeaturesOfEachPSM()
{
    float topScore = 0;
    for (size_t i = 0; i < mSipPSM->isotopicPeakss.size(); i++)
    {
        mSipPSM->isotopicPeakss[i] = findIsotopicPeaks(mSipPSM->precursorScanNumbers[i],
                                                       mSipPSM->parentCharges[i],
                                                       mSipPSM->measuredParentMasses[i],
                                                       mSipPSM->calculatedParentMasses[i]);
        mSipPSM->isotopicPeakNumbers[i] = mSipPSM->isotopicPeakss[i].size();
        const std::string compositionPeptide = peptideSequenceForComposition(mSipPSM->identifiedPeptides[i]);
        mSipPSM->MS1IsotopicAbundances[i] = getSIPelementAbundanceFromMS1(compositionPeptide,
                                                                          mSipPSM->isotopicPeakss[i],
                                                                          mSipPSM->parentCharges[i]);
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

        mSipPSM->MS2IsotopicAbundances[i] = getMS2IsotopicAbundance(mSipPSM->searchNames[i]);
        mSipPSM->isotopicAbundanceDiffs[i] = mSipPSM->MS1IsotopicAbundances[i] - mSipPSM->MS2IsotopicAbundances[i];
    }
}

void PSMfeatureExtractor::writePecorlatorPin(const std::string &fileName, bool doProteinInference)
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    std::ofstream file(fileName);
    if (!file)
    {
        std::cerr << "Unable to open file for writing.\n";
    }
    file << std::fixed << std::setprecision(6);

    const size_t chunkSize = 10000;
    std::stringstream ss;
    ss << "SpecId"
       << "\t"
       << "Label"
       << "\t"
       << "ScanNr"
       << "\t"
       << "ExpMass"
       << "\t"
       << "retentiontime"
       << "\t"
       << "ranks"
       << "\t"
       << "parentCharges"
       << "\t"
       << "massErrors"
       << "\t"
       << "isotopicMassWindowShifts"
       << "\t"
       << "mzShiftFromisolationWindowCenters"
       << "\t"
       << "peptideLengths"
       << "\t"
       << "missCleavageSiteNumbers"
       << "\t"
       << "PTMnumbers"
       << "\t"
       << "isotopicPeakNumbers"
       << "\t"
       << "MS1IsotopicAbundances"
       << "\t"
       << "MS2IsotopicAbundances"
       << "\t"
       << "isotopicAbundanceDiffs"
       << "\t"
       << "WDPscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "MVHscores"
       << "\t"
       << "diffScores"
       << "\t"
       << "log10_precursorIntensities"
       << "\t"
       << "Peptide"
       << "\t"
       << "Proteins"
       << "\n";
    file << ss.str();
    ss.str(std::string());
    ss.clear();
    std::string proteinName;
    std::string peptideSeq;
    std::vector<sipPSM> &psms = sipPSMs;
    for (size_t i = 0; i < psms.size(); i++)
    {
        for (size_t j = 0; j < psms[i].fileNames.size(); j += chunkSize)
        {
            for (size_t k = j; k < std::min(j + chunkSize, psms[i].fileNames.size()); ++k)
            {
                ss << psms[i].fileNames[k] << "." << psms[i].scanNumbers[k] << "." << psms[i].ranks[k] << "\t";
                ss << (psms[i].isDecoys[k] ? -1 : 1) << "\t";
                ss << psms[i].scanNumbers[k] << "\t"
                   << psms[i].calculatedParentMasses[k] << "\t"
                   << psms[i].retentionTimes[k] << "\t"
                   << psms[i].ranks[k] << "\t"
                   << psms[i].parentCharges[k] << "\t"
                   << psms[i].massErrors[k] << "\t"
                   << psms[i].isotopicMassWindowShifts[k] << "\t"
                   << psms[i].mzShiftFromisolationWindowCenters[k] << "\t"
                   << psms[i].peptideLengths[k] << "\t"
                   << psms[i].missCleavageSiteNumbers[k] << "\t"
                   << psms[i].PTMnumbers[k] << "\t"
                   << psms[i].isotopicPeakNumbers[k] << "\t"
                   << psms[i].MS1IsotopicAbundances[k] << "\t"
                   << psms[i].MS2IsotopicAbundances[k] << "\t"
                   << psms[i].isotopicAbundanceDiffs[k] << "\t"
                   << psms[i].WDPscores[k] << "\t"
                   << psms[i].XcorrScores[k] << "\t"
                   << psms[i].MVHscores[k] << "\t"
                   << psms[i].diffScores[k] << "\t"
                   << (psms[i].precursorIntensities[k] > 0
                           ? std::log10(psms[i].precursorIntensities[k])
                           : 0)
                   << "\t";

                if (doProteinInference)
                {
                    // for percolator protein inference format
                    peptideSeq = psms[i].originalPeptides[k];
                    if (peptideSeq.front() == '[')
                        peptideSeq.insert(peptideSeq.begin(), 'n');
                    if (peptideSeq.back() == ']')
                        peptideSeq.push_back('n');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), '[', '.');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), ']', '.');
                }
                else
                    peptideSeq = psms[i].identifiedPeptides[k];
                ss << peptideSeq << "\t";

                proteinName = psms[i].proteinNames[k];
                if (doProteinInference)
                {
                    // for percolator protein inference format
                    proteinName = proteinName.substr(1, psms[i].proteinNames[k].length() - 2);
                    std::replace(proteinName.begin(), proteinName.end(), ',', '\t');
                }
                ss << proteinName;
                ss << "\n";
            }
            file << ss.str();
            ss.str(std::string()); // Clear the stringstream
        }
    }
    file.close();
}