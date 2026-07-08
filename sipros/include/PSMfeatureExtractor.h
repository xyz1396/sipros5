#pragma once
#include "Scan.h"
#include "averagine.h"
#include "isotopicPeak.h"
#include "sipPSM.h"
#include <array>
#include <omp.h>
#include <string>
#include <unordered_map>
#include <vector>

class PSMfeatureExtractor
{
public:
    PSMfeatureExtractor();

    averagine mAveragine;
    std::vector<Scan> MS1Scans;
    std::vector<sipPSM> sipPSMs;
    sipPSM *mSipPSM;
    std::unordered_map<size_t, Scan *> scanNumberMS1ScanMap;
    // N isotopic peaks on each side to consider
    const static int NisotopicPeak = 20;
    std::vector<int> vertexIXs;
    std::array<char, 2> cleavageSites = {'K', 'R'};

    void loadHdf5Ms1(const std::string &hdf5BasePath);
    void initializeFeatureVectors(sipPSM &psm);
    void extractFeaturesForPsm(const std::string &hdf5Path, sipPSM &psm);

    size_t binarySearchPeak(const Scan *mScan, double Mz, int charge);
    // filter isotopic peaks by isotopic envelope shape
    void filterIsotopicPeaks(std::vector<isotopicPeak> &isotopicPeaks, const double calculatedPrecursorMZ);
    void filterIsotopicPeaksTopN(std::vector<isotopicPeak> &isotopicPeaks, const double observedPrecursorMZ,
                                 const size_t topN);
    // return precursor scan number and isotopic peaks
    std::vector<isotopicPeak> findIsotopicPeaks(int &MS1ScanNumber,
                                                const int precursorCharge,
                                                const double observedPrecursorMass,
                                                const double calculatedPrecursorMass);
    double getSIPelementAbundanceFromMS1(const std::string &peptideSeq,
                                         const std::vector<isotopicPeak> &isotopicPeaks, const int precursorCharge);
    std::pair<int, int> getSeqLengthAndMissCleavageSiteNumber(const std::string &peptideSeq);
    int getPTMnumber(const std::string &peptideSeq);
    std::pair<int, double> getMassWindowShiftAndError(const double observedPrecursorMass,
                                                      const double calculatedPrecursorMass);
    double getMS2IsotopicAbundance(const std::string &searchName);
    void extractFeaturesOfEachPSM();
    void writePecorlatorPin(const std::string &fileName, bool doProteinInference);
};
