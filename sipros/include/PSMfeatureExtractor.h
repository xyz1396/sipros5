#pragma once
#include "RaxportHdf5Reader.h"
#include "averagine.h"
#include "sipPSM.h"
#include <array>
#include <functional>
#include <omp.h>
#include <string>
#include <utility>
#include <vector>


class PSMfeatureExtractor
{
public:
    PSMfeatureExtractor();

    averagine mAveragine;
    sipros::RaxportMs1Data ms1Data;
    std::vector<sipPSM> sipPSMs;
    sipPSM *mSipPSM;
    std::array<char, 2> cleavageSites = {'K', 'R'};
    static constexpr int Ms1IsotopeWindow = 10;

    void loadHdf5Ms1(const std::string &hdf5BasePath);
    void initializeFeatureVectors(sipPSM &psm);
    void extractFeaturesForPsm(const std::string &hdf5Path, sipPSM &psm);

    std::pair<int, int> getSeqLengthAndMissCleavageSiteNumber(const std::string &peptideSeq);
    int getPTMnumber(const std::string &peptideSeq);
    std::pair<int, double> getMassWindowShiftAndError(const double observedPrecursorMass,
                                                      const double calculatedPrecursorMass);
    struct Ms1AbundanceResult
    {
        double abundancePct = 0.0;
        int isotopicPeakCount = 0;
    };

    static std::string peptideBodyWithPtms(const std::string &decorated);
    static int countMissCleavage(const std::string &naked);
    static int countPTM(const std::string &decorated);
    static int sipAtomIndex(const std::string &sipAtom);
    static int sipNominalShiftPerAtom(const std::string &sipAtom);
    static double expectedNaturalNominalShiftExceptTarget(const std::array<int, 6> &atomCounts,
                                                          int targetAtomIndex,
                                                          int targetIsotopeIndex);
    static int ms1PeakCharge(const sipros::RaxportMs1Scan &scan, size_t idx);
    static size_t findMs1Peak(const sipros::RaxportMs1Scan &scan,
                              double targetMz,
                              const std::function<double(double)> &mzToleranceDaAt,
                              int requiredCharge = -1);
    static std::vector<isotopicPeak> findMs1IsotopicPeaks(const sipros::RaxportMs1Data *ms1Data,
                                                          int &ms1ScanNumber,
                                                          int precursorCharge,
                                                          double monoPrecursorMz,
                                                          double matchedPrecursorMz,
                                                          int targetNominalShift,
                                                          const std::function<double(double)> &mzToleranceDaAt);
    static Ms1AbundanceResult getSIPelementAbundanceFromMS1Peaks(const std::vector<isotopicPeak> &peaks,
                                                                 double baseMass,
                                                                 const std::string &peptide,
                                                                 int precursorCharge,
                                                                 const std::string &sipAtom);

    void extractFeaturesOfEachPSM();
};
