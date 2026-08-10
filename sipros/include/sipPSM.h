#pragma once
#include <string>
#include <vector>

struct isotopicPeak
{
    double mz = 0.0;
    int charge = 0;
    double intensity = 0.0;
    int isotopeIndex = -1;
};

class alignas(64) sipPSM
{
public:
    std::string fileName;

    // One entry per PSM row. SIP direct search may use pct-specific SpecId stems.
    std::vector<std::string> fileNames;

    std::vector<int> scanNumbers;
    std::vector<int> parentCharges;
    std::vector<double> measuredParentMasses;
    std::vector<double> calculatedParentMasses;
    std::vector<double> precursorNeutronMasses;
    std::string scanType;

    std::string scoringFunction;
    std::vector<int> ranks;
    std::vector<float> scores;
    std::vector<std::string> identifiedPeptides;
    std::vector<std::string> originalPeptides;
    std::vector<std::string> nakePeptides;
    std::vector<std::string> proteinNames;

    // Precursor and score fields used for PIN feature extraction.
    std::vector<int> precursorScanNumbers;
    std::vector<double> isolationWindowCenterMZs;
    std::vector<float> retentionTimes;
    std::vector<float> MVHscores;
    std::vector<float> XcorrScores;
    std::vector<float> WDPscores;
    std::vector<double> precursorRtDiffSeconds;
    std::vector<int> ddaResidualRanks;
    std::vector<float> ddaResidualScores;

    // PIN features.
    std::vector<bool> isDecoys;
    std::vector<int> peptideLengths;
    std::vector<int> missCleavageSiteNumbers;
    std::vector<int> PTMnumbers;
    std::vector<float> diffScores;
    std::vector<int> matchedBIons;
    std::vector<int> matchedYIons;
    std::vector<int> maxConsecutiveBIons;
    std::vector<int> maxConsecutiveYIons;
    std::vector<double> mzShiftFromisolationWindowCenters;
    std::vector<int> isotopicMassWindowShifts;
    std::vector<double> massErrors;
    std::vector<std::vector<isotopicPeak>> isotopicPeakss;
    std::vector<int> isotopicPeakNumbers;
    std::vector<double> MS1IsotopeFitScores;
    std::vector<double> precursorIntensities;
    std::vector<double> MS1IsotopicAbundances;
    std::vector<double> MS2IsotopicAbundances;
    std::vector<double> isotopicAbundanceDiffs;
};
