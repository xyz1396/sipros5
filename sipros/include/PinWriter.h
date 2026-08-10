#pragma once
#include "sipPSM.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class PinWriter
{
public:
    struct SearchSpectraPinRow
    {
        int32_t label = 0;
        int scanNumber = 0;
        int rank = 0;
        int32_t parentCharge = 0;
        int32_t isotopicMassWindowShift = 0;
        int32_t peptideLength = 0;
        int32_t missCleavageSiteNumber = 0;
        int32_t ptmNumber = 0;
        int32_t isotopicPeakNumber = 0;
        double ms1IsotopeFitScore = 0.0;
        int32_t matchedYEnvelope = 0;
        int32_t matchedBEnvelope = 0;
        double expMass = 0.0;
        double observedMass = 0.0;
        double retentionTime = 0.0;
        double massError = 0.0;
        double mzShiftFromIsolationWindowCenter = 0.0;
        double ms1IsotopicAbundance = 0.0;
        double ms2IsotopicAbundance = 0.0;
        double wdpScore = 0.0;
        double xcorrScore = 0.0;
        double mvhScore = 0.0;
        double entropyScore = 0.0;
        double cosineScore = 0.0;
        double deltaRT = 0.0;
        double log10PrecursorIntensity = 0.0;
        double precursorRtDiffSeconds = -1.0;
        std::string peptide;
        std::string proteins;
    };

    static void writePecorlatorPin(const std::string &fileName, const std::vector<sipPSM> &psms, bool doProteinInference);
    static size_t writeSearchSpectraPin(const std::string &fileName,
                                        const std::string &sampleBasename,
                                        const std::vector<SearchSpectraPinRow> &rows);
};
