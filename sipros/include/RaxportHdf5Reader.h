#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class MS2Scan;

namespace sipros
{

enum class PrecursorSource
{
    RaxportCandidates,
    IsolationWindow
};

struct RaxportReadOptions
{
    PrecursorSource precursorSource = PrecursorSource::RaxportCandidates;
    // Regular search creates candidates from the isolation window only. This
    // radius is used later to match already-scored PSMs to nearby MS1 scans.
    int precursorMatchScanRadius = 5;
};

struct RaxportMs1Scan
{
    int scanNumber = 0;
    double retentionTime = 0.0;
    double tic = 0.0;
    std::vector<double> mz;
    std::vector<double> intensity;
    std::vector<int> charge;
};

struct RaxportMs1Data
{
    std::vector<RaxportMs1Scan> scans;
    std::unordered_map<int, size_t> scanNumberToIndex;

    void clear()
    {
        scans.clear();
        scanNumberToIndex.clear();
    }
};

struct RaxportPrecursorMatch
{
    bool found = false;
    int ms1ScanNumber = 0;
    double observedNeutralMass = 0.0;
    double rtDiffSeconds = -1.0;
};

// Match an already-scored theoretical precursor against parent MS1 +/- radius.
// Positive native peak charges must agree; unassigned peaks can support z2-z4.
RaxportPrecursorMatch findRaxportPrecursorMatch(
    const RaxportMs1Data &ms1Data,
    int parentScanNumber,
    double ms2RetentionTimeMinutes,
    int precursorCharge,
    double calculatedNeutralMass,
    double precursorNeutronMass,
    const std::vector<int> &isotopeWindows,
    double neutralMassTolerance,
    int scanRadius = 5);

bool isRaxportHdf5Path(const std::string &path);

bool readRaxportHdf5Scans(const std::string &path,
                          std::vector<MS2Scan *> &ms2Scans,
                          RaxportMs1Data *ms1Data,
                          std::string &error,
                          const std::unordered_set<int> *requestedMs2ScanNumbers = nullptr,
                          const RaxportReadOptions &options = RaxportReadOptions{});

bool readRaxportHdf5Ms1Data(const std::string &path,
                            RaxportMs1Data &ms1Data,
                            std::string &error);

} // namespace sipros
