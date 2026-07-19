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
    Ms1Neighborhood
};

struct RaxportReadOptions
{
    PrecursorSource precursorSource = PrecursorSource::RaxportCandidates;
    int ms1NeighborhoodRadius = 2;
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
