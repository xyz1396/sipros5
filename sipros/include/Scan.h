#pragma once

#include <utility>
#include <vector>

struct alignas(64) Scan
{
    size_t scanNumber = 0;
    float retentionTime = 0;
    double TIC = 0;
    std::vector<double> mz;
    std::vector<double> mass;
    std::vector<double> intensity;

    int precursorScanNumber = 0;
    int precursorCharge = 0;
    double isolationWindowCenterMZ = 0;
    std::vector<int> precursorCharges;
    std::vector<double> precursorMZs;

    std::vector<int> resolution;
    std::vector<float> baseLine;
    std::vector<float> signalToNoise;
    std::vector<int> charge;

    Scan() = default;
    Scan(int mScanNumber, float mRetentionTime, double mTIC)
        : scanNumber(mScanNumber), retentionTime(mRetentionTime), TIC(mTIC) {}
    Scan(int mScanNumber, float mRetentionTime, double mTIC, int mPrecursorScanNumber,
         int mPrecursorCharge, double mIsolationWindowCenterMZ,
         std::vector<int> mPrecursorCharges, std::vector<double> mPrecursorMZs)
        : scanNumber(mScanNumber), retentionTime(mRetentionTime), TIC(mTIC),
          precursorScanNumber(mPrecursorScanNumber), precursorCharge(mPrecursorCharge),
          isolationWindowCenterMZ(mIsolationWindowCenterMZ),
          precursorCharges(std::move(mPrecursorCharges)), precursorMZs(std::move(mPrecursorMZs)) {}
};
