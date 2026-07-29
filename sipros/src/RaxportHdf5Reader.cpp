#include "RaxportHdf5Reader.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <H5Cpp.h>
#include <hdf5.h>

#include "ms2scan.h"
#include "proNovoConfig.h"

namespace fs = std::filesystem;

namespace sipros
{
namespace
{

bool hasObject(hid_t location, const std::string &name)
{
    return H5Lexists(location, name.c_str(), H5P_DEFAULT) > 0;
}

std::string lowerExt(const std::string &path)
{
    std::string ext = fs::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext;
}

void check1D(const H5::DataSet &dataset, const std::string &name, hsize_t &size)
{
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1)
    {
        throw std::runtime_error("Dataset is not one-dimensional: " + name);
    }
    space.getSimpleExtentDims(&size, nullptr);
}

template <typename T>
std::vector<T> read1D(H5::H5File &file, const std::string &name, const H5::PredType &type)
{
    H5::DataSet dataset = file.openDataSet(name);
    hsize_t size = 0;
    check1D(dataset, name, size);
    std::vector<T> values(static_cast<size_t>(size));
    if (size > 0)
    {
        dataset.read(values.data(), type);
    }
    return values;
}

template <typename T>
std::vector<T> readSlice(H5::DataSet &dataset,
                         const H5::PredType &type,
                         hsize_t start,
                         hsize_t count)
{
    std::vector<T> values(static_cast<size_t>(count));
    if (count == 0)
    {
        return values;
    }
    H5::DataSpace fileSpace = dataset.getSpace();
    hsize_t datasetSize = 0;
    fileSpace.getSimpleExtentDims(&datasetSize, nullptr);
    if (start > datasetSize || count > datasetSize - start)
    {
        throw std::runtime_error("HDF5 peak slice is out of bounds");
    }
    H5::DataSpace memSpace(1, &count);
    fileSpace.selectHyperslab(H5S_SELECT_SET, &count, &start);
    dataset.read(values.data(), type, memSpace, fileSpace);
    return values;
}

std::string formatRt(double rt)
{
    std::ostringstream out;
    out << std::setprecision(10) << rt;
    return out.str();
}

void sortPeakTriples(std::vector<double> &mz,
                     std::vector<double> &intensity,
                     std::vector<int> &charge)
{
    if (mz.size() <= 1)
    {
        return;
    }
    std::vector<size_t> order(mz.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return mz[a] < mz[b];
    });

    std::vector<double> sortedMz;
    std::vector<double> sortedIntensity;
    std::vector<int> sortedCharge;
    sortedMz.reserve(order.size());
    sortedIntensity.reserve(order.size());
    sortedCharge.reserve(order.size());
    for (size_t idx : order)
    {
        sortedMz.push_back(mz[idx]);
        sortedIntensity.push_back(idx < intensity.size() ? intensity[idx] : 0.0);
        sortedCharge.push_back(idx < charge.size() ? charge[idx] : 0);
    }
    mz.swap(sortedMz);
    intensity.swap(sortedIntensity);
    charge.swap(sortedCharge);
}

void updateObservedMzRange(const std::vector<double> &mz)
{
    for (double value : mz)
    {
        if (value > ProNovoConfig::maxObservedMz)
        {
            ProNovoConfig::maxObservedMz = value;
        }
        if (value < ProNovoConfig::minObservedMz)
        {
            ProNovoConfig::minObservedMz = value;
        }
    }
}

void finalizeMs2Scan(MS2Scan *scan, bool useReactionChargeForScoring)
{
    scan->isMS2HighRes = (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1);
    scan->isMS1HighRes = true;
    scan->bUseReactionChargeForScoring = useReactionChargeForScoring;
    scan->iMaxCandidateCharge = 0;

    double maxChargedMass = 0.0;
    double maxNeutralMass = 0.0;
    auto considerPrecursor = [&](double mz, int charge) {
        if (mz <= 0.0 || charge <= 0)
        {
            return;
        }
        const double chargedMass = mz * charge;
        const double neutralMass = chargedMass -
                                   static_cast<double>(charge) * ProNovoConfig::getProtonMass();
        maxChargedMass = std::max(maxChargedMass, chargedMass);
        maxNeutralMass = std::max(maxNeutralMass, neutralMass);
    };

    if (useReactionChargeForScoring)
    {
        considerPrecursor(scan->dParentMZ, scan->iParentChargeState);
    }
    const size_t n = std::min(scan->iParentChargeStates.size(), scan->dParentMZs.size());
    for (size_t i = 0; i < n; ++i)
    {
        considerPrecursor(scan->dParentMZs[i], scan->iParentChargeStates[i]);
        scan->iMaxCandidateCharge = std::max(
            scan->iMaxCandidateCharge, scan->iParentChargeStates[i]);
    }
    scan->dParentMass = maxChargedMass;
    scan->dParentNeutralMass = maxNeutralMass;
}

void appendCandidatePrecursors(MS2Scan *scan,
                               long long candidateStart,
                               int candidateCount,
                               const std::vector<int> &candidateCharge,
                               const std::vector<double> &candidateMz)
{
    if (candidateStart < 0 || candidateCount <= 0)
    {
        return;
    }
    const size_t start = static_cast<size_t>(candidateStart);
    if (start >= candidateMz.size())
    {
        return;
    }
    const size_t end = std::min(candidateMz.size(), start + static_cast<size_t>(candidateCount));
    for (size_t i = start; i < end; ++i)
    {
        const int charge = i < candidateCharge.size() ? candidateCharge[i] : 0;
        const double mz = candidateMz[i];
        if (charge <= 0 || mz <= 0.0)
        {
            continue;
        }
        scan->iParentChargeStates.push_back(charge);
        scan->dParentMZs.push_back(mz);
    }
}

size_t appendNeighborhoodMs1Precursors(MS2Scan *scan,
                                       const RaxportMs1Data &ms1Data,
                                       int parentScanNumber,
                                       double isolationCenterMz,
                                       double isolationWidth,
                                       int neighborhoodRadius)
{
    if (scan == nullptr || ms1Data.scans.empty() || neighborhoodRadius < 0 ||
        isolationCenterMz <= 0.0 || isolationWidth <= 0.0 ||
        !std::isfinite(isolationCenterMz) || !std::isfinite(isolationWidth))
    {
        return 0;
    }

    const auto parent = ms1Data.scanNumberToIndex.find(parentScanNumber);
    if (parent == ms1Data.scanNumberToIndex.end())
    {
        return 0;
    }

    const size_t parentIndex = parent->second;
    const size_t radius = static_cast<size_t>(neighborhoodRadius);
    const size_t first = parentIndex > radius ? parentIndex - radius : 0;
    const size_t last = std::min(ms1Data.scans.size() - 1, parentIndex + radius);
    const double lowerMz = isolationCenterMz - isolationWidth / 2.0;
    const double upperMz = isolationCenterMz + isolationWidth / 2.0;
    size_t appended = 0;

    for (size_t scanIndex = first; scanIndex <= last; ++scanIndex)
    {
        const RaxportMs1Scan &ms1Scan = ms1Data.scans[scanIndex];
        auto peak = std::lower_bound(ms1Scan.mz.begin(), ms1Scan.mz.end(), lowerMz);
        for (; peak != ms1Scan.mz.end() && *peak <= upperMz; ++peak)
        {
            const double mz = *peak;
            if (mz <= 0.0 || !std::isfinite(mz))
            {
                continue;
            }
            const size_t peakIndex = static_cast<size_t>(std::distance(ms1Scan.mz.begin(), peak));
            const int nativeCharge = peakIndex < ms1Scan.charge.size()
                                         ? ms1Scan.charge[peakIndex]
                                         : 0;
            if (nativeCharge > 0)
            {
                scan->iParentChargeStates.push_back(nativeCharge);
                scan->dParentMZs.push_back(mz);
                ++appended;
            }
            else if (nativeCharge == 0)
            {
                for (const int guessedCharge : {2, 3, 4})
                {
                    scan->iParentChargeStates.push_back(guessedCharge);
                    scan->dParentMZs.push_back(mz);
                    ++appended;
                }
            }
        }
    }
    return appended;
}

bool requestedScan(const std::unordered_set<int> *requested, int scanNumber)
{
    return requested == nullptr || requested->find(scanNumber) != requested->end();
}

void appendMs1Data(RaxportMs1Data &data,
                   int scanNumber,
                   double retentionTime,
                   double tic,
                   std::vector<double> mz,
                   std::vector<double> intensity,
                   std::vector<int> charge)
{
    sortPeakTriples(mz, intensity, charge);
    RaxportMs1Scan scan;
    scan.scanNumber = scanNumber;
    scan.retentionTime = retentionTime;
    scan.tic = tic;
    scan.mz = std::move(mz);
    scan.intensity = std::move(intensity);
    scan.charge = std::move(charge);
    data.scanNumberToIndex[scan.scanNumber] = data.scans.size();
    data.scans.push_back(std::move(scan));
}

} // namespace

bool isRaxportHdf5Path(const std::string &path)
{
    const std::string ext = lowerExt(path);
    return ext == ".h5" || ext == ".hdf5";
}

bool readRaxportHdf5Scans(const std::string &path,
                          std::vector<MS2Scan *> &ms2Scans,
                          RaxportMs1Data *ms1Data,
                          std::string &error,
                          const std::unordered_set<int> *requestedMs2ScanNumbers,
                          const RaxportReadOptions &options)
{
    error.clear();
    ms2Scans.clear();
    if (ms1Data != nullptr)
    {
        ms1Data->clear();
    }
    if (!isRaxportHdf5Path(path))
    {
        error = "Raxport HDF5 input required (.h5 or .hdf5): " + path;
        return false;
    }
    if (options.ms1NeighborhoodRadius < 0)
    {
        error = "MS1 neighborhood radius must be non-negative";
        return false;
    }

    const bool useMs1Neighborhood =
        options.precursorSource == PrecursorSource::Ms1Neighborhood;

    H5::Exception::dontPrint();
    try
    {
        H5::H5File file(path, H5F_ACC_RDONLY);
        if (!hasObject(file.getId(), "scans") || !hasObject(file.getId(), "peaks") ||
            !hasObject(file.getId(), "reactions"))
        {
            throw std::runtime_error("missing one or more required groups: /scans, /peaks, /reactions");
        }
        if (!useMs1Neighborhood && !hasObject(file.getId(), "precursor_candidates"))
        {
            throw std::runtime_error("missing required group: /precursor_candidates");
        }

        int schemaVersion = 0;
        try
        {
            H5::Attribute attr = file.openAttribute("schema_version");
            attr.read(H5::PredType::NATIVE_INT, &schemaVersion);
        }
        catch (const H5::Exception &)
        {
            throw std::runtime_error("missing root attribute schema_version");
        }
        constexpr int kSupportedSchemaVersion = 6;
        if (schemaVersion != kSupportedSchemaVersion)
        {
            throw std::runtime_error("unsupported Raxport HDF5 schema_version " +
                                     std::to_string(schemaVersion) +
                                     "; expected 6");
        }

        const std::vector<int> scanNumber = read1D<int>(file, "/scans/scan_number", H5::PredType::NATIVE_INT);
        const std::vector<int> msOrder = read1D<int>(file, "/scans/ms_order", H5::PredType::NATIVE_INT);
        const std::vector<double> retentionTime = read1D<double>(file, "/scans/retention_time", H5::PredType::NATIVE_DOUBLE);
        const std::vector<double> tic = read1D<double>(file, "/scans/tic", H5::PredType::NATIVE_DOUBLE);
        const std::vector<int> parentScanNumber = read1D<int>(file, "/scans/parent_scan_number", H5::PredType::NATIVE_INT);
        const std::vector<long long> reactionStart = read1D<long long>(file, "/scans/reaction_start", H5::PredType::NATIVE_LLONG);
        const std::vector<int> reactionCount = read1D<int>(file, "/scans/reaction_count", H5::PredType::NATIVE_INT);
        const std::vector<long long> peakStart = read1D<long long>(file, "/scans/peak_start", H5::PredType::NATIVE_LLONG);
        const std::vector<int> peakCount = read1D<int>(file, "/scans/peak_count", H5::PredType::NATIVE_INT);

        const size_t nScans = scanNumber.size();
        if (msOrder.size() != nScans || retentionTime.size() != nScans || tic.size() != nScans ||
            parentScanNumber.size() != nScans || reactionStart.size() != nScans || reactionCount.size() != nScans ||
            peakStart.size() != nScans || peakCount.size() != nScans)
        {
            throw std::runtime_error("/scans datasets have inconsistent lengths");
        }

        const std::vector<double> reactionPrecursorMass = read1D<double>(file, "/reactions/precursor_mass", H5::PredType::NATIVE_DOUBLE);
        const std::vector<int> reactionChargeState = read1D<int>(file, "/reactions/charge_state", H5::PredType::NATIVE_INT);
        if (reactionChargeState.size() != reactionPrecursorMass.size())
        {
            throw std::runtime_error("/reactions datasets have inconsistent lengths");
        }

        std::vector<long long> reactionCandidateStart;
        std::vector<int> reactionCandidateCount;
        std::vector<double> reactionIsolationWidth;
        std::vector<double> reactionIsolationOffset;
        std::vector<int> candidateCharge;
        std::vector<double> candidateMz;
        if (useMs1Neighborhood)
        {
            reactionIsolationWidth = read1D<double>(file, "/reactions/isolation_width", H5::PredType::NATIVE_DOUBLE);
            reactionIsolationOffset = read1D<double>(file, "/reactions/isolation_width_offset", H5::PredType::NATIVE_DOUBLE);
            if (reactionIsolationWidth.size() != reactionPrecursorMass.size() ||
                reactionIsolationOffset.size() != reactionPrecursorMass.size())
            {
                throw std::runtime_error("/reactions isolation-window datasets have inconsistent lengths");
            }
        }
        else
        {
            reactionCandidateStart = read1D<long long>(file, "/reactions/candidate_start", H5::PredType::NATIVE_LLONG);
            reactionCandidateCount = read1D<int>(file, "/reactions/candidate_count", H5::PredType::NATIVE_INT);
            if (reactionCandidateStart.size() != reactionPrecursorMass.size() ||
                reactionCandidateCount.size() != reactionPrecursorMass.size())
            {
                throw std::runtime_error("/reactions candidate datasets have inconsistent lengths");
            }

            candidateCharge = read1D<int>(file, "/precursor_candidates/charge", H5::PredType::NATIVE_INT);
            candidateMz = read1D<double>(file, "/precursor_candidates/mz", H5::PredType::NATIVE_DOUBLE);
            // These v6 fields describe candidate provenance and confidence.
            // Validate their alignment, but search every candidate Raxport selected.
            const std::vector<int> candidateChargeSource =
                read1D<int>(file, "/precursor_candidates/charge_source", H5::PredType::NATIVE_INT);
            const std::vector<int> candidateIsotopeMatchCount =
                read1D<int>(file, "/precursor_candidates/isotope_match_count", H5::PredType::NATIVE_INT);
            if (candidateCharge.size() != candidateMz.size() ||
                candidateChargeSource.size() != candidateMz.size() ||
                candidateIsotopeMatchCount.size() != candidateMz.size())
            {
                throw std::runtime_error("Raxport schema 6 precursor-candidate datasets have inconsistent lengths");
            }
        }

        H5::DataSet peakMzDataset = file.openDataSet("/peaks/mz");
        H5::DataSet peakIntensityDataset = file.openDataSet("/peaks/intensity");
        hsize_t peakMzSize = 0;
        hsize_t peakIntensitySize = 0;
        check1D(peakMzDataset, "/peaks/mz", peakMzSize);
        check1D(peakIntensityDataset, "/peaks/intensity", peakIntensitySize);
        if (peakMzSize != peakIntensitySize)
        {
            throw std::runtime_error("/peaks/mz and /peaks/intensity have inconsistent lengths");
        }

        bool hasPeakCharge = hasObject(file.getId(), "peaks/charge");
        H5::DataSet peakChargeDataset;
        hsize_t peakChargeSize = 0;
        if (hasPeakCharge)
        {
            peakChargeDataset = file.openDataSet("/peaks/charge");
            check1D(peakChargeDataset, "/peaks/charge", peakChargeSize);
            if (peakChargeSize != peakMzSize)
            {
                hasPeakCharge = false;
            }
        }

        // The neighborhood mode needs both preceding and following MS1 scans,
        // so load the full MS1 series before constructing any MS2 candidates.
        RaxportMs1Data localMs1Data;
        RaxportMs1Data *neighborhoodMs1Data = nullptr;
        if (useMs1Neighborhood)
        {
            neighborhoodMs1Data = ms1Data != nullptr ? ms1Data : &localMs1Data;
            for (size_t i = 0; i < nScans; ++i)
            {
                if (msOrder[i] != 1)
                {
                    continue;
                }
                if (peakStart[i] < 0 || peakCount[i] <= 0)
                {
                    appendMs1Data(*neighborhoodMs1Data, scanNumber[i], retentionTime[i], tic[i],
                                  {}, {}, {});
                    continue;
                }
                const hsize_t start = static_cast<hsize_t>(peakStart[i]);
                const hsize_t count = static_cast<hsize_t>(peakCount[i]);
                if (start > peakMzSize || count > peakMzSize - start)
                {
                    throw std::runtime_error("scan " + std::to_string(scanNumber[i]) +
                                             " has an out-of-bounds peak slice");
                }
                std::vector<double> mz = readSlice<double>(
                    peakMzDataset, H5::PredType::NATIVE_DOUBLE, start, count);
                std::vector<double> intensity = readSlice<double>(
                    peakIntensityDataset, H5::PredType::NATIVE_DOUBLE, start, count);
                std::vector<int> charge;
                if (hasPeakCharge)
                {
                    charge = readSlice<int>(
                        peakChargeDataset, H5::PredType::NATIVE_INT, start, count);
                }
                else
                {
                    charge.assign(static_cast<size_t>(count), 0);
                }
                appendMs1Data(*neighborhoodMs1Data, scanNumber[i], retentionTime[i], tic[i],
                              std::move(mz), std::move(intensity), std::move(charge));
            }
        }

        ms2Scans.reserve(requestedMs2ScanNumbers == nullptr ? nScans / 2 : requestedMs2ScanNumbers->size());
        for (size_t i = 0; i < nScans; ++i)
        {
            const int order = msOrder[i];
            if (order != 1 && order != 2)
            {
                continue;
            }
            if (peakStart[i] < 0 || peakCount[i] <= 0)
            {
                continue;
            }
            const hsize_t start = static_cast<hsize_t>(peakStart[i]);
            const hsize_t count = static_cast<hsize_t>(peakCount[i]);
            if (start > peakMzSize || count > peakMzSize - start)
            {
                throw std::runtime_error("scan " + std::to_string(scanNumber[i]) + " has an out-of-bounds peak slice");
            }

            if (order == 1)
            {
                if (useMs1Neighborhood)
                {
                    continue;
                }
                if (ms1Data == nullptr)
                {
                    continue;
                }
                std::vector<double> mz = readSlice<double>(peakMzDataset, H5::PredType::NATIVE_DOUBLE, start, count);
                std::vector<double> intensity = readSlice<double>(peakIntensityDataset, H5::PredType::NATIVE_DOUBLE, start, count);
                std::vector<int> charge;
                if (hasPeakCharge)
                {
                    charge = readSlice<int>(peakChargeDataset, H5::PredType::NATIVE_INT, start, count);
                }
                else
                {
                    charge.assign(static_cast<size_t>(count), 0);
                }
                appendMs1Data(*ms1Data, scanNumber[i], retentionTime[i], tic[i], std::move(mz), std::move(intensity), std::move(charge));
                continue;
            }

            if (!requestedScan(requestedMs2ScanNumbers, scanNumber[i]))
            {
                continue;
            }

            std::unique_ptr<MS2Scan> scan(new MS2Scan);
            scan->iScanId = scanNumber[i];
            scan->iParentScanID = parentScanNumber[i];
            scan->setRTime(formatRt(retentionTime[i]));
            scan->isMS1HighRes = true;
            scan->isMS2HighRes = (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1);

            const long long rxnStart = reactionStart[i];
            const int rxnCount = reactionCount[i];
            bool havePrimaryReaction = false;
            if (rxnStart >= 0 && rxnCount > 0)
            {
                const size_t startRxn = static_cast<size_t>(rxnStart);
                const size_t endRxn = std::min(reactionPrecursorMass.size(), startRxn + static_cast<size_t>(rxnCount));
                for (size_t r = startRxn; r < endRxn; ++r)
                {
                    if (!havePrimaryReaction)
                    {
                        scan->dParentMZ = reactionPrecursorMass[r];
                        scan->iParentChargeState = r < reactionChargeState.size() ? reactionChargeState[r] : 0;
                        havePrimaryReaction = true;
                    }
                    if (useMs1Neighborhood)
                    {
                        const double isolationCenterMz =
                            reactionPrecursorMass[r] + reactionIsolationOffset[r];
                        appendNeighborhoodMs1Precursors(
                            scan.get(), *neighborhoodMs1Data, parentScanNumber[i],
                            isolationCenterMz, reactionIsolationWidth[r],
                            options.ms1NeighborhoodRadius);
                    }
                    else
                    {
                        appendCandidatePrecursors(scan.get(),
                                                  reactionCandidateStart[r],
                                                  reactionCandidateCount[r],
                                                  candidateCharge,
                                                  candidateMz);
                    }
                }
            }
            if (scan->iParentChargeStates.empty())
            {
                continue;
            }

            scan->vdMZ = readSlice<double>(peakMzDataset, H5::PredType::NATIVE_DOUBLE, start, count);
            scan->vdIntensity = readSlice<double>(peakIntensityDataset, H5::PredType::NATIVE_DOUBLE, start, count);
            if (hasPeakCharge)
            {
                scan->viCharge = readSlice<int>(peakChargeDataset, H5::PredType::NATIVE_INT, start, count);
            }
            else
            {
                scan->viCharge.assign(static_cast<size_t>(count), 0);
            }
            sortPeakTriples(scan->vdMZ, scan->vdIntensity, scan->viCharge);
            updateObservedMzRange(scan->vdMZ);
            finalizeMs2Scan(scan.get(), !useMs1Neighborhood);
            ms2Scans.push_back(scan.release());
        }
        return true;
    }
    catch (const H5::Exception &ex)
    {
        error = "Unable to read Raxport HDF5 file '" + path + "': " + ex.getDetailMsg();
    }
    catch (const std::exception &ex)
    {
        error = "Unable to read Raxport HDF5 file '" + path + "': " + ex.what();
    }

    for (MS2Scan *scan : ms2Scans)
    {
        delete scan;
    }
    ms2Scans.clear();
    if (ms1Data != nullptr)
    {
        ms1Data->clear();
    }
    return false;
}

bool readRaxportHdf5Ms1Data(const std::string &path,
                            RaxportMs1Data &ms1Data,
                            std::string &error)
{
    std::vector<MS2Scan *> ignored;
    // Feature extraction needs only MS1 scans.  Use the Regular reader policy
    // so this helper does not accidentally require the legacy
    // /precursor_candidates group merely to expose /peaks data.
    RaxportReadOptions options;
    options.precursorSource = PrecursorSource::Ms1Neighborhood;
    const bool ok = readRaxportHdf5Scans(
        path, ignored, &ms1Data, error, nullptr, options);
    for (MS2Scan *scan : ignored)
    {
        delete scan;
    }
    ignored.clear();
    return ok;
}

} // namespace sipros
