#ifndef SIPROS_SPECTRA_INDEX_H
#define SIPROS_SPECTRA_INDEX_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sipros
{

struct SpectraIndexMetadata
{
	// A consolidated index stores the exact abundance on every record.
	static constexpr double MixedSipAbundancePct = -1.0;

	std::string chemistryProfileId;
	std::string recordKind;
	double targetSipAbundancePct = 0.0;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double probabilityCutoff = 0.0;
	double generationPpmTolerance = 0.0;
	uint64_t minimumMatchedEnvelopes = 0;
	uint32_t envelopeTopN = 3;
	int label = 1;
};

struct SpectraIndexPrecursorPeak
{
	double mz = 0.0;
	double intensity = 0.0;
};

struct SpectraIndexFragmentPeakInput
{
	double mz = 0.0;
	float theoreticalIntensity = 0.0F;
	float experimentalIntensity = 0.0F;
	uint32_t ionPosition = 0;
	uint8_t ionKind = 0;
	uint8_t reserved[3] = {0, 0, 0};
};

// SFI stores fragment m/z as a fixed-point 0.001-Da bin.  Position occupies
// the upper byte of packedMzPosition; the sign bit of theoreticalBits records
// b/y kind while the remaining bits preserve the original non-negative float
// exactly. Experimental intensities live in a sparse cold sidecar because the
// MVH hot path does not consume them.
struct SpectraIndexFragmentPeak
{
	static constexpr uint32_t MassMask = 0x00ffffffU;
	static constexpr uint32_t KindMask = 0x80000000U;

	uint32_t packedMzPosition = 0;
	uint32_t theoreticalBits = 0;

	uint32_t mzBin() const { return packedMzPosition & MassMask; }
	uint16_t ionPosition() const {
		return static_cast<uint16_t>(packedMzPosition >> 24U);
	}
	uint8_t ionKind() const {
		return static_cast<uint8_t>(
			(theoreticalBits & KindMask) != 0 ? 'y' : 'b');
	}
	float theoreticalIntensity() const {
		const uint32_t magnitudeBits = theoreticalBits & ~KindMask;
		float value = 0.0F;
		std::memcpy(&value, &magnitudeBits, sizeof(value));
		return value;
	}
	double mz() const { return static_cast<double>(mzBin()) * 0.001; }
};

class SpectraIndexExperimentalCursor
{
public:
	float next()
	{
		const uint64_t fragment = fragmentIndex_++;
		const bool present =
			(presenceBits_[fragment >> 6U] &
			 (uint64_t{1} << (fragment & 63U))) != 0;
		return present ? *values_++ : 0.0F;
	}

private:
	friend class SpectraIndex;
	SpectraIndexExperimentalCursor(const uint64_t *presenceBits,
		uint64_t fragmentIndex, const float *values)
		: presenceBits_(presenceBits), fragmentIndex_(fragmentIndex),
		  values_(values) {}

	const uint64_t *presenceBits_ = nullptr;
	uint64_t fragmentIndex_ = 0;
	const float *values_ = nullptr;
};

struct SpectraIndexRecordInput
{
	std::string psmId;
	std::string peptide;
	std::string proteins;
	double retentionMinutes = 0.0;
	double sipAbundancePct = 0.0;
	int charge = 1;
	std::vector<SpectraIndexPrecursorPeak> precursors;
	std::vector<SpectraIndexFragmentPeakInput> fragments;
};

struct SpectraIndexRecord
{
	double topPrecursorMz = 0.0;
	double sumPrecursorIntensity = 0.0;
	double retentionMinutes = 0.0;
	double sipAbundancePct = 0.0;
	uint32_t precursorOffset = 0;
	uint32_t fragmentOffset = 0;
	uint32_t experimentalOffset = 0;
	uint32_t psmIdOffset = 0;
	uint32_t peptideOffset = 0;
	uint32_t proteinsOffset = 0;
	uint32_t generationOrdinal = 0;
	uint16_t precursorCount = 0;
	uint16_t fragmentCount = 0;
	uint16_t psmIdSize = 0;
	uint16_t peptideSize = 0;
	uint16_t proteinsSize = 0;
	int16_t charge = 1;
};

struct SpectraIndexFragmentPosting
{
	uint32_t packed = 0;

	uint32_t massBin() const { return packed >> 8; }
	uint8_t localRecordId() const { return static_cast<uint8_t>(packed & 0xffU); }
};

struct SpectraIndexRtBin
{
	uint32_t rtBin = 0;
	uint32_t postingOffset = 0;
};

struct SpectraIndexBuildStats
{
	uint64_t recordCount = 0;
	uint64_t precursorCount = 0;
	uint64_t fragmentCount = 0;
	uint64_t experimentalValueCount = 0;
	uint64_t productPostingCount = 0;
	uint64_t rtBinCount = 0;
	uint64_t stringBytes = 0;
	uint64_t fileBytes = 0;
	uint32_t blockCount = 0;
	uint32_t threadsUsed = 1;
	double compactValidateSeconds = 0.0;
	double orderSeconds = 0.0;
	double reserveSeconds = 0.0;
	double flattenSeconds = 0.0;
	double productIndexSeconds = 0.0;
	double layoutChecksumSeconds = 0.0;
	double writeSeconds = 0.0;
	double totalSeconds = 0.0;
};

class SpectraIndex
{
public:
	SpectraIndex();
	~SpectraIndex();

	SpectraIndex(const SpectraIndex &) = delete;
	SpectraIndex &operator=(const SpectraIndex &) = delete;
	SpectraIndex(SpectraIndex &&other) noexcept;
	SpectraIndex &operator=(SpectraIndex &&other) noexcept;

	// Consumes per-record strings and peak buffers while flattening them, which
	// keeps peak memory near one compact representation for large indexes.
	static bool write(const std::string &path,
					  const SpectraIndexMetadata &metadata,
					  std::vector<SpectraIndexRecordInput> &records,
					  std::string &error,
					  int threads = 1,
					  SpectraIndexBuildStats *stats = nullptr,
					  bool inputAlreadyCompacted = false);

	bool load(const std::string &path, std::string &error);
	void close();

	const SpectraIndexMetadata &metadata() const { return metadata_; }
	uint64_t recordCount() const { return recordCount_; }
	uint64_t fragmentCount() const { return fragmentCount_; }
	uint64_t rtBinCount() const { return rtBinCount_; }
	uint32_t precursorBlockCount() const { return precursorBlockCount_; }
	static constexpr uint32_t recordBlockCapacity() { return RecordsPerBlock; }
	static constexpr double rtBinWidthMinutes() { return RtBinWidthMinutes; }

	const SpectraIndexRecord &record(uint32_t recordId) const;
	std::string_view psmId(uint32_t recordId) const;
	std::string_view peptide(uint32_t recordId) const;
	std::string_view proteins(uint32_t recordId) const;
	std::pair<const SpectraIndexPrecursorPeak *, const SpectraIndexPrecursorPeak *>
	precursors(uint32_t recordId) const;
	std::pair<const SpectraIndexFragmentPeak *, const SpectraIndexFragmentPeak *>
	fragments(uint32_t recordId) const;
	SpectraIndexExperimentalCursor experimentalIntensities(
		uint32_t recordId) const;

	// Records are sorted by their precursor-envelope apex m/z.
	std::pair<uint32_t, uint32_t> precursorMzRange(double lowerMz,
												 double upperMz) const;
	uint32_t precursorBlockForRecord(uint32_t recordId) const;
	uint32_t precursorBlockRecordBegin(uint32_t block) const;
	uint32_t postingRecordId(uint32_t block,
								 const SpectraIndexFragmentPosting &posting) const;
	std::pair<const SpectraIndexRtBin *, const SpectraIndexRtBin *>
	rtBins(uint32_t block, double lowerRtMinutes, double upperRtMinutes) const;
	std::pair<const SpectraIndexFragmentPosting *,
			  const SpectraIndexFragmentPosting *>
	fragmentRange(uint32_t block, const SpectraIndexRtBin &rtBin,
				  double lowerMz, double upperMz) const;
	std::pair<const SpectraIndexFragmentPosting *,
			  const SpectraIndexFragmentPosting *>
	productPostings(uint32_t block, const SpectraIndexRtBin &rtBin) const;

private:
	struct Header;
	std::string_view stringAt(uint64_t offset, uint32_t size) const;
	void moveFrom(SpectraIndex &other) noexcept;

	static constexpr uint32_t RecordsPerBlock = 256;
	static constexpr double FragmentBinWidth = 0.001;
	static constexpr double RtBinWidthMinutes = 5.0;

	void *mapping_ = nullptr;
	size_t mappingSize_ = 0;
	int mappingFd_ = -1;
	const SpectraIndexRecord *records_ = nullptr;
	const SpectraIndexPrecursorPeak *precursors_ = nullptr;
	const SpectraIndexFragmentPeak *fragments_ = nullptr;
	const uint64_t *experimentalPresenceBits_ = nullptr;
	const float *experimentalValues_ = nullptr;
	const uint64_t *blockRtBinOffsets_ = nullptr;
	const uint64_t *blockProductOffsets_ = nullptr;
	const SpectraIndexRtBin *rtBins_ = nullptr;
	const SpectraIndexFragmentPosting *productPostings_ = nullptr;
	const char *strings_ = nullptr;
	uint64_t recordCount_ = 0;
	uint64_t precursorCount_ = 0;
	uint64_t fragmentCount_ = 0;
	uint64_t experimentalValueCount_ = 0;
	uint64_t rtBinCount_ = 0;
	uint64_t productPostingCount_ = 0;
	uint64_t stringBytes_ = 0;
	uint32_t precursorBlockCount_ = 0;
	SpectraIndexMetadata metadata_;
};

} // namespace sipros

#endif // SIPROS_SPECTRA_INDEX_H
