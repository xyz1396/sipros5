#include "spectraindex.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <new>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef _WIN32
#include "windows_posix_compat.h"
#else
#include <unistd.h>
#include <sys/mman.h>
#endif
#include <omp.h>


namespace sipros
{

static int closeFileDescriptor(int fd)
{
#ifdef _WIN32
	return sipros_close(fd);
#else
	return ::close(fd);
#endif
}

constexpr std::array<char, 8> Magic{{'S', 'I', 'P', 'S', 'F', 'I', '0', '6'}};
constexpr uint32_t Version = 6;
constexpr uint32_t EndianMarker = 0x01020304U;
constexpr uint64_t Alignment = 64;
constexpr uint64_t FnvOffsetBasis = 1469598103934665603ULL;
constexpr uint64_t FnvPrime = 1099511628211ULL;
constexpr double FragmentBinWidth = 0.001;
constexpr double RtBinWidthMinutes = 5.0;
constexpr uint32_t PackedMassBits = 24;
constexpr uint32_t PackedMassMax = (1U << PackedMassBits) - 1U;

template <typename T>
class RawArray
{
public:
	explicit RawArray(size_t size) : size_(size)
	{
		static_assert(std::is_trivially_copyable<T>::value,
			"SFI raw arrays require trivially copyable elements");
		if (size_ > std::numeric_limits<size_t>::max() / sizeof(T))
			throw std::bad_array_new_length();
		if (size_ != 0)
			data_ = static_cast<T *>(::operator new(size_ * sizeof(T)));
	}

	~RawArray() { ::operator delete(data_); }
	RawArray(const RawArray &) = delete;
	RawArray &operator=(const RawArray &) = delete;

	T *data() { return data_; }
	const T *data() const { return data_; }
	size_t size() const { return size_; }
	T &operator[](size_t index) { return data_[index]; }
	const T &operator[](size_t index) const { return data_[index]; }

private:
	T *data_ = nullptr;
	size_t size_ = 0;
};

template <typename Function>
void parallelFor(size_t count, uint32_t requestedThreads, Function function)
{
	if (count == 0)
		return;
	const uint32_t workers = static_cast<uint32_t>(std::min<size_t>(
		count, static_cast<size_t>(std::max<uint32_t>(1, requestedThreads))));
	if (workers == 1)
	{
		for (size_t index = 0; index < count; ++index)
			function(index);
		return;
	}
	std::atomic<size_t> next{0};
	std::vector<std::thread> threads;
	threads.reserve(workers);
	for (uint32_t worker = 0; worker < workers; ++worker)
	{
		threads.emplace_back([&]()
		{
			while (true)
			{
				const size_t index = next.fetch_add(1, std::memory_order_relaxed);
				if (index >= count)
					break;
				function(index);
			}
		});
	}
	for (std::thread &thread : threads)
		thread.join();
}

template <typename Function>
uint64_t parallelXor(size_t count, uint32_t requestedThreads, Function function)
{
	if (count == 0)
		return 0;
	const uint32_t workers = static_cast<uint32_t>(std::min<size_t>(
		count, static_cast<size_t>(std::max<uint32_t>(1, requestedThreads))));
	std::vector<uint64_t> partial(workers, 0);
	std::atomic<size_t> next{0};
	std::vector<std::thread> threads;
	threads.reserve(workers);
	for (uint32_t worker = 0; worker < workers; ++worker)
	{
		threads.emplace_back([&, worker]()
		{
			uint64_t local = 0;
			while (true)
			{
				const size_t index = next.fetch_add(1, std::memory_order_relaxed);
				if (index >= count)
					break;
				local ^= function(index);
			}
			partial[worker] = local;
		});
	}
	for (std::thread &thread : threads)
		thread.join();
	uint64_t combined = 0;
	for (uint64_t value : partial)
		combined ^= value;
	return combined;
}

template <typename T, typename Compare>
void parallelSort(std::vector<T> &values, uint32_t requestedThreads,
				  Compare compare)
{
	if (values.size() < 2)
		return;
	const size_t workers = std::min<size_t>(values.size(),
		static_cast<size_t>(std::max<uint32_t>(1, requestedThreads)));
	if (workers == 1 || values.size() < 4096)
	{
		std::sort(values.begin(), values.end(), compare);
		return;
	}
	const size_t chunkSize = (values.size() + workers - 1) / workers;
	parallelFor(workers, static_cast<uint32_t>(workers), [&](size_t chunk)
	{
		const size_t begin = std::min(values.size(), chunk * chunkSize);
		const size_t end = std::min(values.size(), begin + chunkSize);
		std::sort(values.begin() + static_cast<std::ptrdiff_t>(begin),
			values.begin() + static_cast<std::ptrdiff_t>(end), compare);
	});

	std::vector<T> buffer(values.size());
	bool valuesAreSource = true;
	for (size_t width = chunkSize; width < values.size(); width *= 2)
	{
		const size_t pairCount = (values.size() + 2 * width - 1) / (2 * width);
		parallelFor(pairCount, static_cast<uint32_t>(workers), [&](size_t pair)
		{
			const size_t begin = pair * 2 * width;
			const size_t middle = std::min(values.size(), begin + width);
			const size_t end = std::min(values.size(), begin + 2 * width);
			const T *source = valuesAreSource ? values.data() : buffer.data();
			T *destination = valuesAreSource ? buffer.data() : values.data();
			std::merge(source + begin, source + middle,
				source + middle, source + end, destination + begin, compare);
		});
		valuesAreSource = !valuesAreSource;
		if (width > values.size() / 2)
			break;
	}
	if (!valuesAreSource)
		values.swap(buffer);
}

static uint64_t alignUp(uint64_t value, uint64_t alignment)
{
	return (value + alignment - 1) & ~(alignment - 1);
}

static bool checkedAdd(uint64_t a, uint64_t b, uint64_t &out)
{
	if (a > std::numeric_limits<uint64_t>::max() - b)
		return false;
	out = a + b;
	return true;
}

static bool checkedMultiply(uint64_t a, uint64_t b, uint64_t &out)
{
	if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a)
		return false;
	out = a * b;
	return true;
}

static void hashBytes(uint64_t &hash, const void *data, uint64_t size)
{
	const auto *bytes = static_cast<const unsigned char *>(data);
	for (uint64_t i = 0; i < size; ++i)
	{
		hash ^= bytes[i];
		hash *= FnvPrime;
	}
}

template <typename Peak, typename Score>
void retainTopNInOriginalOrder(std::vector<Peak> &peaks,
							   size_t topN, Score score)
{
	if (topN == 0 || peaks.size() <= topN)
		return;
	std::vector<size_t> order(peaks.size());
	for (size_t i = 0; i < order.size(); ++i)
		order[i] = i;
	std::stable_sort(order.begin(), order.end(),
		[&](size_t left, size_t right)
		{
			const double leftScore = score(peaks[left]);
			const double rightScore = score(peaks[right]);
			if (leftScore != rightScore)
				return leftScore > rightScore;
			return left < right;
		});
	order.resize(topN);
	std::sort(order.begin(), order.end());
	std::vector<Peak> retained;
	retained.reserve(topN);
	for (size_t index : order)
		retained.push_back(std::move(peaks[index]));
	peaks = std::move(retained);
}

void retainCompactEnvelopePeaks(SpectraIndexRecordInput &record,
								uint32_t envelopeTopN)
{
	retainTopNInOriginalOrder(record.precursors, envelopeTopN,
		[](const SpectraIndexPrecursorPeak &peak)
		{ return std::isfinite(peak.intensity) ? peak.intensity :
			-std::numeric_limits<double>::infinity(); });

	std::unordered_map<uint64_t, std::vector<size_t>> envelopeMembers;
	envelopeMembers.reserve(record.fragments.size());
	for (size_t i = 0; i < record.fragments.size(); ++i)
	{
		const auto &fragment = record.fragments[i];
		const uint64_t key = (static_cast<uint64_t>(fragment.ionKind) << 32U) |
			fragment.ionPosition;
		envelopeMembers[key].push_back(i);
	}
	std::vector<char> retained(record.fragments.size(), 0);
	for (auto &entry : envelopeMembers)
	{
		auto &members = entry.second;
		std::stable_sort(members.begin(), members.end(),
			[&](size_t left, size_t right)
			{
				const float leftScore =
					record.fragments[left].theoreticalIntensity;
				const float rightScore =
					record.fragments[right].theoreticalIntensity;
				if (leftScore != rightScore)
					return leftScore > rightScore;
				return left < right;
			});
		const size_t keep = std::min<size_t>(envelopeTopN, members.size());
		for (size_t i = 0; i < keep; ++i)
			retained[members[i]] = 1;
	}
	std::vector<SpectraIndexFragmentPeakInput> compact;
	compact.reserve(record.fragments.size());
	for (size_t i = 0; i < record.fragments.size(); ++i)
		if (retained[i])
			compact.push_back(std::move(record.fragments[i]));
	record.fragments = std::move(compact);
}

static uint64_t mix64(uint64_t value)
{
	value ^= value >> 30;
	value *= 0xbf58476d1ce4e5b9ULL;
	value ^= value >> 27;
	value *= 0x94d049bb133111ebULL;
	return value ^ (value >> 31);
}

uint64_t checksumBytes(const void *data, uint64_t size, uint64_t seed,
					   uint32_t threads)
{
	constexpr uint64_t ChunkBytes = 1ULL << 20;
	const uint64_t chunkCount = size / ChunkBytes + (size % ChunkBytes != 0 ? 1 : 0);
	const auto *bytes = static_cast<const unsigned char *>(data);
	const uint64_t combined = parallelXor(
		static_cast<size_t>(chunkCount), threads, [&](size_t chunkIndex)
	{
		const uint64_t offset = chunkIndex * ChunkBytes;
		const uint64_t count = std::min(ChunkBytes, size - offset);
		uint64_t local = FnvOffsetBasis ^
			mix64(seed + chunkIndex);
		hashBytes(local, bytes + offset, count);
		return mix64(local ^ chunkIndex);
	});
	return mix64(combined ^ size ^ seed);
}

template <typename DataFn, typename SizeFn>
uint64_t checksumBlocks(uint32_t blockCount, uint64_t totalBytes,
	uint64_t seed, uint32_t threads, DataFn dataForBlock, SizeFn bytesForBlock)
{
	const uint64_t combined = parallelXor(
		blockCount, threads, [&](size_t block)
	{
		const uint64_t blockIndex = block;
		const uint64_t bytes = bytesForBlock(static_cast<uint32_t>(blockIndex));
		uint64_t local = FnvOffsetBasis ^ mix64(seed + blockIndex);
		if (bytes != 0)
			hashBytes(local, dataForBlock(static_cast<uint32_t>(blockIndex)), bytes);
		return mix64(local ^ blockIndex);
	});
	return mix64(combined ^ totalBytes ^ seed);
}

uint32_t fragmentBin(double mz)
{
	if (!(mz >= 0.0) || !std::isfinite(mz))
		return std::numeric_limits<uint32_t>::max();
	const double value = std::round(mz / FragmentBinWidth);
	if (value > static_cast<double>(PackedMassMax))
		return std::numeric_limits<uint32_t>::max();
	return static_cast<uint32_t>(value);
}

uint32_t lowerFragmentBin(double mz)
{
	if (!std::isfinite(mz))
		return std::numeric_limits<uint32_t>::max();
	const double value = std::floor(std::max(0.0, mz) / FragmentBinWidth);
	return value <= static_cast<double>(PackedMassMax)
		? static_cast<uint32_t>(value)
		: std::numeric_limits<uint32_t>::max();
}

uint32_t upperFragmentBin(double mz)
{
	if (!(mz >= 0.0) || !std::isfinite(mz))
		return std::numeric_limits<uint32_t>::max();
	const double value = std::ceil(mz / FragmentBinWidth);
	return value <= static_cast<double>(PackedMassMax)
		? static_cast<uint32_t>(value)
		: std::numeric_limits<uint32_t>::max();
}

uint32_t packFragmentPosting(uint32_t massBin, uint8_t localRecordId)
{
	return (massBin << 8) | static_cast<uint32_t>(localRecordId);
}

uint32_t retentionBin(double retentionMinutes)
{
	if (!std::isfinite(retentionMinutes))
		return std::numeric_limits<uint32_t>::max();
	const double value = std::floor(
		std::max(0.0, retentionMinutes) / RtBinWidthMinutes);
	return value <= static_cast<double>(std::numeric_limits<uint32_t>::max())
		? static_cast<uint32_t>(value)
		: std::numeric_limits<uint32_t>::max();
}

template <typename T>
bool sectionFits(uint64_t offset, uint64_t count, uint64_t fileSize)
{
	uint64_t bytes = 0;
	uint64_t end = 0;
	return checkedMultiply(count, sizeof(T), bytes) &&
		checkedAdd(offset, bytes, end) && end <= fileSize;
}

static std::string systemError(const std::string &prefix)
{
	return prefix + ": " +
		std::error_code(errno, std::generic_category()).message();
}


struct SpectraIndex::Header
{
	char magic[8];
	uint32_t version;
	uint32_t endian;
	uint32_t headerSize;
	uint32_t recordSize;
	uint32_t precursorPeakSize;
	uint32_t fragmentPeakSize;
	uint32_t productPostingSize;
	uint32_t rtBinSize;
	uint32_t recordsPerBlock;
	uint32_t blockCount;
	uint32_t sipIsotopeMassNumber;
	int32_t label;
	char sipAtom;
	char reservedChars[7];
	uint64_t fileSize;
	uint64_t recordCount;
	uint64_t precursorCount;
	uint64_t fragmentCount;
	uint64_t rtBinCount;
	uint64_t productPostingCount;
	uint64_t stringBytes;
	uint64_t recordOffset;
	uint64_t precursorOffset;
	uint64_t fragmentOffset;
	uint64_t blockRtBinOffset;
	uint64_t blockProductOffset;
	uint64_t rtBinOffset;
	uint64_t productPostingOffset;
	uint64_t stringOffset;
	double fragmentBinWidth;
	double rtBinWidthMinutes;
	double targetSipAbundancePct;
	double probabilityCutoff;
	double generationPpmTolerance;
	uint64_t minimumMatchedEnvelopes;
	uint64_t payloadChecksum;
	char chemistryProfileId[96];
	char recordKind[24];
	uint32_t envelopeTopN;
	uint32_t reservedEnvelope;
	uint64_t experimentalValueCount;
	uint64_t experimentalPresenceOffset;
	uint64_t experimentalValueOffset;
	uint64_t reserved[4];
};

static_assert(std::is_trivially_copyable<SpectraIndexRecord>::value, "SFI record must be POD");
static_assert(std::is_trivially_copyable<SpectraIndexPrecursorPeak>::value, "SFI precursor must be POD");
static_assert(std::is_trivially_copyable<SpectraIndexFragmentPeak>::value, "SFI fragment must be POD");
static_assert(sizeof(SpectraIndexRecord) == 72, "SFI record layout changed");
static_assert(sizeof(SpectraIndexFragmentPeak) == 8, "SFI fragment layout changed");
static_assert(sizeof(SpectraIndexFragmentPosting) == 4, "SFI packed posting layout changed");
static_assert(sizeof(SpectraIndexRtBin) == 8, "SFI RT-bin layout changed");

SpectraIndex::SpectraIndex() = default;

SpectraIndex::~SpectraIndex()
{
	close();
}

SpectraIndex::SpectraIndex(SpectraIndex &&other) noexcept
{
	moveFrom(other);
}

SpectraIndex &SpectraIndex::operator=(SpectraIndex &&other) noexcept
{
	if (this != &other)
	{
		close();
		moveFrom(other);
	}
	return *this;
}

void SpectraIndex::moveFrom(SpectraIndex &other) noexcept
{
	mapping_ = other.mapping_;
	mappingSize_ = other.mappingSize_;
	mappingFd_ = other.mappingFd_;
	records_ = other.records_;
	precursors_ = other.precursors_;
	fragments_ = other.fragments_;
	experimentalPresenceBits_ = other.experimentalPresenceBits_;
	experimentalValues_ = other.experimentalValues_;
	blockRtBinOffsets_ = other.blockRtBinOffsets_;
	blockProductOffsets_ = other.blockProductOffsets_;
	rtBins_ = other.rtBins_;
	productPostings_ = other.productPostings_;
	strings_ = other.strings_;
	recordCount_ = other.recordCount_;
	precursorCount_ = other.precursorCount_;
	fragmentCount_ = other.fragmentCount_;
	experimentalValueCount_ = other.experimentalValueCount_;
	rtBinCount_ = other.rtBinCount_;
	productPostingCount_ = other.productPostingCount_;
	stringBytes_ = other.stringBytes_;
	precursorBlockCount_ = other.precursorBlockCount_;
	metadata_ = std::move(other.metadata_);
	other.mapping_ = nullptr;
	other.mappingSize_ = 0;
	other.mappingFd_ = -1;
	other.records_ = nullptr;
	other.precursors_ = nullptr;
	other.fragments_ = nullptr;
	other.experimentalPresenceBits_ = nullptr;
	other.experimentalValues_ = nullptr;
	other.blockRtBinOffsets_ = nullptr;
	other.blockProductOffsets_ = nullptr;
	other.rtBins_ = nullptr;
	other.productPostings_ = nullptr;
	other.strings_ = nullptr;
	other.recordCount_ = 0;
	other.precursorCount_ = 0;
	other.fragmentCount_ = 0;
	other.experimentalValueCount_ = 0;
	other.rtBinCount_ = 0;
	other.productPostingCount_ = 0;
	other.stringBytes_ = 0;
	other.precursorBlockCount_ = 0;
}

void SpectraIndex::close()
{
	if (mapping_ != nullptr && mappingSize_ != 0)
		::munmap(mapping_, mappingSize_);
	if (mappingFd_ >= 0)
		closeFileDescriptor(mappingFd_);
	mapping_ = nullptr;
	mappingSize_ = 0;
	mappingFd_ = -1;
	records_ = nullptr;
	precursors_ = nullptr;
	fragments_ = nullptr;
	experimentalPresenceBits_ = nullptr;
	experimentalValues_ = nullptr;
	blockRtBinOffsets_ = nullptr;
	blockProductOffsets_ = nullptr;
	rtBins_ = nullptr;
	productPostings_ = nullptr;
	strings_ = nullptr;
	recordCount_ = 0;
	precursorCount_ = 0;
	fragmentCount_ = 0;
	experimentalValueCount_ = 0;
	rtBinCount_ = 0;
	productPostingCount_ = 0;
	stringBytes_ = 0;
	precursorBlockCount_ = 0;
	metadata_ = SpectraIndexMetadata();
}

bool SpectraIndex::write(const std::string &path,
						 const SpectraIndexMetadata &metadata,
						 std::vector<SpectraIndexRecordInput> &input,
						 std::string &error,
						 int threads,
						 SpectraIndexBuildStats *stats,
						 bool inputAlreadyCompacted)
{
	using Clock = std::chrono::steady_clock;
	const auto totalStart = Clock::now();
	SpectraIndexBuildStats localStats;
	auto elapsed = [](Clock::time_point start)
	{
		return std::chrono::duration<double>(Clock::now() - start).count();
	};
	error.clear();
	if (metadata.chemistryProfileId.empty() ||
		metadata.chemistryProfileId.size() >= sizeof(Header{}.chemistryProfileId) ||
		metadata.recordKind.empty() ||
		metadata.recordKind.size() >= sizeof(Header{}.recordKind) ||
		!std::isfinite(metadata.targetSipAbundancePct) ||
		(metadata.targetSipAbundancePct !=
			 SpectraIndexMetadata::MixedSipAbundancePct &&
		 (metadata.targetSipAbundancePct < 0.0 ||
		  metadata.targetSipAbundancePct > 100.0)) ||
		metadata.sipAtom == '\0' || metadata.envelopeTopN == 0)
	{
		error = "invalid SIP spectra-index metadata";
		return false;
	}
	const uint32_t requestedThreads = static_cast<uint32_t>(std::max(1, threads));
	const auto compactValidateStart = Clock::now();

	struct OrderedInput
	{
		size_t inputIndex = 0;
		double topMz = 0.0;
	};
	std::vector<OrderedInput> candidates(input.size());
	std::vector<char> valid(input.size(), 0);
	std::atomic<bool> invalidInput{false};
	parallelFor(input.size(), requestedThreads, [&](size_t i)
	{
		if (!inputAlreadyCompacted)
			retainCompactEnvelopePeaks(input[i], metadata.envelopeTopN);
		if (!std::isfinite(input[i].retentionMinutes) ||
			retentionBin(input[i].retentionMinutes) == std::numeric_limits<uint32_t>::max() ||
			!std::isfinite(input[i].sipAbundancePct) ||
			input[i].sipAbundancePct < 0.0 ||
			input[i].sipAbundancePct > 100.0)
		{
			invalidInput.store(true, std::memory_order_relaxed);
			return;
		}
		double topMz = 0.0;
		double topIntensity = -1.0;
		for (const auto &peak : input[i].precursors)
		{
			if (std::isfinite(peak.mz) && peak.mz > 0.0 &&
				std::isfinite(peak.intensity) && peak.intensity > topIntensity)
			{
				topIntensity = peak.intensity;
				topMz = peak.mz;
			}
		}
		if (input[i].charge > 0 && topMz > 0.0 &&
			!input[i].peptide.empty() && !input[i].proteins.empty())
		{
			candidates[i] = {i, topMz};
			valid[i] = 1;
		}
	});
	localStats.compactValidateSeconds = elapsed(compactValidateStart);
	if (invalidInput.load(std::memory_order_relaxed))
	{
		error = "invalid retention time or per-record SIP abundance";
		return false;
	}
	const auto orderStart = Clock::now();
	std::vector<OrderedInput> order;
	order.reserve(input.size());
	for (size_t i = 0; i < candidates.size(); ++i)
		if (valid[i])
			order.push_back(candidates[i]);
	std::vector<OrderedInput>().swap(candidates);
	std::vector<char>().swap(valid);
	parallelSort(order, requestedThreads,
		[&](const OrderedInput &a, const OrderedInput &b)
		{
			const auto &left = input[a.inputIndex];
			const auto &right = input[b.inputIndex];
			if (a.topMz != b.topMz)
				return a.topMz < b.topMz;
			if (left.charge != right.charge)
				return left.charge < right.charge;
			return a.inputIndex < b.inputIndex;
		});
	localStats.orderSeconds = elapsed(orderStart);
	if (order.size() > std::numeric_limits<uint32_t>::max())
	{
		error = "SIP spectra index exceeds uint32 record capacity";
		return false;
	}

	const auto reserveStart = Clock::now();
	size_t precursorReserve = 0;
	size_t fragmentReserve = 0;
	size_t experimentalReserve = 0;
	size_t psmStringReserve = 0;
	auto addReserve = [](size_t value, size_t &total)
	{
		if (value > std::numeric_limits<size_t>::max() - total)
			return false;
		total += value;
		return true;
	};
	auto storesExperimentalIntensity = [](float value)
	{
		uint32_t bits = 0;
		std::memcpy(&bits, &value, sizeof(bits));
		return bits != 0;
	};
	for (const OrderedInput &entry : order)
	{
		const auto &source = input[entry.inputIndex];
		if (source.precursors.size() > std::numeric_limits<uint16_t>::max() ||
			source.fragments.size() > std::numeric_limits<uint16_t>::max() ||
			source.psmId.size() > std::numeric_limits<uint16_t>::max() ||
			source.peptide.size() > std::numeric_limits<uint16_t>::max() ||
			source.proteins.size() > std::numeric_limits<uint16_t>::max() ||
			source.charge < std::numeric_limits<int16_t>::min() ||
			source.charge > std::numeric_limits<int16_t>::max())
		{
			error = "one SIP spectra record exceeds compact v6 field capacity";
			return false;
		}
		if (!addReserve(source.precursors.size(), precursorReserve) ||
			!addReserve(source.fragments.size(), fragmentReserve) ||
			!addReserve(source.psmId.size(), psmStringReserve))
		{
			error = "SIP spectra index input size overflows size_t";
			return false;
		}
		for (const auto &peak : source.fragments)
			if (storesExperimentalIntensity(peak.experimentalIntensity) &&
				!addReserve(1, experimentalReserve))
			{
				error = "SIP spectra experimental intensity count overflows size_t";
				return false;
			}
	}
	if (precursorReserve > std::numeric_limits<uint32_t>::max() ||
		fragmentReserve > std::numeric_limits<uint32_t>::max() ||
		experimentalReserve > std::numeric_limits<uint32_t>::max() ||
		psmStringReserve > std::numeric_limits<uint32_t>::max())
	{
		error = "SIP spectra index exceeds compact v6 section capacity";
		return false;
	}

	// PSM IDs are generally unique, while peptide and protein strings repeat for
	// every SIP abundance. Keep the IDs contiguous and intern only the repeated
	// search metadata so MVH peptide reads do not fault the large ID payload.
	size_t stringReserve = psmStringReserve;
	std::unordered_map<std::string, uint32_t> internedStrings;
	internedStrings.reserve(std::min<size_t>(order.size(), 262144));
	auto internString = [&](const std::string &value) -> bool
	{
		if (internedStrings.find(value) != internedStrings.end())
			return true;
		if (value.size() > std::numeric_limits<uint16_t>::max() ||
			value.size() > std::numeric_limits<size_t>::max() - stringReserve ||
			stringReserve + value.size() > std::numeric_limits<uint32_t>::max())
			return false;
		const uint32_t offset = static_cast<uint32_t>(stringReserve);
		internedStrings.emplace(value, offset);
		stringReserve += value.size();
		return true;
	};
	for (const OrderedInput &entry : order)
	{
		const auto &source = input[entry.inputIndex];
		if (!internString(source.peptide) || !internString(source.proteins))
		{
			error = "SIP spectra string dictionary exceeds compact v6 capacity";
			return false;
		}
	}

	RawArray<SpectraIndexRecord> records(order.size());
	RawArray<SpectraIndexPrecursorPeak> precursors(precursorReserve);
	RawArray<SpectraIndexFragmentPeak> fragments(fragmentReserve);
	RawArray<uint64_t> experimentalPresence(
		(fragmentReserve + 63U) / 64U);
	RawArray<float> experimentalValues(experimentalReserve);
	RawArray<char> strings(stringReserve);
	if (experimentalPresence.size() != 0)
		std::memset(experimentalPresence.data(), 0,
			experimentalPresence.size() * sizeof(experimentalPresence[0]));
	for (const auto &entry : internedStrings)
		if (!entry.first.empty())
			std::memcpy(strings.data() + entry.second,
				entry.first.data(), entry.first.size());
	size_t nextPrecursor = 0;
	size_t nextFragment = 0;
	size_t nextExperimental = 0;
	size_t nextPsmString = 0;
	for (size_t ordinal = 0; ordinal < order.size(); ++ordinal)
	{
		const auto &source = input[order[ordinal].inputIndex];
		SpectraIndexRecord record{};
		record.precursorOffset = static_cast<uint32_t>(nextPrecursor);
		record.fragmentOffset = static_cast<uint32_t>(nextFragment);
		record.experimentalOffset = static_cast<uint32_t>(nextExperimental);
		record.precursorCount = static_cast<uint16_t>(source.precursors.size());
		record.fragmentCount = static_cast<uint16_t>(source.fragments.size());
		record.psmIdOffset = static_cast<uint32_t>(nextPsmString);
		record.psmIdSize = static_cast<uint16_t>(source.psmId.size());
		nextPsmString += source.psmId.size();
		record.peptideOffset = internedStrings.find(source.peptide)->second;
		record.peptideSize = static_cast<uint16_t>(source.peptide.size());
		record.proteinsOffset = internedStrings.find(source.proteins)->second;
		record.proteinsSize = static_cast<uint16_t>(source.proteins.size());
		for (size_t index = 0; index < source.fragments.size(); ++index)
		{
			if (!storesExperimentalIntensity(
					source.fragments[index].experimentalIntensity))
				continue;
			const size_t fragmentIndex = nextFragment + index;
			experimentalPresence[fragmentIndex >> 6U] |=
				uint64_t{1} << (fragmentIndex & 63U);
			++nextExperimental;
		}
		nextPrecursor += source.precursors.size();
		nextFragment += source.fragments.size();
		records[ordinal] = record;
	}
	if (nextExperimental != experimentalReserve ||
		nextPsmString != psmStringReserve)
	{
		error = "SIP spectra compact v6 layout accounting mismatch";
		return false;
	}
	localStats.reserveSeconds = elapsed(reserveStart);
	const auto flattenStart = Clock::now();
	std::atomic<bool> flattenFailed{false};
	parallelFor(order.size(), requestedThreads, [&](size_t ordinal)
	{
		if (flattenFailed.load(std::memory_order_relaxed))
			return;
		auto &source = input[order[ordinal].inputIndex];
		SpectraIndexRecord record = records[ordinal];
			record.retentionMinutes = source.retentionMinutes;
			record.sipAbundancePct = source.sipAbundancePct;
			record.charge = static_cast<int16_t>(source.charge);
			record.generationOrdinal = static_cast<uint32_t>(order[ordinal].inputIndex);
			double topPrecursorIntensity = 0.0;
			for (size_t index = 0; index < source.precursors.size(); ++index)
			{
				const auto &peak = source.precursors[index];
				precursors[static_cast<size_t>(record.precursorOffset) + index] = peak;
				if (std::isfinite(peak.intensity) && peak.intensity > 0.0)
					record.sumPrecursorIntensity += peak.intensity;
				if (std::isfinite(peak.mz) && peak.mz > 0.0 &&
					std::isfinite(peak.intensity) &&
					peak.intensity > topPrecursorIntensity)
				{
					record.topPrecursorMz = peak.mz;
					topPrecursorIntensity = peak.intensity;
				}
			}
			size_t experimentalIndex = record.experimentalOffset;
			for (size_t index = 0; index < source.fragments.size(); ++index)
			{
				const auto &peak = source.fragments[index];
				const uint32_t mzBin = fragmentBin(peak.mz);
				if (mzBin == std::numeric_limits<uint32_t>::max() ||
					peak.ionPosition > std::numeric_limits<uint8_t>::max() ||
					!std::isfinite(peak.theoreticalIntensity) ||
					std::signbit(peak.theoreticalIntensity) ||
					!std::isfinite(peak.experimentalIntensity) ||
					(peak.ionKind != static_cast<uint8_t>('b') &&
					 peak.ionKind != static_cast<uint8_t>('B') &&
					 peak.ionKind != static_cast<uint8_t>('y') &&
					 peak.ionKind != static_cast<uint8_t>('Y')))
				{
					flattenFailed.store(true, std::memory_order_relaxed);
					return;
				}
				SpectraIndexFragmentPeak fragment{};
				fragment.packedMzPosition = mzBin |
					(static_cast<uint32_t>(peak.ionPosition) << 24U);
				std::memcpy(&fragment.theoreticalBits,
					&peak.theoreticalIntensity, sizeof(fragment.theoreticalBits));
				if (peak.ionKind == static_cast<uint8_t>('y') ||
					peak.ionKind == static_cast<uint8_t>('Y'))
					fragment.theoreticalBits |= SpectraIndexFragmentPeak::KindMask;
				fragments[static_cast<size_t>(record.fragmentOffset) + index] =
					fragment;
				if (storesExperimentalIntensity(peak.experimentalIntensity))
					experimentalValues[experimentalIndex++] =
						peak.experimentalIntensity;
			}
		if (record.psmIdSize != 0)
			std::memcpy(strings.data() + record.psmIdOffset,
				source.psmId.data(), record.psmIdSize);
		records[ordinal] = record;
		source = SpectraIndexRecordInput{};
	});
	if (flattenFailed.load(std::memory_order_relaxed))
	{
		error = "fragment cannot be represented by the compact SFI layout";
		return false;
	}
	localStats.flattenSeconds = elapsed(flattenStart);

	const uint32_t blockCount = static_cast<uint32_t>(
		(order.size() + RecordsPerBlock - 1) / RecordsPerBlock);
	const uint32_t threadsUsed = blockCount == 0 ? 1U : static_cast<uint32_t>(
		std::min<uint64_t>(static_cast<uint64_t>(requestedThreads), blockCount));
	const auto productIndexStart = Clock::now();
	struct BlockProductData
	{
		std::vector<SpectraIndexRtBin> rtBins;
		std::vector<uint32_t> postings;
	};
	std::vector<BlockProductData> productByBlock(blockCount);
	std::atomic<uint32_t> nextBlock{0};
	std::atomic<bool> productBuildFailed{false};
	auto buildBlocks = [&]()
	{
		while (true)
		{
			const uint32_t block = nextBlock.fetch_add(1, std::memory_order_relaxed);
			if (block >= blockCount)
				break;
			const uint32_t begin = block * RecordsPerBlock;
			const uint32_t end = std::min<uint32_t>(begin + RecordsPerBlock,
				static_cast<uint32_t>(order.size()));
			size_t capacity = 0;
			for (uint32_t recordId = begin; recordId < end; ++recordId)
				capacity += records[recordId].fragmentCount;
			if (capacity > std::numeric_limits<uint32_t>::max())
			{
				productBuildFailed.store(true, std::memory_order_relaxed);
				continue;
			}
			std::vector<uint64_t> entries;
			entries.reserve(capacity);
			for (uint32_t recordId = begin; recordId < end; ++recordId)
			{
				const auto &record = records[recordId];
				const uint32_t rtBin = retentionBin(record.retentionMinutes);
				for (uint64_t i = record.fragmentOffset;
					 i < record.fragmentOffset + record.fragmentCount; ++i)
				{
					const uint32_t posting = packFragmentPosting(
						fragments[static_cast<size_t>(i)].mzBin(),
						static_cast<uint8_t>(recordId - begin));
					entries.push_back((static_cast<uint64_t>(rtBin) << 32) | posting);
				}
			}
			std::sort(entries.begin(), entries.end());
			entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
			auto &output = productByBlock[block];
			output.postings.reserve(entries.size());
			uint32_t currentRtBin = std::numeric_limits<uint32_t>::max();
			for (uint64_t entry : entries)
			{
				const uint32_t rtBin = static_cast<uint32_t>(entry >> 32);
				if (rtBin != currentRtBin)
				{
					output.rtBins.push_back({rtBin,
						static_cast<uint32_t>(output.postings.size())});
					currentRtBin = rtBin;
				}
				output.postings.push_back(static_cast<uint32_t>(entry));
			}
		}
	};
	std::vector<std::thread> indexWorkers;
	indexWorkers.reserve(threadsUsed);
	for (uint32_t worker = 0; worker < threadsUsed; ++worker)
		indexWorkers.emplace_back(buildBlocks);
	for (std::thread &worker : indexWorkers)
		worker.join();
	if (productBuildFailed.load(std::memory_order_relaxed))
	{
		error = "one SFI precursor block exceeds uint32 RT-aware posting capacity";
		return false;
	}
	std::vector<uint64_t> blockRtBinOffsets(static_cast<size_t>(blockCount) + 1, 0);
	std::vector<uint64_t> blockProductOffsets(static_cast<size_t>(blockCount) + 1, 0);
	for (uint32_t block = 0; block < blockCount; ++block)
	{
		const uint64_t rtBinCount = productByBlock[block].rtBins.size();
		const uint64_t postingCount = productByBlock[block].postings.size();
		if (rtBinCount > std::numeric_limits<uint64_t>::max() - blockRtBinOffsets[block] ||
			postingCount > std::numeric_limits<uint64_t>::max() - blockProductOffsets[block])
		{
			error = "SIP spectra RT-aware product index size overflows uint64";
			return false;
		}
		blockRtBinOffsets[block + 1] = blockRtBinOffsets[block] + rtBinCount;
		blockProductOffsets[block + 1] = blockProductOffsets[block] + postingCount;
	}
	localStats.productIndexSeconds = elapsed(productIndexStart);

	Header header{};
	std::memcpy(header.magic, Magic.data(), Magic.size());
	header.version = Version;
	header.endian = EndianMarker;
	header.headerSize = sizeof(Header);
	header.recordSize = sizeof(SpectraIndexRecord);
	header.precursorPeakSize = sizeof(SpectraIndexPrecursorPeak);
	header.fragmentPeakSize = sizeof(SpectraIndexFragmentPeak);
	header.productPostingSize = sizeof(SpectraIndexFragmentPosting);
	header.rtBinSize = sizeof(SpectraIndexRtBin);
	header.recordsPerBlock = RecordsPerBlock;
	header.blockCount = blockCount;
	header.sipIsotopeMassNumber = static_cast<uint32_t>(metadata.sipIsotopeMassNumber);
	header.label = metadata.label;
	header.sipAtom = metadata.sipAtom;
	header.recordCount = records.size();
	header.precursorCount = precursors.size();
	header.fragmentCount = fragments.size();
	header.experimentalValueCount = experimentalValues.size();
	header.rtBinCount = blockRtBinOffsets.back();
	header.productPostingCount = blockProductOffsets.back();
	header.stringBytes = strings.size();
	header.fragmentBinWidth = FragmentBinWidth;
	header.rtBinWidthMinutes = RtBinWidthMinutes;
	header.targetSipAbundancePct = metadata.targetSipAbundancePct;
	header.probabilityCutoff = metadata.probabilityCutoff;
	header.generationPpmTolerance = metadata.generationPpmTolerance;
	header.minimumMatchedEnvelopes = metadata.minimumMatchedEnvelopes;
	header.envelopeTopN = metadata.envelopeTopN;
	std::memcpy(header.chemistryProfileId, metadata.chemistryProfileId.data(),
		metadata.chemistryProfileId.size());
	std::memcpy(header.recordKind, metadata.recordKind.data(), metadata.recordKind.size());

	const auto layoutStart = Clock::now();
	uint64_t end = 0;
	auto place = [&](uint64_t count, uint64_t elementSize, uint64_t &offset) -> bool
	{
		offset = alignUp(end, Alignment);
		uint64_t bytes = 0;
		return checkedMultiply(count, elementSize, bytes) && checkedAdd(offset, bytes, end);
	};
	end = sizeof(Header);
	if (!place(records.size(), sizeof(SpectraIndexRecord), header.recordOffset) ||
		!place(precursors.size(), sizeof(SpectraIndexPrecursorPeak), header.precursorOffset) ||
		!place(fragments.size(), sizeof(SpectraIndexFragmentPeak), header.fragmentOffset) ||
		!place(blockRtBinOffsets.size(), sizeof(uint64_t), header.blockRtBinOffset) ||
		!place(blockProductOffsets.size(), sizeof(uint64_t), header.blockProductOffset) ||
		!place(header.rtBinCount, sizeof(SpectraIndexRtBin), header.rtBinOffset) ||
		!place(header.productPostingCount, sizeof(SpectraIndexFragmentPosting),
			header.productPostingOffset) ||
		!place(experimentalPresence.size(), sizeof(uint64_t),
			header.experimentalPresenceOffset) ||
		!place(experimentalValues.size(), sizeof(float),
			header.experimentalValueOffset) ||
		!place(strings.size(), sizeof(char), header.stringOffset))
	{
		error = "SIP spectra index layout overflows uint64";
		return false;
	}
	header.fileSize = end;
	uint64_t checksum = 0x6a09e667f3bcc909ULL;
	checksum = mix64(checksum ^ checksumBytes(records.data(),
		records.size() * sizeof(records[0]), 0x243f6a8885a308d3ULL,
		threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(precursors.data(),
		precursors.size() * sizeof(precursors[0]), 0x13198a2e03707344ULL,
		threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(fragments.data(),
		fragments.size() * sizeof(fragments[0]), 0xa4093822299f31d0ULL,
		threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(experimentalPresence.data(),
		experimentalPresence.size() * sizeof(experimentalPresence[0]),
		0x3bd39e10cb0ef593ULL, threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(experimentalValues.data(),
		experimentalValues.size() * sizeof(experimentalValues[0]),
		0xc0acf169b5f18a8cULL, threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(blockRtBinOffsets.data(),
		blockRtBinOffsets.size() * sizeof(uint64_t), 0x299f31d0082efa98ULL,
		threadsUsed));
	checksum = mix64(checksum ^ checksumBytes(blockProductOffsets.data(),
		blockProductOffsets.size() * sizeof(uint64_t), 0xec4e6c891be6f7a1ULL,
		threadsUsed));
	checksum = mix64(checksum ^ checksumBlocks(blockCount,
		header.rtBinCount * sizeof(SpectraIndexRtBin), 0x082efa98ec4e6c89ULL,
		threadsUsed,
		[&](uint32_t block) { return productByBlock[block].rtBins.data(); },
		[&](uint32_t block) {
			return static_cast<uint64_t>(productByBlock[block].rtBins.size()) *
				sizeof(SpectraIndexRtBin);
		}));
	checksum = mix64(checksum ^ checksumBlocks(blockCount,
		header.productPostingCount * sizeof(SpectraIndexFragmentPosting),
		0x452821e638d01377ULL, threadsUsed,
		[&](uint32_t block) { return productByBlock[block].postings.data(); },
		[&](uint32_t block) {
			return static_cast<uint64_t>(productByBlock[block].postings.size()) *
				sizeof(SpectraIndexFragmentPosting);
		}));
	checksum = mix64(checksum ^ checksumBytes(strings.data(), strings.size(),
		0xbe5466cf34e90c6cULL, threadsUsed));
	header.payloadChecksum = checksum;
	localStats.layoutChecksumSeconds = elapsed(layoutStart);

	const std::filesystem::path output(path);
	std::error_code ec;
	if (output.has_parent_path())
	{
		std::filesystem::create_directories(output.parent_path(), ec);
		if (ec)
		{
			error = "cannot create SIP spectra index directory: " + ec.message();
			return false;
		}
	}
	const std::filesystem::path temporary = output.string() + ".tmp." + std::to_string(getpid());
	const auto writeStart = Clock::now();
	const int outputFd = ::open(temporary.c_str(),
		O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0600);
	if (outputFd < 0)
	{
		error = systemError("cannot create SIP spectra index " + temporary.string());
		return false;
	}
	if (header.fileSize > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
		::ftruncate(outputFd, static_cast<off_t>(header.fileSize)) != 0)
	{
		error = systemError("cannot size SIP spectra index " + temporary.string());
		closeFileDescriptor(outputFd);
		std::filesystem::remove(temporary, ec);
		return false;
	}
	struct WriteTask
	{
		uint64_t offset = 0;
		const char *data = nullptr;
		size_t bytes = 0;
	};
	std::vector<WriteTask> writeTasks;
	writeTasks.reserve(static_cast<size_t>(blockCount) * 2 + 1024);
	constexpr size_t WriteChunkBytes = 16ULL << 20;
	auto addWriteTasks = [&](uint64_t offset, const void *data, uint64_t bytes)
	{
		const char *cursor = static_cast<const char *>(data);
		while (bytes != 0)
		{
			const size_t count = static_cast<size_t>(
				std::min<uint64_t>(bytes, WriteChunkBytes));
			writeTasks.push_back({offset, cursor, count});
			offset += count;
			cursor += count;
			bytes -= count;
		}
	};
	addWriteTasks(0, &header, sizeof(header));
	addWriteTasks(header.recordOffset, records.data(),
		records.size() * sizeof(records[0]));
	addWriteTasks(header.precursorOffset, precursors.data(),
		precursors.size() * sizeof(precursors[0]));
	addWriteTasks(header.fragmentOffset, fragments.data(),
		fragments.size() * sizeof(fragments[0]));
	addWriteTasks(header.experimentalPresenceOffset,
		experimentalPresence.data(),
		experimentalPresence.size() * sizeof(experimentalPresence[0]));
	addWriteTasks(header.experimentalValueOffset, experimentalValues.data(),
		experimentalValues.size() * sizeof(experimentalValues[0]));
	addWriteTasks(header.blockRtBinOffset, blockRtBinOffsets.data(),
		blockRtBinOffsets.size() * sizeof(uint64_t));
	addWriteTasks(header.blockProductOffset, blockProductOffsets.data(),
		blockProductOffsets.size() * sizeof(uint64_t));
	// The RT bins and product postings remain in per-block vectors after the
	// parallel build, avoiding another multi-gigabyte flattening allocation.
	for (uint32_t block = 0; block < blockCount; ++block)
	{
		const auto &rtBins = productByBlock[block].rtBins;
		addWriteTasks(header.rtBinOffset +
				blockRtBinOffsets[block] * sizeof(SpectraIndexRtBin),
			rtBins.data(), rtBins.size() * sizeof(SpectraIndexRtBin));
	}
	for (uint32_t block = 0; block < blockCount; ++block)
	{
		const auto &blockPostings = productByBlock[block].postings;
		addWriteTasks(header.productPostingOffset +
				blockProductOffsets[block] * sizeof(SpectraIndexFragmentPosting),
			blockPostings.data(),
			blockPostings.size() * sizeof(SpectraIndexFragmentPosting));
	}
	addWriteTasks(header.stringOffset, strings.data(), strings.size());

	std::atomic<bool> writeFailed{false};
	parallelFor(writeTasks.size(), threadsUsed, [&](size_t taskIndex)
	{
		if (writeFailed.load(std::memory_order_relaxed))
			return;
		const WriteTask &task = writeTasks[taskIndex];
		size_t completed = 0;
		while (completed < task.bytes)
		{
			const ssize_t count = ::pwrite(outputFd, task.data + completed,
				task.bytes - completed,
				static_cast<off_t>(task.offset + completed));
			if (count < 0 && errno == EINTR)
				continue;
			if (count <= 0)
			{
				writeFailed.store(true, std::memory_order_relaxed);
				return;
			}
			completed += static_cast<size_t>(count);
		}
	});
	const bool closeOk = closeFileDescriptor(outputFd) == 0;
	const bool writeOk = !writeFailed.load(std::memory_order_relaxed) && closeOk;
	if (!writeOk)
	{
		error = "failed while writing SIP spectra index: " + temporary.string();
		std::filesystem::remove(temporary, ec);
		return false;
	}
	std::filesystem::rename(temporary, output, ec);
	if (ec)
	{
		error = "cannot publish SIP spectra index: " + ec.message();
		std::filesystem::remove(temporary, ec);
		return false;
	}
	localStats.writeSeconds = elapsed(writeStart);
	localStats.recordCount = header.recordCount;
	localStats.precursorCount = header.precursorCount;
	localStats.fragmentCount = header.fragmentCount;
	localStats.experimentalValueCount = header.experimentalValueCount;
	localStats.rtBinCount = header.rtBinCount;
	localStats.productPostingCount = header.productPostingCount;
	localStats.stringBytes = header.stringBytes;
	localStats.fileBytes = header.fileSize;
	localStats.blockCount = header.blockCount;
	localStats.threadsUsed = threadsUsed;
	localStats.totalSeconds = elapsed(totalStart);
	if (stats != nullptr)
		*stats = localStats;
	return true;
}

bool SpectraIndex::load(const std::string &path, std::string &error)
{
	error.clear();
	close();
	const int fd = ::open(path.c_str(), O_RDONLY);
	if (fd < 0)
	{
		error = systemError("cannot open SIP spectra index " + path);
		return false;
	}
	struct stat st{};
	if (::fstat(fd, &st) != 0 || st.st_size < static_cast<off_t>(sizeof(Header)))
	{
		error = "SIP spectra index is truncated: " + path;
		closeFileDescriptor(fd);
		return false;
	}
	const size_t size = static_cast<size_t>(st.st_size);
	void *mapping = ::mmap(nullptr, size, PROT_READ, MAP_SHARED, fd, 0);
	if (mapping == MAP_FAILED)
	{
		error = systemError("cannot mmap SIP spectra index " + path);
		closeFileDescriptor(fd);
		return false;
	}
	const auto *header = static_cast<const Header *>(mapping);
	auto reject = [&](const std::string &reason) -> bool
	{
		error = "invalid SIP spectra index " + path + ": " + reason;
		::munmap(mapping, size);
		closeFileDescriptor(fd);
		return false;
	};
	if (std::memcmp(header->magic, Magic.data(), Magic.size()) != 0)
	{
		if (std::memcmp(header->magic, "SIPSFI", 6) == 0)
		{
			return reject("unsupported SFI version " +
				std::to_string(header->version) +
				"; regenerate the spectra library with this Sipros build");
		}
		return reject("bad magic (HDF5 spectra libraries are not supported)");
	}
	if (header->version != Version || header->endian != EndianMarker ||
		header->headerSize != sizeof(Header) ||
		header->recordSize != sizeof(SpectraIndexRecord) ||
		header->precursorPeakSize != sizeof(SpectraIndexPrecursorPeak) ||
		header->fragmentPeakSize != sizeof(SpectraIndexFragmentPeak) ||
		header->productPostingSize != sizeof(SpectraIndexFragmentPosting) ||
		header->rtBinSize != sizeof(SpectraIndexRtBin) ||
		header->recordsPerBlock != RecordsPerBlock ||
		header->fragmentBinWidth != FragmentBinWidth ||
		header->rtBinWidthMinutes != RtBinWidthMinutes ||
		header->envelopeTopN == 0 ||
		!std::isfinite(header->targetSipAbundancePct) ||
		(header->targetSipAbundancePct !=
			 SpectraIndexMetadata::MixedSipAbundancePct &&
		 (header->targetSipAbundancePct < 0.0 ||
		  header->targetSipAbundancePct > 100.0)) ||
		header->fileSize != size)
		return reject("unsupported or inconsistent layout");
	const uint64_t expectedBlocks =
		(header->recordCount + RecordsPerBlock - 1) / RecordsPerBlock;
	const uint64_t experimentalPresenceWordCount =
		header->fragmentCount / 64U +
		(header->fragmentCount % 64U != 0 ? 1U : 0U);
	if (header->blockCount != expectedBlocks ||
		!sectionFits<SpectraIndexRecord>(header->recordOffset, header->recordCount, size) ||
		!sectionFits<SpectraIndexPrecursorPeak>(header->precursorOffset, header->precursorCount, size) ||
		!sectionFits<SpectraIndexFragmentPeak>(header->fragmentOffset, header->fragmentCount, size) ||
		!sectionFits<uint64_t>(header->blockRtBinOffset, header->blockCount + 1ULL, size) ||
		!sectionFits<uint64_t>(header->blockProductOffset, header->blockCount + 1ULL, size) ||
		!sectionFits<SpectraIndexRtBin>(header->rtBinOffset, header->rtBinCount, size) ||
		!sectionFits<SpectraIndexFragmentPosting>(header->productPostingOffset,
			header->productPostingCount, size) ||
		!sectionFits<uint64_t>(header->experimentalPresenceOffset,
			experimentalPresenceWordCount, size) ||
		!sectionFits<float>(header->experimentalValueOffset,
			header->experimentalValueCount, size) ||
		!sectionFits<char>(header->stringOffset, header->stringBytes, size))
		return reject("section outside file bounds");

	const auto *base = static_cast<const unsigned char *>(mapping);
	const auto *records = reinterpret_cast<const SpectraIndexRecord *>(base + header->recordOffset);
	const auto *precursors = reinterpret_cast<const SpectraIndexPrecursorPeak *>(base + header->precursorOffset);
	const auto *fragments = reinterpret_cast<const SpectraIndexFragmentPeak *>(base + header->fragmentOffset);
	const auto *experimentalPresence = reinterpret_cast<const uint64_t *>(
		base + header->experimentalPresenceOffset);
	const auto *experimentalValues = reinterpret_cast<const float *>(
		base + header->experimentalValueOffset);
	const auto *blockRtBins = reinterpret_cast<const uint64_t *>(base + header->blockRtBinOffset);
	const auto *blockProducts = reinterpret_cast<const uint64_t *>(base + header->blockProductOffset);
	const auto *rtBins = reinterpret_cast<const SpectraIndexRtBin *>(base + header->rtBinOffset);
	const auto *productPostings = reinterpret_cast<const SpectraIndexFragmentPosting *>(
		base + header->productPostingOffset);
	const char *strings = reinterpret_cast<const char *>(base + header->stringOffset);
	if (blockRtBins[header->blockCount] != header->rtBinCount ||
		blockProducts[header->blockCount] != header->productPostingCount)
		return reject("bad block sentinels");
	for (uint32_t block = 0; block < header->blockCount; ++block)
	{
		if (blockRtBins[block] > blockRtBins[block + 1] ||
			blockProducts[block] > blockProducts[block + 1])
			return reject("non-monotonic RT/product block offsets");
		uint32_t previousRtBin = 0;
		bool firstRtBin = true;
		const uint64_t blockPostingCount = blockProducts[block + 1] - blockProducts[block];
		for (uint64_t i = blockRtBins[block]; i < blockRtBins[block + 1]; ++i)
		{
			if ((!firstRtBin && rtBins[i].rtBin <= previousRtBin) ||
				rtBins[i].postingOffset >= blockPostingCount ||
				(i == blockRtBins[block] && rtBins[i].postingOffset != 0))
				return reject("invalid RT-bin directory");
			previousRtBin = rtBins[i].rtBin;
			firstRtBin = false;
		}
	}
	double previousMz = -std::numeric_limits<double>::infinity();
	uint32_t previousExperimentalOffset = 0;
	for (uint64_t i = 0; i < header->recordCount; ++i)
	{
		const auto &record = records[i];
		if (!(record.topPrecursorMz > 0.0) || !std::isfinite(record.topPrecursorMz) ||
			record.topPrecursorMz < previousMz || record.charge <= 0 ||
			!std::isfinite(record.retentionMinutes) ||
			!std::isfinite(record.sipAbundancePct) ||
			record.sipAbundancePct < 0.0 || record.sipAbundancePct > 100.0 ||
			record.precursorOffset > header->precursorCount ||
			record.precursorCount > header->precursorCount - record.precursorOffset ||
			record.fragmentOffset > header->fragmentCount ||
			record.fragmentCount > header->fragmentCount - record.fragmentOffset ||
			record.experimentalOffset > header->experimentalValueCount ||
			(i != 0 && record.experimentalOffset < previousExperimentalOffset) ||
			record.psmIdOffset > header->stringBytes ||
			record.psmIdSize > header->stringBytes - record.psmIdOffset ||
			record.peptideOffset > header->stringBytes ||
			record.peptideSize > header->stringBytes - record.peptideOffset ||
			record.proteinsOffset > header->stringBytes ||
			record.proteinsSize > header->stringBytes - record.proteinsOffset)
			return reject("invalid record data");
		previousMz = record.topPrecursorMz;
		previousExperimentalOffset = record.experimentalOffset;
	}
	uint64_t presentExperimentalValues = 0;
	for (uint64_t i = 0; i < experimentalPresenceWordCount; ++i)
	{
		uint64_t bits = experimentalPresence[i];
		while (bits != 0)
		{
			bits &= bits - 1U;
			++presentExperimentalValues;
		}
	}
	if (header->fragmentCount % 64U != 0 &&
		(experimentalPresence[experimentalPresenceWordCount - 1U] >>
		 (header->fragmentCount % 64U)) != 0)
		return reject("nonzero sparse-intensity padding bits");
	if (presentExperimentalValues != header->experimentalValueCount)
		return reject("sparse experimental-intensity count mismatch");
	const uint32_t checksumThreads = static_cast<uint32_t>(
		std::max(1, omp_get_max_threads()));
	uint64_t checksum = 0x6a09e667f3bcc909ULL;
	checksum = mix64(checksum ^ checksumBytes(records,
		header->recordCount * sizeof(*records), 0x243f6a8885a308d3ULL,
		checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(precursors,
		header->precursorCount * sizeof(*precursors), 0x13198a2e03707344ULL,
		checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(fragments,
		header->fragmentCount * sizeof(*fragments), 0xa4093822299f31d0ULL,
		checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(experimentalPresence,
		experimentalPresenceWordCount * sizeof(*experimentalPresence),
		0x3bd39e10cb0ef593ULL, checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(experimentalValues,
		header->experimentalValueCount * sizeof(*experimentalValues),
		0xc0acf169b5f18a8cULL, checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(blockRtBins,
		(header->blockCount + 1ULL) * sizeof(uint64_t), 0x299f31d0082efa98ULL,
		checksumThreads));
	checksum = mix64(checksum ^ checksumBytes(blockProducts,
		(header->blockCount + 1ULL) * sizeof(uint64_t), 0xec4e6c891be6f7a1ULL,
		checksumThreads));
	checksum = mix64(checksum ^ checksumBlocks(header->blockCount,
		header->rtBinCount * sizeof(*rtBins), 0x082efa98ec4e6c89ULL,
		checksumThreads,
		[&](uint32_t block) { return rtBins + blockRtBins[block]; },
		[&](uint32_t block) {
			return (blockRtBins[block + 1] - blockRtBins[block]) * sizeof(*rtBins);
		}));
	checksum = mix64(checksum ^ checksumBlocks(header->blockCount,
		header->productPostingCount * sizeof(*productPostings),
		0x452821e638d01377ULL, checksumThreads,
		[&](uint32_t block) { return productPostings + blockProducts[block]; },
		[&](uint32_t block) {
			return (blockProducts[block + 1] - blockProducts[block]) *
				sizeof(*productPostings);
		}));
	checksum = mix64(checksum ^ checksumBytes(strings, header->stringBytes,
		0xbe5466cf34e90c6cULL, checksumThreads));
	if (checksum != header->payloadChecksum)
		return reject("payload checksum mismatch");

	// Payload validation deliberately reads every byte. Drop that fully touched
	// view before search so clean checksum-only pages do not remain resident and
	// overlap the scoring caches. The replacement view refers to the same open,
	// read-only file; its pages are faulted back only when the search needs them.
	const Header validatedHeader = *header;
	if (::munmap(mapping, size) != 0)
	{
		error = systemError("cannot release validated SIP spectra index " + path);
		closeFileDescriptor(fd);
		return false;
	}
	mapping = ::mmap(nullptr, size, PROT_READ, MAP_SHARED, fd, 0);
	if (mapping == MAP_FAILED)
	{
		error = systemError("cannot remap validated SIP spectra index " + path);
		closeFileDescriptor(fd);
		return false;
	}
	header = static_cast<const Header *>(mapping);
	if (std::memcmp(header, &validatedHeader, sizeof(Header)) != 0)
		return reject("header changed during validation");
	base = static_cast<const unsigned char *>(mapping);
	records = reinterpret_cast<const SpectraIndexRecord *>(
		base + header->recordOffset);
	precursors = reinterpret_cast<const SpectraIndexPrecursorPeak *>(
		base + header->precursorOffset);
	fragments = reinterpret_cast<const SpectraIndexFragmentPeak *>(
		base + header->fragmentOffset);
	experimentalPresence = reinterpret_cast<const uint64_t *>(
		base + header->experimentalPresenceOffset);
	experimentalValues = reinterpret_cast<const float *>(
		base + header->experimentalValueOffset);
	blockRtBins = reinterpret_cast<const uint64_t *>(
		base + header->blockRtBinOffset);
	blockProducts = reinterpret_cast<const uint64_t *>(
		base + header->blockProductOffset);
	rtBins = reinterpret_cast<const SpectraIndexRtBin *>(
		base + header->rtBinOffset);
	productPostings = reinterpret_cast<const SpectraIndexFragmentPosting *>(
		base + header->productPostingOffset);
	strings = reinterpret_cast<const char *>(base + header->stringOffset);

	mapping_ = mapping;
	mappingSize_ = size;
	mappingFd_ = fd;
	records_ = records;
	precursors_ = precursors;
	fragments_ = fragments;
	experimentalPresenceBits_ = experimentalPresence;
	experimentalValues_ = experimentalValues;
	blockRtBinOffsets_ = blockRtBins;
	blockProductOffsets_ = blockProducts;
	rtBins_ = rtBins;
	productPostings_ = productPostings;
	strings_ = strings;
	recordCount_ = header->recordCount;
	precursorCount_ = header->precursorCount;
	fragmentCount_ = header->fragmentCount;
	experimentalValueCount_ = header->experimentalValueCount;
	rtBinCount_ = header->rtBinCount;
	productPostingCount_ = header->productPostingCount;
	stringBytes_ = header->stringBytes;
	precursorBlockCount_ = header->blockCount;
	metadata_.chemistryProfileId = header->chemistryProfileId;
	metadata_.recordKind = header->recordKind;
	metadata_.targetSipAbundancePct = header->targetSipAbundancePct;
	metadata_.sipAtom = header->sipAtom;
	metadata_.sipIsotopeMassNumber = static_cast<int32_t>(header->sipIsotopeMassNumber);
	metadata_.probabilityCutoff = header->probabilityCutoff;
	metadata_.generationPpmTolerance = header->generationPpmTolerance;
	metadata_.minimumMatchedEnvelopes = header->minimumMatchedEnvelopes;
	metadata_.envelopeTopN = header->envelopeTopN;
	metadata_.label = header->label;
	return true;
}

const SpectraIndexRecord &SpectraIndex::record(uint32_t recordId) const
{
	if (recordId >= recordCount_)
		throw std::out_of_range("SIP spectra-index record id is out of range");
	return records_[recordId];
}

std::string_view SpectraIndex::stringAt(uint64_t offset, uint32_t size) const
{
	if (offset > stringBytes_ || size > stringBytes_ - offset)
		throw std::out_of_range("SIP spectra-index string is out of range");
	return std::string_view(strings_ + offset, size);
}

std::string_view SpectraIndex::psmId(uint32_t id) const
{
	const auto &r = record(id);
	return stringAt(r.psmIdOffset, r.psmIdSize);
}

std::string_view SpectraIndex::peptide(uint32_t id) const
{
	const auto &r = record(id);
	return stringAt(r.peptideOffset, r.peptideSize);
}

std::string_view SpectraIndex::proteins(uint32_t id) const
{
	const auto &r = record(id);
	return stringAt(r.proteinsOffset, r.proteinsSize);
}

std::pair<const SpectraIndexPrecursorPeak *, const SpectraIndexPrecursorPeak *>
SpectraIndex::precursors(uint32_t id) const
{
	const auto &r = record(id);
	return {precursors_ + r.precursorOffset,
		precursors_ + r.precursorOffset + r.precursorCount};
}

std::pair<const SpectraIndexFragmentPeak *, const SpectraIndexFragmentPeak *>
SpectraIndex::fragments(uint32_t id) const
{
	const auto &r = record(id);
	return {fragments_ + r.fragmentOffset,
		fragments_ + r.fragmentOffset + r.fragmentCount};
}

SpectraIndexExperimentalCursor SpectraIndex::experimentalIntensities(
	uint32_t id) const
{
	const auto &r = record(id);
	return SpectraIndexExperimentalCursor(experimentalPresenceBits_,
		r.fragmentOffset, experimentalValues_ + r.experimentalOffset);
}

std::pair<uint32_t, uint32_t> SpectraIndex::precursorMzRange(
	double lowerMz, double upperMz) const
{
	if (recordCount_ == 0 || upperMz < lowerMz)
		return {0, 0};
	const auto *first = std::lower_bound(records_, records_ + recordCount_, lowerMz,
		[](const SpectraIndexRecord &record, double value)
		{ return record.topPrecursorMz < value; });
	const auto *last = std::upper_bound(first, records_ + recordCount_, upperMz,
		[](double value, const SpectraIndexRecord &record)
		{ return value < record.topPrecursorMz; });
	return {static_cast<uint32_t>(first - records_),
		static_cast<uint32_t>(last - records_)};
}

uint32_t SpectraIndex::precursorBlockForRecord(uint32_t recordId) const
{
	if (recordId >= recordCount_)
		throw std::out_of_range("SIP spectra-index record id is out of range");
	return recordId / RecordsPerBlock;
}

uint32_t SpectraIndex::precursorBlockRecordBegin(uint32_t block) const
{
	if (block >= precursorBlockCount_)
		throw std::out_of_range("SIP spectra-index precursor block is out of range");
	return block * RecordsPerBlock;
}

uint32_t SpectraIndex::postingRecordId(
	uint32_t block, const SpectraIndexFragmentPosting &posting) const
{
	const uint32_t id = precursorBlockRecordBegin(block) + posting.localRecordId();
	if (id >= recordCount_)
		throw std::out_of_range("SIP spectra-index posting record is out of range");
	return id;
}

std::pair<const SpectraIndexFragmentPosting *,
			  const SpectraIndexFragmentPosting *>
SpectraIndex::fragmentRange(uint32_t block, const SpectraIndexRtBin &rtBin,
							double lowerMz, double upperMz) const
{
	if (block >= precursorBlockCount_ || lowerMz > upperMz || upperMz < 0.0)
		return {nullptr, nullptr};
	const uint32_t lowerBin = lowerFragmentBin(lowerMz);
	const uint32_t upperBin = upperFragmentBin(upperMz);
	if (lowerBin == std::numeric_limits<uint32_t>::max() ||
		upperBin == std::numeric_limits<uint32_t>::max())
		return {nullptr, nullptr};
	const auto segment = productPostings(block, rtBin);
	const auto *postingBegin = segment.first;
	const auto *postingEnd = segment.second;
	if (postingBegin == nullptr)
		return {nullptr, nullptr};
	const uint32_t lowerPacked = packFragmentPosting(lowerBin, 0);
	const uint32_t upperPacked = packFragmentPosting(upperBin, 0xffU);
	const auto *first = std::lower_bound(postingBegin, postingEnd, lowerPacked,
		[](const SpectraIndexFragmentPosting &posting, uint32_t value)
		{ return posting.packed < value; });
	const auto *last = std::upper_bound(first, postingEnd, upperPacked,
		[](uint32_t value, const SpectraIndexFragmentPosting &posting)
		{ return value < posting.packed; });
	if (first == last)
		return {nullptr, nullptr};
	return {first, last};
}

std::pair<const SpectraIndexFragmentPosting *,
			  const SpectraIndexFragmentPosting *>
SpectraIndex::productPostings(uint32_t block,
							 const SpectraIndexRtBin &rtBin) const
{
	if (block >= precursorBlockCount_)
		return {nullptr, nullptr};
	const uint64_t rtBinIndex = static_cast<uint64_t>(&rtBin - rtBins_);
	if (rtBinIndex < blockRtBinOffsets_[block] ||
		rtBinIndex >= blockRtBinOffsets_[block + 1])
		return {nullptr, nullptr};
	const uint64_t postingBase = blockProductOffsets_[block];
	const uint64_t postingBeginOffset = postingBase + rtBin.postingOffset;
	const uint64_t postingEndOffset = rtBinIndex + 1 < blockRtBinOffsets_[block + 1]
		? postingBase + rtBins_[rtBinIndex + 1].postingOffset
		: blockProductOffsets_[block + 1];
	return {productPostings_ + postingBeginOffset,
		productPostings_ + postingEndOffset};
}

std::pair<const SpectraIndexRtBin *, const SpectraIndexRtBin *>
SpectraIndex::rtBins(uint32_t block, double lowerRtMinutes,
					 double upperRtMinutes) const
{
	if (block >= precursorBlockCount_ || lowerRtMinutes > upperRtMinutes ||
		upperRtMinutes < 0.0 || !std::isfinite(lowerRtMinutes) ||
		!std::isfinite(upperRtMinutes))
		return {nullptr, nullptr};
	const uint32_t lowerBin = retentionBin(std::max(0.0, lowerRtMinutes));
	const uint32_t upperBin = retentionBin(upperRtMinutes);
	if (lowerBin == std::numeric_limits<uint32_t>::max() ||
		upperBin == std::numeric_limits<uint32_t>::max())
		return {nullptr, nullptr};
	const auto *begin = rtBins_ + blockRtBinOffsets_[block];
	const auto *end = rtBins_ + blockRtBinOffsets_[block + 1];
	const auto *first = std::lower_bound(begin, end, lowerBin,
		[](const SpectraIndexRtBin &bin, uint32_t value)
		{ return bin.rtBin < value; });
	const auto *last = std::upper_bound(first, end, upperBin,
		[](uint32_t value, const SpectraIndexRtBin &bin)
		{ return value < bin.rtBin; });
	return first == last
		? std::make_pair(static_cast<const SpectraIndexRtBin *>(nullptr),
			static_cast<const SpectraIndexRtBin *>(nullptr))
		: std::make_pair(first, last);
}

} // namespace sipros
