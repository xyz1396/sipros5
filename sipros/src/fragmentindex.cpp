#include "fragmentindex.h"

#include "MVH.h"
#include "PeptideIsotopeCalculator.h"
#include "peptide.h"
#include "performancelog.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <omp.h>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace sipros
{

namespace
{

constexpr std::array<char, 8> CacheMagic{{'S', 'I', 'P', 'R', 'O', 'S', 'F', 'I'}};
constexpr uint32_t CacheVersion = 5;
constexpr uint32_t EndianMarker = 0x01020304U;
constexpr uint64_t FnvOffsetBasis = 1469598103934665603ULL;
constexpr uint64_t FnvPrime = 1099511628211ULL;
constexpr uint64_t CacheAlignment = 64;

uint64_t alignUp(uint64_t value, uint64_t alignment)
{
	return (value + alignment - 1) / alignment * alignment;
}

bool checkedMultiply(uint64_t left, uint64_t right, uint64_t &result)
{
	if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
	{
		return false;
	}
	result = left * right;
	return true;
}

bool checkedAdd(uint64_t left, uint64_t right, uint64_t &result)
{
	if (right > std::numeric_limits<uint64_t>::max() - left)
	{
		return false;
	}
	result = left + right;
	return true;
}

void hashBytes(uint64_t &hash, const void *data, size_t size)
{
	const auto *bytes = static_cast<const unsigned char *>(data);
	for (size_t i = 0; i < size; ++i)
	{
		hash ^= bytes[i];
		hash *= FnvPrime;
	}
}

template <typename T>
void hashValue(uint64_t &hash, const T &value)
{
	static_assert(std::is_trivially_copyable<T>::value,
		"cache fingerprint values must be trivially copyable");
	hashBytes(hash, &value, sizeof(value));
}

void hashString(uint64_t &hash, const std::string &value)
{
	const uint64_t size = value.size();
	hashValue(hash, size);
	hashBytes(hash, value.data(), value.size());
}

uint64_t mix64(uint64_t value)
{
	value ^= value >> 30;
	value *= 0xbf58476d1ce4e5b9ULL;
	value ^= value >> 27;
	value *= 0x94d049bb133111ebULL;
	return value ^ (value >> 31);
}

uint64_t checksumBytes(const void *data, uint64_t size, uint64_t seed)
{
	constexpr uint64_t ChunkBytes = 1ULL << 20;
	const uint64_t chunkCount = size / ChunkBytes + (size % ChunkBytes != 0 ? 1 : 0);
	const auto *bytes = static_cast<const unsigned char *>(data);
	uint64_t combined = 0;
#pragma omp parallel for schedule(static) reduction(^ : combined)
	for (int64_t chunkIndex = 0;
		 chunkIndex < static_cast<int64_t>(chunkCount);
		 ++chunkIndex)
	{
		const uint64_t offset = static_cast<uint64_t>(chunkIndex) * ChunkBytes;
		const uint64_t count = std::min(ChunkBytes, size - offset);
		uint64_t local = FnvOffsetBasis ^
			mix64(seed + static_cast<uint64_t>(chunkIndex));
		hashBytes(local, bytes + offset, static_cast<size_t>(count));
		combined ^= mix64(local ^ static_cast<uint64_t>(chunkIndex));
	}
	return mix64(combined ^ size ^ seed);
}

uint64_t payloadChecksum(const IndexedPeptideRecord *peptides,
						 uint64_t peptideBytes,
						 const char *strings,
						 uint64_t stringBytes,
						 const uint64_t *blockBinOffsets,
						 uint64_t blockBinBytes,
						 const uint64_t *blockPostingOffsets,
						 uint64_t blockPostingBytes,
						 const FragmentBin *bins,
						 uint64_t binBytes,
						 const FragmentPosting *fragments,
						 uint64_t fragmentBytes)
{
	uint64_t checksum = 0x6a09e667f3bcc909ULL;
	checksum = mix64(checksum ^ checksumBytes(
		peptides, peptideBytes, 0x243f6a8885a308d3ULL));
	checksum = mix64(checksum ^ checksumBytes(
		strings, stringBytes, 0x13198a2e03707344ULL));
	checksum = mix64(checksum ^ checksumBytes(
		blockBinOffsets, blockBinBytes, 0xa4093822299f31d0ULL));
	checksum = mix64(checksum ^ checksumBytes(
		blockPostingOffsets, blockPostingBytes, 0x299f31d0082efa98ULL));
	checksum = mix64(checksum ^ checksumBytes(
		bins, binBytes, 0xec4e6c891be6f7a1ULL));
	checksum = mix64(checksum ^ checksumBytes(
		fragments, fragmentBytes, 0x082efa98ec4e6c89ULL));
	return checksum;
}

struct FastaEntry
{
	std::string name;
	std::string sequence;
};

bool readFastaEntries(const std::string &path,
					  std::vector<FastaEntry> &entries,
					  std::string &error)
{
	std::ifstream input(path);
	if (!input)
	{
		error = "cannot open FASTA for parallel digestion: " + path;
		return false;
	}

	std::string currentName;
	std::string currentSequence;
	std::string line;
	bool sawHeader = false;
	auto flush = [&]()
	{
		if (!currentName.empty() && !currentSequence.empty())
		{
			entries.push_back({currentName, currentSequence});
		}
		currentSequence.clear();
	};
	while (std::getline(input, line))
	{
		if (!line.empty() && line.front() == '>')
		{
			flush();
			currentName = line;
			sawHeader = true;
		}
		else if (!line.empty())
		{
			if (!sawHeader)
			{
				error = "FASTA sequence occurs before the first header: " + path;
				return false;
			}
			currentSequence += line;
		}
	}
	if (!input.eof())
	{
		error = "failed while reading FASTA for parallel digestion: " + path;
		return false;
	}
	flush();
	if (entries.empty())
	{
		error = "FASTA contains no non-empty protein entries: " + path;
		return false;
	}
	return true;
}

std::string systemError(const std::string &prefix)
{
	return prefix + ": " + std::strerror(errno);
}

class CacheBuildLock
{
public:
	~CacheBuildLock()
	{
#if defined(__unix__) || defined(__APPLE__)
		if (fd_ >= 0)
		{
			(void)flock(fd_, LOCK_UN);
			close(fd_);
		}
#endif
	}

	bool lock(const std::string &cachePath, std::string &error)
	{
#if defined(__unix__) || defined(__APPLE__)
		const std::string lockPath = cachePath + ".lock";
		fd_ = open(lockPath.c_str(), O_CREAT | O_RDWR, 0666);
		if (fd_ < 0)
		{
			error = systemError(
				"cannot open fragment-index build lock " + lockPath);
			return false;
		}
		if (flock(fd_, LOCK_EX) != 0)
		{
			error = systemError(
				"cannot acquire fragment-index build lock " + lockPath);
			close(fd_);
			fd_ = -1;
			return false;
		}
		return true;
#else
		(void)cachePath;
		(void)error;
		return true;
#endif
	}

private:
	int fd_ = -1;
};

} // namespace

struct FragmentIndex::CacheHeader
{
	char magic[8];
	uint32_t version;
	uint32_t endian;
	uint32_t headerSize;
	uint32_t peptideRecordSize;
	uint32_t fragmentRecordSize;
	uint32_t fragmentBinRecordSize;
	uint32_t peptidesPerBlock;
	uint32_t precursorBlockCount;
	uint64_t fingerprint;
	uint64_t fileSize;
	uint64_t peptideCount;
	uint64_t stringBytes;
	uint64_t fragmentCount;
	uint64_t fragmentBinCount;
	uint64_t peptideOffset;
	uint64_t stringOffset;
	uint64_t blockBinOffset;
	uint64_t blockPostingOffset;
	uint64_t binOffset;
	uint64_t fragmentOffset;
	double fragmentBinWidth;
	double minimumNeutronMass;
	double maximumNeutronMass;
	double maximumPeptideMass;
	uint64_t payloadChecksum;
	uint64_t reserved[7];
};

static_assert(std::is_trivially_copyable<IndexedPeptideRecord>::value,
	"IndexedPeptideRecord must be safe for binary cache storage");
static_assert(sizeof(IndexedPeptideRecord) == 56,
	"IndexedPeptideRecord must remain the compact v5 cache record");

uint64_t FragmentIndex::computeHeaderChecksum(const CacheHeader &header)
{
	uint64_t hash = FnvOffsetBasis;
	hashBytes(hash, header.magic, sizeof(header.magic));
	hashValue(hash, header.version);
	hashValue(hash, header.endian);
	hashValue(hash, header.headerSize);
	hashValue(hash, header.peptideRecordSize);
	hashValue(hash, header.fragmentRecordSize);
	hashValue(hash, header.fragmentBinRecordSize);
	hashValue(hash, header.peptidesPerBlock);
	hashValue(hash, header.precursorBlockCount);
	hashValue(hash, header.fingerprint);
	hashValue(hash, header.fileSize);
	hashValue(hash, header.peptideCount);
	hashValue(hash, header.stringBytes);
	hashValue(hash, header.fragmentCount);
	hashValue(hash, header.fragmentBinCount);
	hashValue(hash, header.peptideOffset);
	hashValue(hash, header.stringOffset);
	hashValue(hash, header.blockBinOffset);
	hashValue(hash, header.blockPostingOffset);
	hashValue(hash, header.binOffset);
	hashValue(hash, header.fragmentOffset);
	hashValue(hash, header.fragmentBinWidth);
	hashValue(hash, header.minimumNeutronMass);
	hashValue(hash, header.maximumNeutronMass);
	hashValue(hash, header.maximumPeptideMass);
	for (uint64_t value : header.reserved)
	{
		hashValue(hash, value);
	}
	return mix64(hash);
}

FragmentIndex::FragmentIndex() = default;

FragmentIndex::~FragmentIndex()
{
	releaseMapping();
}

void FragmentIndex::releaseMapping()
{
#if defined(__unix__) || defined(__APPLE__)
	if (mapping_ != nullptr && mappingSize_ > 0)
	{
		munmap(mapping_, mappingSize_);
	}
	if (mappingFd_ >= 0)
	{
		close(mappingFd_);
	}
#endif
	mapping_ = nullptr;
	mappingSize_ = 0;
	mappingFd_ = -1;
}

void FragmentIndex::bindOwnedStorage()
{
	peptides_ = ownedPeptides_.empty() ? nullptr : ownedPeptides_.data();
	strings_ = ownedStrings_.empty() ? nullptr : ownedStrings_.data();
	blockBinOffsets_ = ownedBlockBinOffsets_.empty()
		? nullptr
		: ownedBlockBinOffsets_.data();
	blockPostingOffsets_ = ownedBlockPostingOffsets_.empty()
		? nullptr
		: ownedBlockPostingOffsets_.data();
	fragmentBins_ = ownedFragmentBins_.empty()
		? nullptr
		: ownedFragmentBins_.data();
	fragments_ = ownedFragments_.empty() ? nullptr : ownedFragments_.data();
	peptideCount_ = ownedPeptides_.size();
	stringBytes_ = ownedStrings_.size();
	fragmentCount_ = ownedFragments_.size();
	fragmentBinCount_ = ownedFragmentBins_.size();
}

uint64_t FragmentIndex::computeFingerprint(std::string &error) const
{
	error.clear();
	uint64_t hash = FnvOffsetBasis;
	hashValue(hash, CacheVersion);

	const std::string fastaPath = ProNovoConfig::getFASTAfilename();
	std::ifstream fasta(fastaPath, std::ios::binary);
	if (!fasta)
	{
		error = "cannot open FASTA while computing index fingerprint: " + fastaPath;
		return 0;
	}
	std::array<char, 1 << 20> buffer{};
	while (fasta)
	{
		fasta.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
		const std::streamsize count = fasta.gcount();
		if (count > 0)
		{
			hashBytes(hash, buffer.data(), static_cast<size_t>(count));
		}
	}
	if (!fasta.eof())
	{
		error = "failed while reading FASTA for index fingerprint: " + fastaPath;
		return 0;
	}

	hashString(hash, ProNovoConfig::getSearchType());
	hashString(hash, ProNovoConfig::getChemistryProfileId());
	hashString(hash, ProNovoConfig::getCleavageAfterResidues());
	hashString(hash, ProNovoConfig::getCleavageBeforeResidues());
	hashValue(hash, ProNovoConfig::getMinPeptideLength());
	hashValue(hash, ProNovoConfig::getMaxPeptideLength());
	hashValue(hash, ProNovoConfig::getMaxMissedCleavages());
	hashValue(hash, ProNovoConfig::getTestStartRemoval());
	hashValue(hash, ProNovoConfig::getMaxPTMcount());
	hashValue(hash, ProNovoConfig::getTerminusMassN());
	hashValue(hash, ProNovoConfig::getTerminusMassC());

	const auto neutralLosses = ProNovoConfig::getNeutralLossList();
	for (const auto &entry : neutralLosses)
	{
		hashString(hash, entry.first);
		hashString(hash, entry.second);
	}

	std::map<std::string, std::string> ptms;
	ProNovoConfig::getPTMinfo(ptms);
	for (const auto &entry : ptms)
	{
		hashString(hash, entry.first);
		hashString(hash, entry.second);
	}

	hashValue(hash, static_cast<uint64_t>(ProNovoConfig::vsSingleResidueNames.size()));
	for (size_t i = 0; i < ProNovoConfig::vsSingleResidueNames.size(); ++i)
	{
		hashString(hash, ProNovoConfig::vsSingleResidueNames[i]);
		const double mass = i < ProNovoConfig::vdSingleResidueMasses.size()
			? ProNovoConfig::vdSingleResidueMasses[i]
			: 0.0;
		hashValue(hash, mass);
	}

	const unsigned long fragmentMask = MVH::fragmentTypes.to_ulong();
	hashValue(hash, fragmentMask);
	hashValue(hash, PeptidesPerBlock);
	hashValue(hash, FragmentBinWidth);
	return hash;
}

bool FragmentIndex::loadOrBuild(const std::string &cachePath,
								bool forceRebuild,
								std::string &error)
{
	stats_ = FragmentIndexStats{};
	std::string fingerprintError;
	const uint64_t fingerprint = computeFingerprint(fingerprintError);
	if (!fingerprintError.empty())
	{
		error = fingerprintError;
		return false;
	}

	if (!cachePath.empty() && !forceRebuild && fs::exists(cachePath))
	{
		const PerformanceTimer timer;
		std::string loadError;
		if (load(cachePath, fingerprint, loadError))
		{
			const PerformanceTiming timing = timer.elapsed();
			stats_.loadSeconds = timing.wallSeconds;
			stats_.loadCpuSeconds = timing.cpuSeconds;
			stats_.loadedFromCache = true;
			stats_.peptideCount = peptideCount_;
			stats_.fragmentCount = fragmentCount_;
			stats_.stringBytes = stringBytes_;
			stats_.cacheBytes = mappingSize_;
			return true;
		}
		std::cerr << "  Rebuilding invalid v5 peptide cache: "
				  << loadError << std::endl;
	}

	CacheBuildLock buildLock;
	if (!cachePath.empty())
	{
		const fs::path cacheFile(cachePath);
		if (cacheFile.has_parent_path())
		{
			std::error_code directoryError;
			fs::create_directories(cacheFile.parent_path(), directoryError);
			if (directoryError)
			{
				error = "cannot create fragment-index cache directory: " +
					directoryError.message();
				return false;
			}
		}
		if (!buildLock.lock(cachePath, error))
		{
			return false;
		}
		// A sibling process may have published the cache while this process
		// waited. Recheck under the lock before doing any expensive work.
		if (!forceRebuild && fs::exists(cachePath))
		{
			const PerformanceTimer timer;
			std::string loadError;
			if (load(cachePath, fingerprint, loadError))
			{
				const PerformanceTiming timing = timer.elapsed();
				stats_.loadSeconds = timing.wallSeconds;
				stats_.loadCpuSeconds = timing.cpuSeconds;
				stats_.loadedFromCache = true;
				stats_.peptideCount = peptideCount_;
				stats_.fragmentCount = fragmentCount_;
				stats_.stringBytes = stringBytes_;
				stats_.cacheBytes = mappingSize_;
				return true;
			}
		}
	}

	const PerformanceTimer generateTimer;
	if (!build(error))
	{
		return false;
	}
	if (!cachePath.empty())
	{
		const PerformanceTimer saveTimer;
		if (!save(cachePath, fingerprint, error))
		{
			return false;
		}
		const PerformanceTiming saveTiming = saveTimer.elapsed();
		stats_.saveSeconds = saveTiming.wallSeconds;
		stats_.saveCpuSeconds = saveTiming.cpuSeconds;

		// Search from the newly published mmap rather than retaining a private
		// multi-gigabyte owned copy for every scan in this process.
		const FragmentIndexStats buildStats = stats_;
		const PerformanceTimer mapTimer;
		if (!load(cachePath, fingerprint, error))
		{
			return false;
		}
		const PerformanceTiming mapTiming = mapTimer.elapsed();
		stats_ = buildStats;
		stats_.loadedFromCache = false;
		stats_.cacheBytes = mappingSize_;
		stats_.mapSeconds = mapTiming.wallSeconds;
		stats_.mapCpuSeconds = mapTiming.cpuSeconds;
	}
	const PerformanceTiming generateTiming = generateTimer.elapsed();
	stats_.generateSeconds = generateTiming.wallSeconds;
	stats_.generateCpuSeconds = generateTiming.cpuSeconds;
	return true;
}

bool FragmentIndex::build(std::string &error)
{
	error.clear();
	releaseMapping();
	// Discard any storage from an earlier build/load before constructing a new
	// owned index.
	std::vector<IndexedPeptideRecord>().swap(ownedPeptides_);
	std::vector<char>().swap(ownedStrings_);
	std::vector<uint64_t>().swap(ownedBlockBinOffsets_);
	std::vector<uint64_t>().swap(ownedBlockPostingOffsets_);
	std::vector<FragmentBin>().swap(ownedFragmentBins_);
	std::vector<FragmentPosting>().swap(ownedFragments_);
	peptides_ = nullptr;
	strings_ = nullptr;
	blockBinOffsets_ = nullptr;
	blockPostingOffsets_ = nullptr;
	fragmentBins_ = nullptr;
	fragments_ = nullptr;
	peptideCount_ = 0;
	stringBytes_ = 0;
	fragmentCount_ = 0;
	fragmentBinCount_ = 0;
	precursorBlockCount_ = 0;

	if (ProNovoConfig::getSearchType() != "Regular")
	{
		error = "fragment index currently supports only the Regular search profile";
		return false;
	}

	const PerformanceTimer enumerateTimer;
	std::vector<FastaEntry> fastaEntries;
	const PerformanceTimer parseTimer;
	if (!readFastaEntries(
			ProNovoConfig::getFASTAfilename(), fastaEntries, error))
	{
		return false;
	}
	const PerformanceTiming parseTiming = parseTimer.elapsed();
	stats_.fastaParseSeconds = parseTiming.wallSeconds;
	stats_.fastaParseCpuSeconds = parseTiming.cpuSeconds;

	struct ProteinLayout
	{
		uint64_t peptideCount = 0;
		uint64_t stringBytes = 0;
		uint64_t peptideOffset = 0;
		uint64_t stringOffset = 0;
	};
	std::vector<ProteinLayout> layouts(fastaEntries.size());
	int digestFailed = 0;
	std::string digestError;
	const PerformanceTimer countTimer;
#pragma omp parallel for schedule(dynamic, 4) reduction(| : digestFailed)
	for (int64_t entryIndex = 0;
		 entryIndex < static_cast<int64_t>(fastaEntries.size());
		 ++entryIndex)
	{
		try
		{
			const FastaEntry &entry = fastaEntries[static_cast<size_t>(entryIndex)];
			ProteinDatabase database;
			ProteinLayout layout;
			bool firstPeptide = true;
			if (database.setProteinEntry(entry.name, entry.sequence))
			{
				Peptide current;
				while (database.getNextPeptide(&current))
				{
					const std::string &peptide = current.sPeptide;
					const std::string &original = current.sOriginalPeptide;
					const std::string &protein = current.sProteinName;
					const std::string scoring =
						current.neutralLossProcess(peptide);
					const std::array<size_t, 4> sizes{{
						peptide.size(),
						original == peptide ? 0 : original.size(),
						scoring == peptide ? 0 : scoring.size(),
						firstPeptide ? protein.size() : 0}};
					for (size_t size : sizes)
					{
						if (size > std::numeric_limits<uint16_t>::max() ||
							layout.stringBytes >
								std::numeric_limits<uint64_t>::max() - size)
						{
							throw std::overflow_error(
								"indexed peptide metadata size overflow");
						}
						layout.stringBytes += size;
					}
					if (layout.peptideCount ==
						std::numeric_limits<uint64_t>::max())
					{
						throw std::overflow_error(
							"indexed peptide count overflow");
					}
					++layout.peptideCount;
					firstPeptide = false;
				}
			}
			layouts[static_cast<size_t>(entryIndex)] = layout;
		}
		catch (const std::exception &ex)
		{
			digestFailed = 1;
#pragma omp critical(sipros_fragment_index_digest_error)
			{
				if (digestError.empty())
				{
					digestError = ex.what();
				}
			}
		}
	}
	const PerformanceTiming countTiming = countTimer.elapsed();
	stats_.digestCountSeconds = countTiming.wallSeconds;
	stats_.digestCountCpuSeconds = countTiming.cpuSeconds;
	if (digestFailed != 0)
	{
		error = "failed while counting parallel FASTA digest: " + digestError;
		return false;
	}

	uint64_t totalPeptides = 0;
	uint64_t totalStringBytes = 0;
	for (ProteinLayout &layout : layouts)
	{
		layout.peptideOffset = totalPeptides;
		layout.stringOffset = totalStringBytes;
		if (!checkedAdd(totalPeptides, layout.peptideCount, totalPeptides) ||
			!checkedAdd(totalStringBytes, layout.stringBytes, totalStringBytes))
		{
			error = "parallel FASTA digest layout overflows uint64";
			return false;
		}
	}
	if (totalPeptides == 0)
	{
		error = "FASTA produced no peptides for fragment indexing";
		return false;
	}
	if (totalPeptides > std::numeric_limits<uint32_t>::max() ||
		totalPeptides > std::numeric_limits<size_t>::max() ||
		totalStringBytes > std::numeric_limits<uint32_t>::max() ||
		totalStringBytes > std::numeric_limits<size_t>::max())
	{
		error = "fragment index exceeds uint32 peptide/string offset limits";
		return false;
	}
	try
	{
		ownedPeptides_.resize(static_cast<size_t>(totalPeptides));
		ownedStrings_.resize(static_cast<size_t>(totalStringBytes));
	}
	catch (const std::exception &ex)
	{
		error = std::string("cannot allocate parallel FASTA digest: ") + ex.what();
		return false;
	}

	digestFailed = 0;
	digestError.clear();
	const PerformanceTimer fillTimer;
#pragma omp parallel for schedule(dynamic, 4) reduction(| : digestFailed)
	for (int64_t entryIndex = 0;
		 entryIndex < static_cast<int64_t>(fastaEntries.size());
		 ++entryIndex)
	{
		try
		{
			const size_t proteinIndex = static_cast<size_t>(entryIndex);
			const FastaEntry &entry = fastaEntries[proteinIndex];
			const ProteinLayout &layout = layouts[proteinIndex];
			uint64_t localPeptide = 0;
			uint64_t stringCursor = layout.stringOffset;
			uint32_t proteinOffset = 0;
			uint16_t proteinSize = 0;
			ProteinDatabase database;
			if (database.setProteinEntry(entry.name, entry.sequence))
			{
				Peptide current;
				while (database.getNextPeptide(&current))
				{
					const std::string &peptide = current.sPeptide;
					const std::string &original = current.sOriginalPeptide;
					const std::string &protein = current.sProteinName;
					const uint64_t recordIndex =
						layout.peptideOffset + localPeptide;
					if (localPeptide >= layout.peptideCount ||
						recordIndex >= ownedPeptides_.size())
					{
						throw std::runtime_error(
							"parallel digest changed between count and fill passes");
					}
					IndexedPeptideRecord record;
					record.generationOrdinal =
						static_cast<uint32_t>(recordIndex);
					record.beginPosition = current.getBeginPosProtein();
					record.identifyPrefix = current.getIdentifyPrefix();
					record.identifySuffix = current.getIdentifySuffix();
					record.originalPrefix = current.getOriginalPrefix();
					record.originalSuffix = current.getOriginalSuffix();
					for (char symbol : peptide)
					{
						if (std::isalpha(
								static_cast<unsigned char>(symbol)) != 0)
						{
							++record.peptideLength;
						}
					}
					const std::string scoring = current.neutralLossProcess(
						peptide);
					auto copyString = [&](const std::string &value,
									  uint32_t &offset,
									  uint16_t &size)
					{
						if (stringCursor > std::numeric_limits<uint32_t>::max() ||
							value.size() > std::numeric_limits<uint16_t>::max())
						{
							throw std::overflow_error(
								"compact peptide string reference overflow");
						}
						offset = static_cast<uint32_t>(stringCursor);
						size = static_cast<uint16_t>(value.size());
						if (!value.empty())
						{
							std::memcpy(
								ownedStrings_.data() + stringCursor,
								value.data(), value.size());
						}
						stringCursor += value.size();
					};
					if (localPeptide == 0)
					{
						copyString(protein, proteinOffset, proteinSize);
					}
					copyString(peptide, record.peptideOffset,
						record.peptideSize);
					if (original == peptide)
					{
						record.originalOffset = record.peptideOffset;
						record.originalSize = record.peptideSize;
					}
					else
					{
						copyString(original, record.originalOffset,
							record.originalSize);
					}
					if (scoring == peptide)
					{
						record.scoringOffset = record.peptideOffset;
						record.scoringSize = record.peptideSize;
					}
					else
					{
						copyString(scoring, record.scoringOffset,
							record.scoringSize);
					}
					record.proteinOffset = proteinOffset;
					record.proteinSize = proteinSize;
					ownedPeptides_[static_cast<size_t>(recordIndex)] = record;
					++localPeptide;
				}
			}
			if (localPeptide != layout.peptideCount ||
				stringCursor != layout.stringOffset + layout.stringBytes)
			{
				throw std::runtime_error(
					"parallel digest byte layout changed between passes");
			}
		}
		catch (const std::exception &ex)
		{
			digestFailed = 1;
#pragma omp critical(sipros_fragment_index_fill_error)
			{
				if (digestError.empty())
				{
					digestError = ex.what();
				}
			}
		}
	}
	const PerformanceTiming fillTiming = fillTimer.elapsed();
	stats_.digestFillSeconds = fillTiming.wallSeconds;
	stats_.digestFillCpuSeconds = fillTiming.cpuSeconds;
	if (digestFailed != 0)
	{
		error = "failed while filling parallel FASTA digest: " + digestError;
		return false;
	}
	const PerformanceTiming enumerateTiming = enumerateTimer.elapsed();
	stats_.enumerateSeconds = enumerateTiming.wallSeconds;
	stats_.enumerateCpuSeconds = enumerateTiming.cpuSeconds;
	if (ownedPeptides_.empty())
	{
		error = "FASTA produced no peptides for fragment indexing";
		return false;
	}
	if (ownedPeptides_.size() > std::numeric_limits<uint32_t>::max())
	{
		error = "fragment index exceeds the uint32 PeptideId limit";
		return false;
	}

	const PerformanceTimer precursorTimer;
	std::exception_ptr estimateError;
#pragma omp parallel
	{
		PeptideIsotopeCalculator calculator;
#pragma omp for schedule(guided)
		for (int64_t i = 0; i < static_cast<int64_t>(ownedPeptides_.size()); ++i)
		{
			try
			{
				const IndexedPeptideRecord &record = ownedPeptides_[static_cast<size_t>(i)];
				const std::string_view decorated(
					ownedStrings_.data() + record.peptideOffset,
					record.peptideSize);
				std::string composition;
				composition.reserve(decorated.size());
				for (char symbol : decorated)
				{
					if (symbol != '[' && symbol != ']')
					{
						composition.push_back(symbol);
					}
				}
				const auto estimate = calculator.calPrecursorEstimate(composition);
				ownedPeptides_[static_cast<size_t>(i)].precursorMass = estimate.mass;
				ownedPeptides_[static_cast<size_t>(i)].precursorNeutronMass =
					estimate.neutronMass;
			}
			catch (...)
			{
#pragma omp critical(sipros_fragment_index_estimate_error)
				{
					if (!estimateError)
					{
						estimateError = std::current_exception();
					}
				}
			}
		}
	}
	if (estimateError)
	{
		try
		{
			std::rethrow_exception(estimateError);
		}
		catch (const std::exception &ex)
		{
			error = std::string("failed to estimate indexed precursor mass: ") + ex.what();
		}
		return false;
	}

	std::sort(ownedPeptides_.begin(), ownedPeptides_.end(),
		[](const IndexedPeptideRecord &left, const IndexedPeptideRecord &right)
		{
			if (left.precursorMass != right.precursorMass)
			{
				return left.precursorMass < right.precursorMass;
			}
			return left.generationOrdinal < right.generationOrdinal;
		});
	minimumNeutronMass_ = std::numeric_limits<double>::infinity();
	maximumNeutronMass_ = 0.0;
	maximumPeptideMass_ = 0.0;
	for (const IndexedPeptideRecord &record : ownedPeptides_)
	{
		minimumNeutronMass_ = std::min(minimumNeutronMass_,
			record.precursorNeutronMass);
		maximumNeutronMass_ = std::max(maximumNeutronMass_,
			record.precursorNeutronMass);
		maximumPeptideMass_ = std::max(maximumPeptideMass_,
			record.precursorMass);
	}
	const PerformanceTiming precursorTiming = precursorTimer.elapsed();
	stats_.precursorSeconds = precursorTiming.wallSeconds;
	stats_.precursorCpuSeconds = precursorTiming.cpuSeconds;

	const int enabledSeries =
		(MVH::fragmentTypes[FragmentType_B] ? 1 : 0) +
		(MVH::fragmentTypes[FragmentType_Y] ? 1 : 0);
	if (enabledSeries == 0)
	{
		error = "MVH fragment index has neither b nor y ions enabled";
		return false;
	}

	std::vector<uint64_t> peptideFragmentOffsets(ownedPeptides_.size() + 1, 0);
	for (size_t i = 0; i < ownedPeptides_.size(); ++i)
	{
		const uint64_t count = static_cast<uint64_t>(ownedPeptides_[i].peptideLength) *
			static_cast<uint64_t>(enabledSeries);
		if (peptideFragmentOffsets[i] >
			std::numeric_limits<uint64_t>::max() - count)
		{
			error = "fragment posting count overflow";
			return false;
		}
		peptideFragmentOffsets[i + 1] = peptideFragmentOffsets[i] + count;
	}
	const uint64_t totalFragments = peptideFragmentOffsets.back();
	if (totalFragments > static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
	{
		error = "fragment index is too large for this address space";
		return false;
	}

	precursorBlockCount_ = static_cast<uint32_t>(
		(ownedPeptides_.size() + PeptidesPerBlock - 1) / PeptidesPerBlock);
	ownedBlockPostingOffsets_.resize(
		static_cast<size_t>(precursorBlockCount_) + 1);
	for (uint32_t block = 0; block <= precursorBlockCount_; ++block)
	{
		const size_t peptideBegin = std::min(
			static_cast<size_t>(block) * PeptidesPerBlock,
			ownedPeptides_.size());
		ownedBlockPostingOffsets_[block] = peptideFragmentOffsets[peptideBegin];
	}

	const PerformanceTimer fragmentTimer;
	// The temporary 64-bit key keeps generation parallel and sorting simple.
	// Its high bits are the sparse fragment-mass bin and its low byte is the
	// peptide id relative to a bounded precursor block.
	std::vector<uint64_t> encodedFragments;
	try
	{
		encodedFragments.resize(static_cast<size_t>(totalFragments));
	}
	catch (const std::exception &ex)
	{
		error = std::string("cannot allocate fragment build keys: ") + ex.what();
		return false;
	}
	std::exception_ptr fragmentError;
#pragma omp parallel
	{
		std::vector<double> chargedIons;
		std::vector<double> forward;
		std::vector<double> reverse;
		std::vector<char> residues;
#pragma omp for schedule(guided)
		for (int64_t index = 0;
			 index < static_cast<int64_t>(ownedPeptides_.size());
			 ++index)
		{
			try
			{
				const uint32_t peptideId = static_cast<uint32_t>(index);
				const IndexedPeptideRecord &record = ownedPeptides_[peptideId];
				std::string scoring(ownedStrings_.data() + record.scoringOffset,
					record.scoringSize);
				if (!MVH::CalculateSequenceIons(scoring, 2,
						MVH::bUseSmartPlusThreeModel, &chargedIons,
						&forward, &reverse, &residues))
				{
					throw std::runtime_error("MVH rejected indexed scoring sequence");
				}
				uint64_t cursor = peptideFragmentOffsets[peptideId];
				const uint8_t localPeptideId = static_cast<uint8_t>(
					peptideId % PeptidesPerBlock);
				auto appendMass = [&](double mass)
				{
					if (!std::isfinite(mass) || mass <= 0.0 ||
						mass / FragmentBinWidth >
							static_cast<double>(std::numeric_limits<uint32_t>::max()))
					{
						throw std::overflow_error("fragment mass-bin overflow");
					}
					const uint32_t massBin = static_cast<uint32_t>(
						std::floor(mass / FragmentBinWidth));
					encodedFragments[static_cast<size_t>(cursor++)] =
						(static_cast<uint64_t>(massBin) << 8) | localPeptideId;
				};
				if (MVH::fragmentTypes[FragmentType_B])
				{
					for (double mass : forward)
					{
						appendMass(mass);
					}
				}
				if (MVH::fragmentTypes[FragmentType_Y])
				{
					for (double mass : reverse)
					{
						appendMass(mass);
					}
				}
				if (cursor != peptideFragmentOffsets[peptideId + 1])
				{
					throw std::runtime_error("unexpected indexed neutral-fragment count");
				}
			}
			catch (...)
			{
#pragma omp critical(sipros_fragment_index_generation_error)
				{
					if (!fragmentError)
					{
						fragmentError = std::current_exception();
					}
				}
			}
		}
	}
	if (fragmentError)
	{
		try
		{
			std::rethrow_exception(fragmentError);
		}
		catch (const std::exception &ex)
		{
			error = std::string("failed to build fragment postings: ") + ex.what();
		}
		return false;
	}
	const PerformanceTiming fragmentTiming = fragmentTimer.elapsed();
	stats_.fragmentSeconds = fragmentTiming.wallSeconds;
	stats_.fragmentCpuSeconds = fragmentTiming.cpuSeconds;

	const PerformanceTimer sortTimer;
#pragma omp parallel for schedule(dynamic, 1)
	for (int64_t block = 0; block < static_cast<int64_t>(precursorBlockCount_); ++block)
	{
		auto begin = encodedFragments.begin() + static_cast<std::ptrdiff_t>(
			ownedBlockPostingOffsets_[static_cast<size_t>(block)]);
		auto end = encodedFragments.begin() + static_cast<std::ptrdiff_t>(
			ownedBlockPostingOffsets_[static_cast<size_t>(block) + 1]);
		std::sort(begin, end);
	}

	std::vector<uint64_t> blockBinCounts(precursorBlockCount_, 0);
#pragma omp parallel for schedule(static)
	for (int64_t block = 0; block < static_cast<int64_t>(precursorBlockCount_); ++block)
	{
		const uint64_t begin = ownedBlockPostingOffsets_[static_cast<size_t>(block)];
		const uint64_t end = ownedBlockPostingOffsets_[static_cast<size_t>(block) + 1];
		uint64_t count = 0;
		uint32_t previousBin = 0;
		for (uint64_t i = begin; i < end; ++i)
		{
			const uint32_t massBin = static_cast<uint32_t>(encodedFragments[i] >> 8);
			if (i == begin || massBin != previousBin)
			{
				++count;
				previousBin = massBin;
			}
		}
		blockBinCounts[static_cast<size_t>(block)] = count;
	}
	ownedBlockBinOffsets_.assign(
		static_cast<size_t>(precursorBlockCount_) + 1, 0);
	for (uint32_t block = 0; block < precursorBlockCount_; ++block)
	{
		ownedBlockBinOffsets_[block + 1] =
			ownedBlockBinOffsets_[block] + blockBinCounts[block];
	}
	try
	{
		ownedFragmentBins_.resize(
			static_cast<size_t>(ownedBlockBinOffsets_.back()));
		ownedFragments_.resize(static_cast<size_t>(totalFragments));
	}
	catch (const std::exception &ex)
	{
		error = std::string("cannot allocate compact fragment index: ") + ex.what();
		return false;
	}
#pragma omp parallel for schedule(static)
	for (int64_t block = 0; block < static_cast<int64_t>(precursorBlockCount_); ++block)
	{
		const uint64_t postingBegin =
			ownedBlockPostingOffsets_[static_cast<size_t>(block)];
		const uint64_t postingEnd =
			ownedBlockPostingOffsets_[static_cast<size_t>(block) + 1];
		uint64_t binCursor = ownedBlockBinOffsets_[static_cast<size_t>(block)];
		uint32_t previousBin = 0;
		for (uint64_t i = postingBegin; i < postingEnd; ++i)
		{
			const uint32_t massBin = static_cast<uint32_t>(encodedFragments[i] >> 8);
			if (i == postingBegin || massBin != previousBin)
			{
				ownedFragmentBins_[static_cast<size_t>(binCursor++)] = {
					massBin, static_cast<uint32_t>(i - postingBegin)};
				previousBin = massBin;
			}
			ownedFragments_[static_cast<size_t>(i)] = {
				static_cast<uint8_t>(encodedFragments[i] & 0xffU)};
		}
	}
	const PerformanceTiming sortTiming = sortTimer.elapsed();
	stats_.sortSeconds = sortTiming.wallSeconds;
	stats_.sortCpuSeconds = sortTiming.cpuSeconds;

	std::vector<uint64_t>().swap(peptideFragmentOffsets);
	std::vector<uint64_t>().swap(blockBinCounts);
	std::vector<uint64_t>().swap(encodedFragments);
	bindOwnedStorage();
	stats_.peptideCount = peptideCount_;
	stats_.fragmentCount = fragmentCount_;
	stats_.stringBytes = stringBytes_;
	return true;
}

bool FragmentIndex::save(const std::string &path,
						 uint64_t fingerprint,
						 std::string &error)
{
	error.clear();
	if (path.empty())
	{
		return true;
	}
	const fs::path output(path);
	std::error_code ec;
	if (output.has_parent_path())
	{
		fs::create_directories(output.parent_path(), ec);
		if (ec)
		{
			error = "cannot create fragment-index cache directory: " + ec.message();
			return false;
		}
	}
	const fs::path temporary = output.string() + ".tmp." +
		std::to_string(static_cast<unsigned long long>(
#if defined(__unix__) || defined(__APPLE__)
			getpid()
#else
			0
#endif
		));

	CacheHeader header{};
	std::memcpy(header.magic, CacheMagic.data(), CacheMagic.size());
	header.version = CacheVersion;
	header.endian = EndianMarker;
	header.headerSize = sizeof(CacheHeader);
	header.peptideRecordSize = sizeof(IndexedPeptideRecord);
	header.fragmentRecordSize = sizeof(FragmentPosting);
	header.fragmentBinRecordSize = sizeof(FragmentBin);
	header.peptidesPerBlock = PeptidesPerBlock;
	header.precursorBlockCount = precursorBlockCount_;
	header.fingerprint = fingerprint;
	header.peptideCount = peptideCount_;
	header.stringBytes = stringBytes_;
	header.fragmentCount = fragmentCount_;
	header.fragmentBinCount = fragmentBinCount_;
	header.fragmentBinWidth = FragmentBinWidth;
	header.minimumNeutronMass = minimumNeutronMass_;
	header.maximumNeutronMass = maximumNeutronMass_;
	header.maximumPeptideMass = maximumPeptideMass_;
	uint64_t peptideBytes = 0;
	uint64_t blockIndexBytes = 0;
	uint64_t binBytes = 0;
	uint64_t fragmentBytes = 0;
	uint64_t sectionEnd = 0;
	if (!checkedMultiply(peptideCount_, sizeof(IndexedPeptideRecord), peptideBytes) ||
		!checkedMultiply(static_cast<uint64_t>(precursorBlockCount_) + 1,
			sizeof(uint64_t), blockIndexBytes) ||
		!checkedMultiply(fragmentBinCount_, sizeof(FragmentBin), binBytes) ||
		!checkedMultiply(fragmentCount_, sizeof(FragmentPosting), fragmentBytes))
	{
		error = "fragment-index cache layout overflows uint64";
		return false;
	}
	header.peptideOffset = alignUp(sizeof(CacheHeader), CacheAlignment);
	if (!checkedAdd(header.peptideOffset, peptideBytes, sectionEnd))
	{
		error = "fragment-index peptide section overflows uint64";
		return false;
	}
	header.stringOffset = alignUp(sectionEnd, CacheAlignment);
	if (!checkedAdd(header.stringOffset, stringBytes_, sectionEnd))
	{
		error = "fragment-index string section overflows uint64";
		return false;
	}
	header.blockBinOffset = alignUp(sectionEnd, CacheAlignment);
	if (!checkedAdd(header.blockBinOffset, blockIndexBytes, sectionEnd))
	{
		error = "fragment-index block-bin section overflows uint64";
		return false;
	}
	header.blockPostingOffset = alignUp(sectionEnd, CacheAlignment);
	if (!checkedAdd(header.blockPostingOffset, blockIndexBytes, sectionEnd))
	{
		error = "fragment-index block-posting section overflows uint64";
		return false;
	}
	header.binOffset = alignUp(sectionEnd, CacheAlignment);
	if (!checkedAdd(header.binOffset, binBytes, sectionEnd))
	{
		error = "fragment-index sparse-bin section overflows uint64";
		return false;
	}
	header.fragmentOffset = alignUp(sectionEnd, CacheAlignment);
	if (!checkedAdd(header.fragmentOffset, fragmentBytes, header.fileSize))
	{
		error = "fragment-index fragment section overflows uint64";
		return false;
	}
	const uint64_t dataChecksum = payloadChecksum(
		peptides_, peptideBytes,
		strings_, stringBytes_,
		blockBinOffsets_, blockIndexBytes,
		blockPostingOffsets_, blockIndexBytes,
		fragmentBins_, binBytes,
		fragments_, fragmentBytes);
	header.payloadChecksum = mix64(
		dataChecksum ^ computeHeaderChecksum(header));

	std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
	if (!out)
	{
		error = "cannot create fragment-index cache: " + temporary.string();
		return false;
	}
	auto writeAt = [&](uint64_t offset, const void *data, uint64_t bytes) -> bool
	{
		out.seekp(static_cast<std::streamoff>(offset));
		if (!out)
		{
			return false;
		}
		constexpr uint64_t MaxChunk = 1ULL << 30;
		const char *cursor = static_cast<const char *>(data);
		while (bytes > 0)
		{
			const uint64_t chunk = std::min(bytes, MaxChunk);
			out.write(cursor, static_cast<std::streamsize>(chunk));
			if (!out)
			{
				return false;
			}
			cursor += chunk;
			bytes -= chunk;
		}
		return true;
	};

	if (!writeAt(0, &header, sizeof(header)) ||
		!writeAt(header.peptideOffset, peptides_,
			peptideCount_ * sizeof(IndexedPeptideRecord)) ||
		!writeAt(header.stringOffset, strings_, stringBytes_) ||
		!writeAt(header.blockBinOffset, blockBinOffsets_,
			(static_cast<uint64_t>(precursorBlockCount_) + 1) * sizeof(uint64_t)) ||
		!writeAt(header.blockPostingOffset, blockPostingOffsets_,
			(static_cast<uint64_t>(precursorBlockCount_) + 1) * sizeof(uint64_t)) ||
		!writeAt(header.binOffset, fragmentBins_,
			fragmentBinCount_ * sizeof(FragmentBin)) ||
		!writeAt(header.fragmentOffset, fragments_,
			fragmentCount_ * sizeof(FragmentPosting)))
	{
		error = "failed while writing fragment-index cache: " + temporary.string();
		out.close();
		fs::remove(temporary, ec);
		return false;
	}
	out.flush();
	if (!out)
	{
		error = "failed to flush fragment-index cache: " + temporary.string();
		out.close();
		fs::remove(temporary, ec);
		return false;
	}
	out.close();
	if (!out)
	{
		error = "failed to close fragment-index cache: " + temporary.string();
		fs::remove(temporary, ec);
		return false;
	}
#if defined(__unix__) || defined(__APPLE__)
	const int temporaryFd = open(temporary.c_str(), O_RDWR);
	if (temporaryFd < 0 || fsync(temporaryFd) != 0)
	{
		const int syncError = errno;
		if (temporaryFd >= 0)
		{
			close(temporaryFd);
		}
		errno = syncError;
		error = systemError("cannot sync fragment-index cache " + temporary.string());
		fs::remove(temporary, ec);
		return false;
	}
	close(temporaryFd);
#endif
	fs::rename(temporary, output, ec);
	if (ec)
	{
		error = "cannot publish fragment-index cache: " + ec.message();
		fs::remove(temporary, ec);
		return false;
	}
#if defined(__unix__) || defined(__APPLE__)
	const fs::path parent = output.has_parent_path()
		? output.parent_path()
		: fs::path(".");
	const int parentFd = open(parent.c_str(), O_RDONLY);
	if (parentFd >= 0)
	{
		(void)fsync(parentFd);
		close(parentFd);
	}
#endif
	stats_.cacheBytes = header.fileSize;
	return true;
}

bool FragmentIndex::load(const std::string &path,
						 uint64_t fingerprint,
						 std::string &error)
{
	error.clear();
	releaseMapping();
	// After a cold build/save, clear() would retain the complete private index
	// capacity alongside the new read-only mmap.  Swap releases that storage.
	std::vector<IndexedPeptideRecord>().swap(ownedPeptides_);
	std::vector<char>().swap(ownedStrings_);
	std::vector<uint64_t>().swap(ownedBlockBinOffsets_);
	std::vector<uint64_t>().swap(ownedBlockPostingOffsets_);
	std::vector<FragmentBin>().swap(ownedFragmentBins_);
	std::vector<FragmentPosting>().swap(ownedFragments_);
	peptides_ = nullptr;
	strings_ = nullptr;
	blockBinOffsets_ = nullptr;
	blockPostingOffsets_ = nullptr;
	fragmentBins_ = nullptr;
	fragments_ = nullptr;
	peptideCount_ = 0;
	stringBytes_ = 0;
	fragmentCount_ = 0;
	fragmentBinCount_ = 0;
	precursorBlockCount_ = 0;

#if defined(__unix__) || defined(__APPLE__)
	mappingFd_ = open(path.c_str(), O_RDONLY);
	if (mappingFd_ < 0)
	{
		error = systemError("cannot open fragment-index cache " + path);
		return false;
	}
	struct stat fileStat{};
	if (fstat(mappingFd_, &fileStat) != 0)
	{
		error = systemError("cannot stat fragment-index cache " + path);
		releaseMapping();
		return false;
	}
	if (fileStat.st_size < static_cast<off_t>(sizeof(CacheHeader)))
	{
		error = "fragment-index cache is truncated: " + path;
		releaseMapping();
		return false;
	}
	mappingSize_ = static_cast<size_t>(fileStat.st_size);
	mapping_ = mmap(nullptr, mappingSize_, PROT_READ, MAP_SHARED, mappingFd_, 0);
	if (mapping_ == MAP_FAILED)
	{
		mapping_ = nullptr;
		error = systemError("cannot mmap fragment-index cache " + path);
		releaseMapping();
		return false;
	}
	const auto *header = static_cast<const CacheHeader *>(mapping_);
	if (std::memcmp(header->magic, CacheMagic.data(), CacheMagic.size()) != 0 ||
		header->version != CacheVersion ||
		header->endian != EndianMarker ||
		header->headerSize != sizeof(CacheHeader) ||
		header->peptideRecordSize != sizeof(IndexedPeptideRecord) ||
		header->fragmentRecordSize != sizeof(FragmentPosting) ||
		header->fragmentBinRecordSize != sizeof(FragmentBin) ||
		header->peptidesPerBlock != PeptidesPerBlock ||
		header->fragmentBinWidth != FragmentBinWidth)
	{
		error = "file does not use the current v5 peptide-cache schema: " + path;
		releaseMapping();
		return false;
	}
	if (header->fingerprint != fingerprint)
	{
		error = "fragment-index cache fingerprint does not match FASTA/search chemistry";
		releaseMapping();
		return false;
	}
	uint64_t peptideBytes = 0;
	uint64_t blockIndexBytes = 0;
	uint64_t binBytes = 0;
	uint64_t fragmentBytes = 0;
	if (!checkedMultiply(header->peptideCount,
			sizeof(IndexedPeptideRecord), peptideBytes) ||
		!checkedMultiply(static_cast<uint64_t>(header->precursorBlockCount) + 1,
			sizeof(uint64_t), blockIndexBytes) ||
		!checkedMultiply(header->fragmentBinCount,
			sizeof(FragmentBin), binBytes) ||
		!checkedMultiply(header->fragmentCount,
			sizeof(FragmentPosting), fragmentBytes))
	{
		error = "fragment-index cache contains overflowing section sizes: " + path;
		releaseMapping();
		return false;
	}
	auto validRange = [&](uint64_t offset, uint64_t bytes) -> bool
	{
		return offset <= header->fileSize && bytes <= header->fileSize - offset;
	};
	if (header->fileSize != mappingSize_ ||
		header->peptideCount == 0 ||
		header->peptideCount > std::numeric_limits<uint32_t>::max() ||
		header->stringBytes > std::numeric_limits<uint32_t>::max() ||
		header->precursorBlockCount == 0 ||
		header->precursorBlockCount !=
			(header->peptideCount + PeptidesPerBlock - 1) / PeptidesPerBlock ||
		header->fragmentBinCount == 0 ||
		header->fragmentBinCount > header->fragmentCount ||
		!std::isfinite(header->minimumNeutronMass) ||
		!std::isfinite(header->maximumNeutronMass) ||
		header->minimumNeutronMass <= 0.0 ||
		header->maximumNeutronMass < header->minimumNeutronMass ||
		!std::isfinite(header->maximumPeptideMass) ||
		header->maximumPeptideMass <= 0.0 ||
		header->peptideOffset < sizeof(CacheHeader) ||
		header->peptideOffset % CacheAlignment != 0 ||
		header->stringOffset % CacheAlignment != 0 ||
		header->blockBinOffset % CacheAlignment != 0 ||
		header->blockPostingOffset % CacheAlignment != 0 ||
		header->binOffset % CacheAlignment != 0 ||
		header->fragmentOffset % CacheAlignment != 0 ||
		!validRange(header->peptideOffset, peptideBytes) ||
		!validRange(header->stringOffset, header->stringBytes) ||
		!validRange(header->blockBinOffset, blockIndexBytes) ||
		!validRange(header->blockPostingOffset, blockIndexBytes) ||
		!validRange(header->binOffset, binBytes) ||
		!validRange(header->fragmentOffset, fragmentBytes) ||
		header->peptideOffset + peptideBytes > header->stringOffset ||
		header->stringOffset + header->stringBytes > header->blockBinOffset ||
		header->blockBinOffset + blockIndexBytes > header->blockPostingOffset ||
		header->blockPostingOffset + blockIndexBytes > header->binOffset ||
		header->binOffset + binBytes > header->fragmentOffset ||
		header->fragmentOffset + fragmentBytes != header->fileSize)
	{
		error = "fragment-index cache contains invalid offsets: " + path;
		releaseMapping();
		return false;
	}
	const char *base = static_cast<const char *>(mapping_);
	const auto *mappedPeptides = reinterpret_cast<const IndexedPeptideRecord *>(
		base + header->peptideOffset);
	const char *mappedStrings = base + header->stringOffset;
	const auto *mappedBlockBinOffsets = reinterpret_cast<const uint64_t *>(
		base + header->blockBinOffset);
	const auto *mappedBlockPostingOffsets = reinterpret_cast<const uint64_t *>(
		base + header->blockPostingOffset);
	const auto *mappedBins = reinterpret_cast<const FragmentBin *>(
		base + header->binOffset);
	const auto *mappedFragments = reinterpret_cast<const FragmentPosting *>(
		base + header->fragmentOffset);
	const uint64_t dataChecksum = payloadChecksum(
			mappedPeptides, peptideBytes,
			mappedStrings, header->stringBytes,
			mappedBlockBinOffsets, blockIndexBytes,
			mappedBlockPostingOffsets, blockIndexBytes,
			mappedBins, binBytes,
			mappedFragments, fragmentBytes);
	if (mix64(dataChecksum ^ computeHeaderChecksum(*header)) !=
		header->payloadChecksum)
	{
		error = "fragment-index cache payload checksum failed: " + path;
		releaseMapping();
		return false;
	}
	if (mappedBlockBinOffsets[0] != 0 ||
		mappedBlockBinOffsets[header->precursorBlockCount] !=
			header->fragmentBinCount ||
		mappedBlockPostingOffsets[0] != 0 ||
		mappedBlockPostingOffsets[header->precursorBlockCount] !=
			header->fragmentCount)
	{
		error = "fragment-index cache has invalid block endpoints: " + path;
		releaseMapping();
		return false;
	}
	for (uint32_t block = 0; block < header->precursorBlockCount; ++block)
	{
		if (mappedBlockBinOffsets[block] >= mappedBlockBinOffsets[block + 1] ||
			mappedBlockBinOffsets[block + 1] > header->fragmentBinCount ||
			mappedBlockPostingOffsets[block] >=
				mappedBlockPostingOffsets[block + 1] ||
			mappedBlockPostingOffsets[block + 1] > header->fragmentCount)
		{
			error = "fragment-index cache has non-monotonic block offsets: " + path;
			releaseMapping();
			return false;
		}
	}
	int invalidPosting = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(| : invalidPosting)
	for (int64_t block = 0;
		 block < static_cast<int64_t>(header->precursorBlockCount);
		 ++block)
	{
		const uint64_t binBegin =
			mappedBlockBinOffsets[static_cast<size_t>(block)];
		const uint64_t binEnd =
			mappedBlockBinOffsets[static_cast<size_t>(block) + 1];
		const uint64_t postingBegin =
			mappedBlockPostingOffsets[static_cast<size_t>(block)];
		const uint64_t postingEnd =
			mappedBlockPostingOffsets[static_cast<size_t>(block) + 1];
		const uint32_t blockPeptideCount = static_cast<uint32_t>(std::min<uint64_t>(
			PeptidesPerBlock,
			header->peptideCount - static_cast<uint64_t>(block) * PeptidesPerBlock));
		uint32_t previousMassBin = 0;
		for (uint64_t binIndex = binBegin; binIndex < binEnd; ++binIndex)
		{
			const FragmentBin &bin = mappedBins[binIndex];
			const uint64_t relativeEnd = binIndex + 1 < binEnd
				? mappedBins[binIndex + 1].postingOffset
				: postingEnd - postingBegin;
			if ((binIndex > binBegin && bin.massBin <= previousMassBin) ||
				(binIndex == binBegin && bin.postingOffset != 0) ||
				bin.postingOffset >= relativeEnd ||
				relativeEnd > postingEnd - postingBegin)
			{
				invalidPosting = 1;
				break;
			}
			uint8_t previousLocalId = 0;
			for (uint64_t postingIndex = postingBegin + bin.postingOffset;
				 postingIndex < postingBegin + relativeEnd; ++postingIndex)
			{
				const uint8_t localId = mappedFragments[postingIndex].localPeptideId;
				if (localId >= blockPeptideCount ||
					(postingIndex > postingBegin + bin.postingOffset &&
					 localId < previousLocalId))
				{
					invalidPosting = 1;
					break;
				}
				previousLocalId = localId;
			}
			previousMassBin = bin.massBin;
			if (invalidPosting != 0)
				break;
		}
	}
	if (invalidPosting != 0)
	{
		error = "fragment-index cache has invalid or unsorted postings: " + path;
		releaseMapping();
		return false;
	}
	double previousMass = -std::numeric_limits<double>::infinity();
	for (uint64_t peptideId = 0; peptideId < header->peptideCount; ++peptideId)
	{
		const IndexedPeptideRecord &record = mappedPeptides[peptideId];
		auto validString = [&](uint64_t offset, uint32_t size) -> bool
		{
			return offset <= header->stringBytes &&
				size <= header->stringBytes - offset;
		};
		if (!std::isfinite(record.precursorMass) ||
			record.precursorMass <= 0.0 ||
			record.precursorMass < previousMass ||
			!std::isfinite(record.precursorNeutronMass) ||
			record.precursorNeutronMass <= 0.0 ||
			record.generationOrdinal >= header->peptideCount ||
			record.peptideLength == 0 ||
			!validString(record.peptideOffset, record.peptideSize) ||
			!validString(record.originalOffset, record.originalSize) ||
			!validString(record.scoringOffset, record.scoringSize) ||
			!validString(record.proteinOffset, record.proteinSize))
		{
			error = "fragment-index cache has an invalid peptide record: " + path;
			releaseMapping();
			return false;
		}
		previousMass = record.precursorMass;
	}
	peptides_ = mappedPeptides;
	strings_ = mappedStrings;
	blockBinOffsets_ = mappedBlockBinOffsets;
	blockPostingOffsets_ = mappedBlockPostingOffsets;
	fragmentBins_ = mappedBins;
	fragments_ = mappedFragments;
	peptideCount_ = header->peptideCount;
	stringBytes_ = header->stringBytes;
	fragmentCount_ = header->fragmentCount;
	fragmentBinCount_ = header->fragmentBinCount;
	precursorBlockCount_ = header->precursorBlockCount;
	minimumNeutronMass_ = header->minimumNeutronMass;
	maximumNeutronMass_ = header->maximumNeutronMass;
	maximumPeptideMass_ = header->maximumPeptideMass;
	return true;
#else
	(void)path;
	(void)fingerprint;
	error = "memory-mapped fragment-index caches are not implemented on this platform";
	return false;
#endif
}

const IndexedPeptideRecord &FragmentIndex::peptide(uint32_t peptideId) const
{
	if (peptideId >= peptideCount_ || peptides_ == nullptr)
	{
		throw std::out_of_range("fragment-index PeptideId is out of range");
	}
	return peptides_[peptideId];
}

std::string_view FragmentIndex::stringAt(uint64_t offset, uint32_t size) const
{
	if (offset > stringBytes_ || size > stringBytes_ - offset || strings_ == nullptr)
	{
		throw std::out_of_range("fragment-index string offset is out of range");
	}
	return std::string_view(strings_ + offset, size);
}

std::string_view FragmentIndex::peptideSequence(uint32_t peptideId) const
{
	const auto &record = peptide(peptideId);
	return stringAt(record.peptideOffset, record.peptideSize);
}

std::string_view FragmentIndex::originalSequence(uint32_t peptideId) const
{
	const auto &record = peptide(peptideId);
	return stringAt(record.originalOffset, record.originalSize);
}

std::string_view FragmentIndex::scoringSequence(uint32_t peptideId) const
{
	const auto &record = peptide(peptideId);
	return stringAt(record.scoringOffset, record.scoringSize);
}

std::string_view FragmentIndex::proteinNames(uint32_t peptideId) const
{
	const auto &record = peptide(peptideId);
	return stringAt(record.proteinOffset, record.proteinSize);
}

std::pair<uint32_t, uint32_t> FragmentIndex::peptideMassRange(
	double lower, double upper) const
{
	if (peptides_ == nullptr || peptideCount_ == 0 || lower > upper)
	{
		return {0, 0};
	}
	const IndexedPeptideRecord *begin = peptides_;
	const IndexedPeptideRecord *end = peptides_ + peptideCount_;
	const auto first = std::lower_bound(begin, end, lower,
		[](const IndexedPeptideRecord &record, double value)
		{
			return record.precursorMass < value;
		});
	const auto last = std::upper_bound(first, end, upper,
		[](double value, const IndexedPeptideRecord &record)
		{
			return value < record.precursorMass;
		});
	return {
		static_cast<uint32_t>(first - begin),
		static_cast<uint32_t>(last - begin)};
}

uint32_t FragmentIndex::precursorBlockForPeptide(uint32_t peptideId) const
{
	if (peptideId >= peptideCount_ || precursorBlockCount_ == 0)
	{
		throw std::out_of_range("fragment-index peptide block is out of range");
	}
	return peptideId / PeptidesPerBlock;
}

uint32_t FragmentIndex::precursorBlockPeptideBegin(
	uint32_t precursorBlock) const
{
	if (precursorBlock >= precursorBlockCount_)
	{
		throw std::out_of_range("fragment-index precursor block is out of range");
	}
	return precursorBlock * PeptidesPerBlock;
}

uint32_t FragmentIndex::postingPeptideId(
	uint32_t precursorBlock, const FragmentPosting &posting) const
{
	if (precursorBlock >= precursorBlockCount_)
	{
		throw std::out_of_range("fragment-index posting block is out of range");
	}
	const uint32_t peptideId = precursorBlock * PeptidesPerBlock +
		posting.localPeptideId;
	if (peptideId >= peptideCount_)
	{
		throw std::out_of_range("fragment-index local peptide id is out of range");
	}
	return peptideId;
}

std::pair<const FragmentPosting *, const FragmentPosting *>
FragmentIndex::fragmentRange(uint32_t precursorBlock,
							 double lowerNeutralMass,
							 double upperNeutralMass) const
{
	if (fragments_ == nullptr || fragmentBins_ == nullptr ||
		blockBinOffsets_ == nullptr || blockPostingOffsets_ == nullptr ||
		precursorBlock >= precursorBlockCount_ ||
		lowerNeutralMass > upperNeutralMass)
	{
		return {nullptr, nullptr};
	}
	const double scale = std::max(1.0,
		std::max(std::abs(lowerNeutralMass), std::abs(upperNeutralMass)));
	const double guard = 2.0 * std::numeric_limits<float>::epsilon() * scale;
	const double guardedLower = std::max(0.0, lowerNeutralMass - guard);
	const double guardedUpper = std::max(0.0, upperNeutralMass + guard);
	const uint32_t lowerBin = static_cast<uint32_t>(std::min(
		std::floor(guardedLower / FragmentBinWidth),
		static_cast<double>(std::numeric_limits<uint32_t>::max())));
	const uint32_t upperBin = static_cast<uint32_t>(std::min(
		std::floor(guardedUpper / FragmentBinWidth),
		static_cast<double>(std::numeric_limits<uint32_t>::max())));
	const uint64_t binStart = blockBinOffsets_[precursorBlock];
	const uint64_t binStop = blockBinOffsets_[precursorBlock + 1];
	const FragmentBin *binBegin = fragmentBins_ + binStart;
	const FragmentBin *binEnd = fragmentBins_ + binStop;
	const FragmentBin *firstBin = std::lower_bound(binBegin, binEnd, lowerBin,
		[](const FragmentBin &bin, uint32_t value)
		{
			return bin.massBin < value;
		});
	const FragmentBin *lastBin = std::upper_bound(firstBin, binEnd, upperBin,
		[](uint32_t value, const FragmentBin &bin)
		{
			return value < bin.massBin;
		});
	if (firstBin == lastBin)
	{
		return {nullptr, nullptr};
	}
	const uint64_t postingBase = blockPostingOffsets_[precursorBlock];
	const uint64_t firstPosting = postingBase + firstBin->postingOffset;
	const uint64_t lastPosting = lastBin == binEnd
		? blockPostingOffsets_[precursorBlock + 1]
		: postingBase + lastBin->postingOffset;
	return {fragments_ + firstPosting, fragments_ + lastPosting};
}

void FragmentIndex::materializePeptide(uint32_t peptideId, Peptide &result) const
{
	const IndexedPeptideRecord &record = peptide(peptideId);
	result = Peptide();
	result.setPeptide(
		std::string(peptideSequence(peptideId)),
		std::string(originalSequence(peptideId)),
		std::string(proteinNames(peptideId)),
		record.beginPosition,
		record.precursorMass,
		record.identifyPrefix,
		record.identifySuffix,
		record.originalPrefix,
		record.originalSuffix);
	result.setPrecursorNeutronMass(record.precursorNeutronMass);
	result.sNeutralLossPeptide = std::string(scoringSequence(peptideId));
	result.iPeptideLength = record.peptideLength;
}

} // namespace sipros
