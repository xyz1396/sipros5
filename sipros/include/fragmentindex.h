#ifndef SIPROS_FRAGMENT_INDEX_H
#define SIPROS_FRAGMENT_INDEX_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

class Peptide;

namespace sipros
{

// The on-disk records deliberately contain no pointers or STL objects.  This
// lets a completed index be mapped read-only and shared by concurrent searches.
struct IndexedPeptideRecord
{
	double precursorMass = 0.0;
	double precursorNeutronMass = 0.0;
	uint32_t peptideOffset = 0;
	uint32_t originalOffset = 0;
	uint32_t scoringOffset = 0;
	uint32_t proteinOffset = 0;
	uint32_t generationOrdinal = 0;
	int32_t beginPosition = 0;
	uint16_t peptideSize = 0;
	uint16_t originalSize = 0;
	uint16_t scoringSize = 0;
	uint16_t proteinSize = 0;
	uint16_t peptideLength = 0;
	char identifyPrefix = '-';
	char identifySuffix = '-';
	char originalPrefix = '-';
	char originalSuffix = '-';
};

struct FragmentPosting
{
	uint8_t localPeptideId = 0;
};

static_assert(sizeof(FragmentPosting) == 1,
	"FragmentPosting must remain the compact block-local index record");

struct FragmentBin
{
	uint32_t massBin = 0;
	uint32_t postingOffset = 0;
};

static_assert(sizeof(FragmentBin) == 8,
	"FragmentBin must remain the compact sparse-directory record");

struct FragmentIndexStats
{
	uint64_t peptideCount = 0;
	uint64_t fragmentCount = 0;
	uint64_t stringBytes = 0;
	uint64_t cacheBytes = 0;
	double fastaParseSeconds = 0.0;
	double fastaParseCpuSeconds = 0.0;
	double digestCountSeconds = 0.0;
	double digestCountCpuSeconds = 0.0;
	double digestFillSeconds = 0.0;
	double digestFillCpuSeconds = 0.0;
	double enumerateSeconds = 0.0;
	double enumerateCpuSeconds = 0.0;
	double precursorSeconds = 0.0;
	double precursorCpuSeconds = 0.0;
	double fragmentSeconds = 0.0;
	double fragmentCpuSeconds = 0.0;
	double sortSeconds = 0.0;
	double sortCpuSeconds = 0.0;
	double saveSeconds = 0.0;
	double saveCpuSeconds = 0.0;
	double mapSeconds = 0.0;
	double mapCpuSeconds = 0.0;
	double loadSeconds = 0.0;
	double loadCpuSeconds = 0.0;
	double generateSeconds = 0.0;
	double generateCpuSeconds = 0.0;
	bool loadedFromCache = false;
};

class FragmentIndex
{
public:
	FragmentIndex();
	~FragmentIndex();

	FragmentIndex(const FragmentIndex &) = delete;
	FragmentIndex &operator=(const FragmentIndex &) = delete;

	// Open or build the current v5 cache and optionally persist it at cachePath.
	bool loadOrBuild(const std::string &cachePath,
					 bool forceRebuild,
					 std::string &error);

	uint64_t peptideCount() const { return peptideCount_; }
	uint64_t fragmentCount() const { return fragmentCount_; }
	double minimumNeutronMass() const { return minimumNeutronMass_; }
	double maximumNeutronMass() const { return maximumNeutronMass_; }
	double maximumPeptideMass() const { return maximumPeptideMass_; }
	uint32_t precursorBlockCount() const { return precursorBlockCount_; }
	uint64_t fragmentBinCount() const { return fragmentBinCount_; }
	static constexpr size_t peptideBlockCapacity() { return PeptidesPerBlock; }
	const FragmentIndexStats &stats() const { return stats_; }

	const IndexedPeptideRecord &peptide(uint32_t peptideId) const;
	std::string_view peptideSequence(uint32_t peptideId) const;
	std::string_view originalSequence(uint32_t peptideId) const;
	std::string_view scoringSequence(uint32_t peptideId) const;
	std::string_view proteinNames(uint32_t peptideId) const;

	// Half-open PeptideId interval.  PeptideIds are precursor-mass ordered.
	std::pair<uint32_t, uint32_t> peptideMassRange(double lower,
												double upper) const;
	uint32_t precursorBlockForPeptide(uint32_t peptideId) const;
	uint32_t precursorBlockPeptideBegin(uint32_t precursorBlock) const;
	uint32_t postingPeptideId(
		uint32_t precursorBlock, const FragmentPosting &posting) const;

	// Half-open posting interval for one precursor-mass block and one neutral
	// fragment-mass window.  The caller applies exact peptide eligibility.
	std::pair<const FragmentPosting *, const FragmentPosting *>
	fragmentRange(uint32_t precursorBlock,
				  double lowerNeutralMass,
				  double upperNeutralMass) const;

	void materializePeptide(uint32_t peptideId, Peptide &result) const;

private:
	struct CacheHeader;

	bool build(std::string &error);
	bool save(const std::string &path, uint64_t fingerprint,
			  std::string &error);
	bool load(const std::string &path, uint64_t fingerprint,
			  std::string &error);
	uint64_t computeFingerprint(std::string &error) const;
	static uint64_t computeHeaderChecksum(const CacheHeader &header);
	void releaseMapping();
	void bindOwnedStorage();

	std::string_view stringAt(uint64_t offset, uint32_t size) const;

	std::vector<IndexedPeptideRecord> ownedPeptides_;
	std::vector<char> ownedStrings_;
	std::vector<uint64_t> ownedBlockBinOffsets_;
	std::vector<uint64_t> ownedBlockPostingOffsets_;
	std::vector<FragmentBin> ownedFragmentBins_;
	std::vector<FragmentPosting> ownedFragments_;

	const IndexedPeptideRecord *peptides_ = nullptr;
	const char *strings_ = nullptr;
	const uint64_t *blockBinOffsets_ = nullptr;
	const uint64_t *blockPostingOffsets_ = nullptr;
	const FragmentBin *fragmentBins_ = nullptr;
	const FragmentPosting *fragments_ = nullptr;
	uint64_t peptideCount_ = 0;
	uint64_t stringBytes_ = 0;
	uint64_t fragmentCount_ = 0;
	uint64_t fragmentBinCount_ = 0;
	uint32_t precursorBlockCount_ = 0;
	static constexpr uint32_t PeptidesPerBlock = 256;
	static constexpr double FragmentBinWidth = 0.001;
	static_assert(PeptidesPerBlock <=
		static_cast<uint32_t>(std::numeric_limits<uint8_t>::max()) + 1,
		"block-local peptide ids must fit in one byte");
	double minimumNeutronMass_ = 0.0;
	double maximumNeutronMass_ = 0.0;
	double maximumPeptideMass_ = 0.0;

	void *mapping_ = nullptr;
	size_t mappingSize_ = 0;
	int mappingFd_ = -1;
	FragmentIndexStats stats_;
};

} // namespace sipros

#endif // SIPROS_FRAGMENT_INDEX_H
