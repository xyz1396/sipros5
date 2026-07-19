#ifndef SIPROS_INDEXED_DATABASE_SEARCH_H
#define SIPROS_INDEXED_DATABASE_SEARCH_H

#include <cstdint>
#include <utility>
#include <vector>

class MS2Scan;

namespace sipros
{

class FragmentIndex;

struct IndexedCandidate
{
	uint32_t peptideId = 0;
	double measuredMass = 0.0;
	int charge = 0;
	uint32_t hypothesisOrdinal = 0;
};

struct IndexedSearchCounters
{
	uint64_t precursorHypotheses = 0;
	uint64_t exactMassCandidates = 0;
	uint64_t uniquePeptideChargeCandidates = 0;
	uint64_t fragmentPostingsVisited = 0;
	uint64_t fragmentGateSurvivors = 0;
	uint64_t exactMvhCalls = 0;
	uint64_t exactMvhAccepted = 0;

	IndexedSearchCounters &operator+=(const IndexedSearchCounters &other);
};

class IndexedSearchScratch
{
public:
	std::vector<IndexedCandidate> candidates;
};

// Query one already-MVH-preprocessed scan. The preliminary neutral-fragment
// count is deliberately an upper bound; every survivor is returned for exact
// charge-dependent MVH scoring in a separate, independently timed phase.
void queryIndexedScan(
	const FragmentIndex &index,
	MS2Scan &scan,
	const std::vector<std::pair<double, int>> &precursorMassCharges,
	IndexedSearchScratch &scratch,
	std::vector<IndexedCandidate> &survivors,
	IndexedSearchCounters &counters);

void scoreIndexedScanMvh(
	const FragmentIndex &index,
	MS2Scan &scan,
	const std::vector<IndexedCandidate> &survivors,
	std::vector<double> &sequenceIonMasses,
	std::vector<double> &aaForward,
	std::vector<double> &aaReverse,
	std::vector<char> &residues,
	IndexedSearchCounters &counters);

} // namespace sipros

#endif // SIPROS_INDEXED_DATABASE_SEARCH_H
