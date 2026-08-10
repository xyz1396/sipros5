#include "indexeddatabasesearch.h"

#include "MVH.h"
#include "fragmentindex.h"
#include "ms2scan.h"
#include "peptide.h"
#include "proNovoConfig.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <tuple>

namespace sipros
{

namespace
{

// Validated DDA first-gate settings. Square-root transformation does not
// change which experimental peaks belong to the top-N set.
constexpr uint16_t DdaWindowMinMatchedFragments = 4;
constexpr size_t DdaWindowTopPeaks = 200;

bool candidateIdentityLess(const IndexedCandidate &left,
	const IndexedCandidate &right)
{
	return std::tie(left.peptideId, left.charge, left.hypothesisOrdinal) <
		std::tie(right.peptideId, right.charge, right.hypothesisOrdinal);
}

bool mergeIndexedPeptideIfPresent(
	const FragmentIndex &index,
	uint32_t peptideId,
	int precursorCharge,
	std::vector<PeptideUnit *> &topPeptides)
{
	const std::string_view sequence = index.peptideSequence(peptideId);
	const std::string proteins(index.proteinNames(peptideId));
	for (PeptideUnit *top : topPeptides)
	{
		if (top->iMeasuredParentCharge != precursorCharge ||
			std::string_view(top->sIdentifiedPeptide) != sequence)
		{
			continue;
		}
		if (top->sProteinNames != proteins &&
			top->sProteinNames.find("," + proteins) == std::string::npos)
		{
			top->sProteinNames += "," + proteins;
		}
		return true;
	}
	return false;
}

struct ResidualEvidence
{
	PeptideUnit *peptide = nullptr;
	std::vector<size_t> peakIndices;
	double totalPeakWeight = 0.0;
};

ResidualEvidence collectResidualEvidence(
	PeptideUnit *peptide,
	MS2Scan &scan,
	std::vector<double> &sequenceIonMasses,
	std::vector<double> &aaForward,
	std::vector<double> &aaReverse,
	std::vector<char> &residues)
{
	ResidualEvidence evidence;
	evidence.peptide = peptide;
	if (peptide == nullptr || scan.pPeakList == nullptr ||
		scan.pPeakList->pPeaks.empty() ||
		!MVH::CalculateSequenceIons(
			peptide->sPeptideForScoring,
			peptide->iMeasuredParentCharge,
			MVH::bUseSmartPlusThreeModel,
			&sequenceIonMasses,
			&aaForward,
			&aaReverse,
			&residues))
	{
		return evidence;
	}

	const double tolerance = ProNovoConfig::getMassAccuracyFragmentIon();
	const auto &peaks = scan.pPeakList->pPeaks;
	const auto &classes = scan.pPeakList->pClasses;
	for (double ionMz : sequenceIonMasses)
	{
		auto peak = std::lower_bound(peaks.begin(), peaks.end(), ionMz - tolerance);
		size_t best = std::numeric_limits<size_t>::max();
		unsigned char bestClass = std::numeric_limits<unsigned char>::max();
		for (; peak != peaks.end() && *peak <= ionMz + tolerance; ++peak)
		{
			const size_t index = static_cast<size_t>(peak - peaks.begin());
			const unsigned char intensityClass = index < classes.size()
				? static_cast<unsigned char>(classes[index])
				: std::numeric_limits<unsigned char>::max();
			if (intensityClass < bestClass)
			{
				best = index;
				bestClass = intensityClass;
			}
		}
		if (best != std::numeric_limits<size_t>::max())
			evidence.peakIndices.push_back(best);
	}
	std::sort(evidence.peakIndices.begin(), evidence.peakIndices.end());
	evidence.peakIndices.erase(
		std::unique(evidence.peakIndices.begin(), evidence.peakIndices.end()),
		evidence.peakIndices.end());
	for (size_t index : evidence.peakIndices)
	{
		const int intensityClass = static_cast<unsigned char>(classes[index]);
		evidence.totalPeakWeight += std::max(
			1, ProNovoConfig::NumIntensityClasses + 1 - intensityClass);
	}
	return evidence;
}

void assignDdaResidualRanks(
	MS2Scan &scan,
	std::vector<double> &sequenceIonMasses,
	std::vector<double> &aaForward,
	std::vector<double> &aaReverse,
	std::vector<char> &residues)
{
	if (scan.pPeakList == nullptr || scan.vpWeightSumTopPeptides.empty())
		return;
	std::vector<ResidualEvidence> evidence;
	evidence.reserve(scan.vpWeightSumTopPeptides.size());
	for (PeptideUnit *peptide : scan.vpWeightSumTopPeptides)
	{
		evidence.push_back(collectResidualEvidence(
			peptide, scan, sequenceIonMasses, aaForward, aaReverse, residues));
	}
	std::vector<unsigned char> active(scan.pPeakList->pPeaks.size(), 1);
	for (int rank = 1; rank <= PeptideUnit::DdaResidualTopN; ++rank)
	{
		size_t bestIndex = std::numeric_limits<size_t>::max();
		double bestScore = -std::numeric_limits<double>::infinity();
		for (size_t i = 0; i < evidence.size(); ++i)
		{
			ResidualEvidence &candidate = evidence[i];
			if (candidate.peptide == nullptr ||
				candidate.peptide->iDdaResidualRank != 0 ||
				candidate.totalPeakWeight <= 0.0)
			{
				continue;
			}
			size_t activeMatches = 0;
			double activeWeight = 0.0;
			for (size_t peakIndex : candidate.peakIndices)
			{
				if (!active[peakIndex])
					continue;
				++activeMatches;
				const int intensityClass = static_cast<unsigned char>(
					scan.pPeakList->pClasses[peakIndex]);
				activeWeight += std::max(
					1, ProNovoConfig::NumIntensityClasses + 1 - intensityClass);
			}
			if (activeMatches < static_cast<size_t>(
					ProNovoConfig::MinMatchedFragments))
			{
				continue;
			}
			const double residualFraction =
				activeWeight / candidate.totalPeakWeight;
			const double score =
				candidate.peptide->vdScores[2] * residualFraction;
			if (score > bestScore ||
				(score == bestScore &&
				 bestIndex != std::numeric_limits<size_t>::max() &&
				 candidate.peptide->sIdentifiedPeptide <
					evidence[bestIndex].peptide->sIdentifiedPeptide))
			{
				bestIndex = i;
				bestScore = score;
			}
		}
		if (bestIndex == std::numeric_limits<size_t>::max())
			break;
		ResidualEvidence &winner = evidence[bestIndex];
		winner.peptide->iDdaResidualRank = rank;
		winner.peptide->dDdaResidualScore = bestScore;
		for (size_t peakIndex : winner.peakIndices)
			active[peakIndex] = 0;
	}
}

void selectDdaWindowGatePeaks(
	const MS2Scan &scan,
	std::vector<double> &gatePeakMzs)
{
	gatePeakMzs.clear();
	if (scan.pPeakList == nullptr)
		return;
	const std::vector<double> &filteredMzs = scan.pPeakList->pPeaks;
	if (filteredMzs.size() <= DdaWindowTopPeaks)
	{
		gatePeakMzs = filteredMzs;
		return;
	}

	// MVH preprocessing keeps at most 300 peaks. Recover their raw intensities
	// with a monotone merge, retain the top 200, then restore m/z order.
	std::vector<std::pair<double, double>> intensityMz;
	intensityMz.reserve(filteredMzs.size());
	size_t rawIndex = 0;
	for (double mz : filteredMzs)
	{
		while (rawIndex < scan.vdMZ.size() && scan.vdMZ[rawIndex] < mz)
			++rawIndex;
		const double intensity = rawIndex < scan.vdMZ.size() &&
			rawIndex < scan.vdIntensity.size() && scan.vdMZ[rawIndex] == mz
				? scan.vdIntensity[rawIndex]
				: 0.0;
		intensityMz.push_back({intensity, mz});
	}
	std::partial_sort(
		intensityMz.begin(),
		intensityMz.begin() + static_cast<std::ptrdiff_t>(DdaWindowTopPeaks),
		intensityMz.end(),
		[](const auto &left, const auto &right)
		{
			if (left.first != right.first)
				return left.first > right.first;
			return left.second < right.second;
		});
	gatePeakMzs.reserve(DdaWindowTopPeaks);
	for (size_t i = 0; i < DdaWindowTopPeaks; ++i)
		gatePeakMzs.push_back(intensityMz[i].second);
	std::sort(gatePeakMzs.begin(), gatePeakMzs.end());
}

// Gate one contiguous isolation-window range without materializing it, using
// singly charged b/y ions and four matches.
void gateIsolationWindowBlock(
	const FragmentIndex &index,
	uint32_t block,
	uint32_t windowBegin,
	uint32_t windowEnd,
	int precursorCharge,
	uint32_t windowHypothesisOrdinal,
	const std::vector<double> &gatePeakMzs,
	double fragmentTolerance,
	uint16_t hitThreshold,
	std::vector<IndexedCandidate> &survivors,
	IndexedSearchCounters &counters)
{
	const uint32_t peptideBegin = index.precursorBlockPeptideBegin(block);
	const uint32_t peptideEnd = static_cast<uint32_t>(std::min<uint64_t>(
		index.peptideCount(),
		static_cast<uint64_t>(peptideBegin) +
			FragmentIndex::peptideBlockCapacity()));
	windowBegin = std::max(windowBegin, peptideBegin);
	windowEnd = std::min(windowEnd, peptideEnd);
	if (windowBegin >= windowEnd)
		return;

	std::array<uint16_t, FragmentIndex::peptideBlockCapacity()> hitCounts{};
	size_t candidatesBelowThreshold = windowEnd - windowBegin;
	if (hitThreshold > 0)
	{
		for (double observedMz : gatePeakMzs)
		{
			const double neutralMass = observedMz - Proton;
			const auto postings = index.fragmentRange(
				block,
				neutralMass - fragmentTolerance,
				neutralMass + fragmentTolerance);
			if (postings.first == nullptr)
				continue;
			counters.fragmentPostingsVisited +=
				static_cast<uint64_t>(postings.second - postings.first);
			for (const FragmentPosting *posting = postings.first;
				 posting != postings.second; ++posting)
			{
				const uint32_t localId = posting->localPeptideId;
				const uint32_t peptideId = peptideBegin + localId;
				uint16_t &state = hitCounts[localId];
				if (peptideId >= windowBegin && peptideId < windowEnd &&
					state < hitThreshold)
				{
					++state;
					if (state == hitThreshold)
						--candidatesBelowThreshold;
				}
			}
			if (candidatesBelowThreshold == 0)
				break;
		}
	}

	for (uint32_t peptideId = windowBegin; peptideId < windowEnd; ++peptideId)
	{
		const size_t localId = peptideId - peptideBegin;
		if (hitCounts[localId] < hitThreshold)
			continue;
		const IndexedPeptideRecord &record = index.peptide(peptideId);
		survivors.push_back({
			peptideId,
			record.precursorMass,
			precursorCharge,
			windowHypothesisOrdinal});
		++counters.fragmentGateSurvivors;
	}
}

} // namespace

IndexedSearchCounters &IndexedSearchCounters::operator+=(
	const IndexedSearchCounters &other)
{
	isolationWindowRanges += other.isolationWindowRanges;
	isolationWindowCandidates += other.isolationWindowCandidates;
	exactMassCandidates += other.exactMassCandidates;
	uniquePeptideChargeCandidates += other.uniquePeptideChargeCandidates;
	fragmentPostingsVisited += other.fragmentPostingsVisited;
	fragmentGateSurvivors += other.fragmentGateSurvivors;
	exactMvhCalls += other.exactMvhCalls;
	exactMvhAccepted += other.exactMvhAccepted;
	return *this;
}

void queryIndexedScan(
	const FragmentIndex &index,
	MS2Scan &scan,
	IndexedSearchScratch &scratch,
	std::vector<IndexedCandidate> &survivors,
	IndexedSearchCounters &counters)
{
	survivors.clear();
	if (scan.bSkip || scan.pPeakList == nullptr ||
		scan.vIsolationWindowsMz.empty() || index.peptideCount() == 0)
	{
		return;
	}

	selectDdaWindowGatePeaks(scan, scratch.gatePeakMzs);
	if (scratch.gatePeakMzs.empty())
		return;
	const double fragmentTolerance = ProNovoConfig::getMassAccuracyFragmentIon();
	const uint16_t hitThreshold = static_cast<uint16_t>(std::min(
		std::max(ProNovoConfig::MinMatchedFragments,
			static_cast<int>(DdaWindowMinMatchedFragments)),
		static_cast<int>(std::numeric_limits<uint16_t>::max())));
	const double protonMass = ProNovoConfig::getProtonMass();

	// Fragment-first DDA+: the database candidate set is defined only by each
	// acquisition window and the configured charge range, never by an MS1 peak.
	for (int precursorCharge = 1; precursorCharge <= 4; ++precursorCharge)
	{
		for (size_t windowIndex = 0;
			 windowIndex < scan.vIsolationWindowsMz.size(); ++windowIndex)
		{
			const double centerMz = scan.vIsolationWindowsMz[windowIndex].first;
			const double widthMz = scan.vIsolationWindowsMz[windowIndex].second;
			if (!(centerMz > 0.0) || !(widthMz > 0.0) ||
				!std::isfinite(centerMz) || !std::isfinite(widthMz))
			{
				continue;
			}
			const double lower =
				(centerMz - widthMz / 2.0) * precursorCharge -
				static_cast<double>(precursorCharge) * protonMass;
			const double upper =
				(centerMz + widthMz / 2.0) * precursorCharge -
				static_cast<double>(precursorCharge) * protonMass;
			const auto range = index.peptideMassRange(lower, upper);
			++counters.isolationWindowRanges;
			const uint64_t candidateCount =
				static_cast<uint64_t>(range.second - range.first);
			counters.isolationWindowCandidates += candidateCount;
			counters.exactMassCandidates += candidateCount;
			counters.uniquePeptideChargeCandidates += candidateCount;
			if (range.first == range.second)
				continue;

			const uint32_t firstBlock =
				index.precursorBlockForPeptide(range.first);
			const uint32_t lastBlock =
				index.precursorBlockForPeptide(range.second - 1);
			for (uint32_t block = firstBlock; block <= lastBlock; ++block)
			{
				gateIsolationWindowBlock(
					index, block, range.first, range.second,
					precursorCharge, static_cast<uint32_t>(windowIndex),
					scratch.gatePeakMzs, fragmentTolerance, hitThreshold,
					survivors, counters);
			}
		}
	}

	// Multiple reaction records can overlap. Score each peptide/charge once.
	std::sort(survivors.begin(), survivors.end(), candidateIdentityLess);
	survivors.erase(
		std::unique(
			survivors.begin(), survivors.end(),
			[](const IndexedCandidate &left, const IndexedCandidate &right)
			{
				return left.peptideId == right.peptideId &&
					left.charge == right.charge;
			}),
		survivors.end());
	std::sort(survivors.begin(), survivors.end(),
		[&](const IndexedCandidate &left, const IndexedCandidate &right)
		{
			const uint32_t leftOrdinal =
				index.peptide(left.peptideId).generationOrdinal;
			const uint32_t rightOrdinal =
				index.peptide(right.peptideId).generationOrdinal;
			return std::tie(leftOrdinal, left.hypothesisOrdinal,
				left.measuredMass, left.charge) <
				std::tie(rightOrdinal, right.hypothesisOrdinal,
					right.measuredMass, right.charge);
		});
}

void scoreIndexedScanMvh(
	const FragmentIndex &index,
	MS2Scan &scan,
	const std::vector<IndexedCandidate> &survivors,
	std::vector<double> &sequenceIonMasses,
	std::vector<double> &aaForward,
	std::vector<double> &aaReverse,
	std::vector<char> &residues,
	IndexedSearchCounters &counters)
{
	for (const IndexedCandidate &candidate : survivors)
	{
		if (mergeIndexedPeptideIfPresent(
				index, candidate.peptideId, candidate.charge,
				scan.vpWeightSumTopPeptides))
		{
			continue;
		}
		double mvh = 0.0;
		++counters.exactMvhCalls;
		if (MVH::ScoreSequenceVsSpectrum(
				index.scoringSequence(candidate.peptideId),
				candidate.charge,
				&scan,
				&sequenceIonMasses,
				&aaForward,
				&aaReverse,
				mvh,
				&residues))
		{
			++counters.exactMvhAccepted;
			if (scan.vpWeightSumTopPeptides.size() < TOP_N ||
				mvh > scan.vpWeightSumTopPeptides[TOP_N - 1]->dScore)
			{
				Peptide peptide;
				index.materializePeptide(candidate.peptideId, peptide);
				scan.saveScore(
					mvh,
					{candidate.measuredMass, candidate.charge, &peptide},
					scan.vpWeightSumTopPeptides,
					"MVH", 2);
			}
		}
	}
	assignDdaResidualRanks(
		scan, sequenceIonMasses, aaForward, aaReverse, residues);
}

} // namespace sipros
