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
#include <map>
#include <tuple>

namespace sipros
{

namespace
{

bool candidateIdentityLess(const IndexedCandidate &left,
							   const IndexedCandidate &right)
{
	return std::tie(left.peptideId, left.charge, left.measuredMass) <
		std::tie(right.peptideId, right.charge, right.measuredMass);
}

bool mergeIndexedPeptideIfPresent(
	const FragmentIndex &index,
	uint32_t peptideId,
	MS2Scan &scan)
{
	const std::string_view sequence = index.peptideSequence(peptideId);
	for (const PeptideUnit *top : scan.vpWeightSumTopPeptides)
	{
		if (std::string_view(top->sIdentifiedPeptide) == sequence)
		{
			// Preserve the existing protein-name merge behavior, but only create
			// owning strings for the rare candidate that is already in the top N.
			return scan.mergePeptide(
				scan.vpWeightSumTopPeptides,
				std::string(sequence),
				std::string(index.proteinNames(peptideId)));
		}
	}
	return false;
}

} // namespace

IndexedSearchCounters &IndexedSearchCounters::operator+=(
	const IndexedSearchCounters &other)
{
	precursorHypotheses += other.precursorHypotheses;
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
	const std::vector<std::pair<double, int>> &precursorMassCharges,
	IndexedSearchScratch &scratch,
	std::vector<IndexedCandidate> &survivors,
	IndexedSearchCounters &counters)
{
	survivors.clear();
	if (scan.bSkip || scan.pPeakList == nullptr ||
		precursorMassCharges.empty() || index.peptideCount() == 0)
	{
		return;
	}
	struct Hypothesis
	{
		double mass = 0.0;
		uint32_t ordinal = 0;
	};
	std::map<int, std::vector<Hypothesis>> hypothesesByCharge;
	for (size_t i = 0; i < precursorMassCharges.size(); ++i)
	{
		const double mass = precursorMassCharges[i].first;
		const int charge = precursorMassCharges[i].second;
		if (mass > 0.0 && std::isfinite(mass) && charge > 0)
		{
			hypothesesByCharge[charge].push_back(
				{mass, static_cast<uint32_t>(i)});
			++counters.precursorHypotheses;
		}
	}

	const double precursorTolerance = ProNovoConfig::getMassAccuracyParentIon();
	const double fragmentTolerance = ProNovoConfig::getMassAccuracyFragmentIon();
	const auto &isotopeWindows = ProNovoConfig::getParentMassWindows();
	const uint16_t hitThreshold = static_cast<uint16_t>(
		std::min(ProNovoConfig::MinMatchedFragments,
			static_cast<int>(std::numeric_limits<uint16_t>::max())));

	for (const auto &chargeEntry : hypothesesByCharge)
	{
		const int precursorCharge = chargeEntry.first;
		std::vector<IndexedCandidate> &candidates = scratch.candidates;
		candidates.clear();
		if (scratch.candidatePositions.size() != index.peptideCount())
		{
			scratch.candidatePositions.assign(index.peptideCount(), 0);
			scratch.candidateEpoch = 0;
		}
		++scratch.candidateEpoch;
		if (scratch.candidateEpoch == 0)
		{
			std::fill(
				scratch.candidatePositions.begin(),
				scratch.candidatePositions.end(), 0);
			scratch.candidateEpoch = 1;
		}
		const uint64_t epoch =
			static_cast<uint64_t>(scratch.candidateEpoch) << 32;

		for (const Hypothesis &hypothesis : chargeEntry.second)
		{
			for (int isotopeWindow : isotopeWindows)
			{
				const double endpointA = hypothesis.mass -
					static_cast<double>(isotopeWindow) * index.minimumNeutronMass();
				const double endpointB = hypothesis.mass -
					static_cast<double>(isotopeWindow) * index.maximumNeutronMass();
				const double lower = std::min(endpointA, endpointB) - precursorTolerance;
				const double upper = std::max(endpointA, endpointB) + precursorTolerance;
				const auto range = index.peptideMassRange(lower, upper);
				for (uint32_t peptideId = range.first; peptideId < range.second; ++peptideId)
				{
					const IndexedPeptideRecord &record = index.peptide(peptideId);
					const double expected = record.precursorMass +
						static_cast<double>(isotopeWindow) *
							record.precursorNeutronMass;
					if (std::abs(hypothesis.mass - expected) > precursorTolerance)
					{
						continue;
					}
					const IndexedCandidate candidate{
						peptideId, hypothesis.mass, precursorCharge,
						hypothesis.ordinal};
					uint64_t &positionState =
						scratch.candidatePositions[peptideId];
					if ((positionState & 0xffffffff00000000ULL) != epoch)
					{
						const uint64_t position =
							static_cast<uint64_t>(candidates.size());
						positionState = epoch | (position + 1);
						candidates.push_back(candidate);
					}
					else
					{
						const size_t position = static_cast<size_t>(
							(positionState & 0xffffffffULL) - 1);
						IndexedCandidate &representative =
							candidates[position];
						if (std::tie(candidate.hypothesisOrdinal,
								candidate.measuredMass) <
							std::tie(representative.hypothesisOrdinal,
								representative.measuredMass))
						{
							representative = candidate;
						}
					}
				}
			}
		}

		std::sort(candidates.begin(), candidates.end(), candidateIdentityLess);
		counters.exactMassCandidates += candidates.size();
		if (candidates.empty())
		{
			continue;
		}
		counters.uniquePeptideChargeCandidates += candidates.size();

		const int maximumFragmentCharge = precursorCharge <= 2
			? 1
			: precursorCharge - 1;
		// Candidate records are peptide-id ordered, so each bounded precursor
		// block can be scored with a small counter array that remains in L1.
		// UINT32_MAX denotes a peptide that is not an exact precursor candidate.
		for (size_t candidateBegin = 0; candidateBegin < candidates.size();)
		{
			const uint32_t block = index.precursorBlockForPeptide(
				candidates[candidateBegin].peptideId);
			const uint32_t peptideBegin =
				index.precursorBlockPeptideBegin(block);
			size_t candidateEnd = candidateBegin + 1;
			while (candidateEnd < candidates.size() &&
				index.precursorBlockForPeptide(candidates[candidateEnd].peptideId) ==
					block)
			{
				++candidateEnd;
			}
			std::array<uint32_t, FragmentIndex::peptideBlockCapacity()> hitCounts;
			hitCounts.fill(std::numeric_limits<uint32_t>::max());
			for (size_t i = candidateBegin; i < candidateEnd; ++i)
			{
				hitCounts[candidates[i].peptideId - peptideBegin] = 0;
			}
			size_t candidatesBelowThreshold = hitThreshold == 0
				? 0
				: candidateEnd - candidateBegin;
			for (int fragmentCharge = 1;
				 fragmentCharge <= maximumFragmentCharge &&
				 candidatesBelowThreshold != 0;
				 ++fragmentCharge)
			{
				const double tolerance =
					static_cast<double>(fragmentCharge) * fragmentTolerance;
				for (double observedMz : scan.pPeakList->pPeaks)
				{
					const double neutralMass =
						static_cast<double>(fragmentCharge) *
						(observedMz - Proton);
					const auto postings = index.fragmentRange(
						block, neutralMass - tolerance, neutralMass + tolerance);
					if (postings.first == nullptr)
					{
						continue;
					}
					counters.fragmentPostingsVisited +=
						static_cast<uint64_t>(postings.second - postings.first);
					for (const FragmentPosting *posting = postings.first;
						 posting != postings.second;
						 ++posting)
					{
						uint32_t &state = hitCounts[posting->localPeptideId];
						if (state != std::numeric_limits<uint32_t>::max() &&
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
			for (size_t i = candidateBegin; i < candidateEnd; ++i)
			{
				const IndexedCandidate &candidate = candidates[i];
				if (hitCounts[candidate.peptideId - peptideBegin] >= hitThreshold)
				{
					survivors.push_back(candidate);
					++counters.fragmentGateSurvivors;
				}
			}
			candidateBegin = candidateEnd;
		}
	}

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
		if (mergeIndexedPeptideIfPresent(index, candidate.peptideId, scan))
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
					"MVH",
					2);
			}
		}
	}
}

} // namespace sipros
