#include "MVH.h"
#include "PeptideIsotopeCalculator.h"
#include "peptide.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Snapshot
{
	std::string identified;
	std::string original;
	std::string protein;
	int begin = -1;
	char identifyPrefix = '?';
	char identifySuffix = '?';
	char originalPrefix = '?';
	char originalSuffix = '?';
};

void check(bool condition, const std::string &message)
{
	if (!condition)
		throw std::runtime_error(message);
}

std::vector<Snapshot> digest(const std::string &name,
							 const std::string &sequence)
{
	ProteinDatabase database;
	check(database.setProteinEntry(name, sequence),
		"setProteinEntry rejected a valid protein");
	std::vector<Snapshot> result;
	Peptide peptide;
	while (database.getNextPeptide(&peptide))
	{
		result.push_back({
			peptide.getPeptideSeq(),
			peptide.getOriginalPeptideSeq(),
			peptide.getProteinName(),
			peptide.getBeginPosProtein(),
			peptide.getIdentifyPrefix(),
			peptide.getIdentifySuffix(),
			peptide.getOriginalPrefix(),
			peptide.getOriginalSuffix()});
		check(result.size() < 100000,
			"protein peptide iterator did not terminate");
	}
	return result;
}

const Snapshot &findOne(const std::vector<Snapshot> &peptides,
							const std::string &identified)
{
	const auto found = std::find_if(peptides.begin(), peptides.end(),
		[&](const Snapshot &candidate)
		{
			return candidate.identified == identified;
		});
	check(found != peptides.end(), "missing peptide " + identified);
	check(std::count_if(peptides.begin(), peptides.end(),
		[&](const Snapshot &candidate)
		{
			return candidate.identified == identified;
		}) == 1, "duplicate peptide " + identified);
	return *found;
}

bool contains(const std::vector<Snapshot> &peptides,
			  const std::string &identified)
{
	return std::any_of(peptides.begin(), peptides.end(),
		[&](const Snapshot &candidate)
		{
			return candidate.identified == identified;
		});
}

void configure(ProNovoConfig::Profile profile,
			   const std::vector<std::string> &ptms,
			   int maxPtmCount)
{
	check(ProNovoConfig::load(profile), "could not load built-in profile");
	std::string error;
	check(ProNovoConfig::configureVariablePtms(ptms, maxPtmCount, error),
		error);
}

void checkRegularProteoformsAndContext()
{
	configure(ProNovoConfig::Profile::Regular, {"none"}, 1);
	const std::vector<Snapshot> peptides = digest(
		">met_protein descriptive text", "MACDEFGHIKAAAAAAAK");

	const Snapshot &intact = findOne(peptides, "[MACDEFGHIK]");
	check(intact.original == "[MACDEFGHIK]" && intact.protein == "met_protein",
		"retained-M metadata is incorrect");
	check(intact.begin == 0 && intact.identifyPrefix == '-' &&
		intact.originalPrefix == '-' && intact.identifySuffix == 'A' &&
		intact.originalSuffix == 'A',
		"retained-M coordinates or flanks are incorrect");

	const Snapshot &clipped = findOne(peptides, "[ACDEFGHIK]");
	check(clipped.original == "[ACDEFGHIK]" && clipped.begin == 1,
		"Met-excised coordinates are incorrect");
	check(clipped.identifyPrefix == '-' && clipped.originalPrefix == 'M' &&
		clipped.identifySuffix == 'A' && clipped.originalSuffix == 'A',
		"Met-excised identified/original flanks are incorrect");

	const Snapshot &acetylIntact = findOne(peptides, "[%MACDEFGHIK]");
	check(acetylIntact.original == "[MACDEFGHIK]" &&
		acetylIntact.begin == 0 && acetylIntact.identifyPrefix == '-',
		"acetylated retained-M metadata is incorrect");
	const Snapshot &acetylClipped = findOne(peptides, "[%ACDEFGHIK]");
	check(acetylClipped.original == "[ACDEFGHIK]" &&
		acetylClipped.begin == 1 && acetylClipped.identifyPrefix == '-' &&
		acetylClipped.originalPrefix == 'M',
		"acetylated Met-excised metadata is incorrect");

	findOne(peptides, "[AAAAAAAK]");
	check(!contains(peptides, "[%AAAAAAAK]"),
		"an internal peptide was protein-N-terminally acetylated");

	const std::vector<Snapshot> nonMet = digest(
		"non_met", "ACDEFGHIKAAAAAAAK");
	const Snapshot &nonMetTerminus = findOne(nonMet, "[ACDEFGHIK]");
	const Snapshot &acetylNonMet = findOne(nonMet, "[%ACDEFGHIK]");
	check(nonMetTerminus.begin == 0 && acetylNonMet.begin == 0 &&
		nonMetTerminus.originalPrefix == '-' && acetylNonMet.originalPrefix == '-',
		"non-M biological terminus metadata is incorrect");
	check(!contains(nonMet, "[CDEFGHIK]"),
		"non-M protein received an initiator-excision proteoform");
}

void checkTerminalAcetylConsumesPtmLimit()
{
	configure(ProNovoConfig::Profile::Regular, {"oxidation"}, 1);
	const std::vector<Snapshot> maxOne = digest("max_one", "MACDEFGHIK");
	findOne(maxOne, "[%MACDEFGHIK]");
	findOne(maxOne, "[M~ACDEFGHIK]");
	check(!contains(maxOne, "[%M~ACDEFGHIK]"),
		"terminal acetylation did not consume a PTM slot");

	configure(ProNovoConfig::Profile::Regular, {"oxidation"}, 2);
	const std::vector<Snapshot> maxTwo = digest("max_two", "MACDEFGHIK");
	findOne(maxTwo, "[%M~ACDEFGHIK]");

	configure(ProNovoConfig::Profile::Regular, {"none"}, 0);
	const std::vector<Snapshot> maxZero = digest("max_zero", "MACDEFGHIK");
	check(std::none_of(maxZero.begin(), maxZero.end(),
		[](const Snapshot &candidate)
		{
			return candidate.identified.find('%') != std::string::npos;
		}), "terminal acetylation bypassed a zero PTM maximum");
}

void checkLengthBoundaries()
{
	configure(ProNovoConfig::Profile::Regular, {"none"}, 0);
	const std::string clippedSixty = std::string(59, 'A') + "K";
	const std::vector<Snapshot> longProtein = digest(
		"long_boundary", "M" + clippedSixty);
	findOne(longProtein, "[" + clippedSixty + "]");
	check(!contains(longProtein, "[M" + clippedSixty + "]"),
		"over-length retained-M peptide was emitted");

	const std::vector<Snapshot> shortProtein = digest(
		"short_boundary", "MACDEGK");
	findOne(shortProtein, "[MACDEGK]");
	check(!contains(shortProtein, "[ACDEGK]"),
		"under-length Met-excised peptide was emitted");
}

void checkTerminalAcetylChemistry()
{
	configure(ProNovoConfig::Profile::Regular, {"none"}, 1);
	PeptideIsotopeCalculator calculator;
	const auto plain = calculator.calPrecursorEstimate("ACDEFGHIK");
	const auto acetyl = calculator.calPrecursorEstimate("%ACDEFGHIK");
	const double acetylMass = ProNovoConfig::getResidueMass("%");
	check(std::fabs((acetyl.mass - plain.mass) - acetylMass) < 1e-6,
		"terminal acetyl precursor mass is incorrect");
	check(std::fabs(acetylMass - 42.010565) < 1e-5,
		"compiled acetyl mass is incorrect");

	std::string plainSequence = "[ACDEFGHIK]";
	std::string acetylSequence = "[%ACDEFGHIK]";
	std::vector<double> plainCharged, plainForward, plainReverse;
	std::vector<double> acetylCharged, acetylForward, acetylReverse;
	std::vector<char> plainResidues, acetylResidues;
	check(MVH::CalculateSequenceIons(plainSequence, 2,
		MVH::bUseSmartPlusThreeModel, &plainCharged, &plainForward,
		&plainReverse, &plainResidues), "plain MVH ions failed");
	check(MVH::CalculateSequenceIons(acetylSequence, 2,
		MVH::bUseSmartPlusThreeModel, &acetylCharged, &acetylForward,
		&acetylReverse, &acetylResidues), "acetyl MVH ions failed");
	check(plainForward.size() == acetylForward.size() &&
		plainReverse.size() == acetylReverse.size(),
		"terminal acetyl changed fragment counts");
	for (size_t i = 0; i < plainForward.size(); ++i)
	{
		check(std::fabs((acetylForward[i] - plainForward[i]) - acetylMass) < 1e-6,
			"terminal acetyl did not shift a b ion");
	}
	// Conventional y1..y(L-1) ions do not include the protein N terminus.
	for (size_t i = 0; i + 1 < plainReverse.size(); ++i)
	{
		check(std::fabs(acetylReverse[i] - plainReverse[i]) < 1e-8,
			"terminal acetyl unexpectedly shifted a y ion");
	}
}

void checkSipLegacyBehavior()
{
	// A positive PTM limit makes this assertion independent of SIP's usual
	// zero-PTM default: terminal acetylation is disabled by profile, not merely
	// by the count limit.
	configure(ProNovoConfig::Profile::Sip, {"none"}, 1);
	const std::vector<Snapshot> peptides = digest(
		"sip_legacy", "MACDEFGHIKAAAAAAAK");
	const Snapshot &legacy = findOne(peptides, "[ACDEFGHIK]");
	check(legacy.begin == 0 && legacy.identifyPrefix == '-' &&
		legacy.originalPrefix == '-',
		"SIP Met-excised legacy coordinates changed");
	check(!contains(peptides, "[MACDEFGHIK]"),
		"SIP unexpectedly emitted retained initiator methionine");
	check(std::none_of(peptides.begin(), peptides.end(),
		[](const Snapshot &candidate)
		{
			return candidate.identified.find('%') != std::string::npos;
		}), "SIP unexpectedly emitted protein-N-terminal acetylation");
}

} // namespace

int main()
{
	try
	{
		checkRegularProteoformsAndContext();
		checkTerminalAcetylConsumesPtmLimit();
		checkLengthBoundaries();
		checkTerminalAcetylChemistry();
		checkSipLegacyBehavior();
		std::cout << "ok: Regular protein-N-terminal proteoforms and SIP legacy behavior\n";
		return 0;
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
