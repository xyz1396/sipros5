#include "proNovoConfig.h"
#include "isotopologue.h"
#include "SiprosWorkflows.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using PtmDefinition = ProNovoConfig::PtmDefinition;

void require(bool condition, const std::string &message)
{
	if (!condition)
		throw std::runtime_error(message);
}

bool near(double actual, double expected, double tolerance = 1e-10)
{
	return std::abs(actual - expected) <= tolerance;
}

double lowestMass(const IsotopeDistribution &distribution)
{
	require(!distribution.vMass.empty(), "isotope distribution is empty");
	return *std::min_element(
		distribution.vMass.begin(), distribution.vMass.end());
}

double lowestMassForSequence(
	Isotopologue &isotopologue,
	const std::string &sequence)
{
	IsotopeDistribution distribution;
	require(isotopologue.computeIsotopicDistribution(sequence, distribution),
			"failed to compute isotope distribution for " + sequence);
	return lowestMass(distribution);
}

void requireCounts(const sipros::AtomCounts &actual,
					   const sipros::AtomCounts &expected,
				   const std::string &message)
{
	require(actual == expected, message);
}

const ProNovoConfig::PtmDefinition &ptmByToken(const std::string &token)
{
	const auto &catalog = ProNovoConfig::getPtmCatalog();
	const auto found = std::find_if(
		catalog.begin(), catalog.end(),
		[&](const ProNovoConfig::PtmDefinition &definition)
		{ return definition.token == token; });
	require(found != catalog.end(), "missing PTM token " + token);
	return *found;
}

void requireSameComposition(
	const sipros::SourcedComposition &actual,
	const sipros::SourcedComposition &expected,
	const std::string &message)
{
	for (size_t source = 0; source < sipros::IsotopeSourceCount; ++source)
		require(actual.atoms[source] == expected.atoms[source], message);
}

void requireDistributionNear(const IsotopeDistribution &actual,
							 const IsotopeDistribution &expected,
							 const std::string &message)
{
	require(actual.vMass.size() == expected.vMass.size() &&
				actual.vProb.size() == expected.vProb.size(),
			message + " (different sizes)");
	for (size_t index = 0; index < actual.vMass.size(); ++index)
	{
		require(near(actual.vMass[index], expected.vMass[index], 1e-8) &&
					near(actual.vProb[index], expected.vProb[index], 1e-10),
				message + " (different peak)");
	}
}

sipros::SourcedComposition compositionForSymbols(
	const std::string &symbols,
	const std::string &context)
{
	sipros::SourcedComposition composition;
	const auto &chemistry =
		ProNovoConfig::configIsotopologue.mResidueSourcedComposition;
	for (char symbol : symbols)
	{
		const auto found = chemistry.find(std::string(1, symbol));
		require(found != chemistry.end(),
				context + " is missing chemistry symbol " + symbol);
		composition += found->second;
	}
	return composition;
}

void requireNonnegativeComposition(
	const sipros::SourcedComposition &composition,
	const std::string &message)
{
	for (size_t source = 0; source < sipros::IsotopeSourceCount; ++source)
	{
		for (size_t element = 0; element < sipros::ElementCount; ++element)
			require(composition.atoms[source][element] >= 0, message);
	}
}

void requireIonDistributionNear(
	const std::vector<double> &actualMass,
	const std::vector<double> &actualProbability,
	IsotopeDistribution expected,
	double massCorrection,
	const std::string &message)
{
	for (double &mass : expected.vMass)
		mass += massCorrection;
	require(!actualMass.empty() &&
			actualMass.size() == actualProbability.size(),
			message + " (empty or malformed envelope)");
	const double actualTotal = std::accumulate(
		actualProbability.begin(), actualProbability.end(), 0.0);
	const double expectedTotal = std::accumulate(
		expected.vProb.begin(), expected.vProb.end(), 0.0);
	require(near(actualTotal, expectedTotal, 1e-10),
			message + " (different normalization)");
	const auto actualApex = std::max_element(
		actualProbability.begin(), actualProbability.end());
	const auto expectedApex = std::max_element(
		expected.vProb.begin(), expected.vProb.end());
	const size_t actualApexIndex = static_cast<size_t>(
		std::distance(actualProbability.begin(), actualApex));
	const size_t expectedApexIndex = static_cast<size_t>(
		std::distance(expected.vProb.begin(), expectedApex));
	require(near(actualMass[actualApexIndex],
				expected.vMass[expectedApexIndex], 1e-5) &&
				near(*actualApex, *expectedApex, 1e-5),
			message + " (different envelope apex)");
}

void checkProductIonNetCompositionCase(
	const std::string &decoratedPeptide,
	const std::vector<std::string> &residues,
	const std::string &nTermPtms,
	const std::string &cTermPtms,
	char sipAtom,
	double sipAbundance,
	const std::string &label)
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to load SIP profile for " + label);
	require(ProNovoConfig::applySipAbundance(sipAtom, sipAbundance),
			"failed to set partial SIP abundance for " + label);

	Isotopologue &isotopologue = ProNovoConfig::configIsotopologue;
	std::vector<std::vector<double>> yMass;
	std::vector<std::vector<double>> yProbability;
	std::vector<std::vector<double>> bMass;
	std::vector<std::vector<double>> bProbability;
	require(isotopologue.computeProductIon(
			decoratedPeptide, yMass, yProbability, bMass, bProbability),
			"failed to compute product ions for " + label);
	require(yMass.size() == residues.size() - 1 &&
				yProbability.size() == yMass.size() &&
				bMass.size() == residues.size() - 1 &&
				bProbability.size() == bMass.size(),
			"wrong product-ion count for " + label);

	sipros::SourcedComposition bComposition =
		compositionForSymbols(nTermPtms, label + " N terminus");
	for (size_t cleavage = 0; cleavage + 1 < residues.size(); ++cleavage)
	{
		bComposition += compositionForSymbols(
			residues[cleavage], label + " b residue");
		requireNonnegativeComposition(
			bComposition, label + " b-ion net composition is negative");
		IsotopeDistribution expected;
		require(isotopologue.computeIsotopicDistribution(
				bComposition, expected),
			"failed direct b-ion convolution for " + label);
		requireIonDistributionNear(
			bMass[cleavage], bProbability[cleavage], expected,
			0.0,
			label + " b-ion differs from net-composition convolution");
	}

	sipros::SourcedComposition yComposition =
		isotopologue.mResidueSourcedComposition.at("Nterm") +
		isotopologue.mResidueSourcedComposition.at("Cterm") +
		compositionForSymbols(cTermPtms, label + " C terminus");
	for (size_t offset = 0; offset + 1 < residues.size(); ++offset)
	{
		const size_t residue = residues.size() - 1 - offset;
		yComposition += compositionForSymbols(
			residues[residue], label + " y residue");
		requireNonnegativeComposition(
			yComposition, label + " y-ion net composition is negative");
		IsotopeDistribution expected;
		require(isotopologue.computeIsotopicDistribution(
				yComposition, expected),
			"failed direct y-ion convolution for " + label);
		requireIonDistributionNear(
			yMass[offset], yProbability[offset], expected,
			0.0,
			label + " y-ion differs from net-composition convolution");
	}
}

void checkProductIonExactMassesAndPrecursorConservation()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to load SIP profile for exact product-ion masses");
	Isotopologue &isotopologue = ProNovoConfig::configIsotopologue;
	std::vector<std::vector<double>> yMass;
	std::vector<std::vector<double>> yProbability;
	std::vector<std::vector<double>> bMass;
	std::vector<std::vector<double>> bProbability;
	const std::string peptide = "SATPAQAQAVHK";
	require(isotopologue.computeProductIon(
			"[" + peptide + "]", yMass, yProbability, bMass, bProbability),
			"failed to compute exact product-ion masses");

	const double proton = ProNovoConfig::getProtonMass();
	const double expectedB1Mz =
		3.0 * 12.000000 + 5.0 * 1.007825 +
		2.0 * 15.994915 + 14.003074 + proton;
	const double expectedY1Mz =
		6.0 * 12.000000 + 14.0 * 1.007825 +
		2.0 * 15.994915 + 2.0 * 14.003074 + proton;
	require(near(bMass.front().front() + proton, expectedB1Mz, 1e-9),
			"b1 m/z differs from neutral composition plus proton");
	require(near(yMass.front().front() + proton, expectedY1Mz, 1e-9),
			"y1 m/z differs from neutral composition plus proton");

	IsotopeDistribution reconstructedPrecursor = isotopologue.sum(
		IsotopeDistribution(bMass.back(), bProbability.back()),
		IsotopeDistribution(yMass.front(), yProbability.front()));
	IsotopeDistribution directPrecursor;
	require(isotopologue.computePeptideIsotopicDistribution(
			"[" + peptide + "]", directPrecursor),
			"failed to compute direct precursor distribution");
	sipros::SourcedComposition precursorComposition;
	IsotopeDistribution compositionPrecursor;
	require(isotopologue.computeSourcedComposition(
			peptide, precursorComposition) &&
			isotopologue.computeIsotopicDistribution(
				precursorComposition, compositionPrecursor),
			"failed to compute composition-aware precursor distribution");
	requireDistributionNear(
		directPrecursor, compositionPrecursor,
		"decorated direct precursor convolution changed the composition");
	require(near(lowestMass(reconstructedPrecursor),
				lowestMass(directPrecursor), 1e-8),
			"b/y terminal correction changed reconstructed precursor base mass");
	require(near(reconstructedPrecursor.getMostAbundantMass(),
				directPrecursor.getMostAbundantMass(), 1e-5),
			"b/y terminal correction changed reconstructed precursor apex mass");
}

void checkObservedFragmentChargePolicy()
{
	require(sipros::observedPeakChargeMatches(0, 1),
			"unknown raw fragment charge was rejected");
	require(sipros::observedPeakChargeMatches(1, 1),
			"charge-one raw fragment was rejected");
	require(!sipros::observedPeakChargeMatches(2, 1),
			"charge-two raw fragment matched a charge-one theoretical ion");
	require(!sipros::observedPeakChargeMatches(-1, 1),
			"invalid raw fragment charge was accepted");
}

void checkProductIonNetComposition()
{
	checkProductIonExactMassesAndPrecursorConservation();
	checkObservedFragmentChargePolicy();
	checkProductIonNetCompositionCase(
		"[AN!CDEFGK]",
		{"A", "N!", "C", "D", "E", "F", "G", "K"},
		"", "", 'N', 0.43, "deamidation");
	checkProductIonNetCompositionCase(
		"[AC(DEFGHK]",
		{"A", "C(", "D", "E", "F", "G", "H", "K"},
		"", "", 'C', 0.57, "fixed-CAM S-nitrosylation");
	checkProductIonNetCompositionCase(
		"[%AS2CDEFGK]",
		{"A", "S2", "C", "D", "E", "F", "G", "K"},
		"%", "", 'O', 0.61,
		"neutral loss with N-terminal acetylation");
}

void checkNoUninitializedFallback()
{
	require(ProNovoConfig::atomIndex('C') == -1,
			"atom lookup must not fall back before profile initialization");
	std::map<std::string, std::string> ptms;
	require(!ProNovoConfig::getPTMinfo(ptms),
			"PTM lookup must fail before profile initialization");
	bool threw = false;
	try
	{
		(void)ProNovoConfig::getNaturalAtomIsotopeProbabilities(0);
	}
	catch (const std::logic_error &)
	{
		threw = true;
	}
	require(threw,
			"natural isotope lookup must fail before profile initialization");
}

void checkCommonChemistry()
{
	Isotopologue &iso = ProNovoConfig::configIsotopologue;
	require(iso.vAtomIsotopicDistribution.size() == 6,
			"unexpected atom distribution count");
	require(ProNovoConfig::atomIndex('C') == 0 &&
				ProNovoConfig::atomIndex('P') == 4 &&
				ProNovoConfig::atomIndex('S') == 5,
			"unexpected atom order");

	const auto sourcedCysteine =
		iso.mResidueSourcedComposition.find("C");
	require(sourcedCysteine != iso.mResidueSourcedComposition.end(),
			"sourced cysteine composition is missing");
	requireCounts(
		sourcedCysteine->second[sipros::IsotopeSource::Biosynthetic],
		{3, 5, 1, 1, 0, 1},
		"cysteine biosynthetic formula is wrong");
	requireCounts(
		sourcedCysteine->second[sipros::IsotopeSource::ReagentNatural],
		{2, 3, 1, 1, 0, 0},
		"CAM reagent-natural formula is wrong");
	requireCounts(
		sourcedCysteine->second[sipros::IsotopeSource::DigestionSolvent],
		{}, "cysteine unexpectedly contains digestion solvent");

	const auto &phosphorus = iso.vNaturalAtomIsotopicDistribution[4];
	require(phosphorus.vMass == std::vector<double>{30.973762} &&
				phosphorus.vProb == std::vector<double>{1.0},
			"phosphorus is not real monoisotopic P31");
	require(ProNovoConfig::getChemistryProfileId() ==
				"sipros5/source-aware-cam-tryptic-water/v1",
			"unexpected chemistry profile ID");

	std::string error;
	require(ProNovoConfig::validatePreparationChemistry(iso, error), error);

	const double blockedCysteineMass =
		5.0 * 12.0 + 8.0 * 1.007825 + 2.0 * 15.994915 +
		2.0 * 14.003074 + 31.972071;
	require(near(ProNovoConfig::getResidueMass("C"),
				 blockedCysteineMass, 1e-6),
			"natural blocked-cysteine mass changed");
	require(near(blockedCysteineMass, 160.030649, 1e-6),
			"blocked-cysteine formula mass is wrong");
	const double camMass =
		2.0 * 12.0 + 3.0 * 1.007825 + 15.994915 + 14.003074;
	require(near(camMass, 57.021464, 1e-6),
			"CAM formula mass is wrong");

	const auto nTerm = iso.mResidueSourcedComposition.at("Nterm");
	const auto cTerm = iso.mResidueSourcedComposition.at("Cterm");
	requireCounts(nTerm[sipros::IsotopeSource::DigestionSolvent],
				  {0, 1, 0, 0, 0, 0},
				  "N terminus is not digestion-solvent hydrogen");
	requireCounts(cTerm[sipros::IsotopeSource::DigestionSolvent],
				  {0, 1, 1, 0, 0, 0},
				  "C terminus is not digestion-solvent hydroxyl");

	sipros::SourcedComposition peptide;
	require(iso.computeSourcedComposition("C", peptide),
			"failed to compute sourced peptide composition");
	requireCounts(peptide[sipros::IsotopeSource::Biosynthetic],
				  {3, 5, 1, 1, 0, 1},
				  "peptide biosynthetic formula is wrong");
	requireCounts(peptide.naturalSourceTotal(),
				  {2, 5, 2, 1, 0, 0},
				  "peptide natural-source formula is wrong");
	requireCounts(peptide.total(), {5, 10, 3, 2, 0, 1},
				  "peptide total formula is wrong");
	require(near(lowestMassForSequence(iso, "C"),
				 178.041214, 1e-6),
			"natural peptide monoisotopic mass is wrong");

	IsotopeDistribution sourcedDistribution;
	IsotopeDistribution flatDistribution;
	const sipros::SourcedComposition flatComposition =
		sipros::compositionFrom(sipros::IsotopeSource::Biosynthetic,
								{5, 10, 3, 2, 0, 1});
	require(iso.computeIsotopicDistribution(peptide, sourcedDistribution) &&
				iso.computeIsotopicDistribution(
					flatComposition, flatDistribution),
			"failed to compute natural peptide isotope distribution");
	requireDistributionNear(sourcedDistribution, flatDistribution,
						"natural sourced and flat convolutions differ");

	sipros::SourcedComposition phospho = sipros::compositionFrom(
		sipros::IsotopeSource::ReagentNatural,
		{0, 1, 3, 0, 1, 0});
	IsotopeDistribution phosphoDistribution;
	require(iso.computeIsotopicDistribution(phospho, phosphoDistribution) &&
				near(lowestMass(phosphoDistribution), 79.966332, 1e-6),
			"real-phosphorus HPO3 mass is wrong");
}

void checkSourceIsolationAndEndpointMasses()
{
	struct Target
	{
		char atom;
		size_t element;
		double expectedEndpointMass;
	};
	const std::vector<Target> targets{
		{'C', 0, 181.051279},
		{'H', 1, 183.072599},
		{'N', 3, 179.038249},
		{'O', 2, 180.045459},
		{'S', 5, 180.037010}};

	for (const Target &target : targets)
	{
		require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
				"failed to reset SIP profile");
		const auto naturalBefore = ProNovoConfig::configIsotopologue
								   .vNaturalAtomIsotopicDistribution;
		require(ProNovoConfig::applySipAbundance(target.atom, 0.75),
				"failed to apply source-isolation abundance");

		for (sipros::IsotopeSource source :
			 {sipros::IsotopeSource::ReagentNatural,
			  sipros::IsotopeSource::DigestionSolvent})
		{
			sipros::AtomCounts oneAtom{};
			oneAtom[target.element] = 1;
			IsotopeDistribution actual;
			require(ProNovoConfig::configIsotopologue
						.computeIsotopicDistribution(
							sipros::compositionFrom(source, oneAtom), actual),
					"failed natural-source one-atom convolution");
			requireDistributionNear(
				actual, naturalBefore[target.element],
				"natural-source atom followed SIP enrichment");
		}
		sipros::AtomCounts oneBiosyntheticAtom{};
		oneBiosyntheticAtom[target.element] = 1;
		IsotopeDistribution biosyntheticDistribution;
		require(ProNovoConfig::configIsotopologue
					.computeIsotopicDistribution(
						sipros::compositionFrom(
							sipros::IsotopeSource::Biosynthetic,
							oneBiosyntheticAtom),
						biosyntheticDistribution),
			"failed biosynthetic one-atom convolution");
		requireDistributionNear(
			biosyntheticDistribution,
			ProNovoConfig::configIsotopologue
				.vAtomIsotopicDistribution[target.element],
			"biosynthetic atom did not follow SIP enrichment");
		const auto &naturalAfter = ProNovoConfig::configIsotopologue
								  .vNaturalAtomIsotopicDistribution;
		require(naturalAfter.size() == naturalBefore.size(),
			"immutable natural isotope table changed size");
		for (size_t element = 0; element < naturalBefore.size(); ++element)
		{
			require(naturalAfter[element].vMass == naturalBefore[element].vMass &&
						naturalAfter[element].vProb == naturalBefore[element].vProb,
					"immutable natural isotope distributions changed");
		}

		require(ProNovoConfig::applySipAbundance(target.atom, 1.0),
				"failed to apply endpoint abundance");
		require(near(lowestMassForSequence(
						 ProNovoConfig::configIsotopologue, "C"),
				 target.expectedEndpointMass, 2e-6),
				std::string("natural-source atoms shifted at endpoint ") +
					target.atom);
	}
}

void checkSipMassCacheRefreshAndTargetSelection()
{
	std::string error;
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to load SIP profile for mass-cache tests");
	require(ProNovoConfig::configureVariablePtms(
			{"phosphorylation", "deamidation"}, -1, error), error);

	require(ProNovoConfig::applySipAbundance('O', 0.0),
			"failed to apply light-oxygen endpoint");
	const double lightPhospho = ProNovoConfig::getResidueMass("@");
	const double lightFixedCys = ProNovoConfig::getResidueMass("C");
	const double naturalTerminalWater =
		ProNovoConfig::precalcMasses.dCtermOH2;
	require(near(ProNovoConfig::pdAAMassFragment.find('@'), lightPhospho),
			"Comet phosphorylation cache disagrees with compiled residue mass");

	require(ProNovoConfig::applySipAbundance('O', 1.0),
			"failed to apply heavy-oxygen endpoint");
	const double oxygen18Delta = 17.999160 - 15.994915;
	require(near(ProNovoConfig::getResidueMass("@") - lightPhospho,
				 3.0 * oxygen18Delta, 1e-6),
			"phosphorylation cache did not follow all three biosynthetic oxygens");
	require(near(ProNovoConfig::getResidueMass("C") - lightFixedCys,
				 oxygen18Delta, 1e-6),
			"fixed-CAM reagent oxygen incorrectly followed O18 enrichment");
	require(near(ProNovoConfig::pdAAMassFragment.find('@'),
				 ProNovoConfig::getResidueMass("@"), 1e-10),
			"Comet residue cache remained stale after abundance change");
	require(near(ProNovoConfig::precalcMasses.dCtermOH2,
				 naturalTerminalWater, 1e-10),
			"digestion-solvent terminal water followed O18 enrichment");

	require(ProNovoConfig::applySipAbundance('N', 0.0),
			"failed to apply light-nitrogen endpoint");
	const double lightDeamidation = ProNovoConfig::getResidueMass("!");
	require(ProNovoConfig::applySipAbundance('N', 1.0),
			"failed to apply heavy-nitrogen endpoint");
	require(near(ProNovoConfig::getResidueMass("!") - lightDeamidation,
				 -(15.000109 - 14.003074), 1e-6),
			"deamidation did not remove the target-specific N15 mass");

	struct Target
	{
		char atom;
		int massNumber;
		double neutronSpacing;
	};
	const std::vector<Target> targets{
		{'C', 13, 13.003355 - 12.000000},
		{'H', 2, 2.014102 - 1.007825},
		{'N', 15, 15.000109 - 14.003074},
		{'O', 18, (17.999160 - 15.994915) / 2.0},
		{'S', 34, (33.967867 - 31.972071) / 2.0}};
	for (const Target &target : targets)
	{
		require(ProNovoConfig::selectSipTarget(
				target.atom, target.massNumber, error), error);
		require(near(ProNovoConfig::getNeutronMass(),
					 target.neutronSpacing, 1e-9),
				std::string("wrong precursor isotope spacing for ") + target.atom);
	}
	require(!ProNovoConfig::selectSipTarget('O', 17, error) && !error.empty(),
			"unsupported isotope mass number was accepted");
}

void checkRegularProfile()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to load Regular profile");
	require(ProNovoConfig::getSearchType() == "Regular",
			"wrong Regular search type");
	require(ProNovoConfig::getSearchName() == "SE",
			"wrong Regular search name");
	require(ProNovoConfig::getMaxPTMcount() == 3,
			"wrong Regular maximum PTM count");
	require(ProNovoConfig::getMinPeptideLength() == 7 &&
				ProNovoConfig::getMaxPeptideLength() == 60,
			"wrong peptide length defaults");
	require(near(ProNovoConfig::getMassAccuracyParentIon(), 0.01) &&
				near(ProNovoConfig::getMassAccuracyFragmentIon(), 0.01),
			"wrong mass tolerance defaults");

	std::map<std::string, std::string> ptms;
	require(ProNovoConfig::getPTMinfo(ptms), "failed to read Regular PTMs");
	require(ptms == std::map<std::string, std::string>(
					{{"!", "NQ"}, {"~", "M"}}),
			"wrong Regular variable PTMs");
	checkCommonChemistry();

	require(ProNovoConfig::applySipAbundance('C', 0.25),
			"failed to apply C13 abundance");
	require(near(ProNovoConfig::configIsotopologue
					 .vAtomIsotopicDistribution[0]
					 .vProb[1],
				 0.25),
			"biological carbon did not receive C13 abundance");
	require(ProNovoConfig::configIsotopologue
				.vNaturalAtomIsotopicDistribution[0]
				.vProb[1] == 0.0107,
			"natural carbon pool incorrectly followed C13 abundance");

	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to reload Regular profile");
	require(near(ProNovoConfig::configIsotopologue
					 .vAtomIsotopicDistribution[0]
					 .vProb[1],
				 0.0107),
			"profile reload did not restore natural carbon");
}

void checkUnifiedRegularProductIons()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to load Regular profile for product-ion test");
	Isotopologue &isotopologue = ProNovoConfig::configIsotopologue;
	std::vector<std::vector<double>> yMass;
	std::vector<std::vector<double>> yProbability;
	std::vector<std::vector<double>> bMass;
	std::vector<std::vector<double>> bProbability;
	require(isotopologue.computeProductIon(
			"[SATPAQAQAVHK]", yMass, yProbability, bMass, bProbability),
			"failed to compute Regular product ions");

	IsotopeDistribution expectedB1 =
		isotopologue.vResidueIsotopicDistribution.at("S");
	sipros::SourcedComposition expectedYTermComposition =
		isotopologue.mResidueSourcedComposition.at("Nterm") +
		isotopologue.mResidueSourcedComposition.at("Cterm");
	IsotopeDistribution expectedYTerm;
	require(isotopologue.computeIsotopicDistribution(
			expectedYTermComposition, expectedYTerm),
		"failed to build direct Regular y-ion terminus");
	IsotopeDistribution expectedY1 = isotopologue.sum(
		isotopologue.vResidueIsotopicDistribution.at("K"), expectedYTerm);
	requireDistributionNear(
		IsotopeDistribution(bMass.front(), bProbability.front()),
		expectedB1,
		"Regular b1 differs from its neutral composition");
	requireDistributionNear(
		IsotopeDistribution(yMass.front(), yProbability.front()),
		expectedY1,
		"Regular y1 differs from its neutral composition");
}

void checkResiduePtmStateCache()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to load SIP profile for residue-PTM cache test");
	require(ProNovoConfig::applySipAbundance('C', 0.57),
			"failed to apply SIP abundance for residue-PTM cache test");
	Isotopologue &isotopologue = ProNovoConfig::configIsotopologue;
	for (const std::string &state : {"N!", "M~", "C(", "S2"})
	{
		sipros::SourcedComposition cachedComposition;
		IsotopeDistribution cachedDistribution;
		require(isotopologue.getCachedResidueState(
				state, cachedComposition, cachedDistribution),
				"missing cached residue-PTM state " + state);
		sipros::SourcedComposition expectedComposition =
			isotopologue.mResidueSourcedComposition.at(state.substr(0, 1));
		expectedComposition +=
			isotopologue.mResidueSourcedComposition.at(state.substr(1, 1));
		requireSameComposition(
			cachedComposition, expectedComposition,
			"cached sourced composition differs for " + state);
		requireNonnegativeComposition(
			cachedComposition,
			"cached sourced composition is negative for " + state);
		IsotopeDistribution expectedDistribution;
		require(isotopologue.computeIsotopicDistribution(
				cachedComposition, expectedDistribution),
				"failed direct convolution for cached state " + state);
		requireDistributionNear(
			cachedDistribution, expectedDistribution,
			"cached envelope differs for " + state);
	}
}

void checkPtmCatalogAndSelectors()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to load Regular profile for PTM catalog tests");
	const auto &catalog = ProNovoConfig::getPtmCatalog();
	require(catalog.size() == 13, "PTM catalog is incomplete");
	struct ExpectedPtm
	{
		const char *token;
		const char *sites;
		double externalShift;
		double internalShiftWithFixedCam;
	};
	const std::vector<ExpectedPtm> expected{
		{"~", "M", 15.994915, 15.994915},
		{"!", "NQ", 0.984016, 0.984016},
		{"@", "STYHD", 79.966332, 79.966332},
		{">", "STYHD", 79.966332, 79.966332},
		{"<", "ST", 79.966332, 79.966332},
		{"%", "K", 42.010565, 42.010565},
		{"^", "KRED", 14.015650, 14.015650},
		{"&", "KR", 28.031300, 28.031300},
		{"*", "K", 42.046950, 42.046950},
		{"(", "C", 28.990164, -28.031300},
		{")", "Y", 44.985079, 44.985079},
		{"/", "C", 57.021464, 0.0},
		{"$", "D", 45.987721, 45.987721}};
	for (const ExpectedPtm &entry : expected)
	{
		const PtmDefinition definition = ptmByToken(entry.token);
		require(definition.sites == entry.sites,
				std::string("wrong PTM sites for ") + entry.token);
		require(near(definition.externalMonoisotopicShift,
					 entry.externalShift, 1e-6),
				std::string("wrong external PTM shift for ") + entry.token);
		require(near(ProNovoConfig::getResidueMass(entry.token),
					 entry.internalShiftWithFixedCam, 1e-6),
				std::string("wrong internal PTM shift for ") + entry.token);
	}
	require(ptmByToken("~").regularDefault &&
			ptmByToken("!").regularDefault &&
			!ptmByToken("@").regularDefault,
			"Regular PTM default flags are wrong");
	require(!ptmByToken("/").selectable,
			"variable IAA must not be selectable while fixed CAM is active");
	require(near(ProNovoConfig::getResidueMass("1"), 0.0, 1e-8) &&
			near(ProNovoConfig::getResidueMass("2"), -18.010565, 1e-6),
			"phosphorylation neutral-loss helper formulas are wrong");
	require(std::none_of(
			catalog.begin(), catalog.end(),
			[](const PtmDefinition &definition)
			{ return definition.token == "1" || definition.token == "2"; }),
			"neutral-loss helpers leaked into the selectable catalog");

	const auto &sourced =
		ProNovoConfig::configIsotopologue.mResidueSourcedComposition;
	requireCounts(sourced.at("@")[sipros::IsotopeSource::Biosynthetic],
				  {0, 1, 3, 0, 1, 0},
				  "phosphorylation is not biosynthetic");
	requireCounts(sourced.at(")")[sipros::IsotopeSource::Biosynthetic],
				  {0, -1, 0, 0, 0, 0},
				  "nitration did not remove biosynthetic hydrogen");
	requireCounts(sourced.at(")")[sipros::IsotopeSource::ReagentNatural],
				  {0, 0, 2, 1, 0, 0},
				  "nitration is not reagent-natural");
	requireCounts(sourced.at("(")[sipros::IsotopeSource::Biosynthetic],
				  {0, -1, 0, 0, 0, 0},
				  "S-nitrosylation did not remove biosynthetic hydrogen");
	requireCounts(sourced.at("(")[sipros::IsotopeSource::ReagentNatural],
				  {-2, -3, 0, 0, 0, 0},
				  "fixed-CAM S-nitrosylation replacement is wrong");

	std::string error;
	require(ProNovoConfig::configureVariablePtms(
			{"phosphorylation", "~"}, 2, error), error);
	std::map<std::string, std::string> ptms;
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>(
					{{"@", "STYHD"}, {"~", "M"}}) &&
			ProNovoConfig::getMaxPTMcount() == 2,
			"name/token PTM selection is not exact");

	require(ProNovoConfig::configureVariablePtms(
			{">", "phosphorylation-loss-hpo3-h2o"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>(
					{{"<to2", "ST"}, {">to1", "STYHD"}}),
			"neutral-loss search keys are wrong");
	require(ProNovoConfig::getNeutralLossList() ==
			std::vector<std::pair<std::string, std::string>>(
				{{">", "1"}, {"<", "2"}}),
			"neutral-loss replacements were not activated");
	require(ProNovoConfig::getMaxPTMcount() == 3,
			"omitted Regular PTM maximum did not restore profile default");

	const auto beforeUnknownPtms = ptms;
	const auto beforeUnknownLosses = ProNovoConfig::getNeutralLossList();
	const int beforeUnknownMaximum = ProNovoConfig::getMaxPTMcount();
	require(!ProNovoConfig::configureVariablePtms(
			{"not-a-ptm"}, 7, error) && !error.empty(),
			"unknown PTM selector unexpectedly succeeded");
	std::map<std::string, std::string> afterUnknownPtms;
	ProNovoConfig::getPTMinfo(afterUnknownPtms);
	require(afterUnknownPtms == beforeUnknownPtms &&
			ProNovoConfig::getNeutralLossList() == beforeUnknownLosses &&
			ProNovoConfig::getMaxPTMcount() == beforeUnknownMaximum,
			"unknown PTM selector was not atomic");
	require(!ProNovoConfig::configureVariablePtms(
			{"none", "oxidation"}, -1, error),
			"exclusive PTM selector 'none' accepted another selector");
	require(!ProNovoConfig::configureVariablePtms(
			{"oxidation"}, -2, error),
			"invalid PTM maximum unexpectedly succeeded");

	require(ProNovoConfig::configureVariablePtms({"none"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.empty() &&
			ProNovoConfig::getNeutralLossList().empty(),
			"PTM selector 'none' did not close all variable PTMs");
	require(ProNovoConfig::configureVariablePtms({"default"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>(
					{{"!", "NQ"}, {"~", "M"}}),
			"PTM selector 'default' did not restore Regular defaults");
	require(ProNovoConfig::configureVariablePtms({"all"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.size() == 12 &&
			ptms.find("/") == ptms.end(),
			"PTM selector 'all' included unavailable IAA or omitted a PTM");
	require(ProNovoConfig::getNeutralLossList().size() == 2,
			"PTM selector 'all' omitted neutral-loss rules");
	require(!ProNovoConfig::configureVariablePtms({"/"}, -1, error),
			"variable IAA was selectable with fixed CAM active");

	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to reload Regular PTM defaults");
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>(
					{{"!", "NQ"}, {"~", "M"}}) &&
			ProNovoConfig::getMaxPTMcount() == 3 &&
			ProNovoConfig::getNeutralLossList().empty(),
			"profile reload did not reset variable PTMs");
}

void checkFixedPtmConfiguration()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to load Regular profile for fixed PTM tests");
	const auto &fixedCatalog = ProNovoConfig::getFixedPtmCatalog();
	require(fixedCatalog.size() == 1 &&
			fixedCatalog.front().name == "carbamidomethyl" &&
			fixedCatalog.front().sites == "C" &&
			fixedCatalog.front().profileDefault &&
			near(fixedCatalog.front().externalMonoisotopicShift,
				 57.021464, 1e-6),
			"fixed PTM catalog is wrong");
	require(ProNovoConfig::getEnabledFixedPtmNames() ==
			std::vector<std::string>{"carbamidomethyl"},
			"fixed CAM is not enabled by default");
	const auto &fixedCamPtmCatalog = ProNovoConfig::getPtmCatalog();
	const auto fixedCamIaa = std::find_if(
		fixedCamPtmCatalog.begin(), fixedCamPtmCatalog.end(),
		[](const ProNovoConfig::PtmDefinition &definition)
		{ return definition.token == "/"; });
	require(fixedCamIaa != fixedCamPtmCatalog.end() &&
			!fixedCamIaa->selectable,
			"fixed-CAM PTM catalog view made variable IAA selectable");

	sipros::SourcedComposition defaultFixedCysteine;
	require(ProNovoConfig::configIsotopologue.computeSourcedComposition(
			"C", defaultFixedCysteine),
			"failed to read default fixed-CAM cysteine");
	std::string error;
	require(ProNovoConfig::configureFixedPtms({"none"}, error), error);
	require(ProNovoConfig::getEnabledFixedPtmNames().empty(),
			"fixed selector 'none' did not disable CAM");
	const auto &naturalCysteinePtmCatalog = ProNovoConfig::getPtmCatalog();
	const auto naturalCysteineIaa = std::find_if(
		naturalCysteinePtmCatalog.begin(), naturalCysteinePtmCatalog.end(),
		[](const ProNovoConfig::PtmDefinition &definition)
		{ return definition.token == "/"; });
	require(naturalCysteineIaa != naturalCysteinePtmCatalog.end() &&
			naturalCysteineIaa->selectable && !fixedCamIaa->selectable,
			"fixed-chemistry catalog snapshots mutate across selector changes");
	require(ProNovoConfig::getChemistryProfileId() ==
			"sipros5/source-aware-natural-cys-tryptic-water/v1",
			"natural-cysteine chemistry profile ID is wrong");
	require(ProNovoConfig::validatePreparationChemistry(
			ProNovoConfig::configIsotopologue, error), error);
	const auto &naturalCys = ProNovoConfig::configIsotopologue
							 .mResidueSourcedComposition.at("C");
	requireCounts(naturalCys[sipros::IsotopeSource::Biosynthetic],
				  {3, 5, 1, 1, 0, 1},
				  "natural cysteine formula changed when CAM was disabled");
	requireCounts(naturalCys[sipros::IsotopeSource::ReagentNatural], {},
				  "disabled CAM atoms remain on cysteine");
	require(near(ProNovoConfig::getResidueMass("C"), 103.009185, 1e-6),
			"natural cysteine mass is wrong");
	require(ptmByToken("/").selectable &&
			near(ProNovoConfig::getResidueMass("/"), 57.021464, 1e-6),
			"variable IAA did not become available when fixed CAM was disabled");
	requireCounts(
		ProNovoConfig::configIsotopologue.mResidueSourcedComposition
			.at("(")[sipros::IsotopeSource::Biosynthetic],
		{0, -1, 0, 0, 0, 0},
		"S-nitrosylation did not remove natural-cysteine hydrogen");
	requireCounts(
		ProNovoConfig::configIsotopologue.mResidueSourcedComposition
			.at("(")[sipros::IsotopeSource::ReagentNatural],
		{0, 0, 1, 1, 0, 0},
		"S-nitrosylation without fixed CAM is wrong");

	sipros::SourcedComposition naturalSno;
	require(ProNovoConfig::configIsotopologue.computeSourcedComposition(
			"C(", naturalSno),
			"failed to compute natural-cysteine SNO composition");
	const double naturalSnoMass = lowestMassForSequence(
		ProNovoConfig::configIsotopologue, "C(");
	require(ProNovoConfig::configureVariablePtms({"/"}, -1, error), error);
	sipros::SourcedComposition variableCamCysteine;
	require(ProNovoConfig::configIsotopologue.computeSourcedComposition(
			"C/", variableCamCysteine),
			"failed to compute variable-CAM cysteine");
	requireSameComposition(variableCamCysteine, defaultFixedCysteine,
			"variable and fixed CAM cysteine formulas differ");
	require(!ProNovoConfig::configureFixedPtms({"default"}, error),
			"fixed CAM was enabled on top of variable IAA");
	require(ProNovoConfig::getEnabledFixedPtmNames().empty(),
			"failed fixed-PTM request was not atomic");

	require(ProNovoConfig::configureVariablePtms({"none"}, -1, error), error);
	require(ProNovoConfig::configureFixedPtms({"default"}, error), error);
	require(ProNovoConfig::getEnabledFixedPtmNames() ==
			std::vector<std::string>{"carbamidomethyl"},
			"fixed selector 'default' did not restore CAM");
	require(ProNovoConfig::getChemistryProfileId() ==
			"sipros5/source-aware-cam-tryptic-water/v1",
			"fixed-CAM chemistry profile ID was not restored");
	sipros::SourcedComposition fixedSno;
	require(ProNovoConfig::configIsotopologue.computeSourcedComposition(
			"C(", fixedSno),
			"failed to compute fixed-CAM SNO composition");
	requireSameComposition(fixedSno, naturalSno,
			"S-nitrosocysteine depends on starting fixed-CAM state");
	require(near(lowestMassForSequence(
				 ProNovoConfig::configIsotopologue, "C("),
				 naturalSnoMass, 1e-6),
			"S-nitrosocysteine mass depends on starting fixed-CAM state");
	require(!ptmByToken("/").selectable &&
			near(ProNovoConfig::getResidueMass("/"), 0.0, 1e-8),
			"variable IAA did not return to its fixed-CAM compatibility state");

	require(ProNovoConfig::configureChemistryProfileId(
			"sipros5/source-aware-natural-cys-tryptic-water/v1", error),
		error);
	require(ProNovoConfig::getEnabledFixedPtmNames().empty() &&
			near(ProNovoConfig::getResidueMass("C"), 103.009185, 1e-6),
			"natural-Cys library profile did not disable fixed CAM");
	require(ProNovoConfig::configureChemistryProfileId(
			"sipros5/source-aware-cam-tryptic-water/v1", error),
		error);
	require(ProNovoConfig::getEnabledFixedPtmNames() ==
			std::vector<std::string>{"carbamidomethyl"} &&
			near(ProNovoConfig::getResidueMass("C"), 160.030649, 1e-6),
			"fixed-CAM library profile did not restore fixed CAM");
	const auto beforeUnknownProfile =
		ProNovoConfig::getEnabledFixedPtmNames();
	require(!ProNovoConfig::configureChemistryProfileId(
			"sipros5/unknown-chemistry/v99", error) && !error.empty(),
			"unknown library chemistry profile unexpectedly succeeded");
	require(ProNovoConfig::getEnabledFixedPtmNames() == beforeUnknownProfile,
			"unknown library chemistry profile changed fixed PTMs");

	const double beforeInvalidMass = ProNovoConfig::getResidueMass("C");
	const auto beforeInvalidFixed = ProNovoConfig::getEnabledFixedPtmNames();
	require(!ProNovoConfig::configureFixedPtms({"not-a-fixed-ptm"}, error),
			"unknown fixed PTM selector unexpectedly succeeded");
	require(ProNovoConfig::getEnabledFixedPtmNames() == beforeInvalidFixed &&
			near(ProNovoConfig::getResidueMass("C"), beforeInvalidMass, 1e-10),
			"unknown fixed PTM selector was not atomic");
	require(!ProNovoConfig::configureFixedPtms(
			{"none", "default"}, error),
			"exclusive fixed selector 'none' accepted another selector");

	require(ProNovoConfig::configureFixedPtms({"none"}, error), error);
	require(ProNovoConfig::configureFixedPtms({}, error) &&
			ProNovoConfig::getEnabledFixedPtmNames().empty(),
			"empty fixed selector did not preserve the active selection");
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to reload fixed PTM defaults");
	require(ProNovoConfig::getEnabledFixedPtmNames() ==
			std::vector<std::string>{"carbamidomethyl"} &&
			near(ProNovoConfig::getResidueMass("C"), 160.030649, 1e-6),
			"profile reload did not restore fixed CAM");
}

void checkPsmPtmTranslation()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"failed to load chemistry for PSM PTM translation");
	struct TranslationCase
	{
		const char *plain;
		const char *modified;
		const char *expected;
	};
	const std::vector<TranslationCase> cases{
		{"M", "M[147]", "[M~]"},
		{"N", "N[115]", "[N!]"},
		{"Q", "Q[+0.984016]", "[Q!]"},
		{"S", "S[167]", "[S@]"},
		{"K", "K[42.010565]", "[K%]"},
		{"E", "E[14.015650]", "[E^]"},
		{"R", "R[28.031300]", "[R&]"},
		{"K", "K[42.046950]", "[K*]"},
		{"C", "C[28.990164]", "[C(]"},
		{"Y", "Y[44.985079]", "[Y)]"},
		{"D", "D[45.987721]", "[D$]"},
		{"C", "C[160]", "[C]"},
		{"PEP", "n[43]PEP", "[%PEP]"},
		{"PEPTIDE", "", "[PEPTIDE]"}};
	for (const TranslationCase &entry : cases)
	{
		std::string translated;
		std::string error;
		require(ProNovoConfig::translatePsmPeptide(
				entry.plain, entry.modified, translated, error),
			error);
		require(translated == entry.expected,
			std::string("wrong PSM PTM translation for ") + entry.modified +
				": " + translated);
	}

	std::string translated;
	std::string error;
	require(ProNovoConfig::translatePsmPeptide(
			"ACMNK", "n[42.0106]AC[160]M[147]N[115]K[42.04695]",
			translated, error),
		error);
	require(translated == "[%ACM~N!K*]",
			"combined PSM modifications were not converted before mass use");
	require(ProNovoConfig::translatePsmPeptide(
			"K", "K[170]", "1K(42.046950)",
			translated, error),
		error);
	require(translated == "[K*]",
			"Assigned Modifications did not disambiguate nominal K[170] as trimethylation");
	require(ProNovoConfig::translatePsmPeptide(
			"K", "K[170]", "1K(42.010565)",
			translated, error),
		error);
	require(translated == "[K%]",
			"Assigned Modifications did not disambiguate nominal K[170] as acetylation");
	require(ProNovoConfig::translatePsmPeptide(
			"K", "K[170]", "1K(28.031300),1K(14.015650)",
			translated, error),
		error);
	require(translated == "[K&^]",
			"summed nominal modification did not consume both positioned deltas");
	require(ProNovoConfig::translatePsmPeptide(
			"ACMK", "n[43]ACM[147]K[170]",
			"N-term(42.0106),2C(57.0215),3M(15.9949),4K(42.04695)",
			translated, error),
		error);
	require(translated == "[%ACM~K*]",
			"positioned Assigned Modifications were not reconciled with Modified Peptide");
	require(!ProNovoConfig::translatePsmPeptide(
			"A", "A[999]", translated, error) && !error.empty(),
			"unsupported PSM modification unexpectedly translated");
	require(!ProNovoConfig::translatePsmPeptide(
			"M", "M[147foo]", translated, error),
			"partially numeric PSM modification unexpectedly translated");

	require(ProNovoConfig::configureFixedPtms({"none"}, error), error);
	require(ProNovoConfig::translatePsmPeptide(
			"ACD", "ACD", "2C(57.0215)", translated, error),
		error);
	require(translated == "[AC/D]",
			"Assigned fixed CAM omitted from Modified Peptide was not preserved when fixed CAM was disabled");
	require(ProNovoConfig::translatePsmPeptide(
			"C", "C[160]", translated, error),
		error);
	require(translated == "[C/]",
			"FragPipe CAM did not become variable IAA with fixed CAM disabled");
	require(ProNovoConfig::translatePsmPeptide(
			"C", "C[132]", translated, error),
		error);
	require(translated == "[C(]",
			"absolute S-nitrosocysteine mass translated to the wrong token");
	require(!ProNovoConfig::translatePsmPeptide(
			"ACD", "ACD", "2M(15.9949)", translated, error) &&
			!error.empty(),
			"Assigned Modifications residue/position mismatch was accepted");
}

void checkSipProfile()
{
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to load SIP profile");
	require(ProNovoConfig::getSearchType() == "SIP" &&
				ProNovoConfig::getSearchName() == "SIP",
			"wrong SIP identity");
	require(ProNovoConfig::getMaxPTMcount() == 0,
			"wrong SIP maximum PTM count");
	std::map<std::string, std::string> ptms;
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.empty(),
			"SIP profile unexpectedly enables variable PTMs");
	checkCommonChemistry();

	require(ProNovoConfig::applySipAbundance('C', 0.5),
			"failed to apply 50% C13 abundance");
	require(near(ProNovoConfig::getDeductionCoefficient(), -0.005),
			"wrong SIP deduction coefficient at 50% abundance");

	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip),
			"failed to reset SIP profile for PTM selection");
	std::string error;
	require(ProNovoConfig::configureVariablePtms(
			{"oxidation"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>{{"~", "M"}} &&
			ProNovoConfig::getMaxPTMcount() == 3,
			"explicit SIP PTM did not automatically open a three-PTM maximum");
	require(ProNovoConfig::configureVariablePtms(
			{"deamidation"}, 1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) &&
			ptms == std::map<std::string, std::string>{{"!", "NQ"}} &&
			ProNovoConfig::getMaxPTMcount() == 1,
			"explicit SIP PTM maximum override was ignored");
	require(ProNovoConfig::configureVariablePtms({"none"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.empty() &&
			ProNovoConfig::getMaxPTMcount() == 0,
			"SIP selector 'none' did not restore a zero PTM maximum");
	require(ProNovoConfig::configureVariablePtms({"all"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.size() == 12 &&
			ProNovoConfig::getMaxPTMcount() == 3,
			"SIP selector 'all' did not open the selectable catalog");
	require(ProNovoConfig::configureVariablePtms({"default"}, -1, error), error);
	require(ProNovoConfig::getPTMinfo(ptms) && ptms.empty() &&
			ProNovoConfig::getMaxPTMcount() == 0,
			"SIP selector 'default' did not restore profile defaults");
	require(ProNovoConfig::configureVariablePtms({}, 2, error) &&
			ProNovoConfig::getPTMinfo(ptms) && ptms.empty() &&
			ProNovoConfig::getMaxPTMcount() == 2,
			"SIP maximum-only override changed the PTM selection");
	require(ProNovoConfig::load(ProNovoConfig::Profile::Sip) &&
			ProNovoConfig::getMaxPTMcount() == 0,
			"SIP reload did not reset the PTM maximum");
}

} // namespace

int main()
{
	try
	{
		checkNoUninitializedFallback();
		checkRegularProfile();
		checkUnifiedRegularProductIons();
		checkResiduePtmStateCache();
		checkPtmCatalogAndSelectors();
		checkFixedPtmConfiguration();
		checkPsmPtmTranslation();
		checkProductIonNetComposition();
		checkSipProfile();
		checkSourceIsolationAndEndpointMasses();
		checkSipMassCacheRefreshAndTargetSelection();
		std::cout << "Built-in configuration profiles passed\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "Built-in configuration profile test failed: "
				  << error.what() << '\n';
		return 1;
	}
}
