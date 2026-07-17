#include "proNovoConfig.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "PeptideIsotopeCalculator.h"
#include "isotopologue.h"

namespace
{

constexpr const char *kCamChemistryProfileId =
	"sipros5/source-aware-cam-tryptic-water/v1";
constexpr const char *kNaturalCysChemistryProfileId =
	"sipros5/source-aware-natural-cys-tryptic-water/v1";

// The compiled profiles are the single source of truth for Sipros chemistry
// and search defaults. Runtime inputs such as FASTA paths, SIP targets,
// abundances, and optional tolerance overrides do not belong here.
const std::vector<ProNovoConfig::PtmDefinition> &compiledPtmCatalog()
{
	static const std::vector<ProNovoConfig::PtmDefinition> catalog{
		{"oxidation", "~", "M",
		 "Oxidation of methionine; oxygen remains reagent-natural.",
		 15.994915, true, true},
		{"deamidation", "!", "NQ",
		 "Deamidation of asparagine or glutamine; removed H/N are "
		 "biosynthetic and incorporated O is reagent-natural.",
		 0.984016, true, true},
		{"phosphorylation", "@", "STYHD",
		 "Phosphorylation (+HPO3) without forced fragment neutral loss.",
		 79.966332, true, false},
		{"phosphorylation-loss-hpo3", ">", "STYHD",
		 "Phosphorylation whose modified fragments lose HPO3.",
		 79.966332, true, false},
		{"phosphorylation-loss-hpo3-h2o", "<", "ST",
		 "Phosphorylation whose modified fragments lose HPO3 and H2O.",
		 79.966332, true, false},
		{"acetylation", "%", "K",
		 "Acetylation (+C2H2O).",
		 42.010565, true, false},
		{"mono-methylation", "^", "KRED",
		 "Mono-methylation (+CH2).",
		 14.015650, true, false},
		{"di-methylation", "&", "KR",
		 "Di-methylation (+C2H4).",
		 28.031300, true, false},
		{"tri-methylation", "*", "K",
		 "Tri-methylation (+C3H6).",
		 42.046950, true, false},
		{"s-nitrosylation", "(", "C",
		 "S-nitrosylation replaces any fixed CAM, removes biosynthetic H, "
		 "and adds reagent-natural NO.",
		 28.990164, true, false},
		{"nitration", ")", "Y",
		 "Tyrosine nitration removes biosynthetic H and adds reagent-natural NO2.",
		 44.985079, true, false},
		{"iaa-blocking", "/", "C",
		 "Variable CAM/IAA token. It is selectable only when fixed "
		 "carbamidomethylation is disabled; otherwise it is a zero-delta alias.",
		 57.021464, true, false},
		{"beta-methylthiolation", "$", "D",
		 "Beta-methylthiolation (+CH2S).",
		 45.987721, true, false}};
	return catalog;
}

const std::vector<ProNovoConfig::FixedPtmDefinition> &
compiledFixedPtmCatalog()
{
	static const std::vector<ProNovoConfig::FixedPtmDefinition> catalog{
		{"carbamidomethyl", "C",
		 "Carbamidomethylation of cysteine by reagent-natural IAA (+C2H3NO).",
		 57.021464, true}};
	return catalog;
}

std::string trimText(const std::string &text)
{
	const size_t first = text.find_first_not_of(" \t\f\v\n\r");
	if (first == std::string::npos)
		return std::string();
	const size_t last = text.find_last_not_of(" \t\f\v\n\r");
	return text.substr(first, last - first + 1);
}

bool parsePsmModificationMass(const std::string &text, double &mass)
{
	try
	{
		std::string normalized = trimText(text);
		if (!normalized.empty() && normalized.front() == '+')
			normalized.erase(normalized.begin());
		if (normalized.empty())
			return false;
		size_t consumed = 0;
		mass = std::stod(normalized, &consumed);
		return consumed == normalized.size() && std::isfinite(mass);
	}
	catch (const std::exception &)
	{
		return false;
	}
}

double reportedMassMatchError(double reportedMass,
							  double unmodifiedResidueMass,
							  double modificationShift)
{
	const double deltaError = std::abs(reportedMass - modificationShift);
	const double absoluteError = std::abs(
		reportedMass - (unmodifiedResidueMass + modificationShift));
	double best = -1.0;
	// Delta masses in FragPipe's Modified Peptide column retain decimals;
	// absolute residue masses are commonly rounded to an integer.
	if (deltaError <= 0.05)
		best = deltaError;
	if (absoluteError <= 0.5)
		best = best < 0.0 ? absoluteError : std::min(best, absoluteError);
	return best;
}

std::string normalizeSelector(const std::string &selector)
{
	const size_t first = selector.find_first_not_of(" \t\f\v\n\r");
	if (first == std::string::npos)
		return std::string();
	const size_t last = selector.find_last_not_of(" \t\f\v\n\r");
	std::string normalized;
	normalized.reserve(last - first + 1);
	bool previousDash = false;
	for (size_t index = first; index <= last; ++index)
	{
		const unsigned char raw = static_cast<unsigned char>(selector[index]);
		const bool separator = std::isspace(raw) != 0 || raw == '_';
		char value = separator ? '-' : static_cast<char>(std::tolower(raw));
		if (value == '-' && previousDash)
			continue;
		normalized.push_back(value);
		previousDash = value == '-';
	}
	return normalized;
}

std::string ptmSearchKey(const std::string &token)
{
	if (token == ">")
		return ">to1";
	if (token == "<")
		return "<to2";
	return token;
}

std::vector<std::pair<std::string, std::string>> neutralLossesFor(
	const std::map<std::string, std::string> &variablePtms)
{
	std::vector<std::pair<std::string, std::string>> losses;
	if (variablePtms.find(">to1") != variablePtms.end())
		losses.push_back({">", "1"});
	if (variablePtms.find("<to2") != variablePtms.end())
		losses.push_back({"<", "2"});
	return losses;
}

std::map<std::string, std::string> productIonPtmSites()
{
	std::map<std::string, std::string> sites;
	for (const ProNovoConfig::PtmDefinition &definition :
		 ProNovoConfig::getPtmCatalog())
	{
		if (definition.selectable)
			sites[definition.token] = definition.sites;
	}
	// Search-time neutral-loss processing replaces these phosphorylation
	// tokens while retaining the same residue specificity.
	sites["1"] = sites.at(">");
	sites["2"] = sites.at("<");
	return sites;
}

class BuiltInConfig
{
public:
	struct FixedModification
	{
		std::string name;
		std::string sites;
		sipros::SourcedComposition delta;
	};

	std::string searchType;
	std::string searchName;
	int maxPtmCount = 0;
	int defaultMaxPtmCount = 0;
	int minPeptideLength = 7;
	int maxPeptideLength = 60;
	std::string cleavageAfter = "KR";
	std::string cleavageBefore = "ACDEFGHIJKLMNPQRSTVWY";
	int maxMissedCleavages = 2;
	bool removeFirstMethionine = true;
	double parentToleranceDa = 0.01;
	double fragmentToleranceDa = 0.01;
	std::vector<int> parentMassWindows{-2, -1, 0, 1, 2};
	std::map<std::string, std::string> variablePtms;
	std::map<std::string, std::string> defaultVariablePtms;
	std::string atomNames = "CHONPS";
	std::map<std::string, sipros::SourcedComposition> baseResidues;
	std::map<std::string, sipros::SourcedComposition> residues;
	std::vector<FixedModification> fixedModifications;
	std::vector<std::string> defaultFixedPtms;
	std::vector<std::string> enabledFixedPtms;
	std::vector<IsotopeDistribution> atoms;
	double deductionMinValue = 0.0;
	double deductionFold = 0.0;

	bool fixedPtmEnabled(const std::string &name) const
	{
		return std::find(enabledFixedPtms.begin(), enabledFixedPtms.end(), name) !=
			enabledFixedPtms.end();
	}

	const FixedModification &fixedModification(
		const std::string &name) const
	{
		for (const FixedModification &modification : fixedModifications)
		{
			if (modification.name == name)
				return modification;
		}
		throw std::logic_error("Unknown compiled fixed modification.");
	}

	void rebuildResidues()
	{
		residues = baseResidues;
		for (const std::string &name : enabledFixedPtms)
		{
			const FixedModification &modification = fixedModification(name);
			for (char site : modification.sites)
			{
				auto residue = residues.find(std::string(1, site));
				if (residue == residues.end())
					throw std::logic_error(
						"Fixed modification site is absent from built-in residues.");
				residue->second += modification.delta;
			}
		}

		const sipros::SourcedComposition nitrosyl =
			sipros::compositionFrom(
				sipros::IsotopeSource::Biosynthetic,
				{0, -1, 0, 0, 0, 0}) +
			sipros::compositionFrom(
				sipros::IsotopeSource::ReagentNatural,
				{0, 0, 1, 1, 0, 0});
		if (fixedPtmEnabled("carbamidomethyl"))
		{
			const sipros::SourcedComposition &cam =
				fixedModification("carbamidomethyl").delta;
			// Fixed CAM is already present on C. IAA is mass-neutral, while SNO
			// first removes CAM.
			residues["/"] = {};
			residues["("] = nitrosyl - cam;
		}
		else
		{
			residues["/"] = fixedModification("carbamidomethyl").delta;
			residues["("] = nitrosyl;
		}
	}

	static BuiltInConfig create(ProNovoConfig::Profile profile)
	{
		BuiltInConfig config;

		// Six real elements in C,H,O,N,P,S order. Source provenance controls which atoms can follow SIP labeling.
		config.atoms = {
			IsotopeDistribution({12.000000, 13.003355},
								{0.9893, 0.0107}),
			IsotopeDistribution({1.007825, 2.014102},
								{0.999885, 0.000115}),
			IsotopeDistribution({15.994915, 16.999132, 17.999160},
								{0.99757, 0.00038, 0.00205}),
			IsotopeDistribution({14.003074, 15.000109},
								{0.99632, 0.00368}),
			IsotopeDistribution({30.973762}, {1.0}),
			IsotopeDistribution(
				{31.972071, 32.971459, 33.967867, 34.967867, 35.967081},
				{0.9493, 0.0076, 0.0429, 0.0, 0.0002})};

		const auto biosynthetic = [](const sipros::AtomCounts &counts)
		{
			return sipros::compositionFrom(
				sipros::IsotopeSource::Biosynthetic, counts);
		};
		const auto reagentNatural = [](const sipros::AtomCounts &counts)
		{
			return sipros::compositionFrom(
				sipros::IsotopeSource::ReagentNatural, counts);
		};
		const auto digestionSolvent = [](const sipros::AtomCounts &counts)
		{
			return sipros::compositionFrom(
				sipros::IsotopeSource::DigestionSolvent, counts);
		};

		config.baseResidues = {
			{"Nterm", digestionSolvent({0, 1, 0, 0, 0, 0})},
			{"Cterm", digestionSolvent({0, 1, 1, 0, 0, 0})},
			{"J", biosynthetic({6, 11, 1, 1, 0, 0})},
			{"I", biosynthetic({6, 11, 1, 1, 0, 0})},
			{"L", biosynthetic({6, 11, 1, 1, 0, 0})},
			{"A", biosynthetic({3, 5, 1, 1, 0, 0})},
			{"S", biosynthetic({3, 5, 2, 1, 0, 0})},
			{"G", biosynthetic({2, 3, 1, 1, 0, 0})},
			{"V", biosynthetic({5, 9, 1, 1, 0, 0})},
			{"E", biosynthetic({5, 7, 3, 1, 0, 0})},
			{"K", biosynthetic({6, 12, 1, 2, 0, 0})},
			{"T", biosynthetic({4, 7, 2, 1, 0, 0})},
			{"D", biosynthetic({4, 5, 3, 1, 0, 0})},
			{"R", biosynthetic({6, 12, 1, 4, 0, 0})},
			{"P", biosynthetic({5, 7, 1, 1, 0, 0})},
			{"N", biosynthetic({4, 6, 2, 2, 0, 0})},
			{"F", biosynthetic({9, 9, 1, 1, 0, 0})},
			{"Q", biosynthetic({5, 8, 2, 2, 0, 0})},
			{"Y", biosynthetic({9, 9, 2, 1, 0, 0})},
			{"M", biosynthetic({5, 9, 1, 1, 0, 1})},
			{"H", biosynthetic({6, 7, 1, 3, 0, 0})},
			{"C", biosynthetic({3, 5, 1, 1, 0, 1})},
			{"W", biosynthetic({11, 10, 1, 2, 0, 0})},
			// Oxidation adds oxygen from an exogenous reagent/air pool.
			{"~", reagentNatural({0, 0, 1, 0, 0, 0})},
			// Deamidation removes peptide H/N and incorporates natural O.
			{"!", biosynthetic({0, -1, 0, -1, 0, 0}) +
					  reagentNatural({0, 0, 1, 0, 0, 0})},
			{"@", biosynthetic({0, 1, 3, 0, 1, 0})},
			{">", biosynthetic({0, 1, 3, 0, 1, 0})},
			{"<", biosynthetic({0, 1, 3, 0, 1, 0})},
			// Chemistry-only fragment replacements for phospho neutral loss.
			{"1", biosynthetic({0, 0, 0, 0, 0, 0})},
			{"2", biosynthetic({0, -2, -1, 0, 0, 0})},
			{"%", biosynthetic({2, 2, 1, 0, 0, 0})},
			{"^", biosynthetic({1, 2, 0, 0, 0, 0})},
			{"&", biosynthetic({2, 4, 0, 0, 0, 0})},
			{"*", biosynthetic({3, 6, 0, 0, 0, 0})},
			{")", biosynthetic({0, -1, 0, 0, 0, 0}) +
				  reagentNatural({0, 0, 2, 1, 0, 0})},
			{"$", biosynthetic({1, 2, 0, 0, 0, 1})}};

		config.fixedModifications = {
			{"carbamidomethyl", "C",
			 reagentNatural({2, 3, 1, 1, 0, 0})}};
		config.defaultFixedPtms = {"carbamidomethyl"};
		config.enabledFixedPtms = config.defaultFixedPtms;
		config.rebuildResidues();

		switch (profile)
		{
		case ProNovoConfig::Profile::Regular:
			config.searchType = "Regular";
			config.searchName = "SE";
			config.defaultMaxPtmCount = 3;
			config.defaultVariablePtms = {{"!", "NQ"}, {"~", "M"}};
			break;
		case ProNovoConfig::Profile::Sip:
			config.searchType = "SIP";
			config.searchName = "SIP";
			config.defaultMaxPtmCount = 0;
			config.deductionMinValue = 0.005;
			config.deductionFold = 4.0;
			break;
		default:
			throw std::invalid_argument("Unknown built-in Sipros profile.");
		}
		config.maxPtmCount = config.defaultMaxPtmCount;
		config.variablePtms = config.defaultVariablePtms;
		return config;
	}
};

BuiltInConfig activeConfig;
bool configLoaded = false;

bool supportedSipTarget(char sipAtom,
						int &massNumber,
						int &isotopeIndex,
						int &nominalShift)
{
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	switch (atom)
	{
	case 'C':
		massNumber = 13;
		isotopeIndex = 1;
		nominalShift = 1;
		return true;
	case 'H':
		massNumber = 2;
		isotopeIndex = 1;
		nominalShift = 1;
		return true;
	case 'N':
		massNumber = 15;
		isotopeIndex = 1;
		nominalShift = 1;
		return true;
	case 'O':
		massNumber = 18;
		isotopeIndex = 2;
		nominalShift = 2;
		return true;
	case 'S':
		massNumber = 34;
		isotopeIndex = 2;
		nominalShift = 2;
		return true;
	default:
		return false;
	}
}

struct ChemistryBuildState
{
	Isotopologue isotopologue;
	std::vector<std::string> residueNames;
	std::vector<double> residueMasses;
	double terminusMassN = 0.0;
	double terminusMassC = 0.0;
	std::vector<std::vector<double>> naturalProbabilities;
};

bool buildChemistryState(const BuiltInConfig &config,
						 ChemistryBuildState &state,
						 std::string &error)
{
	if (!state.isotopologue.setupIsotopologue(
			config.residues, config.atoms, config.atomNames))
	{
		error = "Failed to build compiled residue isotope distributions.";
		return false;
	}
	if (!state.isotopologue.getSingleResidueMostAbundantMasses(
			state.residueNames,
			state.residueMasses,
			state.terminusMassN,
			state.terminusMassC))
	{
		error = "Failed to derive compiled residue masses.";
		return false;
	}
	state.naturalProbabilities.clear();
	for (const IsotopeDistribution &distribution :
		 state.isotopologue.vNaturalAtomIsotopicDistribution)
	{
		state.naturalProbabilities.push_back(distribution.vProb);
	}
	return true;
}

} // namespace

std::string ProNovoConfig::sFASTAFilename;
std::string ProNovoConfig::sSearchType;
std::string ProNovoConfig::sSearchName;
int ProNovoConfig::iMaxPTMcount = 0;
int ProNovoConfig::iMinPeptideLength = 0;
int ProNovoConfig::iMaxPeptideLength = 0;
std::string ProNovoConfig::sCleavageAfterResidues;
std::string ProNovoConfig::sCleavageBeforeResidues;
int ProNovoConfig::iMaxMissedCleavages = 0;
bool ProNovoConfig::bTestStartRemoval = false;
double ProNovoConfig::dMassAccuracyParentIon = 0.0;
double ProNovoConfig::dMassAccuracyFragmentIon = 0.0;
std::vector<int> ProNovoConfig::viParentMassWindows;
std::vector<std::pair<double, double>> ProNovoConfig::vpPeptideMassWindowOffset;
std::vector<std::pair<std::string, std::string>> ProNovoConfig::vpNeutralLossList;
std::vector<std::string> ProNovoConfig::vsSingleResidueNames;
std::vector<double> ProNovoConfig::vdSingleResidueMasses;
double ProNovoConfig::dTerminusMassN = 0.0;
double ProNovoConfig::dTerminusMassC = 0.0;

Isotopologue ProNovoConfig::configIsotopologue;
std::vector<std::vector<double>>
	ProNovoConfig::naturalAtomIsotopeProbabilities;

Options ProNovoConfig::options;
IonInfo ProNovoConfig::ionInformation;
int ProNovoConfig::iXcorrProcessingOffset = 75;
PrecalcMasses ProNovoConfig::precalcMasses;
double ProNovoConfig::dMaxMS2ScanMass = 0.0;
double ProNovoConfig::dMaxPeptideMass = 0.0;
AminoAcidMasses ProNovoConfig::pdAAMassFragment;
double ProNovoConfig::dHighResInverseBinWidth = 1.0 / 0.02;
double ProNovoConfig::dLowResInverseBinWidth = 1.0 / 1.0005;
double ProNovoConfig::dHighResOneMinusBinOffset = 1.0;
double ProNovoConfig::dLowResOneMinusBinOffset = 0.6;
int ProNovoConfig::iMaxPercusorCharge = 0;

double ProNovoConfig::ClassSizeMultiplier = 2.0;
int ProNovoConfig::NumIntensityClasses = 3;
int ProNovoConfig::minIntensityClassCount = static_cast<int>(
	(std::pow(ClassSizeMultiplier, NumIntensityClasses) - 1.0) /
	(ClassSizeMultiplier - 1.0));
double ProNovoConfig::ticCutoffPercentage = 0.98;
int ProNovoConfig::MaxPeakCount = 300;
int ProNovoConfig::MinMatchedFragments = 5;
double ProNovoConfig::minObservedMz = std::numeric_limits<double>::max();
double ProNovoConfig::maxObservedMz = 0.0;

int ProNovoConfig::INTTOPKEEP = 10;
int ProNovoConfig::iRank = 0;

std::string ProNovoConfig::SIPelement;
double ProNovoConfig::deductionCoefficient = 0.0;
double ProNovoConfig::neutronMass = 0.0;

double AminoAcidMasses::dNULL = -1.0;
double AminoAcidMasses::dERROR = -2.0;

AminoAcidMasses::AminoAcidMasses()
	: vdMasses(AminoAcidMassesSize, -1.0)
{
}

void AminoAcidMasses::clear()
{
	std::fill(vdMasses.begin(), vdMasses.end(), dNULL);
}

double AminoAcidMasses::end()
{
	return dNULL;
}

double AminoAcidMasses::find(char aminoAcid)
{
	return vdMasses.at(static_cast<unsigned char>(aminoAcid));
}

double AminoAcidMasses::operator[](char aminoAcid) const
{
	return vdMasses.at(static_cast<unsigned char>(aminoAcid));
}

double &AminoAcidMasses::operator[](char aminoAcid)
{
	return vdMasses.at(static_cast<unsigned char>(aminoAcid));
}

bool ProNovoConfig::refreshSessionMassCaches()
{
	std::vector<std::string> residueNames;
	std::vector<double> residueMasses;
	double terminusMassN = 0.0;
	double terminusMassC = 0.0;
	if (!configIsotopologue.getSingleResidueMostAbundantMasses(
			residueNames, residueMasses, terminusMassN, terminusMassC) ||
		residueNames.size() != residueMasses.size() ||
		activeConfig.atomNames.size() !=
			configIsotopologue.vAtomIsotopicDistribution.size())
	{
		return false;
	}

	AminoAcidMasses fragmentMasses;
	for (size_t atomIndexValue = 0;
		 atomIndexValue < activeConfig.atomNames.size(); ++atomIndexValue)
	{
		const IsotopeDistribution &distribution =
			configIsotopologue.vAtomIsotopicDistribution[atomIndexValue];
		if (distribution.vMass.empty() ||
			distribution.vMass.size() != distribution.vProb.size())
		{
			return false;
		}
		const auto mostAbundant = std::max_element(
			distribution.vProb.begin(), distribution.vProb.end());
		const size_t isotopeIndex = static_cast<size_t>(
			std::distance(distribution.vProb.begin(), mostAbundant));
		const char atom = static_cast<char>(std::tolower(
			static_cast<unsigned char>(activeConfig.atomNames[atomIndexValue])));
		if (fragmentMasses.find(atom) != fragmentMasses.end())
			return false;
		fragmentMasses[atom] = distribution.vMass[isotopeIndex];
	}
	for (size_t index = 0; index < residueNames.size(); ++index)
	{
		if (residueNames[index].size() != 1)
			return false;
		const char residue = residueNames[index][0];
		if (fragmentMasses.find(residue) != fragmentMasses.end())
			return false;
		fragmentMasses[residue] = residueMasses[index];
	}

	const double hydrogen = fragmentMasses.find('h');
	const double carbon = fragmentMasses.find('c');
	const double oxygen = fragmentMasses.find('o');
	const double nitrogen = fragmentMasses.find('n');
	if (hydrogen == fragmentMasses.end() ||
		carbon == fragmentMasses.end() ||
		oxygen == fragmentMasses.end() ||
		nitrogen == fragmentMasses.end())
	{
		return false;
	}

	PrecalcMasses calculated{};
	const double ammonia = hydrogen * 3.0 + nitrogen;
	const double water = hydrogen * 2.0 + oxygen;
	calculated.iMinus17LowRes = static_cast<int>(
		ammonia * dLowResInverseBinWidth + dLowResOneMinusBinOffset);
	calculated.iMinus17HighRes = static_cast<int>(
		ammonia * dHighResInverseBinWidth + dHighResOneMinusBinOffset);
	calculated.iMinus18LowRes = static_cast<int>(
		water * dLowResInverseBinWidth + dLowResOneMinusBinOffset);
	calculated.iMinus18HighRes = static_cast<int>(
		water * dHighResInverseBinWidth + dHighResOneMinusBinOffset);
	calculated.dNtermProton = PROTON_MASS;
	// Tryptic terminal water is introduced by digestion solvent and therefore
	// stays at its natural isotope distribution even under SIP enrichment.
	calculated.dCtermOH2 = terminusMassN + terminusMassC;
	calculated.dCtermOH2Proton =
		calculated.dCtermOH2 + PROTON_MASS;
	calculated.dCO = oxygen + carbon;
	calculated.dNH2 = nitrogen + hydrogen * 2.0;
	calculated.dNH3 = ammonia;
	calculated.dCOminusH2 = calculated.dCO - hydrogen * 2.0;

	vsSingleResidueNames = std::move(residueNames);
	vdSingleResidueMasses = std::move(residueMasses);
	dTerminusMassN = terminusMassN;
	dTerminusMassC = terminusMassC;
	pdAAMassFragment = std::move(fragmentMasses);
	precalcMasses = calculated;
	return true;
}

bool ProNovoConfig::load(Profile profile)
{
	configLoaded = false;
	activeConfig = BuiltInConfig::create(profile);

	sFASTAFilename.clear();
	sSearchType = activeConfig.searchType;
	sSearchName = activeConfig.searchName;
	iMaxPTMcount = activeConfig.maxPtmCount;
	iMinPeptideLength = activeConfig.minPeptideLength;
	iMaxPeptideLength = activeConfig.maxPeptideLength;
	sCleavageAfterResidues = activeConfig.cleavageAfter;
	sCleavageBeforeResidues = activeConfig.cleavageBefore;
	iMaxMissedCleavages = activeConfig.maxMissedCleavages;
	bTestStartRemoval = activeConfig.removeFirstMethionine;
	dMassAccuracyParentIon = activeConfig.parentToleranceDa;
	dMassAccuracyFragmentIon = activeConfig.fragmentToleranceDa;
	viParentMassWindows = activeConfig.parentMassWindows;
	SIPelement = "C";
	deductionCoefficient = 0.0;

	minObservedMz = std::numeric_limits<double>::max();
	maxObservedMz = 0.0;
	dMaxMS2ScanMass = 0.0;
	dMaxPeptideMass = 0.0;
	iMaxPercusorCharge = 0;
	options = Options{};
	ionInformation = IonInfo{};

	if (!configIsotopologue.setupIsotopologue(
			activeConfig.residues,
			activeConfig.atoms,
			activeConfig.atomNames))
	{
		return false;
	}
	if (!refreshResidueDistributions(configIsotopologue))
		return false;
	if (!refreshSessionMassCaches())
	{
		return false;
	}

	naturalAtomIsotopeProbabilities.clear();
	for (const IsotopeDistribution &distribution :
		 configIsotopologue.vNaturalAtomIsotopicDistribution)
	{
		naturalAtomIsotopeProbabilities.push_back(distribution.vProb);
	}
	vpNeutralLossList = neutralLossesFor(activeConfig.variablePtms);
	configLoaded = true;
	std::string chemistryError;
	if (!validatePreparationChemistry(configIsotopologue, chemistryError))
	{
		std::cerr << "ERROR: " << chemistryError << std::endl;
		configLoaded = false;
		return false;
	}
	setDeductionCoefficient();
	if (!calculatePeptideMassWindowOffset())
	{
		configLoaded = false;
		return false;
	}
	return true;
}

void ProNovoConfig::setFASTAfilename(const std::string &fastaFilename)
{
	sFASTAFilename = fastaFilename;
}

void ProNovoConfig::setMassAccuracy(double parentIonToleranceDa,
									 double fragmentIonToleranceDa)
{
	dMassAccuracyParentIon = parentIonToleranceDa;
	dMassAccuracyFragmentIon = fragmentIonToleranceDa;
	calculatePeptideMassWindowOffset();
}

const std::vector<ProNovoConfig::PtmDefinition> &
ProNovoConfig::getPtmCatalog()
{
	static const std::vector<PtmDefinition> fixedCamView = []
	{
		std::vector<PtmDefinition> view = compiledPtmCatalog();
		for (PtmDefinition &definition : view)
		{
			if (definition.token == "/")
				definition.selectable = false;
		}
		return view;
	}();
	static const std::vector<PtmDefinition> naturalCysteineView = []
	{
		std::vector<PtmDefinition> view = compiledPtmCatalog();
		for (PtmDefinition &definition : view)
		{
			if (definition.token == "/")
				definition.selectable = true;
		}
		return view;
	}();
	const bool camEnabled =
		!configLoaded || activeConfig.fixedPtmEnabled("carbamidomethyl");
	return camEnabled ? fixedCamView : naturalCysteineView;
}

const std::vector<ProNovoConfig::FixedPtmDefinition> &
ProNovoConfig::getFixedPtmCatalog()
{
	return compiledFixedPtmCatalog();
}

std::vector<std::string> ProNovoConfig::getEnabledFixedPtmNames()
{
	return configLoaded ? activeConfig.enabledFixedPtms
						: std::vector<std::string>{};
}

bool ProNovoConfig::translatePsmPeptide(
	const std::string &plainPeptide,
	const std::string &modifiedPeptide,
	std::string &translatedPeptide,
	std::string &error)
{
	return translatePsmPeptide(
		plainPeptide, modifiedPeptide, std::string(),
		translatedPeptide, error);
}

bool ProNovoConfig::translatePsmPeptide(
	const std::string &plainPeptide,
	const std::string &modifiedPeptide,
	const std::string &assignedModifications,
	std::string &translatedPeptide,
	std::string &error)
{
	translatedPeptide.clear();
	error.clear();
	if (!configLoaded)
	{
		error = "Built-in Sipros chemistry is not initialized.";
		return false;
	}

	const std::string modified = trimText(modifiedPeptide);
	const std::string source = modified.empty()
								   ? trimText(plainPeptide)
								   : modified;
	if (source.empty())
	{
		error = "empty peptide";
		return false;
	}

	const auto fixedEnabled = [&](const std::string &name)
	{
		return activeConfig.fixedPtmEnabled(name);
	};
	const auto unmodifiedResidueMass = [&](char residue)
	{
		double mass = getResidueMass(std::string(1, residue));
		for (const FixedPtmDefinition &definition : compiledFixedPtmCatalog())
		{
			if (fixedEnabled(definition.name) &&
				definition.sites.find(residue) != std::string::npos)
			{
				mass -= definition.externalMonoisotopicShift;
			}
		}
		return mass;
	};
	const auto variableTokenForFixed = [&](
		const FixedPtmDefinition &fixed) -> std::string
	{
		for (const PtmDefinition &definition : compiledPtmCatalog())
		{
			if (definition.sites == fixed.sites &&
				std::abs(definition.externalMonoisotopicShift -
						 fixed.externalMonoisotopicShift) <= 1e-6)
			{
				return definition.token;
			}
		}
		return std::string();
	};
	const auto translateResidueModification = [&](
		char residue,
		const std::string &massText,
		std::string &body) -> bool
	{
		double reportedMass = 0.0;
		if (!parsePsmModificationMass(massText, reportedMass))
		{
			error = "invalid modification mass " + massText;
			return false;
		}

		const double baseMass = unmodifiedResidueMass(residue);
		for (const FixedPtmDefinition &definition : compiledFixedPtmCatalog())
		{
			if (definition.sites.find(residue) == std::string::npos)
			{
				continue;
			}
			const double matchError = reportedMassMatchError(
				reportedMass, baseMass,
				definition.externalMonoisotopicShift);
			if (matchError < 0.0)
				continue;
			if (fixedEnabled(definition.name))
				return true;
			const std::string token = variableTokenForFixed(definition);
			if (token.empty())
			{
				error = "fixed modification " + definition.name +
						" has no compiled variable token";
				return false;
			}
			body += token;
			return true;
		}

		const PtmDefinition *best = nullptr;
		double bestError = 0.0;
		for (const PtmDefinition &definition : compiledPtmCatalog())
		{
			if (definition.sites.find(residue) == std::string::npos)
				continue;
			const double matchError = reportedMassMatchError(
				reportedMass, baseMass,
				definition.externalMonoisotopicShift);
			if (matchError >= 0.0 &&
				(best == nullptr || matchError < bestError))
			{
				best = &definition;
				bestError = matchError;
			}
		}
		if (best == nullptr)
		{
			error = std::string("unsupported modification ") + residue +
					"[" + massText + "]";
			return false;
		}
		body += best->token;
		return true;
	};
	const auto translateNTermModification = [&](
		const std::string &massText,
		std::string &body) -> bool
	{
		double reportedMass = 0.0;
		if (!parsePsmModificationMass(massText, reportedMass))
		{
			error = "invalid N-terminal modification mass " + massText;
			return false;
		}
		if (std::abs(reportedMass - 42.010565) > 0.05 &&
			std::abs(reportedMass - 43.0) > 0.5)
		{
			error = "unsupported N-terminal modification n[" + massText + "]";
			return false;
		}
		for (const PtmDefinition &definition : compiledPtmCatalog())
		{
			if (definition.name == "acetylation")
			{
				body += definition.token;
				return true;
			}
		}
		error = "N-terminal acetylation is absent from compiled chemistry";
		return false;
	};

	struct ParsedPeptide
	{
		std::vector<char> residues;
		std::vector<std::string> nTermModifications;
		std::vector<std::vector<std::string>> residueModifications;
	};
	ParsedPeptide parsed;
	for (size_t index = 0; index < source.size();)
	{
		const char raw = source[index];
		if (index == 0 && raw == 'n' && index + 1 < source.size() &&
			source[index + 1] == '[')
		{
			const size_t right = source.find(']', index + 2);
			if (right == std::string::npos)
			{
				error = "unterminated N-terminal modification";
				return false;
			}
			parsed.nTermModifications.push_back(
				source.substr(index + 2, right - index - 2));
			index = right + 1;
			continue;
		}
		if (std::isspace(static_cast<unsigned char>(raw)) != 0)
		{
			++index;
			continue;
		}
		if (!std::isalpha(static_cast<unsigned char>(raw)))
		{
			error = "unsupported peptide character";
			return false;
		}

		const char residue = static_cast<char>(
			std::toupper(static_cast<unsigned char>(raw)));
		parsed.residues.push_back(residue);
		parsed.residueModifications.emplace_back();
		++index;
		while (index < source.size() && source[index] == '[')
		{
			const size_t right = source.find(']', index + 1);
			if (right == std::string::npos)
			{
				error = "unterminated modification";
				return false;
			}
			parsed.residueModifications.back().push_back(
				source.substr(index + 1, right - index - 1));
			index = right + 1;
		}
	}

	if (parsed.residues.empty())
	{
		error = "empty peptide";
		return false;
	}

	std::vector<std::string> assignedNTermModifications;
	std::vector<std::vector<std::string>> assignedResidueModifications(
		parsed.residues.size());
	const std::string assignedText = trimText(assignedModifications);
	for (size_t begin = 0; begin < assignedText.size();)
	{
		const size_t comma = assignedText.find(',', begin);
		const size_t end = comma == std::string::npos
						   ? assignedText.size()
						   : comma;
		const std::string entry = trimText(
			assignedText.substr(begin, end - begin));
		if (entry.empty())
		{
			error = "empty Assigned Modifications entry";
			return false;
		}

		const std::string nTermPrefix = "N-term(";
		const std::string lowerNTermPrefix = "n-term(";
		if ((entry.compare(0, nTermPrefix.size(), nTermPrefix) == 0 ||
			 entry.compare(0, lowerNTermPrefix.size(), lowerNTermPrefix) == 0) &&
			entry.back() == ')')
		{
			const std::string massText = entry.substr(
				nTermPrefix.size(),
				entry.size() - nTermPrefix.size() - 1);
			double mass = 0.0;
			if (!parsePsmModificationMass(massText, mass))
			{
				error = "invalid Assigned Modifications entry: " + entry;
				return false;
			}
			assignedNTermModifications.push_back(massText);
		}
		else
		{
			size_t positionEnd = 0;
			while (positionEnd < entry.size() &&
				   std::isdigit(static_cast<unsigned char>(entry[positionEnd])) != 0)
			{
				++positionEnd;
			}
			if (positionEnd == 0 || positionEnd + 3 > entry.size() ||
				std::isalpha(static_cast<unsigned char>(entry[positionEnd])) == 0 ||
				entry[positionEnd + 1] != '(' || entry.back() != ')')
			{
				error = "invalid Assigned Modifications entry: " + entry;
				return false;
			}

			size_t position = 0;
			try
			{
				position = static_cast<size_t>(
					std::stoull(entry.substr(0, positionEnd)));
			}
			catch (const std::exception &)
			{
				error = "invalid Assigned Modifications position: " + entry;
				return false;
			}
			const char residue = static_cast<char>(std::toupper(
				static_cast<unsigned char>(entry[positionEnd])));
			if (position == 0 || position > parsed.residues.size())
			{
				error = "Assigned Modifications position is outside the peptide: " +
						entry;
				return false;
			}
			if (parsed.residues[position - 1] != residue)
			{
				error = "Assigned Modifications residue does not match peptide at "
						"position " +
						std::to_string(position) + ": " + entry;
				return false;
			}
			const std::string massText = entry.substr(
				positionEnd + 2,
				entry.size() - positionEnd - 3);
			double mass = 0.0;
			if (!parsePsmModificationMass(massText, mass))
			{
				error = "invalid Assigned Modifications entry: " + entry;
				return false;
			}
			assignedResidueModifications[position - 1].push_back(massText);
		}

		if (comma == std::string::npos)
			break;
		begin = comma + 1;
	}

	std::string body;
	std::vector<bool> assignedNTermUsed(
		assignedNTermModifications.size(), false);
	for (const std::string &modifiedMass : parsed.nTermModifications)
	{
		size_t assignedIndex = assignedNTermModifications.size();
		for (size_t index = 0;
			 index < assignedNTermModifications.size(); ++index)
		{
			if (!assignedNTermUsed[index])
			{
				assignedIndex = index;
				break;
			}
		}
		if (assignedIndex < assignedNTermModifications.size())
		{
			if (!translateNTermModification(
					assignedNTermModifications[assignedIndex], body))
				return false;
			assignedNTermUsed[assignedIndex] = true;
		}
		else if (!translateNTermModification(modifiedMass, body))
		{
			return false;
		}
	}
	for (size_t index = 0;
		 index < assignedNTermModifications.size(); ++index)
	{
		if (!assignedNTermUsed[index] &&
			!translateNTermModification(
				assignedNTermModifications[index], body))
		{
			return false;
		}
	}

	for (size_t residueIndex = 0;
		 residueIndex < parsed.residues.size(); ++residueIndex)
	{
		const char residue = parsed.residues[residueIndex];
		body.push_back(residue);
		const auto &assignedAtResidue =
			assignedResidueModifications[residueIndex];
		std::vector<bool> assignedUsed(assignedAtResidue.size(), false);
		for (const std::string &modifiedMassText :
			 parsed.residueModifications[residueIndex])
		{
			double modifiedMass = 0.0;
			if (!parsePsmModificationMass(modifiedMassText, modifiedMass))
			{
				error = "invalid modification mass " + modifiedMassText;
				return false;
			}
			const double baseMass = unmodifiedResidueMass(residue);
			std::vector<std::pair<size_t, double>> availableAssigned;
			for (size_t assignedIndex = 0;
				 assignedIndex < assignedAtResidue.size(); ++assignedIndex)
			{
				if (assignedUsed[assignedIndex])
					continue;
				double assignedMass = 0.0;
				if (!parsePsmModificationMass(
						assignedAtResidue[assignedIndex], assignedMass))
					continue;
				availableAssigned.push_back({assignedIndex, assignedMass});
			}
			if (availableAssigned.size() > 16)
			{
				error = "too many Assigned Modifications at one peptide residue";
				return false;
			}

			// FragPipe may collapse several positioned deltas into one nominal
			// bracketed mass. Match every nonempty subset so fixed CAM plus a
			// variable PTM is reconciled once rather than duplicated or rejected.
			std::vector<size_t> currentAssigned;
			std::vector<size_t> bestAssigned;
			double bestError = std::numeric_limits<double>::infinity();
			const auto considerSubsets = [&](auto &&self,
										 size_t offset,
										 double shiftSum) -> void
			{
				if (offset == availableAssigned.size())
				{
					if (currentAssigned.empty())
						return;
					const double matchError = reportedMassMatchError(
						modifiedMass, baseMass, shiftSum);
					if (matchError >= 0.0 &&
						(matchError + 1e-12 < bestError ||
						 (std::abs(matchError - bestError) <= 1e-12 &&
						  currentAssigned.size() > bestAssigned.size())))
					{
						bestError = matchError;
						bestAssigned = currentAssigned;
					}
					return;
				}
				self(self, offset + 1, shiftSum);
				currentAssigned.push_back(
					availableAssigned[offset].first);
				self(self, offset + 1,
					 shiftSum + availableAssigned[offset].second);
				currentAssigned.pop_back();
			};
			considerSubsets(considerSubsets, 0, 0.0);

			if (!bestAssigned.empty())
			{
				for (size_t assignedIndex : bestAssigned)
				{
					if (!translateResidueModification(
							residue, assignedAtResidue[assignedIndex], body))
						return false;
					assignedUsed[assignedIndex] = true;
				}
			}
			else if (!translateResidueModification(
					 residue, modifiedMassText, body))
			{
				return false;
			}
		}
		for (size_t assignedIndex = 0;
			 assignedIndex < assignedAtResidue.size(); ++assignedIndex)
		{
			if (!assignedUsed[assignedIndex] &&
				!translateResidueModification(
					residue, assignedAtResidue[assignedIndex], body))
			{
				return false;
			}
		}
	}

	translatedPeptide = "[" + body + "]";
	return true;
}

bool ProNovoConfig::configureFixedPtms(
	const std::vector<std::string> &selectors,
	std::string &error)
{
	error.clear();
	if (!configLoaded)
	{
		error = "Built-in Sipros configuration is not initialized.";
		return false;
	}
	if (selectors.empty())
		return true;

	std::vector<std::string> normalizedSelectors;
	normalizedSelectors.reserve(selectors.size());
	for (const std::string &selector : selectors)
		normalizedSelectors.push_back(normalizeSelector(selector));
	if (std::find(normalizedSelectors.begin(), normalizedSelectors.end(),
				  "none") != normalizedSelectors.end() &&
		normalizedSelectors.size() != 1)
	{
		error = "Fixed PTM selector 'none' must be used alone.";
		return false;
	}

	std::vector<std::string> candidateNames;
	auto addName = [&](const std::string &name)
	{
		if (std::find(candidateNames.begin(), candidateNames.end(), name) ==
			candidateNames.end())
		{
			candidateNames.push_back(name);
		}
	};
	for (const std::string &selector : normalizedSelectors)
	{
		if (selector == "none")
			continue;
		if (selector == "default")
		{
			for (const std::string &name : activeConfig.defaultFixedPtms)
				addName(name);
			continue;
		}
		if (selector == "all")
		{
			for (const FixedPtmDefinition &definition : compiledFixedPtmCatalog())
				addName(definition.name);
			continue;
		}

		const FixedPtmDefinition *matched = nullptr;
		for (const FixedPtmDefinition &definition : compiledFixedPtmCatalog())
		{
			if (normalizeSelector(definition.name) == selector)
			{
				matched = &definition;
				break;
			}
		}
		if (matched == nullptr)
		{
			error = "Unknown fixed PTM selector: " + selector;
			return false;
		}
		addName(matched->name);
	}

	const bool candidateCamEnabled =
		std::find(candidateNames.begin(), candidateNames.end(),
				  "carbamidomethyl") != candidateNames.end();
	if (candidateCamEnabled &&
		activeConfig.variablePtms.find("/") != activeConfig.variablePtms.end())
	{
		error = "Cannot enable fixed carbamidomethylation while the variable "
				"iaa-blocking PTM is enabled.";
		return false;
	}

	BuiltInConfig candidate = activeConfig;
	candidate.enabledFixedPtms = candidateNames;
	candidate.rebuildResidues();
	ChemistryBuildState chemistry;
	if (!buildChemistryState(candidate, chemistry, error))
		return false;
	if (!validatePreparationChemistry(chemistry.isotopologue, error))
		return false;

	activeConfig = candidate;
	// Isotopologue has immutable numeric policy members and is intentionally
	// non-assignable. Rebuild the session object after the candidate has already
	// been validated in isolation.
	if (!configIsotopologue.setupIsotopologue(
			activeConfig.residues,
			activeConfig.atoms,
			activeConfig.atomNames))
	{
		error = "Failed to install validated fixed-PTM chemistry.";
		return false;
	}
	if (!refreshResidueDistributions(configIsotopologue))
	{
		error = "Failed to build residue-PTM isotope distributions.";
		return false;
	}
	naturalAtomIsotopeProbabilities = chemistry.naturalProbabilities;
	if (!refreshSessionMassCaches())
	{
		error = "Failed to refresh residue masses for fixed-PTM chemistry.";
		return false;
	}
	setDeductionCoefficient();
	return calculatePeptideMassWindowOffset();
}

bool ProNovoConfig::configureVariablePtms(
	const std::vector<std::string> &selectors,
	int maxPtmCountOverride,
	std::string &error)
{
	error.clear();
	if (!configLoaded)
	{
		error = "Built-in Sipros configuration is not initialized.";
		return false;
	}
	if (maxPtmCountOverride < -1)
	{
		error = "Maximum PTM count must be -1 (unspecified) or nonnegative.";
		return false;
	}

	if (selectors.empty())
	{
		if (maxPtmCountOverride >= 0)
		{
			activeConfig.maxPtmCount = maxPtmCountOverride;
			iMaxPTMcount = maxPtmCountOverride;
		}
		return true;
	}

	std::vector<std::string> normalizedSelectors;
	normalizedSelectors.reserve(selectors.size());
	for (const std::string &selector : selectors)
		normalizedSelectors.push_back(normalizeSelector(selector));
	if (std::find(normalizedSelectors.begin(), normalizedSelectors.end(),
				  "none") != normalizedSelectors.end() &&
		normalizedSelectors.size() != 1)
	{
		error = "Variable PTM selector 'none' must be used alone.";
		return false;
	}

	std::map<std::string, std::string> candidatePtms;
	auto addDefinition = [&](const PtmDefinition &definition) -> bool
	{
		if (!definition.selectable)
		{
			error = "PTM '" + definition.name + "' cannot be selected with "
					"the active fixed chemistry.";
			return false;
		}
		candidatePtms[ptmSearchKey(definition.token)] = definition.sites;
		return true;
	};
	const std::vector<PtmDefinition> activeCatalog = getPtmCatalog();
	for (const std::string &selector : normalizedSelectors)
	{
		if (selector == "none")
			continue;
		if (selector == "default")
		{
			candidatePtms.insert(activeConfig.defaultVariablePtms.begin(),
							 activeConfig.defaultVariablePtms.end());
			continue;
		}
		if (selector == "all")
		{
			for (const PtmDefinition &definition : activeCatalog)
			{
				if (definition.selectable && !addDefinition(definition))
					return false;
			}
			continue;
		}

		const PtmDefinition *matched = nullptr;
		for (const PtmDefinition &definition : activeCatalog)
		{
			if (normalizeSelector(definition.name) == selector ||
				normalizeSelector(definition.token) == selector)
			{
				matched = &definition;
				break;
			}
		}
		if (matched == nullptr)
		{
			error = "Unknown variable PTM selector: " + selector;
			return false;
		}
		if (!addDefinition(*matched))
			return false;
	}

	int candidateMaxPtmCount = activeConfig.defaultMaxPtmCount;
	if (activeConfig.searchType == "SIP" && !candidatePtms.empty())
		candidateMaxPtmCount = 3;
	if (maxPtmCountOverride >= 0)
		candidateMaxPtmCount = maxPtmCountOverride;

	activeConfig.variablePtms = candidatePtms;
	activeConfig.maxPtmCount = candidateMaxPtmCount;
	iMaxPTMcount = candidateMaxPtmCount;
	vpNeutralLossList = neutralLossesFor(candidatePtms);
	return true;
}

char ProNovoConfig::getSeparator()
{
#if _WIN32
	return '\\';
#else
	return '/';
#endif
}

bool ProNovoConfig::getPTMinfo(std::map<std::string, std::string> &ptms)
{
	if (!configLoaded)
		return false;
	ptms = activeConfig.variablePtms;
	return true;
}

std::string ProNovoConfig::getChemistryProfileId()
{
	if (!configLoaded)
		throw std::logic_error(
			"Built-in Sipros configuration is not initialized.");
	if (activeConfig.fixedPtmEnabled("carbamidomethyl"))
	{
		return kCamChemistryProfileId;
	}
	return kNaturalCysChemistryProfileId;
}

bool ProNovoConfig::configureChemistryProfileId(
	const std::string &profileId,
	std::string &error)
{
	error.clear();
	if (!configLoaded)
	{
		error = "Built-in Sipros configuration is not initialized.";
		return false;
	}

	std::vector<std::string> fixedPtmSelectors;
	if (profileId == kCamChemistryProfileId)
	{
		fixedPtmSelectors = {"carbamidomethyl"};
	}
	else if (profileId == kNaturalCysChemistryProfileId)
	{
		fixedPtmSelectors = {"none"};
	}
	else
	{
		error = "Unknown chemistry_profile_id '" + profileId +
			"'. Regenerate the spectra library with the current Sipros build.";
		return false;
	}

	if (!configureFixedPtms(fixedPtmSelectors, error))
		return false;
	if (getChemistryProfileId() != profileId)
	{
		error = "Compiled chemistry did not reproduce chemistry_profile_id '" +
			profileId + "'.";
		return false;
	}
	return true;
}

const std::vector<double> &
ProNovoConfig::getNaturalAtomIsotopeProbabilities(size_t atomIndex)
{
	if (!configLoaded || atomIndex >= naturalAtomIsotopeProbabilities.size())
	{
		throw std::logic_error(
			"Built-in Sipros configuration is not initialized for this atom.");
	}
	return naturalAtomIsotopeProbabilities[atomIndex];
}

int ProNovoConfig::atomIndex(char sipAtom)
{
	if (!configLoaded)
		return -1;
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	const size_t position = activeConfig.atomNames.find(atom);
	return position == std::string::npos ? -1 : static_cast<int>(position);
}

bool ProNovoConfig::validatePreparationChemistry(
	const Isotopologue &iso,
	std::string &error)
{
	if (!configLoaded)
	{
		error = "Built-in Sipros configuration is not initialized.";
		return false;
	}
	const auto cysteine = iso.mResidueSourcedComposition.find("C");
	const auto nTerm = iso.mResidueSourcedComposition.find("Nterm");
	const auto cTerm = iso.mResidueSourcedComposition.find("Cterm");
	const auto variableIaa = iso.mResidueSourcedComposition.find("/");
	const auto sNitrosylation = iso.mResidueSourcedComposition.find("(");
	if (cysteine == iso.mResidueSourcedComposition.end() ||
		nTerm == iso.mResidueSourcedComposition.end() ||
		cTerm == iso.mResidueSourcedComposition.end() ||
		variableIaa == iso.mResidueSourcedComposition.end() ||
		sNitrosylation == iso.mResidueSourcedComposition.end())
	{
		error = "Built-in preparation chemistry lacks cysteine, peptide termini, "
				"or conditional cysteine PTMs.";
		return false;
	}
	const sipros::AtomCounts expectedCysBio{3, 5, 1, 1, 0, 1};
	const sipros::AtomCounts expectedCam{2, 3, 1, 1, 0, 0};
	const sipros::AtomCounts expectedNitrosylBio{0, -1, 0, 0, 0, 0};
	const sipros::AtomCounts expectedNitrosylReagent{0, 0, 1, 1, 0, 0};
	const sipros::AtomCounts noAtoms{};
	const sipros::AtomCounts expectedNTerm{0, 1, 0, 0, 0, 0};
	const sipros::AtomCounts expectedCTerm{0, 1, 1, 0, 0, 0};
	const sipros::AtomCounts cysteineReagent =
		cysteine->second[sipros::IsotopeSource::ReagentNatural];
	const bool fixedCamEnabled = cysteineReagent == expectedCam;
	if (cysteine->second[sipros::IsotopeSource::Biosynthetic] != expectedCysBio ||
		(cysteineReagent != expectedCam && cysteineReagent != noAtoms) ||
		cysteine->second[sipros::IsotopeSource::DigestionSolvent] !=
			noAtoms)
	{
		error = "Built-in preparation chemistry must model cysteine as "
				"biosynthetic C3H5NOS with either zero or one reagent-natural "
				"CAM C2H3NO group.";
		return false;
	}
	const sipros::AtomCounts expectedVariableIaa =
		fixedCamEnabled ? noAtoms : expectedCam;
	sipros::AtomCounts expectedSnoReagent = expectedNitrosylReagent;
	if (fixedCamEnabled)
	{
		for (size_t element = 0; element < expectedSnoReagent.size(); ++element)
			expectedSnoReagent[element] -= expectedCam[element];
	}
	if (variableIaa->second[sipros::IsotopeSource::Biosynthetic] != noAtoms ||
		variableIaa->second[sipros::IsotopeSource::ReagentNatural] !=
			expectedVariableIaa ||
		variableIaa->second[sipros::IsotopeSource::DigestionSolvent] != noAtoms ||
		sNitrosylation->second[sipros::IsotopeSource::Biosynthetic] !=
			expectedNitrosylBio ||
		sNitrosylation->second[sipros::IsotopeSource::ReagentNatural] !=
			expectedSnoReagent ||
		sNitrosylation->second[sipros::IsotopeSource::DigestionSolvent] != noAtoms)
	{
		error = "Conditional variable IAA and S-nitrosylation formulas do not "
				"match the active fixed-CAM state.";
		return false;
	}
	if (nTerm->second[sipros::IsotopeSource::DigestionSolvent] !=
			expectedNTerm ||
		cTerm->second[sipros::IsotopeSource::DigestionSolvent] !=
			expectedCTerm ||
		nTerm->second[sipros::IsotopeSource::Biosynthetic] != noAtoms ||
		nTerm->second[sipros::IsotopeSource::ReagentNatural] != noAtoms ||
		cTerm->second[sipros::IsotopeSource::Biosynthetic] != noAtoms ||
		cTerm->second[sipros::IsotopeSource::ReagentNatural] != noAtoms)
	{
		error = "Built-in preparation chemistry must source peptide-terminal "
				"H2O exclusively from natural-abundance digestion solvent.";
		return false;
	}
	const int phosphorusIndex = atomIndex('P');
	if (phosphorusIndex < 0 ||
		static_cast<size_t>(phosphorusIndex) >=
			iso.vNaturalAtomIsotopicDistribution.size() ||
		iso.vNaturalAtomIsotopicDistribution[static_cast<size_t>(phosphorusIndex)].vMass !=
			std::vector<double>{30.973762} ||
		iso.vNaturalAtomIsotopicDistribution[static_cast<size_t>(phosphorusIndex)].vProb !=
			std::vector<double>{1.0})
	{
		error = "Built-in preparation chemistry must use real monoisotopic "
				"phosphorus; pseudo-carbon element slots are forbidden.";
		return false;
	}
	return true;
}

bool ProNovoConfig::refreshResidueDistributions(Isotopologue &iso)
{
	return iso.refreshResidueStateCache(
		productIonPtmSites(), {"%"}, {});
}

int ProNovoConfig::resolveSipIsotopeIndex(const Isotopologue &iso,
											char sipAtom,
											int isotopeMassNumber)
{
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	int expectedMassNumber = 0;
	int expectedIsotopeIndex = -1;
	int nominalShift = 0;
	if (!supportedSipTarget(atom, expectedMassNumber,
							expectedIsotopeIndex, nominalShift) ||
		(isotopeMassNumber > 0 && isotopeMassNumber != expectedMassNumber))
	{
		throw std::runtime_error(
			"Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
	}

	const int index = atomIndex(atom);
	if (index < 0 ||
		index >= static_cast<int>(iso.vAtomIsotopicDistribution.size()))
	{
		throw std::runtime_error(
			"Built-in SIP atom distribution is unavailable.");
	}
	const auto &distribution =
		iso.vAtomIsotopicDistribution[static_cast<size_t>(index)];
	if (expectedIsotopeIndex >= static_cast<int>(distribution.vProb.size()) ||
		expectedIsotopeIndex >= static_cast<int>(distribution.vMass.size()) ||
		static_cast<int>(std::lround(
			distribution.vMass[static_cast<size_t>(expectedIsotopeIndex)])) !=
			expectedMassNumber)
	{
		throw std::runtime_error(
			"Built-in distribution does not contain the requested SIP isotope.");
	}
	return expectedIsotopeIndex;
}

void ProNovoConfig::setSipAbundance(Isotopologue &iso,
									char sipAtom,
									int isotopeIndex,
									double sipPct)
{
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	int massNumber = 0;
	int expectedIsotopeIndex = -1;
	int nominalShift = 0;
	if (!supportedSipTarget(atom, massNumber,
							expectedIsotopeIndex, nominalShift) ||
		isotopeIndex != expectedIsotopeIndex)
	{
		throw std::runtime_error(
			"Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
	}
	const int index = atomIndex(atom);
	if (index < 0 ||
		index >= static_cast<int>(iso.vAtomIsotopicDistribution.size()) ||
		isotopeIndex >= static_cast<int>(
			iso.vAtomIsotopicDistribution[static_cast<size_t>(index)].vProb.size()))
	{
		throw std::runtime_error(
			"Built-in isotope distribution for SIP target is unusable.");
	}

	auto &probabilities =
		iso.vAtomIsotopicDistribution[static_cast<size_t>(index)].vProb;
	if (!PeptideIsotopeCalculator::changeAtomProbability(
			probabilities, atom, sipPct / 100.0) ||
		!refreshResidueDistributions(iso))
	{
		throw std::runtime_error(
			"Failed to apply SIP abundance to built-in chemistry.");
	}
}

double ProNovoConfig::getIsotopeAbundancePct(
	const Isotopologue &iso,
	char sipAtom,
	int isotopeIndex)
{
	int massNumber = 0;
	int expectedIsotopeIndex = -1;
	int nominalShift = 0;
	if (!supportedSipTarget(sipAtom, massNumber,
							expectedIsotopeIndex, nominalShift) ||
		isotopeIndex != expectedIsotopeIndex)
	{
		throw std::runtime_error(
			"Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
	}
	const int index = atomIndex(sipAtom);
	if (index < 0 ||
		index >= static_cast<int>(iso.vAtomIsotopicDistribution.size()) ||
		isotopeIndex >= static_cast<int>(
			iso.vAtomIsotopicDistribution[static_cast<size_t>(index)].vProb.size()))
	{
		throw std::runtime_error(
			"Built-in isotope abundance is unavailable.");
	}
	return iso.vAtomIsotopicDistribution[static_cast<size_t>(index)]
			   .vProb[static_cast<size_t>(isotopeIndex)] *
		   100.0;
}

bool ProNovoConfig::selectSipTarget(char sipAtom,
									int isotopeMassNumber,
									std::string &error)
{
	error.clear();
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	int expectedMassNumber = 0;
	int isotopeIndex = -1;
	int nominalShift = 0;
	if (!configLoaded ||
		!supportedSipTarget(atom, expectedMassNumber,
							isotopeIndex, nominalShift) ||
		isotopeMassNumber != expectedMassNumber)
	{
		error = "Unsupported SIP isotope. Use C13,H2,N15,O18,S34.";
		return false;
	}
	try
	{
		resolveSipIsotopeIndex(
			configIsotopologue, atom, isotopeMassNumber);
		SIPelement = std::string(1, atom);
		setDeductionCoefficient();
		if (!calculatePeptideMassWindowOffset())
		{
			error = "Cannot calculate precursor windows for the SIP isotope.";
			return false;
		}
	}
	catch (const std::exception &exception)
	{
		error = exception.what();
		return false;
	}
	return true;
}

bool ProNovoConfig::applySipAbundance(char sipAtom, double fraction)
{
	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(sipAtom)));
	int massNumber = 0;
	int isotopeIndex = -1;
	int nominalShift = 0;
	if (!configLoaded ||
		!supportedSipTarget(atom, massNumber, isotopeIndex, nominalShift))
	{
		return false;
	}
	const int index = atomIndex(atom);
	if (index < 0 ||
		index >= static_cast<int>(configIsotopologue.vAtomIsotopicDistribution.size()) ||
		naturalAtomIsotopeProbabilities.size() !=
			configIsotopologue.vAtomIsotopicDistribution.size())
	{
		return false;
	}

	for (size_t atomIndexValue = 0;
		 atomIndexValue < naturalAtomIsotopeProbabilities.size();
		 ++atomIndexValue)
	{
		configIsotopologue.vAtomIsotopicDistribution[atomIndexValue].vProb =
			naturalAtomIsotopeProbabilities[atomIndexValue];
	}
	SIPelement = std::string(1, atom);
	if (!PeptideIsotopeCalculator::changeAtomProbability(
			configIsotopologue.vAtomIsotopicDistribution[static_cast<size_t>(index)].vProb,
			atom,
			fraction) ||
		!refreshResidueDistributions(configIsotopologue))
	{
		return false;
	}
	if (!refreshSessionMassCaches())
	{
		return false;
	}
	setDeductionCoefficient();
	return calculatePeptideMassWindowOffset();
}

double ProNovoConfig::getResidueMass(std::string residue)
{
	if (residue == "|||")
		return 0.0;
	for (size_t index = 0; index < vsSingleResidueNames.size(); ++index)
	{
		if (vsSingleResidueNames[index] == residue)
			return vdSingleResidueMasses[index];
	}
	std::cerr << "ERROR: cannot find residue " << residue << std::endl;
	return 0.0;
}

bool ProNovoConfig::calculatePeptideMassWindowOffset()
{
	vpPeptideMassWindowOffset.clear();
	std::sort(viParentMassWindows.begin(), viParentMassWindows.end());
	double lastLower = 0.0;
	double lastUpper = 0.0;
	bool haveWindow = false;
	for (int massWindow : viParentMassWindows)
	{
		const double currentLower =
			massWindow * getNeutronMass() - dMassAccuracyParentIon;
		const double currentUpper =
			massWindow * getNeutronMass() + dMassAccuracyParentIon;
		if (!haveWindow)
		{
			lastLower = currentLower;
			lastUpper = currentUpper;
			haveWindow = true;
		}
		else if (currentLower <= lastUpper)
		{
			lastUpper = currentUpper;
		}
		else
		{
			vpPeptideMassWindowOffset.push_back({lastLower, lastUpper});
			lastLower = currentLower;
			lastUpper = currentUpper;
		}
	}
	if (haveWindow)
		vpPeptideMassWindowOffset.push_back({lastLower, lastUpper});
	return haveWindow;
}

bool ProNovoConfig::getPeptideMassWindows(
	double peptideMass,
	std::vector<std::pair<double, double>> &peptideMassWindows)
{
	for (const auto &offset : vpPeptideMassWindowOffset)
	{
		peptideMassWindows.push_back(
			{peptideMass + offset.first, peptideMass + offset.second});
	}
	return !vpPeptideMassWindowOffset.empty();
}

bool ProNovoConfig::getPeptideMassWindows(
	double peptideMass,
	double precursorNeutronMass,
	std::vector<std::pair<double, double>> &peptideMassWindows)
{
	if (!(precursorNeutronMass > 0.0) ||
		!std::isfinite(precursorNeutronMass))
		return false;

	double lastLower = 0.0;
	double lastUpper = 0.0;
	bool haveWindow = false;
	for (int massWindow : viParentMassWindows)
	{
		const double currentLower = peptideMass +
			massWindow * precursorNeutronMass - dMassAccuracyParentIon;
		const double currentUpper = peptideMass +
			massWindow * precursorNeutronMass + dMassAccuracyParentIon;
		if (!haveWindow)
		{
			lastLower = currentLower;
			lastUpper = currentUpper;
			haveWindow = true;
		}
		else if (currentLower <= lastUpper)
		{
			lastUpper = currentUpper;
		}
		else
		{
			peptideMassWindows.push_back({lastLower, lastUpper});
			lastLower = currentLower;
			lastUpper = currentUpper;
		}
	}
	if (haveWindow)
		peptideMassWindows.push_back({lastLower, lastUpper});
	return haveWindow;
}

void ProNovoConfig::setDeductionCoefficient()
{
	if (!configLoaded)
		throw std::logic_error("Built-in Sipros configuration is not initialized.");
	if (SIPelement.size() != 1)
		throw std::runtime_error("Unsupported SIP target.");

	const char atom = static_cast<char>(
		std::toupper(static_cast<unsigned char>(SIPelement[0])));
	int massNumber = 0;
	int isotopeIndex = -1;
	int nominalShift = 0;
	if (!supportedSipTarget(atom, massNumber, isotopeIndex, nominalShift))
	{
		throw std::runtime_error(
			"Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
	}
	const int index = atomIndex(atom);
	if (index < 0 ||
		index >= static_cast<int>(configIsotopologue.vAtomIsotopicDistribution.size()))
	{
		throw std::runtime_error("Built-in SIP target distribution is unavailable.");
	}
	const auto &distribution =
		configIsotopologue.vAtomIsotopicDistribution[static_cast<size_t>(index)];
	if (isotopeIndex >= static_cast<int>(distribution.vMass.size()) ||
		isotopeIndex >= static_cast<int>(distribution.vProb.size()))
	{
		throw std::runtime_error("Built-in SIP target distribution is incomplete.");
	}
	SIPelement = std::string(1, atom);
	neutronMass =
		(distribution.vMass[static_cast<size_t>(isotopeIndex)] -
		 distribution.vMass[0]) /
		static_cast<double>(nominalShift);
	deductionCoefficient =
		-(activeConfig.deductionMinValue +
		  activeConfig.deductionFold *
			  std::pow(
				  distribution.vProb[static_cast<size_t>(isotopeIndex)] - 0.5,
				  8));
}
