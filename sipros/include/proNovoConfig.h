#ifndef PRONOVOCONFIG_H_
#define PRONOVOCONFIG_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "omp.h"

namespace sipros
{

constexpr std::size_t ElementCount = 6;
constexpr std::size_t IsotopeSourceCount = 3;

enum class Element : std::size_t
{
	Carbon = 0,
	Hydrogen = 1,
	Oxygen = 2,
	Nitrogen = 3,
	Phosphorus = 4,
	Sulfur = 5
};

// Atom provenance is part of the compiled ProNovo chemistry. Only
// Biosynthetic atoms follow the selected SIP abundance. ReagentNatural and
// DigestionSolvent atoms retain the natural isotope distribution of the real
// element when isotope envelopes are convolved.
enum class IsotopeSource : std::size_t
{
	Biosynthetic = 0,
	ReagentNatural = 1,
	DigestionSolvent = 2
};

using AtomCounts = std::array<int, ElementCount>;

struct SourcedComposition
{
	std::array<AtomCounts, IsotopeSourceCount> atoms{};

	AtomCounts &operator[](IsotopeSource source)
	{
		return atoms[static_cast<std::size_t>(source)];
	}

	const AtomCounts &operator[](IsotopeSource source) const
	{
		return atoms[static_cast<std::size_t>(source)];
	}

	SourcedComposition &operator+=(const SourcedComposition &other)
	{
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				atoms[source][element] += other.atoms[source][element];
			}
		}
		return *this;
	}

	SourcedComposition &operator-=(const SourcedComposition &other)
	{
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				atoms[source][element] -= other.atoms[source][element];
			}
		}
		return *this;
	}

	AtomCounts total() const
	{
		AtomCounts result{};
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				result[element] += atoms[source][element];
			}
		}
		return result;
	}

	AtomCounts naturalSourceTotal() const
	{
		AtomCounts result{};
		for (std::size_t source =
				 static_cast<std::size_t>(IsotopeSource::ReagentNatural);
			 source < IsotopeSourceCount;
			 ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				result[element] += atoms[source][element];
			}
		}
		return result;
	}
};

inline SourcedComposition operator+(SourcedComposition left,
									const SourcedComposition &right)
{
	left += right;
	return left;
}

inline SourcedComposition operator-(SourcedComposition left,
									const SourcedComposition &right)
{
	left -= right;
	return left;
}

inline SourcedComposition compositionFrom(IsotopeSource source,
										const AtomCounts &counts)
{
	SourcedComposition composition;
	composition[source] = counts;
	return composition;
}

} // namespace sipros

class Isotopologue;

using namespace std;

typedef long long INT64;

// To keep time information of functions.
#define CLOCKSTART                        \
	INT64 mem_start = checkMemoryUsage(); \
	double begin = omp_get_wtime();       \
	cout << "  " << __FUNCTION__ << " started" << endl
#define CLOCKSTOP                                      \
	INT64 mem_end = checkMemoryUsage();                 \
	double end = omp_get_wtime();                       \
	cout << "  " << __FUNCTION__ << " finished in " << (end - begin) << "s"; \
	INT64 mem_delta = mem_end - mem_start;              \
	if (mem_delta != 0)                                 \
	{                                                   \
		cout << ", memory delta " << mem_delta << " MB"; \
	}                                                   \
	cout << endl

// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
	// get KB memory into count
	unsigned int count = 0;

#if defined(__linux__)
	ifstream f("/proc/self/status"); // read the linux file
	while (!f.eof())
	{
		string key;
		f >> key;
		if (key == "VmData:")
		{ // size of data
			f >> count;
			break;
		}
	}
	f.close();
#endif

	// return MBs memory (size of data)
	return (count / 1024);
};

//--------------Comet Begin------------
#define PROTON_MASS 1.00727646688
#define NUM_ION_SERIES 9
#define NUM_SP_IONS 200 // num ions for preliminary scoring

#define ION_SERIES_A 0
#define ION_SERIES_B 1
#define ION_SERIES_C 2
#define ION_SERIES_X 3
#define ION_SERIES_Y 4
#define ION_SERIES_Z 5

#define SPARSE_MATRIX_SIZE 100
#define FLOAT_ZERO 1e-6 // 0.000001

#define MAX_FRAGMENT_CHARGE 5
#define MAX_PEPTIDE_LEN 150 // max # of AA for a peptide

struct Options // output parameters
{
	int iRemovePrecursor; // 0=no, 1=yes, 2=ETD precursors
	double dMinIntensity;
	double dRemovePrecursorTol;

	Options()
	{
		iRemovePrecursor = 0; // 0=no, 1=yes, 2=ETD precursors
		dMinIntensity = 0;
		dRemovePrecursorTol = 0;
	}
};

struct IonInfo
{
	int iNumIonSeriesUsed;
	int piSelectedIonSeries[NUM_ION_SERIES];
	int bUseNeutralLoss;
	int iIonVal[NUM_ION_SERIES];
	IonInfo()
	{
		bUseNeutralLoss = 0;
		iIonVal[ION_SERIES_A] = 0;
		iIonVal[ION_SERIES_B] = 1;
		iIonVal[ION_SERIES_C] = 0;
		iIonVal[ION_SERIES_X] = 0;
		iIonVal[ION_SERIES_Y] = 1;
		iIonVal[ION_SERIES_Z] = 0;
		iNumIonSeriesUsed = 2;
		piSelectedIonSeries[0] = 1;
		piSelectedIonSeries[1] = 4;
	}
};

struct PrecalcMasses
{
	double dNtermProton;		 // dAddNterminusPeptide + PROTON_MASS
	double dCtermOH2Proton;		 // dAddCterminusPeptide + dOH2fragment + PROTON_MASS
	double dCtermOH2;			 // dAddCterminusPeptide + dOHfragment
	int iMinus17HighRes;		 // BIN'd value of mass(NH3)
	int iMinus17LowRes;
	int iMinus18HighRes; // BIN'd value of mass(H2O)
	int iMinus18LowRes;
	double dCO;
	double dNH3;
	double dNH2;
	double dCOminusH2;
};

#define AminoAcidMassesSize 256
// store the mass for different amino acids
class AminoAcidMasses
{
public:
	static double dNULL;
	static double dERROR;
	vector<double> vdMasses;
	// double vdMasses[AminoAcidMassesSize];

	// construct function
	AminoAcidMasses();
	// clear vdMasses
	void clear();
	// reach an empty spot
	double end();
	// return the mass for the given amino acid
	double find(char _cAminoAcid);

	double operator[](char _cAminoAcid) const;

	double &operator[](char _cAminoAcid);
};
//--------------Comet End------------
class ProNovoConfig
{
public:
	enum class Profile
	{
		Regular,
		Sip
	};

	struct PtmDefinition
	{
		string name;
		string token;
		string sites;
		string description;
		double externalMonoisotopicShift;
		bool selectable;
		bool regularDefault;
	};

	struct FixedPtmDefinition
	{
		string name;
		string sites;
		string description;
		double externalMonoisotopicShift;
		bool profileDefault;
	};

	// Initialize session-wide state from a compiled profile.
	static bool load(Profile profile);

	// Return every compiled PTM. The selectable flag reflects whether each token
	// is valid for the active fixed chemistry.
	static const vector<PtmDefinition> &getPtmCatalog();
	static const vector<FixedPtmDefinition> &getFixedPtmCatalog();
	static vector<string> getEnabledFixedPtmNames();
	// Convert FragPipe-style bracketed masses to the compiled one-character
	// chemistry tokens before any theoretical mass or spectrum calculation.
	static bool translatePsmPeptide(
		const string &plainPeptide,
		const string &modifiedPeptide,
		string &translatedPeptide,
		string &error);
	static bool translatePsmPeptide(
		const string &plainPeptide,
		const string &modifiedPeptide,
		const string &assignedModifications,
		string &translatedPeptide,
		string &error);
	static bool configureFixedPtms(
		const vector<string> &selectors,
		string &error);

	/*
	 * Configure the exact variable-PTM set after loading a profile. An empty
	 * selector list preserves the loaded profile defaults; -1 leaves the maximum
	 * PTM count unspecified. Selectors accept catalog names, one-character
	 * peptide tokens, and the special values default, none, and all.
	 */
	static bool configureVariablePtms(
		const vector<string> &selectors,
		int maxPtmCountOverride,
		string &error);

	static vector<pair<string, string>> getNeutralLossList()
	{
		return vpNeutralLossList;
	}

	static string getSearchName()
	{
		return sSearchName;
	}

	static void setSearchName(const string &searchName)
	{
		sSearchName = searchName;
	}

	static string getSearchType()
	{
		return sSearchType;
	}

	static string getChemistryProfileId();
	// Apply the exact compiled fixed-PTM state named by spectra-library
	// metadata. Unknown or obsolete profile IDs are rejected; there is no
	// compatibility fallback.
	static bool configureChemistryProfileId(
		const string &profileId,
		string &error);

	static char getSeparator();

	// Set the active FASTA database path.
	static void setFASTAfilename(const string &fastaFilename);
	// Return the active FASTA database path.
	static string getFASTAfilename()
	{
		return sFASTAFilename;
	}

	// retrieve the Minimum length of a peptide
	static int getMinPeptideLength()
	{
		return iMinPeptideLength;
	}

	// retrieve the max length of a peptide
	static int getMaxPeptideLength()
	{
		return iMaxPeptideLength;
	}

	// Return the parent-ion mass tolerance in Da.
	static double getMassAccuracyParentIon()
	{
		return dMassAccuracyParentIon;
	}

	// Return the fragment-ion mass tolerance in Da.
	static double getMassAccuracyFragmentIon()
	{
		return dMassAccuracyFragmentIon;
	}

	static void setMassAccuracy(double parentIonToleranceDa, double fragmentIonToleranceDa);

	// Expand a peptide mass into the compiled parent-mass windows.
	static bool getPeptideMassWindows(double dPeptideMass,
									  vector<pair<double, double>> &vpPeptideMassWindows);
	// Expand a peptide mass with its source-aware, composition-weighted
	// nominal-neutron spacing. FASTA search uses this overload exclusively.
	static bool getPeptideMassWindows(double dPeptideMass,
									  double precursorNeutronMass,
									  vector<pair<double, double>> &vpPeptideMassWindows);
	static const vector<int> &getParentMassWindows()
	{
		return viParentMassWindows;
	}

	// Return the active maximum number of variable PTMs per peptide.
	static int getMaxPTMcount()
	{
		return iMaxPTMcount;
	}

	// Return the active compiled cleavage rules.
	static string getCleavageAfterResidues()
	{
		return sCleavageAfterResidues;
	}
	static string getCleavageBeforeResidues()
	{
		return sCleavageBeforeResidues;
	}
	static int getMaxMissedCleavages()
	{
		return iMaxMissedCleavages;
	}
	static bool getTestStartRemoval()
	{
		return bTestStartRemoval;
	}

	static bool getPTMinfo(map<string, string> &mPTMinfo);

	static Isotopologue configIsotopologue;
	static vector<vector<double>> naturalAtomIsotopeProbabilities;
	static const vector<double> &getNaturalAtomIsotopeProbabilities(
		size_t atomIndex);
	static vector<string> vsSingleResidueNames;
	static vector<double> vdSingleResidueMasses;

	static double getResidueMass(string sResidue);

	static double getTerminusMassN()
	{
		return dTerminusMassN;
	}
	static double getTerminusMassC()
	{
		return dTerminusMassC;
	}

	static double getProtonMass()
	{
		return 1.007276466;
	}

	static double getNeutronMass()
	{
		return neutronMass;
	}

	static double pnorm(double dMean, double dStandardDeviation,
						double dRandomVariable)
	{
		double dZScore = (dRandomVariable - dMean) / dStandardDeviation;
		double dProbability = 0.5 * erfc(-dZScore / sqrt(2.0));
		return dProbability;
	}

	static double scoreError(double dMassError)
	{

		//	pnorm function
		return (1.0 - pnorm(0, (getMassAccuracyFragmentIon() / 2), fabs(dMassError))) * 2.0;

		//  sigmoid function
		//	return ( 1/(1+exp(dMassError*600-3)));
	}
	static string &getSetSIPelement() { return SIPelement; }
	// Compute the score deduction coefficient for the active SIP target.
	static void setDeductionCoefficient();
	static int atomIndex(char sipAtom);
	static bool validatePreparationChemistry(const Isotopologue &iso,
										 std::string &error);
	static bool refreshResidueDistributions(Isotopologue &iso);
	static int resolveSipIsotopeIndex(const Isotopologue &iso, char sipAtom, int isotopeMassNumber);
	static void setSipAbundance(Isotopologue &iso, char sipAtom, int isotopeIndex, double sipPct);
	static double getIsotopeAbundancePct(const Isotopologue &iso, char sipAtom, int isotopeIndex);
	// Select the isotope spacing used by precursor-window and library matching
	// without changing the current isotope abundances.
	static bool selectSipTarget(char sipAtom, int isotopeMassNumber,
								std::string &error);
	static bool applySipAbundance(char sipAtom, double fraction);
	// get deduction coefficient in score function
	static double getDeductionCoefficient() { return deductionCoefficient; }

	//---------------Comet Begin---------------------
	static Options options;
	static IonInfo ionInformation;
	static int iXcorrProcessingOffset;
	static PrecalcMasses precalcMasses;
	static double dMaxMS2ScanMass;
	static double dMaxPeptideMass;
	// static map<char, double> pdAAMassFragment;
	static AminoAcidMasses pdAAMassFragment;
	static double dHighResInverseBinWidth;
	static double dLowResInverseBinWidth;
	static double dHighResOneMinusBinOffset;
	static double dLowResOneMinusBinOffset;
	static int iMaxPercusorCharge;
	//---------------Comet End-----------------------

	//---------------Myrimatch Begin-----------------
	static double ClassSizeMultiplier;
	static int NumIntensityClasses;
	static int minIntensityClassCount;
	static double ticCutoffPercentage;
	static int MaxPeakCount;
	static int MinMatchedFragments;
	static double minObservedMz;
	static double maxObservedMz;
	//---------------Myrimatch End-------------------

	//---------------Sipros Score Begin--------------
	static int INTTOPKEEP; // the top n PSM for calculation of other two scores
	static int iRank;
	//---------------Sipros Score End----------------
	static string sCleavageAfterResidues;
	static string sCleavageBeforeResidues;

private:
	// Active search-session state.
	static string sFASTAFilename;
	static string sSearchType;
	static string sSearchName;

	static int iMaxPTMcount;

	static int iMinPeptideLength;
	static int iMaxPeptideLength;

	static int iMaxMissedCleavages;
	static bool bTestStartRemoval;

	static double dMassAccuracyParentIon;
	static double dMassAccuracyFragmentIon;
	static vector<int> viParentMassWindows;

	static vector<pair<double, double>> vpPeptideMassWindowOffset;

	static vector<pair<string, string>> vpNeutralLossList;

	static double dTerminusMassN;
	static double dTerminusMassC;

	static string SIPelement;
	// for deductionCoefficient compute in SIP search
	static double neutronMass, deductionCoefficient;

	static bool refreshSessionMassCaches();
	static bool calculatePeptideMassWindowOffset();
};

#endif /*PRONOVOCONFIG_H_*/
