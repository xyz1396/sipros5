#ifndef MS2SCAN_H
#define MS2SCAN_H

#include <vector>
#include <string>
#include <array>

#include "peptide.h"
#include "MVH.h"
#include "CometSearchMod.h"
// #include "TandemMassSpectrum.h"

#define BIN_RES 1000
#define LOW_BIN_RES 10
#define TOP_N 50
#define TOP_N_SIP 10
#define TOPPICKNUM 20
#define ZERO 0.00000001
#define SMALLINCREMENT 0.00000001 // to solve the incorrect cut off

using namespace std;

class ProductIon // for each ion
{
	char cIonType;				// y or b
	int iIonNumber;				// 1: y1 or b1; 2: y2 or b2; ...
	int iCharge;				// charge state
	int iMostAbundantPeakIndex; // related to the item with highest intensity
	double dMostAbundantMass;	// related to the item with highest intensity
	double dMostAbundantMZ;		// related to the item with highest intensity
	double dMZError;
	double dMassError;					 // calculate based on dMZError
	double dScoreWeight;				 // related to mass error and intensity
	bool bComplementaryFragmentObserved; // if y_k and b_n-k co-exist, where n is # of residues
public:
	ProductIon();
	~ProductIon();
	// set the basic info
	void setProductIon(char cIonTypeInput, int iIonNumberInput, int iChargeInput);
	// set observed info
	void setObservedInfo(double dMZErrorInput, double dWeightInput, double dMostAbundantMZInput,
						 int iMostAbundantPeakIndexInput);
	void setComplementaryFragmentObserved(bool bComplementaryFragmentObservedInput);
	char getIonType()
	{
		return cIonType;
	};
	int getIonNumber()
	{
		return iIonNumber;
	};
	int getCharge()
	{
		return iCharge;
	};
	double getMZError()
	{
		return dMZError;
	};
	double getMassError()
	{
		return dMassError;
	};
	double getScoreWeight()
	{
		return dScoreWeight;
	};
	double getMostAbundantMass()
	{
		return dMostAbundantMass;
	};
	double getMostAbundantMZ()
	{
		return dMostAbundantMZ;
	};
	int getMostAbundantPeakIndex()
	{
		return iMostAbundantPeakIndex;
	};
	bool getComplementaryFragmentObserved()
	{
		return bComplementaryFragmentObserved;
	};
};

/*
 class ScanUnit
 {
 public:
 double intensity;
 bool match;
 ScanUnit(double inten, bool mat){intensity = inten; match=mat;};
 };*/

class PeptideUnit
{
public:
	static constexpr int DdaResidualTopN = 5;

	// the mass of the matched precursor peak
	double dMeasuredParentMass;
	// the charge of the matched precursor peak
	int iMeasuredParentCharge;
	double dCalculatedParentMass;
	double dScore;
	string sIdentifiedPeptide;
	string sOriginalPeptide;
	string sProteinNames;
	string sScoringFunction;
	char cIdentifyPrefix;
	char cIdentifySuffix;
	char cOriginalPrefix;
	char cOriginalSuffix;

	// Sipros Ensemble
	string sPeptideForScoring;
	// WDPScore, XcorrScore, MVHScore
	std::array<double, 3> vdScores;
	// ranks of the 3 scores
	std::array<double, 3> vdRank;
	int iDdaResidualRank = 0;
	double dDdaResidualScore = 0.0;
	double dPepNeutralMass;
	double dPrecursorNeutronMass;
	double iPepLength;
	static int iNumScores;

	// SIP
	vector<vector<double>> vvdYionMass;
	vector<vector<double>> vvdYionProb;
	vector<vector<double>> vvdBionMass;
	vector<vector<double>> vvdBionProb;

	// void setPeptideUnitInfo(const Peptide *currentPeptide, const double &dScore, string sScoringFunction);
	// scoreIX=0: WDP; scoreIX=1: Xcorr; scoreIX=2: MVH
	void setPeptideUnitInfo(const tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple, 
	const double &dScore, string sScoringFunction, const int scoreIX);
	void setIonMassProb(const Peptide *currentPeptide);
	// get topN peaks in each isotopic envolope
	void setIonMassProb(const Peptide *currentPeptide, int topN);
};

class PeakList
{
public:
	int iLowestMass;
	// excluded
	int iHighestMass;
	vector<short> pMassHub;
	vector<double> pPeaks;
	vector<char> pClasses;
	int iPeakSize;
	int iMassHubSize;
	int iMassHubPairSizeMinusOne;
	static char iNULL;

	PeakList(map<double, char> *_peakData);
	~PeakList();

	char end() const
	{
		return iNULL;
	}

	char findNear(double mz, double tolerance) const
	{
		if (iPeakSize == 0)
		{
			return iNULL;
		}

		const int mzUpper = static_cast<int>(mz + tolerance);
		const int mzLower = static_cast<int>(mz - tolerance);
		if (mzUpper < iLowestMass || mzLower > iHighestMass)
		{
			return iNULL;
		}

		int binBegin = 0;
		int binEnd = iMassHubPairSizeMinusOne;
		if (mzLower >= iLowestMass)
		{
			binBegin = mzLower - iLowestMass;
		}
		if (mzUpper <= iHighestMass)
		{
			binEnd = mzUpper - iLowestMass;
		}

		double minimumDifference = 1000000.0;
		char intensityClass = iNULL;
		const short *massHub = pMassHub.data();
		const double *peaks = pPeaks.data();
		const char *classes = pClasses.data();
		for (int bin = binBegin; bin <= binEnd; ++bin)
		{
			const int peakBegin = massHub[bin * 2];
			if (peakBegin == -1)
			{
				continue;
			}
			const int peakEnd = massHub[bin * 2 + 1];
			for (int peak = peakBegin; peak < peakEnd; ++peak)
			{
				const double difference = std::fabs(mz - peaks[peak]);
				if (difference < minimumDifference)
				{
					minimumDifference = difference;
					intensityClass = classes[peak];
				}
			}
		}
		return minimumDifference < tolerance ? intensityClass : iNULL;
	}

	int size() const
	{
		return iPeakSize;
	}
};

class MS2Scan
{
public:
	int bin_res, iMaxMZ, iMinMZ;

	double dMassTolerance; // fragment-ion mass tolerance
	double dProtonMass;	   // proton mass

	string sScanType; // format: FT-MS1/FT-MS2@CID

	vector<double> vdpreprocessedMZ;
	vector<double> vdpreprocessedIntensity;
	vector<int> vipreprocessedCharge; // the value of vipreprocessedCharge is zero for low-resolution MS2Scan

	// vdMaxMzIntensity[i] is the maximum intensity at M/Z window of i plus and minus iMzRange
	vector<double> vdMaxMzIntensity;
	vector<double> vdMzIntensity;
	vector<double> vdHighIntensity; // thresholds

	vector<int> vbPeakPresenceBins;
	vector<pair<int, int>> vbPeakPresenceBins2D; // first is lower bounder, second is upper bound
	vector<int> viIntensityRank;				 // ranks of preprocessed intensities staring with 0

	void preprocessLowMS2();
	void preprocessHighMS2();
	void initialPreprocess();
	void sortPeakList(); // bubble sort the peak list by MZ in ascending order
	int getMaxValueIndex(const vector<double> &vdData);
	void normalizeMS2scan();
	void setIntensityThreshold();
	void filterMS2scan();
	static bool mygreater(double i, double j);
	void binCalculation();
	void binCalculation2D(); // replace binCalculation();
	// scoreIX=0: WDP; scoreIX=1: Xcorr; scoreIX=2: MVH
	void saveScore(const double &dScore,
				   const tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple,
				   vector<PeptideUnit *> &vpTopPeptides, string sScoreFunction, const int scoreIX);
	// scoreIX=0: WDP; scoreIX=1: Xcorr; scoreIX=2: MVH
	void saveScoreSIP(const double &dScore,
					  const tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple,
					  vector<PeptideUnit *> &vpTopPeptides,
					  string sScoreFunction, const int scoreIX);
	static bool GreaterScore(PeptideUnit *p1, PeptideUnit *p2);
	void WeightCompare(const string &sPeptide, vector<bool> &vbFragmentZ2);
	bool searchMZ(const double &dTarget, int &iIndex4Found);   // corresponds to binCalculation()
	bool searchMZ2D(const double &dTarget, int &iIndex4Found); // corresponds to binCalculation2D()
	bool searchMZ2D(const double &dTarget, const double &dErrRange, int &iIndex4Found);
	// merge same peptide. If no same peptide return false, otherwise, return true
	bool mergePeptide(vector<PeptideUnit *> &vpTopPeptides, const string &sPeptide, const string &sProteinName);
	void sortPreprocessedIntensity(); // sort preprocessed Intensity;

	void cleanup();

	// static bool mySUGreater (ScanUnit s1, ScanUnit s2);

	bool binarySearch(const double &dTarget, const vector<double> &vdList, const double &dTolerance,
					  vector<int> &viIndex4Found);

	// build the map between y or b ions for observed intensity and related mass (only for one ion)
	bool findProductIon(const vector<double> &vdIonMass, // expected mass
						const vector<double> &vdIonProb,
						// expected intensity based on the summation of all related intensities
						const int &iCharge, double &dScoreWeight, double &dAverageMZError, double &dMostAbundantObservedMZ, // with highest intensity
						int &iMostAbundantPeakIndex);																		// start with 0

	bool findProductIonSIP(const vector<double> &vdIonMass, // expected mass
						   const vector<double> &vdIonProb,
						   // expected intensity based on the summation of all related intensities
						   const int &iCharge, double &dScoreWeight, double &dAverageMZError, double &dMostAbundantObservedMZ, // with highest intensity
						   int &iMostAbundantPeakIndex);																	   // start with 0
	MS2Scan();
	//    MS2Scan(const MS2Scan *& cMS2Scan);
	~MS2Scan();

	int iParentChargeState;	   // Parent ion charge state
	double dParentMZ;		   // Parent ion M/Z
	// The reaction charge belongs to the isolation-window center.  When precursor
	// hypotheses come from neighboring MS1 peaks, it must remain metadata only;
	// scoring then uses the charge attached to each individual MS1 peak instead.
	bool bUseReactionChargeForScoring;
	int iMaxCandidateCharge;
	double dParentNeutralMass; // Parent neutral mass
	// Parent mass with charge
	double dParentMass;
	int iScanId;													   // Product ions in the scan
	int iParentScanID; 												   // Parent scan ID for DIA
	vector<int> iParentChargeStates;								   // Parent ion charge states for DIA
	vector<double> dParentMZs;										   // Parent ion M/Z for DIA
	// Keep the acquisition window independent from detected MS1 peaks.  Regular
	// search can then discover a weak co-isolated precursor from its fragments
	// before targeted MS1 feature extraction validates it.
	vector<pair<double, double>> vIsolationWindowsMz; // center m/z, full width
	vector<tuple<double, int, Peptide *>> vMassChargePeptidePtrTuples; // current set of peptides to be scored

	vector<PeptideUnit *> vpWeightSumTopPeptides;

	// the scores of all scored peptides
	vector<double> vdWeightSumAllScores;

	// int inumberofWeightSumScore;
	// double dsumofWeightScore;
	// double dsumofSquareWeightSumScore;

	vector<double> vdMZ;
	vector<double> vdIntensity;
	vector<int> viCharge; // the value of viCharge is zero for low-resolution MS2

	bool isMS1HighRes; // Is the MS1 scan a high-resolution scan?
	bool isMS2HighRes; // Is this MS2 scan a high-resolution scan?
	// vector<Peptide *> vpPeptides; // current set of peptides to be scored
	bool bSetMS2Flag; // false when preprocess fail on bad data

	// preprocess this scan, including
	// (1) remove noise peaks
	// (2) normalize intensity
	void preprocess();
	void preprocessMvh(multimap<double, double> *pIntenSortedPeakPreData);
	void preprocessXcorr();
	// score all peptides matched to this scan
	void scorePeptides();
	void scorePeptidesMVH(vector<double> *sequenceIonMasses, vector<double> *pdAAforward,
						  vector<double> *pdAAreverse, vector<char> *Seqs);
	void scorePeptidesXcorr(bool *pbDuplFragment, double *_pdAAforward, double *_pdAAreverse,
							unsigned int ***_uiBinnedIonMasses);
	void scorePeptidesLowMS2();
	void scorePeptidesHighMS2();

	// scoring functions
	void scoreWeightSum(tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple);
	void scoreRankSum(tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple);
	// void scoreRankSum_test(Peptide * currentPeptide);
	double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity);
	void scoreWeightSumHighMS2(tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple);
	void scoreRankSumHighMS2(tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple);
	// normalize raw scores
	void postprocess();
	double CalculateRankSum(double r1, double n1, double n2);

	void setScanType(string sScanType)
	{
		this->sScanType = sScanType;
	};
	string getScanType()
	{
		return this->sScanType;
	};

	//-----------Comet Begin-------------
	struct Query *pQuery;
	//-----------Comet End---------------
	//-----------Myrimatch Begin-------------
	bool bSkip;
	map<double, char> *peakData;
	PeakList *pPeakList;
	vector<int> *intenClassCounts;
	int totalPeakBins;
	double mzLowerBound;
	double mzUpperBound;
	//-----------Myrimatch End---------------
	//-----------Features--------------------
	double dSumIntensity;
	double dMaxIntensity;
	void sumIntensity();
	void scoreFeatureCalculation();
	void scoreFeatureCalculationWDPSip();
	int iNumPeptideAssigned;
	int getMaxNumProteinPsm();
	string sRTime;
	void setRTime(string _sRTime)
	{
		sRTime = _sRTime;
	};
	string getRTime() const
	{
		return sRTime;
	};
	bool isAnyScoreInTopN(int _iIndexPeptide, int _iRankThreshold);
	//-----------Features End----------------
	//-----------WDP End---------------------
	// for adding features after first time scoring
	double scoreWeightSum(string *currentPeptide, int measuredCharge,
					  vector<double> *pvdYionMass, vector<double> *pvdBionMass);
	// void scoreRankSum_test(Peptide * currentPeptide);
	// for adding features after first time scoring
	double scoreWeightSumHighMS2(const string *currentPeptide, const int measuredCharge,
								 const vector<vector<double>> *vvdYionMass,
								 const vector<vector<double>> *vvdYionProb,
								 const vector<vector<double>> *vvdBionMass,
								 const vector<vector<double>> *vvdBionProb);
	//-----------WDP End---------------------
};

#endif // MS2SCAN_H
