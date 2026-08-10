#ifndef MS2SCANVECTOR_H
#define MS2SCANVECTOR_H

#include <string>
#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>
#include <filesystem>

#include "ms2scan.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"
#include "MVH.h"
#include "CometSearchMod.h"
#include "RaxportHdf5Reader.h"
#include "performancelog.h"

#define ZERO            0.00000001
#define PEPTIDE_ARRAY_SIP_SIZE  200000

using namespace std;

namespace sipros
{
class FragmentIndex;
}

struct ScoredPsmRow
{
	int scanNumber = 0;
	int parentCharge = 0;
	int precursorScanNumber = 0;
	double isolationWindowCenterMZ = 0.0;
	double measuredParentMass = 0.0;
	double calculatedParentMass = 0.0;
	double precursorNeutronMass = 0.0;
	string scanType;
	string searchName;
	double ms2IsotopicAbundancePct = 1.07;
	double retentionTime = 0.0;
	float mvhScore = 0.0f;
	float xcorrScore = 0.0f;
	float wdpScore = 0.0f;
	int rank = 0;
	// Absolute RT difference in seconds to the precursor matched after scoring.
	// -1 is the explicit missing value when no peak matches in parent MS1 +/-5.
	double precursorRtDiffSeconds = -1.0;
	int ddaResidualRank = 0;
	float ddaResidualScore = 0.0f;
	int matchedBIons = 0;
	int matchedYIons = 0;
	int maxConsecutiveBIons = 0;
	int maxConsecutiveYIons = 0;
	string identifiedPeptide;
	string originalPeptide;
	string nakedPeptide;
	string proteinNames;
	bool isDecoy = false;
};

struct RegularSearchStatistics
{
	sipros::PerformanceTiming preprocess;
	sipros::PerformanceTiming queryIndex;
	sipros::PerformanceTiming mvhScoring;
	sipros::PerformanceTiming xcorrScoring;
	sipros::PerformanceTiming wdpScoring;
	sipros::PerformanceTiming psmFeatures;
	uint64_t scanCount = 0;
	uint64_t skippedScans = 0;
	uint64_t isolationWindowRanges = 0;
	uint64_t isolationWindowCandidates = 0;
	uint64_t exactMassCandidates = 0;
	uint64_t fragmentPostingsVisited = 0;
	uint64_t fragmentGateSurvivors = 0;
	uint64_t exactMvhCalls = 0;
	uint64_t exactMvhAccepted = 0;
	uint64_t candidatePsms = 0;
};

class MS2ScanVector {
	// All MS2 scans to be scored
	// the MS2 scans are sorted by their precursor masses
	vector<MS2Scan *> vpAllMS2Scans;
	// the precursor mass of these MS2 scans
	// this is used for quickly inserting a peptide into corrent MS2 scans
	// vpAllMS2Scans and vpPrecursorMasses are in the same order
	vector<double> vpPrecursorMasses;
	// for peptide assignment
	vector<tuple<double, int, MS2Scan *>> vAllPrecursorMassChargeMS2ScanPtrTuples;
	//vector <int> mass_w; // mass window
	string sScanFilename;    // the scan filename
	string sOutputFile; // the output file name
	sipros::RaxportReadOptions raxportReadOptions;
	sipros::RaxportMs1Data raxportMs1Data;
	shared_ptr<const sipros::FragmentIndex> fragmentIndex;
	RegularSearchStatistics regularSearchStatistics;

	// this should be moved to Peptide or shared through ProNovoConfig
	map<char, double> mapResidueMass; // mass except N and C termini;

	bool loadRaxportHdf5File();
	static bool mygreater(double i, double j);
	static bool myless(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);
	static bool mylessScanId(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);

	// find every MS2 scan whose precursor mass matches peptide mass
	bool assignPeptides2Scans(Peptide * currentPeptide);
	void estimateAndAssignPeptides(vector<Peptide *> &vpPeptideArray);
	pair<int, int> GetRangeFromMass(double lb, double ub);

	// regular search functions
	void preProcessAllMs2Mvh(); // pre-process all MS2 scans before mvh
	void searchDatabaseMvh(); // lossless precursor/fragment indexed regular search
	void postProcessAllMs2WdpXcorr(); // post processing after mvh scoring
	void postProcessAllMs2Wdp(); // post processing using wdp scoring
	void postProcessAllMs2Xcorr(); // post processing using xcorr scoring

	// sip search functions
	void preProcessAllMs2WdpSip();  // pre-process all MS2 scans before wdp for sip search
	void searchDatabaseWdpSip(); // search all ms2 scans against the protein list using wdp for sip search
	void processPeptideArrayWdpSip(vector<Peptide*>& vpPeptideArray); // process peptide array using wdp score for sip search
	void postProcessAllMs2MvhXcorr(); // post processing after wdp scoring
	void postProcessAllMs2XcorrSip(); // post processing using xcorr scoring for sip search
	void postProcessAllMs2MvhSip(); // post processing using mvh scoring for sip search

	void setOutputFile(const string & sScanFilenameInput, const string & sOutputDirectory);

	void GetAllRangeFromMass(double dPeptideMass,
						 double precursorNeutronMass,
						 vector<pair<int, int> > & vpPeptideMassRanges);
	string ParsePath(string sPath);

public:
	MS2ScanVector(const string &sScanFilenameInput,
				  const string &sOutputDirectory,
				  const sipros::RaxportReadOptions &readOptions = sipros::RaxportReadOptions{});
	~MS2ScanVector();

	// Populate vpAllMS2Scans from a Raxport schema v6 HDF5 file.
	// Return false if there is a problem with the file.
	bool loadMassData();
	size_t scanCount() const { return vpAllMS2Scans.size(); }
	size_t precursorHypothesisCount() const
	{
		return vAllPrecursorMassChargeMS2ScanPtrTuples.size();
	}
	sipros::RaxportMs1Data releaseMs1Data()
	{
		return std::move(raxportMs1Data);
	}
	void configureFragmentIndex(
		const shared_ptr<const sipros::FragmentIndex> &sharedIndex)
	{
		fragmentIndex = sharedIndex;
	}
	const RegularSearchStatistics &regularStatistics() const
	{
		return regularSearchStatistics;
	}

	// begin the database searching
	void startProcessingMvh(); // start functions to process the loaded HDF5 file using mvh as the prime score
	void startProcessingWdpSip(); // start functions to process the loaded HDF5 file with WDP as prime score without tasking
	void clearSearchResults();
	void appendScoredPsmRows(vector<ScoredPsmRow> &rows, bool isDecoy, int topKeep, double ms2IsotopicAbundancePct = 1.07) const;
	size_t matchScoredPsmPrecursors(vector<ScoredPsmRow> &rows) const;

	// variables for the MVH thread
	vector<vector<double> *> _ppdAAforward;
	vector<vector<double> *> _ppdAAreverse;
	vector<vector<double> *> psequenceIonMasses;
	vector<vector<char> *> pSeqs;
	int num_max_threads;
	void preMvh();
	void postMvh();
	// variables for the Xcorr thread
	vector<bool *> vpbDuplFragmentGlobal;
	vector<double *> v_pdAAforwardGlobal;
	vector<double *> v_pdAAreverseGlobal;
	vector<unsigned int ***> v_uiBinnedIonMassesGlobal;
	vector<vector<unsigned char> > vvpbDuplFragmentGlobal;
	vector<vector<double> > vvdBinnedIonMassesGlobal;
	vector<vector<int> > vvdBinGlobal;
	void preXcorr();
	void postXcorr();
};

#endif // MS2SCANVECTOR_H
