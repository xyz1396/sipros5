#include "ms2scanvector.h"
#include "PeptideIsotopeCalculator.h"
#include "RaxportHdf5Reader.h"
#include "fragmentindex.h"
#include "indexeddatabasesearch.h"
#include "performancelog.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

MS2ScanVector::MS2ScanVector(const string &sScanFilenameInput,
							 const string &sOutputDirectory,
							 const sipros::RaxportReadOptions &readOptions)
	: raxportReadOptions(readOptions)
{
	unsigned int n;
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;
	sScanFilename = sScanFilenameInput;
	setOutputFile(sScanFilenameInput, sOutputDirectory);
	for (n = 0; n < vsSingleResidueNames.size(); ++n)
		mapResidueMass[vsSingleResidueNames[n][0]] = vdSingleResidueMasses[n];

}

MS2ScanVector::~MS2ScanVector()
{
	// the destructors will free memory from vpAllMS2Scans
	vector<MS2Scan *>::iterator it;
	for (it = vpAllMS2Scans.begin(); it != vpAllMS2Scans.end(); ++it)
	{
		delete (*it);
	}
	vpAllMS2Scans.clear();
	vpPrecursorMasses.clear();
}

void MS2ScanVector::setOutputFile(const string &sScanFilenameInput, const string &sOutputDirectory)
{
	std::filesystem::path inPath = sScanFilenameInput;
	std::filesystem::path outPath = sOutputDirectory;
	// get file name without extension and path
	std::string baseName = inPath.stem().string();
	string searchName;
	if ((ProNovoConfig::getSearchName() == "") || (ProNovoConfig::getSearchName() == "Null") || (ProNovoConfig::getSearchName() == "NULL") || (ProNovoConfig::getSearchName() == "null"))
		searchName = "Null";
	else
	{
		searchName = ProNovoConfig::getSearchName();
		std::replace(searchName.begin(), searchName.end(), '.', '_');
	}
	std::filesystem::path outFileName = baseName + "." + searchName + ".pin";
	outPath = outPath / outFileName;
	sOutputFile = outPath.string();
}

bool MS2ScanVector::loadRaxportHdf5File()
{
	std::string error;
	if (!sipros::readRaxportHdf5Scans(
			sScanFilename, vpAllMS2Scans, &raxportMs1Data, error, nullptr,
			raxportReadOptions))
	{
		std::cerr << error << std::endl;
		return false;
	}
	if (vpAllMS2Scans.empty())
	{
		std::cerr << "No ms_order == 2 scans found in Raxport HDF5 file: " << sScanFilename << std::endl;
		return false;
	}

	vAllPrecursorMassChargeMS2ScanPtrTuples.clear();
	vpPrecursorMasses.clear();
	vAllPrecursorMassChargeMS2ScanPtrTuples.reserve(vpAllMS2Scans.size());
	vpPrecursorMasses.reserve(vpAllMS2Scans.size());

	for (MS2Scan *scan : vpAllMS2Scans)
	{
		if (scan == nullptr)
		{
			continue;
		}
		std::vector<std::pair<double, int>> scanPrecursorMassCharges;
		const bool useReactionCharge =
			raxportReadOptions.precursorSource == sipros::PrecursorSource::RaxportCandidates;
		scan->bUseReactionChargeForScoring = useReactionCharge;
		scan->iMaxCandidateCharge = 0;
		if (!useReactionCharge)
		{
			// Defensively recompute these limits solely from the selected MS1 peak
			// hypotheses below; the reaction center remains output metadata only.
			scan->dParentNeutralMass = 0.0;
			scan->dParentMass = 0.0;
		}

		auto appendPrecursorMz = [&](double mz, int charge)
		{
			if (charge <= 0 || mz <= 0.0)
			{
				return;
			}
			const double chargedMass = mz * charge;
			const double parentNeutralMass = chargedMass -
				static_cast<double>(charge) * ProNovoConfig::getProtonMass();
			const double deduplicationTolerance = useReactionCharge
				? ProNovoConfig::getMassAccuracyParentIon()
				: 0.0;
			for (const auto &existing : scanPrecursorMassCharges)
			{
				if (existing.second == charge &&
					std::abs(existing.first - parentNeutralMass) <= deduplicationTolerance)
				{
					return;
				}
			}
			scanPrecursorMassCharges.push_back({parentNeutralMass, charge});
			scan->iMaxCandidateCharge = std::max(scan->iMaxCandidateCharge, charge);
			vAllPrecursorMassChargeMS2ScanPtrTuples.push_back({parentNeutralMass, charge, scan});
			if (parentNeutralMass > scan->dParentNeutralMass)
			{
				scan->dParentNeutralMass = parentNeutralMass;
			}
			if (chargedMass > scan->dParentMass)
			{
				scan->dParentMass = chargedMass;
			}
		};

		const size_t nCandidates = std::min(
			scan->dParentMZs.size(), scan->iParentChargeStates.size());
		for (size_t j = 0; j < nCandidates; ++j)
		{
			appendPrecursorMz(
				scan->dParentMZs[j], scan->iParentChargeStates[j]);
		}
		if (!useReactionCharge)
		{
			for (const auto &window : scan->vIsolationWindowsMz)
			{
				const double upperMz = window.first + window.second / 2.0;
				for (int charge = 1; charge <= 4; ++charge)
				{
					const double chargedMass = upperMz * charge;
					const double neutralMass = chargedMass -
						static_cast<double>(charge) *
							ProNovoConfig::getProtonMass();
					scan->dParentMass =
						std::max(scan->dParentMass, chargedMass);
					scan->dParentNeutralMass =
						std::max(scan->dParentNeutralMass, neutralMass);
					scan->iMaxCandidateCharge =
						std::max(scan->iMaxCandidateCharge, charge);
				}
			}
		}

		if (useReactionCharge && scan->iParentChargeState <= 0 && nCandidates > 0)
		{
			// Downstream direct-SIP preprocessing requires one representative charge;
			// scored candidates retain their native precursor charges.
			scan->iParentChargeState = 3;
		}
	}

	std::sort(vAllPrecursorMassChargeMS2ScanPtrTuples.begin(), vAllPrecursorMassChargeMS2ScanPtrTuples.end(),
			  [](const std::tuple<double, int, MS2Scan *> &a, const std::tuple<double, int, MS2Scan *> &b)
			  { return std::get<0>(a) < std::get<0>(b); });
	for (const auto &entry : vAllPrecursorMassChargeMS2ScanPtrTuples)
	{
		vpPrecursorMasses.push_back(std::get<0>(entry));
	}
	return true;
}

bool MS2ScanVector::loadMassData()
{
	bool bReVal = false;
	std::string filenameStr = this->sScanFilename;
	std::transform(filenameStr.begin(), filenameStr.end(), filenameStr.begin(),
				   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	std::string fileNameSuffix;
	const size_t dot = filenameStr.rfind('.');
	if (dot != std::string::npos)
	{
		fileNameSuffix = filenameStr.substr(dot + 1);
	}

	if (fileNameSuffix == "h5" || fileNameSuffix == "hdf5")
	{
		bReVal = loadRaxportHdf5File();
	}
	else
	{
		std::cerr << "Raxport HDF5 scan input required (.h5 or .hdf5); unsupported file: "
			  << this->sScanFilename << std::endl;
	}

	double mass = 0;
	int charge = 0;
	MS2Scan *pMS2Scan = nullptr;
	for (const auto &tuple : vAllPrecursorMassChargeMS2ScanPtrTuples)
	{
		std::tie(mass, charge, pMS2Scan) = tuple;
		if (pMS2Scan->dParentNeutralMass < mass)
			pMS2Scan->dParentNeutralMass = mass;
		mass = mass + charge * ProNovoConfig::getProtonMass();
		if (pMS2Scan->dParentMass < mass)
			pMS2Scan->dParentMass = mass;
		if (ProNovoConfig::dMaxMS2ScanMass < mass)
			ProNovoConfig::dMaxMS2ScanMass = mass;
		if (ProNovoConfig::iMaxPercusorCharge < charge)
			ProNovoConfig::iMaxPercusorCharge = charge;
	}
	for (const MS2Scan *scan : vpAllMS2Scans)
	{
		if (scan == nullptr)
			continue;
		ProNovoConfig::dMaxMS2ScanMass =
			std::max(ProNovoConfig::dMaxMS2ScanMass, scan->dParentMass);
		ProNovoConfig::iMaxPercusorCharge =
			std::max(ProNovoConfig::iMaxPercusorCharge, scan->iMaxCandidateCharge);
	}
	return bReVal;
}

bool MS2ScanVector::mygreater(double i, double j)
{
	return (i > j);
}

bool MS2ScanVector::myless(MS2Scan *pMS2Scan1, MS2Scan *pMS2Scan2)
{
	return (pMS2Scan1->dParentNeutralMass < pMS2Scan2->dParentNeutralMass);
}

void MS2ScanVector::preProcessAllMs2Mvh()
{
	const sipros::PerformanceTimer timer;
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	int num_threads = omp_get_max_threads();
	vector<multimap<double, double> *> vpIntenSortedPeakPreData;
	for (int j = 0; j < num_threads; ++j)
	{
		vpIntenSortedPeakPreData.push_back(new multimap<double, double>());
	}

#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
	{
		// vpAllMS2Scans.at(i)->preprocess();
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans.at(i)->preprocessMvh(vpIntenSortedPeakPreData.at(iThreadId));
	}

	for (int j = 0; j < num_threads; ++j)
	{
		delete vpIntenSortedPeakPreData.at(j);
	}
	vpIntenSortedPeakPreData.clear();
	int maxPeakBins = 0;
	int iNumSkippedScans = 0;
	for (i = 0; i < iScanSize; i++)
	{
		if (vpAllMS2Scans.at(i)->bSkip)
		{
			iNumSkippedScans++;
		}
		if (vpAllMS2Scans.at(i)->totalPeakBins > maxPeakBins)
		{
			maxPeakBins = vpAllMS2Scans.at(i)->totalPeakBins;
		}
	}
	MVH::initialLnTable(maxPeakBins);

	{
#pragma omp parallel for schedule(guided)
		for (i = 0; i < iScanSize; i++)
		{
			vpAllMS2Scans.at(i)->sumIntensity();
		}
	}

	regularSearchStatistics.preprocess = timer.elapsed();
	regularSearchStatistics.scanCount = vpAllMS2Scans.size();
	regularSearchStatistics.skippedScans =
		static_cast<uint64_t>(iNumSkippedScans);
}

void MS2ScanVector::preProcessAllMs2WdpSip()
{
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();
	cout << "  Preprocessing scans: " << vpAllMS2Scans.size() << endl;

#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
		vpAllMS2Scans.at(i)->preprocess();
}

void MS2ScanVector::GetAllRangeFromMass(
	double dPeptideMass,
	double precursorNeutronMass,
	vector<std::pair<int, int>> &vpPeptideMassRanges)
// all ranges of MS2 scans are stored in  vpPeptideMassWindows
{
	int i;
	pair<int, int> pairMS2Range;
	pair<int, int> lastPairRange(-100, -100);
	vector<pair<double, double>> vpPeptideMassWindows;
	vpPeptideMassWindows.clear();
	vpPeptideMassRanges.clear();
	if (!ProNovoConfig::getPeptideMassWindows(
			dPeptideMass, precursorNeutronMass, vpPeptideMassWindows))
		throw std::runtime_error(
			"FASTA peptide is missing a valid composition-weighted precursor neutron mass.");
	for (i = 0; i < (int)vpPeptideMassWindows.size(); i++)
	{
		pairMS2Range = GetRangeFromMass(vpPeptideMassWindows.at(i).first, vpPeptideMassWindows.at(i).second);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
		{
			if ((lastPairRange.first < 0) || (lastPairRange.second < 0))
				lastPairRange = pairMS2Range;
			else
			{
				if (lastPairRange.second > pairMS2Range.first)
					lastPairRange.second = pairMS2Range.second;
				else
				{
					vpPeptideMassRanges.push_back(lastPairRange);
					lastPairRange = pairMS2Range;
				}
			}
		}
	}
	if ((lastPairRange.first > -1) && (lastPairRange.second > -1))
		vpPeptideMassRanges.push_back(lastPairRange);
}

pair<int, int> MS2ScanVector::GetRangeFromMass(double lb, double ub)
// lb and ub are lower and upper bounds of acceptable parent mass values
{
	pair<int, int> p;
	int low = 0, high = static_cast<int>(vpPrecursorMasses.size()) - 1, mid;
	double target;
	target = (lb + ub) / 2.0;
	// double lb, ub;  // lower and upper bounds on acceptable parent mass values
	// ub = target + error;
	// lb = target - error;
	while ((high - low) > 1)
	{
		mid = (high + low) / 2;
		if (vpPrecursorMasses[mid] > target)
			high = mid;
		else
			low = mid;
	}

	// Iterate till we get to the first element > than the lower bound
	int ndx = low;
	// cout<<scan_mass_list_.size()<<endl;
	if (vpPrecursorMasses[ndx] >= lb)
	{
		while (ndx >= 0 && vpPrecursorMasses[ndx] >= lb)
			ndx--;
		ndx++;
	}
	else
		while (ndx < (int)vpPrecursorMasses.size() && vpPrecursorMasses[ndx] < lb)
			ndx++;
	if (ndx == (int)vpPrecursorMasses.size() || vpPrecursorMasses[ndx] > ub)
		p = make_pair(-1, -1);
	else
	{
		low = ndx;
		while (ndx < (int)vpPrecursorMasses.size() && vpPrecursorMasses[ndx] <= ub)
			ndx++;
		high = ndx - 1;
		p = make_pair(low, high);
	}

	if (p.first == -1)
		return p;
	if (vpPrecursorMasses[p.first] < lb)
		cerr << "ERROR L " << vpPrecursorMasses[p.first] << " " << lb << endl;
	if (vpPrecursorMasses[p.second] > ub)
		cerr << "ERROR U " << vpPrecursorMasses[p.second] << " " << ub << endl;
	return p;
}

/*
 void MS2ScanVector::assignPeptides2Scans(Peptide * currentPeptide)
 {
 int i;
 pair<int, int> pairMS2Range;
 //   cout<<currentPeptide->getPeptideMass()<<endl;
 pairMS2Range = GetRangeFromMass(currentPeptide->getPeptideMass(),
 ProNovoConfig::getMassAccuracyParentIon());
 if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
 for (i= pairMS2Range.first; i<= pairMS2Range.second; i++)
 vpAllMS2Scans.at(i)->vpPeptides.push_back(currentPeptide);
 }*/

bool MS2ScanVector::assignPeptides2Scans(Peptide *currentPeptide)
{
	int i, j;
	bool bAssigned = false;
	vector<pair<int, int>> vpPeptideMassRanges;
	pair<int, int> pairMS2Range;

	GetAllRangeFromMass(currentPeptide->getPeptideMass(),
						currentPeptide->getPrecursorNeutronMass(),
						vpPeptideMassRanges);

	for (j = 0; j < (int)vpPeptideMassRanges.size(); j++)
	{
		pairMS2Range = vpPeptideMassRanges.at(j);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
		{
			for (i = pairMS2Range.first; i <= pairMS2Range.second; i++)
			{ // vpAllMS2ScanPtrs.at(i)->vpPeptides.push_back(currentPeptide);
				// for DIA and DDA with large isolation window
				tuple<double, int, Peptide *> currentMassChargePeptidePtrTuple =
					{
						get<0>(vAllPrecursorMassChargeMS2ScanPtrTuples[i]),
						get<1>(vAllPrecursorMassChargeMS2ScanPtrTuples[i]),
						currentPeptide};
				get<2>(vAllPrecursorMassChargeMS2ScanPtrTuples[i])
					->vMassChargePeptidePtrTuples.push_back(currentMassChargePeptidePtrTuple);
			}
			bAssigned = true;
		}
	}
	return bAssigned;
}

void MS2ScanVector::estimateAndAssignPeptides(
	vector<Peptide *> &vpPeptideArray)
{
	std::exception_ptr estimateError;
	const int peptideCount = static_cast<int>(vpPeptideArray.size());

#pragma omp parallel
	{
		PeptideIsotopeCalculator calculator;
#pragma omp for schedule(guided)
		for (int i = 0; i < peptideCount; ++i)
		{
			try
			{
				const string &decorated = vpPeptideArray[i]->getPeptideSeq();
				string compositionSequence;
				compositionSequence.reserve(decorated.size());
				for (char symbol : decorated)
					if (symbol != '[' && symbol != ']')
						compositionSequence.push_back(symbol);

				const auto estimate =
					calculator.calPrecursorEstimate(compositionSequence);
				vpPeptideArray[i]->setPeptideMass(estimate.mass);
				vpPeptideArray[i]->setPrecursorNeutronMass(
					estimate.neutronMass);
			}
			catch (...)
			{
#pragma omp critical(sipros_precursor_estimate_error)
				{
					if (!estimateError)
						estimateError = std::current_exception();
				}
			}
		}
	}

	if (estimateError)
	{
		for (Peptide *peptide : vpPeptideArray)
			delete peptide;
		vpPeptideArray.clear();
		std::rethrow_exception(estimateError);
	}

	// Scan candidate vectors are mutable, so assignment remains serial after
	// every candidate has received its final modal precursor estimate.
	vector<Peptide *> assigned;
	assigned.reserve(vpPeptideArray.size());
	for (Peptide *peptide : vpPeptideArray)
	{
		if (assignPeptides2Scans(peptide))
		{
			assigned.push_back(peptide);
			if (peptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass)
				ProNovoConfig::dMaxPeptideMass = peptide->getPeptideMass();
		}
		else
		{
			delete peptide;
		}
	}
	vpPeptideArray.swap(assigned);
}

void MS2ScanVector::processPeptideArrayWdpSip(vector<Peptide *> &vpPeptideArray)
{
	int iPeptideArraySize, iScanSize;
	estimateAndAssignPeptides(vpPeptideArray);
	iPeptideArraySize = (int)vpPeptideArray.size();
	if (iPeptideArraySize == 0)
		return;

	//    for (int i=0; i< (int) vpPeptideArray.size(); i++)
	//      cout<<vpPeptideArray.at(i)->getPeptideSeq() <<"\t"
	//	  <<vpPeptideArray.at(i)->getOriginalPeptideSeq()<<"\t"
	//	  <<vpPeptideArray.at(i)->getProteinName()<<"\t"
	//	  <<vpPeptideArray.at(i)->getBeginPosProtein()<<"\t"
	//	  <<vpPeptideArray.at(i)->getPeptideMass()<<endl;

	// cout<<"calculating fragments of "<<vpPeptideArray.size()<<"  peptides"<<endl;
#pragma omp parallel for shared(vpPeptideArray) \
	schedule(guided)

	for (int i = 0; i < iPeptideArraySize; i++)
		vpPeptideArray[i]->preprocessing(vpAllMS2Scans.at(0)->isMS2HighRes, mapResidueMass);
	// vpPeptideArray[i]->calculateExpectedFragments(mapResidueMass);

	// cout<<"scoring "<<vpPeptideArray.size()<<"  peptides"<<endl;
	//  every MS2 scans scores their matched peptides

	iScanSize = (int)vpAllMS2Scans.size();
#pragma omp parallel for schedule(guided)

	for (int i = 0; i < iScanSize; i++)
		vpAllMS2Scans[i]->scorePeptides();

	// free memory of all peptide objects
	for (int i = 0; i < static_cast<int>(vpPeptideArray.size()); i++)
		delete vpPeptideArray[i];

	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::searchDatabaseMvh()
{
	this->preMvh();
	if (!fragmentIndex)
	{
		this->postMvh();
		throw std::runtime_error(
			"fragment-index search was started without a prepared index");
	}
	const sipros::FragmentIndex &index = *fragmentIndex;
	if (index.maximumPeptideMass() > ProNovoConfig::dMaxPeptideMass)
	{
		ProNovoConfig::dMaxPeptideMass = index.maximumPeptideMass();
	}

	const sipros::PerformanceTimer querySetupTimer;
	const int threadCount = omp_get_max_threads();
	std::vector<sipros::IndexedSearchScratch> scratch(
		static_cast<size_t>(threadCount));
	std::vector<sipros::IndexedSearchCounters> threadCounters(
		static_cast<size_t>(threadCount));

	sipros::IndexedSearchCounters queryTotals;
	sipros::IndexedSearchCounters mvhTotals;
	sipros::PerformanceTiming queryTiming = querySetupTimer.elapsed();
	sipros::PerformanceTiming mvhTiming;
	// Full-isolation-window candidate vectors are broader than point-MS1
	// queries.  Compact batches cap transient survivor memory before scoring.
	constexpr size_t ScanBatchSize = 64;
	for (size_t batchBegin = 0; batchBegin < vpAllMS2Scans.size();
		 batchBegin += ScanBatchSize)
	{
		const size_t batchEnd = std::min(
			vpAllMS2Scans.size(), batchBegin + ScanBatchSize);
		std::vector<std::vector<sipros::IndexedCandidate>> batchCandidates(
			batchEnd - batchBegin);
		threadCounters.assign(
			static_cast<size_t>(threadCount), sipros::IndexedSearchCounters{});
		const sipros::PerformanceTimer queryTimer;
#pragma omp parallel for schedule(guided)
		for (int64_t i = static_cast<int64_t>(batchBegin);
			 i < static_cast<int64_t>(batchEnd); ++i)
		{
			const int threadId = omp_get_thread_num();
			sipros::queryIndexedScan(
				index,
				*vpAllMS2Scans[static_cast<size_t>(i)],
				scratch[static_cast<size_t>(threadId)],
				batchCandidates[static_cast<size_t>(i) - batchBegin],
				threadCounters[static_cast<size_t>(threadId)]);
		}
		queryTiming += queryTimer.elapsed();
		for (const auto &counter : threadCounters)
		{
			queryTotals += counter;
		}

		threadCounters.assign(
			static_cast<size_t>(threadCount), sipros::IndexedSearchCounters{});
		const sipros::PerformanceTimer mvhTimer;
#pragma omp parallel for schedule(guided)
		for (int64_t i = static_cast<int64_t>(batchBegin);
			 i < static_cast<int64_t>(batchEnd); ++i)
		{
			const int threadId = omp_get_thread_num();
			sipros::scoreIndexedScanMvh(
				index,
				*vpAllMS2Scans[static_cast<size_t>(i)],
				batchCandidates[static_cast<size_t>(i) - batchBegin],
				*psequenceIonMasses[static_cast<size_t>(threadId)],
				*_ppdAAforward[static_cast<size_t>(threadId)],
				*_ppdAAreverse[static_cast<size_t>(threadId)],
				*pSeqs[static_cast<size_t>(threadId)],
				threadCounters[static_cast<size_t>(threadId)]);
		}
		mvhTiming += mvhTimer.elapsed();
		for (const auto &counter : threadCounters)
		{
			mvhTotals += counter;
		}
	}
	std::vector<sipros::IndexedSearchScratch>().swap(scratch);
	regularSearchStatistics.queryIndex = queryTiming;
	regularSearchStatistics.mvhScoring = mvhTiming;
	regularSearchStatistics.isolationWindowRanges =
		queryTotals.isolationWindowRanges;
	regularSearchStatistics.isolationWindowCandidates =
		queryTotals.isolationWindowCandidates;
	regularSearchStatistics.exactMassCandidates =
		queryTotals.exactMassCandidates;
	regularSearchStatistics.fragmentPostingsVisited =
		queryTotals.fragmentPostingsVisited;
	regularSearchStatistics.fragmentGateSurvivors =
		queryTotals.fragmentGateSurvivors;
	regularSearchStatistics.exactMvhCalls = mvhTotals.exactMvhCalls;
	regularSearchStatistics.exactMvhAccepted = mvhTotals.exactMvhAccepted;

	this->postMvh();
	MVH::destroyLnTable();
	PeptideUnit::iNumScores = 1;
}

void MS2ScanVector::searchDatabaseWdpSip()
{
	CLOCKSTART;
	ProteinDatabase myProteinDatabase;
	vector<Peptide *> vpPeptideArray;
	Peptide *currentPeptide;
	myProteinDatabase.loadDatabase();

	if (myProteinDatabase.getFirstProtein())
	{
		currentPeptide = new Peptide;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide))
		{
			vpPeptideArray.push_back(currentPeptide);
			// create a new peptide for the next iteration
			currentPeptide = new Peptide;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIP_SIZE)
			{
				processPeptideArrayWdpSip(vpPeptideArray);
			}
		}
		// the last peptide object is an empty object and need to be deleted
		delete currentPeptide;
		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// Process the final partial SIP peptide batch.
		if (!vpPeptideArray.empty())
		{
			processPeptideArrayWdpSip(vpPeptideArray);
		}

		PeptideUnit::iNumScores = 1;
	}
	cout << "  WDP search done" << endl;
	CLOCKSTOP;
}

void MS2ScanVector::startProcessingMvh()
{
	regularSearchStatistics = RegularSearchStatistics{};
	// Preprocessing all MS2 scans by mult-threading
	preProcessAllMs2Mvh();

	// Search all MS2 scans with the prepared precursor/fragment index.
	searchDatabaseMvh();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2WdpXcorr();

}

void MS2ScanVector::startProcessingWdpSip()
{
	// SIP abundance changes the modal masses of biosynthetic residues and PTMs.
	// Refresh the low-resolution fragment lookup for every abundance iteration.
	mapResidueMass.clear();
	for (size_t index = 0;
		 index < ProNovoConfig::vsSingleResidueNames.size();
		 ++index)
	{
		const std::string &name = ProNovoConfig::vsSingleResidueNames[index];
		if (name.size() == 1 &&
			index < ProNovoConfig::vdSingleResidueMasses.size())
		{
			mapResidueMass[name.front()] =
				ProNovoConfig::vdSingleResidueMasses[index];
		}
	}

	// Preprocessing all MS2 scans by multi-threading
	preProcessAllMs2WdpSip();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseWdpSip();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2MvhXcorr();

}

void MS2ScanVector::postProcessAllMs2WdpXcorr()
{
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();
	uint64_t candidateCount = 0;
	for (const MS2Scan *scan : vpAllMS2Scans)
	{
		candidateCount += scan->vpWeightSumTopPeptides.size();
	}
	regularSearchStatistics.candidatePsms = candidateCount;

	const sipros::PerformanceTimer xcorrTimer;
	postProcessAllMs2Xcorr();
	regularSearchStatistics.xcorrScoring = xcorrTimer.elapsed();

	const sipros::PerformanceTimer wdpTimer;
	postProcessAllMs2Wdp();
	regularSearchStatistics.wdpScoring = wdpTimer.elapsed();

	const sipros::PerformanceTimer featureTimer;
#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
	{
		vpAllMS2Scans.at(i)->scoreFeatureCalculation();
	}
	regularSearchStatistics.psmFeatures = featureTimer.elapsed();
}

void MS2ScanVector::postProcessAllMs2MvhXcorr()
{
	CLOCKSTART;
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	postProcessAllMs2XcorrSip();
	cout << "  Xcorr search done" << endl;

	postProcessAllMs2MvhSip();
	cout << "  MVH search done" << endl;

#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
	{
		vpAllMS2Scans.at(i)->scoreFeatureCalculationWDPSip();
	}
	CLOCKSTOP;
}

void MS2ScanVector::postProcessAllMs2Wdp()
{
	vector<vector<vector<double>> *> vpvvdYionMass;
	vector<vector<vector<double>> *> vpvvdYionProb;
	vector<vector<vector<double>> *> vpvvdBionMass;
	vector<vector<vector<double>> *> vpvvdBionProb;
	vector<vector<double> *> vpvdYionMass;
	vector<vector<double> *> vpvdBionMass;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i)
	{
		vpvvdYionMass.push_back(new vector<vector<double>>());
		vpvvdYionProb.push_back(new vector<vector<double>>());
		vpvvdBionMass.push_back(new vector<vector<double>>());
		vpvvdBionProb.push_back(new vector<vector<double>>());
		vpvdYionMass.push_back(new vector<double>());
		vpvdBionMass.push_back(new vector<double>());
	}

	int iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanSize; i++)
	{
		vpAllMS2Scans.at(i)->preprocess();
		bool bHighRes = vpAllMS2Scans.at(i)->isMS2HighRes;
		int iThreadId = omp_get_thread_num();
		for (int j = 0; j < ((int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size()); j++)
		{
			Peptide::preprocessing(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring, bHighRes,
								   mapResidueMass, vpvvdYionMass.at(iThreadId), vpvvdYionProb.at(iThreadId),
								   vpvvdBionMass.at(iThreadId), vpvvdBionProb.at(iThreadId), vpvdYionMass.at(iThreadId),
								   vpvdBionMass.at(iThreadId));
			double dWeightSum = 0;
			if (bHighRes)
			{
				dWeightSum = vpAllMS2Scans.at(i)->scoreWeightSumHighMS2(
					&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
					vpAllMS2Scans[i]->vpWeightSumTopPeptides[j]->iMeasuredParentCharge,
					vpvvdYionMass.at(iThreadId), vpvvdYionProb.at(iThreadId), vpvvdBionMass.at(iThreadId),
					vpvvdBionProb.at(iThreadId));
			}
			else
			{
				const int wdpCharge = vpAllMS2Scans.at(i)->bUseReactionChargeForScoring
					? vpAllMS2Scans.at(i)->iParentChargeState
					: vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->iMeasuredParentCharge;
				dWeightSum = vpAllMS2Scans.at(i)->scoreWeightSum(
					&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
					wdpCharge,
					vpvdYionMass.at(iThreadId), vpvdBionMass.at(iThreadId));
			}
			vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores[0] = dWeightSum;
		}
	}

	for (int i = 0; i < num_max_threads; ++i)
	{
		delete vpvvdYionMass.at(i);
		delete vpvvdYionProb.at(i);
		delete vpvvdBionMass.at(i);
		delete vpvvdBionProb.at(i);
		delete vpvdYionMass.at(i);
		delete vpvdBionMass.at(i);
	}
	++PeptideUnit::iNumScores;
}

void MS2ScanVector::postProcessAllMs2Xcorr()
{
	// MH: Must be equal to largest possible array
	int iArraySizePreprocess = (int)((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0) * ProNovoConfig::dHighResInverseBinWidth);
	// MH: Must be equal to largest possible array
	int iArraySizeScore = (int)((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	vector<double *> vpdTmpRawData;
	vector<double *> vpdTmpFastXcorrData;
	vector<double *> vpdTmpCorrelationData;
	vector<double *> vpdTmpSmoothedSpectrum;
	vector<double *> vpdTmpPeakExtracted;
	vector<bool *> vpbDuplFragment;
	vector<double *> v_pdAAforward;
	vector<double *> v_pdAAreverse;
	vector<unsigned int ***> v_uiBinnedIonMasses;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i)
	{
		// CometSearchMod::Preprocess clears the active range before every use.
		// Avoid eagerly touching the full global-maximum allocation here.
		vpdTmpRawData.push_back(new double[iArraySizePreprocess]);
		vpdTmpFastXcorrData.push_back(new double[iArraySizePreprocess]);
		vpdTmpCorrelationData.push_back(new double[iArraySizePreprocess]);
		vpdTmpSmoothedSpectrum.push_back(new double[iArraySizePreprocess]);
		vpdTmpPeakExtracted.push_back(new double[iArraySizePreprocess]);
		vpbDuplFragment.push_back(new bool[iArraySizeScore]());
		v_pdAAforward.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverse.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int ***_uiBinnedIonMasses = new unsigned int **[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			_uiBinnedIonMasses[ii] = new unsigned int *[9]();
			for (int j = 0; j < 9; j++)
			{
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMasses.push_back(_uiBinnedIonMasses);
	}

	int iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanSize; i++)
	{
		int iThreadId = omp_get_thread_num();
		struct Query *pQuery = new Query();
		if (!CometSearchMod::Preprocess(pQuery, vpAllMS2Scans.at(i), vpdTmpRawData.at(iThreadId),
										vpdTmpFastXcorrData.at(iThreadId), vpdTmpCorrelationData.at(iThreadId),
										vpdTmpSmoothedSpectrum.at(iThreadId), vpdTmpPeakExtracted.at(iThreadId)))
		{
			cout << "Error Post Xcorr." << endl;
			exit(1);
		}
		else
		{
			if (vpAllMS2Scans.at(i)->pQuery != NULL)
			{
				delete vpAllMS2Scans.at(i)->pQuery;
				vpAllMS2Scans.at(i)->pQuery = NULL;
			}
			vpAllMS2Scans.at(i)->pQuery = pQuery;
			for (int j = 0; j < (int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++)
			{
				double dXcorr = 0;
				CometSearchMod::ScorePeptides(&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
										  vpbDuplFragment.at(iThreadId), v_pdAAforward.at(iThreadId), v_pdAAreverse.at(iThreadId),
										  vpAllMS2Scans.at(i),
										  vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->iMeasuredParentCharge,
										  v_uiBinnedIonMasses.at(iThreadId), dXcorr);
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores[1] = dXcorr;
			}
			delete pQuery;
			vpAllMS2Scans.at(i)->pQuery = NULL;
		}
	}

	for (int i = 0; i < num_max_threads; ++i)
	{
		delete[] vpdTmpRawData.at(i);
		delete[] vpdTmpFastXcorrData.at(i);
		delete[] vpdTmpCorrelationData.at(i);
		delete[] vpdTmpSmoothedSpectrum.at(i);
		delete[] vpdTmpPeakExtracted.at(i);
		delete[] vpbDuplFragment.at(i);
		delete[] v_pdAAforward.at(i);
		delete[] v_pdAAreverse.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			for (int j = 0; j < 9; j++)
			{
				delete[] v_uiBinnedIonMasses.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMasses.at(i)[ii];
		}
		delete[] v_uiBinnedIonMasses.at(i);
	}
	++PeptideUnit::iNumScores;
}

void MS2ScanVector::postProcessAllMs2XcorrSip()
{
	// MH: Must be equal to largest possible array
	int iArraySizePreprocess = (int)((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0) * ProNovoConfig::dHighResInverseBinWidth);
	// MH: Must be equal to largest possible array
	int iArraySizeScore = (int)((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	vector<double *> vpdTmpRawData;
	vector<double *> vpdTmpFastXcorrData;
	vector<double *> vpdTmpCorrelationData;
	vector<double *> vpdTmpSmoothedSpectrum;
	vector<double *> vpdTmpPeakExtracted;
	vector<bool *> vpbDuplFragment;
	vector<double *> v_pdAAforward;
	vector<double *> v_pdAAreverse;
	vector<unsigned int ***> v_uiBinnedIonMasses;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i)
	{
		vpdTmpRawData.push_back(new double[iArraySizePreprocess]());
		vpdTmpFastXcorrData.push_back(new double[iArraySizePreprocess]());
		vpdTmpCorrelationData.push_back(new double[iArraySizePreprocess]());
		vpdTmpSmoothedSpectrum.push_back(new double[iArraySizePreprocess]());
		vpdTmpPeakExtracted.push_back(new double[iArraySizePreprocess]());
		vpbDuplFragment.push_back(new bool[iArraySizeScore]());
		v_pdAAforward.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverse.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int ***_uiBinnedIonMasses = new unsigned int **[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			_uiBinnedIonMasses[ii] = new unsigned int *[9]();
			for (int j = 0; j < 9; j++)
			{
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMasses.push_back(_uiBinnedIonMasses);
	}

	int iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	preXcorr();

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanSize; i++)
	{
		int iThreadId = omp_get_thread_num();
		struct Query *pQuery = new Query();
		if (!CometSearchMod::Preprocess(pQuery, vpAllMS2Scans.at(i), vpdTmpRawData.at(iThreadId),
										vpdTmpFastXcorrData.at(iThreadId), vpdTmpCorrelationData.at(iThreadId),
										vpdTmpSmoothedSpectrum.at(iThreadId), vpdTmpPeakExtracted.at(iThreadId)))
		{
			cout << "Error Post Xcorr." << endl;
			exit(1);
		}
		else
		{
			if (vpAllMS2Scans.at(i)->pQuery != NULL)
			{
				delete vpAllMS2Scans.at(i)->pQuery;
				vpAllMS2Scans.at(i)->pQuery = NULL;
			}
			vpAllMS2Scans.at(i)->pQuery = pQuery;
			for (int j = 0; j < (int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++)
			{
				double dXcorr = 0;
				PeptideUnit *pepUnit = vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j);
				CometSearchMod::ScorePeptidesSIPNoCancelOut(pepUnit->vvdYionMass, pepUnit->vvdYionProb,
															pepUnit->vvdBionMass, pepUnit->vvdBionProb, vpAllMS2Scans.at(i),
															vvpbDuplFragmentGlobal.at(iThreadId), vvdBinnedIonMassesGlobal.at(iThreadId),
															vvdBinGlobal.at(iThreadId), dXcorr);
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores[1] = dXcorr;
			}
			delete pQuery;
			vpAllMS2Scans.at(i)->pQuery = NULL;
		}
	}

	postXcorr();

	for (int i = 0; i < num_max_threads; ++i)
	{
		delete[] vpdTmpRawData.at(i);
		delete[] vpdTmpFastXcorrData.at(i);
		delete[] vpdTmpCorrelationData.at(i);
		delete[] vpdTmpSmoothedSpectrum.at(i);
		delete[] vpdTmpPeakExtracted.at(i);
		delete[] vpbDuplFragment.at(i);
		delete[] v_pdAAforward.at(i);
		delete[] v_pdAAreverse.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			for (int j = 0; j < 9; j++)
			{
				delete[] v_uiBinnedIonMasses.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMasses.at(i)[ii];
		}
		delete[] v_uiBinnedIonMasses.at(i);
	}
	++PeptideUnit::iNumScores;
}

double roundMy(double f, int precision)
{
	if (f == 0.0f)
		return +0.0f;

	double multiplier = pow(10.0, (double)precision); // moves f over <precision> decimal places
	f *= multiplier;
	f = floor(f + 0.5f);
	return f / multiplier;
}

void MS2ScanVector::postProcessAllMs2MvhSip()
{
	int num_threads = omp_get_max_threads();
	vector<multimap<double, double> *> vpIntenSortedPeakPreData;
	for (int j = 0; j < num_threads; ++j)
	{
		vpIntenSortedPeakPreData.push_back(new multimap<double, double>());
	}
	this->preMvh();
	int iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	double totalPeakSpace = ProNovoConfig::maxObservedMz - ProNovoConfig::minObservedMz;
	int totalPeakBins = (int)roundMy(totalPeakSpace / (ProNovoConfig::getMassAccuracyFragmentIon() * 2.0f), 0);
	MVH::initialLnTable(totalPeakBins);

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanSize; i++)
	{
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans.at(i)->preprocessMvh(vpIntenSortedPeakPreData.at(iThreadId));
		if (!vpAllMS2Scans.at(i)->bSkip)
		{
			for (int j = 0; j < (int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++)
			{
				double dMvh = 0;
				PeptideUnit *pepUnit = vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j);
				MVH::ScoreSequenceVsSpectrumSIP(pepUnit->sPeptideForScoring, pepUnit->iMeasuredParentCharge, vpAllMS2Scans.at(i),
												psequenceIonMasses.at(iThreadId), pepUnit->vvdYionMass, pepUnit->vvdYionProb,
												pepUnit->vvdBionMass, pepUnit->vvdBionProb, dMvh, pSeqs.at(iThreadId));
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores[2] = dMvh;
			}
		}
		else
		{
			for (int j = 0; j < (int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++)
			{
				double dMvh = 0;
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores[2] = dMvh;
			}
		}
		if (vpAllMS2Scans.at(i)->pPeakList != NULL)
		{
			delete vpAllMS2Scans.at(i)->pPeakList;
			vpAllMS2Scans.at(i)->pPeakList = NULL;
		}
		if (vpAllMS2Scans.at(i)->intenClassCounts != NULL)
		{
			delete vpAllMS2Scans.at(i)->intenClassCounts;
			vpAllMS2Scans.at(i)->intenClassCounts = NULL;
		}
	}

	for (int j = 0; j < num_threads; ++j)
	{
		delete vpIntenSortedPeakPreData.at(j);
	}
	vpIntenSortedPeakPreData.clear();
	this->postMvh();
	MVH::destroyLnTable();
	++PeptideUnit::iNumScores;
}

bool MS2ScanVector::mylessScanId(MS2Scan *pMS2Scan1, MS2Scan *pMS2Scan2)
{
	return (pMS2Scan1->iScanId < pMS2Scan2->iScanId);
}

string MS2ScanVector::ParsePath(string sPath)
// return the filename without path
{
	string sTailFileName;
	size_t iPosition;
	iPosition = sPath.rfind(ProNovoConfig::getSeparator());
	if (iPosition == string::npos)
		sTailFileName = sPath;
	else
		sTailFileName = sPath.substr(iPosition + 1);
	return sTailFileName;
}

void MS2ScanVector::clearSearchResults()
{
	for (MS2Scan *scan : vpAllMS2Scans)
	{
		if (scan == nullptr)
		{
			continue;
		}
		for (PeptideUnit *peptide : scan->vpWeightSumTopPeptides)
		{
			delete peptide;
		}
		scan->vpWeightSumTopPeptides.clear();
		scan->vdWeightSumAllScores.clear();
		scan->vMassChargePeptidePtrTuples.clear();
	}
}

static string formatContextPeptide(char prefix, const string &peptide, char suffix)
{
	string value;
	if (prefix != '-')
	{
		value.push_back(prefix);
	}
	if (peptide.find('[') != string::npos || peptide.find(']') != string::npos)
	{
		value += peptide;
	}
	else
	{
		value.push_back('[');
		value += peptide;
		value.push_back(']');
	}
	if (suffix != '-')
	{
		value.push_back(suffix);
	}
	return value;
}

static string stripPeptideForFeatures(const string &peptide)
{
	string stripped;
	stripped.reserve(peptide.size());
	for (char ch : peptide)
	{
		if (std::isalpha(static_cast<unsigned char>(ch)))
		{
			stripped.push_back(ch);
		}
	}
	return stripped;
}

struct FragmentSeriesFeatures
{
	int matchedB = 0;
	int matchedY = 0;
	int maxConsecutiveB = 0;
	int maxConsecutiveY = 0;
};

static bool matchesMvhPeak(const MS2Scan *scan, double mz)
{
	if (scan == nullptr ||
		mz < scan->mzLowerBound || mz > scan->mzUpperBound)
	{
		return false;
	}
	const double tolerance = ProNovoConfig::getMassAccuracyFragmentIon();
	if (scan->pPeakList != nullptr)
	{
		const char peakClass = scan->pPeakList->findNear(mz, tolerance);
		return peakClass != scan->pPeakList->end() && peakClass > 0;
	}
	// MVH releases its temporary PeakList before PIN rows are collected.
	// Raxport keeps vdMZ sorted, so use the retained experimental peaks rather
	// than silently reporting every matched-ion feature as zero.
	const auto peak = lower_bound(
		scan->vdMZ.begin(), scan->vdMZ.end(), mz - tolerance);
	return peak != scan->vdMZ.end() && *peak < mz + tolerance;
}

static bool matchesIsotopeEnvelope(
	const MS2Scan *scan, const vector<double> &neutralMasses, int charge)
{
	if (charge <= 0)
	{
		return false;
	}
	for (double neutralMass : neutralMasses)
	{
		const double mz =
			(neutralMass + Proton * static_cast<double>(charge)) /
			static_cast<double>(charge);
		if (matchesMvhPeak(scan, mz))
		{
			return true;
		}
	}
	return false;
}

static int longestMatchedSeries(const vector<bool> &matches)
{
	int longest = 0;
	int current = 0;
	for (bool matched : matches)
	{
		current = matched ? current + 1 : 0;
		longest = max(longest, current);
	}
	return longest;
}

static FragmentSeriesFeatures calculateFragmentSeriesFeatures(
	const PeptideUnit *peptide, MS2Scan *scan)
{
	FragmentSeriesFeatures features;
	if (peptide == nullptr || scan == nullptr)
	{
		return features;
	}

	vector<double> sequenceIons;
	vector<double> aaForward;
	vector<double> aaReverse;
	vector<char> residues;
	if (!MVH::CalculateSequenceIons(
			peptide->sPeptideForScoring,
			peptide->iMeasuredParentCharge,
			MVH::bUseSmartPlusThreeModel,
			&sequenceIons,
			&aaForward,
			&aaReverse,
			&residues) ||
		residues.size() < 2)
	{
		return features;
	}

	const int peptideLength = static_cast<int>(residues.size());
	vector<bool> bMatches(static_cast<size_t>(peptideLength - 1), false);
	vector<bool> yMatches(static_cast<size_t>(peptideLength - 1), false);
	const bool sipEnvelopes =
		ProNovoConfig::getSearchType() == "SIP" &&
		peptide->vvdBionMass.size() == bMatches.size() &&
		peptide->vvdYionMass.size() == yMatches.size();
	const auto matchB = [&](size_t index, int charge, double regularMz) {
		return sipEnvelopes
			? matchesIsotopeEnvelope(
				scan, peptide->vvdBionMass[index], charge)
			: matchesMvhPeak(scan, regularMz);
	};
	const auto matchY = [&](size_t index, int charge, double regularMz) {
		return sipEnvelopes
			? matchesIsotopeEnvelope(
				scan, peptide->vvdYionMass[index], charge)
			: matchesMvhPeak(scan, regularMz);
	};

	if (peptide->iMeasuredParentCharge > 2 && MVH::bUseSmartPlusThreeModel)
	{
		int totalStrongBasic = 0;
		int totalWeakBasic = 0;
		for (char residue : residues)
		{
			if (residue == 'R' || residue == 'K' || residue == 'H')
				++totalStrongBasic;
			else if (residue == 'Q' || residue == 'N')
				++totalWeakBasic;
		}
		const int totalBasicity =
			totalStrongBasic * 4 + totalWeakBasic * 2 + peptideLength - 2;
		const double chargeDenominator =
			static_cast<double>(peptide->iMeasuredParentCharge - 1);
		int bStrongBasic = 0;
		int bWeakBasic = 0;
		for (int cleavage = 1; cleavage < peptideLength; ++cleavage)
		{
			const char previousResidue =
				residues[static_cast<size_t>(cleavage - 1)];
			if (previousResidue == 'R' || previousResidue == 'K' ||
				previousResidue == 'H')
			{
				++bStrongBasic;
			}
			else if (previousResidue == 'Q' || previousResidue == 'N')
			{
				++bWeakBasic;
			}

			const int bScore =
				bStrongBasic * 4 + bWeakBasic * 2 + cleavage;
			const double basicityRatio =
				static_cast<double>(bScore) /
				static_cast<double>(totalBasicity);
			int first = 1;
			int last = peptide->iMeasuredParentCharge - 1;
			while (first < last)
			{
				const int middle = first + (last - first) / 2;
				if (static_cast<double>(middle) / chargeDenominator <=
					basicityRatio)
				{
					first = middle + 1;
				}
				else
				{
					last = middle;
				}
			}
			const int bCharge = first;
			const int yCharge =
				peptide->iMeasuredParentCharge - bCharge;
			const size_t bIndex =
				static_cast<size_t>(cleavage - 1);
			const size_t yIndex =
				static_cast<size_t>(peptideLength - 1 - cleavage);
			const double bMz =
				(aaForward[bIndex] + Proton * bCharge) / bCharge;
			const double yMz =
				(aaReverse[yIndex] + Proton * yCharge) / yCharge;
			bMatches[bIndex] = matchB(bIndex, bCharge, bMz);
			yMatches[yIndex] = matchY(yIndex, yCharge, yMz);
		}
	}
	else
	{
		for (size_t ionIndex = 0;
			ionIndex < bMatches.size(); ++ionIndex)
		{
			const int maxCharge =
				peptide->iMeasuredParentCharge > 2
				? peptide->iMeasuredParentCharge - 1 : 1;
			for (int charge = 1; charge <= maxCharge; ++charge)
			{
				const double bMz =
					(aaForward[ionIndex] + Proton * charge) / charge;
				const double yMz =
					(aaReverse[ionIndex] + Proton * charge) / charge;
				bMatches[ionIndex] =
					bMatches[ionIndex] ||
					matchB(ionIndex, charge, bMz);
				yMatches[ionIndex] =
					yMatches[ionIndex] ||
					matchY(ionIndex, charge, yMz);
			}
		}
	}

	features.matchedB = static_cast<int>(
		count(bMatches.begin(), bMatches.end(), true));
	features.matchedY = static_cast<int>(
		count(yMatches.begin(), yMatches.end(), true));
	features.maxConsecutiveB = longestMatchedSeries(bMatches);
	features.maxConsecutiveY = longestMatchedSeries(yMatches);
	return features;
}

void MS2ScanVector::appendScoredPsmRows(vector<ScoredPsmRow> &rows, bool isDecoy, int topKeep, double ms2IsotopicAbundancePct) const
{
	vector<MS2Scan *> orderedScans = vpAllMS2Scans;
	sort(orderedScans.begin(), orderedScans.end(), mylessScanId);
	const int keep = topKeep > 0 ? topKeep : ProNovoConfig::INTTOPKEEP;
	vector<vector<ScoredPsmRow>> rowsByScan(orderedScans.size());
	const long long scanCount = static_cast<long long>(orderedScans.size());
#pragma omp parallel for schedule(guided)
	for (long long scanIndex = 0; scanIndex < scanCount; ++scanIndex)
	{
		MS2Scan *scan = orderedScans[static_cast<size_t>(scanIndex)];
		if (scan == nullptr || scan->vpWeightSumTopPeptides.empty())
		{
			continue;
		}
		vector<ScoredPsmRow> &scanRows =
			rowsByScan[static_cast<size_t>(scanIndex)];
		scanRows.reserve(scan->vpWeightSumTopPeptides.size());
		for (int j = 0; j < static_cast<int>(scan->vpWeightSumTopPeptides.size()); ++j)
		{
			if (!scan->isAnyScoreInTopN(j, keep))
			{
				continue;
			}
			PeptideUnit *peptide = scan->vpWeightSumTopPeptides.at(j);
			if (peptide == nullptr)
			{
				continue;
			}
			ScoredPsmRow row;
			row.scanNumber = scan->iScanId;
			row.parentCharge = peptide->iMeasuredParentCharge;
			row.precursorScanNumber = scan->iParentScanID;
			row.isolationWindowCenterMZ = scan->dParentMZ;
			row.measuredParentMass = peptide->dMeasuredParentMass;
			row.calculatedParentMass = peptide->dPepNeutralMass;
			row.precursorNeutronMass = peptide->dPrecursorNeutronMass;
			row.scanType = scan->getScanType();
			row.searchName = ProNovoConfig::getSearchName();
			row.ms2IsotopicAbundancePct = ms2IsotopicAbundancePct;
			row.retentionTime = std::atof(scan->getRTime().c_str());
			row.wdpScore = static_cast<float>(peptide->vdScores[0]);
			row.xcorrScore = static_cast<float>(peptide->vdScores[1]);
			row.mvhScore = static_cast<float>(peptide->vdScores[2]);
			row.rank = static_cast<int>(peptide->vdRank[0]);
			row.ddaResidualRank = peptide->iDdaResidualRank;
			row.ddaResidualScore =
				static_cast<float>(peptide->dDdaResidualScore);
			const FragmentSeriesFeatures fragmentFeatures =
				calculateFragmentSeriesFeatures(peptide, scan);
			row.matchedBIons = fragmentFeatures.matchedB;
			row.matchedYIons = fragmentFeatures.matchedY;
			row.maxConsecutiveBIons =
				fragmentFeatures.maxConsecutiveB;
			row.maxConsecutiveYIons =
				fragmentFeatures.maxConsecutiveY;
			row.identifiedPeptide = formatContextPeptide(peptide->cIdentifyPrefix, peptide->sIdentifiedPeptide, peptide->cIdentifySuffix);
			row.originalPeptide = formatContextPeptide(peptide->cOriginalPrefix, peptide->sOriginalPeptide, peptide->cOriginalSuffix);
			row.nakedPeptide = stripPeptideForFeatures(peptide->sIdentifiedPeptide);
			row.proteinNames = "{" + peptide->sProteinNames + "}";
			row.isDecoy = isDecoy;
			scanRows.push_back(std::move(row));
		}
	}

	size_t appendedCount = 0;
	for (const vector<ScoredPsmRow> &scanRows : rowsByScan)
		appendedCount += scanRows.size();
	rows.reserve(rows.size() + appendedCount);
	for (vector<ScoredPsmRow> &scanRows : rowsByScan)
	{
		for (ScoredPsmRow &row : scanRows)
			rows.push_back(std::move(row));
	}
}

size_t MS2ScanVector::matchScoredPsmPrecursors(
	vector<ScoredPsmRow> &rows) const
{
	long long matchedCount = 0;
#pragma omp parallel for schedule(guided) reduction(+ : matchedCount)
	for (long long i = 0; i < static_cast<long long>(rows.size()); ++i)
	{
		ScoredPsmRow &row = rows[static_cast<size_t>(i)];
		const int parentScanNumber = row.precursorScanNumber;
		const sipros::RaxportPrecursorMatch match =
			sipros::findRaxportPrecursorMatch(
				raxportMs1Data,
				parentScanNumber,
				row.retentionTime,
				row.parentCharge,
				row.calculatedParentMass,
				row.precursorNeutronMass,
				ProNovoConfig::getParentMassWindows(),
				ProNovoConfig::getMassAccuracyParentIon(),
				raxportReadOptions.precursorMatchScanRadius);
		row.precursorRtDiffSeconds = match.rtDiffSeconds;
		if (!match.found)
		{
			continue;
		}
		row.precursorScanNumber = match.ms1ScanNumber;
		row.measuredParentMass = match.observedNeutralMass;
		++matchedCount;
	}
	return static_cast<size_t>(matchedCount);
}

void MS2ScanVector::preMvh()
{
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i)
	{
		_ppdAAforward.push_back(new vector<double>());
		_ppdAAreverse.push_back(new vector<double>());
		psequenceIonMasses.push_back(new vector<double>());
		pSeqs.push_back(new vector<char>());
	}
}

void MS2ScanVector::postMvh()
{
	for (int i = 0; i < num_max_threads; ++i)
	{
		delete _ppdAAforward.at(i);
		delete _ppdAAreverse.at(i);
		delete psequenceIonMasses.at(i);
		delete pSeqs.at(i);
	}
	_ppdAAforward.clear();
	_ppdAAreverse.clear();
	psequenceIonMasses.clear();
	pSeqs.clear();
}

void MS2ScanVector::preXcorr()
{
	// MH: Must be equal to largest possible array
	int iArraySizeScore = (int)((ProNovoConfig::dMaxMS2ScanMass * 2 + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i)
	{
		vpbDuplFragmentGlobal.push_back(new bool[iArraySizeScore]());
		v_pdAAforwardGlobal.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverseGlobal.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int ***_uiBinnedIonMasses = new unsigned int **[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			_uiBinnedIonMasses[ii] = new unsigned int *[9]();
			for (int j = 0; j < 9; j++)
			{
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMassesGlobal.push_back(_uiBinnedIonMasses);
	}
	// SIP mode variable
	for (int i = 0; i < num_max_threads; ++i)
	{
		vvpbDuplFragmentGlobal.push_back(vector<unsigned char>());
		vvpbDuplFragmentGlobal.back().clear();
		vvpbDuplFragmentGlobal.back().resize(iArraySizeScore, false);
		vvdBinnedIonMassesGlobal.push_back(vector<double>());
		vvdBinnedIonMassesGlobal.back().clear();
		vvdBinnedIonMassesGlobal.back().resize(iArraySizeScore, 0);
		vvdBinGlobal.push_back(vector<int>());
	}
}

void MS2ScanVector::postXcorr()
{
	for (int i = 0; i < num_max_threads; ++i)
	{
		delete[] vpbDuplFragmentGlobal.at(i);
		delete[] v_pdAAforwardGlobal.at(i);
		delete[] v_pdAAreverseGlobal.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++)
		{
			for (int j = 0; j < 9; j++)
			{
				delete[] v_uiBinnedIonMassesGlobal.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMassesGlobal.at(i)[ii];
		}
		delete[] v_uiBinnedIonMassesGlobal.at(i);
	}
	vpbDuplFragmentGlobal.clear();
	v_pdAAforwardGlobal.clear();
	v_pdAAreverseGlobal.clear();
	v_uiBinnedIonMassesGlobal.clear();
	// SIP mode variable
	vvpbDuplFragmentGlobal.clear();
	vvdBinnedIonMassesGlobal.clear();
	vvdBinGlobal.clear();
}
