#include "ms2scanvector.h"
#include "RaxportHdf5Reader.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <utility>

MS2ScanVector::MS2ScanVector(const string &sScanFilenameInput, const string &sOutputDirectory,
							 const string &sConfigFilename, bool bScreenOutput)
{
	unsigned int n;
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;
	sScanFilename = sScanFilenameInput;
	sConfigFile = sConfigFilename;
	// mass_w = ProNovoConfig::getParentMassWindows();
	setOutputFile(sScanFilenameInput, sOutputDirectory);
	this->bScreenOutput = bScreenOutput;
	for (n = 0; n < vsSingleResidueNames.size(); ++n)
		mapResidueMass[vsSingleResidueNames[n][0]] = vdSingleResidueMasses[n];

	iOpenMPTaskNum = 0;
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
	if (!sipros::readRaxportHdf5Scans(sScanFilename, vpAllMS2Scans, nullptr, error))
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
		const double primaryChargedMass = scan->dParentMZ * scan->iParentChargeState;
		if (scan->iParentChargeState > 0 && scan->dParentMZ > 0.0)
		{
			const double parentNeutralMass = primaryChargedMass -
				static_cast<double>(scan->iParentChargeState) * ProNovoConfig::getProtonMass();
			scan->dParentNeutralMass = parentNeutralMass;
			scan->dParentMass = primaryChargedMass;
			vAllPrecursorMassChargeMS2ScanPtrTuples.push_back({parentNeutralMass, scan->iParentChargeState, scan});
		}

		const size_t nCandidates = std::min(scan->dParentMZs.size(), scan->iParentChargeStates.size());
		for (size_t j = 0; j < nCandidates; ++j)
		{
			const int charge = scan->iParentChargeStates[j];
			const double mz = scan->dParentMZs[j];
			if (charge <= 0 || mz <= 0.0)
			{
				continue;
			}
			const double candidateChargedMass = mz * charge;
			if (scan->iParentChargeState > 0 &&
				std::abs(candidateChargedMass - primaryChargedMass) <= ProNovoConfig::getMassAccuracyParentIon())
			{
				continue;
			}
			const double parentNeutralMass = candidateChargedMass -
				static_cast<double>(charge) * ProNovoConfig::getProtonMass();
			vAllPrecursorMassChargeMS2ScanPtrTuples.push_back({parentNeutralMass, charge, scan});
			if (parentNeutralMass > scan->dParentNeutralMass)
			{
				scan->dParentNeutralMass = parentNeutralMass;
			}
			if (candidateChargedMass > scan->dParentMass)
			{
				scan->dParentMass = candidateChargedMass;
			}
		}

		if (scan->iParentChargeState <= 0 && nCandidates > 0)
		{
			// Keep the historical wide-window fallback charge for downstream preprocessing.
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
#ifdef Ticktock
	TOCK1ST(loadHdf5file);
#endif
	return true;
}

bool MS2ScanVector::loadMassData()
{
	CLOCKSTART;

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
	ProNovoConfig::getSetFileNameSuffix() = fileNameSuffix;

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
	if (bReVal)
	{
		std::cout << "\nload Raxport HDF5 mass data done.\n" << std::endl;
	}
	CLOCKSTOP;
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
	CLOCKSTART;

	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();
	if (bScreenOutput)
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans " << endl;

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

	if (bScreenOutput)
	{
		cout << "Preprocessing Done." << endl;
	}

	CLOCKSTOP;
}

void MS2ScanVector::preProcessAllMs2WdpSip()
{
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();
	if (bScreenOutput)
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans " << endl;

#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
		vpAllMS2Scans.at(i)->preprocess();
}

void MS2ScanVector::GetAllRangeFromMass(double dPeptideMass, vector<std::pair<int, int>> &vpPeptideMassRanges)
// all ranges of MS2 scans are stored in  vpPeptideMassWindows
{
	int i;
	pair<int, int> pairMS2Range;
	pair<int, int> lastPairRange(-100, -100);
	vector<pair<double, double>> vpPeptideMassWindows;
	vpPeptideMassWindows.clear();
	vpPeptideMassRanges.clear();
	ProNovoConfig::getPeptideMassWindows(dPeptideMass, vpPeptideMassWindows);
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
	int low = 0, high = vpPrecursorMasses.size() - 1, mid;
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

	GetAllRangeFromMass(currentPeptide->getPeptideMass(), vpPeptideMassRanges);

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

void MS2ScanVector::processPeptideArrayWdpSip(vector<Peptide *> &vpPeptideArray)
{
	int i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = (int)vpPeptideArray.size();

	//    for (int i=0; i< (int) vpPeptideArray.size(); i++)
	//      cout<<vpPeptideArray.at(i)->getPeptideSeq() <<"\t"
	//	  <<vpPeptideArray.at(i)->getOriginalPeptideSeq()<<"\t"
	//	  <<vpPeptideArray.at(i)->getProteinName()<<"\t"
	//	  <<vpPeptideArray.at(i)->getBeginPosProtein()<<"\t"
	//	  <<vpPeptideArray.at(i)->getPeptideMass()<<endl;

	// cout<<"calculating fragments of "<<vpPeptideArray.size()<<"  peptides"<<endl;
#pragma omp parallel for shared(vpPeptideArray) private(i) \
	schedule(guided)

	for (i = 0; i < iPeptideArraySize; i++)
		vpPeptideArray[i]->preprocessing(vpAllMS2Scans.at(0)->isMS2HighRes, mapResidueMass);
	// vpPeptideArray[i]->calculateExpectedFragments(mapResidueMass);

	// cout<<"scoring "<<vpPeptideArray.size()<<"  peptides"<<endl;
	//  every MS2 scans scores their matched peptides

	iScanSize = (int)vpAllMS2Scans.size();
#pragma omp parallel for schedule(guided)

	for (i = 0; i < iScanSize; i++)
		vpAllMS2Scans[i]->scorePeptides();

	// free memory of all peptide objects
	for (i = 0; i < (int)vpPeptideArray.size(); i++)
		delete vpPeptideArray[i];

	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::processPeptideArrayMvh(vector<Peptide *> &vpPeptideArray)
{
	int i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = (int)vpPeptideArray.size();

#pragma omp parallel for shared(vpPeptideArray) private(i) \
	schedule(guided)

	for (i = 0; i < iPeptideArraySize; i++)
	{
		vpPeptideArray.at(i)->preprocessingMVH();
	}

	// every MS2 scans scores their matched peptides

	iScanSize = (int)vpAllMS2Scans.size();
#pragma omp parallel for schedule(guided)

	for (i = 0; i < iScanSize; i++)
	{
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans[i]->scorePeptidesMVH(psequenceIonMasses.at(iThreadId), _ppdAAforward.at(iThreadId),
										   _ppdAAreverse.at(iThreadId), pSeqs.at(iThreadId));
	}

	// free memory of all peptide objects
	for (i = 0; i < (int)vpPeptideArray.size(); i++)
		delete vpPeptideArray.at(i);

	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::processPeptideArrayMvhTask(vector<Peptide *> &vpPeptideArray, omp_lock_t *pLck)
{
	int i, j, k, iPeptideArraySize;
	vector<pair<int, int>> vpPeptideMassRanges;
	pair<int, int> pairMS2Range;
	bool bAssigned = false;
	bool bMerged = false, bScored = false, bProcessed = false;
	double dMvh = 0;
	int iThreadId = omp_get_thread_num();

	iPeptideArraySize = vpPeptideArray.size();

	MS2Scan *scanPtr;
	int precursorCharge;
	double precursorMass;
	for (i = 0; i < iPeptideArraySize; i++)
	{
		bAssigned = false;
		bProcessed = false;
		// assign peptide to scan
		GetAllRangeFromMass(vpPeptideArray.at(i)->getPeptideMass(), vpPeptideMassRanges);
		for (j = 0; j < (int)vpPeptideMassRanges.size(); j++)
		{
			pairMS2Range = vpPeptideMassRanges.at(j);
			if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
			{
				bAssigned = true;
				if (bAssigned && !bProcessed)
				{
					vpPeptideArray.at(i)->preprocessingMVH();
					bProcessed = true;
				}
				// calculate the score
				for (k = pairMS2Range.first; k <= pairMS2Range.second; ++k)
				{
					tie(precursorMass, precursorCharge, scanPtr) = vAllPrecursorMassChargeMS2ScanPtrTuples[k];
					if (scanPtr->bSkip)
					{
						continue;
					}
					omp_set_lock(&(pLck[k]));
					bMerged = scanPtr->mergePeptide(scanPtr->vpWeightSumTopPeptides,
													vpPeptideArray.at(i)->getPeptideForScoring(), vpPeptideArray.at(i)->getProteinName());
					omp_unset_lock(&(pLck[k]));
					// not merged then score
					if (!bMerged)
					{
						bScored = MVH::ScoreSequenceVsSpectrum(vpPeptideArray.at(i)->sNeutralLossPeptide, precursorCharge,
															   scanPtr, psequenceIonMasses.at(iThreadId), _ppdAAforward.at(iThreadId),
															   _ppdAAreverse.at(iThreadId), dMvh, pSeqs.at(iThreadId));
						if (bScored)
						{
							omp_set_lock(&(pLck[k]));
							scanPtr->saveScore(dMvh, {precursorMass, precursorCharge, vpPeptideArray.at(i)},
											   scanPtr->vpWeightSumTopPeptides,
											   "MVH", 2);
							omp_unset_lock(&(pLck[k]));
						}
					}
				}
			}
		}
	}

	// free memory of all peptide objects
	for (i = 0; i < iPeptideArraySize; i++)
	{
		delete vpPeptideArray.at(i);
	}

	// empty peptide array
	vpPeptideArray.clear();

	omp_set_lock(&lckOpenMpTaskNum);
	iOpenMPTaskNum--;
	if (iOpenMPTaskNum < TASKRESUME_SIZE)
	{
		if (omp_test_lock(&lckOpenMpTaskNumHalfed) == 0)
		{
			// lock is set'
			omp_unset_lock(&lckOpenMpTaskNumHalfed);
			// cout << endl << "Resume" << endl;
		}
		else
		{
			omp_unset_lock(&lckOpenMpTaskNumHalfed);
		}
	}
	omp_unset_lock(&lckOpenMpTaskNum);
}

void MS2ScanVector::searchDatabaseMvh()
{
	CLOCKSTART;

	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide *currentPeptide;
	myProteinDatabase.loadDatabase();
	this->preMvh();
	if (myProteinDatabase.getFirstProtein())
	{
		currentPeptide = new Peptide;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide))
		{
			// assign the pointers of peptides to appropriete MS2Scans
			if (assignPeptides2Scans(currentPeptide))
			{
				// save the new peptide to the array
				vpPeptideArray.push_back(currentPeptide);
				if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass)
				{
					ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
				}
			}
			else
			{
				delete currentPeptide;
			}
			// create a new peptide for the next iteration
			currentPeptide = new Peptide;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE)
				processPeptideArrayMvh(vpPeptideArray);
		}
		// the last peptide object is an empty object and need to be deleted
		delete currentPeptide;
		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
		//    cout<<vpPeptideArray.size()<<endl;
		if (!vpPeptideArray.empty())
			processPeptideArrayMvh(vpPeptideArray);
	}
	CLOCKSTOP;

	this->postMvh();
	MVH::destroyLnTable();
	PeptideUnit::iNumScores = 1;
	cout << "MVH search done.\n"
		 << endl;
}

void MS2ScanVector::searchDatabaseWdpSip()
{
	CLOCKSTART;
	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide *currentPeptide;
	myProteinDatabase.loadDatabase();

	if (myProteinDatabase.getFirstProtein())
	{
		currentPeptide = new Peptide;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide))
		{

			// assign the pointers of peptides to appropriete MS2Scans
			if (assignPeptides2Scans(currentPeptide))
			{
				// save the new peptide to the array
				vpPeptideArray.push_back(currentPeptide);
			}
			else
			{
				// not assigned, delete
				delete currentPeptide;
			}
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
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
		if (!vpPeptideArray.empty())
		{
			processPeptideArrayWdpSip(vpPeptideArray);
		}

		PeptideUnit::iNumScores = 1;
	}
	cout << "\nWDP search done.\n"
		 << endl;
	CLOCKSTOP;
}

void MS2ScanVector::searchDatabaseMvhTask()
{
	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide *currentPeptide;
	myProteinDatabase.loadDatabase();
	this->preMvh();

	int iScanNum = this->vpAllMS2Scans.size();
	omp_lock_t *pLck = new omp_lock_t[iScanNum];
	iOpenMPTaskNum = 0;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanNum; ++i)
	{
		omp_init_lock(&(pLck[i]));
	}
	omp_init_lock(&lckOpenMpTaskNum);
	omp_init_lock(&lckOpenMpTaskNumHalfed);
	// omp_set_lock(&lckOpenMpTaskNumHalfed);
#pragma omp parallel
	{
#pragma omp single
		{
			if (myProteinDatabase.getFirstProtein())
			{
				currentPeptide = new Peptide();
				// get one peptide from the database at a time, until there is no more peptide
				while (myProteinDatabase.getNextPeptide(currentPeptide))
				{
					if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass)
					{
						ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
					}
					vpPeptideArray.push_back(currentPeptide);
					if (vpPeptideArray.size() >= TASKPEPTIDE_ARRAY_SIZE)
					{
						omp_set_lock(&lckOpenMpTaskNum);
						iOpenMPTaskNum++;
						omp_unset_lock(&lckOpenMpTaskNum);
#pragma omp task firstprivate(vpPeptideArray)
						{
							processPeptideArrayMvhTask(vpPeptideArray, pLck);
						}

						vpPeptideArray.clear();

						if (iOpenMPTaskNum >= TASKWAIT_SIZE)
						{
							// cout << endl << "Stop" << endl;
							omp_set_lock(&lckOpenMpTaskNumHalfed);
						}
					}
					currentPeptide = new Peptide();
				}
				delete currentPeptide;
			}
			if (!vpPeptideArray.empty())
			{
#pragma omp task firstprivate(vpPeptideArray)
				{
					processPeptideArrayMvhTask(vpPeptideArray, pLck);
				}
				vpPeptideArray.clear();
			}
#pragma omp taskwait
		}
	}

	omp_destroy_lock(&lckOpenMpTaskNumHalfed);
	omp_destroy_lock(&lckOpenMpTaskNum);

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanNum; ++i)
	{
		omp_destroy_lock(&(pLck[i]));
	}

	this->postMvh();
	MVH::destroyLnTable();
	PeptideUnit::iNumScores = 1;
	cout << "MVH search done.\n"
		 << endl;
}

void MS2ScanVector::startProcessingMvh()
{
	// Preprocessing all MS2 scans by mult-threading
	preProcessAllMs2Mvh();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseMvh();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2WdpXcorr();

}

void MS2ScanVector::startProcessingMvhTask()
{
	// Preprocessing all MS2 scans by mult-threading
	preProcessAllMs2Mvh();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseMvhTask();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2WdpXcorr();

}

void MS2ScanVector::startProcessingWdpSip()
{

	// Preprocessing all MS2 scans by multi-threading
	preProcessAllMs2WdpSip();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseWdpSip();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2MvhXcorr();

}

void MS2ScanVector::postProcessAllMs2WdpXcorr()
{
	CLOCKSTART;
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	postProcessAllMs2Xcorr();
	cout << "\nXcorr search done." << endl;

	postProcessAllMs2Wdp();
	cout << "\nWDP search done.\n"
		 << endl;

#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
	{
		vpAllMS2Scans.at(i)->scoreFeatureCalculation();
	}
	CLOCKSTOP;
}

void MS2ScanVector::postProcessAllMs2MvhXcorr()
{
	CLOCKSTART;
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();

	postProcessAllMs2XcorrSip();
	cout << "\nXcorr search done.\n" << endl;

	postProcessAllMs2MvhSip();
	cout << "MVH search done.\n" << endl;

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
				dWeightSum = vpAllMS2Scans.at(i)->scoreWeightSum(
					&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
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
											  vpAllMS2Scans.at(i), v_uiBinnedIonMasses.at(iThreadId), dXcorr);
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

void MS2ScanVector::appendScoredPsmRows(vector<ScoredPsmRow> &rows, bool isDecoy, int topKeep) const
{
	vector<MS2Scan *> orderedScans = vpAllMS2Scans;
	sort(orderedScans.begin(), orderedScans.end(), mylessScanId);
	const int keep = topKeep > 0 ? topKeep : ProNovoConfig::INTTOPKEEP;
	for (MS2Scan *scan : orderedScans)
	{
		if (scan == nullptr || scan->vpWeightSumTopPeptides.empty())
		{
			continue;
		}
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
			row.scanType = scan->getScanType();
			row.searchName = ProNovoConfig::getSearchName();
			row.retentionTime = static_cast<float>(std::atof(scan->getRTime().c_str()));
			row.wdpScore = static_cast<float>(peptide->vdScores[0]);
			row.xcorrScore = static_cast<float>(peptide->vdScores[1]);
			row.mvhScore = static_cast<float>(peptide->vdScores[2]);
			row.rank = static_cast<int>(peptide->vdRank[0]);
			row.identifiedPeptide = formatContextPeptide(peptide->cIdentifyPrefix, peptide->sIdentifiedPeptide, peptide->cIdentifySuffix);
			row.originalPeptide = formatContextPeptide(peptide->cOriginalPrefix, peptide->sOriginalPeptide, peptide->cOriginalSuffix);
			row.nakedPeptide = stripPeptideForFeatures(peptide->sIdentifiedPeptide);
			row.proteinNames = "{" + peptide->sProteinNames + "}";
			row.isDecoy = isDecoy;
			rows.push_back(std::move(row));
		}
	}
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
		vvpbDuplFragmentGlobal.push_back(vector<bool>());
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
