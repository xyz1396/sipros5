#include "proNovoConfig.h"

#include <cctype>
#include <cmath>
#include <stdexcept>

#include "averagine.h"

string ProNovoConfig::sFilename = "SiprosConfig.cfg";

void ProNovoConfig::setFASTAfilename(const string &fastaFilename)
{
	sFASTAFilename = fastaFilename;
}

void ProNovoConfig::setMassAccuracy(double parentIonToleranceDa, double fragmentIonToleranceDa)
{
	dMassAccuracyParentIon = parentIonToleranceDa;
	dMassAccuracyFragmentIon = fragmentIonToleranceDa;
}

#if _WIN32
string ProNovoConfig::sWorkingDirectory = ".\\";
#else
string ProNovoConfig::sWorkingDirectory = ".";
#endif

// variables from the Peptide_Identification element
string ProNovoConfig::sFASTAFilename = "";
string ProNovoConfig::sSearchType = "Regular";
string ProNovoConfig::sSearchName = "Null";
string ProNovoConfig::sFragmentationMethod = "CID";

int ProNovoConfig::iMaxPTMcount = 0;

int ProNovoConfig::iMinPeptideLength = 6;
int ProNovoConfig::iMaxPeptideLength = 60;

string ProNovoConfig::sCleavageAfterResidues = "KR";
string ProNovoConfig::sCleavageBeforeResidues = "ACDEFGHIJKLMNPQRSTVWXY";
int ProNovoConfig::iMaxMissedCleavages = 2;
bool ProNovoConfig::bTestStartRemoval = false;

double ProNovoConfig::dMassAccuracyParentIon = 0.05;
double ProNovoConfig::dMassAccuracyFragmentIon = 0.05;
vector<int> ProNovoConfig::viParentMassWindows;

ProNovoConfig *ProNovoConfig::ProNovoConfigSingleton = NULL;

vector<string> ProNovoConfig::vsSingleResidueNames;
vector<double> ProNovoConfig::vdSingleResidueMasses;

double ProNovoConfig::dTerminusMassN = 1.0072765;
double ProNovoConfig::dTerminusMassC = 17.00265;

string ProNovoConfig::sElementList = "";

map<string, string> ProNovoConfig::mapConfigKeyValues;
Isotopologue ProNovoConfig::configIsotopologue;
vector<vector<double>>
    ProNovoConfig::naturalAtomIsotopeProbabilities;

const vector<double> &
ProNovoConfig::getNaturalAtomIsotopeProbabilities(size_t atomIndex)
{
    static const vector<double> empty;
    if (atomIndex < naturalAtomIsotopeProbabilities.size())
        return naturalAtomIsotopeProbabilities[atomIndex];
    if (atomIndex <
        configIsotopologue.vAtomIsotopicDistribution.size())
    {
        return configIsotopologue
            .vAtomIsotopicDistribution[atomIndex].vProb;
    }
    return empty;
}

vector<pair<double, double>> ProNovoConfig::vpPeptideMassWindowOffset;

vector<pair<string, string>> ProNovoConfig::vpNeutralLossList;

//---------------Comet Begin---------------------
bool ProNovoConfig::bXcorrEnable = false;
Options ProNovoConfig::options;
double ProNovoConfig::dInverseBinWidth = 0;	  // this is used in BIN() many times so use inverse binWidth to do multiply vs. divide
double ProNovoConfig::dOneMinusBinOffset = 0; // this is used in BIN() many times so calculate once
IonInfo ProNovoConfig::ionInformation;
int ProNovoConfig::iXcorrProcessingOffset = 75;
PrecalcMasses ProNovoConfig::precalcMasses;
double ProNovoConfig::dMaxMS2ScanMass = 0;
double ProNovoConfig::dMaxPeptideMass = 0;
// map<char, double> ProNovoConfig::pdAAMassFragment;
AminoAcidMasses ProNovoConfig::pdAAMassFragment;
double ProNovoConfig::dHighResFragmentBinSize = 0.02;
double ProNovoConfig::dHighResFragmentBinStartOffset = 0;
double ProNovoConfig::dLowResFragmentBinSize = 1.0005;
double ProNovoConfig::dLowResFragmentBinStartOffset = 0.4;
double ProNovoConfig::dHighResInverseBinWidth = 1.0 / ProNovoConfig::dHighResFragmentBinSize;
double ProNovoConfig::dLowResInverseBinWidth = 1.0 / ProNovoConfig::dLowResFragmentBinSize;
double ProNovoConfig::dHighResOneMinusBinOffset = 1.0 - ProNovoConfig::dHighResFragmentBinStartOffset;
double ProNovoConfig::dLowResOneMinusBinOffset = 1.0 - ProNovoConfig::dLowResFragmentBinStartOffset;
int ProNovoConfig::iMaxPercusorCharge = 0;
//---------------Comet End-----------------------

//---------------Myrimatch Begin-----------------
bool ProNovoConfig::bMvhEnable = true;
double ProNovoConfig::ClassSizeMultiplier = 2;
int ProNovoConfig::NumIntensityClasses = 3;
int ProNovoConfig::minIntensityClassCount = int((pow(ClassSizeMultiplier, NumIntensityClasses) - 1) / (ClassSizeMultiplier - 1));
double ProNovoConfig::ticCutoffPercentage = 0.98;
int ProNovoConfig::MaxPeakCount = 300;
int ProNovoConfig::MinMatchedFragments = 5;
double ProNovoConfig::minObservedMz = numeric_limits<double>::max();
double ProNovoConfig::maxObservedMz = 0;
//---------------Myrimatch End-------------------

//---------------Sipros Score--------------------
bool ProNovoConfig::bWeightDotSumEnable = false;
bool ProNovoConfig::bLessIsotopicDistribution = false;
bool ProNovoConfig::bMultiScores = true;
string ProNovoConfig::sDecoyPrefix = "";
int ProNovoConfig::INTTOPKEEP = 10;
int ProNovoConfig::iRank = 0;
//---------------Sipros Score--------------------

double AminoAcidMasses::dNULL = -1;
double AminoAcidMasses::dERROR = -2;

string ProNovoConfig::SIPelement = "C";
int ProNovoConfig::SIPisotopeMassNumber = 13;
double ProNovoConfig::minValue = 0;
double ProNovoConfig::fold = 0;
double ProNovoConfig::deductionCoefficient = 0;
string ProNovoConfig::fileNameSuffix = "h5";
// carbon isotopic delta mass in default
double ProNovoConfig::neutronMass = 1.003355;

int ProNovoConfig::atomIndex(char sipAtom)
{
	const string atoms = "CHONPS";
	const char atom = static_cast<char>(std::toupper(static_cast<unsigned char>(sipAtom)));
	const size_t pos = atoms.find(atom);
	return pos == string::npos ? -1 : static_cast<int>(pos);
}

static bool supportedSipTarget(char sipAtom,
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

bool ProNovoConfig::refreshResidueDistributions(Isotopologue &iso)
{
	for (const auto &kv : iso.mResidueAtomicComposition)
	{
		IsotopeDistribution dist;
		if (!iso.computeIsotopicDistribution(kv.second, dist))
		{
			return false;
		}
		iso.vResidueIsotopicDistribution[kv.first] = dist;
	}
	return true;
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
    if (!supportedSipTarget(
            atom, expectedMassNumber, expectedIsotopeIndex, nominalShift) ||
        (isotopeMassNumber > 0 &&
         isotopeMassNumber != expectedMassNumber))
    {
        throw std::runtime_error(
            "Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
    }

    const int atomIx = atomIndex(atom);
    if (atomIx < 0 ||
        atomIx >= static_cast<int>(
            iso.vAtomIsotopicDistribution.size()))
    {
        throw std::runtime_error(
            "Configured SIP atom distribution is not available.");
    }
    const auto &distribution =
        iso.vAtomIsotopicDistribution[static_cast<size_t>(atomIx)];
    if (expectedIsotopeIndex >=
            static_cast<int>(distribution.vProb.size()) ||
        expectedIsotopeIndex >=
            static_cast<int>(distribution.vMass.size()) ||
        static_cast<int>(std::lround(
            distribution.vMass[
                static_cast<size_t>(expectedIsotopeIndex)])) !=
            expectedMassNumber)
    {
        throw std::runtime_error(
            "Configured distribution does not contain the requested "
            "supported SIP isotope.");
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
    if (!supportedSipTarget(
            atom, massNumber, expectedIsotopeIndex, nominalShift) ||
        isotopeIndex != expectedIsotopeIndex)
    {
        throw std::runtime_error(
            "Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
    }

    const int atomIx = atomIndex(atom);
    if (atomIx < 0 ||
        atomIx >= static_cast<int>(
            iso.vAtomIsotopicDistribution.size()) ||
        isotopeIndex >= static_cast<int>(
            iso.vAtomIsotopicDistribution[
                static_cast<size_t>(atomIx)].vProb.size()))
    {
        throw std::runtime_error(
            "Configured isotopic distribution for SIP target is not usable.");
    }

    auto &probabilities =
        iso.vAtomIsotopicDistribution[
            static_cast<size_t>(atomIx)].vProb;
    if (!averagine::changeAtomProbability(
            probabilities, atom, sipPct / 100.0))
    {
        throw std::runtime_error(
            "Configured isotopic distribution for SIP target is not usable.");
    }
    if (!refreshResidueDistributions(iso))
        throw std::runtime_error(
            "Failed to refresh residue isotopic distributions.");
}

double ProNovoConfig::getIsotopeAbundancePct(
    const Isotopologue &iso,
    char sipAtom,
    int isotopeIndex)
{
    int massNumber = 0;
    int expectedIsotopeIndex = -1;
    int nominalShift = 0;
    if (!supportedSipTarget(
            sipAtom, massNumber, expectedIsotopeIndex, nominalShift) ||
        isotopeIndex != expectedIsotopeIndex)
    {
        throw std::runtime_error(
            "Unsupported SIP isotope. Use C13,H2,N15,O18,S34.");
    }
    const int atomIx = atomIndex(sipAtom);
    if (atomIx < 0 ||
        atomIx >= static_cast<int>(
            iso.vAtomIsotopicDistribution.size()) ||
        isotopeIndex >= static_cast<int>(
            iso.vAtomIsotopicDistribution[
                static_cast<size_t>(atomIx)].vProb.size()))
    {
        throw std::runtime_error(
            "Configured isotope abundance is not available.");
    }
    return iso.vAtomIsotopicDistribution[
               static_cast<size_t>(atomIx)]
               .vProb[static_cast<size_t>(isotopeIndex)] *
           100.0;
}

bool ProNovoConfig::applySipAbundance(char sipAtom,
                                        double fraction)
{
    const char atom = static_cast<char>(
        std::toupper(static_cast<unsigned char>(sipAtom)));
    int massNumber = 0;
    int isotopeIndex = -1;
    int nominalShift = 0;
    if (!supportedSipTarget(
            atom, massNumber, isotopeIndex, nominalShift))
        return false;

    const int atomIx = atomIndex(atom);
    if (atomIx < 0 ||
        atomIx >= static_cast<int>(
            configIsotopologue.vAtomIsotopicDistribution.size()))
        return false;

    if (naturalAtomIsotopeProbabilities.size() ==
        configIsotopologue.vAtomIsotopicDistribution.size())
    {
        for (size_t index = 0;
             index < naturalAtomIsotopeProbabilities.size();
             ++index)
        {
            configIsotopologue.vAtomIsotopicDistribution[index].vProb =
                naturalAtomIsotopeProbabilities[index];
        }
    }

    getSetSIPelement() = std::string(1, atom);
    const bool updated = averagine::changeAtomProbability(
        configIsotopologue.vAtomIsotopicDistribution[
            static_cast<size_t>(atomIx)].vProb,
        atom,
        fraction);
    if (!updated ||
        !refreshResidueDistributions(configIsotopologue))
        return false;

    setDeductionCoefficient(false);
    if (!calculatePeptideMassWindowOffset())
        return false;
    getSetSIPelement() = std::string(1, atom);
    return true;
}

AminoAcidMasses::AminoAcidMasses()
{
	vdMasses.clear();
	vdMasses.resize(AminoAcidMassesSize, 0);
	for (int i = 0; i < AminoAcidMassesSize; i++)
	{
		vdMasses.at(i) = dNULL;
	}
}

void AminoAcidMasses::clear()
{
	for (int i = 0; i < AminoAcidMassesSize; i++)
	{
		vdMasses.at(i) = dNULL;
	}
}

double AminoAcidMasses::end()
{
	return dNULL;
}

double AminoAcidMasses::find(char _cAminoAcid)
{
	/*
	if (((int)_cAminoAcid) >= AminoAcidMassesSize ) {
		cerr << "error AminoAcidMasses. " << endl;
		exit(1);
		return dERROR;
	}*/
	return vdMasses.at((int)_cAminoAcid);
}

double AminoAcidMasses::operator[](char _cAminoAcid) const
{
	/*
	if (((int)_cAminoAcid) >= AminoAcidMassesSize) {
		cerr << "error AminoAcidMasses. " << endl;
		exit(1);
		return dERROR;
	}*/
	return vdMasses.at((int)_cAminoAcid);
}

double &AminoAcidMasses::operator[](char _cAminoAcid)
{
	/*
	if (_cAminoAcid >= AminoAcidMassesSize || _cAminoAcid < 0) {
		cerr << "error AminoAcidMasses. " << endl;
		exit(1);
		return dERROR;
	}*/
	return vdMasses.at((int)_cAminoAcid);
}

ProNovoConfig::ProNovoConfig()
{
}

bool ProNovoConfig::setFilename(const string &sConfigFileName)
{
	if (ProNovoConfigSingleton == NULL)
	{
		ProNovoConfigSingleton = new ProNovoConfig;
	}

	sFilename = sConfigFileName;

	// Try loading the file.
	if (!ProNovoConfigSingleton->parseConfigKeyValues())
	{
		cerr << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	if (!ProNovoConfigSingleton->getParameters())
	{
		return false;
	}

	naturalAtomIsotopeProbabilities.clear();
	for (const IsotopeDistribution &distribution :
		 configIsotopologue.vAtomIsotopicDistribution)
	{
		naturalAtomIsotopeProbabilities.push_back(
			distribution.vProb);
	}

	// parse neutral loss
	ProNovoConfigSingleton->NeutralLoss();

	// Compute the target-specific score coefficient and validate the
	// configured element/isotope pair before any workflow starts.
	try
	{
		ProNovoConfigSingleton->setDeductionCoefficient();
	}
	catch (const std::exception &error)
	{
		cerr << "ERROR: " << error.what() << endl;
		return false;
	}
	if (!ProNovoConfigSingleton->calculatePeptideMassWindowOffset())
	{
		return false;
	}

	// If everything goes fine return 0.
	return true;
}

bool ProNovoConfig::setWorkingDirectory(const string &sDirectoryName)
{
	if (sDirectoryName[sDirectoryName.size() - 1] == ProNovoConfig::getSeparator())
	{
		sWorkingDirectory = sDirectoryName;
	}
	else
	{
		sWorkingDirectory = sDirectoryName + ProNovoConfig::getSeparator();
	}

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

bool ProNovoConfig::getAtomIsotopicComposition(char cAtom, vector<double> &vdAtomicMass, vector<double> &vdComposition)
{

	// clear the input vectors
	vdAtomicMass.clear();
	vdComposition.clear();

	string sData;
	istringstream issStream;
	double dValue;
	string sAtom = "X";
	sAtom[0] = cAtom;

	map<string, string> mapElementMasses;
	if (!getConfigMasterKeyValue("[Peptide_Identification]Element_Masses", mapElementMasses))
	{
		cerr << "Error: cannot retrieve Element Masses." << endl;
		return false;
	}

	map<string, string> mapElementPercent;
	if (!getConfigMasterKeyValue("[Peptide_Identification]Element_Percent", mapElementPercent))
	{
		cerr << "Error: cannot retrieve Element Percent." << endl;
		return false;
	}

	map<string, string>::iterator iterMass = mapElementMasses.find(sAtom);
	if (iterMass == mapElementMasses.end())
	{
		cerr << "Error: cannot find element masses for element " << sAtom << endl;
		return false;
	}
	sData = iterMass->second;
	replaceDelimitor(sData, ',', '\t');
	// clear end of file state
	issStream.clear();
	// re-set the string associated with issStream
	issStream.str(sData);
	while (!(issStream.eof()))
	{
		issStream >> dValue;
		vdAtomicMass.push_back(dValue);
	}

	map<string, string>::iterator iterPercent = mapElementPercent.find(sAtom);
	if (iterPercent == mapElementPercent.end())
	{
		cerr << "Error: cannot find element percent for element " << sAtom << endl;
		return false;
	}
	sData = iterPercent->second;
	replaceDelimitor(sData, ',', '\t');
	// clear end of file state
	issStream.clear();
	// re-set the string associated with issStream
	issStream.str(sData);
	while (!(issStream.eof()))
	{
		issStream >> dValue;
		vdComposition.push_back(dValue);
	}

	return true;
}

bool ProNovoConfig::getResidueElementalComposition(string &sAtomicCompositionTable)
{
	sAtomicCompositionTable = "";

	map<string, string> mapResidueTable;
	if (!getConfigMasterKeyValue("[Peptide_Identification]Residue", mapResidueTable))
	{
		cerr << "Error: cannot retrieve Elemental composition of amino acid residues." << endl;
		return false;
	}

	map<string, string>::iterator iter;

	for (iter = mapResidueTable.begin(); iter != mapResidueTable.end(); ++iter)
	{
		sAtomicCompositionTable.append(iter->first);
		sAtomicCompositionTable.append(",\t");
		sAtomicCompositionTable.append(iter->second);
		sAtomicCompositionTable.append("\n");
	}
	replaceDelimitor(sAtomicCompositionTable, ',', '\t');
	return true;
}

bool ProNovoConfig::getPTMinfo(map<string, string> &mPTMinfo)
{
	mPTMinfo.clear();
	if (!getConfigMasterKeyValue("[Peptide_Identification]PTM", mPTMinfo))
	{
		cerr << "Error: cannot retrieve PTM information." << endl;
		return false;
	}
	return true;
}

bool ProNovoConfig::getParameters()
{

	string sTemp;
	istringstream issStream;

	// Extract the elements inside <Peptide_Identification>
	getConfigValue("[Peptide_Identification]Search_Type", sSearchType);
	getConfigValue("[Peptide_Identification]Search_Name", sSearchName);

	getConfigValue("[Peptide_Identification]FASTA_Database", sFASTAFilename);
	getConfigValue("[Peptide_Identification]Fragmentation_Method", sFragmentationMethod);

	if (sSearchType == "Regular")
	{
		getConfigValue("[Peptide_Identification]Max_PTM_Count", sTemp);
		issStream.clear();
		issStream.str(sTemp);
		issStream >> iMaxPTMcount;
	}

	getConfigValue("[Peptide_Identification]Search_Mass_Tolerance_Parent_Ion", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	issStream >> dMassAccuracyParentIon;

	getConfigValue("[Peptide_Identification]Mass_Tolerance_Fragment_Ions", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	issStream >> dMassAccuracyFragmentIon;

	getConfigValue("[Peptide_Identification]Parent_Mass_Windows", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	string sField;
	viParentMassWindows.clear();
	while (getline(issStream, sField, ','))
	{
		istringstream issField(sField);
		int iWindow;
		issField >> iWindow;
		viParentMassWindows.push_back(iWindow);
	}

	// read Peptide_Length
	getConfigValue("[Peptide_Identification]Minimum_Peptide_Length", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	issStream >> iMinPeptideLength;

	getConfigValue("[Peptide_Identification]Maximum_Peptide_Length", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	issStream >> iMaxPeptideLength;

	// read Cleavage_Rules
	getConfigValue("[Peptide_Identification]Cleave_After_Residues", sCleavageAfterResidues);
	getConfigValue("[Peptide_Identification]Cleave_Before_Residues", sCleavageBeforeResidues);

	getConfigValue("[Peptide_Identification]Maximum_Missed_Cleavages", sTemp);
	issStream.clear();
	issStream.str(sTemp);
	issStream >> iMaxMissedCleavages;

	getConfigValue("[Peptide_Identification]Try_First_Methionine", sTemp);
	if (sTemp == "TRUE" || sTemp == "True" || sTemp == "true" || sTemp == "T")
		bTestStartRemoval = true;
	else
		bTestStartRemoval = false;

	sElementList = "";
	getConfigValue("[Peptide_Identification]Element_List", sTemp);
	replaceDelimitor(sTemp, ',', '\t');
	issStream.clear();
	issStream.str(sTemp);
	while (!(issStream.eof()))
	{
		string sAtom;
		issStream >> sAtom;
		if (sAtom.size() == 1)
			sElementList.append(sAtom);
		else
			cerr << "Warning: Ignore an invalid element at Element_List " << sAtom << endl;
	}

	// populate vpPeptideMassWindowOffset
	calculatePeptideMassWindowOffset();

	// setup configIsotopologue, this must be done after the other parameters have been configured.
	string sResidueElementalComposition;
	getResidueElementalComposition(sResidueElementalComposition);
	configIsotopologue.setupIsotopologue(sResidueElementalComposition, sElementList);
	configIsotopologue.getSingleResidueMostAbundantMasses(vsSingleResidueNames, vdSingleResidueMasses, dTerminusMassN, dTerminusMassC);

	return true;
}

void ProNovoConfig::NeutralLoss()
{
	map<string, string> mPTMinfo;
	map<string, string>::iterator iter;
	pair<string, string> pCurrentPair;
	string sCurrentWholePTM, sOriginalPTM, sChangedPTM;
	vpNeutralLossList.clear();
	getPTMinfo(mPTMinfo);
	for (iter = mPTMinfo.begin(); iter != mPTMinfo.end(); iter++)
	{
		sCurrentWholePTM = iter->first;
		if (sCurrentWholePTM.length() > 1)
		{
			if (sCurrentWholePTM.substr(1, 2) == "to")
			{
				// consider neutral loss
				sOriginalPTM = sCurrentWholePTM.substr(0, 1);
				// if it is like PTM{@to}, sChangedPTM is ""
				sChangedPTM = (sCurrentWholePTM.length() == 3) ? "" : sCurrentWholePTM.substr(3, 1);
				pCurrentPair = make_pair(sOriginalPTM, sChangedPTM);
				vpNeutralLossList.push_back(pCurrentPair);
			}
			else
			{
				cerr << "illeagal ptm: " << sCurrentWholePTM << endl;
				exit(0);
			}
		}
	}
}

double ProNovoConfig::getResidueMass(string sResidue)
{
	unsigned int i;
	double dResidueMass = 0.0;
	if (sResidue == "|||")
	{
		dResidueMass = 0.0;
		return dResidueMass;
	}
	for (i = 0; i < vsSingleResidueNames.size(); ++i)
	{
		if (vsSingleResidueNames[i] == sResidue)
		{
			dResidueMass = vdSingleResidueMasses[i];
			return dResidueMass;
		}
	}

	cerr << "ERROR: cannot find residue " << sResidue << endl;
	return dResidueMass;
}

void ProNovoConfig::replaceDelimitor(string &sLine, char cOldDelimitor, char cNewDelimitor)
{
	int iLength = sLine.length();
	for (int i = 0; i < iLength; ++i)
	{
		if (sLine[i] == cOldDelimitor)
			sLine[i] = cNewDelimitor;
	}
	return;
}

// parse the cfg file to populate mapConfigKeyValues
bool ProNovoConfig::parseConfigKeyValues()
{
	bool bReVal = true;
	string sline, sWhiteSpaces(" \t\f\v\n\r");
	size_t poundPos, whitespacePos;
	//    map<string,string>::iterator it;

	ifstream config_stream(sFilename.c_str());
	bReVal = config_stream.is_open();
	mapConfigKeyValues.clear();
	if (bReVal)
	{
		while (!config_stream.eof())
		{
			sline.clear();
			getline(config_stream, sline);
			poundPos = sline.find("#");
			if (poundPos != string::npos)
				sline.erase(poundPos);
			whitespacePos = sline.find_last_not_of(sWhiteSpaces);
			if (whitespacePos != string::npos)
				sline.erase(whitespacePos + 1);
			else
				// if no character is non-whitespace, make string clear
				// the previous version is sline.erase(0), but it can't be accepted by pgCC 11.9-0
				// sline.erase(0);
				sline.clear();
			whitespacePos = sline.find_first_not_of(sWhiteSpaces);
			if (whitespacePos != string::npos)
				if (whitespacePos != 0)
					sline.erase(0, whitespacePos);
			if (sline.length() > 0)
				parseConfigLine(sline);
		}

		//	for ( it=mapConfigKeyValues.begin() ; it != mapConfigKeyValues.end(); it++ )
		//	    cout << (*it).first << " => " << (*it).second << endl;
	}
	else
		cerr << "Can't open configure file " << sFilename << endl;
	config_stream.clear();
	config_stream.close();
	return bReVal;
}

// get the value of a key;
bool ProNovoConfig::getConfigValue(string sConfigKey, string &sConfigValue)
{
	sConfigValue = "";

	map<string, string>::iterator iter = mapConfigKeyValues.find(sConfigKey);
	if (iter != mapConfigKeyValues.end())
	{
		sConfigValue = iter->second;
		return true;
	}
	else
	{
		sConfigValue = "";
		cerr << "Warning: Cannot find parameter " << sConfigKey << " in the Config file." << endl;
		return false;
	}
}

// get a set of key-value pairs, given a master key
bool ProNovoConfig::getConfigMasterKeyValue(string sMasterKey, map<string, string> &mapKeyValueSet)
{
	bool bReVal = true;
	map<string, string>::iterator iter;
	size_t iKeyLength;
	string sCurrentCoreKey;
	string sCurrentKey;

	mapKeyValueSet.clear();
	iKeyLength = sMasterKey.length();
	for (iter = mapConfigKeyValues.begin(); iter != mapConfigKeyValues.end(); iter++)
	{
		// cout << (*it).first << " => " << (*it).second << endl;
		sCurrentKey = (*iter).first;
		if ((sCurrentKey.substr(0, iKeyLength + 1) == (sMasterKey + "{")) && (sCurrentKey.at(sCurrentKey.length() - 1) == '}') && (sCurrentKey.length() > (iKeyLength + 2)))
		{
			sCurrentCoreKey = sCurrentKey.substr(iKeyLength + 1, sCurrentKey.length() - iKeyLength - 2);
			mapKeyValueSet.insert(pair<string, string>(sCurrentCoreKey, (*iter).second));
		}
	}

	return bReVal;
}

// parse one line in Configfile
bool ProNovoConfig::parseConfigLine(const std::string &sLine)
{
	bool bReVal = true;
	size_t equalPos, leftendPos, rightBeginPos; // position of "=", last nonwhitespace before "=", first nonwhitespace after "="
	string sKey, sValue;
	pair<map<string, string>::iterator, bool> ret; // if ret.second == false, key is not unique
												   //    cout<<"beg!"<<sLine<<"!end"<<endl;
	if ((sLine.at(0) == '[') && (sLine.at(sLine.length() - 1) == ']'))
		sSectionName = sLine;
	else
	{
		if (sSectionName == "")
		{
			cerr << "can't find the section name" << endl;
			bReVal = false;
		}
		else
		{
			equalPos = sLine.find("=");
			if (equalPos == string::npos)
			{
				cerr << "can't find = " << endl;
				bReVal = false;
			}
			else
			{
				if ((equalPos == 0) || (equalPos == (sLine.length() - 1)))
				{
					cerr << "can't find key or value" << endl;
					bReVal = false;
				}
				else
				{
					leftendPos = sLine.find_last_not_of(" \t\f\v\n\r", equalPos - 1);
					rightBeginPos = sLine.find_first_not_of(" \t\f\v\n\r", equalPos + 1);
					sKey = sLine.substr(0, leftendPos + 1);
					sValue = sLine.substr(rightBeginPos);
					// cout<<"beg!"<<sSectionName+sKey<<"!"<<sValue<<"!end"<<endl;
					ret = mapConfigKeyValues.insert(pair<string, string>(sSectionName + sKey, sValue));
					if (ret.second == false)
					{
						cerr << "Key " << sSectionName + sKey << " has existed with value of " << ret.first->second << endl;
						bReVal = false;
					}
				}
			}
		}
	}
	return bReVal;
}

bool ProNovoConfig::calculatePeptideMassWindowOffset()
{
	bool bReVal = true;
	int i;
	double dLastUpperBound = -1000, dLastLowerBound = -1000; // last range of acceptable parent mass
	double dCurrentUpperBound, dCurrentLowerBound;

	vpPeptideMassWindowOffset.clear();
	sort(viParentMassWindows.begin(), viParentMassWindows.end());
	for (i = 0; i < (int)viParentMassWindows.size(); i++)
	{
		dCurrentLowerBound = viParentMassWindows.at(i) * getNeutronMass() - dMassAccuracyParentIon;
		dCurrentUpperBound = viParentMassWindows.at(i) * getNeutronMass() + dMassAccuracyParentIon;
		if (dLastUpperBound < -100)
		{
			dLastLowerBound = dCurrentLowerBound;
			dLastUpperBound = dCurrentUpperBound;
		}
		else
		{
			if (dCurrentLowerBound <= dLastUpperBound)
				dLastUpperBound = dCurrentUpperBound;
			else
			{
				vpPeptideMassWindowOffset.push_back(pair<double, double>(dLastLowerBound, dLastUpperBound));
				dLastLowerBound = dCurrentLowerBound;
				dLastUpperBound = dCurrentUpperBound;
			}
		}
	}
	if (dLastUpperBound > -100)
		vpPeptideMassWindowOffset.push_back(pair<double, double>(dLastLowerBound, dLastUpperBound));

	return bReVal;
}

bool ProNovoConfig::getPeptideMassWindows(double dPeptideMass, vector<pair<double, double>> &vpPeptideMassWindows)
{
	bool bReVal = true;
	double dCurrentLowerBound, dCurrentUpperBound;
	int i;
	for (i = 0; i < (int)vpPeptideMassWindowOffset.size(); i++)
	{
		dCurrentLowerBound = dPeptideMass + vpPeptideMassWindowOffset.at(i).first;
		dCurrentUpperBound = dPeptideMass + vpPeptideMassWindowOffset.at(i).second;
		vpPeptideMassWindows.push_back(pair<double, double>(dCurrentLowerBound, dCurrentUpperBound));
	}
	return bReVal;
}

// compute deduction coefficient in score function
// only suitbale for carbon and nitrogen SIP now
void ProNovoConfig::setDeductionCoefficient(bool readConfigElement)
{
    if (readConfigElement)
    {
        std::string configuredAtom;
        if (!getConfigValue(
                "[Stable_Isotope_Probing]SIP_Element",
                configuredAtom))
        {
            throw std::runtime_error(
                "Missing SIP_Element in configuration.");
        }
        getSetSIPelement() = configuredAtom;
    }

    const std::string selected = getSetSIPelement();
    if (selected.size() != 1)
    {
        throw std::runtime_error(
            "Unsupported SIP isotope in configuration. "
            "Use C13,H2,N15,O18,S34.");
    }
    const char atom = static_cast<char>(std::toupper(
        static_cast<unsigned char>(selected[0])));
    int massNumber = 0;
    int isotopeIndex = -1;
    int nominalShift = 0;
    if (!supportedSipTarget(
            atom, massNumber, isotopeIndex, nominalShift))
    {
        throw std::runtime_error(
            "Unsupported SIP isotope in configuration. "
            "Use C13,H2,N15,O18,S34.");
    }

    if (readConfigElement)
    {
        std::string configuredMassText;
        if (!getConfigValue(
                "[Stable_Isotope_Probing]SIP_Element_Isotope",
                configuredMassText))
        {
            throw std::runtime_error(
                "Missing SIP_Element_Isotope in configuration.");
        }
        size_t consumed = 0;
        int configuredMassNumber = 0;
        try
        {
            configuredMassNumber =
                std::stoi(configuredMassText, &consumed);
        }
        catch (const std::exception &)
        {
            throw std::runtime_error(
                "Invalid SIP_Element_Isotope in configuration: " +
                configuredMassText);
        }
        if (consumed != configuredMassText.size() ||
            configuredMassNumber != massNumber)
        {
            throw std::runtime_error(
                "Unsupported SIP isotope in configuration: " +
                std::string(1, atom) + configuredMassText +
                ". Use C13,H2,N15,O18,S34.");
        }
        SIPisotopeMassNumber = configuredMassNumber;
    }
    else
    {
        SIPisotopeMassNumber = massNumber;
    }
    getSetSIPelement() = std::string(1, atom);

    const int atomIx = atomIndex(atom);
    const auto &distribution =
        configIsotopologue.vAtomIsotopicDistribution[
            static_cast<size_t>(atomIx)];
    neutronMass =
        (distribution.vMass[static_cast<size_t>(isotopeIndex)] -
         distribution.vMass[0]) /
        nominalShift;

    string minValueStr;
    string foldStr;
    getConfigValue(
        "[Stable_Isotope_Probing]minValue", minValueStr);
    getSetMinValue() = stod(minValueStr);
    getConfigValue(
        "[Stable_Isotope_Probing]fold", foldStr);
    getSetFold() = stod(foldStr);
    deductionCoefficient =
        -(getSetMinValue() +
          getSetFold() *
              std::pow(
                  distribution.vProb[
                      static_cast<size_t>(isotopeIndex)] -
                      0.5,
                  8));
}
