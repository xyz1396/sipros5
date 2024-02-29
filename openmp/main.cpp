#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include "proNovoConfig.h"
#include "ms2scanvector.h"

using namespace std;
namespace fs = std::filesystem;

static void searchFT2Files(std::vector<std::string> &vsFT2Filenames, const std::string &sWorkingDirectory) {
    const std::vector<std::string> patterns = {".ft2", ".FT2", ".ms2", ".MS2", ".mzml", ".mzML"};    
    for (const auto &entry : fs::directory_iterator(sWorkingDirectory)) {
        if (entry.is_regular_file()) {
            for (const auto &pattern : patterns) {
                if (entry.path().extension() == pattern) {
                    vsFT2Filenames.push_back(entry.path().string());
                    break;
                }
            }
        }
    }
    if (vsFT2Filenames.empty()) {
        std::cerr << "no scan file in the working directory" << std::endl;
        exit(1);
    }
    char separator = ProNovoConfig::getSeparator();
    for (auto &filename : vsFT2Filenames) {
        filename = sWorkingDirectory + separator + filename;
    }
}

static void searchConfigureFiles(vector<string> &vsConfigureFilenames, const string &sConfigFileDirectory)
{
	const std::vector<std::string> patterns = { ".CFG", ".cfg"};
	for (const auto& entry : fs::directory_iterator(sConfigFileDirectory)) {
		if (entry.is_regular_file()) {
			for (const auto& pattern : patterns) {
				if (entry.path().extension() == pattern) {
					vsConfigureFilenames.push_back(entry.path().string());
					break;
				}
			}
		}
	}
	if (vsConfigureFilenames.empty()) {
		std::cerr << "no scan file in the working directory" << std::endl;
		exit(1);
	}
	char separator = ProNovoConfig::getSeparator();
	for (auto& filename : vsConfigureFilenames) {
		filename = sConfigFileDirectory + separator + filename;
	}
}

/*
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

static void initializeArguments(int argc, char **argv, vector<string> &vsFT2Filenames, string &sWorkingDirectory, vector<string> &vsConfigureFilenames,
						 string &sSingleWorkingFile, string &sFastaFile, string &sOutputDirectory, bool &bScreenOutput)
{
	int i;
	string sConfigFileDirectory, sConfigFilename;

	// Grab command line arguments
	vector<string> vsArguments;

	sWorkingDirectory = "";
	sConfigFilename = "";
	sSingleWorkingFile = "";
	sFastaFile = "";
	sOutputDirectory = "";
	sConfigFileDirectory = "";
	bScreenOutput = true;

	while (argc--)
		vsArguments.push_back(*argv++);
	for (i = 1; i <= (int)vsArguments.size() - 1; i++)
		if (vsArguments[i] == "-w")
			sWorkingDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-c")
			sConfigFilename = vsArguments[++i];
		else if (vsArguments[i] == "-f")
			sSingleWorkingFile = vsArguments[++i];
		else if (vsArguments[i] == "-fasta")
			sFastaFile = vsArguments[++i];
		else if (vsArguments[i] == "-o")
			sOutputDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-g")
			sConfigFileDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-s")
			bScreenOutput = false;
		else if ((vsArguments[i] == "-h") || (vsArguments[i] == "--help"))
		{
			cout << "Usage 1: " << endl;
			cout << "-w WorkingDirectory -c ConfigurationFile -o output directory" << endl;
			cout << "Usage 2: " << endl;
			cout << "-f single_ms2_file -c ConfigurationFile -o output directory" << endl;
			cout << "Usage 3: " << endl;
			// set FASTA file out of cfg file
			cout << "-f single_ms2_file -fasta proteins.fasta -c ConfigurationFile -o output directory" << endl;
			cout << "Other options: " << endl;
			cout << "-s silence all standard output." << endl;
			exit(0);
		}
		else if (vsArguments[i] == "-p")
		{
			MVH::ProbabilityCutOff = atof(vsArguments[++i].c_str());
			CometSearchMod::ProbabilityCutOff = MVH::ProbabilityCutOff;
		}
		else
		{
			cerr << "Unknown option " << vsArguments[i] << endl
				 << endl;
			exit(1);
		}
	if ((sWorkingDirectory == "") && (sSingleWorkingFile == ""))
		sWorkingDirectory = ".";

	if ((sWorkingDirectory != "") && (sSingleWorkingFile != ""))
	{
		cerr << "Either a input scan file or the directory of input scan files needs to be specified" << endl;
		exit(1);
	}

	if ((sConfigFilename == "") && (sConfigFileDirectory == ""))
		// sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "SiprosConfig.cfg";
		//  Without specifying configure file and configure file directory, the default configure file directory
		sConfigFileDirectory = sWorkingDirectory;

	if (sConfigFileDirectory != "")
		searchConfigureFiles(vsConfigureFilenames, sConfigFileDirectory);
	else
		vsConfigureFilenames.push_back(sConfigFilename);

	if (sSingleWorkingFile != "")
		vsFT2Filenames.push_back(sSingleWorkingFile);
	else
		searchFT2Files(vsFT2Filenames, sWorkingDirectory);
	if ((sOutputDirectory == "") && (sWorkingDirectory != ""))
		sOutputDirectory = sWorkingDirectory;
}

static void handleScan(const string &sFT2filename, const string &sOutputDirectory, const string &sConfigFilename, bool bScreenOutput)
{
	MS2ScanVector *pMainMS2ScanVector = new MS2ScanVector(sFT2filename, sOutputDirectory, sConfigFilename, bScreenOutput);

	if (bScreenOutput)
	{
		cout << "Reading MS2 scan file " << sFT2filename << endl;
		cout << "Using Configuration file " << sConfigFilename << endl;
	}

	if (!pMainMS2ScanVector->loadMassData())
		cerr << "Error: Failed to load file: " << sFT2filename << endl;
	else
	{
		if (ProNovoConfig::getSearchType() == "SIP")
		{
			pMainMS2ScanVector->startProcessingWdpSip();
		}
		else
		{
			// search all MS2 scans and write output to a file
			pMainMS2ScanVector->startProcessingMvh();
		}
	}
	delete pMainMS2ScanVector; // free memory of vpAllMS2Scans
}

int main(int argc, char **argv)
{
	// record the start time point
	double begin = omp_get_wtime();
	// A list of FT2/MS2 files to be searched
	vector<string> vsFT2Filenames;
	// A list of configure files
	vector<string> vsConfigureFilenames;
	bool bScreenOutput;
	string sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sFastaFile, sOutputDirectory;
	initializeArguments(argc, argv, vsFT2Filenames, sWorkingDirectory, vsConfigureFilenames, sSingleWorkingFile, sFastaFile, sOutputDirectory, bScreenOutput);
	// omp_set_num_threads(32);

	// Process one FT2 file at a time
	for (int i = 0; i < (int)vsFT2Filenames.size(); i++)
	{
		for (int j = 0; j < (int)vsConfigureFilenames.size(); ++j)
		{
			// Load config file
			sConfigFilename = vsConfigureFilenames.at(j);
			if (!ProNovoConfig::setFilename(sConfigFilename))
			{
				cerr << "Could not load config file " << sConfigFilename << endl;
				exit(1);
			}
			// set FASTA file out of cfg file
			if (sFastaFile != "")
				ProNovoConfig::setFASTAfilename(sFastaFile);
			handleScan(vsFT2Filenames.at(i), sOutputDirectory, sConfigFilename, bScreenOutput);
		}
	}
	// record the end time point
	double end = omp_get_wtime();
	cout << "\nSipros finished in :" << double(end - begin) << " Seconds." << endl
		 << endl;
	return 0;
}
