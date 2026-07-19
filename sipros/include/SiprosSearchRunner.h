#ifndef SIPROS_SEARCH_RUNNER_H
#define SIPROS_SEARCH_RUNNER_H

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "RaxportHdf5Reader.h"

namespace sipros
{

class FragmentIndex;

class TextUtils
{
public:
	static std::string toLower(std::string value);
	static std::string trim(const std::string &value);
	static std::vector<std::string> splitTab(const std::string &line);
	static bool parseSipAtomSpec(const std::string &spec, char &sipAtom, int &sipIsotopeMassNumber);
};

struct DatabaseSearchArguments
{
	std::vector<std::string> scanFiles;
	std::string workingDirectory;
	std::string fastaFile;
	std::string outputDirectory;
	std::string pinOutputFile;
	std::string fragmentIndexCache;
	std::string sipElementSpec;
	std::string sipRangeSpec;
	std::vector<std::string> fixedPtmSelectors;
	std::vector<std::string> ptmSelectors;
	RaxportReadOptions raxportReadOptions;
	double sipStepPct = 1.0;
	double toleranceMs1Da = 0.0;
	double toleranceMs2Da = 0.0;
	int pinLabel = 1;
	int topPsmsPerScan = 0;
	int maxPtmCountOverride = -1;
	bool sipStepProvided = false;
	bool toleranceMs1Provided = false;
	bool toleranceMs2Provided = false;
	bool listPtms = false;
	bool rebuildFragmentIndex = false;
	bool prepareOnly = false;
	bool precursorSourceProvided = false;
	bool showHelp = false;
};

class SiprosSearchRunner
{
public:
	static void printUsage(std::ostream &out, const std::string &prog = "sipros search-fasta");

	bool initializeArguments(int argc, char **argv,
						 DatabaseSearchArguments &args,
						 std::ostream &out,
						 std::ostream &err) const;

	// Configure the profile once and build/map the Regular peptide index once.
	// The same prepared runner can then search every scan file in one invocation.
	int prepare(const DatabaseSearchArguments &args);

	int runScan(const std::string &scanFile,
				const DatabaseSearchArguments &args);

private:
	static std::vector<std::string> listFilesWithExtensions(const std::string &directory,
										 const std::vector<std::string> &extensions);

	std::shared_ptr<FragmentIndex> fragmentIndex_;
	bool prepared_ = false;
	bool preparedSipMode_ = false;
};

} // namespace sipros

#endif // SIPROS_SEARCH_RUNNER_H
