#ifndef SIPROS_SEARCH_RUNNER_H
#define SIPROS_SEARCH_RUNNER_H

#include <iosfwd>
#include <string>
#include <vector>

namespace sipros
{

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
	std::string configFile;
	std::string workingDirectory;
	std::string singleWorkingFile;
	std::string fastaFile;
	std::string outputDirectory;
	std::string sipElementSpec;
	std::string sipRangeSpec;
	double sipStepPct = 1.0;
	int pinLabel = 1;
	int topPsmsPerScan = 0;
	bool screenOutput = true;
	bool showHelp = false;
};

class SiprosSearchRunner
{
public:
	static void printUsage(std::ostream &out);

	bool initializeArguments(int argc, char **argv,
						 DatabaseSearchArguments &args,
						 std::ostream &out,
						 std::ostream &err) const;

	int runScan(const std::string &scanFile,
				const DatabaseSearchArguments &args) const;

private:
	static std::vector<std::string> listFilesWithExtensions(const std::string &directory,
										 const std::vector<std::string> &extensions);
};

} // namespace sipros

#endif // SIPROS_SEARCH_RUNNER_H
