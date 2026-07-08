#include "SiprosWorkflows.h"

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace
{

void printUsage(const char *prog)
{
	std::cerr << "Usage:\n"
			  << "  " << prog << " [database-search options]\n"
			  << "  " << prog << " theoretical-spectra [options]\n"
			  << "  " << prog << " experimental-spectra [options]\n"
			  << "  " << prog << " search-spectra [options]\n\n"
			  << "Run a subcommand with --help for command-specific options.\n";
}

template <typename Workflow>
int runShifted(int argc, char **argv, int skip, const std::string &displayName, Workflow &workflow)
{
	std::vector<char *> shifted;
	std::string programName = displayName;
	shifted.push_back(programName.data());
	for (int i = skip; i < argc; ++i)
	{
		shifted.push_back(argv[i]);
	}
	return workflow.run(static_cast<int>(shifted.size()), shifted.data());
}

bool isDatabaseOption(const char *value)
{
	return value != nullptr && value[0] == '-';
}

} // namespace

int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		DatabaseSearchWorkflow workflow;
		return workflow.run(argc, argv);
	}

	const std::string command = argv[1];
	if (command == "-h" || command == "--help")
	{
		printUsage(argv[0]);
		return 0;
	}
	if (command == "theoretical-spectra")
	{
		TheoreticalSpectraWorkflow workflow;
		return runShifted(argc, argv, 2, std::string(argv[0]) + " theoretical-spectra", workflow);
	}
	if (command == "experimental-spectra")
	{
		ExperimentalSpectraWorkflow workflow;
		return runShifted(argc, argv, 2, std::string(argv[0]) + " experimental-spectra", workflow);
	}
	if (command == "search-spectra")
	{
		SearchSpectraWorkflow workflow;
		return runShifted(argc, argv, 2, std::string(argv[0]) + " search-spectra", workflow);
	}
	if (isDatabaseOption(argv[1]))
	{
		DatabaseSearchWorkflow workflow;
		return workflow.run(argc, argv);
	}

	std::cerr << "Unknown sipros subcommand: " << command << "\n\n";
	printUsage(argv[0]);
	return 1;
}
