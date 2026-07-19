#include "SiprosWorkflows.h"

#include <iostream>
#include <string>

#include <omp.h>

#include "SiprosSearchRunner.h"
#include "performancelog.h"

int DatabaseSearchWorkflow::run(int argc, char **argv)
{
	sipros::SiprosSearchRunner runner;
	sipros::DatabaseSearchArguments args;
	if (!runner.initializeArguments(argc, argv, args, std::cout, std::cerr))
	{
		return 1;
	}
	if (args.showHelp)
	{
		return 0;
	}
	const int threadCount = omp_get_max_threads();
	std::cout << "Sipros FASTA search\n"
			  << "  Input files: " << args.scanFiles.size() << "\n"
			  << "  Mode       : "
			  << (args.prepareOnly ? "prepare peptide cache" : "search H5 scans") << "\n"
			  << "  Threads    : " << threadCount << "\n";
	const sipros::PerformanceTimer runTimer;
	if (runner.prepare(args) != 0)
	{
		return 1;
	}
	if (args.prepareOnly)
	{
		sipros::printPerformanceHeader(std::cout, "Run summary", threadCount);
		sipros::printPerformanceStage(
			std::cout, "Prepare peptide cache", runTimer.elapsed());
		sipros::printPerformanceFooter(std::cout);
		std::cout << "\nSipros peptide-cache preparation complete.\n";
		return 0;
	}

	int result = 0;
	for (const std::string &scanFile : args.scanFiles)
	{
		const int scanResult = runner.runScan(scanFile, args);
		if (scanResult != 0)
		{
			result = scanResult;
		}
	}

	sipros::printPerformanceHeader(std::cout, "Run summary", threadCount);
	sipros::printPerformanceStage(
		std::cout, "Complete FASTA search", runTimer.elapsed(),
		sipros::formatPerformanceCount(args.scanFiles.size()) +
			" input file(s)");
	sipros::printPerformanceFooter(std::cout);
	std::cout << "\nSipros search complete.\n";
	return result;
}
