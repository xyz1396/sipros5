#include "SiprosWorkflows.h"

#include <iostream>
#include <string>

#include <omp.h>

#include "SiprosSearchRunner.h"

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

#pragma omp parallel
	{
		const int tid = omp_get_thread_num();
		if (tid == 0)
		{
			std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
		}
	}

	const double begin = omp_get_wtime();

	int result = 0;
	for (const std::string &scanFile : args.scanFiles)
	{
		const int scanResult = runner.runScan(scanFile, args);
		if (scanResult != 0)
		{
			result = scanResult;
		}
	}

	const double end = omp_get_wtime();
	std::cout << "Sipros finished in " << double(end - begin) << " Seconds." << std::endl
			  << std::endl;
	return result;
}
