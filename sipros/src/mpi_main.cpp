#include <mpi.h>

#include <iostream>
#include <ostream>
#include <streambuf>
#include <string>
#include <vector>

#include "proNovoConfig.h"
#include "SiprosSearchRunner.h"

#define WORKTAG 1
#define DIETAG 2


struct NullBuffer : std::streambuf
{
	int overflow(int c) override { return c; }
};

struct UnitOfWork
{
	std::string scanFile;
};

void printUsage(const char *prog)
{
	std::cerr << "Usage:\n"
			  << "  " << prog << " search-fasta [options]\n\n"
			  << "Run the search-fasta subcommand with --help for command-specific options.\n";
}

std::vector<UnitOfWork> buildWorkload(const sipros::DatabaseSearchArguments &args)
{
	std::vector<UnitOfWork> workload;
	for (const std::string &scanFile : args.scanFiles)
	{
		workload.push_back({scanFile});
	}
	return workload;
}

void masterProcess(const std::vector<UnitOfWork> &workload)
{
	size_t i = 0;
	int currentWorkId = 0;
	int processorCount = 0;
	int result = 0;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &processorCount);
	const int slaveCount = processorCount - 1;
	const size_t initialWorkers = workload.size() <= static_cast<size_t>(slaveCount)
								? workload.size()
								: static_cast<size_t>(slaveCount);

	for (i = 1; i <= initialWorkers; ++i)
	{
		currentWorkId = static_cast<int>(i - 1);
		MPI_Send(&currentWorkId, 1, MPI_INT, static_cast<int>(i), WORKTAG, MPI_COMM_WORLD);
	}
	if (static_cast<int>(workload.size()) > slaveCount)
	{
		currentWorkId = slaveCount;
		while (currentWorkId < static_cast<int>(workload.size()))
		{
			MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send(&currentWorkId, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
			++currentWorkId;
		}
	}
	for (i = 1; i <= static_cast<size_t>(slaveCount); ++i)
	{
		MPI_Send(nullptr, 0, MPI_INT, static_cast<int>(i), DIETAG, MPI_COMM_WORLD);
	}
	std::cout << "Master process done" << std::endl;
}

void slaveProcess(const std::vector<UnitOfWork> &workload, const sipros::DatabaseSearchArguments &args)
{
	MPI_Status status;
	int currentWorkId = 0;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	sipros::SiprosSearchRunner runner;
	while (true)
	{
		MPI_Recv(&currentWorkId, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == DIETAG)
		{
			break;
		}
		const UnitOfWork &work = workload.at(static_cast<size_t>(currentWorkId));
		std::cout << "Slave process " << rank << " started " << work.scanFile
				  << std::endl;
		ProNovoConfig::iRank = rank;
		const int result = runner.runScan(work.scanFile, args);
		std::cout << "Slave process " << rank << " finished " << work.scanFile << std::endl;
		MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	std::cout << "Slave process " << rank << " done" << std::endl;
}


int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	NullBuffer nullBuffer;
	std::ostream nullStream(&nullBuffer);

	sipros::SiprosSearchRunner runner;
	sipros::DatabaseSearchArguments args;

	if (argc <= 1)
	{
		if (rank == 0)
		{
			printUsage(argv[0]);
		}
		MPI_Finalize();
		return 1;
	}
	const std::string command = argv[1];
	if (command == "-h" || command == "--help")
	{
		if (rank == 0)
		{
			printUsage(argv[0]);
		}
		MPI_Finalize();
		return 0;
	}
	if (command != "search-fasta")
	{
		if (rank == 0)
		{
			std::cerr << "Unknown siprosMPI subcommand: " << command << "\n\n";
			printUsage(argv[0]);
		}
		MPI_Finalize();
		return 1;
	}

	std::vector<char *> shifted;
	std::string programName = std::string(argv[0]) + " search-fasta";
	shifted.push_back(programName.data());
	for (int i = 2; i < argc; ++i)
	{
		shifted.push_back(argv[i]);
	}
	const bool parsed = runner.initializeArguments(
		static_cast<int>(shifted.size()), shifted.data(), args, rank == 0 ? std::cout : nullStream, rank == 0 ? std::cerr : nullStream);
	if (!parsed || args.showHelp)
	{
		MPI_Finalize();
		return parsed ? 0 : 1;
	}

	const std::vector<UnitOfWork> workload = buildWorkload(args);
	if (rank == 0)
	{
		masterProcess(workload);
	}
	else
	{
		slaveProcess(workload, args);
	}
	MPI_Finalize();
	return 0;
}
