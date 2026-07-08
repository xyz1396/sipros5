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

namespace
{

struct NullBuffer : std::streambuf
{
	int overflow(int c) override { return c; }
};

struct UnitOfWork
{
	std::string scanFile;
};

std::vector<UnitOfWork> buildWorkload(const sipros::DatabaseSearchArguments &args)
{
	std::vector<UnitOfWork> workload;
	for (const std::string &scanFile : args.scanFiles)
	{
		workload.push_back({scanFile});
	}
	return workload;
}

void masterProcess(const std::vector<UnitOfWork> &workload, bool screenOutput)
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
	if (screenOutput)
	{
		std::cout << "Master process is done." << std::endl;
	}
}

void slaveProcess(const std::vector<UnitOfWork> &workload, const sipros::DatabaseSearchArguments &args, bool screenOutput)
{
	MPI_Status status;
	int currentWorkId = 0;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	while (true)
	{
		MPI_Recv(&currentWorkId, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == DIETAG)
		{
			break;
		}
		const UnitOfWork &work = workload.at(static_cast<size_t>(currentWorkId));
		std::cout << "slave Rank:\t" << rank << "\tscan:\t" << work.scanFile
				  << "\tCfg:\t" << args.configFile << std::endl;
		ProNovoConfig::iRank = rank;
		sipros::SiprosSearchRunner runner;
		const int result = runner.runScan(work.scanFile, args);
		if (screenOutput)
		{
			std::cout << work.scanFile << " and " << args.configFile
					  << " is done by Slave process " << rank << std::endl;
		}
		MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if (screenOutput)
	{
		std::cout << "Slave process " << rank << " is done." << std::endl;
	}
}

} // namespace

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	NullBuffer nullBuffer;
	std::ostream nullStream(&nullBuffer);

	sipros::SiprosSearchRunner runner;
	sipros::DatabaseSearchArguments args;
	const bool parsed = runner.initializeArguments(
		argc, argv, args, rank == 0 ? std::cout : nullStream, rank == 0 ? std::cerr : nullStream);
	if (!parsed || args.showHelp)
	{
		MPI_Finalize();
		return parsed ? 0 : 1;
	}

	const std::vector<UnitOfWork> workload = buildWorkload(args);
	if (rank == 0)
	{
		masterProcess(workload, args.screenOutput);
	}
	else
	{
		slaveProcess(workload, args, args.screenOutput);
	}
	MPI_Finalize();
	return 0;
}
