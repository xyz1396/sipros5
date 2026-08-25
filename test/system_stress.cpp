// Standalone CPU and memory integrity stress test.
//
// This intentionally does not link any Sipros, HDF5, BLAS, or Percolator code.
// Build example:
//   x86_64-conda-linux-gnu-g++ -O2 -std=c++17 -fopenmp test/system_stress.cpp -o /tmp/sipros_system_stress
// Run example:
//   /tmp/sipros_system_stress --threads 16 --memory-gib 12 --seconds 300

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>


constexpr std::uint64_t kRoundSalt = UINT64_C(0x9e3779b97f4a7c15);
constexpr std::uint64_t kThreadSalt = UINT64_C(0xd1b54a32d192ed03);

struct Options
{
	int threads = std::max(1, std::min(8, omp_get_max_threads()));
	double memoryGiB = 1.0;
	int seconds = 60;
	bool allowHighMemory = false;
};

std::uint64_t mix64(std::uint64_t value) noexcept
{
	value ^= value >> 30;
	value *= UINT64_C(0xbf58476d1ce4e5b9);
	value ^= value >> 27;
	value *= UINT64_C(0x94d049bb133111eb);
	return value ^ (value >> 31);
}

std::uint64_t expectedValue(std::size_t index, std::uint64_t round, int threadId) noexcept
{
	return mix64(static_cast<std::uint64_t>(index) ^
		((round + 1) * kRoundSalt) ^
		((static_cast<std::uint64_t>(threadId) + 1) * kThreadSalt));
}

std::size_t parseAvailableMemoryBytes()
{
	std::ifstream input("/proc/meminfo");
	std::string key;
	std::uint64_t valueKiB = 0;
	std::string unit;
	while (input >> key >> valueKiB >> unit)
	{
		if (key == "MemAvailable:")
		{
			const std::uint64_t bytes = valueKiB * UINT64_C(1024);
			return bytes > std::numeric_limits<std::size_t>::max()
				? std::numeric_limits<std::size_t>::max()
				: static_cast<std::size_t>(bytes);
		}
	}
	return 0;
}

std::string formatGiB(std::size_t bytes)
{
	std::ostringstream output;
	output << std::fixed << std::setprecision(2)
		   << static_cast<long double>(bytes) /
				  static_cast<long double>(UINT64_C(1024) * 1024 * 1024);
	return output.str();
}

int parsePositiveInt(const std::string &text, const char *option)
{
	std::size_t consumed = 0;
	const long long value = std::stoll(text, &consumed);
	if (consumed != text.size() || value <= 0 || value > std::numeric_limits<int>::max())
	{
		throw std::invalid_argument(std::string(option) + " requires a positive integer");
	}
	return static_cast<int>(value);
}

double parsePositiveDouble(const std::string &text, const char *option)
{
	std::size_t consumed = 0;
	const double value = std::stod(text, &consumed);
	if (consumed != text.size() || !std::isfinite(value) || value <= 0.0)
	{
		throw std::invalid_argument(std::string(option) + " requires a positive number");
	}
	return value;
}

void printUsage(const char *program)
{
	std::cout
		<< "Usage: " << program << " [options]\n"
		<< "  --threads N            OpenMP worker threads (default: up to 8)\n"
		<< "  --memory-gib N         Total allocated memory in GiB (default: 1)\n"
		<< "  --seconds N            Minimum test duration in seconds (default: 60)\n"
		<< "  --allow-high-memory    Permit allocation above 70% of MemAvailable\n"
		<< "  -h, --help             Show this help\n";
}

Options parseOptions(int argc, char **argv)
{
	Options options;
	for (int i = 1; i < argc; ++i)
	{
		const std::string argument(argv[i]);
		auto requireValue = [&](const char *option) -> std::string
		{
			if (i + 1 >= argc)
			{
				throw std::invalid_argument(std::string(option) + " requires a value");
			}
			return argv[++i];
		};

		if (argument == "--threads")
		{
			options.threads = parsePositiveInt(requireValue("--threads"), "--threads");
		}
		else if (argument == "--memory-gib")
		{
			options.memoryGiB = parsePositiveDouble(requireValue("--memory-gib"), "--memory-gib");
		}
		else if (argument == "--seconds")
		{
			options.seconds = parsePositiveInt(requireValue("--seconds"), "--seconds");
		}
		else if (argument == "--allow-high-memory")
		{
			options.allowHighMemory = true;
		}
		else if (argument == "-h" || argument == "--help")
		{
			printUsage(argv[0]);
			std::exit(0);
		}
		else
		{
			throw std::invalid_argument("unknown option: " + argument);
		}
	}
	return options;
}


int main(int argc, char **argv)
{
	Options options;
	try
	{
		options = parseOptions(argc, argv);
	}
	catch (const std::exception &error)
	{
		std::cerr << "Error: " << error.what() << "\n\n";
		printUsage(argv[0]);
		return 2;
	}

	const long double requestedBytesExact =
		static_cast<long double>(options.memoryGiB) *
		static_cast<long double>(UINT64_C(1024) * 1024 * 1024);
	if (requestedBytesExact > static_cast<long double>(std::numeric_limits<std::size_t>::max()))
	{
		std::cerr << "Requested memory is too large for this process.\n";
		return 2;
	}

	std::size_t requestedBytes = static_cast<std::size_t>(requestedBytesExact);
	requestedBytes -= requestedBytes % sizeof(std::uint64_t);
	const std::size_t elementCount = requestedBytes / sizeof(std::uint64_t);
	if (elementCount < static_cast<std::size_t>(options.threads))
	{
		std::cerr << "Requested memory is too small for the thread count.\n";
		return 2;
	}

	const std::size_t availableBytes = parseAvailableMemoryBytes();
	if (!options.allowHighMemory && availableBytes > 0 &&
		static_cast<long double>(requestedBytes) >
			static_cast<long double>(availableBytes) * 0.7L)
	{
		std::cerr << "Refusing to allocate " << formatGiB(requestedBytes)
				  << " GiB: this exceeds 70% of MemAvailable ("
				  << formatGiB(availableBytes) << " GiB).\n"
				  << "Use a smaller --memory-gib value or explicitly pass --allow-high-memory.\n";
		return 2;
	}

	std::cout << "Standalone OpenMP system stress test\n"
			  << "  threads requested : " << options.threads << '\n'
			  << "  memory requested  : " << formatGiB(requestedBytes) << " GiB\n"
			  << "  MemAvailable      : "
			  << (availableBytes > 0 ? formatGiB(availableBytes) + " GiB" : "unknown") << '\n'
			  << "  duration           : at least " << options.seconds << " seconds\n";

	std::unique_ptr<std::uint64_t[]> memory(new (std::nothrow) std::uint64_t[elementCount]);
	if (!memory)
	{
		std::cerr << "Allocation failed before the stress test started.\n";
		return 3;
	}

	omp_set_dynamic(0);
	omp_set_num_threads(options.threads);
	const auto testStart = std::chrono::steady_clock::now();
	auto lastReport = testStart - std::chrono::seconds(1);
	std::uint64_t round = 0;
	std::uint64_t cumulativeChecksum = 0;

	do
	{
		std::size_t errorCount = 0;
		std::uint64_t roundChecksum = 0;
		std::atomic<std::size_t> firstBadIndex(std::numeric_limits<std::size_t>::max());
		std::uint64_t firstExpected = 0;
		std::uint64_t firstObserved = 0;
		int actualThreads = 0;
		const double roundStart = omp_get_wtime();

#pragma omp parallel num_threads(options.threads) reduction(+ : errorCount) reduction(^ : roundChecksum)
		{
			const int threadId = omp_get_thread_num();
			const int teamSize = omp_get_num_threads();
#pragma omp single
			actualThreads = teamSize;

			const std::size_t baseCount = elementCount / static_cast<std::size_t>(teamSize);
			const std::size_t remainder = elementCount % static_cast<std::size_t>(teamSize);
			const std::size_t begin = static_cast<std::size_t>(threadId) * baseCount +
				std::min(static_cast<std::size_t>(threadId), remainder);
			const std::size_t count = baseCount +
				(static_cast<std::size_t>(threadId) < remainder ? 1 : 0);
			const std::size_t end = begin + count;

			for (std::size_t index = begin; index < end; ++index)
			{
				memory[index] = expectedValue(index, round, threadId);
			}

			for (std::size_t index = begin; index < end; ++index)
			{
				const std::uint64_t expected = expectedValue(index, round, threadId);
				const std::uint64_t observed = memory[index];
				roundChecksum ^= mix64(observed + static_cast<std::uint64_t>(index));
				if (observed != expected)
				{
					++errorCount;
					std::size_t noError = std::numeric_limits<std::size_t>::max();
					if (firstBadIndex.compare_exchange_strong(noError, index))
					{
						firstExpected = expected;
						firstObserved = observed;
					}
				}
			}
		}

		if (actualThreads != options.threads)
		{
			std::cerr << "OpenMP created " << actualThreads << " threads instead of "
					  << options.threads << ".\n";
			return 4;
		}
		if (errorCount != 0)
		{
			std::cerr << "MEMORY ERROR: " << errorCount << " mismatched words; first index "
					  << firstBadIndex.load() << ", expected 0x" << std::hex << firstExpected
					  << ", observed 0x" << firstObserved << std::dec << ".\n";
			return 5;
		}

		cumulativeChecksum ^= roundChecksum;
		++round;
		const double roundSeconds = omp_get_wtime() - roundStart;
		const double movedGiB = static_cast<double>(requestedBytes) * 2.0 /
			static_cast<double>(UINT64_C(1024) * 1024 * 1024);
		const auto now = std::chrono::steady_clock::now();
		if (round == 1 || now - lastReport >= std::chrono::seconds(1))
		{
			std::cout << "  round " << round << " passed in " << std::fixed << std::setprecision(2)
					  << roundSeconds << " s (" << movedGiB / roundSeconds
					  << " GiB/s read+write), checksum 0x" << std::hex << roundChecksum
					  << std::dec << '\n'
					  << std::flush;
			lastReport = now;
		}
	}
	while (std::chrono::duration_cast<std::chrono::seconds>(
			   std::chrono::steady_clock::now() - testStart)
			   .count() < options.seconds);

	std::cout << "PASS: " << round << " complete round(s), no memory mismatches; final checksum 0x"
			  << std::hex << cumulativeChecksum << std::dec << '\n';
	return 0;
}
