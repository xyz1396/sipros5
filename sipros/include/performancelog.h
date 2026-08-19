#ifndef SIPROS_PERFORMANCE_LOG_H
#define SIPROS_PERFORMANCE_LOG_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

#include <omp.h>

namespace sipros
{

struct PerformanceTiming
{
	double wallSeconds = 0.0;
	double cpuSeconds = 0.0;

	double speedup() const
	{
		return wallSeconds > 0.0 ? cpuSeconds / wallSeconds : 0.0;
	}
};

double processCpuSeconds();

inline PerformanceTiming &operator+=(PerformanceTiming &left,
										 const PerformanceTiming &right)
{
	left.wallSeconds += right.wallSeconds;
	left.cpuSeconds += right.cpuSeconds;
	return left;
}

class PerformanceTimer
{
public:
	PerformanceTimer()
		: wallStart_(omp_get_wtime()), cpuStart_(processCpuSeconds())
	{
	}

	PerformanceTiming elapsed() const
	{
		const double cpuEnd = processCpuSeconds();
		PerformanceTiming result;
		result.wallSeconds = std::max(0.0, omp_get_wtime() - wallStart_);
		if (cpuStart_ >= 0.0 && cpuEnd >= 0.0)
		{
			result.cpuSeconds = std::max(0.0, cpuEnd - cpuStart_);
		}
		return result;
	}

private:
	double wallStart_ = 0.0;
	double cpuStart_ = 0.0;
};

inline std::string formatPerformanceSeconds(double seconds)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(seconds < 1.0 ? 3 : 2)
		<< seconds << "s";
	return out.str();
}

inline std::string formatPerformanceSpeedup(const PerformanceTiming &timing)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(2) << timing.speedup() << "x";
	return out.str();
}

inline std::string formatPerformanceCount(uint64_t value)
{
	std::string result = std::to_string(value);
	for (std::ptrdiff_t position =
			 static_cast<std::ptrdiff_t>(result.size()) - 3;
		 position > 0; position -= 3)
	{
		result.insert(static_cast<size_t>(position), 1, ',');
	}
	return result;
}

inline void printPerformanceRule(std::ostream &out)
{
	out << "  " << std::string(66, '-') << "\n";
}

inline void printPerformanceHeader(std::ostream &out,
								   const std::string &title,
								   int threadCount)
{
	out << "\n" << title << "\n"
		<< "  Parallel speedup = process CPU time / wall time"
		<< "; 1.00x is one fully busy CPU thread\n"
		<< "  OpenMP threads: " << threadCount
		<< " (ideal upper bound: " << threadCount << ".00x)\n"
		<< "  " << std::left << std::setw(33) << "Stage"
		<< std::right << std::setw(11) << "Wall time"
		<< std::setw(11) << "CPU time"
		<< std::setw(11) << "Speedup" << "\n"
		<< "  " << std::string(66, '-') << "\n";
}

inline void printPerformanceStage(std::ostream &out,
								  const std::string &stage,
								  const PerformanceTiming &timing,
								  const std::string &detail = std::string())
{
	out << "  " << std::left << std::setw(33) << stage
		<< std::right << std::setw(11) << formatPerformanceSeconds(timing.wallSeconds)
		<< std::setw(11) << formatPerformanceSeconds(timing.cpuSeconds)
		<< std::setw(11) << formatPerformanceSpeedup(timing);
	if (!detail.empty())
	{
		out << "  " << detail;
	}
	out << "\n";
}

inline void printPerformanceFooter(std::ostream &out)
{
	printPerformanceRule(out);
	out.flush();
}

} // namespace sipros

#endif // SIPROS_PERFORMANCE_LOG_H
