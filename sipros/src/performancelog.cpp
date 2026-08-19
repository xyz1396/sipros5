#include "performancelog.h"

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace sipros
{

double processCpuSeconds()
{
#if defined(_WIN32)
	FILETIME creationTime{};
	FILETIME exitTime{};
	FILETIME kernelTime{};
	FILETIME userTime{};
	if (!GetProcessTimes(GetCurrentProcess(), &creationTime, &exitTime,
						 &kernelTime, &userTime))
	{
		return -1.0;
	}
	ULARGE_INTEGER kernel{};
	ULARGE_INTEGER user{};
	kernel.LowPart = kernelTime.dwLowDateTime;
	kernel.HighPart = kernelTime.dwHighDateTime;
	user.LowPart = userTime.dwLowDateTime;
	user.HighPart = userTime.dwHighDateTime;
	return static_cast<double>(kernel.QuadPart + user.QuadPart) / 1.0e7;
#else
	const std::clock_t value = std::clock();
	return value == static_cast<std::clock_t>(-1)
		? -1.0
		: static_cast<double>(value) / static_cast<double>(CLOCKS_PER_SEC);
#endif
}

} // namespace sipros
