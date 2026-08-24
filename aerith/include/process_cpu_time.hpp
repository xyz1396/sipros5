#pragma once

#include <ctime>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace aerith {

// Keep the same units as std::clock() so existing timing calculations can
// divide elapsed ticks by CLOCKS_PER_SEC. On Windows, std::clock() measures
// wall time, so query the process user and kernel times directly instead.
using ProcessCpuTick = double;

inline ProcessCpuTick process_cpu_time_ticks() {
#ifdef _WIN32
    FILETIME creation_time{};
    FILETIME exit_time{};
    FILETIME kernel_time{};
    FILETIME user_time{};
    if (!GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time,
                         &kernel_time, &user_time)) {
        return 0.0;
    }
    ULARGE_INTEGER kernel{};
    kernel.LowPart = kernel_time.dwLowDateTime;
    kernel.HighPart = kernel_time.dwHighDateTime;
    ULARGE_INTEGER user{};
    user.LowPart = user_time.dwLowDateTime;
    user.HighPart = user_time.dwHighDateTime;
    constexpr double filetime_ticks_per_second = 10000000.0;
    return (static_cast<double>(kernel.QuadPart) +
            static_cast<double>(user.QuadPart)) *
           (static_cast<double>(CLOCKS_PER_SEC) /
            filetime_ticks_per_second);
#else
    return static_cast<ProcessCpuTick>(std::clock());
#endif
}

}  // namespace aerith
