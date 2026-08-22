#include "workflow.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#ifdef _WIN32
#  define NOMINMAX
#  include <windows.h>
#else
#  include <cerrno>
#  include <csignal>
#  include <fcntl.h>
#  include <sched.h>
#  include <sys/types.h>
#  include <sys/wait.h>
#  include <unistd.h>
extern char** environ;
#endif

namespace siproswf {

static std::string lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

#ifdef _WIN32
static std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}
#endif

int parse_integer(const std::string& option, const std::string& value) {
    try {
        std::size_t consumed = 0;
        const int result = std::stoi(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("trailing text");
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(option + " requires an integer, got '" + value + "'");
    }
}

double parse_number(const std::string& option, const std::string& value) {
    try {
        std::size_t consumed = 0;
        const double result = std::stod(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("trailing text");
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(option + " requires a number, got '" + value + "'");
    }
}

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t when = std::chrono::system_clock::to_time_t(now);
    std::tm local{};
#ifdef _WIN32
    localtime_s(&local, &when);
#else
    localtime_r(&when, &local);
#endif
    std::ostringstream output;
    output << std::put_time(&local, "%Y-%m-%d %H:%M:%S");
    return output.str();
}

#ifdef _WIN32
std::wstring widen(const std::string& value) {
    if (value.empty()) return {};
    const int size = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
                                         static_cast<int>(value.size()), nullptr, 0);
    if (size <= 0) throw std::runtime_error("Invalid UTF-8 text in process argument");
    std::wstring result(static_cast<std::size_t>(size), L'\0');
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
                        static_cast<int>(value.size()), result.data(), size);
    return result;
}

std::string narrow(const std::wstring& value) {
    if (value.empty()) return {};
    const int size = WideCharToMultiByte(CP_UTF8, 0, value.data(),
                                         static_cast<int>(value.size()), nullptr, 0,
                                         nullptr, nullptr);
    std::string result(static_cast<std::size_t>(size), '\0');
    WideCharToMultiByte(CP_UTF8, 0, value.data(), static_cast<int>(value.size()),
                        result.data(), size, nullptr, nullptr);
    return result;
}

std::wstring quote_windows(const std::wstring& value) {
    if (!value.empty() && value.find_first_of(L" \t\n\v\"") == std::wstring::npos) {
        return value;
    }
    std::wstring result = L"\"";
    std::size_t backslashes = 0;
    for (const wchar_t c : value) {
        if (c == L'\\') {
            ++backslashes;
        } else if (c == L'\"') {
            result.append(backslashes * 2 + 1, L'\\');
            result.push_back(L'\"');
            backslashes = 0;
        } else {
            result.append(backslashes, L'\\');
            backslashes = 0;
            result.push_back(c);
        }
    }
    result.append(backslashes * 2, L'\\');
    result.push_back(L'\"');
    return result;
}

std::vector<wchar_t> windows_environment(
    const std::map<std::string, std::string>& updates,
    const std::wstring& path_prefix = {}) {
    struct Entry { std::wstring name; std::wstring value; };
    std::map<std::wstring, Entry> entries;
    auto canonical = [](std::wstring name) {
        std::transform(name.begin(), name.end(), name.begin(), ::towupper);
        return name;
    };
    LPWCH block = GetEnvironmentStringsW();
    if (block == nullptr) throw std::runtime_error("GetEnvironmentStringsW failed");
    for (const wchar_t* cursor = block; *cursor != L'\0'; cursor += std::wcslen(cursor) + 1) {
        const std::wstring item(cursor);
        const auto equals = item.find(L'=', item.front() == L'=' ? 1 : 0);
        if (equals == std::wstring::npos) continue;
        Entry entry{item.substr(0, equals), item.substr(equals + 1)};
        entries[canonical(entry.name)] = std::move(entry);
    }
    FreeEnvironmentStringsW(block);
    for (const auto& [key, value] : updates) {
        Entry entry{widen(key), widen(value)};
        entries[canonical(entry.name)] = std::move(entry);
    }
    if (!path_prefix.empty()) {
        const auto found = entries.find(L"PATH");
        if (found == entries.end()) {
            entries.emplace(L"PATH", Entry{L"PATH", path_prefix});
        } else {
            found->second.value = path_prefix + L';' + found->second.value;
        }
    }
    std::vector<wchar_t> result;
    for (const auto& [_, entry] : entries) {
        result.insert(result.end(), entry.name.begin(), entry.name.end());
        result.push_back(L'=');
        result.insert(result.end(), entry.value.begin(), entry.value.end());
        result.push_back(L'\0');
    }
    result.push_back(L'\0');
    return result;
}

std::string windows_error(const std::string& context) {
    const DWORD code = GetLastError();
    LPWSTR buffer = nullptr;
    FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                       FORMAT_MESSAGE_IGNORE_INSERTS,
                   nullptr, code, 0, reinterpret_cast<LPWSTR>(&buffer), 0, nullptr);
    const std::string message = buffer ? narrow(buffer) : std::string("unknown error");
    if (buffer) LocalFree(buffer);
    return context + " (Windows error " + std::to_string(code) + "): " + trim_copy(message);
}
#endif

void log_process_chunk(Logger& logger, std::string& pending,
                       const char* data, std::size_t size) {
    pending.append(data, size);
    std::size_t begin = 0;
    while (begin < pending.size()) {
        const std::size_t end = pending.find_first_of("\r\n", begin);
        if (end == std::string::npos) break;
        if (end > begin) logger.info(pending.substr(begin, end - begin));
        begin = end + 1;
        if (pending[end] == '\r' && begin < pending.size() &&
            pending[begin] == '\n') {
            ++begin;
        }
    }
    pending.erase(0, begin);
}

void flush_process_output(Logger& logger, std::string& pending) {
    if (!pending.empty()) logger.info(pending);
    pending.clear();
}

constexpr const char kLogRule[] =
    "================================================================================";

std::string clean_log_message(std::string message) {
    std::replace(message.begin(), message.end(), '\r', '\n');
    return message;
}

Logger::Logger(const std::filesystem::path& log_path, Sink sink) : sink_(std::move(sink)) {
    if (!log_path.empty()) {
        const std::filesystem::path parent = log_path.parent_path();
        if (!parent.empty()) {
            std::error_code error;
            std::filesystem::create_directories(parent, error);
            if (error) {
                throw std::runtime_error("Unable to create workflow log directory " +
                                         parent.string() + ": " + error.message());
            }
        }
        stream_.open(log_path, std::ios::out | std::ios::trunc);
        if (!stream_) throw std::runtime_error("Unable to create workflow log: " + log_path.string());
        stream_ << kLogRule << '\n'
                << "SIPROS WORKFLOW LOG" << '\n'
                << "Session started: " << timestamp() << '\n'
                << kLogRule << '\n';
        stream_.flush();
    }
}

Logger::~Logger() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (stream_) {
        stream_ << kLogRule << '\n'
                << "Session closed : " << timestamp() << '\n'
                << kLogRule << '\n';
        stream_.flush();
    }
}

void Logger::info(const std::string& message) { write("INFO", message); }
void Logger::warning(const std::string& message) { write("WARN", message); }
void Logger::error(const std::string& message) { write("ERROR", message); }
void Logger::debug(const std::string& message) { write("DEBUG", message); }

void Logger::write(const char* level, const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::istringstream lines(clean_log_message(message));
    std::string line;
    bool wrote = false;
    const std::string when = timestamp();
    while (std::getline(lines, line)) {
        const auto last = line.find_last_not_of(" \t");
        if (last == std::string::npos) continue;
        line.erase(last + 1);
        std::ostringstream formatted_line;
        formatted_line << when << " | " << std::left << std::setw(5) << level
                       << " | " << line;
        const std::string formatted = formatted_line.str();
        if (stream_) stream_ << formatted << '\n';
        std::cerr << formatted << '\n';
        if (sink_) sink_(formatted);
        wrote = true;
    }
    if (stream_ && wrote) stream_.flush();
}

int ThreadAllocation::peak_threads() const {
    int total = 0;
    for (int i = 0; i < worker_count && i < static_cast<int>(task_threads.size()); ++i) {
        total += task_threads[static_cast<std::size_t>(i)];
    }
    return total;
}

int physical_cpu_count() {
#ifdef _WIN32
    DWORD length = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);
    if (length > 0 && GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
        std::vector<unsigned char> data(length);
        auto* info = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(data.data());
        if (GetLogicalProcessorInformationEx(RelationProcessorCore, info, &length)) {
            int cores = 0;
            DWORD offset = 0;
            while (offset < length) {
                auto* current = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(data.data() + offset);
                if (current->Relationship == RelationProcessorCore) ++cores;
                if (current->Size == 0) break;
                offset += current->Size;
            }
            if (cores > 0) return cores;
        }
    }
#else
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    const bool affinity_known = sched_getaffinity(0, sizeof(allowed), &allowed) == 0;
    std::set<std::pair<std::string, std::string>> cores;
    const std::filesystem::path cpu_root("/sys/devices/system/cpu");
    std::error_code error;
    if (std::filesystem::is_directory(cpu_root, error)) {
        for (const auto& entry : std::filesystem::directory_iterator(cpu_root, error)) {
            const std::string name = entry.path().filename().string();
            if (name.size() <= 3 || name.rfind("cpu", 0) != 0 ||
                !std::all_of(name.begin() + 3, name.end(), [](unsigned char c) {
                    return std::isdigit(c) != 0;
                })) continue;
            const int cpu = std::stoi(name.substr(3));
            if (affinity_known && (cpu >= CPU_SETSIZE || !CPU_ISSET(cpu, &allowed))) continue;
            std::ifstream package(entry.path() / "topology/physical_package_id");
            std::ifstream core(entry.path() / "topology/core_id");
            std::string package_id, core_id;
            if (package >> package_id && core >> core_id) cores.emplace(package_id, core_id);
        }
    }
    if (!cores.empty()) return static_cast<int>(cores.size());
#endif
    const unsigned int fallback = std::thread::hardware_concurrency();
    return std::max(1, static_cast<int>(fallback == 0 ? 1 : fallback));
}

int effective_thread_count(int requested, int available) {
    if (requested < 0) throw std::invalid_argument("Thread count must be non-negative");
    const int cpu_count = available > 0 ? available : physical_cpu_count();
    return requested == 0 ? cpu_count : std::min(requested, cpu_count);
}

ThreadAllocation allocate_threads(int total_threads, int task_count,
                                  int minimum_threads_per_task) {
    if (total_threads <= 0) throw std::invalid_argument("total_threads must be positive");
    if (task_count < 0) throw std::invalid_argument("task_count must be non-negative");
    if (minimum_threads_per_task <= 0) {
        throw std::invalid_argument("minimum_threads_per_task must be positive");
    }
    ThreadAllocation allocation;
    if (task_count == 0) return allocation;
    const int effective_minimum = std::min(total_threads, minimum_threads_per_task);
    allocation.worker_count = std::min(task_count, std::max(1, total_threads / effective_minimum));
    const int per_task = total_threads / allocation.worker_count;
    const int remainder = total_threads % allocation.worker_count;
    allocation.task_threads.reserve(static_cast<std::size_t>(task_count));
    for (int i = 0; i < allocation.worker_count; ++i) {
        allocation.task_threads.push_back(per_task + (i < remainder ? 1 : 0));
    }
    while (static_cast<int>(allocation.task_threads.size()) < task_count) {
        allocation.task_threads.push_back(per_task);
    }
    return allocation;
}

std::map<std::string, std::string> thread_environment(int thread_count) {
    if (thread_count <= 0) throw std::invalid_argument("thread_count must be positive");
    const std::string value = std::to_string(thread_count);
    return {
        {"OMP_NUM_THREADS", value}, {"OMP_THREAD_LIMIT", value},
        {"OMP_MAX_ACTIVE_LEVELS", "1"}, {"OMP_DYNAMIC", "FALSE"},
        {"OPENBLAS_NUM_THREADS", "1"}, {"MKL_NUM_THREADS", "1"},
        {"NUMEXPR_NUM_THREADS", "1"}, {"VECLIB_MAXIMUM_THREADS", "1"},
        {"BLIS_NUM_THREADS", "1"}, {"GOMAXPROCS", value},
        {"DOTNET_PROCESSOR_COUNT", value}, {"COMPlus_ProcessorCount", value},
    };
}

std::string quote_display_argument(const std::string& value) {
    if (!value.empty() && value.find_first_of(" \t\n\r\"'") == std::string::npos) return value;
    std::string result = "\"";
    for (const char c : value) {
        if (c == '\\' || c == '\"') result.push_back('\\');
        result.push_back(c);
    }
    result.push_back('\"');
    return result;
}

std::string display_command(const Command& command) {
    std::ostringstream result;
    result << quote_display_argument(command.executable.string());
    for (const auto& argument : command.arguments) result << ' ' << quote_display_argument(argument);
    return result.str();
}

int run_process(const Command& command, Logger& logger, std::atomic_bool& cancelled) {
    if (cancelled.load()) throw std::runtime_error("Workflow cancelled");
    logger.info("Running process (" + std::to_string(command.cpu_threads) + " CPU " +
                (command.cpu_threads == 1 ? "thread): " : "threads): ") + display_command(command));
    const auto started = std::chrono::steady_clock::now();
#ifdef _WIN32
    SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    HANDLE read_pipe = nullptr;
    HANDLE write_pipe = nullptr;
    if (!CreatePipe(&read_pipe, &write_pipe, &security, 0)) {
        throw std::runtime_error(windows_error("CreatePipe failed"));
    }
    SetHandleInformation(read_pipe, HANDLE_FLAG_INHERIT, 0);
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = write_pipe;
    startup.hStdError = write_pipe;
    startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);

    std::wstring executable = command.executable.wstring();
    std::wstring command_line = quote_windows(executable);
    for (const auto& argument : command.arguments) {
        command_line.push_back(L' ');
        command_line += quote_windows(widen(argument));
    }
    std::vector<wchar_t> mutable_command(command_line.begin(), command_line.end());
    mutable_command.push_back(L'\0');
    std::wstring packaged_library_directory;
    std::error_code library_error;
    const std::filesystem::path library_candidate = command.executable.parent_path() / L"lib";
    if (std::filesystem::is_directory(library_candidate, library_error)) {
        packaged_library_directory = library_candidate.wstring();
    }
    auto environment = windows_environment(command.environment,
                                           packaged_library_directory);
    std::wstring working_directory;
    if (command.working_directory) working_directory = command.working_directory->wstring();

    PROCESS_INFORMATION process{};
    const BOOL created = CreateProcessW(
        executable.c_str(), mutable_command.data(), nullptr, nullptr, TRUE,
        CREATE_UNICODE_ENVIRONMENT | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW,
        environment.data(), command.working_directory ? working_directory.c_str() : nullptr,
        &startup, &process);
    CloseHandle(write_pipe);
    if (!created) {
        CloseHandle(read_pipe);
        throw std::runtime_error(windows_error("Unable to start " + command.executable.string()));
    }
    HANDLE job = CreateJobObjectW(nullptr, nullptr);
    if (job != nullptr) {
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
        limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        SetInformationJobObject(job, JobObjectExtendedLimitInformation, &limits, sizeof(limits));
        if (!AssignProcessToJobObject(job, process.hProcess)) {
            CloseHandle(job);
            job = nullptr;
        }
    }
    CloseHandle(process.hThread);
    char buffer[8192];
    std::string pending_output;
    bool running = true;
    while (running) {
        DWORD available = 0;
        while (PeekNamedPipe(read_pipe, nullptr, 0, nullptr, &available, nullptr) && available > 0) {
            DWORD read = 0;
            const DWORD amount = std::min<DWORD>(available, sizeof(buffer));
            if (!ReadFile(read_pipe, buffer, amount, &read, nullptr) || read == 0) break;
            log_process_chunk(logger, pending_output, buffer, read);
        }
        if (cancelled.load()) {
            if (job) TerminateJobObject(job, 130);
            else TerminateProcess(process.hProcess, 130);
        }
        const DWORD wait = WaitForSingleObject(process.hProcess, 20);
        running = wait == WAIT_TIMEOUT;
        if (wait == WAIT_FAILED) break;
    }
    DWORD available = 0;
    while (PeekNamedPipe(read_pipe, nullptr, 0, nullptr, &available, nullptr) && available > 0) {
        DWORD read = 0;
        if (!ReadFile(read_pipe, buffer, std::min<DWORD>(available, sizeof(buffer)), &read, nullptr) || read == 0) break;
        log_process_chunk(logger, pending_output, buffer, read);
    }
    flush_process_output(logger, pending_output);
    DWORD exit_code = 1;
    GetExitCodeProcess(process.hProcess, &exit_code);
    CloseHandle(read_pipe);
    CloseHandle(process.hProcess);
    if (job) CloseHandle(job);
    const int result = static_cast<int>(exit_code);
#else
    int pipes[2];
    if (pipe(pipes) != 0) throw std::runtime_error("pipe failed: " + std::string(std::strerror(errno)));
    std::map<std::string, std::string> environment;
    for (char** item = environ; item && *item; ++item) {
        std::string entry(*item);
        const auto equals = entry.find('=');
        if (equals != std::string::npos) environment[entry.substr(0, equals)] = entry.substr(equals + 1);
    }
    for (const auto& update : command.environment) environment[update.first] = update.second;
    std::vector<std::string> argument_storage;
    argument_storage.reserve(command.arguments.size() + 1);
    argument_storage.push_back(command.executable.string());
    argument_storage.insert(argument_storage.end(), command.arguments.begin(), command.arguments.end());
    std::vector<char*> arguments;
    for (auto& value : argument_storage) arguments.push_back(value.data());
    arguments.push_back(nullptr);
    std::vector<std::string> environment_storage;
    for (const auto& entry : environment) environment_storage.push_back(entry.first + "=" + entry.second);
    std::vector<char*> environment_values;
    for (auto& value : environment_storage) environment_values.push_back(value.data());
    environment_values.push_back(nullptr);

    const pid_t pid = fork();
    if (pid < 0) {
        close(pipes[0]); close(pipes[1]);
        throw std::runtime_error("fork failed: " + std::string(std::strerror(errno)));
    }
    if (pid == 0) {
        setpgid(0, 0);
        dup2(pipes[1], STDOUT_FILENO);
        dup2(pipes[1], STDERR_FILENO);
        close(pipes[0]); close(pipes[1]);
        if (command.working_directory && chdir(command.working_directory->c_str()) != 0) {
            dprintf(STDERR_FILENO, "Unable to change working directory: %s\n", std::strerror(errno));
            _exit(127);
        }
        execve(command.executable.c_str(), arguments.data(), environment_values.data());
        dprintf(STDERR_FILENO, "Unable to start %s: %s\n", command.executable.c_str(), std::strerror(errno));
        _exit(127);
    }
    close(pipes[1]);
    fcntl(pipes[0], F_SETFL, fcntl(pipes[0], F_GETFL) | O_NONBLOCK);
    int status = 0;
    bool sent_term = false;
    auto cancel_started = std::chrono::steady_clock::time_point{};
    char buffer[8192];
    std::string pending_output;
    while (true) {
        const ssize_t count = read(pipes[0], buffer, sizeof(buffer));
        if (count > 0) {
            log_process_chunk(logger, pending_output, buffer,
                              static_cast<std::size_t>(count));
        }
        const pid_t waited = waitpid(pid, &status, WNOHANG);
        if (waited == pid) break;
        if (waited < 0 && errno != EINTR) break;
        if (cancelled.load() && !sent_term) {
            kill(-pid, SIGTERM);
            sent_term = true;
            cancel_started = std::chrono::steady_clock::now();
        } else if (sent_term && std::chrono::steady_clock::now() - cancel_started > std::chrono::seconds(2)) {
            kill(-pid, SIGKILL);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    while (true) {
        const ssize_t count = read(pipes[0], buffer, sizeof(buffer));
        if (count <= 0) break;
        log_process_chunk(logger, pending_output, buffer,
                          static_cast<std::size_t>(count));
    }
    flush_process_output(logger, pending_output);
    close(pipes[0]);
    const int result = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
#endif
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (cancelled.load()) throw std::runtime_error("Workflow cancelled");
    if (result != 0) {
        throw std::runtime_error("Process exited with status " + std::to_string(result) + ": " +
                                 display_command(command));
    }
    std::ostringstream timing;
    timing << "Process completed in " << std::fixed << std::setprecision(3) << seconds << " s";
    logger.info(timing.str());
    return result;
}

std::string mode_name(WorkflowMode mode) {
    switch (mode) {
        case WorkflowMode::RegularFasta: return "regular FASTA";
        case WorkflowMode::SipFasta: return "SIP FASTA";
        case WorkflowMode::FastSip: return "fast SIP";
        case WorkflowMode::Spectra: return "spectra search";
        default: return "automatic";
    }
}

ParseResult parse_arguments(const std::vector<std::string>& arguments) {
    ParseResult result;
    if (arguments.empty()) { result.launch_gui = true; return result; }
    int explicit_modes = 0;
    auto set_mode = [&](WorkflowMode mode) { result.options.mode = mode; ++explicit_modes; };
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        const std::string& option = arguments[i];
        auto value = [&]() -> std::string {
            if (++i >= arguments.size()) throw std::runtime_error(option + " requires a value");
            return arguments[i];
        };
        if (option == "-h" || option == "--help") result.show_help = true;
        else if (option == "--gui") result.launch_gui = true;
        else if (option == "--regular-fasta-search") set_mode(WorkflowMode::RegularFasta);
        else if (option == "--sip-fasta-search") set_mode(WorkflowMode::SipFasta);
        else if (option == "--fast-sip-search") set_mode(WorkflowMode::FastSip);
        else if (option == "-i" || option == "--input") result.options.input = value();
        else if (option == "-f" || option == "--fasta") result.options.fasta = value();
        else if (option == "-o" || option == "--output") result.options.output = value();
        else if (option == "--ptm") result.options.ptms.push_back(value());
        else if (option == "--fixed-ptm") result.options.fixed_ptms.push_back(value());
        else if (option == "--max-ptm-count") result.options.max_ptm_count = parse_integer(option, value());
        else if (option == "--toleranceMS1") result.options.tolerance_ms1 = parse_number(option, value());
        else if (option == "--toleranceMS2") result.options.tolerance_ms2 = parse_number(option, value());
        else if (option == "-e" || option == "--element") result.options.element = value();
        else if (option == "-r" || option == "--range") result.options.sip_range = value();
        else if (option == "-p" || option == "--precision") result.options.precision = value();
        else if (option == "--psm-tsv") result.options.psm_tsv = value();
        else if (option == "--unlabeled-input") result.options.unlabeled_input = value();
        else if (option == "--spectra-dir") result.options.spectra_dir = value();
        else if (option == "-n" || option == "--nPrecursor") result.options.n_precursor = parse_integer(option, value());
        else if (option == "--product-top-isotopes") result.options.product_top_isotopes = parse_integer(option, value());
        else if (option == "-t" || option == "--thread") result.options.threads = parse_integer(option, value());
        else if (option == "--aerith-sample-parallelism") result.options.aerith_sample_parallelism = parse_integer(option, value());
        else if (option == "--topN" || option == "--top-psms-per-scan") result.options.top_psms_per_scan = parse_integer(option, value());
        else if (option == "--rt-tolerance") result.options.rt_tolerance = parse_number(option, value());
        else if (option == "--sfi-envelope-top-n") result.options.sfi_envelope_top_n = parse_integer(option, value());
        else if (option == "--mvh-cascade-top-n") result.options.mvh_cascade_top_n = parse_integer(option, value());
        else if (option == "--ignorePCT") result.options.ignore_pct = true;
        else if (option == "--negative_control") result.options.negative_control = value();
        else if (option == "--label_threshold") result.options.label_threshold = parse_number(option, value());
        else if (option == "--dryrun") result.options.dry_run = true;
        else throw std::runtime_error("Unknown option: " + option);
    }
    if (explicit_modes > 1) throw std::runtime_error("Search-mode switches are mutually exclusive");
    if (result.show_help) return result;
    if (result.launch_gui) {
        if (arguments.size() != 1) {
            throw std::runtime_error("--gui cannot be combined with headless workflow arguments");
        }
        return result;
    }
    validate_options(result.options, &result.warnings);
    return result;
}

void validate_options(WorkflowOptions& options, std::vector<std::string>* warnings) {
    if (options.input.empty()) throw std::runtime_error("--input is required");
    if (options.fasta.empty()) throw std::runtime_error("--fasta is required");
    if (options.output.empty()) throw std::runtime_error("--output is required");
    if (options.threads < 0) throw std::runtime_error("Thread number must be non-negative");
    if (options.aerith_sample_parallelism <= 0) throw std::runtime_error("--aerith-sample-parallelism must be positive");
    if (options.top_psms_per_scan <= 0) throw std::runtime_error("--topN must be positive");
    if (options.n_precursor <= 0) throw std::runtime_error("--nPrecursor must be positive");
    if (options.product_top_isotopes <= 0) throw std::runtime_error("--product-top-isotopes must be positive");
    if (options.rt_tolerance < 0) throw std::runtime_error("--rt-tolerance must be non-negative");
    if (options.sfi_envelope_top_n <= 0) throw std::runtime_error("--sfi-envelope-top-n must be positive");
    if (options.mvh_cascade_top_n <= 0) throw std::runtime_error("--mvh-cascade-top-n must be positive");
    if (options.max_ptm_count && *options.max_ptm_count < 0) throw std::runtime_error("--max-ptm-count must be non-negative");

    const bool spectra_arguments = options.psm_tsv || options.unlabeled_input || options.spectra_dir;
    if (spectra_arguments && options.mode != WorkflowMode::Auto && options.mode != WorkflowMode::Spectra) {
        throw std::runtime_error("Standalone spectra options cannot be combined with a FASTA search-mode switch");
    }
    if (spectra_arguments) options.mode = WorkflowMode::Spectra;
    if (options.mode == WorkflowMode::Spectra && !options.spectra_dir && (!options.psm_tsv || !options.unlabeled_input)) {
        throw std::runtime_error("Spectra search requires --psm-tsv and --unlabeled-input unless --spectra-dir is provided");
    }
    if (options.mode == WorkflowMode::Spectra && (!options.ptms.empty() || options.max_ptm_count)) {
        throw std::runtime_error("--ptm and --max-ptm-count apply only to FASTA search");
    }
    if (options.spectra_dir && !options.fixed_ptms.empty()) {
        throw std::runtime_error("--fixed-ptm cannot be used with --spectra-dir");
    }

    if (!options.element.empty() && options.element != "R" && options.element != "r") {
        // Normalize the supported isotope spellings explicitly.
        const std::string lowered = lower_copy(options.element);
        if (lowered == "c13") options.element = "C13";
        else if (lowered == "h2") options.element = "H2";
        else if (lowered == "n15") options.element = "N15";
        else if (lowered == "o18") options.element = "O18";
        else if (lowered == "s34") options.element = "S34";
        else throw std::runtime_error("--element must be one of C13, H2, N15, O18, or S34");
    } else if (!options.element.empty()) {
        options.element = "R";
    }

    if (options.mode == WorkflowMode::RegularFasta) {
        if (!options.negative_control.empty()) {
            throw std::runtime_error(
                "--negative_control applies only to SIP FASTA and fast SIP searches");
        }
        if (options.ignore_pct) {
            throw std::runtime_error(
                "--ignorePCT applies only to SIP FASTA and fast SIP searches");
        }
        if (!options.element.empty() && options.element != "R") {
            throw std::runtime_error("--regular-fasta-search cannot use a SIP --element");
        }
        if (options.sip_range || options.precision) {
            throw std::runtime_error("--range and --precision cannot be used with --regular-fasta-search");
        }
        options.element = "R";
    } else if (options.mode == WorkflowMode::SipFasta || options.mode == WorkflowMode::FastSip ||
               options.mode == WorkflowMode::Spectra) {
        if (options.element.empty()) options.element = "C13";
        if (options.element == "R") throw std::runtime_error(mode_name(options.mode) + " requires a SIP isotope");
    } else {
        if (options.element.empty() || options.element == "R") {
            options.element = "R";
            options.mode = WorkflowMode::RegularFasta;
        } else {
            options.mode = WorkflowMode::SipFasta;
        }
    }

    const int available = physical_cpu_count();
    if (options.threads > available && warnings) {
        warnings->push_back("Requested " + std::to_string(options.threads) + " threads, but only " +
                            std::to_string(available) + " physical CPU cores are available; using " +
                            std::to_string(available));
    }
    options.threads = effective_thread_count(options.threads, available);
    if (options.threads < kMinimumSiprosThreads && warnings) {
        warnings->push_back("The 8-thread Sipros minimum cannot fit within this workflow budget; jobs will run serially");
    }
}

void print_help(std::ostream& output) {
    output <<
        "Sipros workflow (siproswf)\n\n"
        "Usage:\n"
        "  siproswf                                      Open the ImGui interface\n"
        "  siproswf [MODE] -i INPUT -f FASTA -o OUTPUT [OPTIONS]\n\n"
        "Modes (mutually exclusive):\n"
        "  --regular-fasta-search   Regular target/decoy FASTA search\n"
        "  --sip-fasta-search       SIP target/decoy FASTA search\n"
        "  --fast-sip-search        Regular search, filtering, SFI generation, SIP search\n\n"
        "Core options:\n"
        "  -i, --input PATH         raw/.d/.d.zip/HDF5 path, directory, or comma list\n"
        "  -f, --fasta PATH         Target FASTA path\n"
        "  -o, --output PATH        Output directory\n"
        "  -e, --element ISOTOPE    C13, H2, N15, O18, or S34\n"
        "  -r, --range VALUE        SIP percentage/range/list (default 0-100)\n"
        "  -p, --precision VALUE    SIP percentage precision (default 1)\n"
        "  -t, --thread N           Total physical-core budget (0 = all physical cores)\n"
        "  -n, --nPrecursor N       Raxport precursor-envelope apex count (default 6)\n"
        "  --ptm VALUE              Variable PTM; repeat for multiple selections\n"
        "  --fixed-ptm VALUE        Fixed PTM; repeat for multiple selections\n"
        "  --max-ptm-count N        Maximum variable PTMs per peptide\n"
        "  --toleranceMS1 DA        MS1 tolerance (default 0.01)\n"
        "  --toleranceMS2 DA        MS2 tolerance (default 0.01)\n"
        "  --topN N                 PSM rows retained per scan (default 20)\n"
        "  --product-top-isotopes N Product isotope peaks (default 5)\n"
        "  --aerith-sample-parallelism N  Concurrent Aerith samples (default 3)\n"
        "  --negative_control LIST  Negative-control sample names\n"
        "  --label_threshold PCT    Negative-control label threshold (default 2)\n"
        "  --ignorePCT              Ignore SIP abundance features during filtering\n"
        "  --dryrun                 Validate and show workflow setup without child tools\n\n"
        "Standalone spectra search:\n"
        "  --psm-tsv PATH --unlabeled-input PATH, or --spectra-dir PATH\n"
        "  --rt-tolerance MINUTES --sfi-envelope-top-n N --mvh-cascade-top-n N\n\n"
        "Compatibility: without a mode switch, a SIP --element selects SIP FASTA; otherwise\n"
        "regular FASTA is selected. Use --gui to open the interface explicitly.\n";
}

ToolPaths locate_tools(const std::filesystem::path& executable_path) {
#ifdef _WIN32
    const std::string raxport_name = "Raxport-win-x64.exe";
    const std::string sipros_name = "sipros.exe";
    const std::string aerith_name = "aerith.exe";
#else
    const std::string raxport_name = "Raxport-linux-x64";
    const std::string sipros_name = "sipros";
    const std::string aerith_name = "aerith";
#endif
    std::error_code error;
    const std::filesystem::path executable_dir = std::filesystem::absolute(executable_path, error).parent_path();
    const std::vector<std::filesystem::path> roots = {
        executable_dir / "lib",
        executable_dir,
    };
    auto locate = [&](const char* override_name, const std::string& filename) {
        if (const char* override_value = std::getenv(override_name); override_value && *override_value) {
            return std::filesystem::path(override_value);
        }
        for (const auto& root : roots) {
            const std::filesystem::path candidate = root / filename;
            if (std::filesystem::exists(candidate, error)) return candidate;
        }
        return roots.front() / filename;
    };
    return {
        locate("SIPROSWF_RAXPORT", raxport_name),
        locate("SIPROSWF_SIPROS", sipros_name),
        locate("SIPROSWF_AERITH", aerith_name),
    };
}

} // namespace siproswf
