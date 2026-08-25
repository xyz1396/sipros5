#include "workflow.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <future>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace siproswf {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> values;
    std::istringstream input(value);
    std::string item;
    while (std::getline(input, item, delimiter)) {
        item = trim(item);
        if (!item.empty()) values.push_back(item);
    }
    return values;
}

std::vector<std::string> split_preserving_empty(const std::string& value, char delimiter) {
    std::vector<std::string> values;
    std::size_t start = 0;
    while (true) {
        const std::size_t end = value.find(delimiter, start);
        values.push_back(value.substr(start, end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return values;
}

bool raw_input(const std::filesystem::path& path) {
    const std::string text = lower(path.string());
    return text.size() >= 4 &&
           (text.rfind(".raw") == text.size() - 4 ||
            text.rfind(".d") == text.size() - 2 ||
            (text.size() >= 6 && text.rfind(".d.zip") == text.size() - 6));
}

bool hdf5_input(const std::filesystem::path& path) {
    const std::string text = lower(path.string());
    return (text.size() >= 3 && text.rfind(".h5") == text.size() - 3) ||
           (text.size() >= 5 && text.rfind(".hdf5") == text.size() - 5);
}

std::string sample_base_name(const std::filesystem::path& path) {
    const std::string name = path.filename().string();
    const std::string lowered = lower(name);
    if (lowered.size() >= 6 && lowered.rfind(".d.zip") == lowered.size() - 6) {
        return name.substr(0, name.size() - 6);
    }
    return path.stem().string();
}

std::string number(double value) {
    std::ostringstream output;
    output << std::setprecision(12) << value;
    return output.str();
}

bool path_entry_exists(const std::filesystem::path& path) {
    std::error_code error;
    const auto status = std::filesystem::symlink_status(path, error);
    return !error && status.type() != std::filesystem::file_type::not_found;
}

std::string environment_value(const char* name) {
    const char* value = std::getenv(name);
    return value ? value : "";
}

Command child_command(const std::filesystem::path& executable, std::vector<std::string> arguments,
                      int threads, std::map<std::string, std::string> extra = {}) {
    Command command;
    command.executable = executable;
    command.arguments = std::move(arguments);
    command.cpu_threads = threads;
    command.environment = thread_environment(threads);
    command.environment.insert(extra.begin(), extra.end());
    return command;
}

template <typename State>
void check_cancelled(State& state) {
    if (state.cancelled->load()) throw std::runtime_error("Workflow cancelled");
}

template <typename State>
void run_command(State& state, Command command) {
    check_cancelled(state);
    run_process(command, *state.logger, *state.cancelled);
}

template <typename State>
void log_allocation(State& state, const std::string& phase,
                    const ThreadAllocation& allocation) {
    if (allocation.worker_count == 0) {
        state.logger->info(phase + ": no jobs");
        return;
    }
    const auto [minimum, maximum] = std::minmax_element(
        allocation.task_threads.begin(), allocation.task_threads.end());
    const std::string per_job = *minimum == *maximum
        ? std::to_string(*minimum)
        : std::to_string(*minimum) + "-" + std::to_string(*maximum);
    state.logger->info(
        phase + ": up to " + std::to_string(allocation.worker_count) +
        " concurrent processes; " + per_job + " CPU " +
        ((*minimum == 1 && *maximum == 1) ? "thread" : "threads") +
        " per process; " + std::to_string(allocation.peak_threads()) + "/" +
        std::to_string(state.options.threads) + " thread budget at peak");
}

std::string gibibytes(std::uint64_t bytes) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(2)
           << static_cast<long double>(bytes) /
                  static_cast<long double>(std::uint64_t{1024} * 1024 * 1024);
    return output.str();
}

template <typename State>
ThreadAllocation sipros_search_allocation(
    State& state, int task_count, int minimum_threads_per_task,
    const std::string& phase) {
    const std::uint64_t usable_memory = usable_memory_bytes();
    const bool serial = serialize_sipros_searches(
        state.options.threads, usable_memory);
    if (serial && task_count > 1) {
        std::vector<std::string> reasons;
        if (state.options.threads <= kMinimumSiprosThreads) {
            reasons.push_back(
                "thread budget " + std::to_string(state.options.threads) +
                " <= " + std::to_string(kMinimumSiprosThreads));
        }
        if (usable_memory > 0 &&
            usable_memory <= kParallelSiprosMinimumUsableMemoryBytes) {
            reasons.push_back(
                "usable RAM " + gibibytes(usable_memory) + " GiB <= 32 GiB");
        }
        std::ostringstream message;
        message << phase << ": target and decoy searches will run serially because ";
        for (std::size_t index = 0; index < reasons.size(); ++index) {
            if (index != 0) message << " and ";
            message << reasons[index];
        }
        state.logger->info(message.str());
    }
    return allocate_sipros_search_threads(
        state.options.threads, task_count, minimum_threads_per_task,
        usable_memory);
}

template <typename State>
void run_parallel(State& state, std::vector<Command> commands,
                  const ThreadAllocation& allocation) {
    if (allocation.worker_count == 0) return;
    for (std::size_t offset = 0; offset < commands.size();
         offset += static_cast<std::size_t>(allocation.worker_count)) {
        check_cancelled(state);
        const std::size_t count = std::min<std::size_t>(
            static_cast<std::size_t>(allocation.worker_count), commands.size() - offset);
        std::vector<std::future<void>> futures;
        futures.reserve(count);
        for (std::size_t index = 0; index < count; ++index) {
            Command command = std::move(commands[offset + index]);
            const std::size_t allocation_index = offset + index;
            command.cpu_threads = allocation.task_threads[allocation_index];
            command.environment = thread_environment(command.cpu_threads);
            futures.push_back(std::async(std::launch::async, [&state, command = std::move(command)]() mutable {
                run_command(state, std::move(command));
            }));
        }
        std::exception_ptr failure;
        for (auto& future : futures) {
            try { future.get(); }
            catch (...) { if (!failure) failure = std::current_exception(); }
        }
        if (failure) std::rethrow_exception(failure);
    }
}

void add_ptm_arguments(std::vector<std::string>& arguments,
                       const WorkflowOptions& options, bool include_variable = true) {
    if (include_variable) {
        for (const auto& ptm : options.ptms) {
            arguments.push_back("--ptm"); arguments.push_back(ptm);
        }
    }
    for (const auto& ptm : options.fixed_ptms) {
        arguments.push_back("--fixed-ptm"); arguments.push_back(ptm);
    }
    if (include_variable && options.max_ptm_count) {
        arguments.push_back("--max-ptm-count");
        arguments.push_back(std::to_string(*options.max_ptm_count));
    }
}

void add_tolerance_arguments(std::vector<std::string>& arguments,
                             const WorkflowOptions& options) {
    arguments.push_back("--tolerance-ms1"); arguments.push_back(number(options.tolerance_ms1));
    arguments.push_back("--tolerance-ms2"); arguments.push_back(number(options.tolerance_ms2));
}

struct Workflow::SearchState {
    WorkflowOptions options;
    WorkflowMode mode = WorkflowMode::RegularFasta;
    ToolPaths tools;
    Logger* logger = nullptr;
    std::atomic_bool* cancelled = nullptr;
    std::filesystem::path output;
    std::filesystem::path fasta;
    std::filesystem::path decoy;
    std::optional<std::filesystem::path> psm_tsv;
    std::optional<std::filesystem::path> unlabeled_input;
    std::optional<std::filesystem::path> spectra_dir;
    std::filesystem::path generated_spectra_dir;
    bool stage_hdf5_copies = true;
    std::string decoy_prefix;
    std::vector<std::filesystem::path> raw_files;
    std::vector<std::filesystem::path> hdf5_input_files;
    std::vector<std::string> base_names;
    std::vector<std::string> raw_base_names;
    std::vector<std::string> hdf5_base_names;
    std::map<std::string, std::filesystem::path> hdf5_paths;
};

template <typename State>
std::filesystem::path expected_hdf5(const State& state, const std::string& base) {
    return state.output / base / (base + ".h5");
}

bool complete_hdf5(const std::filesystem::path& path) {
    std::error_code error;
    return std::filesystem::exists(path, error) && std::filesystem::file_size(path, error) >= 1024ULL * 1024ULL;
}

template <typename State>
std::vector<std::filesystem::path> input_entries(State& state, const std::string& input) {
    if (input.find(',') != std::string::npos) {
        const auto entries = split(input, ',');
        state.logger->info("Input is a comma-separated file list with " +
                           std::to_string(entries.size()) + " entries");
        std::vector<std::filesystem::path> result;
        for (const auto& entry : entries) result.emplace_back(entry);
        return result;
    }
    const std::filesystem::path path(input);
    if (raw_input(path) || hdf5_input(path)) return {path};
    std::error_code error;
    if (std::filesystem::is_directory(path, error)) {
        state.logger->info(input + " is a directory");
        std::vector<std::filesystem::path> result;
        for (const auto& entry : std::filesystem::directory_iterator(path)) result.push_back(entry.path());
        std::sort(result.begin(), result.end());
        return result;
    }
    return {path};
}

template <typename State>
void discover_inputs(State& state) {
    state.raw_files.clear(); state.hdf5_input_files.clear(); state.base_names.clear();
    state.raw_base_names.clear(); state.hdf5_base_names.clear();
    std::error_code error;
    for (const auto& entry : input_entries(state, state.options.input)) {
        const std::filesystem::path resolved = std::filesystem::absolute(entry, error).lexically_normal();
        if (!std::filesystem::exists(resolved, error)) {
            throw std::runtime_error(resolved.string() + " does not exist");
        }
        const std::string base = sample_base_name(resolved);
        if (raw_input(resolved)) {
            state.raw_files.push_back(resolved); state.raw_base_names.push_back(base);
        } else if (hdf5_input(resolved)) {
            state.hdf5_input_files.push_back(resolved); state.hdf5_base_names.push_back(base);
        } else {
            state.logger->warning("Skipping unsupported input: " + resolved.string());
            continue;
        }
        state.base_names.push_back(base);
    }
    if (state.raw_files.empty() && state.hdf5_input_files.empty()) {
        throw std::runtime_error("No raw/.d/.d.zip or Raxport HDF5 files found in " + state.options.input);
    }
    auto log_files = [&](const char* label, const std::vector<std::filesystem::path>& files) {
        std::ostringstream message;
        message << label << " files (" << files.size() << "):";
        if (files.empty()) message << " none";
        for (const auto& path : files) message << "\n  " << path.string();
        state.logger->info(message.str());
    };
    log_files("RAW", state.raw_files);
    log_files("HDF5", state.hdf5_input_files);
    const auto controls = split(state.options.negative_control, ',');
    if (!controls.empty()) {
        std::ostringstream message;
        message << "Negative control samples (" << controls.size() << "):";
        for (const auto& control : controls) message << "\n  " << control;
        state.logger->info(message.str());
    }
}

template <typename State>
void create_sample_directories(State& state) {
    for (const auto& base : state.base_names) std::filesystem::create_directories(state.output / base);
}

template <typename State>
void reverse_fasta(State& state) {
    state.logger->info("Reversing fasta sequences to " + state.decoy.string());
    if (!std::filesystem::exists(state.fasta)) throw std::runtime_error("Fasta file does not exist: " + state.fasta.string());
    std::filesystem::create_directories(state.decoy.parent_path());
    std::ifstream input(state.fasta);
    std::ofstream output(state.decoy, std::ios::trunc);
    if (!input || !output) throw std::runtime_error("Unable to create decoy FASTA: " + state.decoy.string());
    std::string header, sequence, line;
    auto write_record = [&]() {
        if (header.empty()) return;
        std::string decoy_sequence = sequence;
        if (state.mode == WorkflowMode::RegularFasta && !sequence.empty()) {
            std::reverse(decoy_sequence.begin() + 1, decoy_sequence.end());
        } else {
            std::reverse(decoy_sequence.begin(), decoy_sequence.end());
        }
        output << header << '\n' << decoy_sequence << '\n';
    };
    while (std::getline(input, line)) {
        if (!line.empty() && line.front() == '>') {
            write_record();
            header = ">" + state.decoy_prefix + trim(line.substr(1));
            sequence.clear();
        } else {
            sequence += trim(line);
        }
    }
    write_record();
}

template <typename State>
void run_raxport(State& state, const std::filesystem::path& raw, const std::filesystem::path& directory,
                 const std::filesystem::path& expected, int threads) {
    std::error_code error;
    if (std::filesystem::exists(expected, error)) {
        const auto size = std::filesystem::file_size(expected, error);
        if (!error && size >= 1024ULL * 1024ULL) {
            state.logger->info("HDF5 file already exists, skipping conversion: " + expected.string());
            return;
        }
        state.logger->warning("Removing incomplete HDF5 conversion output: " + expected.string() +
                              " (" + std::to_string(error ? 0 : size) + " bytes)");
        std::filesystem::remove(expected, error);
    }
    std::filesystem::create_directories(directory);
    std::string heap_limit = environment_value("SIPROS_RAXPORT_GC_HEAP_LIMIT");
    if (heap_limit.empty()) heap_limit = environment_value("DOTNET_GCHeapHardLimit");
    if (heap_limit.empty()) heap_limit = std::to_string(128ULL * 1024ULL * 1024ULL);
    Command command = child_command(
        state.tools.raxport,
        {"-f", raw.string(), "-o", directory.string(), "--format", "hdf5",
         "-n", std::to_string(state.options.n_precursor)},
        threads, {{"DOTNET_GCHeapHardLimit", heap_limit}});
    run_command(state, std::move(command));
    if (!std::filesystem::exists(expected)) {
        throw std::runtime_error("Raxport did not create expected HDF5 file: " + expected.string());
    }
}

template <typename State>
void prepare_hdf5(State& state) {
    state.logger->info("Preparing Raxport HDF5 scan inputs");
    state.hdf5_paths.clear();
    for (std::size_t i = 0; i < state.hdf5_input_files.size(); ++i) {
        const std::filesystem::path source = state.hdf5_input_files[i];
        const std::string base = state.hdf5_base_names[i];
        const std::filesystem::path expected = expected_hdf5(state, base);
        std::filesystem::create_directories(expected.parent_path());
        std::error_code error;
        const bool same = std::filesystem::equivalent(source, expected, error);
        if (!same && !path_entry_exists(expected)) {
            if (!state.stage_hdf5_copies) {
                std::filesystem::create_symlink(std::filesystem::absolute(source), expected, error);
                if (error) {
                    error.clear(); std::filesystem::create_hard_link(source, expected, error);
                }
                if (error) {
                    error.clear(); std::filesystem::copy_file(source, expected, std::filesystem::copy_options::overwrite_existing, error);
                }
                if (error) throw std::runtime_error("Unable to stage HDF5 input: " + error.message());
                state.logger->info("Reusing HDF5 input through lightweight link: " + expected.string());
            } else {
                state.logger->info("Staging HDF5 input " + source.string() + " to " + expected.string());
                std::filesystem::copy_file(source, expected, std::filesystem::copy_options::overwrite_existing, error);
                if (error) throw std::runtime_error("Unable to copy HDF5 input: " + error.message());
            }
        }
        state.hdf5_paths[base] = expected;
    }
    struct Job { std::filesystem::path raw, directory, expected; std::string base; };
    std::vector<Job> jobs;
    for (std::size_t i = 0; i < state.raw_files.size(); ++i) {
        const std::string& base = state.raw_base_names[i];
        const std::filesystem::path expected = expected_hdf5(state, base);
        if (complete_hdf5(expected)) {
            state.logger->info("HDF5 file already exists, skipping conversion: " + expected.string());
            state.hdf5_paths[base] = expected;
        } else {
            jobs.push_back({state.raw_files[i], expected.parent_path(), expected, base});
        }
    }
    const ThreadAllocation allocation = allocate_threads(state.options.threads, static_cast<int>(jobs.size()));
    log_allocation(state, "Raxport conversion", allocation);
    for (std::size_t offset = 0; offset < jobs.size(); offset += static_cast<std::size_t>(std::max(1, allocation.worker_count))) {
        const std::size_t count = std::min<std::size_t>(allocation.worker_count, jobs.size() - offset);
        std::vector<std::future<void>> futures;
        for (std::size_t index = 0; index < count; ++index) {
            const Job job = jobs[offset + index];
            const int threads = allocation.task_threads[offset + index];
            futures.push_back(std::async(std::launch::async, [&state, job, threads] {
                run_raxport(state, job.raw, job.directory, job.expected, threads);
            }));
        }
        std::exception_ptr failure;
        for (auto& future : futures) {
            try { future.get(); } catch (...) { if (!failure) failure = std::current_exception(); }
        }
        if (failure) std::rethrow_exception(failure);
    }
    for (const auto& job : jobs) state.hdf5_paths[job.base] = job.expected;
    std::ostringstream paths;
    paths << "HDF5 scan files:";
    for (const auto& [base, path] : state.hdf5_paths) paths << "\n  " << base << ": " << path.string();
    state.logger->info(paths.str());
}

template <typename State>
void regular_fasta_search(State& state) {
    std::filesystem::create_directories(state.output);
    const std::filesystem::path target_cache = state.output / "target.sfi";
    const std::filesystem::path decoy_cache = state.output / "decoy.sfi";
    auto preparation = [&](const std::filesystem::path& fasta, const std::filesystem::path& cache) {
        std::vector<std::string> args{"search-fasta", "-fasta", fasta.string(), "--prepare-only"};
        add_tolerance_arguments(args, state.options);
        add_ptm_arguments(args, state.options);
        args.insert(args.end(), {"--fragment-index-cache", cache.string()});
        return child_command(state.tools.sipros, std::move(args), state.options.threads);
    };
    state.logger->info("Sipros Regular FASTA cache preparation: target then guarded decoy; 1 process x " +
                       std::to_string(state.options.threads) + " threads");
    run_command(state, preparation(state.fasta, target_cache));
    run_command(state, preparation(state.decoy, decoy_cache));

    const std::string phase =
        "Sipros Regular FASTA paired target/decoy cache-H5 search";
    const ThreadAllocation pair_allocation = sipros_search_allocation(
        state, 2, 1, phase);
    log_allocation(state, phase, pair_allocation);
    for (const auto& base : state.base_names) {
        const std::filesystem::path sample_dir = state.output / base;
        std::filesystem::create_directories(sample_dir);
        auto search = [&](const std::filesystem::path& fasta, const std::filesystem::path& cache,
                          const std::string& pin, const std::string& label) {
            std::vector<std::string> args{
                "search-fasta", "-fasta", fasta.string(), "-f", state.hdf5_paths.at(base).string(),
                "-o", sample_dir.string(), "--pin-output", pin, "--pin-label", label,
                "--top-psms-per-scan", std::to_string(state.options.top_psms_per_scan)};
            add_tolerance_arguments(args, state.options);
            add_ptm_arguments(args, state.options);
            args.insert(args.end(), {"--fragment-index-cache", cache.string()});
            return child_command(state.tools.sipros, std::move(args), 1);
        };
        std::vector<Command> pair;
        pair.push_back(search(state.fasta, target_cache, base + "_target.pin", "1"));
        pair.push_back(search(state.decoy, decoy_cache, base + "_decoy.pin", "-1"));
        run_parallel(state, std::move(pair), pair_allocation);
    }
}

template <typename State>
void sip_fasta_search(State& state) {
    const std::string sip_range = state.options.sip_range.value_or("0-100");
    const std::string precision = state.options.precision.value_or("1");
    std::vector<Command> commands;
    for (const auto& base : state.base_names) {
        const std::filesystem::path sample_dir = state.output / base;
        auto search = [&](const std::filesystem::path& fasta, const std::string& pin, const std::string& label) {
            std::vector<std::string> args{
                "search-fasta", "-fasta", fasta.string(), "-f", state.hdf5_paths.at(base).string(),
                "-o", sample_dir.string(), "--pin-output", pin,
                "-a", state.options.element, "-b", sip_range, "-s", precision,
                "--pin-label", label, "--top-psms-per-scan",
                std::to_string(state.options.top_psms_per_scan)};
            add_tolerance_arguments(args, state.options);
            add_ptm_arguments(args, state.options);
            return child_command(state.tools.sipros, std::move(args), 1);
        };
        commands.push_back(search(state.fasta, base + "_target.pin", "1"));
        commands.push_back(search(state.decoy, base + "_decoy.pin", "-1"));
    }
    const ThreadAllocation allocation = sipros_search_allocation(
        state, static_cast<int>(commands.size()), kMinimumSiprosThreads,
        "Sipros FASTA search");
    log_allocation(state, "Sipros FASTA search", allocation);
    run_parallel(state, std::move(commands), allocation);
}

void require_sfi_pair(const std::filesystem::path& directory) {
    int targets = 0, decoys = 0;
    std::error_code error;
    for (const auto& entry : std::filesystem::directory_iterator(directory, error)) {
        if (!entry.is_regular_file() || lower(entry.path().extension().string()) != ".sfi") continue;
        if (lower(entry.path().filename().string()).find("decoy") == std::string::npos) ++targets;
        else ++decoys;
    }
    if (error || targets != 1 || decoys != 1) {
        throw std::runtime_error("SFI spectra library must contain exactly one target and one generated-decoy index; found " +
                                 std::to_string(targets) + " target and " + std::to_string(decoys) +
                                 " decoy files in " + directory.string());
    }
}

template <typename State>
std::filesystem::path resolve_unlabeled_hdf5(State& state) {
    if (!state.unlabeled_input || state.unlabeled_input->empty()) {
        throw std::runtime_error("--unlabeled-input is required when --spectra-dir is not provided");
    }
    std::error_code error;
    const std::filesystem::path path = std::filesystem::absolute(*state.unlabeled_input, error).lexically_normal();
    if (std::filesystem::is_directory(path, error)) {
        int count = 0;
        for (const auto& entry : std::filesystem::recursive_directory_iterator(path, error)) {
            if (entry.is_regular_file() && hdf5_input(entry.path())) ++count;
        }
        if (count == 0) throw std::runtime_error("No HDF5 files found under unlabeled input directory: " + path.string());
        state.logger->info("Using " + std::to_string(count) + " HDF5 files under " + path.string() +
                           " to match filtered regular-search PSMs");
        return path;
    }
    if (!std::filesystem::exists(path)) throw std::runtime_error("Unlabeled input does not exist: " + path.string());
    if (hdf5_input(path)) return path;
    if (!raw_input(path)) throw std::runtime_error("Unsupported unlabeled input; use raw/.d/.d.zip or Raxport HDF5: " + path.string());
    const std::string base = sample_base_name(path);
    const std::filesystem::path directory = state.output / "unlabeled_hdf5";
    const std::filesystem::path expected = directory / (base + ".h5");
    run_raxport(state, path, directory, expected, state.options.threads);
    return expected;
}

template <typename State>
std::filesystem::path generate_or_reuse_spectra(State& state) {
    if (state.spectra_dir) {
        std::error_code error;
        const std::filesystem::path directory = std::filesystem::absolute(*state.spectra_dir, error).lexically_normal();
        if (!std::filesystem::is_directory(directory)) throw std::runtime_error("--spectra-dir does not exist: " + directory.string());
        require_sfi_pair(directory);
        state.logger->info("Reusing SFI spectra library from " + directory.string());
        return directory;
    }
    if (!state.psm_tsv) throw std::runtime_error("--psm-tsv is required to generate a spectra library");
    const std::filesystem::path unlabeled = resolve_unlabeled_hdf5(state);
    std::filesystem::create_directories(state.generated_spectra_dir);
    std::vector<std::string> args{
        "experimental-spectra", "-i", state.psm_tsv->string(), "-f", unlabeled.string(),
        "-o", state.generated_spectra_dir.string(), "-a", state.options.element,
        "-b", state.options.sip_range.value_or("0-100"), "-s", state.options.precision.value_or("1"),
        "--decoy", "-t", std::to_string(state.options.threads),
        "--envelope-top-n", std::to_string(state.options.sfi_envelope_top_n)};
    add_ptm_arguments(args, state.options, false);
    run_command(state, child_command(state.tools.sipros, std::move(args), state.options.threads));
    require_sfi_pair(state.generated_spectra_dir);
    return state.generated_spectra_dir;
}

template <typename State>
void search_spectra_samples(State& state, const std::filesystem::path& spectra) {
    const std::string phase =
        "Sipros spectra paired target/decoy SFI-H5 search";
    const ThreadAllocation allocation = sipros_search_allocation(
        state, 2, 1, phase);
    log_allocation(state, phase, allocation);
    for (const auto& base : state.base_names) {
        const std::filesystem::path sample_dir = state.output / base;
        std::filesystem::create_directories(sample_dir);
        auto search = [&](const std::string& label) {
            std::vector<std::string> args{
                "search-spectra", "-f", state.hdf5_paths.at(base).string(), "--sfi", spectra.string(),
                "-o", sample_dir.string(), "--tolerance-ms1", number(state.options.tolerance_ms1),
                "--tolerance-ms1-unit", "da", "--tolerance-ms2", number(state.options.tolerance_ms2),
                "--tolerance-ms2-unit", "da", "--rt-tolerance", number(state.options.rt_tolerance),
                "--mvh-cascade-top-n", std::to_string(state.options.mvh_cascade_top_n),
                "--top-psms-per-scan", std::to_string(state.options.top_psms_per_scan),
                "--sfi-label", label};
            return child_command(state.tools.sipros, std::move(args), 1);
        };
        std::vector<Command> pair{search("target"), search("decoy")};
        // Each paired spectra command receives its allocated thread quota explicitly.
        for (std::size_t i = 0; i < pair.size(); ++i) {
            pair[i].arguments.push_back("-t");
            pair[i].arguments.push_back(std::to_string(allocation.task_threads[i]));
        }
        run_parallel(state, std::move(pair), allocation);
    }
}

Workflow::Workflow(WorkflowOptions options, ToolPaths tools, Logger& logger,
                   std::atomic_bool& cancelled)
    : options_(std::move(options)), tools_(std::move(tools)), logger_(logger), cancelled_(cancelled) {}

Workflow::SearchState Workflow::make_search(
    WorkflowMode mode, const std::filesystem::path& output, std::optional<std::filesystem::path> psm_tsv,
    std::optional<std::filesystem::path> unlabeled_input, std::optional<std::filesystem::path> spectra_dir,
    bool stage_hdf5_copies) const {
    SearchState state;
    state.options = options_;
    state.mode = mode;
    state.tools = tools_;
    state.logger = &logger_;
    state.cancelled = &cancelled_;
    state.output = output;
    state.fasta = options_.fasta;
    if (mode != WorkflowMode::Spectra) state.decoy = output / "decoy.faa";
    state.psm_tsv = std::move(psm_tsv);
    state.unlabeled_input = std::move(unlabeled_input);
    state.spectra_dir = std::move(spectra_dir);
    state.generated_spectra_dir = output / "spectra";
    state.stage_hdf5_copies = stage_hdf5_copies;
    const bool spectra_mode = state.psm_tsv || state.unlabeled_input || state.spectra_dir;
    state.decoy_prefix = spectra_mode ? "DECOY_" : "Decoy_";
    state.options.element = mode == WorkflowMode::RegularFasta ? "R" : options_.element;
    return state;
}

void Workflow::run_fasta(SearchState& state) {
    reverse_fasta(state);
    state.logger->info("Workflow CPU thread budget: " + std::to_string(state.options.threads) +
                       " threads (" + std::to_string(physical_cpu_count()) + " physical cores available)");
    discover_inputs(state);
    create_sample_directories(state);
    if (state.options.dry_run) return;
    prepare_hdf5(state);
    if (state.mode == WorkflowMode::RegularFasta) regular_fasta_search(state);
    else sip_fasta_search(state);
}

void Workflow::run_spectra(SearchState& state) {
    if (state.options.element == "R") throw std::runtime_error("Spectra search requires a SIP element such as C13");
    state.logger->info("Workflow CPU thread budget: " + std::to_string(state.options.threads) +
                       " threads (" + std::to_string(physical_cpu_count()) + " physical cores available)");
    discover_inputs(state);
    create_sample_directories(state);
    if (state.options.dry_run) return;
    prepare_hdf5(state);
    const std::filesystem::path spectra = generate_or_reuse_spectra(state);
    search_spectra_samples(state, spectra);
}

void Workflow::run_filter(const SearchState& state, const std::filesystem::path& output,
                          bool assemble_proteins, const std::string& sip_isotope,
                          const std::string& negative_control, bool include_spectra,
                          const std::filesystem::path& prediction_cache) {
    if (state.base_names.empty()) {
        logger_.info("Aerith: no jobs");
        return;
    }
    std::vector<std::string> args{
        "--decoy-prefix", state.decoy_prefix,
        "--sample-parallelism", std::to_string(options_.aerith_sample_parallelism)};
    if (!prediction_cache.empty()) {
        args.insert(args.end(), {"--prediction-cache", prediction_cache.string()});
    }
    if (assemble_proteins) {
        args.insert(args.end(), {"--database", state.fasta.string()});
        if (!state.decoy.empty()) {
            args.insert(args.end(), {"--decoy-database", state.decoy.string()});
        }
        args.insert(args.end(), {"--protein-output-dir", output.string()});
    } else {
        args.insert(args.end(), {"--no-protein-assembly", "--filtered-only"});
    }
    if (options_.ignore_pct) args.push_back("--ignore-pct");
    if (std::any_of(options_.fixed_ptms.begin(), options_.fixed_ptms.end(),
                    [](const std::string& value) { return lower(trim(value)) == "none"; })) {
        args.push_back("--no-fixed-cam");
    }
    if (!sip_isotope.empty() && sip_isotope != "R") {
        args.insert(args.end(), {"--sip-isotope", sip_isotope});
    }
    for (const auto& ptm : options_.fixed_ptms) args.insert(args.end(), {"--fixed-ptm", ptm});
    for (const auto& ptm : options_.ptms) args.insert(args.end(), {"--ptm", ptm});
    if (options_.max_ptm_count) {
        args.insert(args.end(), {"--max-ptm-count", std::to_string(*options_.max_ptm_count)});
    }
    args.insert(args.end(), {
        "--product-top-isotopes", std::to_string(options_.product_top_isotopes),
        "--quant-top-isotopes", std::to_string(options_.n_precursor)});
    if (!negative_control.empty()) {
        args.insert(args.end(), {"--negative-control", negative_control,
                                 "--label-threshold", number(options_.label_threshold)});
    }
    for (const auto& base : state.base_names) {
        const std::filesystem::path prefix = output / base / base;
        args.insert(args.end(), {
            "--target-pin", prefix.string() + "_target.pin",
            "--decoy-pin", prefix.string() + "_decoy.pin",
            "--output-prefix", prefix.string()});
        if (include_spectra) {
            const auto found = state.hdf5_paths.find(base);
            const std::filesystem::path spectra = found == state.hdf5_paths.end() ? expected_hdf5(state, base) : found->second;
            args.insert(args.end(), {"--spectra", spectra.string()});
        }
    }
    std::string operation = assemble_proteins ? "filtering and protein assembly" : "filtering only";
    if (assemble_proteins && include_spectra) operation = "filtering, chromatographic quantification, and protein assembly";
    logger_.info("Running Aerith cross-sample " + operation + " with " +
                 std::to_string(options_.threads) + " CPU threads");
    if (include_spectra) {
        logger_.info("Aerith DIA-NN device policy: CUDA preferred; automatic CPU fallback when CUDA is unavailable or fails");
    }
    Command command = child_command(tools_.aerith, std::move(args), options_.threads);
    command.environment["MKL_NUM_THREADS"] = std::to_string(options_.threads);
    logger_.info(display_command(command));
    if (options_.dry_run) return;
    run_process(command, logger_, cancelled_);
    if (assemble_proteins) {
        const std::filesystem::path protein = output / "combined_protein.tsv";
        if (!std::filesystem::exists(protein)) throw std::runtime_error("Aerith did not create native protein report: " + protein.string());
        if (!negative_control.empty()) {
            for (const auto* name : {"SIP_filtered_psms.tsv", "SIP_target_psms.tsv", "SIP_decoy_psms.tsv",
                                     "combined_protein_with_SIP_filtered_PSM.tsv"}) {
                if (!std::filesystem::exists(output / name)) {
                    throw std::runtime_error("Aerith did not create native negative-control output: " + (output / name).string());
                }
            }
        }
    }
}

void Workflow::run_fast_sip() {
    const std::filesystem::path root(options_.output);
    const std::filesystem::path regular_output = root / "regular";
    const std::filesystem::path spectra_output = root / "spectra_search";
    const std::filesystem::path prediction_cache = regular_output / "regular_search_predictions";
    std::filesystem::create_directories(regular_output); std::filesystem::create_directories(spectra_output);
    std::vector<std::pair<std::string, double>> timings;
    auto phase = std::chrono::steady_clock::now();

    logger_.info("Fast SIP phase 1/4: regular target/decoy FASTA search");
    SearchState regular = make_search(WorkflowMode::RegularFasta, regular_output,
                                      std::nullopt, std::nullopt, std::nullopt, false);
    run_fasta(regular);
    timings.push_back({"1/4 Regular FASTA search", std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()});

    phase = std::chrono::steady_clock::now();
    logger_.info("Fast SIP phase 2/4: Aerith filtering of regular PSMs");
    if (!options_.dry_run) {
        for (const auto* suffix : {".bin", ".spectrum", ".rt"}) {
            const std::filesystem::path cache_file = prediction_cache.string() + suffix;
            std::error_code error;
            const bool removed = std::filesystem::remove(cache_file, error);
            if (error) {
                throw std::runtime_error("Unable to remove stale prediction cache " +
                                         cache_file.string() + ": " + error.message());
            }
            if (removed) logger_.info("Removed stale prediction cache: " + cache_file.string());
        }
    }
    logger_.info("Aerith regular filtering predicts DIA-NN spectra and RT once and caches unique target forms");
    run_filter(regular, regular_output, true, "", "", true, prediction_cache);
    if (!options_.dry_run) {
        for (const auto& base : regular.base_names) {
            const std::filesystem::path path = regular_output / base / (base + "_filtered_psms.tsv");
            std::ifstream input(path);
            std::string header;
            int label_column = -1, accepted = 0;
            if (std::getline(input, header)) {
                const auto columns = split_preserving_empty(header, '\t');
                for (std::size_t i = 0; i < columns.size(); ++i) if (columns[i] == "Label") label_column = static_cast<int>(i);
                std::string row;
                while (std::getline(input, row)) {
                    const auto values = split_preserving_empty(row, '\t');
                    if (label_column < 0 || label_column >= static_cast<int>(values.size()) || trim(values[label_column]) == "1") ++accepted;
                }
            }
            if (accepted == 0) logger_.warning("Aerith regular-search filtering retained 0 target PSMs for " + base);
            else logger_.info("Aerith regular-search filtering retained " + std::to_string(accepted) + " target PSMs for " + base);
        }
    }
    timings.push_back({"2/4 Regular Aerith filter", std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()});

    phase = std::chrono::steady_clock::now();
    logger_.info("Fast SIP phase 3/4: filtered-PSM SFI generation and spectra search");
    SearchState spectra = make_search(WorkflowMode::Spectra, spectra_output,
                                      regular_output, regular_output, std::nullopt, false);
    spectra.base_names = regular.base_names;
    spectra.hdf5_paths = regular.hdf5_paths;
    spectra.generated_spectra_dir = spectra_output;
    if (!options_.dry_run) {
        const std::filesystem::path generated = generate_or_reuse_spectra(spectra);
        search_spectra_samples(spectra, generated);
    }
    timings.push_back({"3/4 SFI + spectra search", std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()});

    phase = std::chrono::steady_clock::now();
    logger_.info("Fast SIP phase 4/4: Aerith filtering and reporting");
    logger_.info("Aerith spectra-search filtering reuses cached target predictions");
    run_filter(spectra, spectra_output, true, options_.element,
               options_.negative_control, true, prediction_cache);
    timings.push_back({"4/4 SIP Aerith + reports", std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()});

    double total = 0.0;
    std::ostringstream report;
    report << "Fast SIP phase timing (wall clock)";
    for (const auto& [label, seconds] : timings) {
        total += seconds;
        report << "\n  " << std::left << std::setw(28) << label << " : " << std::right
               << std::fixed << std::setprecision(3) << seconds << " s";
    }
    report << "\n  " << std::left << std::setw(28) << "Fast SIP total" << " : " << std::right << total << " s";
    logger_.info(report.str());
}

void Workflow::run() {
    const auto started = std::chrono::steady_clock::now();
    std::filesystem::create_directories(options_.output);
    logger_.info("Workflow started: " + mode_name(options_.mode));
    if (options_.mode == WorkflowMode::FastSip) {
        run_fast_sip();
    } else if (options_.mode == WorkflowMode::Spectra) {
        SearchState state = make_search(
            WorkflowMode::Spectra, options_.output,
            options_.psm_tsv ? std::optional<std::filesystem::path>(*options_.psm_tsv) : std::nullopt,
            options_.unlabeled_input ? std::optional<std::filesystem::path>(*options_.unlabeled_input) : std::nullopt,
            options_.spectra_dir ? std::optional<std::filesystem::path>(*options_.spectra_dir) : std::nullopt);
        run_spectra(state);
        run_filter(state, options_.output, false, options_.element,
                   options_.negative_control, false);
    } else {
        SearchState state = make_search(options_.mode, options_.output);
        run_fasta(state);
        run_filter(state, options_.output, true, options_.element,
                   options_.negative_control, true);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    logger_.info("Workflow completed successfully");
    logger_.info("Results directory: " + options_.output);
    logger_.info("Total wall time: " + number(seconds) + " s");
}

} // namespace siproswf
