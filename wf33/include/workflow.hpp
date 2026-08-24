#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <fstream>
#include <map>
#include <mutex>
#include <optional>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

namespace siproswf {

inline constexpr const char kWorkflowLogFilename[] = "sipros_workflow.log";

constexpr int kMinimumSiprosThreads = 8;
inline constexpr std::uint64_t kParallelSiprosMinimumUsableMemoryBytes =
    std::uint64_t{32} * 1024 * 1024 * 1024;

enum class WorkflowMode {
    Auto,
    RegularFasta,
    SipFasta,
    FastSip,
    Spectra,
};

struct WorkflowOptions {
    WorkflowMode mode = WorkflowMode::Auto;
    std::string input;
    std::string fasta;
    std::string output;
    std::vector<std::string> ptms;
    std::vector<std::string> fixed_ptms;
    std::optional<int> max_ptm_count;
    double tolerance_ms1 = 0.01;
    double tolerance_ms2 = 0.01;
    std::string element;
    std::optional<std::string> sip_range;
    std::optional<std::string> precision;
    std::optional<std::string> psm_tsv;
    std::optional<std::string> unlabeled_input;
    std::optional<std::string> spectra_dir;
    int n_precursor = 6;
    int product_top_isotopes = 5;
    int threads = 0;
    int aerith_sample_parallelism = 3;
    int top_psms_per_scan = 20;
    double rt_tolerance = 5.0;
    int sfi_envelope_top_n = 3;
    int mvh_cascade_top_n = 150;
    bool ignore_pct = false;
    std::string negative_control;
    double label_threshold = 2.0;
    bool dry_run = false;
};

struct ParseResult {
    WorkflowOptions options;
    bool show_help = false;
    bool launch_gui = false;
    std::vector<std::string> warnings;
};

struct ThreadAllocation {
    int worker_count = 0;
    std::vector<int> task_threads;

    [[nodiscard]] int peak_threads() const;
};

struct ToolPaths {
    std::filesystem::path raxport;
    std::filesystem::path sipros;
    std::filesystem::path aerith;
};

struct Command {
    std::filesystem::path executable;
    std::vector<std::string> arguments;
    std::map<std::string, std::string> environment;
    std::optional<std::filesystem::path> working_directory;
    int cpu_threads = 1;
};

class Logger {
public:
    using Sink = std::function<void(const std::string&)>;

    Logger(const std::filesystem::path& log_path, Sink sink = {});
    ~Logger();
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);
    void debug(const std::string& message);

private:
    void write(const char* level, const std::string& message);
    std::mutex mutex_;
    std::ofstream stream_;
    Sink sink_;
};

int physical_cpu_count();
std::uint64_t usable_memory_bytes();
int effective_thread_count(int requested, int available = 0);
ThreadAllocation allocate_threads(int total_threads, int task_count,
                                  int minimum_threads_per_task = 1);
bool serialize_sipros_searches(int total_threads,
                               std::uint64_t usable_memory);
ThreadAllocation allocate_sipros_search_threads(
    int total_threads, int task_count, int minimum_threads_per_task,
    std::uint64_t usable_memory);
std::map<std::string, std::string> thread_environment(int thread_count);

ParseResult parse_arguments(const std::vector<std::string>& arguments);
void validate_options(WorkflowOptions& options, std::vector<std::string>* warnings = nullptr);
void print_help(std::ostream& output);
std::string mode_name(WorkflowMode mode);

std::string quote_display_argument(const std::string& value);
std::string display_command(const Command& command);
int run_process(const Command& command, Logger& logger, std::atomic_bool& cancelled);

ToolPaths locate_tools(const std::filesystem::path& executable_path);

class Workflow {
public:
    Workflow(WorkflowOptions options, ToolPaths tools, Logger& logger,
             std::atomic_bool& cancelled);
    void run();

private:
    struct SearchState;

    WorkflowOptions options_;
    ToolPaths tools_;
    Logger& logger_;
    std::atomic_bool& cancelled_;

    SearchState make_search(WorkflowMode mode, const std::filesystem::path& output,
                            std::optional<std::filesystem::path> psm_tsv = std::nullopt,
                            std::optional<std::filesystem::path> unlabeled_input = std::nullopt,
                            std::optional<std::filesystem::path> spectra_dir = std::nullopt,
                            bool stage_hdf5_copies = true) const;
    void run_fasta(SearchState& state);
    void run_spectra(SearchState& state);
    void run_fast_sip();
    void run_filter(const SearchState& state, const std::filesystem::path& output,
                    bool assemble_proteins, const std::string& sip_isotope,
                    const std::string& negative_control = {},
                    bool include_spectra = false,
                    const std::filesystem::path& prediction_cache = {});
};

int run_gui(const std::filesystem::path& executable_path);

} // namespace siproswf
