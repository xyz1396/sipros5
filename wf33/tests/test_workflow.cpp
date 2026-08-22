#include "workflow.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <future>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace siproswf {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void require_throws(Function function, const std::string& message) {
    try { function(); }
    catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}

void test_thread_allocation() {
    const auto exact = siproswf::allocate_threads(9, 2);
    require(exact.worker_count == 2, "nine threads should run a pair concurrently");
    require(exact.task_threads == std::vector<int>({5, 4}), "odd budget should favor first task");
    require(exact.peak_threads() == 9, "peak must equal total budget");

    const auto minimum = siproswf::allocate_threads(
        15, 6, siproswf::kMinimumSiprosThreads);
    require(minimum.worker_count == 1, "15 threads cannot host two eight-thread jobs");
    require(minimum.peak_threads() <= 15, "allocation may not oversubscribe");
    require_throws([] { siproswf::allocate_threads(0, 1); },
                   "zero total must be rejected");
}

void test_cli_modes() {
    auto regular = siproswf::parse_arguments({
        "--regular-fasta-search", "-i", "input.h5", "-f", "db.faa", "-o", "out"});
    require(regular.options.mode == siproswf::WorkflowMode::RegularFasta,
            "regular switch");
    require(regular.options.element == "R", "regular chemistry");

    auto sip = siproswf::parse_arguments({
        "--sip-fasta-search", "-i", "input.h5", "-f", "db.faa", "-o", "out"});
    require(sip.options.mode == siproswf::WorkflowMode::SipFasta, "SIP switch");
    require(sip.options.element == "C13", "SIP switch defaults to C13");

    auto sip_control = siproswf::parse_arguments({
        "--sip-fasta-search", "-i", "input.h5", "-f", "db.faa", "-o", "out",
        "--negative_control", "control_a,control_b"});
    require(sip_control.options.negative_control == "control_a,control_b",
            "SIP negative-control samples");

    auto fast = siproswf::parse_arguments({
        "--fast-sip-search", "-i", "input.h5", "-f", "db.faa", "-o", "out"});
    require(fast.options.mode == siproswf::WorkflowMode::FastSip, "fast SIP switch");

    auto legacy = siproswf::parse_arguments({
        "-i", "input.h5", "-f", "db.faa", "-e", "n15", "-o", "out"});
    require(legacy.options.mode == siproswf::WorkflowMode::SipFasta,
            "legacy isotope inference");
    require(legacy.options.element == "N15", "isotope normalization");

    auto spectra = siproswf::parse_arguments({
        "-i", "input.h5", "-f", "db.faa", "--spectra-dir", "sfi", "-o", "out"});
    require(spectra.options.mode == siproswf::WorkflowMode::Spectra,
            "spectra inference");

    require_throws([] {
        siproswf::parse_arguments({"--regular-fasta-search", "--sip-fasta-search",
                               "-i", "x", "-f", "x", "-o", "x"});
    }, "multiple modes must fail");
    require_throws([] {
        siproswf::parse_arguments({"--regular-fasta-search", "-e", "C13",
                               "-i", "x", "-f", "x", "-o", "x"});
    }, "regular/SIP contradiction must fail");
    require_throws([] {
        siproswf::parse_arguments({"--regular-fasta-search", "-i", "x", "-f", "x",
                               "-o", "x", "--negative_control", "control"});
    }, "regular search must reject negative-control samples");
    require_throws([] {
        siproswf::parse_arguments({"--gui", "-i", "x"});
    }, "GUI and headless arguments must not be mixed");
}

void test_environment_limits() {
    const auto environment = siproswf::thread_environment(7);
    require(environment.at("OMP_NUM_THREADS") == "7", "OpenMP quota");
    require(environment.at("DOTNET_PROCESSOR_COUNT") == "7", ".NET quota");
    require(environment.at("OPENBLAS_NUM_THREADS") == "1", "nested BLAS quota");
    require(siproswf::effective_thread_count(999999) <=
                siproswf::physical_cpu_count(),
            "explicit thread count must be capped to physical cores");
}

void test_persistent_workflow_log() {
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("siproswf-log-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    const std::filesystem::path log_path = root / "output" / kWorkflowLogFilename;
    try {
        {
            Logger logger(log_path);
            logger.info("\n\npersistent-info\r\n\r\n  \t\r\npersistent-detail\n");
            logger.error("persistent-error");
        }
        std::ifstream input(log_path);
        const std::string contents((std::istreambuf_iterator<char>(input)), {});
        require(input.good() || input.eof(), "workflow log must be readable");
        require(contents.find("SIPROS WORKFLOW LOG") != std::string::npos,
                "workflow log must contain a readable session header");
        require(contents.find("| INFO  | persistent-info") !=
                    std::string::npos,
                "workflow INFO message must be written to the output log");
        require(contents.find("| INFO  | persistent-detail") !=
                    std::string::npos,
                "multiline workflow messages must retain meaningful lines");
        require(contents.find("| ERROR | persistent-error") !=
                    std::string::npos,
                "workflow terminal error must be written to the output log");
        require(contents.find("\n\n") == std::string::npos,
                "workflow log must not contain redundant empty lines");
        require(contents.find("Session closed : ") != std::string::npos,
                "workflow log must contain a clean session footer");
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove_all(root, ignored);
        throw;
    }
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
}

void test_packaged_tool_layout() {
    if (std::getenv("SIPROSWF_RAXPORT") != nullptr ||
        std::getenv("SIPROSWF_SIPROS") != nullptr ||
        std::getenv("SIPROSWF_AERITH") != nullptr) {
        return;
    }
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("siproswf-tools-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    const std::filesystem::path library = root / "lib";
#ifdef _WIN32
    const std::filesystem::path workflow = root / "siproswf.exe";
    const std::filesystem::path raxport = library / "Raxport-win-x64.exe";
    const std::filesystem::path sipros = library / "sipros.exe";
    const std::filesystem::path aerith = library / "aerith.exe";
#else
    const std::filesystem::path workflow = root / "siproswf";
    const std::filesystem::path raxport = library / "Raxport-linux-x64";
    const std::filesystem::path sipros = library / "sipros";
    const std::filesystem::path aerith = library / "aerith";
#endif
    std::filesystem::create_directories(library);
    try {
        for (const auto& tool : {raxport, sipros, aerith}) {
            std::ofstream output(tool);
            output << "test";
        }
        const siproswf::ToolPaths tools = siproswf::locate_tools(workflow);
        require(tools.raxport == raxport, "packaged Raxport must resolve from lib");
        require(tools.sipros == sipros, "packaged Sipros must resolve from lib");
        require(tools.aerith == aerith, "packaged Aerith must resolve from lib");
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove_all(root, ignored);
        throw;
    }
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
}

void test_bin_tool_layout() {
    if (std::getenv("SIPROSWF_RAXPORT") != nullptr ||
        std::getenv("SIPROSWF_SIPROS") != nullptr ||
        std::getenv("SIPROSWF_AERITH") != nullptr) {
        return;
    }
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("siproswf-bin-tools-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
#ifdef _WIN32
    const std::filesystem::path workflow = root / "siproswf.exe";
    const std::filesystem::path raxport = root / "Raxport-win-x64.exe";
    const std::filesystem::path sipros = root / "sipros.exe";
    const std::filesystem::path aerith = root / "aerith.exe";
#else
    const std::filesystem::path workflow = root / "siproswf";
    const std::filesystem::path raxport = root / "Raxport-linux-x64";
    const std::filesystem::path sipros = root / "sipros";
    const std::filesystem::path aerith = root / "aerith";
#endif
    std::filesystem::create_directories(root);
    try {
        for (const auto& tool : {raxport, sipros, aerith}) {
            std::ofstream output(tool);
            output << "test";
        }
        const siproswf::ToolPaths tools = siproswf::locate_tools(workflow);
        require(tools.raxport == raxport, "bin Raxport must resolve beside siproswf");
        require(tools.sipros == sipros, "bin Sipros must resolve beside siproswf");
        require(tools.aerith == aerith, "bin Aerith must resolve beside siproswf");
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove_all(root, ignored);
        throw;
    }
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
}

void test_decoy_rules() {
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("siproswf-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    try {
        const std::filesystem::path fasta = root / "target.faa";
        const std::filesystem::path hdf5 = root / "sample.h5";
        { std::ofstream output(fasta); output << ">protein\nABCDE\n"; }
        { std::ofstream output(hdf5); output << "test"; }
        siproswf::ToolPaths tools{root / "raxport", root / "sipros", root / "aerith"};
        std::atomic_bool cancelled{false};

        siproswf::WorkflowOptions regular;
        regular.mode = siproswf::WorkflowMode::RegularFasta;
        regular.input = hdf5.string(); regular.fasta = fasta.string();
        regular.output = (root / "regular").string(); regular.element = "R";
        regular.threads = 1; regular.dry_run = true;
        std::filesystem::create_directories(regular.output);
        siproswf::Logger regular_log(std::filesystem::path(regular.output) / "test.log");
        siproswf::Workflow regular_workflow(regular, tools, regular_log, cancelled);
        regular_workflow.run();
        std::ifstream regular_decoy(std::filesystem::path(regular.output) / "decoy.faa");
        const std::string regular_text((std::istreambuf_iterator<char>(regular_decoy)), {});
        require(regular_text.find(">Decoy_protein\nAEDCB") != std::string::npos,
                "regular decoy must preserve the protein N-terminal residue");

        siproswf::WorkflowOptions sip = regular;
        sip.mode = siproswf::WorkflowMode::SipFasta; sip.element = "C13";
        sip.output = (root / "sip").string(); std::filesystem::create_directories(sip.output);
        siproswf::Logger sip_log(std::filesystem::path(sip.output) / "test.log");
        siproswf::Workflow sip_workflow(sip, tools, sip_log, cancelled);
        sip_workflow.run();
        std::ifstream sip_decoy(std::filesystem::path(sip.output) / "decoy.faa");
        const std::string sip_text((std::istreambuf_iterator<char>(sip_decoy)), {});
        require(sip_text.find(">Decoy_protein\nEDCBA") != std::string::npos,
                "SIP decoy must retain full reversal");
    } catch (...) {
        std::error_code ignored; std::filesystem::remove_all(root, ignored); throw;
    }
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
}

void test_process_runner(const std::filesystem::path& fake_tool) {
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("siproswf-process-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    try {
        std::vector<std::string> logs;
        siproswf::Logger logger(root / "process.log", [&](const std::string& line) {
            logs.push_back(line);
        });
        std::atomic_bool cancelled{false};
        siproswf::Command command;
        command.executable = fake_tool;
        command.arguments = {"space value", "literal&shell-token"};
        command.environment = siproswf::thread_environment(2);
        command.environment["SIPROSWF_TEST_VALUE"] = "kept";
        command.cpu_threads = 2;
        require(siproswf::run_process(command, logger, cancelled) == 0,
                "fake child should succeed");
        std::string combined;
        for (const auto& line : logs) combined += line + "\n";
        require(combined.find("[space value]") != std::string::npos,
                "argument containing spaces must remain one argument");
        require(combined.find("[literal&shell-token]") != std::string::npos,
                "shell metacharacters must remain literal");
        require(combined.find("environment:kept") != std::string::npos,
                "child environment override must be visible");
        require(combined.find("fake-tool-stderr") != std::string::npos,
                "stderr must stream into the workflow log");

        logs.clear();
        command.arguments = {"--split-line"};
        require(siproswf::run_process(command, logger, cancelled) == 0,
                "fake child with split output should succeed");
        combined.clear();
        for (const auto& line : logs) combined += line + "\n";
        require(combined.find(std::string(88, '-')) != std::string::npos,
                "one child output line split across pipe reads must remain one log line");
        const auto contains_payload = [&](const std::string& payload) {
            const std::string marker = "| INFO  | ";
            return std::any_of(logs.begin(), logs.end(), [&](const std::string& line) {
                const auto found = line.find(marker);
                return found != std::string::npos &&
                       line.substr(found + marker.size()) == payload;
            });
        };
        require(contains_payload(std::string(88, '-')),
                "the reassembled child output must be one complete log payload");
        require(!contains_payload(std::string(66, '-')) &&
                    !contains_payload(std::string(22, '-')),
                "a partial pipe-read fragment must not be logged as a complete line");

        command.arguments = {"--exit", "7"};
        require_throws([&] { siproswf::run_process(command, logger, cancelled); },
                       "nonzero child status must fail the workflow");

        command.arguments = {"--sleep", "5000"};
        auto running = std::async(std::launch::async, [&] {
            require_throws([&] { siproswf::run_process(command, logger, cancelled); },
                           "cancelled child must fail the workflow");
        });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        cancelled.store(true);
        require(running.wait_for(std::chrono::seconds(3)) == std::future_status::ready,
                "cancellation must terminate the process tree promptly");
        running.get();
    } catch (...) {
        std::error_code ignored; std::filesystem::remove_all(root, ignored); throw;
    }
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
}

} // namespace siproswf

int main(int argc, char** argv) {
    try {
        siproswf::test_thread_allocation();
        siproswf::test_cli_modes();
        siproswf::test_environment_limits();
        siproswf::test_persistent_workflow_log();
        siproswf::test_packaged_tool_layout();
        siproswf::test_bin_tool_layout();
        siproswf::test_decoy_rules();
        siproswf::require(argc == 2, "fake-tool path argument is required");
        siproswf::test_process_runner(argv[1]);
        std::cout << "siproswf tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "siproswf test failure: " << error.what() << '\n';
        return 1;
    }
}
