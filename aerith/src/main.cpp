#include "filter.hpp"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

class CommandLine final {
public:
    static bool parse(
        int argc, char** argv, aerith::Config& config, int& exit_status);

private:
    static void usage(std::ostream& out);
    static double number(const std::string& value, const char* option);
    static unsigned int unsigned_number(
        const std::string& value, const char* option);
    static void validate(const aerith::Config& config, int& exit_status);
};

void CommandLine::usage(std::ostream& out) {
    out << "Aerith: DIA-NN features and Percolator-style SVM filtering for Sipros PIN files\n\n"
        << "Usage: aerith --target-pin TARGET.pin --decoy-pin DECOY.pin\n"
        << "              --output-prefix PATH [repeat the three options per sample] [options]\n"
        << "       aerith --input MERGED.pin --output-prefix PATH [options]\n\n"
        << "Options:\n"
        << "  --max-iterations INT     Semi-supervised SVM iterations (default 10)\n"
        << "  --c-pos FLOAT            Positive SVM class cost (default 1)\n"
        << "  --c-neg FLOAT            Negative SVM class cost (default 1)\n"
        << "  --initial-score NAME     Score used to seed RT training (default WDPscores)\n"
        << "  --q-threshold FLOAT      Reporting threshold (default 0.01)\n"
        << "  --train-fdr FLOAT        RT/SVM selection threshold (default 0.01)\n"
        << "  --rt-ridge FLOAT         RT OLS regularization (default 1e-4)\n"
        << "  --database FILE          Target FASTA used for protein assembly\n"
        << "  --decoy-database FILE    Decoy FASTA used for protein assembly\n"
        << "  --protein-output-dir DIR Write combined protein.tsv/protein.fas here\n"
        << "  --no-protein-assembly    Skip native picked-FDR/razor protein assembly\n"
        << "  --filtered-only          Write only PREFIX_filtered_psms.tsv\n"
        << "  --no-fixed-cam           Do not report fixed C carbamidomethylation\n"
        << "  --spectra FILE           Raxport HDF5 MS2 file; repeat once per sample\n"
        << "  --spectrum-model FILE    DIA-NN TorchScript model (default: beside aerith)\n"
        << "  --rt-model FILE          DIA-NN RT TorchScript model (default: beside aerith)\n"
        << "  --no-predicted-rt        Use legacy Aerith RT model instead of DIA-NN RT\n"
        << "  --fragment-ppm FLOAT     Fragment matching tolerance (default 20)\n"
        << "  --decoy-prefix TEXT      Target-decoy FASTA prefix (default Decoy_)\n"
        << "  --ignore-pct             Exclude SIP abundance columns from the SVM\n"
        << "  -h, --help               Show this help\n\n"
        << "Outputs per sample are PREFIX_target_psms.tsv, PREFIX_decoy_psms.tsv,\n"
        << "and PREFIX_filtered_psms.tsv. With target and decoy databases, Aerith writes\n"
        << "native protein.tsv/protein.fas reports without intermediate pepXML.\n";
}

double CommandLine::number(const std::string& value, const char* option) {
    std::size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size()) {
        throw std::runtime_error(std::string("Invalid value for ") + option + ": " + value);
    }
    return parsed;
}

unsigned int CommandLine::unsigned_number(
    const std::string& value, const char* option) {
    std::size_t used = 0;
    const auto parsed = std::stoul(value, &used);
    if (used != value.size() || parsed > 1000) {
        throw std::runtime_error(std::string("Invalid value for ") + option + ": " + value);
    }
    return static_cast<unsigned int>(parsed);
}

void CommandLine::validate(const aerith::Config& config, int& exit_status) {
    const bool paired = !config.target_pins.empty() || !config.decoy_pins.empty();
    if ((config.inputs.empty() && !paired) || config.output_prefixes.empty()) {
        usage(std::cerr);
        exit_status = EXIT_FAILURE;
        return;
    }
    if (paired && !config.inputs.empty()) {
        throw std::runtime_error("Use either --input or target/decoy PIN pairs, not both");
    }
    if (paired && (config.target_pins.size() != config.decoy_pins.size() ||
                   config.target_pins.size() != config.output_prefixes.size())) {
        throw std::runtime_error(
            "Repeat --target-pin, --decoy-pin, and --output-prefix equally");
    }
    const std::size_t samples = paired ? config.target_pins.size() : config.inputs.size();
    if (!config.spectrum_paths.empty() && config.spectrum_paths.size() != samples) {
        throw std::runtime_error("Repeat --spectra exactly once per input sample");
    }
    if (!paired && config.output_prefixes.size() != 1 &&
        config.output_prefixes.size() != config.inputs.size()) {
        throw std::runtime_error(
            "Use one aggregate --output-prefix or repeat it once per --input");
    }
    if (!(config.q_threshold > 0.0 && config.q_threshold <= 1.0) ||
        !(config.train_fdr > 0.0 && config.train_fdr <= 1.0) ||
        config.rt_ridge < 0.0 || config.svm_c_pos <= 0.0 ||
        config.svm_c_neg <= 0.0 || config.fragment_ppm <= 0.0) {
        throw std::runtime_error(
            "FDR thresholds must be in (0,1], SVM costs and fragment ppm positive, "
            "and RT ridge non-negative");
    }
    if (!config.protein_output_dir.empty() &&
        (config.database_path.empty() || config.decoy_database_path.empty())) {
        throw std::runtime_error(
            "--protein-output-dir requires --database and --decoy-database");
    }
    if (config.assemble_proteins &&
        (config.database_path.empty() != config.decoy_database_path.empty())) {
        throw std::runtime_error(
            "Protein assembly requires both --database and --decoy-database");
    }
    if (config.assemble_proteins) {
        for (const auto& database : {
                 config.database_path, config.decoy_database_path}) {
            if (!database.empty() && !std::filesystem::is_regular_file(database)) {
                throw std::runtime_error(
                    "Protein database does not exist: " + database);
            }
        }
    }
    exit_status = EXIT_SUCCESS;
}

bool CommandLine::parse(
    int argc, char** argv, aerith::Config& config, int& exit_status) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&](const char* option) -> std::string {
            if (++i >= argc) {
                throw std::runtime_error(std::string("Missing value for ") + option);
            }
            return argv[i];
        };
        if (arg == "-h" || arg == "--help") {
            usage(std::cout);
            exit_status = EXIT_SUCCESS;
            return false;
        } else if (arg == "--input") {
            config.inputs.push_back(value("--input"));
        } else if (arg == "--target-pin") {
            config.target_pins.push_back(value("--target-pin"));
        } else if (arg == "--decoy-pin") {
            config.decoy_pins.push_back(value("--decoy-pin"));
        } else if (arg == "--output-prefix") {
            config.output_prefixes.push_back(value("--output-prefix"));
        } else if (arg == "--database") {
            config.database_path = value("--database");
        } else if (arg == "--decoy-database") {
            config.decoy_database_path = value("--decoy-database");
        } else if (arg == "--protein-output-dir") {
            config.protein_output_dir = value("--protein-output-dir");
        } else if (arg == "--no-protein-assembly") {
            config.assemble_proteins = false;
        } else if (arg == "--filtered-only") {
            config.filtered_only = true;
        } else if (arg == "--no-fixed-cam") {
            config.fixed_cam = false;
        } else if (arg == "--spectra") {
            config.spectrum_paths.push_back(value("--spectra"));
        } else if (arg == "--spectrum-model") {
            config.spectrum_model_path = value("--spectrum-model");
        } else if (arg == "--rt-model") {
            config.rt_model_path = value("--rt-model");
        } else if (arg == "--no-predicted-rt") {
            config.predict_rt = false;
        } else if (arg == "--fragment-ppm") {
            config.fragment_ppm = number(value("--fragment-ppm"), "--fragment-ppm");
        } else if (arg == "--decoy-prefix") {
            config.decoy_prefix = value("--decoy-prefix");
        } else if (arg == "--ignore-pct") {
            config.ignore_pct = true;
        } else if (arg == "--initial-score") {
            config.initial_score = value("--initial-score");
        } else if (arg == "--q-threshold") {
            config.q_threshold = number(value("--q-threshold"), "--q-threshold");
        } else if (arg == "--train-fdr") {
            config.train_fdr = number(value("--train-fdr"), "--train-fdr");
        } else if (arg == "--rt-ridge" || arg == "--ridge") {
            config.rt_ridge = number(value(arg.c_str()), arg.c_str());
        } else if (arg == "--c-pos") {
            config.svm_c_pos = number(value("--c-pos"), "--c-pos");
        } else if (arg == "--c-neg") {
            config.svm_c_neg = number(value("--c-neg"), "--c-neg");
        } else if (arg == "--max-iterations") {
            config.max_iterations = unsigned_number(
                value("--max-iterations"), "--max-iterations");
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }
    exit_status = EXIT_SUCCESS;
    validate(config, exit_status);
    return exit_status == EXIT_SUCCESS;
}

int main(int argc, char** argv) {
    try {
        aerith::Config config;
        std::error_code executable_error;
        auto executable_path =
            std::filesystem::read_symlink("/proc/self/exe", executable_error);
        if (executable_error) executable_path = std::filesystem::absolute(argv[0]);
        const auto executable = executable_path.parent_path();
        const auto locate_model = [&](const char* name) {
            const auto sibling = executable / name;
            if (std::filesystem::exists(sibling)) return sibling.string();
            return std::string{};
        };
        config.spectrum_model_path = locate_model("diann-2.6.1-fragmentation.pt");
        config.rt_model_path = locate_model("diann-2.6.1-retention-time.pt");
        int exit_status = EXIT_SUCCESS;
        if (!CommandLine::parse(argc, argv, config, exit_status)) {
            return exit_status;
        }
        const auto summary = aerith::run(config);
        std::ostringstream report;
        aerith::print_summary(report, summary);
        const auto report_text = report.str();
        std::cout << report_text;
        if (!summary.protein_assembly_stages.empty() &&
            !summary.protein_output_dir.empty()) {
            const auto log_path =
                std::filesystem::path(summary.protein_output_dir) / "aerith.log";
            std::ofstream log(log_path);
            if (!log) {
                throw std::runtime_error(
                    "Cannot create Aerith log: " + log_path.string());
            }
            log << report_text;
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "aerith: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
