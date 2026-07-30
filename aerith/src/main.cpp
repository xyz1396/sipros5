#include "filter.hpp"

#include <cstdlib>
#include <cmath>
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
        << "  --protein-output-dir DIR Write combined_*.tsv reports here\n"
        << "  --no-protein-assembly    Skip native picked-FDR/razor protein assembly\n"
        << "  --filtered-only          Write only PREFIX_filtered_psms.tsv\n"
        << "  --no-fixed-cam           Do not report fixed C carbamidomethylation\n"
        << "  --spectra FILE           Raxport HDF5 MS2 file; repeat once per sample\n"
        << "  --spectrum-model FILE    DIA-NN TorchScript model (default: beside aerith)\n"
        << "  --rt-model FILE          DIA-NN RT TorchScript model (default: beside aerith)\n"
        << "  --prediction-cache FILE  Reuse/save generated spectrum and RT features\n"
        << "  --no-streaming           Keep all sample PSMs in RAM (diagnostic)\n"
        << "  --sample-parallelism INT Concurrent streamed samples (default 3; raises RAM)\n"
        << "  --no-predicted-rt        Use legacy Aerith RT model instead of DIA-NN RT\n"
        << "  --sip-isotope NAME       Shift predicted/quant ions for C13,H2,N15,O18,S34\n"
        << "  --fixed-ptm NAME         Fixed-PTM chemistry used by the FASTA search; repeatable\n"
        << "  --ptm NAME               Variable-PTM chemistry used by the FASTA search; repeatable\n"
        << "  --max-ptm-count INT      FASTA-search variable-PTM limit\n"
        << "  --fragment-ppm FLOAT     Fragment matching tolerance (default 20)\n"
        << "  --product-top-isotopes INT Top SIP isotopes per predicted product ion (default 5)\n"
        << "  --quant-mz-ppm FLOAT     MS1 XIC tolerance (default 10)\n"
        << "  --quant-rt-window FLOAT  MS1 XIC RT window in minutes (default 0.4)\n"
        << "  --quant-min-isotopes INT Minimum isotope traces (default 2)\n"
        << "  --quant-top-isotopes INT Top theoretical SIP precursor peaks (default 6)\n"
        << "  --quant-min-scans INT    Minimum MS1 scans in a feature (default 3)\n"
        << "  --quant-intensity-mode INT Ion intensity: 0 apex, 1 area, 2 auto (default 2)\n"
        << "  --no-quant-normalization Keep combined intensities on their raw scale\n"
        << "  --no-mbr                 Disable match-between-runs transfer\n"
        << "  --mbr-rt-window FLOAT    MBR local RT alignment window in minutes (default 1)\n"
        << "  --mbr-top-runs INT       Maximum donor runs per acceptor (default 10)\n"
        << "  --mbr-min-correlation FLOAT Minimum overlap-weighted donor correlation (default 0)\n"
        << "  --mbr-ion-fdr FLOAT      Transferred-ion FDR threshold (default 0.01)\n"
        << "  --mbr-sip-bin-width FLOAT SIP abundance-bin width in percent (default 10)\n"
        << "  --decoy-prefix TEXT      Target-decoy FASTA prefix (default Decoy_)\n"
        << "  --protein-reference FILE Existing combined_protein.tsv for SIP PSM mapping\n"
        << "  --sip-protein-output FILE Native SIP protein/PSM mapping output\n"
        << "  --negative-control NAME  Input sample basename used as an unlabeled control; comma-separated\n"
        << "  --label-threshold FLOAT  Minimum MS2 SIP abundance for negative-control filtering (default 2)\n"
        << "  --ignore-pct             Exclude SIP abundance columns from the SVM\n"
        << "  -h, --help               Show this help\n\n"
        << "Outputs per sample are PREFIX_target_psms.tsv, PREFIX_decoy_psms.tsv,\n"
        << "and PREFIX_filtered_psms.tsv. With target and decoy databases, Aerith writes\n"
        << "native sample and combined reports without intermediate pepXML.\n";
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
    if (config.sample_parallelism == 0) {
        throw std::runtime_error(
            "--sample-parallelism must be a positive integer");
    }
    if (!(config.q_threshold > 0.0 && config.q_threshold <= 1.0) ||
        !(config.train_fdr > 0.0 && config.train_fdr <= 1.0) ||
        config.rt_ridge < 0.0 || config.svm_c_pos <= 0.0 ||
        config.svm_c_neg <= 0.0 || config.fragment_ppm <= 0.0 ||
        config.product_top_isotopes == 0 ||
        config.quant_mz_ppm <= 0.0 || config.quant_rt_window <= 0.0 ||
        config.quant_min_isotopes == 0 ||
        config.quant_top_isotopes == 0 ||
        config.quant_min_isotopes > config.quant_top_isotopes ||
        config.quant_min_scans == 0 || config.quant_intensity_mode > 2 ||
        config.mbr_rt_window <= 0.0 ||
        config.mbr_top_runs == 0 || config.mbr_min_correlation < -1.0 ||
        config.mbr_min_correlation > 1.0 || config.mbr_ion_fdr <= 0.0 ||
        config.mbr_ion_fdr > 1.0 ||
        !std::isfinite(config.mbr_sip_bin_width) ||
        config.mbr_sip_bin_width <= 0.0 ||
        config.mbr_sip_bin_width > 100.0 ||
        !std::isfinite(config.label_threshold) ||
        config.label_threshold < 0.0 || config.label_threshold > 100.0) {
        throw std::runtime_error(
            "FDR thresholds must be in (0,1], SVM costs and fragment ppm positive, "
            "quantification tolerances/minima valid, and RT ridge non-negative");
    }
    if (!config.protein_output_dir.empty() &&
        (config.database_path.empty() || config.decoy_database_path.empty())) {
        throw std::runtime_error(
            "--protein-output-dir requires --database and --decoy-database");
    }
    if (config.protein_reference_path.empty() !=
        config.sip_protein_output_path.empty()) {
        throw std::runtime_error(
            "--protein-reference and --sip-protein-output must be provided together");
    }
    if (!config.protein_reference_path.empty() &&
        !std::filesystem::is_regular_file(config.protein_reference_path)) {
        throw std::runtime_error(
            "Protein reference does not exist: " +
            config.protein_reference_path);
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
    if (!config.negative_control_samples.empty() &&
        (!paired || config.protein_output_dir.empty() ||
         !config.assemble_proteins || config.sip_isotope.empty())) {
        throw std::runtime_error(
            "--negative-control requires paired target/decoy inputs and "
            "native SIP protein assembly with --protein-output-dir and "
            "--sip-isotope");
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
        } else if (arg == "--prediction-cache") {
            config.prediction_cache_path = value("--prediction-cache");
        } else if (arg == "--no-streaming") {
            config.stream_samples = false;
        } else if (arg == "--sample-parallelism") {
            config.sample_parallelism = unsigned_number(
                value("--sample-parallelism"), "--sample-parallelism");
        } else if (arg == "--no-predicted-rt") {
            config.predict_rt = false;
        } else if (arg == "--sip-isotope") {
            config.sip_isotope = value("--sip-isotope");
        } else if (arg == "--fixed-ptm") {
            config.fixed_ptm_selectors.push_back(value("--fixed-ptm"));
        } else if (arg == "--ptm") {
            config.ptm_selectors.push_back(value("--ptm"));
        } else if (arg == "--max-ptm-count") {
            config.max_ptm_count = static_cast<int>(
                unsigned_number(value("--max-ptm-count"), "--max-ptm-count"));
        } else if (arg == "--fragment-ppm") {
            config.fragment_ppm = number(value("--fragment-ppm"), "--fragment-ppm");
        } else if (arg == "--product-top-isotopes") {
            config.product_top_isotopes =
                unsigned_number(value("--product-top-isotopes"),
                                "--product-top-isotopes");
        } else if (arg == "--quant-mz-ppm") {
            config.quant_mz_ppm =
                number(value("--quant-mz-ppm"), "--quant-mz-ppm");
        } else if (arg == "--quant-rt-window") {
            config.quant_rt_window =
                number(value("--quant-rt-window"), "--quant-rt-window");
        } else if (arg == "--quant-min-isotopes") {
            config.quant_min_isotopes =
                unsigned_number(value("--quant-min-isotopes"),
                                "--quant-min-isotopes");
        } else if (arg == "--quant-top-isotopes") {
            config.quant_top_isotopes =
                unsigned_number(value("--quant-top-isotopes"),
                                "--quant-top-isotopes");
        } else if (arg == "--quant-min-scans") {
            config.quant_min_scans =
                unsigned_number(value("--quant-min-scans"),
                                "--quant-min-scans");
        } else if (arg == "--quant-intensity-mode") {
            config.quant_intensity_mode =
                unsigned_number(value("--quant-intensity-mode"),
                                "--quant-intensity-mode");
        } else if (arg == "--no-quant-normalization") {
            config.quant_normalize = false;
        } else if (arg == "--no-mbr") {
            config.mbr = false;
        } else if (arg == "--mbr-rt-window") {
            config.mbr_rt_window =
                number(value("--mbr-rt-window"), "--mbr-rt-window");
        } else if (arg == "--mbr-top-runs") {
            config.mbr_top_runs =
                unsigned_number(value("--mbr-top-runs"), "--mbr-top-runs");
        } else if (arg == "--mbr-min-correlation") {
            config.mbr_min_correlation = number(
                value("--mbr-min-correlation"), "--mbr-min-correlation");
        } else if (arg == "--mbr-ion-fdr") {
            config.mbr_ion_fdr =
                number(value("--mbr-ion-fdr"), "--mbr-ion-fdr");
        } else if (arg == "--mbr-sip-bin-width") {
            config.mbr_sip_bin_width = number(
                value("--mbr-sip-bin-width"), "--mbr-sip-bin-width");
        } else if (arg == "--decoy-prefix") {
            config.decoy_prefix = value("--decoy-prefix");
        } else if (arg == "--protein-reference") {
            config.protein_reference_path = value("--protein-reference");
        } else if (arg == "--sip-protein-output") {
            config.sip_protein_output_path = value("--sip-protein-output");
        } else if (arg == "--negative-control" ||
                   arg == "--negative_control") {
            std::istringstream names(value(arg.c_str()));
            std::string name;
            while (std::getline(names, name, ',')) {
                const auto begin = name.find_first_not_of(" \t\r\n");
                const auto end = name.find_last_not_of(" \t\r\n");
                if (begin == std::string::npos) {
                    throw std::runtime_error(
                        "Empty sample in " + arg);
                }
                config.negative_control_samples.push_back(
                    name.substr(begin, end - begin + 1));
            }
        } else if (arg == "--label-threshold") {
            config.label_threshold = number(
                value("--label-threshold"), "--label-threshold");
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
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "aerith: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
