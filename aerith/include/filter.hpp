#pragma once

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace aerith {

struct SampleModelSummary {
    std::string name;
    std::size_t psms = 0;
    std::size_t target_ids = 0;
    std::size_t distinct_target_peptides = 0;
    double pi0 = 1.0;
    std::array<unsigned int, 3> score_iterations{};
    std::array<std::vector<double>, 3> feature_weights;
};

struct StageTiming {
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
};

struct AccelerationTiming {
    std::string name;
    StageTiming timing;
    bool uses_omp = false;
    bool uses_simd = false;
};

struct Config {
    std::vector<std::string> inputs;
    std::vector<std::string> target_pins;
    std::vector<std::string> decoy_pins;
    std::vector<std::string> output_prefixes;
    std::vector<std::string> spectrum_paths;
    std::string database_path;
    std::string decoy_database_path;
    std::string protein_output_dir;
    std::string spectrum_model_path;
    std::string rt_model_path;
    std::string prediction_cache_path;
    std::string sip_isotope;
    std::vector<std::string> fixed_ptm_selectors;
    std::vector<std::string> ptm_selectors;
    int max_ptm_count = -1;
    std::string protein_reference_path;
    std::string sip_protein_output_path;
    std::vector<std::string> negative_control_samples;
    double label_threshold = 2.0;
    bool predict_rt = true;
    std::string decoy_prefix = "Decoy_";
    std::string initial_score = "WDPscores";
    double q_threshold = 0.01;
    double train_fdr = 0.01;
    double rt_ridge = 1e-4;
    double svm_c_pos = 1.0;
    double svm_c_neg = 1.0;
    unsigned int max_iterations = 10;
    double fragment_ppm = 20.0;
    unsigned int product_top_isotopes = 5;
    bool ignore_pct = false;
    bool assemble_proteins = true;
    bool filtered_only = false;
    bool fixed_cam = true;
    double quant_mz_ppm = 10.0;
    double quant_rt_window = 0.4;
    unsigned int quant_min_isotopes = 2;
    unsigned int quant_top_isotopes = 6;
    unsigned int quant_min_scans = 3;
    unsigned int quant_intensity_mode = 2;
    bool quant_normalize = true;
    bool mbr = true;
    bool stream_samples = true;
    unsigned int sample_parallelism = 3;
    double mbr_rt_window = 1.0;
    unsigned int mbr_top_runs = 10;
    double mbr_min_correlation = 0.0;
    double mbr_ion_fdr = 0.01;
    double mbr_sip_bin_width = 10.0;
};

struct Summary {
    std::vector<std::string> input_paths;
    std::vector<std::string> output_prefixes;
    std::size_t files = 0;
    std::size_t psms = 0;
    std::size_t targets = 0;
    std::size_t decoys = 0;
    std::size_t removed_decoy_peptide_collisions = 0;
    std::size_t rt_training_targets = 0;
    std::size_t target_ids = 0;
    std::size_t distinct_target_peptides = 0;
    std::size_t distinct_target_peptide_forms = 0;
    std::size_t distinct_target_ptm_peptides = 0;
    std::size_t target_ptm_psms = 0;
    std::size_t protein_ids = 0;
    std::size_t mbr_ions = 0;
    std::size_t negative_control_candidates = 0;
    std::size_t negative_control_input_psms = 0;
    std::size_t negative_control_threshold_filtered_psms = 0;
    std::size_t negative_control_targets = 0;
    std::size_t negative_control_decoys = 0;
    std::size_t negative_control_target_ids = 0;
    double negative_control_label_threshold = 0.0;
    std::string negative_control_output_path;
    std::string negative_control_target_output_path;
    std::string negative_control_decoy_output_path;
    std::string negative_control_protein_output_path;
    StageTiming negative_control_timing;
    std::vector<AccelerationTiming> negative_control_stages;
    std::vector<std::string> negative_control_feature_names;
    SampleModelSummary negative_control_model;
    std::string protein_output_dir;
    bool filtered_only = false;
    double reporting_fdr = 0.01;
    double rt_r2 = 0.0;
    bool used_internal_rt_model = false;
    std::string spectrum_prediction_device;
    std::string rt_prediction_device;
    StageTiming read_timing;
    StageTiming spectrum_prediction_timing;
    StageTiming spectrum_entropy_timing;
    StageTiming rt_prediction_inference_timing;
    StageTiming predicted_rt_timing;
    StageTiming fold_setup_timing;
    StageTiming rt_model_timing;
    StageTiming svm_model_timing;
    StageTiming statistics_timing;
    StageTiming quantification_timing;
    std::vector<AccelerationTiming> quantification_stages;
    StageTiming write_timing;
    StageTiming protein_assembly_timing;
    std::vector<AccelerationTiming> protein_assembly_stages;
    StageTiming total_timing;
    double omp_speedup_ratio = 1.0;
    double omp_parallel_efficiency = 1.0;
    std::string score_model = "global_rt_samplewise_omp_simd_l2_svm_3fold";
    unsigned int threads = 1;
    unsigned int sample_parallelism = 1;
    std::vector<std::string> feature_names;
    std::vector<SampleModelSummary> sample_models;
};

Summary run(const Config& config);
void print_summary(std::ostream& output, const Summary& summary);

std::vector<double> target_decoy_qvalues(const std::vector<double>& scores,
                                         const std::vector<int>& labels);
std::vector<double> mixmax_qvalues(const std::vector<double>& scores,
                                   const std::vector<int>& labels,
                                   double* pi0 = nullptr);
std::string stripped_peptide(const std::string& peptide);

} // namespace aerith
