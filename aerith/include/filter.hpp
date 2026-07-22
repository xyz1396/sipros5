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

struct Config {
    std::vector<std::string> inputs;
    std::vector<std::string> target_pins;
    std::vector<std::string> decoy_pins;
    std::vector<std::string> output_prefixes;
    std::vector<std::string> spectrum_paths;
    std::string database_path;
    std::string spectrum_model_path;
    std::string rt_model_path;
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
    bool ignore_pct = false;
};

struct Summary {
    std::vector<std::string> input_paths;
    std::vector<std::string> output_prefixes;
    std::size_t files = 0;
    std::size_t psms = 0;
    std::size_t targets = 0;
    std::size_t decoys = 0;
    std::size_t rt_training_targets = 0;
    std::size_t target_ids = 0;
    std::size_t distinct_target_peptides = 0;
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
    StageTiming write_timing;
    StageTiming total_timing;
    double omp_speedup_ratio = 1.0;
    double omp_parallel_efficiency = 1.0;
    std::string score_model = "global_rt_samplewise_omp_simd_l2_svm_3fold";
    unsigned int threads = 1;
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
