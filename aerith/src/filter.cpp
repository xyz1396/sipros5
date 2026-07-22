#include "filter.hpp"
#include "pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>

namespace aerith {

Summary run(const Config& config) {
    using Clock = std::chrono::steady_clock;
    const auto elapsed_seconds = [](Clock::time_point begin, Clock::time_point end) {
        return std::chrono::duration<double>(end - begin).count();
    };
    const auto cpu_seconds = [](std::clock_t begin, std::clock_t end) {
        return static_cast<double>(end - begin) / CLOCKS_PER_SEC;
    };
    const auto total_begin = Clock::now();
    const std::clock_t cpu_begin = std::clock();
    const auto read_begin = Clock::now();
    const std::clock_t read_cpu_begin = std::clock();
    auto data = PinReader::read(config);
    const std::clock_t read_cpu_end = std::clock();
    const auto read_end = Clock::now();

    // This stage intentionally covers the whole generated entropy feature:
    // unique peptide-charge preparation, Torch spectrum prediction, HDF5
    // spectrum loading, fragment matching, and entropy calculation.
    const auto spectrum_entropy_begin = Clock::now();
    const std::clock_t spectrum_entropy_cpu_begin = std::clock();
    SpectralEntropyFeature::add(config, data);
    const std::clock_t spectrum_entropy_cpu_end = std::clock();
    const auto spectrum_entropy_end = Clock::now();

    // DIA-NN RT inference and the sample-specific MSBooster-compatible LOESS
    // calibration are timed together as one generated-feature stage.
    const auto predicted_rt_begin = Clock::now();
    const std::clock_t predicted_rt_cpu_begin = std::clock();
    PredictedRetentionTimeFeature::add(config, data);
    const std::clock_t predicted_rt_cpu_end = std::clock();
    const auto predicted_rt_end = Clock::now();
    const bool has_diann_rt = std::find(
        data.feature_names.begin(), data.feature_names.end(),
        "delta_RT_loess") != data.feature_names.end();
    const bool use_internal_rt_model = !has_diann_rt;

    Summary summary;
    summary.input_paths = data.input_paths;
    summary.output_prefixes = config.output_prefixes;
    summary.threads = static_cast<unsigned int>(omp_get_max_threads());
    summary.feature_names = data.feature_names;
    summary.used_internal_rt_model = use_internal_rt_model;
    summary.spectrum_prediction_device = data.spectrum_prediction_device;
    summary.rt_prediction_device = data.rt_prediction_device;
    summary.spectrum_prediction_timing = data.spectrum_prediction_timing;
    summary.rt_prediction_inference_timing = data.rt_prediction_timing;
    if (use_internal_rt_model) {
        summary.feature_names.push_back("sqrtAbsDeltaRT");
    } else {
        summary.score_model =
            "global_diann_rt_samplewise_omp_simd_l2_svm_3fold";
    }
    summary.files = data.input_paths.size();
    summary.psms = data.rows.size();
    for (const auto& row : data.rows) {
        row.label == 1 ? ++summary.targets : ++summary.decoys;
    }

    std::vector<int> labels;
    std::vector<double> initial_scores;
    labels.reserve(data.rows.size());
    initial_scores.reserve(data.rows.size());
    for (const auto& row : data.rows) {
        labels.push_back(row.label);
        initial_scores.push_back(row.initial_score);
    }

    const auto fold_begin = Clock::now();
    const std::clock_t fold_cpu_begin = std::clock();
    const auto outer_folds = SvmRescorer::assign_folds(data);
    std::vector<double> diagnostic_q;
    if (use_internal_rt_model) {
        // Diagnostic only: every fitted value is recomputed inside training folds.
        diagnostic_q = target_decoy_qvalues(initial_scores, labels);
    }
    const std::clock_t fold_cpu_end = std::clock();
    const auto fold_end = Clock::now();

    const auto rt_begin = Clock::now();
    const std::clock_t rt_cpu_begin = std::clock();
    RtResult rt;
    if (use_internal_rt_model) {
        rt = RetentionTimeModel::fit(
            data, outer_folds, initial_scores, labels, diagnostic_q,
            config.train_fdr, config.rt_ridge);
        const bool complete_rt = std::all_of(
            rt.residuals.begin(), rt.residuals.end(),
            [&](const auto& residuals) {
                return residuals.size() == data.rows.size();
            });
        if (!complete_rt || !std::isfinite(rt.r2)) {
            throw std::runtime_error("Unable to fit the nested chemical RT models");
        }
        summary.rt_training_targets = rt.training_count;
        summary.rt_r2 = rt.r2;
    }
    const std::clock_t rt_cpu_end = std::clock();
    const auto rt_end = Clock::now();

    const auto svm_begin = Clock::now();
    const std::clock_t svm_cpu_begin = std::clock();
    auto fitted = SvmRescorer::fit(
        data, rt.residuals, outer_folds, config.train_fdr,
        config.max_iterations, config.svm_c_pos, config.svm_c_neg);
    std::vector<double> scores = std::move(fitted.scores);
    const std::clock_t svm_cpu_end = std::clock();
    const auto svm_end = Clock::now();

    const auto statistics_begin = Clock::now();
    const std::clock_t statistics_cpu_begin = std::clock();
    std::vector<double> q(data.rows.size(), 1.0), pep(data.rows.size(), 1.0);
    summary.sample_models.resize(data.input_paths.size());
    std::unordered_set<std::string> peptides;
    for (std::size_t file = 0; file < data.input_paths.size(); ++file) {
        auto& sample = summary.sample_models[file];
        sample.name = config.output_prefixes.size() == data.input_paths.size()
            ? std::filesystem::path(config.output_prefixes[file]).filename().string()
            : std::filesystem::path(data.input_paths[file]).stem().string();
        sample.score_iterations = fitted.iterations[file];
        sample.feature_weights = std::move(fitted.calibrated_weights[file]);
        std::vector<std::size_t> rows;
        std::vector<double> sample_scores;
        std::vector<int> sample_labels;
        for (std::size_t i = 0; i < data.rows.size(); ++i) {
            if (data.rows[i].file_id != file) continue;
            rows.push_back(i);
            sample_scores.push_back(scores[i]);
            sample_labels.push_back(labels[i]);
        }
        sample.psms = rows.size();
        const auto sample_q = mixmax_qvalues(
            sample_scores, sample_labels, &sample.pi0);
        const auto sample_pep = SvmRescorer::local_error_probabilities(
            sample_scores, sample_labels);
        std::unordered_set<std::string> sample_peptides;
        for (std::size_t local = 0; local < rows.size(); ++local) {
            const auto i = rows[local];
            q[i] = sample_q[local];
            pep[i] = sample_pep[local];
            if (data.rows[i].label == 1 && q[i] <= config.q_threshold) {
                ++sample.target_ids;
                const auto peptide = stripped_peptide(data.rows[i].peptide);
                sample_peptides.insert(peptide);
                peptides.insert(peptide);
            }
        }
        sample.distinct_target_peptides = sample_peptides.size();
        summary.target_ids += sample.target_ids;
    }
    summary.distinct_target_peptides = peptides.size();
    const std::clock_t statistics_cpu_end = std::clock();
    const auto statistics_end = Clock::now();

    const auto write_begin = Clock::now();
    const std::clock_t write_cpu_begin = std::clock();
    ResultWriter::write(config, data, scores, q, pep, rt, outer_folds);
    const std::clock_t write_cpu_end = std::clock();
    const auto write_end = Clock::now();
    const std::clock_t cpu_end = std::clock();

    summary.read_timing = {
        elapsed_seconds(read_begin, read_end),
        cpu_seconds(read_cpu_begin, read_cpu_end)};
    summary.spectrum_entropy_timing = {
        elapsed_seconds(spectrum_entropy_begin, spectrum_entropy_end),
        cpu_seconds(spectrum_entropy_cpu_begin, spectrum_entropy_cpu_end)};
    summary.predicted_rt_timing = {
        elapsed_seconds(predicted_rt_begin, predicted_rt_end),
        cpu_seconds(predicted_rt_cpu_begin, predicted_rt_cpu_end)};
    summary.fold_setup_timing = {
        elapsed_seconds(fold_begin, fold_end),
        cpu_seconds(fold_cpu_begin, fold_cpu_end)};
    summary.rt_model_timing = {
        elapsed_seconds(rt_begin, rt_end),
        cpu_seconds(rt_cpu_begin, rt_cpu_end)};
    summary.svm_model_timing = {
        elapsed_seconds(svm_begin, svm_end),
        cpu_seconds(svm_cpu_begin, svm_cpu_end)};
    summary.statistics_timing = {
        elapsed_seconds(statistics_begin, statistics_end),
        cpu_seconds(statistics_cpu_begin, statistics_cpu_end)};
    summary.write_timing = {
        elapsed_seconds(write_begin, write_end),
        cpu_seconds(write_cpu_begin, write_cpu_end)};
    summary.total_timing = {
        elapsed_seconds(total_begin, write_end),
        cpu_seconds(cpu_begin, cpu_end)};
    if (summary.total_timing.wall_seconds > 0.0) {
        summary.omp_speedup_ratio =
            summary.total_timing.cpu_seconds / summary.total_timing.wall_seconds;
        summary.omp_parallel_efficiency =
            summary.omp_speedup_ratio / std::max(1u, summary.threads);
    }
    return summary;
}

} // namespace aerith
