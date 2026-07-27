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

std::string peptide_form(const std::string& peptide) {
    const auto open = peptide.find('[');
    const auto close = peptide.rfind(']');
    if (open != std::string::npos && close != std::string::npos &&
        close > open) {
        return peptide.substr(open + 1, close - open - 1);
    }
    return peptide;
}

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
    summary.filtered_only = config.filtered_only;
    summary.reporting_fdr = config.q_threshold;
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
    summary.removed_decoy_peptide_collisions =
        data.removed_decoy_peptide_collisions;
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
        for (std::size_t local = 0; local < rows.size(); ++local) {
            const auto i = rows[local];
            q[i] = sample_q[local];
            pep[i] = sample_pep[local];
        }
    }
    if (config.assemble_proteins && !config.database_path.empty() &&
        !config.decoy_database_path.empty()) {
        ProteinAssembler::sequential_filter(config, data, scores, pep, q);
    }

    std::unordered_set<std::string> peptides;
    std::unordered_set<std::string> peptide_forms;
    std::unordered_set<std::string> ptm_peptides;
    std::vector<std::unordered_set<std::string>> sample_peptides(
        data.input_paths.size());
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto& psm = data.rows[row];
        if (psm.label != 1 || q[row] > config.q_threshold ||
            psm.file_id >= summary.sample_models.size()) {
            continue;
        }
        auto& sample = summary.sample_models[psm.file_id];
        ++sample.target_ids;
        const auto peptide = stripped_peptide(psm.peptide);
        sample_peptides[psm.file_id].insert(peptide);
        peptides.insert(peptide);
        const auto form = peptide_form(psm.peptide);
        peptide_forms.insert(form);
        if (psm.ptm_count > 0) {
            ptm_peptides.insert(form);
            ++summary.target_ptm_psms;
        }
    }
    for (std::size_t file = 0; file < summary.sample_models.size(); ++file) {
        summary.sample_models[file].distinct_target_peptides =
            sample_peptides[file].size();
        summary.target_ids += summary.sample_models[file].target_ids;
    }
    summary.distinct_target_peptides = peptides.size();
    summary.distinct_target_peptide_forms = peptide_forms.size();
    summary.distinct_target_ptm_peptides = ptm_peptides.size();
    const std::clock_t statistics_cpu_end = std::clock();
    const auto statistics_end = Clock::now();

    const auto quantification_begin = Clock::now();
    const std::clock_t quantification_cpu_begin = std::clock();
    auto quantification = ChromatographicQuantifier::add(config, data, q);
    summary.mbr_ions = data.transferred_ions.size();
    summary.quantification_stages = std::move(quantification.stages);
    const std::clock_t quantification_cpu_end = std::clock();
    const auto quantification_end = Clock::now();

    const auto write_begin = Clock::now();
    const std::clock_t write_cpu_begin = std::clock();
    ResultWriter::write(config, data, scores, q, pep, rt, outer_folds);
    const std::clock_t write_cpu_end = std::clock();
    const auto write_end = Clock::now();

    const auto protein_begin = Clock::now();
    const std::clock_t protein_cpu_begin = std::clock();
    if (config.assemble_proteins && !config.database_path.empty() &&
        !config.decoy_database_path.empty()) {
        const auto assembly =
            ProteinAssembler::write(config, data, scores, q, pep);
        summary.protein_ids = assembly.proteins;
        summary.protein_output_dir = assembly.output_dir;
        summary.protein_assembly_stages = assembly.stages;
    }
    const std::clock_t protein_cpu_end = std::clock();
    const auto protein_end = Clock::now();
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
    summary.quantification_timing = {
        elapsed_seconds(quantification_begin, quantification_end),
        cpu_seconds(
            quantification_cpu_begin, quantification_cpu_end)};
    summary.write_timing = {
        elapsed_seconds(write_begin, write_end),
        cpu_seconds(write_cpu_begin, write_cpu_end)};
    summary.protein_assembly_timing = {
        elapsed_seconds(protein_begin, protein_end),
        cpu_seconds(protein_cpu_begin, protein_cpu_end)};
    summary.total_timing = {
        elapsed_seconds(total_begin, protein_end),
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
