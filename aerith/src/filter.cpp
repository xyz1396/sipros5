#include "filter.hpp"
#include "isotope.hpp"
#include "pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace aerith {

namespace {

using Clock = std::chrono::steady_clock;

constexpr std::uint64_t kPredictionCacheMagic = 0x4145525052454431ULL;
constexpr std::uint32_t kPredictionCacheVersion = 2;

void prediction_hash_bytes(
    std::uint64_t& hash, const void* bytes, std::size_t size) {
    const auto* data = static_cast<const unsigned char*>(bytes);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= data[index];
        hash *= 1099511628211ULL;
    }
}

template <typename Value>
void prediction_hash_value(std::uint64_t& hash, const Value& value) {
    prediction_hash_bytes(hash, &value, sizeof(value));
}

void prediction_hash_string(std::uint64_t& hash, const std::string& value) {
    prediction_hash_bytes(hash, value.data(), value.size());
    constexpr unsigned char separator = 0xff;
    prediction_hash_bytes(hash, &separator, 1);
}

std::uint64_t prediction_cache_fingerprint(
    const Config& config, const Dataset& data) {
    std::uint64_t hash = 1469598103934665603ULL;
    prediction_hash_value(hash, config.fragment_ppm);
    prediction_hash_value(hash, config.product_top_isotopes);
    prediction_hash_value(hash, config.predict_rt);
    prediction_hash_string(hash, config.spectrum_model_path);
    prediction_hash_string(hash, config.rt_model_path);
    prediction_hash_string(hash, config.sip_isotope);
    for (const auto& path : config.spectrum_paths) {
        prediction_hash_string(hash, path);
    }
    for (const auto& row : data.rows) {
        prediction_hash_string(hash, row.id);
        prediction_hash_string(hash, row.peptide);
        prediction_hash_value(hash, row.file_id);
        prediction_hash_value(hash, row.scan);
        prediction_hash_value(hash, row.charge);
        prediction_hash_value(hash, row.retention);
    }
    return hash;
}

template <typename Value>
bool read_prediction_value(std::istream& input, Value& value) {
    return static_cast<bool>(input.read(
        reinterpret_cast<char*>(&value), sizeof(value)));
}

template <typename Value>
void write_prediction_value(std::ostream& output, const Value& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

bool load_prediction_cache(
    const Config& config, Dataset& data, std::size_t base_feature_count) {
    if (config.prediction_cache_path.empty() ||
        !std::filesystem::is_regular_file(config.prediction_cache_path)) {
        return false;
    }
    std::ifstream input(config.prediction_cache_path, std::ios::binary);
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    std::uint64_t fingerprint = 0;
    std::uint64_t row_count = 0;
    std::uint64_t cached_base_features = 0;
    std::uint32_t generated_count = 0;
    std::uint8_t diagnostics = 0;
    if (!read_prediction_value(input, magic) ||
        !read_prediction_value(input, version) ||
        !read_prediction_value(input, fingerprint) ||
        !read_prediction_value(input, row_count) ||
        !read_prediction_value(input, cached_base_features) ||
        !read_prediction_value(input, generated_count) ||
        !read_prediction_value(input, diagnostics) ||
        magic != kPredictionCacheMagic ||
        version != kPredictionCacheVersion ||
        fingerprint != prediction_cache_fingerprint(config, data) ||
        row_count != data.rows.size() ||
        cached_base_features != base_feature_count ||
        generated_count == 0 || generated_count > 16) {
        return false;
    }
    std::vector<std::string> names;
    names.reserve(generated_count);
    for (std::uint32_t index = 0; index < generated_count; ++index) {
        std::uint32_t length = 0;
        if (!read_prediction_value(input, length) || length > 1024) {
            return false;
        }
        std::string name(length, '\0');
        if (!input.read(name.data(), length)) return false;
        names.push_back(std::move(name));
    }
    const std::size_t values_per_row = generated_count + 2;
    if (data.rows.size() >
        std::numeric_limits<std::size_t>::max() / values_per_row) {
        return false;
    }
    std::vector<float> values(data.rows.size() * values_per_row);
    if (!input.read(
            reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(float)))) {
        return false;
    }
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        auto& psm = data.rows[row];
        if (psm.features.size() != base_feature_count) return false;
        const auto offset = row * values_per_row;
        psm.features.insert(
            psm.features.end(), values.begin() + offset,
            values.begin() + offset + generated_count);
        psm.delta_rt_loess_real = values[offset + generated_count];
        psm.predicted_rt_real_units = values[offset + generated_count + 1];
    }
    data.feature_names.insert(
        data.feature_names.end(), names.begin(), names.end());
    data.generated_feature_names.insert(
        data.generated_feature_names.end(), names.begin(), names.end());
    data.has_predicted_rt_diagnostics = diagnostics != 0;
    data.spectrum_prediction_device = "cache";
    data.rt_prediction_device = "cache";
    return true;
}

void save_prediction_cache(
    const Config& config, const Dataset& data,
    std::size_t base_feature_count) {
    if (config.prediction_cache_path.empty()) return;
    const std::size_t generated_count = data.generated_feature_names.size();
    if (generated_count == 0 || generated_count > 16) return;
    for (const auto& row : data.rows) {
        if (row.features.size() != base_feature_count + generated_count) {
            throw std::runtime_error(
                "Cannot cache inconsistent generated prediction features");
        }
    }
    const std::filesystem::path path(config.prediction_cache_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const auto temporary = path.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "Cannot write prediction cache: " + temporary);
    }
    write_prediction_value(output, kPredictionCacheMagic);
    write_prediction_value(output, kPredictionCacheVersion);
    write_prediction_value(
        output, prediction_cache_fingerprint(config, data));
    write_prediction_value(
        output, static_cast<std::uint64_t>(data.rows.size()));
    write_prediction_value(
        output, static_cast<std::uint64_t>(base_feature_count));
    write_prediction_value(
        output, static_cast<std::uint32_t>(generated_count));
    write_prediction_value(
        output, static_cast<std::uint8_t>(
            data.has_predicted_rt_diagnostics));
    for (const auto& name : data.generated_feature_names) {
        write_prediction_value(
            output, static_cast<std::uint32_t>(name.size()));
        output.write(name.data(), static_cast<std::streamsize>(name.size()));
    }
    for (const auto& row : data.rows) {
        output.write(
            reinterpret_cast<const char*>(
                row.features.data() + base_feature_count),
            static_cast<std::streamsize>(generated_count * sizeof(float)));
        write_prediction_value(output, row.delta_rt_loess_real);
        write_prediction_value(output, row.predicted_rt_real_units);
    }
    output.close();
    if (!output) {
        throw std::runtime_error(
            "Failed while writing prediction cache: " + temporary);
    }
    std::error_code rename_error;
    std::filesystem::rename(temporary, path, rename_error);
    if (rename_error) {
        std::filesystem::remove(path, rename_error);
        rename_error.clear();
        std::filesystem::rename(temporary, path, rename_error);
    }
    if (rename_error) {
        throw std::runtime_error(
            "Cannot install prediction cache: " + path.string());
    }
}

} // namespace

std::string peptide_form(const std::string& peptide) {
    const auto open = peptide.find('[');
    const auto close = peptide.rfind(']');
    if (open != std::string::npos && close != std::string::npos &&
        close > open) {
        return peptide.substr(open + 1, close - open - 1);
    }
    return peptide;
}

NegativeControlResult NegativeControlFilter::run(
    const Config& config, const Dataset& source,
    const std::vector<double>& source_q) {
    using Clock = std::chrono::steady_clock;
    NegativeControlResult result;
    const auto total_begin = Clock::now();
    const auto total_cpu_begin = std::clock();
    const auto record_stage = [&](
        std::string name, Clock::time_point wall_begin,
        std::clock_t cpu_begin, bool uses_omp = false) {
        const auto wall_end = Clock::now();
        const auto cpu_end = std::clock();
        result.stages.push_back({
            std::move(name),
            {std::chrono::duration<double>(
                 wall_end - wall_begin).count(),
             static_cast<double>(cpu_end - cpu_begin) / CLOCKS_PER_SEC},
            uses_omp,
            false});
    };
    if (source.rows.size() != source_q.size()) {
        throw std::runtime_error(
            "Native negative-control filtering received inconsistent PSM arrays");
    }
    if (config.negative_control_samples.empty()) return {};
    if (config.output_prefixes.size() != source.input_paths.size()) {
        throw std::runtime_error(
            "Native negative-control filtering requires one output prefix per sample");
    }

    const auto selection_begin = Clock::now();
    const auto selection_cpu_begin = std::clock();
    std::vector<std::string> sample_names;
    sample_names.reserve(config.output_prefixes.size());
    std::unordered_map<std::string, std::size_t> sample_indices;
    for (std::size_t file = 0; file < config.output_prefixes.size(); ++file) {
        const auto name = std::filesystem::path(
            config.output_prefixes[file]).filename().string();
        if (name.empty() || !sample_indices.emplace(name, file).second) {
            throw std::runtime_error(
                "Negative-control sample basenames must be nonempty and unique");
        }
        sample_names.push_back(name);
    }
    std::unordered_set<std::size_t> controls;
    for (const auto& name : config.negative_control_samples) {
        const auto found = sample_indices.find(name);
        if (found == sample_indices.end()) {
            throw std::runtime_error(
                "Negative-control sample not found in input basenames: " + name);
        }
        controls.insert(found->second);
    }
    if (controls.size() == sample_names.size()) {
        throw std::runtime_error(
            "At least one non-control sample is required for negative-control filtering");
    }

    Dataset data;
    data.input_paths = {"native negative-control PSMs"};
    data.headers = source.headers;
    data.columns = source.columns;
    data.numeric_columns = source.numeric_columns;
    data.feature_names = source.feature_names;
    data.generated_feature_names = source.generated_feature_names;
    data.has_predicted_rt_diagnostics =
        source.has_predicted_rt_diagnostics;
    const bool append_sample_name =
        data.columns.find("SampleName") == data.columns.end();
    if (append_sample_name) {
        data.columns.emplace("SampleName", data.headers.size());
        data.headers.push_back("SampleName");
    }
    for (std::size_t row = 0; row < source.rows.size(); ++row) {
        const auto& psm = source.rows[row];
        if (psm.label != 1 || source_q[row] > config.q_threshold ||
            psm.file_id >= sample_names.size()) {
            continue;
        }
        ++result.input_psms;
        if (psm.ms2_isotopic_abundance < config.label_threshold) {
            ++result.threshold_filtered_psms;
            continue;
        }
        data.rows.push_back(psm);
        auto& selected = data.rows.back();
        selected.sample_name = sample_names[psm.file_id];
        selected.label = controls.count(psm.file_id) == 0 ? 1 : -1;
        selected.file_id = 0;
        selected.ms1_isotopic_abundance =
            std::clamp(selected.ms1_isotopic_abundance, 0.0, 100.0);
    }
    result.candidates = data.rows.size();
    record_stage(
        "Select and relabel accepted in-memory PSMs",
        selection_begin, selection_cpu_begin);
    if (data.rows.empty()) {
        throw std::runtime_error(
            "No accepted SIP PSMs meet the negative-control label threshold");
    }

    for (const auto& psm : data.rows) {
        psm.label == 1 ? ++result.targets : ++result.decoys;
    }
    if (result.targets == 0 || result.decoys == 0) {
        throw std::runtime_error(
            "Native negative-control filtering requires target and control PSMs "
            "after primary target/decoy filtering");
    }

    const auto rt_feature = std::find(
        data.feature_names.begin(), data.feature_names.end(),
        "delta_RT_loess");
    if (rt_feature == data.feature_names.end()) {
        throw std::runtime_error(
            "Native negative-control filtering requires the primary "
            "delta_RT_loess feature; RT recomputation is disabled");
    }

    const auto fold_begin = Clock::now();
    const auto fold_cpu_begin = std::clock();
    std::vector<int> labels;
    labels.reserve(data.rows.size());
    for (const auto& psm : data.rows) {
        labels.push_back(psm.label);
    }
    const auto outer_folds = SvmRescorer::assign_folds(data);
    record_stage(
        "Assign SIP-Negative-control folds",
        fold_begin, fold_cpu_begin, true);
    RtResult rt;
    const auto svm_begin = Clock::now();
    const auto svm_cpu_begin = std::clock();
    auto fitted = SvmRescorer::fit(
        data, rt.residuals, outer_folds, config.train_fdr,
        config.max_iterations, config.svm_c_pos, config.svm_c_neg);
    auto scores = std::move(fitted.scores);
    record_stage(
        "Fit and score SIP-Negative-control SVM folds",
        svm_begin, svm_cpu_begin, true);
    const auto statistics_begin = Clock::now();
    const auto statistics_cpu_begin = std::clock();
    double pi0 = 1.0;
    const auto q = mixmax_qvalues(scores, labels, &pi0);
    const auto pep =
        SvmRescorer::local_error_probabilities(scores, labels);
    std::unordered_set<std::string> accepted_peptides;
    std::unordered_set<std::string> accepted_donor_psm_ids;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        if (data.rows[row].label == 1 && q[row] <= config.q_threshold) {
            ++result.target_ids;
            accepted_peptides.insert(
                stripped_peptide(data.rows[row].peptide));
            accepted_donor_psm_ids.insert(data.rows[row].id);
        }
    }
    for (const auto& transfer : source.transferred_ions) {
        if (accepted_donor_psm_ids.count(
                transfer.donor_psm_id) == 0 ||
            transfer.psm.file_id >= sample_names.size() ||
            controls.count(transfer.psm.file_id) != 0) {
            continue;
        }
        data.transferred_ions.push_back(transfer);
        data.transferred_ions.back().psm.sample_name =
            sample_names[transfer.psm.file_id];
    }
    result.feature_names = data.feature_names;
    result.model.name = "SIP-Negative-control";
    result.model.psms = data.rows.size();
    result.model.target_ids = result.target_ids;
    result.model.distinct_target_peptides = accepted_peptides.size();
    result.model.pi0 = pi0;
    result.model.score_iterations = fitted.iterations.at(0);
    result.model.feature_weights =
        std::move(fitted.calibrated_weights.at(0));
    record_stage(
        "Compute SIP-Negative-control q-values and PEPs",
        statistics_begin, statistics_cpu_begin);

    const auto write_begin = Clock::now();
    const auto write_cpu_begin = std::clock();
    const std::filesystem::path output_dir(config.protein_output_dir);
    Config output_config = config;
    output_config.inputs.clear();
    output_config.target_pins.clear();
    output_config.decoy_pins.clear();
    output_config.spectrum_paths.clear();
    output_config.output_prefixes = {(output_dir / "SIP").string()};
    output_config.filtered_only = false;
    output_config.assemble_proteins = false;
    output_config.protein_reference_path =
        (output_dir / "combined_protein.tsv").string();
    output_config.sip_protein_output_path =
        (output_dir / "combined_protein_with_SIP_filtered_PSM.tsv").string();
    ResultWriter::write(
        output_config, data, scores, q, pep, rt, outer_folds);
    ProteinAssembler::write_sip_psm_mapping(
        output_config, data, q);
    result.output_path =
        output_config.output_prefixes.front() + "_filtered_psms.tsv";
    result.target_output_path =
        output_config.output_prefixes.front() + "_target_psms.tsv";
    result.decoy_output_path =
        output_config.output_prefixes.front() + "_decoy_psms.tsv";
    result.protein_output_path =
        output_config.sip_protein_output_path;
    record_stage(
        "Write native SIP-Negative-control reports",
        write_begin, write_cpu_begin, true);
    const auto total_end = Clock::now();
    const auto total_cpu_end = std::clock();
    result.timing = {
        std::chrono::duration<double>(total_end - total_begin).count(),
        static_cast<double>(total_cpu_end - total_cpu_begin) /
            CLOCKS_PER_SEC};
    return result;
}

Summary run_monolithic(const Config& config) {
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
    initialize_sip_isotope_model(config);
    const std::clock_t read_cpu_end = std::clock();
    const auto read_end = Clock::now();

    // This stage intentionally covers the whole generated entropy feature:
    // unique peptide-charge preparation, Torch spectrum prediction, HDF5
    // spectrum loading, fragment matching, and entropy calculation.
    const auto spectrum_entropy_begin = Clock::now();
    const std::clock_t spectrum_entropy_cpu_begin = std::clock();
    const std::size_t base_feature_count = data.feature_names.size();
    const bool predictions_cached = load_prediction_cache(
        config, data, base_feature_count);
    if (!predictions_cached) {
        SpectralEntropyFeature::add(config, data);
    }
    const std::clock_t spectrum_entropy_cpu_end = std::clock();
    const auto spectrum_entropy_end = Clock::now();

    // DIA-NN RT inference and the sample-specific MSBooster-compatible LOESS
    // calibration are timed together as one generated-feature stage.
    const auto predicted_rt_begin = Clock::now();
    const std::clock_t predicted_rt_cpu_begin = std::clock();
    if (!predictions_cached) {
        PredictedRetentionTimeFeature::add(config, data);
        save_prediction_cache(config, data, base_feature_count);
    }
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
    if (!config.sip_protein_output_path.empty()) {
        ProteinAssembler::write_sip_psm_mapping(config, data, q);
    }
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
    if (!config.negative_control_samples.empty()) {
        const auto negative =
            NegativeControlFilter::run(config, data, q);
        summary.negative_control_candidates = negative.candidates;
        summary.negative_control_input_psms = negative.input_psms;
        summary.negative_control_threshold_filtered_psms =
            negative.threshold_filtered_psms;
        summary.negative_control_targets = negative.targets;
        summary.negative_control_decoys = negative.decoys;
        summary.negative_control_target_ids = negative.target_ids;
        summary.negative_control_label_threshold =
            config.label_threshold;
        summary.negative_control_output_path = negative.output_path;
        summary.negative_control_target_output_path =
            negative.target_output_path;
        summary.negative_control_decoy_output_path =
            negative.decoy_output_path;
        summary.negative_control_protein_output_path =
            negative.protein_output_path;
        summary.negative_control_timing = negative.timing;
        summary.negative_control_stages = negative.stages;
        summary.negative_control_feature_names =
            negative.feature_names;
        summary.negative_control_model = negative.model;
    }
    const auto completion_end = Clock::now();
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
        elapsed_seconds(total_begin, completion_end),
        cpu_seconds(cpu_begin, cpu_end)};
    if (summary.total_timing.wall_seconds > 0.0) {
        summary.omp_speedup_ratio =
            summary.total_timing.cpu_seconds / summary.total_timing.wall_seconds;
        summary.omp_parallel_efficiency =
            summary.omp_speedup_ratio / std::max(1u, summary.threads);
    }
    return summary;
}

namespace {

Config streamed_sample_config(const Config& config, std::size_t sample) {
    Config result = config;
    result.inputs.clear();
    result.target_pins = {config.target_pins.at(sample)};
    result.decoy_pins = {config.decoy_pins.at(sample)};
    result.output_prefixes = {config.output_prefixes.at(sample)};
    result.spectrum_paths = {config.spectrum_paths.at(sample)};
    result.prediction_cache_path.clear();
    return result;
}

struct StreamSpool {
    std::filesystem::path target;
    std::filesystem::path decoy;
};

struct StreamTemporaryDirectory {
    std::filesystem::path path;

    StreamTemporaryDirectory() {
        std::uint64_t discriminator = static_cast<std::uint64_t>(
            std::chrono::steady_clock::now().time_since_epoch().count());
#if defined(__unix__) || defined(__APPLE__)
        discriminator ^= static_cast<std::uint64_t>(::getpid()) << 32;
#endif
        path = std::filesystem::path("/dev/shm") /
            ("aerith-stream-" + std::to_string(discriminator));
        std::filesystem::create_directories(path);
    }

    ~StreamTemporaryDirectory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

struct StreamedSampleState {
    std::size_t file = 0;
    Config config;
    Dataset data;
    std::vector<std::size_t> folds;
    SvmFit fitted;
    std::vector<double> scores;
    std::vector<double> q;
    std::vector<double> pep;
    std::vector<std::size_t> compact_indices;
    SampleModelSummary model;
};

template <typename Function>
void run_stream_sample_batch(
    std::size_t count, unsigned int total_threads, Function&& function) {
    if (count == 0) return;
    if (count == 1) {
        function(0);
        return;
    }
    std::vector<std::thread> workers;
    workers.reserve(count);
    std::mutex failure_mutex;
    std::exception_ptr failure;
    for (std::size_t sample = 0; sample < count; ++sample) {
        const auto quotient = total_threads / count;
        const auto remainder = total_threads % count;
        const auto thread_budget = static_cast<unsigned int>(
            std::max<std::size_t>(1, quotient + (sample < remainder)));
        workers.emplace_back([&, sample, thread_budget] {
            omp_set_dynamic(0);
            omp_set_num_threads(static_cast<int>(thread_budget));
            try {
                function(sample);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        });
    }
    for (auto& worker : workers) worker.join();
    if (failure) std::rethrow_exception(failure);
}

void write_stream_score_spool(
    const std::filesystem::path& path, int label,
    const Dataset& data, const std::vector<double>& scores,
    const std::vector<double>& primary_q, const std::vector<double>& pep,
    const std::vector<std::size_t>& compact_indices) {
    std::vector<std::size_t> order;
    order.reserve(data.rows.size() / 2);
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        if (data.rows[row].label == label) order.push_back(row);
    }
    std::stable_sort(
        order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
            return scores[left] > scores[right];
        });
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "Cannot create streamed score spool: " + path.string());
    }
    output << std::setprecision(17);
    for (const auto row : order) {
        const auto compact = compact_indices[row];
        if (compact == std::numeric_limits<std::size_t>::max()) {
            output << '-';
        } else {
            output << compact;
        }
        const auto& psm = data.rows[row];
        output << '\t' << psm.id << '\t' << scores[row] << '\t'
               << primary_q[row] << '\t' << pep[row] << '\t'
               << psm.peptide << '\t' << psm.proteins << '\n';
    }
}

void write_stream_score_table(
    const std::filesystem::path& spool, const std::string& output_path,
    const std::vector<double>& final_q, bool protein_filtered) {
    std::ifstream input(spool);
    if (!input) {
        throw std::runtime_error(
            "Cannot read streamed score spool: " + spool.string());
    }
    const std::filesystem::path output_file(output_path);
    if (output_file.has_parent_path()) {
        std::filesystem::create_directories(output_file.parent_path());
    }
    std::ofstream output(output_file);
    if (!output) {
        throw std::runtime_error(
            "Cannot create output: " + output_file.string());
    }
    output << "PSMId\tSVMscore\tq-value\tposterior_error_prob"
           << "\tpeptide\tproteinIds\n" << std::setprecision(10);
    std::string line;
    while (std::getline(input, line)) {
        std::array<std::string_view, 7> fields{};
        std::size_t begin = 0;
        for (std::size_t field = 0; field < fields.size(); ++field) {
            const auto tab = line.find('\t', begin);
            const auto end = tab == std::string::npos ? line.size() : tab;
            fields[field] = std::string_view(line).substr(begin, end - begin);
            begin = tab == std::string::npos ? line.size() : tab + 1;
        }
        double qvalue = std::stod(std::string(fields[3]));
        if (fields[0] != "-") {
            const auto compact = static_cast<std::size_t>(
                std::stoull(std::string(fields[0])));
            if (compact >= final_q.size()) {
                throw std::runtime_error(
                    "Invalid compact PSM index in streamed score spool");
            }
            qvalue = final_q[compact];
        } else if (protein_filtered) {
            qvalue = 1.0;
        }
        output << fields[1] << '\t'
               << std::stod(std::string(fields[2])) << '\t' << qvalue
               << '\t' << std::stod(std::string(fields[4])) << '\t'
               << fields[5] << '\t' << fields[6] << '\n';
    }
}

void add_timing(StageTiming& total, Clock::time_point wall_begin,
                std::clock_t cpu_begin) {
    total.wall_seconds += std::chrono::duration<double>(
        Clock::now() - wall_begin).count();
    total.cpu_seconds += static_cast<double>(
        std::clock() - cpu_begin) / CLOCKS_PER_SEC;
}

Summary run_streamed(const Config& config) {
    const auto total_begin = Clock::now();
    const auto total_cpu_begin = std::clock();
    initialize_sip_isotope_model(config);

    Summary summary;
    summary.files = config.target_pins.size();
    summary.output_prefixes = config.output_prefixes;
    summary.filtered_only = config.filtered_only;
    summary.reporting_fdr = config.q_threshold;
    summary.threads = static_cast<unsigned int>(omp_get_max_threads());
    summary.sample_parallelism = std::max(
        1u, std::min({config.sample_parallelism, summary.threads,
                      static_cast<unsigned int>(summary.files)}));
    summary.used_internal_rt_model = false;
    summary.score_model =
        "streamed_diann_rt_samplewise_omp_simd_l2_svm_3fold";
    summary.sample_models.resize(summary.files);

    StageTiming read_timing;
    std::unordered_set<std::string> global_target_peptides;
    global_target_peptides.reserve(100000);

    // The discovery pass parses only peptide and charge columns and retains
    // one exemplar per key. It does not instantiate complete PSM features.
    auto discovery_begin = Clock::now();
    auto discovery_cpu_begin = std::clock();
    auto prediction_catalog = PinReader::discover_predictions(
        config, global_target_peptides);
    add_timing(read_timing, discovery_begin, discovery_cpu_begin);

    const auto spectrum_stage_begin = Clock::now();
    const auto spectrum_stage_cpu_begin = std::clock();
    auto spectrum_predictions =
        SpectralEntropyFeature::predict(config, prediction_catalog);
    summary.spectrum_prediction_device = spectrum_predictions.device();
    summary.spectrum_prediction_timing = spectrum_predictions.timing();
    const StageTiming spectrum_inference_outer{
        std::chrono::duration<double>(
            Clock::now() - spectrum_stage_begin).count(),
        static_cast<double>(std::clock() - spectrum_stage_cpu_begin) /
            CLOCKS_PER_SEC};

    const auto rt_prediction_begin = Clock::now();
    const auto rt_prediction_cpu_begin = std::clock();
    auto rt_predictions =
        PredictedRetentionTimeFeature::predict(config, prediction_catalog);
    summary.rt_prediction_device = rt_predictions.device();
    summary.rt_prediction_inference_timing = rt_predictions.timing();
    const StageTiming rt_inference_outer{
        std::chrono::duration<double>(
            Clock::now() - rt_prediction_begin).count(),
        static_cast<double>(std::clock() - rt_prediction_cpu_begin) /
            CLOCKS_PER_SEC};

    // Prediction maps now own the canonical strings and model outputs.
    // Release discovery-only exemplars before any complete sample is read.
    prediction_catalog = {};

    Dataset data;
    data.input_paths.reserve(summary.files);
    std::vector<double> scores;
    std::vector<double> q;
    std::vector<double> pep;
    StreamTemporaryDirectory temporary;
    std::vector<StreamSpool> spools(summary.files);
    StageTiming entropy_scoring;
    StageTiming rt_calibration;
    StageTiming fold_setup;
    StageTiming svm_model;
    StageTiming primary_statistics;

    for (std::size_t batch_begin = 0; batch_begin < summary.files;
         batch_begin += summary.sample_parallelism) {
        // At most N complete feature matrices are live. Each external sample
        // worker receives a share of the global OpenMP budget, and batches are
        // merged in input order to keep results deterministic.
        const std::size_t batch_count = std::min<std::size_t>(
            summary.sample_parallelism, summary.files - batch_begin);
        std::vector<StreamedSampleState> batch(batch_count);
        for (std::size_t local = 0; local < batch_count; ++local) {
            auto& state = batch[local];
            state.file = batch_begin + local;
            state.config = streamed_sample_config(config, state.file);
        }
        const auto run_phase = [&](StageTiming& timing, auto&& function) {
            const auto begin = Clock::now();
            const auto cpu_begin = std::clock();
            run_stream_sample_batch(
                batch_count, summary.threads,
                [&](std::size_t local) { function(batch[local]); });
            add_timing(timing, begin, cpu_begin);
        };

        run_phase(read_timing, [&](StreamedSampleState& state) {
            state.data = PinReader::read(
                state.config, &global_target_peptides);
        });
        run_phase(entropy_scoring, [&](StreamedSampleState& state) {
            SpectralEntropyFeature::add(
                state.config, state.data, spectrum_predictions);
        });
        run_phase(rt_calibration, [&](StreamedSampleState& state) {
            PredictedRetentionTimeFeature::add(
                state.config, state.data, rt_predictions);
        });
        run_phase(fold_setup, [&](StreamedSampleState& state) {
            state.folds = SvmRescorer::assign_folds(state.data);
        });
        run_phase(svm_model, [&](StreamedSampleState& state) {
            RtResult rt;
            state.fitted = SvmRescorer::fit(
                state.data, rt.residuals, state.folds, config.train_fdr,
                config.max_iterations, config.svm_c_pos, config.svm_c_neg);
            state.scores = std::move(state.fitted.scores);
            state.folds.clear();
            state.folds.shrink_to_fit();
        });
        run_phase(primary_statistics, [&](StreamedSampleState& state) {
            std::vector<int> labels;
            labels.reserve(state.data.rows.size());
            for (const auto& psm : state.data.rows) {
                labels.push_back(psm.label);
            }
            auto& model = state.model;
            model.name = std::filesystem::path(
                config.output_prefixes[state.file]).filename().string();
            model.psms = state.data.rows.size();
            model.score_iterations = state.fitted.iterations.at(0);
            model.feature_weights =
                std::move(state.fitted.calibrated_weights.at(0));
            state.q = mixmax_qvalues(
                state.scores, labels, &model.pi0);
            state.pep = SvmRescorer::local_error_probabilities(
                state.scores, labels);
            state.compact_indices.assign(
                state.data.rows.size(),
                std::numeric_limits<std::size_t>::max());
            std::size_t compact = 0;
            for (std::size_t row = 0; row < state.data.rows.size(); ++row) {
                if (!(state.pep[row] < 0.5) &&
                    state.q[row] > config.q_threshold) {
                    continue;
                }
                state.compact_indices[row] = compact++;
            }
        });

        auto merge_begin = Clock::now();
        auto merge_cpu_begin = std::clock();
        std::size_t next_compact = data.rows.size();
        for (auto& state : batch) {
            auto& sample = state.data;
            summary.psms += sample.rows.size();
            summary.removed_decoy_peptide_collisions +=
                sample.removed_decoy_peptide_collisions;
            for (const auto& psm : sample.rows) {
                psm.label == 1 ? ++summary.targets : ++summary.decoys;
            }
            if (data.feature_names.empty()) {
                data.headers = sample.headers;
                data.columns = sample.columns;
                data.numeric_columns = sample.numeric_columns;
                data.feature_names = sample.feature_names;
                data.generated_feature_names = sample.generated_feature_names;
                data.has_predicted_rt_diagnostics =
                    sample.has_predicted_rt_diagnostics;
                data.spectrum_prediction_device =
                    sample.spectrum_prediction_device;
                data.rt_prediction_device = sample.rt_prediction_device;
            } else if (data.feature_names != sample.feature_names ||
                       data.headers != sample.headers) {
                throw std::runtime_error(
                    "Streamed sample feature schemas differ");
            }
            data.input_paths.push_back(sample.input_paths.front());
            std::size_t retained = 0;
            for (auto& compact : state.compact_indices) {
                if (compact == std::numeric_limits<std::size_t>::max()) {
                    continue;
                }
                compact += next_compact;
                ++retained;
            }
            next_compact += retained;
            summary.sample_models[state.file] = std::move(state.model);
            if (!config.filtered_only) {
                spools[state.file] = {
                    temporary.path /
                        (std::to_string(state.file) + ".target.tsv"),
                    temporary.path /
                        (std::to_string(state.file) + ".decoy.tsv")};
            }
        }
        add_timing(primary_statistics, merge_begin, merge_cpu_begin);

        if (!config.filtered_only) {
            run_phase(primary_statistics, [&](StreamedSampleState& state) {
                write_stream_score_spool(
                    spools[state.file].target, 1, state.data, state.scores,
                    state.q, state.pep, state.compact_indices);
                write_stream_score_spool(
                    spools[state.file].decoy, -1, state.data, state.scores,
                    state.q, state.pep, state.compact_indices);
            });
        }

        merge_begin = Clock::now();
        merge_cpu_begin = std::clock();
        data.rows.reserve(next_compact);
        scores.reserve(next_compact);
        q.reserve(next_compact);
        pep.reserve(next_compact);
        for (auto& state : batch) {
            for (std::size_t row = 0; row < state.data.rows.size(); ++row) {
                if (state.compact_indices[row] ==
                    std::numeric_limits<std::size_t>::max()) {
                    continue;
                }
                state.data.rows[row].file_id = state.file;
                data.rows.push_back(std::move(state.data.rows[row]));
                scores.push_back(state.scores[row]);
                q.push_back(state.q[row]);
                pep.push_back(state.pep[row]);
            }
        }
        add_timing(primary_statistics, merge_begin, merge_cpu_begin);
    }

    summary.input_paths = data.input_paths;
    summary.feature_names = data.feature_names;
    summary.read_timing = read_timing;
    summary.spectrum_entropy_timing = {
        spectrum_inference_outer.wall_seconds +
            entropy_scoring.wall_seconds,
        spectrum_inference_outer.cpu_seconds +
            entropy_scoring.cpu_seconds};
    summary.predicted_rt_timing = {
        rt_inference_outer.wall_seconds + rt_calibration.wall_seconds,
        rt_inference_outer.cpu_seconds + rt_calibration.cpu_seconds};
    summary.fold_setup_timing = fold_setup;
    summary.svm_model_timing = svm_model;
    summary.statistics_timing = primary_statistics;

    if (config.assemble_proteins && !config.database_path.empty() &&
        !config.decoy_database_path.empty()) {
        ProteinAssembler::sequential_filter(config, data, scores, pep, q);
    }

    std::unordered_set<std::string> peptides;
    std::unordered_set<std::string> peptide_forms;
    std::unordered_set<std::string> ptm_peptides;
    std::vector<std::unordered_set<std::string>> sample_peptides(
        summary.files);
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto& psm = data.rows[row];
        if (psm.label != 1 || q[row] > config.q_threshold) continue;
        auto& model = summary.sample_models[psm.file_id];
        ++model.target_ids;
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
    for (std::size_t file = 0; file < summary.files; ++file) {
        summary.sample_models[file].distinct_target_peptides =
            sample_peptides[file].size();
        summary.target_ids += summary.sample_models[file].target_ids;
    }
    summary.distinct_target_peptides = peptides.size();
    summary.distinct_target_peptide_forms = peptide_forms.size();
    summary.distinct_target_ptm_peptides = ptm_peptides.size();

    auto begin = Clock::now();
    auto cpu_begin = std::clock();
    auto quantification = ChromatographicQuantifier::add(config, data, q);
    summary.mbr_ions = data.transferred_ions.size();
    summary.quantification_stages = std::move(quantification.stages);
    add_timing(summary.quantification_timing, begin, cpu_begin);

    begin = Clock::now();
    cpu_begin = std::clock();
    Config filtered_config = config;
    filtered_config.filtered_only = true;
    RtResult no_internal_rt;
    ResultWriter::write(
        filtered_config, data, scores, q, pep,
        no_internal_rt, {});
    if (!config.filtered_only) {
        const bool protein_filtered = config.assemble_proteins &&
            !config.database_path.empty() &&
            !config.decoy_database_path.empty();
        for (std::size_t batch_begin = 0; batch_begin < summary.files;
             batch_begin += summary.sample_parallelism) {
            const std::size_t batch_count = std::min<std::size_t>(
                summary.sample_parallelism, summary.files - batch_begin);
            run_stream_sample_batch(
                batch_count, summary.threads, [&](std::size_t local) {
                    const auto file = batch_begin + local;
                    const auto& prefix = config.output_prefixes[file];
                    write_stream_score_table(
                        spools[file].target, prefix + "_target_psms.tsv",
                        q, protein_filtered);
                    write_stream_score_table(
                        spools[file].decoy, prefix + "_decoy_psms.tsv",
                        q, protein_filtered);
                });
        }
    }
    if (!config.sip_protein_output_path.empty()) {
        ProteinAssembler::write_sip_psm_mapping(config, data, q);
    }
    add_timing(summary.write_timing, begin, cpu_begin);

    begin = Clock::now();
    cpu_begin = std::clock();
    if (config.assemble_proteins && !config.database_path.empty() &&
        !config.decoy_database_path.empty()) {
        const auto assembly =
            ProteinAssembler::write(config, data, scores, q, pep);
        summary.protein_ids = assembly.proteins;
        summary.protein_output_dir = assembly.output_dir;
        summary.protein_assembly_stages = assembly.stages;
    }
    add_timing(summary.protein_assembly_timing, begin, cpu_begin);

    if (!config.negative_control_samples.empty()) {
        const auto negative = NegativeControlFilter::run(config, data, q);
        summary.negative_control_candidates = negative.candidates;
        summary.negative_control_input_psms = negative.input_psms;
        summary.negative_control_threshold_filtered_psms =
            negative.threshold_filtered_psms;
        summary.negative_control_targets = negative.targets;
        summary.negative_control_decoys = negative.decoys;
        summary.negative_control_target_ids = negative.target_ids;
        summary.negative_control_label_threshold = config.label_threshold;
        summary.negative_control_output_path = negative.output_path;
        summary.negative_control_target_output_path =
            negative.target_output_path;
        summary.negative_control_decoy_output_path =
            negative.decoy_output_path;
        summary.negative_control_protein_output_path =
            negative.protein_output_path;
        summary.negative_control_timing = negative.timing;
        summary.negative_control_stages = negative.stages;
        summary.negative_control_feature_names = negative.feature_names;
        summary.negative_control_model = negative.model;
    }

    summary.total_timing = {
        std::chrono::duration<double>(Clock::now() - total_begin).count(),
        static_cast<double>(std::clock() - total_cpu_begin) /
            CLOCKS_PER_SEC};
    if (summary.total_timing.wall_seconds > 0.0) {
        summary.omp_speedup_ratio = summary.total_timing.cpu_seconds /
            summary.total_timing.wall_seconds;
        summary.omp_parallel_efficiency = summary.omp_speedup_ratio /
            std::max(1u, summary.threads);
    }
    return summary;
}

} // namespace

Summary run(const Config& config) {
    const bool can_stream = config.stream_samples &&
        config.target_pins.size() > 1 &&
        config.target_pins.size() == config.decoy_pins.size() &&
        config.target_pins.size() == config.output_prefixes.size() &&
        config.target_pins.size() == config.spectrum_paths.size() &&
        config.predict_rt;
    return can_stream ? run_streamed(config) : run_monolithic(config);
}

} // namespace aerith
