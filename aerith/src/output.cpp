#include "pipeline.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <omp.h>

namespace aerith {

template <typename Number>
std::string compact_number(Number value) {
    char buffer[64];
    const auto result = std::to_chars(
        std::begin(buffer), std::end(buffer), value);
    if (result.ec != std::errc{}) {
        throw std::runtime_error("Cannot format stored PIN numeric value");
    }
    return std::string(buffer, result.ptr);
}

std::string ResultWriter::formatted_number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(10) << value;
    return stream.str();
}

std::string assigned_modifications_text(
    const std::vector<std::string>& modifications) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < modifications.size(); ++index) {
        if (index != 0) stream << ", ";
        stream << modifications[index];
    }
    return stream.str();
}

void write_ptm_annotation(
    std::ostream& stream, const Psm& psm, bool fixed_cam) {
    const auto modifications = modification_info(psm.peptide, fixed_cam);
    stream << '\t'
           << (modifications.modified_peptide.empty()
                   ? modifications.sequence
                   : modifications.modified_peptide)
           << '\t' << assigned_modifications_text(modifications.assigned);
}

void ResultWriter::write_original_field(
    std::ostream& stream, const Dataset& data, const Psm& psm,
    std::size_t column) {
    const auto& name = data.headers[column];
    if (name == "SpecId") stream << psm.id;
    else if (name == "Label") stream << psm.label;
    else if (name == "ScanNr") stream << psm.scan;
    else if (name == "retentiontime") {
        stream << compact_number(psm.retention);
    } else if (name == "ExpMass") {
        stream << compact_number(psm.exp_mass);
    } else if (name == "ObservedMass") {
        stream << compact_number(psm.observed_mass);
    }
    else if (name == "ranks") stream << psm.rank;
    else if (name == "diffScores") stream << formatted_number(psm.diff_score);
    else if (name == "MS1IsotopicAbundances") {
        stream << formatted_number(psm.ms1_isotopic_abundance);
    } else if (name == "MS2IsotopicAbundances") {
        stream << formatted_number(psm.ms2_isotopic_abundance);
    } else if (name == "SampleName") {
        stream << psm.sample_name;
    } else if (name == "Peptide") {
        stream << psm.peptide;
    } else if (name == "Proteins") {
        stream << psm.proteins;
    } else {
        const auto feature = std::find(
            data.numeric_columns.begin(),
            data.numeric_columns.end(), column);
        if (feature == data.numeric_columns.end()) {
            throw std::runtime_error(
                "Cannot reconstruct stored PIN column: " + name);
        }
        const auto index = static_cast<std::size_t>(
            feature - data.numeric_columns.begin());
        if (index >= psm.features.size()) {
            throw std::runtime_error(
                "Stored PSM has no value for PIN column: " + name);
        }
        stream << compact_number(psm.features[index]);
    }
}

std::vector<std::size_t> ResultWriter::selected_rows(
    const Dataset& data, std::size_t file, int label) {
    std::vector<std::size_t> rows;
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        if ((file == std::numeric_limits<std::size_t>::max() ||
             data.rows[i].file_id == file) &&
            (label == 0 || data.rows[i].label == label)) {
            rows.push_back(i);
        }
    }
    return rows;
}

void ResultWriter::write_results(
    const std::string& path, int label, std::size_t file, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep, bool fixed_cam) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }
    std::ofstream stream(output);
    if (!stream) {
        throw std::runtime_error("Cannot create output: " + path);
    }
    stream << "PSMId\tSVMscore\tq-value\tposterior_error_prob"
           << "\tpeptide\tmodifiedPeptide\tassignedModifications"
           << "\tproteinIds\n";
    auto order = selected_rows(data, file, label);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b];
    });
    stream << std::setprecision(10);
    for (const auto i : order) {
        const auto& psm = data.rows[i];
        stream << psm.id << '\t' << scores[i] << '\t' << q[i] << '\t'
               << pep[i] << '\t' << psm.peptide;
        write_ptm_annotation(stream, psm, fixed_cam);
        stream << '\t' << psm.proteins << '\n';
    }
}

void ResultWriter::write_filtered_results(
    const std::string& path, std::size_t file, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep, const std::vector<double>& rt_residual,
    double threshold, bool sip_output, bool fixed_cam) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("Cannot create output: " + path);
    stream << "PSMId\tSVMscore\tq-value\tposterior_error_prob";
    for (const auto& header : data.headers) {
        if (header == "SpecId" ||
            (sip_output && (header == "Label" ||
                            header == "diffScores"))) {
            continue;
        }
        if (header == "Peptide") {
            for (const auto& generated : data.generated_feature_names) {
                if (generated == "delta_RT_loess") continue;
                stream << '\t' << generated;
            }
        }
        stream << '\t' << header;
        if (header == "Peptide") {
            stream << "\tModifiedPeptide\tAssignedModifications";
        }
        if (header == "retentiontime" && data.has_predicted_rt_diagnostics) {
            stream << "\tpred_RT_real_units\tdelta_RT_loess_real";
        }
        if (header == "parentCharges" && !rt_residual.empty()) {
            stream << "\tsqrtAbsDeltaRT";
        }
    }
    stream << '\n' << std::setprecision(10);
    auto rows = selected_rows(data, file, sip_output ? 1 : 0);
    rows.erase(std::remove_if(rows.begin(), rows.end(), [&](std::size_t i) {
        return q[i] > threshold;
    }), rows.end());
    std::stable_sort(rows.begin(), rows.end(), [&](std::size_t a, std::size_t b) {
        if (data.rows[a].scan != data.rows[b].scan) return data.rows[a].scan < data.rows[b].scan;
        return data.rows[a].rank < data.rows[b].rank;
    });
    for (const auto i : rows) {
        const auto& psm = data.rows[i];
        stream << psm.id << '\t' << scores[i] << '\t' << q[i] << '\t' << pep[i];
        for (std::size_t column = 0; column < data.headers.size(); ++column) {
            if (data.headers[column] == "SpecId" ||
                (sip_output && (data.headers[column] == "Label" ||
                                data.headers[column] == "diffScores"))) {
                continue;
            }
            if (data.headers[column] == "Peptide") {
                const std::size_t first_generated =
                    data.feature_names.size() - data.generated_feature_names.size();
                for (std::size_t generated = 0;
                     generated < data.generated_feature_names.size(); ++generated) {
                    if (data.generated_feature_names[generated] ==
                        "delta_RT_loess") {
                        continue;
                    }
                    stream << '\t' << formatted_number(
                        psm.features[first_generated + generated]);
                }
            }
            stream << '\t';
            write_original_field(stream, data, psm, column);
            if (data.headers[column] == "Peptide") {
                write_ptm_annotation(stream, psm, fixed_cam);
            }
            if (data.headers[column] == "retentiontime" &&
                data.has_predicted_rt_diagnostics) {
                stream << '\t' << formatted_number(psm.predicted_rt_real_units)
                       << '\t' << formatted_number(psm.delta_rt_loess_real);
            }
            if (data.headers[column] == "parentCharges" && !rt_residual.empty()) {
                stream << '\t' << rt_residual[i];
            }
        }
        stream << '\n';
    }
}

void ResultWriter::write(
    const Config& config, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep, const RtResult& rt,
    const std::vector<std::size_t>& outer_folds) {
    const bool has_internal_rt = std::all_of(
        rt.residuals.begin(), rt.residuals.end(),
        [&](const auto& residuals) {
            return residuals.size() == data.rows.size();
        });
    const bool has_partial_internal_rt = std::any_of(
        rt.residuals.begin(), rt.residuals.end(),
        [](const auto& residuals) { return !residuals.empty(); });
    if (!has_internal_rt && has_partial_internal_rt) {
        throw std::runtime_error("Incomplete internal RT residuals");
    }
    std::vector<double> held_out_rt_residuals;
    if (has_internal_rt) {
        held_out_rt_residuals.resize(data.rows.size());
        #pragma omp parallel for simd schedule(static)
        for (std::ptrdiff_t row = 0;
             row < static_cast<std::ptrdiff_t>(data.rows.size()); ++row) {
            const auto i = static_cast<std::size_t>(row);
            held_out_rt_residuals[i] = rt.residuals[outer_folds[i]][i];
        }
    }
    const bool aggregate =
        config.output_prefixes.size() == 1 && data.input_paths.size() > 1;
    const std::size_t output_count = aggregate ? 1 : data.input_paths.size();
    std::exception_ptr failure;
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t output_index = 0;
         output_index < static_cast<std::ptrdiff_t>(output_count); ++output_index) {
        try {
            const auto file = aggregate ? std::numeric_limits<std::size_t>::max()
                                        : static_cast<std::size_t>(output_index);
            const auto& prefix =
                config.output_prefixes[static_cast<std::size_t>(output_index)];
            if (!config.filtered_only) {
                write_results(prefix + "_target_psms.tsv", 1, file,
                              data, scores, q, pep, config.fixed_cam);
                write_results(prefix + "_decoy_psms.tsv", -1, file,
                              data, scores, q, pep, config.fixed_cam);
            }
            write_filtered_results(prefix + "_filtered_psms.tsv", file, data,
                                   scores, q, pep, held_out_rt_residuals,
                                   config.q_threshold,
                                   !config.sip_isotope.empty(),
                                   config.fixed_cam);
        } catch (...) {
            #pragma omp critical(aerith_write_failure)
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
}

void print_summary(std::ostream& output, const Summary& summary) {
    const auto old_flags = output.flags();
    const auto old_precision = output.precision();
    const auto old_fill = output.fill();
    const auto print_summary_value = [&](
        const std::string& label, const auto& value,
        const char* suffix = "") {
        output << "  " << std::left << std::setw(32) << label
               << std::right << value << suffix << '\n';
    };
    output << "Aerith SVM+RT report\n"
           << "========================================================================\n"
           << "Inputs\n";
    for (const auto& path : summary.input_paths) output << "  " << path << '\n';
    output << "\nDataset\n";
    print_summary_value("Files", summary.files);
    print_summary_value("Input PSMs", summary.psms);
    print_summary_value("Targets", summary.targets);
    print_summary_value("Decoys", summary.decoys);
    print_summary_value(
        "Removed colliding decoy PSMs",
        summary.removed_decoy_peptide_collisions);
    print_summary_value("OpenMP threads", summary.threads);
    print_summary_value("Concurrent samples", summary.sample_parallelism);
    print_summary_value("Score model", summary.score_model);

    output << "\nResults\n";
    print_summary_value(
        "Reporting FDR", summary.reporting_fdr * 100.0, "%");
    print_summary_value("Target PSMs", summary.target_ids);
    print_summary_value(
        "Distinct naked peptides", summary.distinct_target_peptides);
    print_summary_value(
        "Distinct peptide forms", summary.distinct_target_peptide_forms);
    print_summary_value(
        "Distinct PTM peptide forms", summary.distinct_target_ptm_peptides);
    print_summary_value("PSMs carrying PTMs", summary.target_ptm_psms);
    if (summary.mbr_ions > 0) {
        print_summary_value("MBR transferred ions", summary.mbr_ions);
    }
    if (summary.protein_ids > 0) {
        print_summary_value("Protein IDs at 1% FDR", summary.protein_ids);
    }
    if (summary.negative_control_candidates > 0) {
        print_summary_value(
            "SIP-Negative-control threshold",
            summary.negative_control_label_threshold, "%");
        print_summary_value(
            "Primary-filtered target PSMs",
            summary.negative_control_input_psms);
        print_summary_value(
            "PSMs below SIP threshold",
            summary.negative_control_threshold_filtered_psms);
        print_summary_value(
            "SIP-Negative-control candidates",
            summary.negative_control_candidates);
        print_summary_value(
            "SIP-Negative-control targets",
            summary.negative_control_targets);
        print_summary_value(
            "SIP-Negative-control decoys",
            summary.negative_control_decoys);
        print_summary_value(
            "SIP-Negative-control target IDs",
            summary.negative_control_target_ids);
    }
    output << std::fixed << std::setprecision(6);
    if (summary.used_internal_rt_model) {
        print_summary_value(
            "RT training targets/fold", summary.rt_training_targets);
        print_summary_value("Held-out RT R2", summary.rt_r2);
    } else {
        print_summary_value("RT feature source", "DIA-NN prediction");
        print_summary_value("Aerith internal RT model", "skipped");
    }
    output << "\nSample-specific SVM results\n";
    output << std::left << std::setw(30) << "Sample"
           << std::right << std::setw(14) << "Input PSMs"
           << std::setw(12) << "IDs"
           << std::setw(12) << "Peptides"
           << std::setw(12) << "pi0"
           << std::setw(15) << "Iterations" << '\n'
           << std::string(95, '-') << '\n';
    for (const auto& sample : summary.sample_models) {
        const std::string iterations = sample.model_valid
            ? std::to_string(sample.score_iterations[0]) + "/" +
                std::to_string(sample.score_iterations[1]) + "/" +
                std::to_string(sample.score_iterations[2])
            : "skipped";
        output << std::left << std::setw(30) << sample.name
               << std::right << std::setw(14) << sample.psms
               << std::setw(12) << sample.target_ids
               << std::setw(12) << sample.distinct_target_peptides
               << std::setw(12) << sample.pi0
               << std::setw(15) << iterations << '\n';
        if (!sample.model_valid) {
            output << "  WARNING: " << sample.name
                   << " accepted 0 PSMs because SVM initialization failed: "
                   << sample.warning << '\n';
        }
    }
    output << "\nTiming by stage (seconds)\n"
           << std::left << std::setw(64) << "Stage"
           << std::right << std::setw(14) << "Wall time"
           << std::setw(14) << "CPU time"
           << std::setw(12) << "Speedup" << '\n'
           << std::string(104, '-') << '\n';
    const auto print_wrapped_timing_detail = [&](const std::string& detail) {
        constexpr std::size_t line_width = 104;
        const std::string first_prefix = "      Detail: ";
        const std::string continuation(first_prefix.size(), ' ');
        std::size_t begin = 0;
        bool first = true;
        while (begin < detail.size()) {
            const std::string& prefix = first ? first_prefix : continuation;
            const std::size_t available = line_width - prefix.size();
            std::size_t end = std::min(detail.size(), begin + available);
            if (end < detail.size()) {
                const std::size_t space = detail.rfind(' ', end);
                if (space != std::string::npos && space > begin) end = space;
            }
            output << prefix << detail.substr(begin, end - begin) << '\n';
            begin = end;
            while (begin < detail.size() && detail[begin] == ' ') ++begin;
            first = false;
        }
    };
    const auto print_timing = [&](const std::string& label,
                                  const StageTiming& timing) {
        const double speedup = timing.wall_seconds > 0.0
            ? timing.cpu_seconds / timing.wall_seconds : 0.0;
        const bool overlong = label.size() > 64;
        const std::string display = overlong
            ? label.substr(0, 61) + "..." : label;
        output << std::left << std::setw(64) << display
               << std::right << std::setw(14) << timing.wall_seconds
               << std::setw(14) << timing.cpu_seconds
               << std::setw(11) << speedup << "x\n";
        if (overlong) print_wrapped_timing_detail("Full stage: " + label);
    };
    const auto displayed_seconds = [](double seconds) {
        constexpr double scale = 1000000.0;
        return std::round(seconds * scale) / scale;
    };
    StageTiming accounted_timing;
    const auto print_top_level_timing = [&](const std::string& label,
                                            const StageTiming& timing) {
        print_timing(label, timing);
        accounted_timing.wall_seconds += displayed_seconds(timing.wall_seconds);
        accounted_timing.cpu_seconds += displayed_seconds(timing.cpu_seconds);
    };
    const bool has_spectrum_entropy =
        std::find(summary.feature_names.begin(), summary.feature_names.end(),
                  "unweighted_spectral_entropy") != summary.feature_names.end();
    const bool has_predicted_rt =
        std::find(summary.feature_names.begin(), summary.feature_names.end(),
                  "delta_RT_loess") != summary.feature_names.end();
    print_top_level_timing(
        "Read, merge, and rerank PIN files", summary.read_timing);
    if (has_spectrum_entropy) {
        print_top_level_timing(
            "Predict spectra and compute entropy",
            summary.spectrum_entropy_timing);
        if (!summary.spectrum_prediction_device.empty()) {
            print_timing(
                "  DIA-NN spectra prediction (" +
                    summary.spectrum_prediction_device + ")",
                summary.spectrum_prediction_timing);
        }
        print_timing(
            "  Spectrum prediction cache .bin file read",
            summary.spectrum_cache_read_timing);
        print_timing(
            "  Spectrum prediction cache .bin file merge/write",
            summary.spectrum_cache_write_timing);
    }
    if (has_predicted_rt) {
        print_top_level_timing(
            "Predict RT and compute delta-RT",
            summary.predicted_rt_timing);
        if (!summary.rt_prediction_device.empty()) {
            print_timing(
                "  DIA-NN RT prediction (" + summary.rt_prediction_device + ")",
                summary.rt_prediction_inference_timing);
        }
        print_timing(
            "  RT prediction cache .bin file read",
            summary.rt_cache_read_timing);
        print_timing(
            "  RT prediction cache .bin file merge/write",
            summary.rt_cache_write_timing);
    }
    print_top_level_timing(
        "Assign folds and seed q-values", summary.fold_setup_timing);
    if (summary.used_internal_rt_model) {
        print_top_level_timing("Fit nested RT models", summary.rt_model_timing);
    }
    print_top_level_timing("Fit and score SVM folds", summary.svm_model_timing);
    print_top_level_timing(
        "Compute q-values and PEPs", summary.statistics_timing);
    print_top_level_timing("Quantification total", summary.quantification_timing);
    for (const auto& stage : summary.quantification_stages) {
        print_timing("  " + stage.name, stage.timing);
    }
    print_top_level_timing("Write result files", summary.write_timing);
    print_top_level_timing(
        "Protein assembly total", summary.protein_assembly_timing);
    for (const auto& stage : summary.protein_assembly_stages) {
        print_timing("  " + stage.name, stage.timing);
    }
    if (summary.negative_control_candidates > 0) {
        print_top_level_timing(
            "SIP-Negative-control filtering total",
            summary.negative_control_timing);
        for (const auto& stage : summary.negative_control_stages) {
            print_timing("  " + stage.name, stage.timing);
        }
    }
    const StageTiming coordination_timing{
        displayed_seconds(summary.total_timing.wall_seconds) -
            accounted_timing.wall_seconds,
        displayed_seconds(summary.total_timing.cpu_seconds) -
            accounted_timing.cpu_seconds};
    print_timing(
        "Workflow coordination and timing overhead", coordination_timing);
    print_timing("Total", summary.total_timing);
    output << "  Overall OpenMP efficiency: "
           << summary.omp_parallel_efficiency * 100.0 << "%\n";
    bool wrote_calibration_header = false;
    for (const auto& stage : summary.quantification_stages) {
        if (stage.detail.empty()) continue;
        if (!wrote_calibration_header) {
            output << "\nMBR calibration audit\n";
            wrote_calibration_header = true;
        }
        output << "  " << stage.name << '\n';
        constexpr std::size_t line_width = 104;
        const std::string prefix = "    ";
        std::size_t begin = 0;
        while (begin < stage.detail.size()) {
            const std::size_t available = line_width - prefix.size();
            std::size_t end = std::min(
                stage.detail.size(), begin + available);
            if (end < stage.detail.size()) {
                const std::size_t space = stage.detail.rfind(' ', end);
                if (space != std::string::npos && space > begin) end = space;
            }
            output << prefix
                   << stage.detail.substr(begin, end - begin) << '\n';
            begin = end;
            while (begin < stage.detail.size() &&
                   stage.detail[begin] == ' ') ++begin;
        }
    }
    output << "\nSVM feature weights\n"
           << "  Positive raises the target score; negative lowers it.\n"
           << "  Values are calibrated weights on each fold's standardized scale.\n";
    const auto print_feature_weights = [&](
        const std::vector<std::string>& feature_names,
        const SampleModelSummary& sample,
        bool print_sample_name = true) {
        const auto mean_weight = [&](std::size_t index) {
            return (sample.feature_weights[0][index] +
                    sample.feature_weights[1][index] +
                    sample.feature_weights[2][index]) / 3.0;
        };
        std::vector<std::size_t> order(feature_names.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return std::abs(mean_weight(a)) > std::abs(mean_weight(b));
        });
        output << '\n';
        if (print_sample_name) {
            output << "  Sample: " << sample.name << '\n';
        }
        output << std::left << std::setw(36) << "Feature"
               << std::right << std::setw(13) << "Mean"
               << std::setw(13) << "Fold 1"
               << std::setw(13) << "Fold 2" << std::setw(13) << "Fold 3" << '\n'
               << std::string(88, '-') << '\n'
               << std::fixed << std::showpos << std::setprecision(6);
        const auto print_weight = [&](const std::string& name, std::size_t index) {
            output << std::left << std::setw(36) << name << std::right
                   << std::setw(13) << mean_weight(index)
                   << std::setw(13) << sample.feature_weights[0][index]
                   << std::setw(13) << sample.feature_weights[1][index]
                   << std::setw(13) << sample.feature_weights[2][index] << '\n';
        };
        for (const auto index : order) {
            print_weight(feature_names[index], index);
        }
        output << std::noshowpos;
        print_weight("(intercept)", feature_names.size());
    };
    for (const auto& sample : summary.sample_models) {
        if (sample.model_valid) {
            print_feature_weights(summary.feature_names, sample);
        }
    }
    if (summary.negative_control_candidates > 0 &&
        summary.negative_control_model.model_valid) {
        output
            << "\nSIP-Negative-control SVM feature weights\n"
            << "  These weights are fitted during the native in-memory "
               "secondary control pass.\n";
        print_feature_weights(
            summary.negative_control_feature_names,
            summary.negative_control_model,
            false);
    } else if (summary.negative_control_candidates > 0) {
        output << "\nSIP-Negative-control SVM warning\n  "
               << summary.negative_control_model.warning
               << "; 0 PSMs accepted.\n";
    }
    output << "\nOutputs\n";
    for (const auto& prefix : summary.output_prefixes) {
        if (!summary.filtered_only) {
            output << "  " << prefix << "_target_psms.tsv\n"
                   << "  " << prefix << "_decoy_psms.tsv\n";
        }
        output << "  " << prefix << "_filtered_psms.tsv\n";
    }
    if (!summary.protein_output_dir.empty()) {
        output << "  " << summary.protein_output_dir << "/combined_psm.tsv\n"
               << "  " << summary.protein_output_dir << "/combined_ion.tsv\n"
               << "  " << summary.protein_output_dir
               << "/combined_modified_peptide.tsv\n"
               << "  " << summary.protein_output_dir
               << "/combined_peptide.tsv\n"
               << "  " << summary.protein_output_dir
               << "/combined_protein.tsv\n"
               << "  " << summary.protein_output_dir
               << "/combined_protein.fas\n";
    }
    if (!summary.negative_control_target_output_path.empty()) {
        output << "  " << summary.negative_control_target_output_path << '\n';
    }
    if (!summary.negative_control_decoy_output_path.empty()) {
        output << "  " << summary.negative_control_decoy_output_path << '\n';
    }
    if (!summary.negative_control_output_path.empty()) {
        output << "  " << summary.negative_control_output_path << '\n';
    }
    if (!summary.negative_control_protein_output_path.empty()) {
        output << "  " << summary.negative_control_protein_output_path << '\n';
    }
    output.flags(old_flags);
    output.precision(old_precision);
    output.fill(old_fill);
}

} // namespace aerith
