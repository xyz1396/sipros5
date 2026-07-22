#include "pipeline.hpp"

#include <algorithm>
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

std::vector<std::string_view> ResultWriter::split_tabs(const std::string& line) {
    std::vector<std::string_view> fields;
    std::size_t begin = 0;
    while (true) {
        const auto tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.emplace_back(line.data() + begin, line.size() - begin);
            break;
        }
        fields.emplace_back(line.data() + begin, tab - begin);
        begin = tab + 1;
    }
    return fields;
}

std::size_t ResultWriter::required_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("PIN file is missing required column: " + name);
    }
    return found->second;
}

std::string ResultWriter::formatted_number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(10) << value;
    return stream.str();
}

std::vector<std::string_view> ResultWriter::row_fields(
    const Dataset& data, const Psm& psm) {
    auto fields = split_tabs(psm.raw_line);
    if (fields.size() != data.headers.size()) {
        throw std::runtime_error("Stored PIN row no longer matches its header");
    }
    return fields;
}

void ResultWriter::write_original_field(
    std::ostream& stream, const Dataset& data, const Psm& psm,
    const std::vector<std::string_view>& fields, std::size_t column) {
    const auto& name = data.headers[column];
    if (name == "SpecId") stream << psm.id;
    else if (name == "ranks") stream << psm.rank;
    else if (name == "diffScores") stream << formatted_number(psm.diff_score);
    else stream << fields[column];
}

std::string ResultWriter::original_field(
    const Dataset& data, const Psm& psm, const std::string& name) {
    const auto column = required_column(data.columns, name);
    if (name == "SpecId") return psm.id;
    if (name == "ranks") return std::to_string(psm.rank);
    if (name == "diffScores") return formatted_number(psm.diff_score);
    const auto fields = row_fields(data, psm);
    return std::string(fields[column]);
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
    const std::vector<double>& pep) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }
    std::ofstream stream(output);
    if (!stream) {
        throw std::runtime_error("Cannot create output: " + path);
    }
    stream << "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";
    auto order = selected_rows(data, file, label);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b];
    });
    stream << std::setprecision(10);
    for (const auto i : order) {
        const auto& psm = data.rows[i];
        stream << psm.id << '\t' << scores[i] << '\t' << q[i] << '\t' << pep[i] << '\t'
               << psm.peptide << '\t' << psm.proteins << '\n';
    }
}

void ResultWriter::write_filtered_results(
    const std::string& path, std::size_t file, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep, const std::vector<double>& rt_residual,
    double threshold) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("Cannot create output: " + path);
    stream << "PSMId\tscore\tq-value\tposterior_error_prob";
    for (const auto& header : data.headers) {
        if (header == "SpecId") continue;
        if (header == "Peptide") {
            for (const auto& generated : data.generated_feature_names) {
                if (generated == "delta_RT_loess") continue;
                stream << '\t' << generated;
            }
        }
        stream << '\t' << header;
        if (header == "retentiontime" && data.has_predicted_rt_diagnostics) {
            stream << "\tpred_RT_real_units\tdelta_RT_loess_real";
        }
        if (header == "parentCharges" && !rt_residual.empty()) {
            stream << "\tsqrtAbsDeltaRT";
        }
    }
    stream << '\n' << std::setprecision(10);
    auto rows = selected_rows(data, file);
    rows.erase(std::remove_if(rows.begin(), rows.end(), [&](std::size_t i) {
        return q[i] > threshold;
    }), rows.end());
    std::stable_sort(rows.begin(), rows.end(), [&](std::size_t a, std::size_t b) {
        if (data.rows[a].scan != data.rows[b].scan) return data.rows[a].scan < data.rows[b].scan;
        return data.rows[a].rank < data.rows[b].rank;
    });
    for (const auto i : rows) {
        const auto& psm = data.rows[i];
        const auto fields = row_fields(data, psm);
        stream << psm.id << '\t' << scores[i] << '\t' << q[i] << '\t' << pep[i];
        for (std::size_t column = 0; column < data.headers.size(); ++column) {
            if (data.headers[column] == "SpecId") continue;
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
            write_original_field(stream, data, psm, fields, column);
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

std::string ResultWriter::xml_escape(std::string_view value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char ch : value) {
        switch (ch) {
        case '&': escaped += "&amp;"; break;
        case '<': escaped += "&lt;"; break;
        case '>': escaped += "&gt;"; break;
        case '\"': escaped += "&quot;"; break;
        case '\'': escaped += "&apos;"; break;
        default: escaped += ch; break;
        }
    }
    return escaped;
}

std::string ResultWriter::trim(std::string value) {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

std::string ResultWriter::sanitized_protein(
    std::string protein, const std::string& decoy_prefix) {
    protein = trim(std::move(protein));
    std::string active_prefix;
    for (const auto& prefix : {decoy_prefix, std::string("DECOY_"), std::string("Decoy_")}) {
        if (!prefix.empty() && protein.rfind(prefix, 0) == 0) {
            active_prefix = decoy_prefix;
            protein.erase(0, prefix.size());
            break;
        }
    }
    if (protein.find('|') == std::string::npos) {
        std::replace_if(protein.begin(), protein.end(), [](char ch) {
            return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n';
        }, '_');
        protein = "sp|" + protein + "|" + protein;
    }
    return active_prefix + protein;
}

std::vector<std::string> ResultWriter::protein_ids(
    const Psm& psm, const std::string& decoy_prefix) {
    std::string proteins = psm.proteins;
    if (!proteins.empty() && proteins.front() == '{') proteins.erase(proteins.begin());
    if (!proteins.empty() && proteins.back() == '}') proteins.pop_back();
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= proteins.size()) {
        const auto comma = proteins.find(',', begin);
        const auto end = comma == std::string::npos ? proteins.size() : comma;
        auto protein = trim(proteins.substr(begin, end - begin));
        if (!protein.empty()) result.push_back(sanitized_protein(std::move(protein), decoy_prefix));
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    return result;
}

void ResultWriter::write_pepxml(
    const std::string& path, std::size_t file, const Config& config,
    const Dataset& data, const std::vector<double>& scores,
    const std::vector<double>& pep) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("Cannot create output: " + path);
    std::string database = config.database_path;
    if (database.empty()) {
        auto root = output.parent_path().parent_path();
        database = (root / "targetDecoy.faa").string();
    }
    database = std::filesystem::absolute(database).lexically_normal().string();
    auto rows = selected_rows(data, file);
    rows.erase(std::remove_if(rows.begin(), rows.end(), [&](std::size_t i) {
        return pep[i] >= 0.5 || (config.decoy_prefix == "DECOY_" && data.rows[i].label == -1);
    }), rows.end());
    std::stable_sort(rows.begin(), rows.end(), [&](std::size_t a, std::size_t b) {
        if (data.rows[a].scan != data.rows[b].scan) return data.rows[a].scan < data.rows[b].scan;
        return data.rows[a].rank < data.rows[b].rank;
    });
    stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           << "<msms_pipeline_analysis>\n  <analysis_summary/>\n  <msms_run_summary>\n"
           << "    <search_summary precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_engine=\"X! Tandem\" search_engine_version=\"Sipros\">\n"
           << "      <search_database local_path=\"" << xml_escape(database) << "\" type=\"AA\"/>\n"
           << "      <terminal_modification terminus=\"N\" massdiff=\"42.010565\" mass=\"43.018390\" variable=\"Y\" protein_terminus=\"Y\" description=\"Protein N-terminal acetylation\"/>\n"
           << "      <parameter name=\"database_name\" value=\"" << xml_escape(database) << "\"/>\n"
           << "    </search_summary>\n";
    std::size_t index = 0;
    for (const auto i : rows) {
        const auto& psm = data.rows[i];
        const auto proteins = protein_ids(psm, config.decoy_prefix);
        if (proteins.empty()) continue;
        const auto open = psm.peptide.find('[');
        const auto close = psm.peptide.find(']', open == std::string::npos ? 0 : open + 1);
        std::string previous = open == std::string::npos ? "-" : psm.peptide.substr(0, open);
        std::string next = close == std::string::npos ? "-" : psm.peptide.substr(close + 1);
        if (previous.empty()) previous = "-";
        if (next.empty()) next = "-";
        const bool nterm_acetyl =
            open != std::string::npos && open + 1 < psm.peptide.size() &&
            psm.peptide[open + 1] == '%';
        const std::string sequence = stripped_peptide(psm.peptide);
        const std::string probability = formatted_number(1.0 - pep[i]);
        stream << "    <spectrum_query start_scan=\"" << psm.scan
               << "\" assumed_charge=\"" << xml_escape(original_field(data, psm, "parentCharges"))
               << "\" spectrum=\"" << xml_escape(psm.id) << "\" end_scan=\"" << psm.scan
               << "\" index=\"" << index++ << "\" precursor_neutral_mass=\"" << psm.exp_mass
               << "\" retention_time_sec=\"" << psm.retention << "\">\n"
               << "      <search_result>\n        <search_hit peptide=\"" << xml_escape(sequence)
               << "\" massdiff=\"" << xml_escape(original_field(data, psm, "massErrors"))
               << "\" calc_neutral_pep_mass=\"" << psm.exp_mass
               << "\" num_missed_cleavages=\"" << xml_escape(original_field(data, psm, "missCleavageSiteNumbers"))
               << "\" num_tol_term=\"2\" protein_descr=\"" << xml_escape(proteins.front())
               << "\" num_tot_proteins=\"" << proteins.size() << "\" hit_rank=\"" << psm.rank
               << "\" protein=\"" << xml_escape(proteins.front())
               << "\" peptide_prev_aa=\"" << xml_escape(previous)
               << "\" peptide_next_aa=\"" << xml_escape(next) << "\" is_rejected=\"0\">\n";
        for (std::size_t protein = 1; protein < proteins.size(); ++protein) {
            stream << "          <alternative_protein protein_descr=\"" << xml_escape(proteins[protein])
                   << "\" protein=\"" << xml_escape(proteins[protein])
                   << "\" peptide_prev_aa=\"" << xml_escape(previous)
                   << "\" peptide_next_aa=\"" << xml_escape(next)
                   << "\" num_tol_term=\"2\"/>\n";
        }
        if (nterm_acetyl) {
            stream << "          <modification_info modified_peptide=\"n[43]" << xml_escape(sequence)
                   << "\" mod_nterm_mass=\"43.018390\"/>\n";
        }
        stream << "          <analysis_result analysis=\"peptideprophet\">\n"
               << "            <peptideprophet_result probability=\"" << probability
               << "\" all_ntt_prob=\"(" << probability << ',' << probability << ',' << probability
               << ")\"/>\n          </analysis_result>\n"
               << "          <hyperscore value=\"" << scores[i] << "\"/>\n"
               << "          <nextscore value=\"0\"/>\n          <expect value=\"0\"/>\n"
               << "        </search_hit>\n      </search_result>\n    </spectrum_query>\n";
    }
    stream << "  </msms_run_summary>\n</msms_pipeline_analysis>\n";
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
            write_results(prefix + "_target_psms.tsv", 1, file,
                          data, scores, q, pep);
            write_results(prefix + "_decoy_psms.tsv", -1, file,
                          data, scores, q, pep);
            write_filtered_results(prefix + "_filtered_psms.tsv", file, data,
                                   scores, q, pep, held_out_rt_residuals,
                                   config.q_threshold);
            write_pepxml(prefix + ".pep.xml", file, config, data, scores, pep);
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
    output << "Aerith SVM+RT report\n"
           << "========================================================================\n"
           << "Inputs\n";
    for (const auto& path : summary.input_paths) output << "  " << path << '\n';
    output << "\nDataset\n"
           << "  Files                         " << summary.files << '\n'
           << "  Input PSMs                    " << summary.psms << '\n'
           << "  Targets                       " << summary.targets << '\n'
           << "  Decoys                        " << summary.decoys << '\n'
           << "  OpenMP threads                " << summary.threads << '\n'
           << "  Score model                   " << summary.score_model << "\n\n"
           << "Results\n"
           << "  Target PSMs at threshold      " << summary.target_ids << '\n'
           << "  Distinct target peptides      " << summary.distinct_target_peptides << '\n';
    output << std::fixed << std::setprecision(6);
    if (summary.used_internal_rt_model) {
        output << "  RT training targets/fold      " << summary.rt_training_targets << '\n'
               << "  Held-out RT R2                " << summary.rt_r2 << '\n';
    } else {
        output << "  RT feature source             DIA-NN prediction\n"
               << "  Aerith internal RT model      skipped\n";
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
        const std::string iterations = std::to_string(sample.score_iterations[0]) + "/" +
            std::to_string(sample.score_iterations[1]) + "/" +
            std::to_string(sample.score_iterations[2]);
        output << std::left << std::setw(30) << sample.name
               << std::right << std::setw(14) << sample.psms
               << std::setw(12) << sample.target_ids
               << std::setw(12) << sample.distinct_target_peptides
               << std::setw(12) << sample.pi0
               << std::setw(15) << iterations << '\n';
    }
    output << "\nTiming (seconds)\n"
           << std::left << std::setw(38) << "Stage"
           << std::right << std::setw(13) << "Wall time"
           << std::setw(13) << "CPU time"
           << std::setw(11) << "Speedup" << '\n'
           << std::string(75, '-') << '\n';
    const auto print_timing = [&](const std::string& label,
                                  const StageTiming& timing) {
        const double speedup = timing.wall_seconds > 0.0
            ? timing.cpu_seconds / timing.wall_seconds : 0.0;
        output << std::left << std::setw(38) << label
               << std::right << std::setw(13) << timing.wall_seconds
               << std::setw(13) << timing.cpu_seconds
               << std::setw(10) << speedup << "x\n";
    };
    print_timing("Read, merge, and rerank PIN files", summary.read_timing);
    if (std::find(summary.feature_names.begin(), summary.feature_names.end(),
                  "unweighted_spectral_entropy") != summary.feature_names.end()) {
        print_timing("Predict spectra and compute entropy",
                     summary.spectrum_entropy_timing);
        if (!summary.spectrum_prediction_device.empty()) {
            print_timing(
                "  DIA-NN spectra prediction (" +
                    summary.spectrum_prediction_device + ")",
                summary.spectrum_prediction_timing);
        }
    }
    if (std::find(summary.feature_names.begin(), summary.feature_names.end(),
                  "delta_RT_loess") != summary.feature_names.end()) {
        print_timing("Predict RT and compute delta-RT",
                     summary.predicted_rt_timing);
        if (!summary.rt_prediction_device.empty()) {
            print_timing(
                "  DIA-NN RT prediction (" + summary.rt_prediction_device + ")",
                summary.rt_prediction_inference_timing);
        }
    }
    print_timing("Assign folds and seed q-values", summary.fold_setup_timing);
    if (summary.used_internal_rt_model) {
        print_timing("Fit nested RT models", summary.rt_model_timing);
    }
    print_timing("Fit and score SVM folds", summary.svm_model_timing);
    print_timing("Compute q-values and PEPs", summary.statistics_timing);
    print_timing("Write result files", summary.write_timing);
    print_timing("Total", summary.total_timing);
    output << "  Overall OpenMP efficiency: "
           << summary.omp_parallel_efficiency * 100.0 << "%\n";
    output << "\nSVM feature weights\n"
           << "  Positive raises the target score; negative lowers it.\n"
           << "  Values are calibrated weights on each fold's standardized scale.\n";
    for (const auto& sample : summary.sample_models) {
        const auto mean_weight = [&](std::size_t index) {
            return (sample.feature_weights[0][index] +
                    sample.feature_weights[1][index] +
                    sample.feature_weights[2][index]) / 3.0;
        };
        std::vector<std::size_t> order(summary.feature_names.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return std::abs(mean_weight(a)) > std::abs(mean_weight(b));
        });
        output << "\n  Sample: " << sample.name << '\n'
               << std::left << std::setw(36) << "Feature"
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
        for (const auto index : order) print_weight(summary.feature_names[index], index);
        output << std::noshowpos;
        print_weight("(intercept)", summary.feature_names.size());
    }
    output << "\nOutputs\n";
    for (const auto& prefix : summary.output_prefixes) {
        output << "  " << prefix << "_target_psms.tsv\n"
               << "  " << prefix << "_decoy_psms.tsv\n"
               << "  " << prefix << "_filtered_psms.tsv\n"
               << "  " << prefix << ".pep.xml\n";
    }
    output.flags(old_flags);
    output.precision(old_precision);
    output.fill(old_fill);
}

} // namespace aerith
