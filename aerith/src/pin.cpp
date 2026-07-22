#include "pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

#include <omp.h>

namespace aerith {

std::vector<std::string_view> PinReader::split_tabs(const std::string& line) {
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

double PinReader::parse_number(std::string_view field, const std::string& column,
                               std::size_t line_number) {
    if (field.empty()) {
        return 0.0;
    }
    std::string copy(field);
    std::size_t used = 0;
    double value = 0.0;
    try {
        value = std::stod(copy, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("Non-numeric value in column " + column + " at line " +
                                 std::to_string(line_number));
    }
    if (used != copy.size() || !std::isfinite(value)) {
        throw std::runtime_error("Invalid numeric value in column " + column + " at line " +
                                 std::to_string(line_number));
    }
    return value;
}

std::size_t PinReader::required_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("PIN file is missing required column: " + name);
    }
    return found->second;
}

Dataset PinReader::read_file(const Config& config, const std::string& input_path,
                             std::size_t file_id, int expected_label) {
    std::ifstream input(input_path);
    if (!input) {
        throw std::runtime_error("Cannot open input PIN: " + input_path);
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Input PIN is empty: " + input_path);
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    const auto header_views = split_tabs(line);
    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header_views.size(); ++i) {
        headers.emplace_back(header_views[i]);
        columns.emplace(headers.back(), i);
    }

    const auto id_col = required_column(columns, "SpecId");
    const auto label_col = required_column(columns, "Label");
    const auto scan_col = required_column(columns, "ScanNr");
    const auto peptide_col = required_column(columns, "Peptide");
    const auto proteins_col = required_column(columns, "Proteins");
    const auto rt_col = required_column(columns, "retentiontime");
    const auto mass_col = required_column(columns, "ExpMass");
    const auto initial_col = required_column(columns, config.initial_score);
    const auto merge_col = required_column(columns, "WDPscores");
    const auto rank_col = required_column(columns, "ranks");
    const auto charge_col = required_column(columns, "parentCharges");
    const auto diff_col = columns.find("diffScores");

    const std::unordered_set<std::string> excluded{
        "SpecId", "Label", "ScanNr", "retentiontime", "ExpMass",
        "Peptide", "Proteins", "SampleName"};
    std::vector<std::size_t> numeric_columns;
    Dataset data;
    data.input_paths.push_back(input_path);
    data.headers = headers;
    data.columns = columns;
    for (std::size_t i = 0; i < headers.size(); ++i) {
        const bool pct = headers[i] == "MS1IsotopeFitScore" ||
            headers[i] == "MS1IsotopicAbundances" ||
            headers[i] == "MS2IsotopicAbundances" ||
            headers[i] == "isotopicAbundanceDiffs";
        if (excluded.count(headers[i]) == 0 && !(config.ignore_pct && pct)) {
            numeric_columns.push_back(i);
            data.feature_names.push_back(headers[i]);
        }
    }
    data.numeric_columns = numeric_columns;
    if (numeric_columns.empty()) {
        throw std::runtime_error("PIN file has no numeric model features");
    }

    std::size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        const auto fields = split_tabs(line);
        if (fields.size() != headers.size()) {
            throw std::runtime_error("PIN line " + std::to_string(line_number) + " has " +
                                     std::to_string(fields.size()) + " fields; expected " +
                                     std::to_string(headers.size()));
        }
        Psm row;
        row.raw_line = line;
        row.id = std::string(fields[id_col]);
        row.peptide = std::string(fields[peptide_col]);
        row.proteins = std::string(fields[proteins_col]);
        row.label = static_cast<int>(parse_number(fields[label_col], "Label", line_number));
        if (row.label != 1 && row.label != -1) {
            throw std::runtime_error("Label must be 1 or -1 at line " +
                                     std::to_string(line_number));
        }
        if (expected_label != 0 && row.label != expected_label) {
            throw std::runtime_error(
                "PIN label does not match its target/decoy input at line " +
                std::to_string(line_number) + " in " + input_path);
        }
        row.retention = parse_number(fields[rt_col], "retentiontime", line_number);
        row.exp_mass = parse_number(fields[mass_col], "ExpMass", line_number);
        row.file_id = file_id;
        row.scan = static_cast<std::uint64_t>(
            parse_number(fields[scan_col], "ScanNr", line_number));
        row.initial_score = parse_number(fields[initial_col], config.initial_score, line_number);
        row.merge_score = parse_number(fields[merge_col], "WDPscores", line_number);
        row.rank = static_cast<std::size_t>(
            parse_number(fields[rank_col], "ranks", line_number));
        row.charge = static_cast<int>(
            parse_number(fields[charge_col], "parentCharges", line_number));
        if (row.charge <= 0) {
            throw std::runtime_error("parentCharges must be positive at line " +
                                     std::to_string(line_number));
        }
        if (diff_col != columns.end()) {
            row.diff_score = parse_number(fields[diff_col->second], "diffScores", line_number);
        }
        row.features.reserve(numeric_columns.size());
        for (const auto column : numeric_columns) {
            row.features.push_back(static_cast<float>(
                parse_number(fields[column], headers[column], line_number)));
        }
        data.rows.push_back(std::move(row));
    }
    if (data.rows.empty()) {
        throw std::runtime_error("Input PIN has no PSM rows");
    }
    return data;
}

Dataset PinReader::read(const Config& config) {
    Dataset combined;
    const bool paired = !config.target_pins.empty();
    const std::size_t samples =
        paired ? config.target_pins.size() : config.inputs.size();
    combined.input_paths.resize(samples);
    std::vector<Dataset> parts(paired ? samples * 2 : samples);
    std::exception_ptr failure;
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t task = 0;
         task < static_cast<std::ptrdiff_t>(parts.size()); ++task) {
        try {
            const auto index = static_cast<std::size_t>(task);
            const auto file = paired ? index / 2 : index;
            const int label = paired ? (index % 2 == 0 ? 1 : -1) : 0;
            const auto& path = paired
                ? (label == 1 ? config.target_pins[file] : config.decoy_pins[file])
                : config.inputs[file];
            parts[index] = read_file(config, path, file, label);
        } catch (...) {
            #pragma omp critical(aerith_read_failure)
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
    for (std::size_t file_id = 0; file_id < samples; ++file_id) {
        combined.input_paths[file_id] = paired
            ? config.target_pins[file_id] + " + " + config.decoy_pins[file_id]
            : config.inputs[file_id];
    }
    for (auto& part : parts) {
        if (combined.feature_names.empty()) {
            combined.feature_names = part.feature_names;
            combined.headers = part.headers;
            combined.columns = part.columns;
            combined.numeric_columns = part.numeric_columns;
        } else if (combined.feature_names != part.feature_names) {
            throw std::runtime_error("PIN numeric feature schemas differ between inputs");
        } else if (combined.headers != part.headers) {
            throw std::runtime_error("PIN column schemas differ between inputs");
        }
        for (auto& row : part.rows) {
            combined.rows.push_back(std::move(row));
        }
    }
    if (paired) {
        const auto rank_feature = std::find(
            combined.feature_names.begin(), combined.feature_names.end(), "ranks");
        const auto diff_feature = std::find(
            combined.feature_names.begin(), combined.feature_names.end(), "diffScores");
        #pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t file_index = 0;
             file_index < static_cast<std::ptrdiff_t>(samples); ++file_index) {
            const auto file = static_cast<std::size_t>(file_index);
            std::unordered_map<std::uint64_t, std::vector<std::size_t>> scans;
            for (std::size_t i = 0; i < combined.rows.size(); ++i) {
                if (combined.rows[i].file_id == file) {
                    scans[combined.rows[i].scan].push_back(i);
                }
            }
            for (auto& entry : scans) {
                auto& rows = entry.second;
                std::stable_sort(rows.begin(), rows.end(), [&](std::size_t a, std::size_t b) {
                    return combined.rows[a].merge_score > combined.rows[b].merge_score;
                });
                const double top = combined.rows[rows.front()].merge_score;
                for (std::size_t position = 0; position < rows.size(); ++position) {
                    auto& row = combined.rows[rows[position]];
                    row.rank = position + 1;
                    row.diff_score = top - row.merge_score;
                    if (config.initial_score == "ranks") {
                        row.initial_score = static_cast<double>(row.rank);
                    } else if (config.initial_score == "diffScores") {
                        row.initial_score = row.diff_score;
                    }
                    const auto dot = row.id.rfind('.');
                    if (dot != std::string::npos) {
                        row.id.resize(dot + 1);
                        row.id += std::to_string(row.rank);
                    }
                    if (rank_feature != combined.feature_names.end()) {
                        row.features[static_cast<std::size_t>(
                            rank_feature - combined.feature_names.begin())] =
                            static_cast<float>(row.rank);
                    }
                    if (diff_feature != combined.feature_names.end()) {
                        row.features[static_cast<std::size_t>(
                            diff_feature - combined.feature_names.begin())] =
                            static_cast<float>(row.diff_score);
                    }
                }
            }
        }
        // Match the historical merged-PIN row order without materializing that
        // intermediate file: sample, numeric scan, then merged WDP rank.
        std::vector<std::size_t> sample_offsets(samples + 1, 0);
        for (const auto& row : combined.rows) ++sample_offsets[row.file_id + 1];
        std::partial_sum(
            sample_offsets.begin(), sample_offsets.end(), sample_offsets.begin());
        #pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t file_index = 0;
             file_index < static_cast<std::ptrdiff_t>(samples); ++file_index) {
            const auto file = static_cast<std::size_t>(file_index);
            std::stable_sort(combined.rows.begin() + sample_offsets[file],
                             combined.rows.begin() + sample_offsets[file + 1],
                             [](const Psm& a, const Psm& b) {
                if (a.scan != b.scan) return a.scan < b.scan;
                return a.rank < b.rank;
            });
        }
    }
    return combined;
}

} // namespace aerith
