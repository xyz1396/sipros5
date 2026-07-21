#pragma once

// Internal interfaces shared by Aerith's implementation units.

#include "filter.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace aerith {

constexpr std::size_t kRtFolds = 3;

struct Psm {
    std::string id;
    std::string peptide;
    std::string proteins;
    int label = 0;
    double initial_score = 0.0;
    double retention = 0.0;
    double exp_mass = 0.0;
    std::size_t file_id = 0;
    std::uint64_t scan = 0;
    std::size_t rank = 0;
    double merge_score = 0.0;
    double diff_score = 0.0;
    std::string raw_line;
    std::vector<float> features;
};

struct Dataset {
    std::vector<Psm> rows;
    std::vector<std::string> feature_names;
    std::vector<std::string> input_paths;
    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> columns;
    std::vector<std::size_t> numeric_columns;
};

struct RtResult {
    std::array<std::vector<double>, kRtFolds> residuals;
    double r2 = 0.0;
    std::size_t training_count = 0;
};

struct SvmFit {
    std::vector<double> scores;
    std::vector<std::array<unsigned int, kRtFolds>> iterations;
    std::vector<std::array<std::vector<double>, kRtFolds>> calibrated_weights;
};

class PinReader {
public:
    static Dataset read(const Config& config);

private:
    static std::vector<std::string_view> split_tabs(const std::string& line);
    static double parse_number(std::string_view field, const std::string& column,
                               std::size_t line_number);
    static std::size_t required_column(
        const std::unordered_map<std::string, std::size_t>& columns,
        const std::string& name);
    static Dataset read_file(const Config& config, const std::string& input_path,
                             std::size_t file_id, int expected_label);
};

class RetentionTimeModel {
public:
    static RtResult fit(const Dataset& data,
                        const std::vector<std::size_t>& outer_folds,
                        const std::vector<double>& initial_scores,
                        const std::vector<int>& labels,
                        const std::vector<double>& diagnostic_q,
                        double train_fdr, double ridge);
    static std::string peptide_sequence(const std::string& peptide);
};

class SvmRescorer {
public:
    static std::vector<std::size_t> assign_folds(const Dataset& data);
    static SvmFit fit(
        const Dataset& data,
        const std::array<std::vector<double>, kRtFolds>& fold_rt_features,
        const std::vector<std::size_t>& folds, double train_fdr,
        unsigned int max_iterations, double c_pos, double c_neg);
    static std::vector<double> local_error_probabilities(
        const std::vector<double>& scores, const std::vector<int>& labels);
};

class ResultWriter {
public:
    static void write(
        const Config& config, const Dataset& data,
        const std::vector<double>& scores, const std::vector<double>& q,
        const std::vector<double>& pep, const RtResult& rt,
        const std::vector<std::size_t>& outer_folds);

private:
    static std::vector<std::string_view> split_tabs(const std::string& line);
    static std::size_t required_column(
        const std::unordered_map<std::string, std::size_t>& columns,
        const std::string& name);
    static std::string formatted_number(double value);
    static std::vector<std::string_view> row_fields(
        const Dataset& data, const Psm& psm);
    static void write_original_field(
        std::ostream& stream, const Dataset& data, const Psm& psm,
        const std::vector<std::string_view>& fields, std::size_t column);
    static std::string original_field(
        const Dataset& data, const Psm& psm, const std::string& name);
    static std::vector<std::size_t> selected_rows(
        const Dataset& data, std::size_t file, int label = 0);
    static void write_results(
        const std::string& path, int label, std::size_t file,
        const Dataset& data, const std::vector<double>& scores,
        const std::vector<double>& q, const std::vector<double>& pep);
    static void write_filtered_results(
        const std::string& path, std::size_t file, const Dataset& data,
        const std::vector<double>& scores, const std::vector<double>& q,
        const std::vector<double>& pep,
        const std::vector<double>& rt_residual, double threshold);
    static std::string xml_escape(std::string_view value);
    static std::string trim(std::string value);
    static std::string sanitized_protein(
        std::string protein, const std::string& decoy_prefix);
    static std::vector<std::string> protein_ids(
        const Psm& psm, const std::string& decoy_prefix);
    static void write_pepxml(
        const std::string& path, std::size_t file, const Config& config,
        const Dataset& data, const std::vector<double>& scores,
        const std::vector<double>& pep);
};

} // namespace aerith
