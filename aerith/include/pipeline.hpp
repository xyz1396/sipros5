#pragma once

// Internal interfaces shared by Aerith's implementation units.

#include "filter.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
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
    double observed_mass = 0.0;
    double calculated_mass = 0.0;
    double calculated_mz = 0.0;
    std::size_t file_id = 0;
    std::uint64_t scan = 0;
    std::size_t rank = 0;
    int charge = 0;
    int missed_cleavages = 0;
    int ptm_count = 0;
    double merge_score = 0.0;
    double diff_score = 0.0;
    double mass_error = 0.0;
    double log10_precursor_intensity = 0.0;
    double ms1_isotopic_abundance = 0.0;
    double ms2_isotopic_abundance = 0.0;
    int sip_abundance_bin = -1;
    std::string sample_name;
    double quantified_intensity = 0.0;
    double apex_retention = 0.0;
    double retention_start = 0.0;
    double retention_end = 0.0;
    double retention_fwhm = 0.0;
    double quant_mass_error_ppm = 0.0;
    double quant_isotope_kl = 0.0;
    double quant_isotope_correlation = 0.0;
    double quant_isotope_fraction = 0.0;
    double quant_isotope_apex_spread = 0.0;
    std::uint64_t parent_scan = 0;
    std::uint64_t apex_scan = 0;
    std::size_t traced_scans = 0;
    bool quantification_attempted = false;
    bool has_chromatographic_feature = false;
    float delta_rt_loess_real = 0.0f;
    float predicted_rt_real_units = 0.0f;
    std::vector<float> features;
};

// Canonical, human-readable rendering of Sipros' compact peptide PTM tokens.
// Protein assembly and the flat Aerith score tables share this conversion so
// a peptide is annotated identically in every output.
struct ModificationInfo {
    std::string sequence;
    std::string modified_peptide;
    std::vector<std::string> assigned;
    std::string localization;
};

ModificationInfo modification_info(
    const std::string& decorated, bool fixed_cam);

struct TransferredIon {
    Psm psm;
    double score = 0.0;
    double qvalue = 1.0;
    std::size_t donor_row = std::numeric_limits<std::size_t>::max();
    std::string donor_psm_id;
};

struct Dataset {
    std::vector<Psm> rows;
    std::vector<std::string> feature_names;
    std::vector<std::string> input_paths;
    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> columns;
    std::vector<std::size_t> numeric_columns;
    std::vector<std::string> generated_feature_names;
    bool has_predicted_rt_diagnostics = false;
    std::string spectrum_prediction_device;
    std::string rt_prediction_device;
    StageTiming spectrum_prediction_timing;
    StageTiming rt_prediction_timing;
    std::vector<TransferredIon> transferred_ions;
    std::size_t removed_decoy_peptide_collisions = 0;
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
    std::vector<char> sample_valid;
    std::vector<std::string> sample_warnings;
};

class PinReader {
public:
    static Dataset read(
        const Config& config,
        const std::unordered_set<std::string>* global_target_peptides =
            nullptr);
    static Dataset discover_predictions(
        const Config& config,
        std::unordered_set<std::string>& global_target_peptides);

private:
    static std::vector<std::string_view> split_tabs(const std::string& line);
    static double parse_number(std::string_view field, const std::string& column,
                               std::size_t line_number);
    static std::size_t required_column(
        const std::unordered_map<std::string, std::size_t>& columns,
        const std::string& name);
    static Dataset read_file(
        const Config& config, const std::string& input_path,
        std::size_t file_id, int expected_label,
        Psm* destination, std::size_t expected_rows);
};

// Prediction discovery may encounter a decoy form before the same normalized
// peptide/charge key appears as a target in a later sample.  Keep one
// exemplar per model key, but always promote the target exemplar so a
// targets-only cache cannot accidentally omit a real target prediction.
void upsert_prediction_exemplar(
    Dataset& catalog,
    std::unordered_map<std::string, std::size_t>& key_to_row,
    const std::string& key,
    const Psm& exemplar);

class SpectrumPredictionLibrary {
public:
    struct Impl;
    SpectrumPredictionLibrary();
    ~SpectrumPredictionLibrary();
    SpectrumPredictionLibrary(SpectrumPredictionLibrary&&) noexcept;
    SpectrumPredictionLibrary& operator=(
        SpectrumPredictionLibrary&&) noexcept;
    SpectrumPredictionLibrary(const SpectrumPredictionLibrary&) = delete;
    SpectrumPredictionLibrary& operator=(
        const SpectrumPredictionLibrary&) = delete;

    std::string device() const;
    StageTiming timing() const;

private:
    std::unique_ptr<Impl> impl_;
    friend class SpectralEntropyFeature;
};

class SpectralEntropyFeature {
public:
    static void add(const Config& config, Dataset& data);
    static SpectrumPredictionLibrary predict(
        const Config& config, const Dataset& unique_peptides);
    static void add(
        const Config& config, Dataset& data,
        const SpectrumPredictionLibrary& predictions);
};

class RtPredictionLibrary {
public:
    struct Impl;
    RtPredictionLibrary();
    ~RtPredictionLibrary();
    RtPredictionLibrary(RtPredictionLibrary&&) noexcept;
    RtPredictionLibrary& operator=(RtPredictionLibrary&&) noexcept;
    RtPredictionLibrary(const RtPredictionLibrary&) = delete;
    RtPredictionLibrary& operator=(const RtPredictionLibrary&) = delete;

    std::string device() const;
    StageTiming timing() const;

private:
    std::unique_ptr<Impl> impl_;
    friend class PredictedRetentionTimeFeature;
};

class PredictedRetentionTimeFeature {
public:
    static void add(const Config& config, Dataset& data);
    static RtPredictionLibrary predict(
        const Config& config, const Dataset& unique_peptides);
    static void add(
        const Config& config, Dataset& data,
        const RtPredictionLibrary& predictions);
};

#if defined(AERITH_WITH_TORCH) || defined(AERITH_TEST_WITH_TORCH)
// Exposes the shared DIA-NN tokenizer to the regression test. Regular and SIP
// searches both call the same encoder through PredictedRetentionTimeFeature.
std::vector<std::int64_t> diann_rt_tokens_for_testing(const Psm& psm);
// Scores one synthetic experimental SIP spectrum with the matching and an
// intentionally wrong abundance envelope through the production entropy path.
std::array<float, 4> sip_entropy_abundance_scores_for_testing();
#endif

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
    static std::string formatted_number(double value);
    static void write_original_field(
        std::ostream& stream, const Dataset& data, const Psm& psm,
        std::size_t column);
    static std::vector<std::size_t> selected_rows(
        const Dataset& data, std::size_t file, int label = 0);
    static void write_results(
        const std::string& path, int label, std::size_t file,
        const Dataset& data, const std::vector<double>& scores,
        const std::vector<double>& q, const std::vector<double>& pep,
        bool fixed_cam);
    static void write_filtered_results(
        const std::string& path, std::size_t file, const Dataset& data,
        const std::vector<double>& scores, const std::vector<double>& q,
        const std::vector<double>& pep,
        const std::vector<double>& rt_residual, double threshold,
        bool sip_output, bool fixed_cam);
};

struct QuantificationResult {
    std::vector<AccelerationTiming> stages;
};

class ChromatographicQuantifier {
public:
    static QuantificationResult add(
        const Config& config, Dataset& data, const std::vector<double>& q);
};

struct ProteinAssemblyResult {
    std::size_t proteins = 0;
    std::string output_dir;
    std::vector<AccelerationTiming> stages;
};

struct NegativeControlResult {
    std::size_t candidates = 0;
    std::size_t input_psms = 0;
    std::size_t threshold_filtered_psms = 0;
    std::size_t targets = 0;
    std::size_t decoys = 0;
    std::size_t target_ids = 0;
    std::string output_path;
    std::string target_output_path;
    std::string decoy_output_path;
    std::string protein_output_path;
    StageTiming timing;
    std::vector<AccelerationTiming> stages;
    std::vector<std::string> feature_names;
    SampleModelSummary model;
};

class NegativeControlFilter {
public:
    static NegativeControlResult run(
        const Config& config, const Dataset& source,
        const std::vector<double>& source_q);
};

class ProteinAssembler {
public:
    static void sequential_filter(
        const Config& config, const Dataset& data,
        const std::vector<double>& scores, const std::vector<double>& pep,
        std::vector<double>& q);
    static ProteinAssemblyResult write(
        const Config& config, const Dataset& data,
        const std::vector<double>& scores, const std::vector<double>& q,
        const std::vector<double>& pep);
    static void write_sip_psm_mapping(
        const Config& config, const Dataset& data,
        const std::vector<double>& q);
};

} // namespace aerith
