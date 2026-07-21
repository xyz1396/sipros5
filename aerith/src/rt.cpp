#include "pipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <omp.h>

namespace aerith {

class RtImplementation final {
public:
    static RtResult fit(
        const Dataset& data, const std::vector<std::size_t>& outer_folds,
        const std::vector<double>& initial_scores,
        const std::vector<int>& labels,
        const std::vector<double>& diagnostic_q,
        double train_fdr, double ridge);
    static std::string peptide_sequence(const std::string& peptide);

private:
    static constexpr std::size_t kAminoAcids = 20;
    static constexpr std::size_t kRtFeatures = kAminoAcids * 3 + 3;
    static constexpr std::size_t kChemicalFeatures = 29;
    static constexpr std::size_t kEnhancedRtFeatures =
        kRtFeatures + kChemicalFeatures;
    inline static constexpr std::array<char, kAminoAcids> kResidues{
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

    struct ParsedPeptide {
        std::string sequence;
        double oxidation = 0.0;
        double deamidation = 0.0;
        double acetylation = 0.0;
        double phosphorylation = 0.0;
        double methylation = 0.0;
        double other_modification = 0.0;
        std::array<double, 5> atom_delta{};
    };

    struct RtAlignment {
        double slope = 1.0;
        double intercept = 0.0;
    };

    struct RtAlignmentModel {
        std::vector<double> max_rt;
        std::vector<RtAlignment> runs;
    };

    static std::vector<double> solve(
        std::vector<double> matrix, std::vector<double> rhs,
        std::size_t size);
    static std::array<int, 26> residue_map();
    static std::string_view peptide_body(const std::string& peptide);
    static ParsedPeptide parse_peptide(const std::string& peptide);
    static RtAlignmentModel fit_retention_alignment(
        const Dataset& data, const std::vector<double>& training_q,
        const std::vector<std::size_t>& outer_folds,
        std::size_t held_out_fold, double train_fdr);
    static double aligned_retention(
        const Psm& psm, const RtAlignmentModel& model);
    static std::array<double, kRtFeatures> base_rt_embedding(const Psm& psm);
    static double residue_value(
        char aa, const std::array<double, kAminoAcids>& values);
    static double peptide_charge(const std::string& sequence, double ph);
    static double estimated_pi(const std::string& sequence);
    static std::array<double, kEnhancedRtFeatures> enhanced_rt_embedding(
        const Psm& psm);
    static std::size_t peptide_fold(const std::string& sequence);

    template <std::size_t Dimensions, typename Embed>
    static RtResult fit_nested_rt(
        const Dataset& data, const std::vector<std::size_t>& outer_folds,
        const std::vector<double>& initial_scores,
        const std::vector<int>& labels,
        const std::vector<double>& diagnostic_q,
        double train_fdr, double ridge, Embed embed);
};

std::vector<double> RtImplementation::solve(
    std::vector<double> matrix, std::vector<double> rhs, std::size_t size) {
    for (std::size_t col = 0; col < size; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < size; ++row) {
            if (std::abs(matrix[row * size + col]) >
                std::abs(matrix[pivot * size + col])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot * size + col]) < 1e-14) {
            throw std::runtime_error("Regularized model matrix is singular");
        }
        if (pivot != col) {
            for (std::size_t j = col; j < size; ++j) {
                std::swap(matrix[col * size + j], matrix[pivot * size + j]);
            }
            std::swap(rhs[col], rhs[pivot]);
        }
        const double diagonal = matrix[col * size + col];
        for (std::size_t j = col; j < size; ++j) {
            matrix[col * size + j] /= diagonal;
        }
        rhs[col] /= diagonal;
        for (std::size_t row = 0; row < size; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = matrix[row * size + col];
            if (factor == 0.0) {
                continue;
            }
            for (std::size_t j = col; j < size; ++j) {
                matrix[row * size + j] -= factor * matrix[col * size + j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    return rhs;
}

std::array<int, 26> RtImplementation::residue_map() {
    std::array<int, 26> map{};
    map.fill(-1);
    for (std::size_t i = 0; i < kResidues.size(); ++i) {
        map[static_cast<std::size_t>(kResidues[i] - 'A')] = static_cast<int>(i);
    }
    return map;
}

std::string_view RtImplementation::peptide_body(const std::string& peptide) {
    const auto open = peptide.find('[');
    const auto close = peptide.rfind(']');
    if (open != std::string::npos && close > open && open <= 1 &&
        (close + 1 == peptide.size() || close + 2 == peptide.size())) {
        return std::string_view(peptide).substr(open + 1, close - open - 1);
    }
    const auto first_dot = peptide.find('.');
    const auto last_dot = peptide.rfind('.');
    if (first_dot != std::string::npos && last_dot > first_dot) {
        return std::string_view(peptide).substr(first_dot + 1, last_dot - first_dot - 1);
    }
    return peptide;
}

RtImplementation::ParsedPeptide RtImplementation::parse_peptide(
    const std::string& peptide) {
    ParsedPeptide parsed;
    const auto body = peptide_body(peptide);
    bool bracket_annotation = false;
    for (const char ch : body) {
        if (ch == '[') {
            bracket_annotation = true;
            continue;
        }
        if (ch == ']') {
            bracket_annotation = false;
            continue;
        }
        if (bracket_annotation) {
            continue;
        }
        if (ch >= 'A' && ch <= 'Z') {
            parsed.sequence.push_back(ch);
            continue;
        }
        // Atom order is C, H, N, O, S. Values match Sipros' compiled PTM catalog.
        switch (ch) {
        case '~':
            ++parsed.oxidation;
            parsed.atom_delta[3] += 1.0;
            break;
        case '!':
            ++parsed.deamidation;
            parsed.atom_delta[1] -= 1.0;
            parsed.atom_delta[2] -= 1.0;
            parsed.atom_delta[3] += 1.0;
            break;
        case '%':
            ++parsed.acetylation;
            parsed.atom_delta[0] += 2.0;
            parsed.atom_delta[1] += 2.0;
            parsed.atom_delta[3] += 1.0;
            break;
        case '@': case '>': case '<':
            ++parsed.phosphorylation;
            parsed.atom_delta[1] += 1.0;
            parsed.atom_delta[3] += 3.0;
            break;
        case '^':
            ++parsed.methylation;
            parsed.atom_delta[0] += 1.0;
            parsed.atom_delta[1] += 2.0;
            break;
        case '&':
            ++parsed.methylation;
            parsed.atom_delta[0] += 2.0;
            parsed.atom_delta[1] += 4.0;
            break;
        case '*':
            ++parsed.methylation;
            parsed.atom_delta[0] += 3.0;
            parsed.atom_delta[1] += 6.0;
            break;
        case '(':
            ++parsed.other_modification;
            parsed.atom_delta[1] -= 1.0;
            parsed.atom_delta[2] += 1.0;
            parsed.atom_delta[3] += 1.0;
            break;
        case ')':
            ++parsed.other_modification;
            parsed.atom_delta[1] -= 1.0;
            parsed.atom_delta[2] += 1.0;
            parsed.atom_delta[3] += 2.0;
            break;
        case '/':
            ++parsed.other_modification;
            parsed.atom_delta[0] += 2.0;
            parsed.atom_delta[1] += 3.0;
            parsed.atom_delta[2] += 1.0;
            parsed.atom_delta[3] += 1.0;
            break;
        case '$':
            ++parsed.other_modification;
            parsed.atom_delta[0] += 1.0;
            parsed.atom_delta[1] += 2.0;
            parsed.atom_delta[4] += 1.0;
            break;
        default:
            break;
        }
    }
    return parsed;
}

RtImplementation::RtAlignmentModel RtImplementation::fit_retention_alignment(
    const Dataset& data, const std::vector<double>& training_q,
    const std::vector<std::size_t>& outer_folds, std::size_t held_out_fold,
    double train_fdr) {
    const std::size_t files = data.input_paths.size();
    RtAlignmentModel model{std::vector<double>(files, 0.0),
                           std::vector<RtAlignment>(files)};
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        if (outer_folds[i] != held_out_fold) {
            model.max_rt[data.rows[i].file_id] =
                std::max(model.max_rt[data.rows[i].file_id], data.rows[i].retention);
        }
    }
    for (double& maximum : model.max_rt) {
        if (!(maximum > 0.0)) maximum = 1.0;
        maximum = std::ceil(maximum);
    }
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::unordered_map<std::string, std::vector<double>> peptide_rt;
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        const auto& psm = data.rows[i];
        if (outer_folds[i] == held_out_fold || psm.label != 1 ||
            training_q[i] > train_fdr) continue;
        auto entry = peptide_rt.try_emplace(std::string(peptide_body(psm.peptide)),
                                            files, missing).first;
        double& rt = entry->second[psm.file_id];
        if (!std::isfinite(rt) || psm.retention < rt) rt = psm.retention;
    }
    struct PeptideTimes { std::vector<double> normalized; double mean = 0.0; };
    std::vector<PeptideTimes> times;
    times.reserve(peptide_rt.size());
    for (const auto& entry : peptide_rt) {
        PeptideTimes peptide{std::vector<double>(files, missing), 0.0};
        std::size_t observed = 0;
        for (std::size_t file = 0; file < files; ++file) {
            if (std::isfinite(entry.second[file])) {
                peptide.normalized[file] = entry.second[file] / model.max_rt[file];
                peptide.mean += peptide.normalized[file];
                ++observed;
            }
        }
        if (observed) {
            peptide.mean /= observed;
            if (std::isfinite(peptide.mean)) times.push_back(std::move(peptide));
        }
    }
    for (std::size_t file = 0; file < files; ++file) {
        std::size_t count = 0;
        double dot = 0.0, sum_x = 0.0, sum_y = 0.0;
        for (const auto& peptide : times) {
            const double x = peptide.normalized[file];
            if (std::isfinite(x)) {
                ++count; dot += x * peptide.mean; sum_x += x; sum_y += peptide.mean;
            }
        }
        if (!count) continue;
        const double x_mean = sum_x / count, y_mean = sum_y / count;
        double sx2 = 1e-8;
        for (const auto& peptide : times) {
            const double x = peptide.normalized[file];
            if (std::isfinite(x)) sx2 += (x - x_mean) * (x - x_mean);
        }
        const double slope = (dot - count * x_mean * y_mean) / sx2;
        const double intercept = y_mean - slope * x_mean;
        if (std::isfinite(slope)) model.runs[file].slope = slope;
        if (std::isfinite(intercept)) model.runs[file].intercept = intercept;
    }
    return model;
}

double RtImplementation::aligned_retention(
    const Psm& psm, const RtAlignmentModel& model) {
    const auto alignment = model.runs[psm.file_id];
    return (psm.retention / model.max_rt[psm.file_id]) * alignment.slope +
           alignment.intercept;
}
std::array<double, RtImplementation::kRtFeatures>
RtImplementation::base_rt_embedding(const Psm& psm) {
    static const auto map = residue_map();
    std::array<double, kRtFeatures> values{};
    const auto sequence = parse_peptide(psm.peptide).sequence;
    const std::size_t cterm = sequence.size() >= 3 ? sequence.size() - 3 : 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const char aa = sequence[i];
        if (aa < 'A' || aa > 'Z') {
            continue;
        }
        const int index = map[static_cast<std::size_t>(aa - 'A')];
        if (index < 0) {
            continue;
        }
        values[static_cast<std::size_t>(index)] += 1.0;
        if (i == 0 || i == 1) {
            values[kAminoAcids + static_cast<std::size_t>(index)] += 1.0;
        }
        if (i == cterm || i == cterm + 1) {
            values[kAminoAcids * 2 + static_cast<std::size_t>(index)] += 1.0;
        }
    }
    values[kRtFeatures - 3] = static_cast<double>(sequence.size());
    values[kRtFeatures - 2] = std::log1p(std::max(0.0, psm.exp_mass));
    values[kRtFeatures - 1] = 1.0;
    return values;
}

double RtImplementation::residue_value(
    char aa, const std::array<double, kAminoAcids>& values) {
    static const auto map = residue_map();
    if (aa < 'A' || aa > 'Z') {
        return 0.0;
    }
    const int index = map[static_cast<std::size_t>(aa - 'A')];
    return index >= 0 ? values[static_cast<std::size_t>(index)] : 0.0;
}

double RtImplementation::peptide_charge(const std::string& sequence, double ph) {
    const auto count = [&](char aa) {
        return static_cast<double>(std::count(sequence.begin(), sequence.end(), aa));
    };
    const auto positive = [ph](double pka) { return 1.0 / (1.0 + std::pow(10.0, ph - pka)); };
    const auto negative = [ph](double pka) { return 1.0 / (1.0 + std::pow(10.0, pka - ph)); };
    return positive(9.69) + count('H') * positive(6.00) + count('K') * positive(10.50) +
           count('R') * positive(12.40) - negative(2.34) - count('D') * negative(3.86) -
           count('E') * negative(4.25) - count('C') * negative(8.33) -
           count('Y') * negative(10.07);
}

double RtImplementation::estimated_pi(const std::string& sequence) {
    double low = 0.0;
    double high = 14.0;
    for (int iteration = 0; iteration < 40; ++iteration) {
        const double mid = (low + high) / 2.0;
        if (peptide_charge(sequence, mid) > 0.0) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return (low + high) / 2.0;
}

std::array<double, RtImplementation::kEnhancedRtFeatures>
RtImplementation::enhanced_rt_embedding(const Psm& psm) {
    std::array<double, kEnhancedRtFeatures> values{};
    const auto base = base_rt_embedding(psm);
    std::copy(base.begin(), base.end(), values.begin());
    const auto parsed = parse_peptide(psm.peptide);
    const auto& sequence = parsed.sequence;
    const double length = std::max<std::size_t>(1, sequence.size());
    constexpr std::array<double, kAminoAcids> hydrophobicity{
        1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
        1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3};
    // Polymer residue compositions in C, H, N, O, S order.
    constexpr std::array<std::array<double, 5>, kAminoAcids> atoms{{
        {{3,5,1,1,0}}, {{3,5,1,1,1}}, {{4,5,1,3,0}}, {{5,7,1,3,0}},
        {{9,9,1,1,0}}, {{2,3,1,1,0}}, {{6,7,3,1,0}}, {{6,11,1,1,0}},
        {{6,12,2,1,0}}, {{6,11,1,1,0}}, {{5,9,1,1,1}}, {{4,6,2,2,0}},
        {{5,7,1,1,0}}, {{5,8,2,2,0}}, {{6,12,4,1,0}}, {{3,5,1,2,0}},
        {{4,7,1,2,0}}, {{5,9,1,1,0}}, {{11,10,2,1,0}}, {{9,9,1,2,0}}
    }};
    static const auto map = residue_map();
    std::vector<double> hydro;
    hydro.reserve(sequence.size());
    std::array<double, 5> atom_count{0.0, 2.0, 0.0, 1.0, 0.0}; // terminal H2O
    double hydrophobic = 0.0;
    double aromatic = 0.0;
    double polar = 0.0;
    double acidic = 0.0;
    double basic = 0.0;
    for (const char aa : sequence) {
        const double h = residue_value(aa, hydrophobicity);
        hydro.push_back(h);
        hydrophobic += std::string_view("AVILMFWY").find(aa) != std::string_view::npos ? 1.0 : 0.0;
        aromatic += std::string_view("FWY").find(aa) != std::string_view::npos ? 1.0 : 0.0;
        polar += std::string_view("STNQC").find(aa) != std::string_view::npos ? 1.0 : 0.0;
        acidic += aa == 'D' || aa == 'E' ? 1.0 : 0.0;
        basic += aa == 'K' || aa == 'R' || aa == 'H' ? 1.0 : 0.0;
        if (aa >= 'A' && aa <= 'Z') {
            const int index = map[static_cast<std::size_t>(aa - 'A')];
            if (index >= 0) {
                for (std::size_t atom = 0; atom < atom_count.size(); ++atom) {
                    atom_count[atom] += atoms[static_cast<std::size_t>(index)][atom];
                }
            }
        }
    }
    for (std::size_t atom = 0; atom < atom_count.size(); ++atom) {
        atom_count[atom] += parsed.atom_delta[atom];
    }
    const double hydro_sum = std::accumulate(hydro.begin(), hydro.end(), 0.0);
    const double hydro_mean = hydro_sum / length;
    double hydro_variance = 0.0;
    double hydro_abs = 0.0;
    for (const double h : hydro) {
        hydro_variance += (h - hydro_mean) * (h - hydro_mean);
        hydro_abs += std::abs(h);
    }
    const auto autocorrelation = [&](std::size_t lag) {
        if (hydro.size() <= lag) {
            return 0.0;
        }
        double sum = 0.0;
        for (std::size_t i = lag; i < hydro.size(); ++i) {
            sum += hydro[i] * hydro[i - lag];
        }
        return sum / static_cast<double>(hydro.size() - lag);
    };
    const auto terminal_hydro = [&](bool n_terminal) {
        const std::size_t count = std::min<std::size_t>(3, hydro.size());
        if (count == 0) {
            return 0.0;
        }
        const auto begin = n_terminal ? hydro.begin() : hydro.end() - static_cast<std::ptrdiff_t>(count);
        return std::accumulate(begin, begin + static_cast<std::ptrdiff_t>(count), 0.0) / count;
    };
    const double total_modifications = parsed.oxidation + parsed.deamidation +
        parsed.acetylation + parsed.phosphorylation + parsed.methylation +
        parsed.other_modification;
    std::size_t out = kRtFeatures;
    values[out++] = hydro_mean;
    values[out++] = std::sqrt(hydro_variance / length);
    values[out++] = hydro_abs / length;
    values[out++] = autocorrelation(1);
    values[out++] = autocorrelation(2);
    values[out++] = hydrophobic / length;
    values[out++] = aromatic / length;
    values[out++] = polar / length;
    values[out++] = acidic / length;
    values[out++] = basic / length;
    values[out++] = static_cast<double>(std::count(sequence.begin(), sequence.end(), 'P')) / length;
    values[out++] = static_cast<double>(std::count(sequence.begin(), sequence.end(), 'G')) / length;
    values[out++] = estimated_pi(sequence);
    values[out++] = peptide_charge(sequence, 2.5);
    values[out++] = psm.exp_mass / length;
    values[out++] = parsed.oxidation;
    values[out++] = parsed.deamidation;
    values[out++] = parsed.acetylation;
    values[out++] = parsed.phosphorylation;
    values[out++] = parsed.methylation;
    values[out++] = parsed.other_modification;
    values[out++] = total_modifications / length;
    values[out++] = terminal_hydro(true);
    values[out++] = terminal_hydro(false);
    for (const double count : atom_count) {
        values[out++] = count / length;
    }
    return values;
}


std::size_t RtImplementation::peptide_fold(const std::string& sequence) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (const unsigned char ch : sequence) {
        hash ^= ch;
        hash *= 1099511628211ULL;
    }
    return static_cast<std::size_t>(hash % kRtFolds);
}

template <std::size_t Dimensions, typename Embed>
RtResult RtImplementation::fit_nested_rt(
    const Dataset& data, const std::vector<std::size_t>& outer_folds,
    const std::vector<double>& initial_scores, const std::vector<int>& labels,
    const std::vector<double>& diagnostic_q, double train_fdr, double ridge,
    Embed embed) {
    struct Accumulator {
        std::vector<double> xtx = std::vector<double>(Dimensions * Dimensions);
        std::vector<double> xty = std::vector<double>(Dimensions);
        std::size_t count = 0;
    };
    struct FoldModel {
        std::vector<double> training_q;
        RtAlignmentModel alignment;
        std::array<Accumulator, kRtFolds + 1> statistics;
        std::vector<double> full_beta;
        std::array<std::vector<double>, kRtFolds> inner_beta;
        bool valid = false;
    };
    struct Metric {
        double sum = 0.0;
        double sum2 = 0.0;
        double sse = 0.0;
        std::size_t count = 0;
    };
    auto coefficients = [&](const Accumulator& all, const Accumulator* excluded) {
        std::vector<double> matrix(Dimensions * Dimensions), rhs(Dimensions);
        for (std::size_t j = 0; j < Dimensions; ++j) {
            rhs[j] = all.xty[j] - (excluded ? excluded->xty[j] : 0.0);
            for (std::size_t k = 0; k <= j; ++k) {
                const double value = all.xtx[j * Dimensions + k] -
                    (excluded ? excluded->xtx[j * Dimensions + k] : 0.0);
                matrix[j * Dimensions + k] = value;
                matrix[k * Dimensions + j] = value;
            }
        }
        double trace = 0.0;
        for (std::size_t j = 0; j < Dimensions; ++j)
            trace += matrix[j * Dimensions + j];
        const double mean_diagonal = trace / Dimensions;
        for (std::size_t j = 0; j < Dimensions; ++j) {
            if (j != kRtFeatures - 1) {
                const double feature_scale = std::max(
                    std::abs(matrix[j * Dimensions + j]), mean_diagonal * 1e-12);
                matrix[j * Dimensions + j] +=
                    std::max(1e-12, ridge * feature_scale);
            }
        }
        return solve(std::move(matrix), std::move(rhs), Dimensions);
    };

    using CachedEmbedding = std::array<float, Dimensions>;
    std::vector<CachedEmbedding> embeddings(data.rows.size());
    std::vector<std::uint8_t> peptide_folds(data.rows.size());
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t row = 0;
         row < static_cast<std::ptrdiff_t>(data.rows.size()); ++row) {
        const auto i = static_cast<std::size_t>(row);
        const auto values = embed(data.rows[i]);
        for (std::size_t feature = 0; feature < Dimensions; ++feature)
            embeddings[i][feature] = static_cast<float>(values[feature]);
        peptide_folds[i] = static_cast<std::uint8_t>(
            peptide_fold(parse_peptide(data.rows[i].peptide).sequence));
    }

    std::array<FoldModel, kRtFolds> models;
    #pragma omp parallel for schedule(static) num_threads(kRtFolds)
    for (int outer_index = 0; outer_index < static_cast<int>(kRtFolds);
         ++outer_index) {
        const auto outer = static_cast<std::size_t>(outer_index);
        std::vector<std::size_t> training_rows;
        std::vector<double> subset_scores;
        std::vector<int> subset_labels;
        training_rows.reserve(data.rows.size() * (kRtFolds - 1) / kRtFolds);
        subset_scores.reserve(training_rows.capacity());
        subset_labels.reserve(training_rows.capacity());
        for (std::size_t i = 0; i < data.rows.size(); ++i) {
            if (outer_folds[i] == outer) continue;
            training_rows.push_back(i);
            subset_scores.push_back(initial_scores[i]);
            subset_labels.push_back(labels[i]);
        }
        const auto subset_q = target_decoy_qvalues(subset_scores, subset_labels);
        models[outer].training_q.assign(data.rows.size(), 1.0);
        for (std::size_t i = 0; i < training_rows.size(); ++i)
            models[outer].training_q[training_rows[i]] = subset_q[i];
        models[outer].alignment = fit_retention_alignment(
            data, models[outer].training_q, outer_folds, outer, train_fdr);
    }

    RtResult result;
    for (auto& residuals : result.residuals) residuals.resize(data.rows.size());
    using FoldStatistics = std::array<Accumulator, kRtFolds + 1>;
    using ThreadStatistics = std::array<FoldStatistics, kRtFolds>;
    // Fixed shards make floating-point reductions reproducible at different
    // OpenMP thread counts while still exposing enough work for large teams.
    constexpr std::size_t reduction_shards = 48;
    constexpr std::ptrdiff_t reduction_shard_count =
        static_cast<std::ptrdiff_t>(reduction_shards);
    std::array<ThreadStatistics, reduction_shards> private_statistics;
    std::array<std::array<Metric, kRtFolds>, reduction_shards> private_metrics{};
    const auto work_items = static_cast<std::ptrdiff_t>(
        data.rows.size() * kRtFolds);

    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for (std::ptrdiff_t shard_index = 0;
             shard_index < static_cast<std::ptrdiff_t>(reduction_shards);
             ++shard_index) {
            const auto shard = static_cast<std::size_t>(shard_index);
            const auto begin = work_items * shard_index / reduction_shard_count;
            const auto end = work_items * (shard_index + 1) / reduction_shard_count;
            for (std::ptrdiff_t item = begin; item < end; ++item) {
                const auto outer = static_cast<std::size_t>(item) / data.rows.size();
                const auto i = static_cast<std::size_t>(item) % data.rows.size();
                const auto& psm = data.rows[i];
                if (outer_folds[i] == outer || psm.label != 1 ||
                    models[outer].training_q[i] > train_fdr) continue;
                const auto& x = embeddings[i];
                const double y = aligned_retention(psm, models[outer].alignment);
                const auto peptide_bucket =
                    static_cast<std::size_t>(peptide_folds[i]) + 1;
                for (const auto bucket : {std::size_t{0}, peptide_bucket}) {
                    auto& accumulator = private_statistics[shard][outer][bucket];
                    for (std::size_t j = 0; j < Dimensions; ++j) {
                        accumulator.xty[j] += static_cast<double>(x[j]) * y;
                        for (std::size_t k = 0; k <= j; ++k) {
                            accumulator.xtx[j * Dimensions + k] +=
                                static_cast<double>(x[j]) * x[k];
                        }
                    }
                    ++accumulator.count;
                }
            }
        }

        #pragma omp single
        {
            for (std::size_t outer = 0; outer < kRtFolds; ++outer) {
                for (std::size_t bucket = 0; bucket <= kRtFolds; ++bucket) {
                    auto& destination = models[outer].statistics[bucket];
                    for (const auto& thread_statistics : private_statistics) {
                        const auto& source = thread_statistics[outer][bucket];
                        destination.count += source.count;
                        for (std::size_t j = 0; j < Dimensions; ++j)
                            destination.xty[j] += source.xty[j];
                        for (std::size_t index = 0;
                             index < Dimensions * Dimensions; ++index)
                            destination.xtx[index] += source.xtx[index];
                    }
                }
            }
        }

        #pragma omp for schedule(static)
        for (int outer_index = 0; outer_index < static_cast<int>(kRtFolds);
             ++outer_index) {
            auto& model = models[static_cast<std::size_t>(outer_index)];
            const auto& all = model.statistics[0];
            bool enough = all.count >= Dimensions * 2;
            for (std::size_t fold = 0; fold < kRtFolds; ++fold) {
                enough = enough &&
                    all.count - model.statistics[fold + 1].count >= Dimensions * 2;
            }
            if (!enough) continue;
            model.full_beta = coefficients(all, nullptr);
            for (std::size_t fold = 0; fold < kRtFolds; ++fold) {
                model.inner_beta[fold] = coefficients(
                    all, &model.statistics[fold + 1]);
            }
            model.valid = true;
        }

        #pragma omp for schedule(static)
        for (std::ptrdiff_t shard_index = 0;
             shard_index < static_cast<std::ptrdiff_t>(reduction_shards);
             ++shard_index) {
            const auto shard = static_cast<std::size_t>(shard_index);
            const auto begin = work_items * shard_index / reduction_shard_count;
            const auto end = work_items * (shard_index + 1) / reduction_shard_count;
            for (std::ptrdiff_t item = begin; item < end; ++item) {
                const auto outer = static_cast<std::size_t>(item) / data.rows.size();
                const auto i = static_cast<std::size_t>(item) % data.rows.size();
                const auto& model = models[outer];
                if (!model.valid) continue;
                const auto& psm = data.rows[i];
                const auto& beta = outer_folds[i] == outer
                    ? model.full_beta
                    : model.inner_beta[static_cast<std::size_t>(peptide_folds[i])];
                const auto& x = embeddings[i];
                double raw_prediction = 0.0;
                #pragma omp simd reduction(+:raw_prediction)
                for (std::size_t feature = 0; feature < Dimensions; ++feature)
                    raw_prediction += static_cast<double>(x[feature]) * beta[feature];
                const double observed = aligned_retention(psm, model.alignment);
                const double prediction = std::clamp(raw_prediction, 0.0, 1.0);
                const double delta =
                    std::clamp(std::abs(observed - prediction), 0.001, 0.999);
                result.residuals[outer][i] = std::sqrt(delta);
                if (outer_folds[i] == outer && psm.label == 1 &&
                    diagnostic_q[i] <= train_fdr) {
                    auto& metric = private_metrics[shard][outer];
                    metric.sum += observed;
                    metric.sum2 += observed * observed;
                    metric.sse +=
                        (raw_prediction - observed) * (raw_prediction - observed);
                    ++metric.count;
                }
            }
        }
    }

    std::array<double, kRtFolds> metric_sum{}, metric_sum2{}, metric_sse{};
    std::array<std::size_t, kRtFolds> metric_count{};
    for (std::size_t outer = 0; outer < kRtFolds; ++outer) {
        if (!models[outer].valid) return RtResult{};
        result.training_count += models[outer].statistics[0].count;
        for (const auto& metrics : private_metrics) {
            metric_sum[outer] += metrics[outer].sum;
            metric_sum2[outer] += metrics[outer].sum2;
            metric_sse[outer] += metrics[outer].sse;
            metric_count[outer] += metrics[outer].count;
        }
    }
    for (std::size_t outer = 0; outer < kRtFolds; ++outer) {
        models[outer].training_q.clear();
        models[outer].training_q.shrink_to_fit();
    }
    result.training_count /= kRtFolds;
    const double sum = std::accumulate(metric_sum.begin(), metric_sum.end(), 0.0);
    const double sum2 = std::accumulate(metric_sum2.begin(), metric_sum2.end(), 0.0);
    const double sse = std::accumulate(metric_sse.begin(), metric_sse.end(), 0.0);
    const std::size_t count =
        std::accumulate(metric_count.begin(), metric_count.end(), std::size_t{0});
    const double centered = count ? sum2 - sum * sum / count : 0.0;
    result.r2 = centered > 0.0 ? 1.0 - sse / centered : 0.0;
    return result;
}

RtResult RtImplementation::fit(
    const Dataset& data, const std::vector<std::size_t>& outer_folds,
    const std::vector<double>& initial_scores, const std::vector<int>& labels,
    const std::vector<double>& diagnostic_q, double train_fdr, double ridge) {
    return fit_nested_rt<kEnhancedRtFeatures>(
        data, outer_folds, initial_scores, labels, diagnostic_q,
        train_fdr, ridge, enhanced_rt_embedding);
}

std::string RtImplementation::peptide_sequence(const std::string& peptide) {
    return parse_peptide(peptide).sequence;
}

RtResult RetentionTimeModel::fit(
    const Dataset& data, const std::vector<std::size_t>& outer_folds,
    const std::vector<double>& initial_scores, const std::vector<int>& labels,
    const std::vector<double>& diagnostic_q, double train_fdr, double ridge) {
    return RtImplementation::fit(
        data, outer_folds, initial_scores, labels, diagnostic_q,
        train_fdr, ridge);
}

std::string RetentionTimeModel::peptide_sequence(const std::string& peptide) {
    return RtImplementation::peptide_sequence(peptide);
}

std::string stripped_peptide(const std::string& peptide) {
    return RetentionTimeModel::peptide_sequence(peptide);
}

} // namespace aerith
