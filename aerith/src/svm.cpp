#include "pipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace aerith {

class SvmImplementation final {
public:
    static std::vector<std::size_t> assign_folds(const Dataset& data);
    static SvmFit fit(
        const Dataset& data,
        const std::array<std::vector<double>, kRtFolds>& fold_extra,
        const std::vector<std::size_t>& folds, double train_fdr,
        unsigned int max_iterations, double c_pos, double c_neg);
    static std::vector<double> local_error_probabilities(
        const std::vector<double>& scores, const std::vector<int>& labels);

private:
    struct SvmMatrix {
        std::size_t row_begin = 0;
        std::size_t rows = 0;
        std::size_t columns = 0;
        std::vector<float> values;

        const float* row(std::size_t i) const {
            return values.data() + (i - row_begin) * columns;
        }
    };

    static SvmMatrix make_matrix(
        const Dataset& data, const std::vector<double>& extra,
        const std::vector<std::size_t>& training_rows,
        std::size_t row_begin, std::size_t row_end);
    static double score(
        const SvmMatrix& matrix, std::size_t row,
        const std::vector<double>& weights);
    static std::vector<double> training_qvalues(
        const std::vector<double>& scores, const std::vector<int>& labels);
    static std::size_t count_confident(
        const std::vector<double>& scores, const std::vector<int>& labels,
        double fdr);
    static std::vector<double> initial_direction(
        const SvmMatrix& matrix, const Dataset& data,
        const std::vector<std::size_t>& training_rows, double train_fdr);
    static std::vector<double> fit_svm(
        const SvmMatrix& matrix, const Dataset& data,
        const std::vector<std::size_t>& training_rows,
        const std::vector<double>& training_q, double train_fdr,
        double c_pos, double c_neg);
    static std::pair<double, double> training_score_calibration(
        const std::vector<double>& scores, const std::vector<int>& labels,
        double fdr);
};

SvmImplementation::SvmMatrix SvmImplementation::make_matrix(
    const Dataset& data, const std::vector<double>& extra,
    const std::vector<std::size_t>& training_rows,
    std::size_t row_begin, std::size_t row_end) {
    SvmMatrix matrix;
    matrix.row_begin = row_begin;
    matrix.rows = row_end - row_begin;
    matrix.columns = data.feature_names.size() + 1;
    auto raw = [&](std::size_t i, std::size_t j) {
        return j < data.feature_names.size()
                   ? static_cast<double>(data.rows[i].features[j]) : extra[i];
    };
    std::vector<double> mean(matrix.columns), sum2(matrix.columns), scale(matrix.columns, 1.0);
    for (const auto i : training_rows) {
        for (std::size_t j = 0; j < matrix.columns; ++j) {
            const double value = raw(i, j);
            mean[j] += value;
            sum2[j] += value * value;
        }
    }
    for (std::size_t j = 0; j < matrix.columns; ++j) {
        mean[j] /= training_rows.size();
        const double variance = sum2[j] / training_rows.size() - mean[j] * mean[j];
        if (variance > 1e-20) scale[j] = std::sqrt(variance);
    }
    matrix.values.resize(matrix.rows * matrix.columns);
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t row_index = 0;
         row_index < static_cast<std::ptrdiff_t>(matrix.rows); ++row_index) {
        const auto local_row = static_cast<std::size_t>(row_index);
        const auto i = row_begin + local_row;
        #pragma omp simd
        for (std::size_t j = 0; j < matrix.columns; ++j) {
            matrix.values[local_row * matrix.columns + j] =
                static_cast<float>((raw(i, j) - mean[j]) / scale[j]);
        }
    }
    return matrix;
}

double SvmImplementation::score(const SvmMatrix& matrix, std::size_t row,
                                const std::vector<double>& weights) {
    const float* values = matrix.row(row);
    double score = weights.back();
    #pragma omp simd reduction(+:score)
    for (std::size_t j = 0; j < matrix.columns; ++j) {
        score += static_cast<double>(values[j]) * weights[j];
    }
    return score;
}

std::vector<double> SvmImplementation::training_qvalues(
    const std::vector<double>& scores, const std::vector<int>& labels) {
    return mixmax_qvalues(scores, labels);
}

std::size_t SvmImplementation::count_confident(
    const std::vector<double>& scores, const std::vector<int>& labels,
    double fdr) {
    const auto q = training_qvalues(scores, labels);
    std::size_t count = 0;
    for (std::size_t i = 0; i < q.size(); ++i)
        count += labels[i] == 1 && q[i] <= fdr ? 1u : 0u;
    return count;
}

std::vector<std::size_t> SvmImplementation::assign_folds(const Dataset& data) {
    std::vector<std::size_t> folds(data.rows.size());
    for (std::size_t file = 0; file < data.input_paths.size(); ++file) {
        std::vector<std::size_t> order;
        for (std::size_t i = 0; i < data.rows.size(); ++i)
            if (data.rows[i].file_id == file) order.push_back(i);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return data.rows[a].scan < data.rows[b].scan;
        });
        std::array<long long, kRtFolds> remaining{};
        std::size_t left = order.size();
        for (std::size_t n = kRtFolds; n > 0; --n) {
            remaining[n - 1] = static_cast<long long>(left / n);
            left -= static_cast<std::size_t>(remaining[n - 1]);
        }
        std::uint64_t seed = 1;
        auto draw = [&]() {
            seed = (seed * 279470273u) % 4294967291u;
            return static_cast<std::size_t>(seed % kRtFolds);
        };
        std::size_t groups = 0;
        for (std::size_t begin = 0; begin < order.size();) {
            std::size_t end = begin + 1;
            while (end < order.size() &&
                   data.rows[order[end]].scan == data.rows[order[begin]].scan) ++end;
            std::size_t fold = draw();
            while (remaining[fold] <= 0) fold = draw();
            for (std::size_t i = begin; i < end; ++i) {
                folds[order[i]] = fold;
                --remaining[fold];
            }
            ++groups;
            begin = end;
        }
        if (groups < kRtFolds) {
            throw std::runtime_error(
                "Each sample requires at least three spectra for SVM filtering");
        }
    }
    return folds;
}

std::vector<double> SvmImplementation::initial_direction(
    const SvmMatrix& matrix, const Dataset& data,
    const std::vector<std::size_t>& rows, double fdr) {
    std::vector<int> labels(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) labels[i] = data.rows[rows[i]].label;
    std::size_t best_count = 0, best_feature = 0;
    double best_sign = 1.0;
    for (std::size_t feature = 0; feature < matrix.columns; ++feature) {
        for (double sign : {1.0, -1.0}) {
            std::vector<double> scores(rows.size());
            for (std::size_t i = 0; i < rows.size(); ++i)
                scores[i] = sign * matrix.row(rows[i])[feature];
            const auto found = count_confident(scores, labels, fdr);
            if (found > best_count) {
                best_count = found;
                best_feature = feature;
                best_sign = sign;
            }
        }
    }
    if (best_count == 0)
        throw std::runtime_error("Cannot initialize SVM: no confident target PSMs");
    std::vector<double> weights(matrix.columns + 1, 0.0);
    weights[best_feature] = best_sign;
    return weights;
}

std::vector<double> SvmImplementation::fit_svm(
    const SvmMatrix& matrix, const Dataset& data,
    const std::vector<std::size_t>& rows, const std::vector<double>& q,
    double fdr, double c_pos, double c_neg) {
    std::vector<std::size_t> selected;
    std::size_t positives = 0, negatives = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (data.rows[rows[i]].label == -1) {
            selected.push_back(rows[i]); ++negatives;
        } else if (q[i] <= fdr) {
            selected.push_back(rows[i]); ++positives;
        }
    }
    if (positives == 0 || negatives == 0)
        throw std::runtime_error("SVM training requires confident targets and decoys");
    const std::size_t dimensions = matrix.columns + 1;
    std::vector<double> alpha(selected.size(), 0.0);
    std::vector<double> diagonal(selected.size());
    std::vector<double> quadratic(selected.size());
    std::vector<std::size_t> order(selected.size());
    std::iota(order.begin(), order.end(), 0);
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const double cost = data.rows[selected[i]].label == 1 ? c_pos : c_neg;
        diagonal[i] = 1.0 / cost;
        const float* values = matrix.row(selected[i]);
        double norm = 1.0; // regularized intercept feature
        #pragma omp simd reduction(+:norm)
        for (std::size_t j = 0; j < matrix.columns; ++j) {
            norm += static_cast<double>(values[j]) * values[j];
        }
        quadratic[i] = norm + diagonal[i];
    }

    std::vector<double> weights(dimensions, 0.0);
    std::size_t active = selected.size();
    double pg_max_old = std::numeric_limits<double>::infinity();
    std::uint64_t random_state = 1;
    constexpr double tolerance = 0.1;
    constexpr unsigned int max_solver_iterations = 300;
    for (unsigned int iteration = 0; iteration < max_solver_iterations; ++iteration) {
        double pg_max = -std::numeric_limits<double>::infinity();
        double pg_min = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < active; ++i) {
            random_state = (random_state * 279470273u) % 4294967291u;
            const std::size_t swap_index = i + random_state % (active - i);
            std::swap(order[i], order[swap_index]);
        }
        std::size_t position = 0;
        while (position < active) {
            const std::size_t i = order[position];
            const int label = data.rows[selected[i]].label;
            const double gradient = label * score(matrix, selected[i], weights) - 1.0 +
                                    alpha[i] * diagonal[i];
            double projected = 0.0;
            if (alpha[i] == 0.0) {
                if (gradient > pg_max_old) {
                    --active;
                    std::swap(order[position], order[active]);
                    continue;
                }
                if (gradient < 0.0) projected = gradient;
            } else {
                projected = gradient;
            }
            pg_max = std::max(pg_max, projected);
            pg_min = std::min(pg_min, projected);
            if (std::abs(projected) > 1e-12) {
                const double old = alpha[i];
                alpha[i] = std::max(0.0, old - gradient / quadratic[i]);
                const double update = (alpha[i] - old) * label;
                const float* values = matrix.row(selected[i]);
                #pragma omp simd
                for (std::size_t j = 0; j < matrix.columns; ++j) {
                    weights[j] += update * values[j];
                }
                weights.back() += update;
            }
            ++position;
        }
        if (pg_max - pg_min <= tolerance && std::abs(pg_max) <= tolerance &&
            std::abs(pg_min) <= tolerance) {
            if (active == selected.size()) break;
            active = selected.size();
            pg_max_old = std::numeric_limits<double>::infinity();
            continue;
        }
        pg_max_old = pg_max <= 0.0 ? std::numeric_limits<double>::infinity() : pg_max;
    }
    return weights;
}

std::pair<double, double> SvmImplementation::training_score_calibration(
    const std::vector<double>& scores, const std::vector<int>& labels, double fdr) {
    const auto q = training_qvalues(scores, labels);
    std::vector<std::size_t> order(scores.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b];
    });
    const std::size_t median = std::max<std::size_t>(1,
        static_cast<std::size_t>(std::count(labels.begin(), labels.end(), -1)) / 2);
    std::size_t decoys = 0;
    double fdr_score = scores[order.front()], median_decoy = fdr_score + 1.0;
    for (const auto i : order) {
        if (q[i] < fdr) fdr_score = scores[i];
        if (labels[i] == -1 && ++decoys == median) {
            median_decoy = scores[i];
            break;
        }
    }
    double scale = fdr_score - median_decoy;
    if (!(scale > 0.0)) {
        scale = 1.0;
        fdr_score += 1.0;
    }
    return {fdr_score, scale};
}
SvmFit SvmImplementation::fit(
    const Dataset& data,
    const std::array<std::vector<double>, kRtFolds>& fold_extra,
    const std::vector<std::size_t>& folds, double train_fdr,
    unsigned int max_iterations, double c_pos, double c_neg) {
    std::vector<std::pair<std::size_t, std::size_t>> ranges(
        data.input_paths.size());
    std::size_t cursor = 0;
    for (std::size_t file = 0; file < ranges.size(); ++file) {
        const std::size_t begin = cursor;
        while (cursor < data.rows.size() && data.rows[cursor].file_id == file) ++cursor;
        ranges[file] = {begin, cursor};
    }

    SvmFit result;
    result.scores.resize(data.rows.size());
    result.iterations.resize(data.input_paths.size());
    result.calibrated_weights.resize(data.input_paths.size());
    const int jobs = static_cast<int>(data.input_paths.size() * kRtFolds);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int job = 0; job < jobs; ++job) {
        const auto file = static_cast<std::size_t>(job) / kRtFolds;
        const auto fold = static_cast<std::size_t>(job) % kRtFolds;
        const auto [begin, end] = ranges[file];
        std::vector<std::size_t> rows;
        std::vector<int> labels;
        rows.reserve((end - begin) * (kRtFolds - 1) / kRtFolds);
        labels.reserve(rows.capacity());
        for (std::size_t i = begin; i < end; ++i) {
            if (folds[i] == fold) continue;
            rows.push_back(i);
            labels.push_back(data.rows[i].label);
        }
        auto matrix = make_matrix(
            data, fold_extra[fold], rows, begin, end);
        auto weights = initial_direction(matrix, data, rows, train_fdr);
        std::array<std::size_t, 2> previous{};
        for (unsigned int iteration = 0; iteration < max_iterations; ++iteration) {
            std::vector<double> scores(rows.size());
            for (std::size_t i = 0; i < rows.size(); ++i)
                scores[i] = score(matrix, rows[i], weights);
            const auto q = training_qvalues(scores, labels);
            weights = fit_svm(matrix, data, rows, q, train_fdr, c_pos, c_neg);
            for (std::size_t i = 0; i < rows.size(); ++i)
                scores[i] = score(matrix, rows[i], weights);
            const auto found = count_confident(scores, labels, train_fdr);
            ++result.iterations[file][fold];
            if (iteration >= 2 &&
                static_cast<double>(found) <= 1.01 * previous[iteration % 2]) {
                break;
            }
            previous[iteration % 2] = found;
        }
        std::vector<double> training_scores(rows.size());
        for (std::size_t i = 0; i < rows.size(); ++i) {
            training_scores[i] = score(matrix, rows[i], weights);
        }
        const auto [offset, scale] = training_score_calibration(
            training_scores, labels, train_fdr);
        auto& calibrated = result.calibrated_weights[file][fold];
        calibrated.resize(weights.size());
        for (std::size_t feature = 0; feature < matrix.columns; ++feature) {
            calibrated[feature] = weights[feature] / scale;
        }
        calibrated.back() = (weights.back() - offset) / scale;
        for (std::size_t i = begin; i < end; ++i) {
            if (folds[i] == fold) {
                result.scores[i] = (score(matrix, i, weights) - offset) / scale;
            }
        }
    }
    return result;
}

std::vector<double> SvmImplementation::local_error_probabilities(
    const std::vector<double>& scores, const std::vector<int>& labels) {
    std::vector<std::size_t> order(scores.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b];
    });
    std::vector<std::size_t> decoys(scores.size() + 1, 0);
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
        decoys[rank + 1] = decoys[rank] + (labels[order[rank]] == -1 ? 1 : 0);
    }
    std::vector<double> result(scores.size(), 1.0);
    constexpr std::size_t half_window = 250;
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
        const std::size_t lo = rank > half_window ? rank - half_window : 0;
        const std::size_t hi = std::min(order.size(), rank + half_window + 1);
        const double d = static_cast<double>(decoys[hi] - decoys[lo]);
        const double t = static_cast<double>((hi - lo) - (decoys[hi] - decoys[lo]));
        result[order[rank]] = std::min(1.0, (d + 1.0) / std::max(1.0, t));
    }
    return result;
}

std::vector<std::size_t> SvmRescorer::assign_folds(const Dataset& data) {
    return SvmImplementation::assign_folds(data);
}

SvmFit SvmRescorer::fit(
    const Dataset& data,
    const std::array<std::vector<double>, kRtFolds>& fold_extra,
    const std::vector<std::size_t>& folds, double train_fdr,
    unsigned int max_iterations, double c_pos, double c_neg) {
    return SvmImplementation::fit(
        data, fold_extra, folds, train_fdr, max_iterations, c_pos, c_neg);
}

std::vector<double> SvmRescorer::local_error_probabilities(
    const std::vector<double>& scores, const std::vector<int>& labels) {
    return SvmImplementation::local_error_probabilities(scores, labels);
}

} // namespace aerith
