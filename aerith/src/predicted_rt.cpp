#include "pipeline.hpp"
#include "prediction_cache.hpp"
#include "torch_device.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef AERITH_WITH_TORCH
#include <torch/script.h>
#endif

namespace aerith {

#ifndef AERITH_WITH_TORCH

struct RtPredictionLibrary::Impl {};

RtPredictionLibrary::RtPredictionLibrary()
    : impl_(std::make_unique<Impl>()) {}
RtPredictionLibrary::~RtPredictionLibrary() = default;
RtPredictionLibrary::RtPredictionLibrary(
    RtPredictionLibrary&&) noexcept = default;
RtPredictionLibrary& RtPredictionLibrary::operator=(
    RtPredictionLibrary&&) noexcept = default;
std::string RtPredictionLibrary::device() const { return {}; }
StageTiming RtPredictionLibrary::timing() const { return {}; }
StageTiming RtPredictionLibrary::cache_read_timing() const { return {}; }
StageTiming RtPredictionLibrary::cache_write_timing() const { return {}; }

RtPredictionLibrary PredictedRetentionTimeFeature::predict(
    const Config& config, const Dataset&) {
    if (config.predict_rt && !config.rt_model_path.empty()) {
        throw std::runtime_error(
            "Aerith was built without LibTorch; rebuild with Torch_DIR to predict RT");
    }
    return {};
}

void PredictedRetentionTimeFeature::add(
    const Config& config, Dataset& data,
    const RtPredictionLibrary&) {
    add(config, data);
}

void PredictedRetentionTimeFeature::add(const Config& config, Dataset&) {
    if (config.predict_rt && !config.rt_model_path.empty()) {
        throw std::runtime_error(
            "Aerith was built without LibTorch; rebuild with Torch_DIR to predict RT");
    }
}

#else

constexpr std::size_t kPredictionBatch = 4096;
constexpr std::size_t kRegressionSize = 5000;
constexpr std::size_t kGlobalRegressionSize = 4000;
constexpr std::size_t kRtSelectionBins = 50;
constexpr std::size_t kPerBinRegressionSize = 20;
constexpr std::size_t kMinimumLinearPoints = 10;
constexpr std::size_t kMinimumLoessPoints = 100;
constexpr std::size_t kRegressionSplits = 5;
constexpr int kRobustIterations = 2;
constexpr std::array<double, 4> kBandwidths{0.01, 0.05, 0.1, 0.2};

struct EncodedPeptide {
    std::string key;
    std::vector<std::int64_t> tokens;
};

struct RtPredictionLibrary::Impl {
    std::unordered_map<std::string, float> predictions;
    std::string device;
    StageTiming timing;
    StageTiming cache_read_timing;
    StageTiming cache_write_timing;
};

RtPredictionLibrary::RtPredictionLibrary()
    : impl_(std::make_unique<Impl>()) {}
RtPredictionLibrary::~RtPredictionLibrary() = default;
RtPredictionLibrary::RtPredictionLibrary(
    RtPredictionLibrary&&) noexcept = default;
RtPredictionLibrary& RtPredictionLibrary::operator=(
    RtPredictionLibrary&&) noexcept = default;
std::string RtPredictionLibrary::device() const {
    return impl_->device;
}
StageTiming RtPredictionLibrary::timing() const {
    return impl_->timing;
}
StageTiming RtPredictionLibrary::cache_read_timing() const {
    return impl_->cache_read_timing;
}
StageTiming RtPredictionLibrary::cache_write_timing() const {
    return impl_->cache_write_timing;
}

struct TrainingPoint {
    std::size_t row = 0;
    double observed = 0.0;
    double predicted = 0.0;
    double wdp = 0.0;
};

struct MonotoneCurve {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> tangent;

    double operator()(double value) const {
        if (x.empty()) throw std::runtime_error("Cannot evaluate an empty RT curve");
        if (x.size() == 1) return y.front();
        if (value <= x.front()) return y.front() + tangent.front() * (value - x.front());
        if (value >= x.back()) return y.back() + tangent.back() * (value - x.back());
        const auto upper = std::upper_bound(x.begin(), x.end(), value);
        const std::size_t right = static_cast<std::size_t>(upper - x.begin());
        const std::size_t left = right - 1;
        const double width = x[right] - x[left];
        const double t = (value - x[left]) / width;
        const double t2 = t * t;
        const double t3 = t2 * t;
        return (2.0 * t3 - 3.0 * t2 + 1.0) * y[left] +
               (t3 - 2.0 * t2 + t) * width * tangent[left] +
               (-2.0 * t3 + 3.0 * t2) * y[right] +
               (t3 - t2) * width * tangent[right];
    }
};

std::string rt_peptide_body(const std::string& peptide) {
    const auto first_bracket = peptide.find('[');
    const auto last_bracket = peptide.rfind(']');
    if (first_bracket != std::string::npos && last_bracket > first_bracket &&
        peptide.find('[', first_bracket + 1) == std::string::npos &&
        first_bracket <= 2 && peptide.size() - last_bracket <= 3) {
        return peptide.substr(first_bracket + 1, last_bracket - first_bracket - 1);
    }
    const auto first_dot = peptide.find('.');
    const auto last_dot = peptide.rfind('.');
    if (first_dot != std::string::npos && last_dot > first_dot) {
        return peptide.substr(first_dot + 1, last_dot - first_dot - 1);
    }
    return peptide;
}

std::int64_t residue_token(char residue) {
    switch (residue) {
    case 'G': return 3;  case 'A': return 4;  case 'V': return 5;
    case 'I': return 6;  case 'L': return 7;  case 'P': return 8;
    case 'F': return 9;  case 'W': return 10; case 'M': case 'X': return 11;
    case 'S': return 13; case 'T': return 14; case 'Y': return 15;
    case 'Q': return 16; case 'E': return 17; case 'N': return 18;
    case 'D': return 19; case 'K': case 'O': return 20;
    case 'R': return 21; case 'H': return 22; case 'C': case 'U': return 24;
    default: throw std::runtime_error(std::string("Unsupported DIA-NN RT residue: ") + residue);
    }
}

bool modification_symbol(char value) {
    switch (value) {
    case '~': case '!': case '@': case '>': case '<': case '%':
    case '^': case '&': case '*': case '(': case ')': case '/': case '$':
        return true;
    default:
        return false;
    }
}

bool numeric_modification(const std::string& body, std::size_t& position,
                          double& mass_shift) {
    if (position >= body.size() || body[position] != '[') return false;
    const auto close = body.find(']', position + 1);
    if (close == std::string::npos) return false;
    try {
        std::size_t consumed = 0;
        const std::string value = body.substr(position + 1, close - position - 1);
        mass_shift = std::stod(value, &consumed);
        if (consumed != value.size()) return false;
    } catch (const std::exception&) {
        return false;
    }
    position = close + 1;
    return true;
}

bool close_mass(double observed, double expected) {
    return std::abs(observed - expected) < 0.02;
}

EncodedPeptide encode_peptide(const Psm& psm) {
    const std::string body = rt_peptide_body(psm.peptide);
    EncodedPeptide encoded;
    encoded.key = body + '\x1f' + std::to_string(psm.charge);
    std::size_t position = 0;
    bool acetylated_nterm = !body.empty() && body.front() == '%';
    if (acetylated_nterm) {
        ++position;
    } else {
        double nterm_shift = 0.0;
        const std::size_t original = position;
        if (numeric_modification(body, position, nterm_shift)) {
            acetylated_nterm = close_mass(nterm_shift, 42.0106);
            if (!acetylated_nterm) position = original;
        }
    }
    encoded.tokens.push_back(acetylated_nterm ? 29 : 1);
    std::size_t residues = 0;
    while (position < body.size()) {
        const char residue = body[position++];
        if (!(residue >= 'A' && residue <= 'Z')) continue;
        char modification = '\0';
        if (position < body.size() && modification_symbol(body[position])) {
            modification = body[position++];
        }
        double numeric_shift = 0.0;
        const bool numeric =
            numeric_modification(body, position, numeric_shift);
        std::int64_t token = residue_token(residue);
        // This is the exact DIA-NN 2.6 dict.txt vocabulary used by rt.d0.pt.
        // Unknown search modifications, including deamidation, are stripped as
        // they are by DIA-NN/MSBooster with --strip-unknown-mods.
        if (residue == 'C' && modification != '(' &&
            (!numeric || close_mass(numeric_shift, 57.0215))) token = 25;
        else if (residue == 'M' &&
                 (modification == '~' || (numeric && close_mass(numeric_shift, 15.9949)))) {
            token = 26;
        } else if ((residue == 'S' || residue == 'T' || residue == 'Y') &&
                   (modification == '@' || modification == '>' || modification == '<' ||
                    (numeric && close_mass(numeric_shift, 79.9663)))) {
            token = residue == 'S' ? 31 : (residue == 'T' ? 32 : 33);
        } else if (residue == 'K' &&
                   (modification == '&' || (numeric && close_mass(numeric_shift, 28.0313)))) {
            token = 39;
        }
        encoded.tokens.push_back(token);
        ++residues;
    }
    encoded.tokens.push_back(2);
    if (residues < 5) {
        throw std::runtime_error("DIA-NN RT prediction requires at least five residues: " +
                                 psm.peptide);
    }
    return encoded;
}

std::vector<std::int64_t> diann_rt_tokens_for_testing(const Psm& psm) {
    return encode_peptide(psm).tokens;
}

std::string token_key(const std::vector<std::int64_t>& tokens) {
    std::string key;
    key.reserve(tokens.size());
    for (const auto token : tokens) key.push_back(static_cast<char>(token));
    return key;
}

#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 15
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
std::vector<float> predict_irt(const std::filesystem::path& model_path,
                               const Dataset& data,
                               std::string& selected_device) {
    struct UniquePeptide {
        std::vector<std::int64_t> tokens;
        std::vector<std::size_t> rows;
    };
    std::unordered_map<std::string, std::size_t> unique_index;
    std::unordered_map<std::string_view, std::size_t> peptide_index;
    std::vector<UniquePeptide> unique;
    unique.reserve(data.rows.size() / 2);
    peptide_index.reserve(data.rows.size() / 2);
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const std::string_view peptide = data.rows[row].peptide;
        const auto cached = peptide_index.find(peptide);
        if (cached != peptide_index.end()) {
            unique[cached->second].rows.push_back(row);
            continue;
        }
        auto encoded = encode_peptide(data.rows[row]);
        const std::string key = token_key(encoded.tokens);
        const auto inserted = unique_index.emplace(key, unique.size());
        if (inserted.second) unique.push_back({std::move(encoded.tokens), {}});
        const std::size_t index = inserted.first->second;
        peptide_index.emplace(peptide, index);
        unique[index].rows.push_back(row);
    }
    std::unordered_map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < unique.size(); ++i) {
        groups[unique[i].tokens.size()].push_back(i);
    }

    // DIA-NN 2.6.1 Academia models/rt.d0.pt. The checkpoint was renamed to
    // diann-2.6.1-retention-time.pt. Source release and license:
    // https://github.com/vdemichev/DiaNN/releases/tag/2.0
    return run_torch_prefer_cuda(
        "DIA-NN RT prediction", selected_device,
        [&](const torch::Device& device) {
        auto model = load_torch_model_on_device(model_path, device);
        c10::InferenceMode inference_mode;
        std::vector<float> predictions(data.rows.size());
        for (auto& group_entry : groups) {
            const std::size_t token_count = group_entry.first;
            auto& group = group_entry.second;
            for (std::size_t begin = 0; begin < group.size(); begin += kPredictionBatch) {
                const std::size_t count = std::min(kPredictionBatch, group.size() - begin);
                auto cpu_input = torch::empty(
                    {static_cast<long>(count), static_cast<long>(token_count + 1)},
                    torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU));
                auto accessor = cpu_input.accessor<std::int64_t, 2>();
                for (std::size_t local = 0; local < count; ++local) {
                    const auto& peptide = unique[group[begin + local]];
                    // The RT graph discards this leading charge column, but DIA-NN
                    // still supplies it to the shared peptide-model interface.
                    accessor[local][0] = 0;
                    for (std::size_t column = 0; column < token_count; ++column) {
                        accessor[local][column + 1] = peptide.tokens[column];
                    }
                }
                const auto input = move_torch_input(cpu_input, device);
                const auto output =
                    model.forward({input}).toTensor().to(torch::kCPU).contiguous();
                if (static_cast<std::size_t>(output.numel()) != count) {
                    throw std::runtime_error("Unexpected DIA-NN RT output shape");
                }
                const float* values = output.data_ptr<float>();
                for (std::size_t local = 0; local < count; ++local) {
                    // DIA-NN performs this affine conversion outside Torch. The
                    // double intermediate and single float cast match its binary.
                    const float irt = static_cast<float>(
                        static_cast<double>(values[local]) * 250.0 - 60.0);
                    for (const auto row : unique[group[begin + local]].rows) {
                        predictions[row] = irt;
                    }
                }
            }
        }
        return predictions;
    });
}
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 15
#pragma GCC diagnostic pop
#endif

struct RtCacheLoad {
    bool compatible = false;
    std::size_t hits = 0;
    std::unordered_set<std::string> exact_keys;
};

RtCacheLoad load_rt_library_cache(
    const Config& config, const Dataset& unique_peptides,
    RtPredictionLibrary::Impl& library) {
    const auto timing_begin = std::chrono::steady_clock::now();
    const auto timing_cpu_begin = std::clock();
    RtCacheLoad result;
    auto cache = read_prediction_cache(config, false);
    result.compatible = cache.compatible;
    if (cache.compatible) {
        result.exact_keys.reserve(unique_peptides.rows.size());
        for (const auto& psm : unique_peptides.rows) {
            const auto exact = prediction_cache_key(psm.peptide, psm.charge);
            const auto found = cache.entries.find(exact);
            if (found == cache.entries.end() || !found->second.has_rt ||
                !result.exact_keys.insert(exact).second) {
                continue;
            }
            library.predictions.emplace(
                token_key(encode_peptide(psm).tokens), found->second.rt);
            ++result.hits;
        }
    }
    if (result.hits != 0) library.device = "cache";
    library.timing = {};
    library.cache_read_timing = {
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - timing_begin).count(),
        static_cast<double>(std::clock() - timing_cpu_begin) /
            CLOCKS_PER_SEC};
    return result;
}

StageTiming append_rt_library_cache(
    const Config& config,
    const Dataset& unique_peptides,
    const std::unordered_map<std::string, float>& predictions,
    const RtCacheLoad& loaded) {
    const auto timing_begin = std::chrono::steady_clock::now();
    const auto timing_cpu_begin = std::clock();
    std::unordered_map<std::string, PredictionCacheEntry> updates;
    updates.reserve(unique_peptides.rows.size());
    for (const auto& psm : unique_peptides.rows) {
        if (psm.label != 1) continue;
        const auto exact = prediction_cache_key(psm.peptide, psm.charge);
        if (loaded.exact_keys.find(exact) != loaded.exact_keys.end()) continue;
        const auto prediction = predictions.find(
            token_key(encode_peptide(psm).tokens));
        if (prediction == predictions.end()) continue;
        auto& update = updates[exact];
        update.has_rt = true;
        update.rt = prediction->second;
    }
    update_prediction_cache(config, updates);
    return {
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - timing_begin).count(),
        static_cast<double>(std::clock() - timing_cpu_begin) /
            CLOCKS_PER_SEC};
}

double median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 == 0
        ? (values[middle - 1] + values[middle]) / 2.0 : values[middle];
}

double tricube(double value) {
    const double x = std::min(1.0, std::abs(value));
    const double cube = x * x * x;
    const double one_minus = 1.0 - cube;
    return one_minus * one_minus * one_minus;
}

std::vector<double> loess_smooth(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const std::vector<double>& user_weights,
                                 double bandwidth, int robust_iterations) {
    const std::size_t n = x.size();
    if (n != y.size() || n != user_weights.size() || n < 2) {
        throw std::runtime_error("Invalid LOESS input");
    }
    const std::size_t bandwidth_points = static_cast<std::size_t>(bandwidth * n);
    if (bandwidth_points < 2) throw std::runtime_error("LOESS bandwidth is too small");
    std::vector<double> fitted(n), residual(n), robustness(n, 1.0);
    constexpr double accuracy = 1e-12;
    for (int iteration = 0; iteration <= robust_iterations; ++iteration) {
        std::size_t left = 0;
        std::size_t right = bandwidth_points - 1;
        for (std::size_t i = 0; i < n; ++i) {
            if (i > 0 && right + 1 < n &&
                x[right + 1] - x[i] < x[i] - x[left]) {
                ++left;
                ++right;
            }
            const std::size_t edge = x[i] - x[left] > x[right] - x[i] ? left : right;
            const double edge_distance = std::abs(x[edge] - x[i]);
            const double denominator = edge_distance > 0.0 ? 1.0 / edge_distance : 0.0;
            double sum_weights = 0.0, sum_x = 0.0, sum_x2 = 0.0;
            double sum_y = 0.0, sum_xy = 0.0;
            for (std::size_t k = left; k <= right; ++k) {
                const double distance = std::abs(x[k] - x[i]);
                const double weight = tricube(distance * denominator) *
                                      robustness[k] * user_weights[k];
                const double xw = x[k] * weight;
                sum_weights += weight;
                sum_x += xw;
                sum_x2 += x[k] * xw;
                sum_y += y[k] * weight;
                sum_xy += y[k] * xw;
            }
            if (!(sum_weights > 0.0)) throw std::runtime_error("Zero LOESS weight");
            const double mean_x = sum_x / sum_weights;
            const double mean_y = sum_y / sum_weights;
            const double mean_x2 = sum_x2 / sum_weights;
            const double mean_xy = sum_xy / sum_weights;
            const double variance = mean_x2 - mean_x * mean_x;
            const double beta = std::sqrt(std::abs(variance)) < accuracy
                ? 0.0 : (mean_xy - mean_x * mean_y) / variance;
            const double alpha = mean_y - beta * mean_x;
            fitted[i] = beta * x[i] + alpha;
            residual[i] = std::abs(y[i] - fitted[i]);
        }
        if (iteration == robust_iterations) break;
        const double median_residual = median(residual);
        if (std::abs(median_residual) < accuracy) break;
        for (std::size_t i = 0; i < n; ++i) {
            const double argument = residual[i] / (6.0 * median_residual);
            if (argument >= 1.0) robustness[i] = 0.0;
            else {
                const double weight = 1.0 - argument * argument;
                robustness[i] = weight * weight;
            }
        }
    }
    return fitted;
}

MonotoneCurve isotonic_curve(const std::vector<double>& x,
                              const std::vector<double>& y) {
    struct Block { double x_sum = 0.0; double y_sum = 0.0; std::size_t count = 0; };
    std::vector<Block> blocks;
    blocks.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        blocks.push_back({x[i], y[i], 1});
        while (blocks.size() > 1) {
            const auto& a = blocks[blocks.size() - 2];
            const auto& b = blocks.back();
            // pairAdjacentViolators 1.4.16 merges equal adjacent levels too.
            if (a.y_sum / a.count < b.y_sum / b.count) break;
            Block merged{a.x_sum + b.x_sum, a.y_sum + b.y_sum, a.count + b.count};
            blocks.pop_back();
            blocks.back() = merged;
        }
    }
    MonotoneCurve curve;
    curve.x.reserve(blocks.size());
    curve.y.reserve(blocks.size());
    for (const auto& block : blocks) {
        curve.x.push_back(block.x_sum / block.count);
        curve.y.push_back(block.y_sum / block.count);
    }
    if (curve.x.size() == 1) {
        curve.tangent.push_back(0.0);
        return curve;
    }
    const std::size_t segments = curve.x.size() - 1;
    std::vector<double> secant(segments);
    for (std::size_t i = 0; i < segments; ++i) {
        secant[i] = (curve.y[i + 1] - curve.y[i]) / (curve.x[i + 1] - curve.x[i]);
    }
    curve.tangent.resize(curve.x.size());
    curve.tangent.front() = secant.front();
    curve.tangent.back() = secant.back();
    for (std::size_t i = 1; i < segments; ++i) {
        curve.tangent[i] = (secant[i - 1] + secant[i]) / 2.0;
    }
    for (std::size_t i = 0; i < segments; ++i) {
        if (secant[i] == 0.0) {
            curve.tangent[i] = 0.0;
            curve.tangent[i + 1] = 0.0;
            continue;
        }
        const double alpha = curve.tangent[i] / secant[i];
        const double beta = curve.tangent[i + 1] / secant[i];
        const double norm = alpha * alpha + beta * beta;
        if (norm > 9.0) {
            const double scale = 3.0 / std::sqrt(norm);
            curve.tangent[i] = scale * alpha * secant[i];
            curve.tangent[i + 1] = scale * beta * secant[i];
        }
    }
    return curve;
}

MonotoneCurve fit_curve(std::vector<double> x, std::vector<double> y,
                        double bandwidth) {
    std::vector<std::size_t> order(x.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return x[a] < x[b];
    });
    std::vector<double> sorted_x(x.size()), sorted_y(y.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        sorted_x[i] = x[order[i]];
        sorted_y[i] = y[order[i]];
    }
    double previous_original = -1.0;
    for (std::size_t i = 0; i < sorted_x.size(); ++i) {
        const double original = sorted_x[i];
        if (i != 0 && original == previous_original) sorted_x[i] = sorted_x[i - 1] + 1e-8;
        else previous_original = original;
    }
    if (sorted_x.size() < kMinimumLoessPoints || bandwidth > 1.0) bandwidth = 1.0;
    else if (bandwidth < 2.0 / sorted_x.size()) bandwidth = 3.0 / sorted_x.size();

    std::vector<double> unit_weights(sorted_x.size(), 1.0);
    auto fitted = loess_smooth(sorted_x, sorted_y, unit_weights,
                               bandwidth, kRobustIterations);
    if (fitted.size() > 100) {
        std::vector<double> weights(fitted.size());
        const std::size_t radius = fitted.size() / 100;
        for (std::size_t i = 0; i < fitted.size(); ++i) {
            const std::size_t begin = i > radius ? i - radius : 0;
            const std::size_t end = std::min(fitted.size(), i + radius);
            weights[i] = std::abs(median(std::vector<double>(
                fitted.begin() + static_cast<std::ptrdiff_t>(begin),
                fitted.begin() + static_cast<std::ptrdiff_t>(end))));
        }
        fitted = loess_smooth(sorted_x, sorted_y, weights,
                              bandwidth, kRobustIterations);
    }
    auto isotonic = isotonic_curve(sorted_x, fitted);
    for (std::size_t i = 0; i < fitted.size(); ++i) fitted[i] = isotonic(sorted_x[i]);

    if (fitted.size() > 100) {
        std::vector<double> residual(fitted.size());
        for (std::size_t i = 0; i < fitted.size(); ++i) residual[i] = fitted[i] - sorted_y[i];
        const double mean = std::accumulate(residual.begin(), residual.end(), 0.0) /
                            residual.size();
        double variance = 0.0;
        for (const double value : residual) variance += (value - mean) * (value - mean);
        variance /= std::max<std::size_t>(1, residual.size() - 1);
        const double standard_deviation = std::sqrt(variance);
        std::vector<double> kept_x, kept_y, kept_residual;
        kept_x.reserve(sorted_x.size());
        kept_y.reserve(sorted_y.size());
        kept_residual.reserve(residual.size());
        for (std::size_t i = 0; i < residual.size(); ++i) {
            if (standard_deviation > 0.0 &&
                std::abs((residual[i] - mean) / standard_deviation) > 2.0) continue;
            kept_x.push_back(sorted_x[i]);
            kept_y.push_back(sorted_y[i]);
            kept_residual.push_back(residual[i]);
        }
        std::vector<double> weights(kept_x.size());
        const std::size_t radius = kept_x.size() / 100;
        for (std::size_t i = 0; i < kept_x.size(); ++i) {
            const std::size_t begin = i > radius ? i - radius : 0;
            const std::size_t end = std::min(kept_x.size(), i + radius);
            weights[i] = std::abs(median(std::vector<double>(
                kept_residual.begin() + static_cast<std::ptrdiff_t>(begin),
                kept_residual.begin() + static_cast<std::ptrdiff_t>(end))));
        }
        fitted = loess_smooth(kept_x, kept_y, weights,
                              bandwidth, kRobustIterations);
        return isotonic_curve(kept_x, fitted);
    }
    return isotonic;
}

MonotoneCurve select_curve(const std::vector<TrainingPoint>& points) {
    if (points.size() < kMinimumLinearPoints) {
        throw std::runtime_error("DIA-NN delta-RT requires at least ten calibration PSMs");
    }
    if (points.size() < kMinimumLoessPoints) {
        double mean_x = 0.0, mean_y = 0.0;
        for (const auto& point : points) { mean_x += point.observed; mean_y += point.predicted; }
        mean_x /= points.size();
        mean_y /= points.size();
        double numerator = 0.0, denominator = 0.0;
        for (const auto& point : points) {
            numerator += (point.observed - mean_x) * (point.predicted - mean_y);
            denominator += (point.observed - mean_x) * (point.observed - mean_x);
        }
        const double slope = denominator > 0.0 ? numerator / denominator : 0.0;
        MonotoneCurve curve;
        curve.x = {points.front().observed, points.front().observed + 1.0};
        curve.y = {mean_y + slope * (curve.x[0] - mean_x),
                   mean_y + slope * (curve.x[1] - mean_x)};
        curve.tangent = {slope, slope};
        return curve;
    }

    constexpr std::size_t candidate_count =
        kRegressionSplits * kBandwidths.size();
    std::array<double, candidate_count> candidate_mse;
    candidate_mse.fill(std::numeric_limits<double>::infinity());
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t task = 0;
         task < static_cast<std::ptrdiff_t>(candidate_count); ++task) {
        const std::size_t candidate = static_cast<std::size_t>(task);
        const std::size_t split = candidate / kBandwidths.size();
        const double bandwidth = kBandwidths[candidate % kBandwidths.size()];
        std::vector<double> train_x, train_y;
        train_x.reserve(points.size() - points.size() / kRegressionSplits);
        train_y.reserve(points.size() - points.size() / kRegressionSplits);
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i % kRegressionSplits == split) continue;
            train_x.push_back(points[i].observed);
            train_y.push_back(points[i].predicted);
        }
        try {
            const auto curve = fit_curve(train_x, train_y, bandwidth);
            double mse = 0.0;
            std::size_t count = 0;
            for (std::size_t i = split; i < points.size();
                 i += kRegressionSplits) {
                const double difference =
                    curve(points[i].observed) - points[i].predicted;
                mse += difference * difference;
                ++count;
            }
            candidate_mse[candidate] = mse / count;
        } catch (const std::exception&) {
        }
    }
    std::array<double, kRegressionSplits> best{};
    for (std::size_t split = 0; split < kRegressionSplits; ++split) {
        double best_mse = std::numeric_limits<double>::infinity();
        double best_bandwidth = 1.0;
        for (std::size_t bandwidth = 0; bandwidth < kBandwidths.size();
             ++bandwidth) {
            const double mse =
                candidate_mse[split * kBandwidths.size() + bandwidth];
            if (mse < best_mse) {
                best_mse = mse;
                best_bandwidth = kBandwidths[bandwidth];
            }
        }
        best[split] = best_bandwidth;
    }
    double bandwidth = std::accumulate(best.begin(), best.end(), 0.0) / best.size();
    bandwidth = std::round(bandwidth * 10000.0) / 10000.0;
    std::vector<double> x, y;
    x.reserve(points.size());
    y.reserve(points.size());
    for (const auto& point : points) { x.push_back(point.observed); y.push_back(point.predicted); }
    while (true) {
        try { return fit_curve(x, y, bandwidth); }
        catch (const std::exception&) {
            if (bandwidth >= 1.0) throw;
            bandwidth = std::min(1.0, bandwidth * 2.0);
        }
    }
}

std::vector<TrainingPoint> calibration_points(const Dataset& data,
                                               const std::vector<float>& predicted,
                                               const std::vector<std::size_t>& file_rows) {
    const auto wdp_column = std::find(
        data.feature_names.begin(), data.feature_names.end(), "WDPscores");
    if (wdp_column == data.feature_names.end()) {
        throw std::runtime_error("DIA-NN delta-RT calibration requires WDPscores");
    }
    const std::size_t wdp_index = static_cast<std::size_t>(
        wdp_column - data.feature_names.begin());

    std::vector<std::size_t> order;
    double minimum_rt = std::numeric_limits<double>::infinity();
    double maximum_rt = -std::numeric_limits<double>::infinity();
    order.reserve(file_rows.size());
    for (const auto i : file_rows) {
        if (!std::isfinite(data.rows[i].retention)) continue;
        order.push_back(i);
        minimum_rt = std::min(minimum_rt, data.rows[i].retention);
        maximum_rt = std::max(maximum_rt, data.rows[i].retention);
    }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (data.rows[a].scan != data.rows[b].scan) return data.rows[a].scan < data.rows[b].scan;
        return data.rows[a].rank < data.rows[b].rank;
    });

    std::unordered_map<std::string, std::vector<TrainingPoint>> by_precursor;
    for (const auto i : order) {
        if (predicted[i] == 0.0f) continue;
        const std::string key = rt_peptide_body(data.rows[i].peptide) + '\x1f' +
                                std::to_string(data.rows[i].charge);
        by_precursor[key].push_back({
            i, data.rows[i].retention, predicted[i], data.rows[i].features[wdp_index]});
    }

    std::vector<TrainingPoint> candidates;
    candidates.reserve(by_precursor.size());
    for (auto& entry : by_precursor) {
        auto& matches = entry.second;
        const double best = std::max_element(
            matches.begin(), matches.end(), [](const auto& a, const auto& b) {
                return a.wdp < b.wdp;
            })->wdp;
        std::vector<TrainingPoint> tied;
        for (const auto& point : matches) if (point.wdp == best) tied.push_back(point);
        candidates.push_back(tied[tied.size() / 2]);
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        if (a.wdp != b.wdp) return a.wdp > b.wdp;
        return a.row < b.row;
    });
    if (candidates.size() <= kRegressionSize) return candidates;

    std::vector<bool> selected(candidates.size(), false);
    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(kRegressionSize);
    for (std::size_t i = 0; i < kGlobalRegressionSize; ++i) {
        selected[i] = true;
        selected_indices.push_back(i);
    }

    std::array<std::vector<std::size_t>, kRtSelectionBins> bins;
    const double rt_width = maximum_rt - minimum_rt;
    for (std::size_t i = kGlobalRegressionSize; i < candidates.size(); ++i) {
        std::size_t bin = 0;
        if (rt_width > 0.0) {
            bin = std::min(
                kRtSelectionBins - 1,
                static_cast<std::size_t>(
                    (candidates[i].observed - minimum_rt) / rt_width * kRtSelectionBins));
        }
        bins[bin].push_back(i);
    }
    for (const auto& bin : bins) {
        const std::size_t count = std::min(kPerBinRegressionSize, bin.size());
        for (std::size_t i = 0; i < count; ++i) {
            selected[bin[i]] = true;
            selected_indices.push_back(bin[i]);
        }
    }

    // Sparse RT bins may contain fewer than 20 remaining precursors. Preserve
    // the requested 5,000-point regression size by taking the best remaining
    // WDP candidates after every bin has had its opportunity to contribute.
    for (std::size_t i = kGlobalRegressionSize;
         selected_indices.size() < kRegressionSize && i < candidates.size(); ++i) {
        if (selected[i]) continue;
        selected[i] = true;
        selected_indices.push_back(i);
    }
    std::sort(selected_indices.begin(), selected_indices.end());
    std::vector<TrainingPoint> result;
    result.reserve(selected_indices.size());
    for (const auto index : selected_indices) result.push_back(candidates[index]);
    return result;
}

double inverse_lookup(const std::vector<double>& predicted,
                      const std::vector<double>& experimental,
                      double target) {
    const auto upper = std::lower_bound(predicted.begin(), predicted.end(), target);
    if (upper == predicted.begin()) return experimental.front();
    if (upper == predicted.end()) return experimental.back();
    const std::size_t right =
        static_cast<std::size_t>(upper - predicted.begin());
    if (*upper == target) return experimental[right];
    const std::size_t left = right - 1;
    const double ratio = (target - predicted[left]) /
                         (predicted[right] - predicted[left]);
    return experimental[left] +
           ratio * (experimental[right] - experimental[left]);
}

float rounded_feature(double value) {
    return static_cast<float>(std::round(value * 10000.0) / 10000.0);
}

RtPredictionLibrary PredictedRetentionTimeFeature::predict(
    const Config& config, const Dataset& unique_peptides) {
    RtPredictionLibrary library;
    if (!config.predict_rt) return library;
    if (config.rt_model_path.empty()) {
        throw std::runtime_error(
            "Cannot locate the DIA-NN RT model; use --rt-model");
    }
    const std::filesystem::path model_path(config.rt_model_path);
    if (!std::filesystem::is_regular_file(model_path)) {
        throw std::runtime_error("DIA-NN RT model does not exist: " + model_path.string());
    }
    const auto loaded = load_rt_library_cache(
        config, unique_peptides, *library.impl_);
    Dataset missing;
    std::unordered_map<std::string, std::size_t> missing_keys;
    missing_keys.reserve(unique_peptides.rows.size());
    for (const auto& psm : unique_peptides.rows) {
        const auto key = token_key(encode_peptide(psm).tokens);
        if (library.impl_->predictions.find(key) ==
                library.impl_->predictions.end()) {
            upsert_prediction_exemplar(missing, missing_keys, key, psm);
        }
    }
    if (!config.prediction_cache_path.empty()) {
        std::unordered_map<std::string, int> requested_labels;
        requested_labels.reserve(unique_peptides.rows.size());
        for (const auto& psm : unique_peptides.rows) {
            const auto key = prediction_cache_key(
                psm.peptide, psm.charge);
            const auto [found, inserted] = requested_labels.emplace(
                key, psm.label);
            if (!inserted && psm.label == 1) found->second = 1;
        }
        std::size_t target_reused = 0;
        std::size_t target_predict = 0;
        std::size_t decoy_reused = 0;
        std::size_t decoy_predict = 0;
        for (const auto& [key, label] : requested_labels) {
            const bool reused =
                loaded.exact_keys.find(key) != loaded.exact_keys.end();
            if (label == 1) {
                reused ? ++target_reused : ++target_predict;
            } else {
                reused ? ++decoy_reused : ++decoy_predict;
            }
        }
        std::cerr << "DIA-NN RT cache: target reused="
                  << target_reused << ", predict=" << target_predict
                  << "; decoy reused=" << decoy_reused
                  << ", predict=" << decoy_predict << '\n';
    }
    if (missing.rows.empty()) {
        library.impl_->cache_write_timing = append_rt_library_cache(
            config, unique_peptides,
            library.impl_->predictions, loaded);
        return library;
    }
    const auto prediction_begin = std::chrono::steady_clock::now();
    const std::clock_t prediction_cpu_begin = std::clock();
    std::string prediction_device;
    const auto predicted = predict_irt(
        model_path, missing, prediction_device);
    const std::clock_t prediction_cpu_end = std::clock();
    const auto prediction_end = std::chrono::steady_clock::now();
    library.impl_->timing = {
        std::chrono::duration<double>(prediction_end - prediction_begin).count(),
        static_cast<double>(prediction_cpu_end - prediction_cpu_begin) /
            CLOCKS_PER_SEC};
    std::unordered_map<std::string, float> inferred;
    inferred.reserve(missing.rows.size());
    for (std::size_t row = 0; row < missing.rows.size(); ++row) {
        const auto encoded = encode_peptide(missing.rows[row]);
        inferred.emplace(
            token_key(encoded.tokens), predicted[row]);
    }
    library.impl_->predictions.merge(inferred);
    library.impl_->cache_write_timing = append_rt_library_cache(
        config, unique_peptides, library.impl_->predictions, loaded);
    library.impl_->device = loaded.hits == 0
        ? prediction_device : "cache+" + prediction_device;
    return library;
}

void PredictedRetentionTimeFeature::add(
    const Config& config, Dataset& data,
    const RtPredictionLibrary& library) {
    if (!config.predict_rt) return;
    constexpr std::array<const char*, 3> names{
        "delta_RT_loess", "delta_RT_loess_real", "pred_RT_real_units"};
    for (const char* name : names) {
        if (std::find(data.feature_names.begin(), data.feature_names.end(), name) !=
            data.feature_names.end()) {
            throw std::runtime_error(std::string("PIN already contains ") + name);
        }
    }
    std::vector<float> predicted(data.rows.size());
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto encoded = encode_peptide(data.rows[row]);
        const auto found = library.impl_->predictions.find(
            token_key(encoded.tokens));
        if (found == library.impl_->predictions.end()) {
            throw std::runtime_error(
                "Missing DIA-NN RT prediction for " +
                data.rows[row].peptide);
        }
        predicted[row] = found->second;
    }
    data.rt_prediction_device = library.impl_->device;
    data.rt_prediction_timing = library.impl_->timing;
    data.rt_cache_read_timing = library.impl_->cache_read_timing;
    data.rt_cache_write_timing = library.impl_->cache_write_timing;
    std::vector<std::array<float, 3>> features(data.rows.size());
    std::vector<std::vector<std::size_t>> rows_by_file(data.input_paths.size());
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        rows_by_file[data.rows[i].file_id].push_back(i);
    }
    for (std::size_t file = 0; file < data.input_paths.size(); ++file) {
        auto points = calibration_points(data, predicted, rows_by_file[file]);
        const auto curve = select_curve(points);
        double minimum_rt = std::numeric_limits<double>::infinity();
        double maximum_rt = -std::numeric_limits<double>::infinity();
        for (const auto index : rows_by_file[file]) {
            minimum_rt = std::min(minimum_rt, data.rows[index].retention);
            maximum_rt = std::max(maximum_rt, data.rows[index].retention);
        }
        constexpr std::size_t increments = 10000;
        const double increment = (maximum_rt - minimum_rt) / increments;
        std::vector<double> inverse_predicted;
        std::vector<double> inverse_experimental;
        inverse_predicted.reserve(increments);
        inverse_experimental.reserve(increments);
        for (std::size_t i = 0; i < increments; ++i) {
            const double experimental = minimum_rt + i * increment;
            const double predicted_value = curve(experimental);
            if (!inverse_predicted.empty() &&
                predicted_value == inverse_predicted.back()) {
                inverse_experimental.back() = experimental;
            } else {
                inverse_predicted.push_back(predicted_value);
                inverse_experimental.push_back(experimental);
            }
        }
        #pragma omp parallel for schedule(static)
        for (std::ptrdiff_t index = 0;
             index < static_cast<std::ptrdiff_t>(rows_by_file[file].size()); ++index) {
            const auto i =
                rows_by_file[file][static_cast<std::size_t>(index)];
            const double calibrated = curve(data.rows[i].retention);
            const double predicted_minutes = inverse_lookup(
                inverse_predicted, inverse_experimental, predicted[i]);
            features[i] = {
                rounded_feature(std::abs(calibrated - predicted[i])),
                rounded_feature(std::abs(data.rows[i].retention - predicted_minutes)),
                rounded_feature(predicted_minutes)};
        }
    }
    data.feature_names.emplace_back(names[0]);
    data.generated_feature_names.emplace_back(names[0]);
    data.has_predicted_rt_diagnostics = true;
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        data.rows[i].features.push_back(features[i][0]);
        data.rows[i].delta_rt_loess_real = features[i][1];
        data.rows[i].predicted_rt_real_units = features[i][2];
    }
}

void PredictedRetentionTimeFeature::add(
    const Config& config, Dataset& data) {
    if (!config.predict_rt) return;
    auto predictions = predict(config, data);
    add(config, data, predictions);
}

#endif

} // namespace aerith
