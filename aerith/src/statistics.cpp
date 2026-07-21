#include "filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace aerith {

std::vector<double> target_decoy_qvalues(const std::vector<double>& scores,
                                         const std::vector<int>& labels) {
    if (scores.size() != labels.size()) {
        throw std::runtime_error("Score and label arrays have different sizes");
    }
    std::vector<std::size_t> order(scores.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b];
    });
    std::vector<double> q(scores.size(), 1.0);
    std::size_t target = 0;
    std::size_t decoy = 1;
    for (const auto index : order) {
        if (labels[index] == -1) {
            ++decoy;
        } else {
            ++target;
        }
        q[index] = std::min(1.0, static_cast<double>(decoy) /
                                    static_cast<double>(std::max<std::size_t>(1, target)));
    }
    double minimum = 1.0;
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        minimum = std::min(minimum, q[*it]);
        q[*it] = minimum;
    }
    return q;
}

class MixMaxEstimator final {
public:
    static std::vector<double> calculate(
        const std::vector<double>& scores, const std::vector<int>& labels,
        double* pi0_output);

private:
    struct RankedScore {
        double score = 0.0;
        bool target = false;
        std::size_t row = 0;
    };

    static std::vector<RankedScore> rank_scores(
        const std::vector<double>& scores, const std::vector<int>& labels);
    static std::vector<double> percolator_pvalues(
        const std::vector<RankedScore>& ranked);
    static double estimate_pi0(const std::vector<double>& pvalues);
    static void mixmax_tail_counts(
        const std::vector<RankedScore>& ranked,
        std::vector<double>& target_counts,
        std::vector<double>& decoy_counts);
};

std::vector<MixMaxEstimator::RankedScore> MixMaxEstimator::rank_scores(
    const std::vector<double>& scores, const std::vector<int>& labels) {
    std::vector<RankedScore> ranked;
    ranked.reserve(scores.size());
    for (std::size_t i = 0; i < scores.size(); ++i) {
        ranked.push_back({scores[i], labels[i] == 1, i});
    }
    std::stable_sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) {
        return a.score > b.score;
    });
    return ranked;
}

std::vector<double> MixMaxEstimator::percolator_pvalues(
    const std::vector<RankedScore>& ranked) {
    std::vector<double> pvalues;
    std::size_t targets = 0;
    for (const auto& item : ranked) {
        targets += item.target ? 1u : 0u;
    }
    pvalues.reserve(targets);
    std::size_t decoys = 1;
    for (std::size_t begin = 0; begin < ranked.size();) {
        std::size_t end = begin + 1;
        while (end < ranked.size() && ranked[end].score == ranked[begin].score) {
            ++end;
        }
        std::size_t tied_targets = 0;
        std::size_t tied_decoys = 0;
        for (std::size_t i = begin; i < end; ++i) {
            ranked[i].target ? ++tied_targets : ++tied_decoys;
        }
        for (std::size_t i = 0; i < tied_targets; ++i) {
            pvalues.push_back(static_cast<double>(decoys) +
                              static_cast<double>(tied_decoys * (i + 1)) /
                                  static_cast<double>(tied_targets + 1));
        }
        decoys += tied_decoys;
        begin = end;
    }
    const double denominator = static_cast<double>(decoys);
    for (double& value : pvalues) {
        value /= denominator;
    }
    return pvalues;
}

// Storey's bootstrap estimator as implemented by Percolator's
// PosteriorEstimator::estimatePi0. Sampling is deterministic for reproducibility.
double MixMaxEstimator::estimate_pi0(const std::vector<double>& pvalues) {
    if (pvalues.empty()) {
        return 1.0;
    }
    constexpr unsigned lambdas_count = 100;
    constexpr double max_lambda = 0.5;
    std::vector<double> lambdas;
    std::vector<double> estimates;
    for (unsigned i = 0; i <= lambdas_count; ++i) {
        const double lambda = (static_cast<double>(i + 1) / lambdas_count) * max_lambda;
        const auto start = std::lower_bound(pvalues.begin(), pvalues.end(), lambda);
        const double remaining = static_cast<double>(pvalues.end() - start);
        const double value = remaining / pvalues.size() / (1.0 - lambda);
        if (value > 0.0) {
            lambdas.push_back(lambda);
            estimates.push_back(value);
        }
    }
    if (estimates.empty()) {
        return 1.0;
    }
    const double minimum = *std::min_element(estimates.begin(), estimates.end());
    std::vector<double> mse(estimates.size(), 0.0);
    std::uint64_t state = 1;
    auto random_index = [&](std::size_t size) {
        state = (state * 279470273u) % 4294967291u;
        const double uniform = static_cast<double>(state) / 4294967292.0;
        return std::min(size - 1, static_cast<std::size_t>(uniform * size));
    };
    const std::size_t draws = std::min<std::size_t>(pvalues.size(), 1000);
    std::vector<double> sample(draws);
    for (unsigned boot = 0; boot < 100; ++boot) {
        for (double& value : sample) {
            value = pvalues[random_index(pvalues.size())];
        }
        std::sort(sample.begin(), sample.end());
        for (std::size_t i = 0; i < lambdas.size(); ++i) {
            const auto start = std::lower_bound(sample.begin(), sample.end(), lambdas[i]);
            const double value = static_cast<double>(sample.end() - start) /
                                 sample.size() / (1.0 - lambdas[i]);
            const double delta = value - minimum;
            mse[i] += delta * delta;
        }
    }
    const auto best = std::min_element(mse.begin(), mse.end()) - mse.begin();
    return std::clamp(estimates[static_cast<std::size_t>(best)], 0.0, 1.0);
}

void MixMaxEstimator::mixmax_tail_counts(
    const std::vector<RankedScore>& ranked,
    std::vector<double>& target_counts,
    std::vector<double>& decoy_counts) {
    std::size_t targets = 0;
    std::size_t decoys = 0;
    std::size_t queued_decoys = 0;
    for (auto it = ranked.rbegin(); it != ranked.rend(); ++it) {
        if (it->target) {
            ++targets;
        } else {
            ++decoys;
            ++queued_decoys;
        }
        const auto next = it + 1;
        if (next == ranked.rend() || next->score != it->score) {
            for (std::size_t i = 0; i < queued_decoys; ++i) {
                target_counts.push_back(static_cast<double>(targets));
                decoy_counts.push_back(static_cast<double>(decoys));
            }
            queued_decoys = 0;
        }
    }
}

std::vector<double> MixMaxEstimator::calculate(
    const std::vector<double>& scores, const std::vector<int>& labels,
    double* pi0_output) {
    if (scores.size() != labels.size()) {
        throw std::runtime_error("Score and label arrays have different sizes");
    }
    const auto ranked = rank_scores(scores, labels);
    const auto pvalues = percolator_pvalues(ranked);
    const double pi0 = estimate_pi0(pvalues);
    if (pi0_output != nullptr) {
        *pi0_output = pi0;
    }

    std::vector<double> tail_targets;
    std::vector<double> tail_decoys;
    if (pi0 < 1.0) {
        mixmax_tail_counts(ranked, tail_targets, tail_decoys);
    }

    std::vector<double> ranked_q(ranked.size(), 1.0);
    std::size_t decoys_seen = 1;
    std::size_t targets_seen = 0;
    double alternative_nulls = 0.0;
    for (std::size_t begin = 0; begin < ranked.size();) {
        std::size_t end = begin + 1;
        while (end < ranked.size() && ranked[end].score == ranked[begin].score) {
            ++end;
        }
        std::size_t tied_targets = 0;
        std::size_t tied_decoys = 0;
        for (std::size_t i = begin; i < end; ++i) {
            ranked[i].target ? ++tied_targets : ++tied_decoys;
        }
        targets_seen += tied_targets;
        decoys_seen += tied_decoys;
        if (pi0 < 1.0 && tied_decoys > 0) {
            const std::size_t j = tail_targets.size() - (decoys_seen - 1);
            const double count_targets = tail_targets.at(j);
            const double count_decoys = tail_decoys.at(j);
            double probability = (count_targets - pi0 * count_decoys) /
                                 ((1.0 - pi0) * count_decoys);
            probability = std::clamp(probability, 0.0, 1.0);
            alternative_nulls += tied_decoys * probability * (1.0 - pi0);
        }
        const double fdr = std::min(1.0, (decoys_seen * pi0 + alternative_nulls) /
                                                std::max<std::size_t>(1, targets_seen));
        for (std::size_t i = begin; i < end; ++i) {
            ranked_q[i] = fdr;
        }
        begin = end;
    }
    double minimum = 1.0;
    for (auto it = ranked_q.rbegin(); it != ranked_q.rend(); ++it) {
        minimum = std::min(minimum, *it);
        *it = minimum;
    }
    std::vector<double> result(scores.size(), 1.0);
    for (std::size_t rank = 0; rank < ranked.size(); ++rank) {
        result[ranked[rank].row] = ranked_q[rank];
    }
    return result;
}

std::vector<double> mixmax_qvalues(const std::vector<double>& scores,
                                   const std::vector<int>& labels,
                                   double* pi0_output) {
    return MixMaxEstimator::calculate(scores, labels, pi0_output);
}


} // namespace aerith
