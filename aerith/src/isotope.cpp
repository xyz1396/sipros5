#include "isotope.hpp"

#include "filter.hpp"
#include "isotopologue.h"
#include "proNovoConfig.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

namespace aerith {

bool initialized = false;
char sip_atom = '\0';
int sip_mass_number = 0;
int sip_isotope_index = -1;
std::unique_ptr<Isotopologue> pristine_isotopologue;
std::atomic<std::uint64_t> model_generation{0};

std::pair<char, int> parse_sip_isotope(std::string value) {
    value.erase(
        std::remove_if(value.begin(), value.end(), [](unsigned char ch) {
            return std::isspace(ch) != 0;
        }),
        value.end());
    std::transform(
        value.begin(), value.end(), value.begin(), [](unsigned char ch) {
            return static_cast<char>(std::toupper(ch));
        });
    if (value == "C13") return {'C', 13};
    if (value == "H2") return {'H', 2};
    if (value == "N15") return {'N', 15};
    if (value == "O18") return {'O', 18};
    if (value == "S34") return {'S', 34};
    throw std::runtime_error(
        "Unsupported Aerith SIP isotope '" + value +
        "'; use C13,H2,N15,O18,S34");
}

std::string isotope_peptide(const std::string& peptide) {
    const auto open = peptide.find('[');
    const auto close = peptide.find(']', open == std::string::npos
        ? 0 : open + 1);
    if (open != std::string::npos && close != std::string::npos &&
        close > open) {
        return peptide.substr(open, close - open + 1);
    }
    return '[' + peptide + ']';
}

Isotopologue& enriched_isotopologue(double abundance_pct) {
    if (!initialized) {
        throw std::runtime_error(
            "SIP isotope model was used before initialization");
    }
    if (!std::isfinite(abundance_pct) ||
        abundance_pct < 0.0 || abundance_pct > 100.0) {
        throw std::runtime_error(
            "SIP isotope abundance must be within [0,100]");
    }
    const auto key = static_cast<long long>(
        std::llround(abundance_pct * 1000000.0));
    thread_local std::uint64_t local_generation = 0;
    thread_local std::unordered_map<
        long long, std::unique_ptr<Isotopologue>> models;
    const auto current_generation = model_generation.load();
    if (local_generation != current_generation) {
        models.clear();
        local_generation = current_generation;
    }
    const auto found = models.find(key);
    if (found != models.end()) return *found->second;
    auto model = std::make_unique<Isotopologue>(*pristine_isotopologue);
    ProNovoConfig::setSipAbundance(
        *model, sip_atom, sip_isotope_index, abundance_pct);
    auto inserted = models.emplace(key, std::move(model));
    return *inserted.first->second;
}

void initialize_sip_isotope_model(const Config& config) {
    initialized = false;
    pristine_isotopologue.reset();
    sip_atom = '\0';
    sip_mass_number = 0;
    sip_isotope_index = -1;
    if (config.sip_isotope.empty()) return;

    const auto parsed = parse_sip_isotope(config.sip_isotope);
    if (!ProNovoConfig::load(ProNovoConfig::Profile::Sip)) {
        throw std::runtime_error(
            "Could not initialize Sipros SIP chemistry in Aerith");
    }
    std::string error;
    std::vector<std::string> fixed = config.fixed_ptm_selectors;
    if (fixed.empty() && !config.fixed_cam) fixed.push_back("none");
    if (!ProNovoConfig::configureFixedPtms(fixed, error)) {
        throw std::runtime_error(
            "Could not configure Aerith fixed PTMs: " + error);
    }
    if (!ProNovoConfig::configureVariablePtms(
            config.ptm_selectors, config.max_ptm_count, error)) {
        throw std::runtime_error(
            "Could not configure Aerith variable PTMs: " + error);
    }
    if (!ProNovoConfig::selectSipTarget(
            parsed.first, parsed.second, error)) {
        throw std::runtime_error(
            "Could not configure Aerith SIP isotope: " + error);
    }

    sip_atom = parsed.first;
    sip_mass_number = parsed.second;
    sip_isotope_index = ProNovoConfig::resolveSipIsotopeIndex(
        ProNovoConfig::configIsotopologue, sip_atom, sip_mass_number);
    pristine_isotopologue = std::make_unique<Isotopologue>(
        ProNovoConfig::configIsotopologue);
    ++model_generation;
    initialized = true;
}

bool sip_isotope_model_enabled() {
    return initialized;
}

std::vector<TheoreticalIsotopePeak> precursor_isotope_peaks(
    const std::string& peptide, double abundance_pct, std::size_t top_n) {
    if (top_n == 0) {
        throw std::runtime_error(
            "At least one theoretical precursor isotope is required");
    }
    auto& model = enriched_isotopologue(abundance_pct);
    IsotopeDistribution distribution;
    const auto decorated = isotope_peptide(peptide);
    if (!model.computePeptideIsotopicDistribution(
            decorated, distribution)) {
        throw std::runtime_error(
            "Could not compute theoretical precursor isotope envelope for " +
            peptide);
    }
    const std::size_t count = std::min(
        distribution.vMass.size(), distribution.vProb.size());
    std::vector<TheoreticalIsotopePeak> peaks;
    peaks.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (distribution.vProb[index] > 0.0 &&
            std::isfinite(distribution.vProb[index]) &&
            std::isfinite(distribution.vMass[index])) {
            peaks.push_back({
                distribution.vMass[index], distribution.vProb[index]});
        }
    }
    std::stable_sort(
        peaks.begin(), peaks.end(), [](const auto& left, const auto& right) {
            if (left.probability != right.probability) {
                return left.probability > right.probability;
            }
            return left.neutral_mass < right.neutral_mass;
        });
    if (peaks.size() > top_n) peaks.resize(top_n);
    if (peaks.empty()) {
        throw std::runtime_error(
            "Theoretical precursor isotope envelope is empty for " + peptide);
    }
    return peaks;
}

ProductIsotopeEnvelopes product_isotope_envelopes(
    const std::string& peptide, double abundance_pct) {
    auto& model = enriched_isotopologue(abundance_pct);
    ProductIsotopeEnvelopes result;
    const auto decorated = isotope_peptide(peptide);
    if (!model.computeProductIon(
            decorated, result.y_mass, result.y_probability,
            result.b_mass, result.b_probability)) {
        throw std::runtime_error(
            "Could not compute theoretical product-ion isotope envelopes for " +
            peptide);
    }
    return result;
}

} // namespace aerith
