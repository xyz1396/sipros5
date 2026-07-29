#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace aerith {

struct Config;

struct TheoreticalIsotopePeak {
    double neutral_mass = 0.0;
    double probability = 0.0;
};

struct ProductIsotopeEnvelopes {
    std::vector<std::vector<double>> y_mass;
    std::vector<std::vector<double>> y_probability;
    std::vector<std::vector<double>> b_mass;
    std::vector<std::vector<double>> b_probability;
};

void initialize_sip_isotope_model(const Config& config);
bool sip_isotope_model_enabled();
std::vector<TheoreticalIsotopePeak> precursor_isotope_peaks(
    const std::string& peptide, double abundance_pct, std::size_t top_n);
ProductIsotopeEnvelopes product_isotope_envelopes(
    const std::string& peptide, double abundance_pct);

} // namespace aerith
