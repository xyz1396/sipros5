#pragma once

#include "pipeline.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace aerith {

using IonIntensityMap = std::unordered_map<std::string, double>;

double psm_intensity(const Psm& psm);
std::string ion_form(const Psm& psm);
void update_ion_intensity(
    IonIntensityMap& intensities, const std::string& ion, double intensity);
double top_three_intensity(const IonIntensityMap& intensities);
double summed_intensity(const IonIntensityMap& intensities);
double summed_isotope_apex_intensity(
    const std::array<double, 3>& isotope_apices);
double summed_isotope_apex_intensity(
    const std::vector<double>& isotope_apices);
std::vector<double> mbr_posterior_probabilities(
    const std::vector<double>& target_scores,
    const std::vector<double>& decoy_scores,
    const std::vector<double>& identified_target_scores,
    const std::vector<double>& identified_decoy_scores,
    double* false_prior = nullptr);

struct IntensityNormalizer {
    static constexpr std::size_t bins = 10;

    std::unordered_map<std::string, std::size_t> ion_bins;
    std::vector<std::array<double, bins>> factors;
    std::array<double, bins - 1> upper_mz{};

    double factor(std::size_t sample, const std::string& ion) const;
};

IntensityNormalizer build_intensity_normalizer(
    const std::vector<IonIntensityMap>& sample_ions,
    std::vector<std::pair<double, std::string>> ion_mz,
    bool enabled);
double normalized_intensity(
    const IonIntensityMap& ions, std::size_t sample,
    const IntensityNormalizer& normalizer);
double normalized_top_three(
    const IonIntensityMap& ions, std::size_t sample,
    const IntensityNormalizer& normalizer);
std::vector<double> maxlfq(
    const std::vector<IonIntensityMap>& sample_ions,
    const IntensityNormalizer& normalizer);

} // namespace aerith
