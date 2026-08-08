#pragma once

#include "filter.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace aerith {

struct PredictionCacheFragment {
    float mz = 0.0f;
    float intensity = 0.0f;
    char ion_kind = '\0';
    std::uint32_t ion_position = 0;
    std::int32_t charge = 1;
};

struct PredictionCacheEntry {
    bool has_spectrum = false;
    bool has_rt = false;
    float rt = 0.0f;
    std::vector<PredictionCacheFragment> fragments;
};

struct PredictionCacheData {
    bool compatible = false;
    std::unordered_map<std::string, PredictionCacheEntry> entries;
};

std::filesystem::path prediction_cache_file_path(const Config& config);
std::string prediction_cache_peptide_body(const std::string& peptide);
std::string prediction_cache_key(const std::string& peptide, int charge);
PredictionCacheData read_prediction_cache(
    const Config& config, bool load_fragments);
void update_prediction_cache(
    const Config& config,
    const std::unordered_map<std::string, PredictionCacheEntry>& updates);

} // namespace aerith
