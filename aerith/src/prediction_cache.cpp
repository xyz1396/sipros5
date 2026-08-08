#include "prediction_cache.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <string_view>

namespace aerith {
namespace {

constexpr std::uint64_t kPredictionCacheMagic = 0x4145525052454433ULL;
constexpr std::uint8_t kHasSpectrum = 1;
constexpr std::uint8_t kHasRt = 2;
constexpr std::uint64_t kMaximumEntries = 1000000000ULL;
constexpr std::uint32_t kMaximumKeySize = 65536;
constexpr std::uint32_t kMaximumFragments = 10000;

template <typename Value>
bool read_value(std::istream& input, Value& value) {
    return static_cast<bool>(input.read(
        reinterpret_cast<char*>(&value), sizeof(value)));
}

template <typename Value>
void write_value(std::ostream& output, const Value& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void hash_bytes(std::uint64_t& hash, const void* bytes, std::size_t size) {
    const auto* value = static_cast<const unsigned char*>(bytes);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= value[index];
        hash *= 1099511628211ULL;
    }
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_bytes(hash, value.data(), value.size());
    constexpr unsigned char separator = 0xff;
    hash_bytes(hash, &separator, 1);
}

std::uint64_t prediction_cache_fingerprint(const Config& config) {
    std::uint64_t hash = 1469598103934665603ULL;
    hash_string(hash, config.spectrum_model_path);
    hash_string(hash, config.rt_model_path);
    return hash;
}

bool read_entry(
    std::istream& input, bool load_fragments,
    std::string& key, PredictionCacheEntry& entry) {
    std::uint32_t key_size = 0;
    std::uint8_t flags = 0;
    std::uint32_t fragment_count = 0;
    if (!read_value(input, key_size) || key_size > kMaximumKeySize) {
        return false;
    }
    key.assign(key_size, '\0');
    if (!input.read(key.data(), key.size()) ||
        !read_value(input, flags) ||
        !read_value(input, entry.rt) ||
        !read_value(input, fragment_count) ||
        fragment_count > kMaximumFragments) {
        return false;
    }
    entry.has_spectrum = (flags & kHasSpectrum) != 0;
    entry.has_rt = (flags & kHasRt) != 0;
    if (!entry.has_spectrum && fragment_count != 0) return false;
    if (!load_fragments) {
        constexpr std::streamoff fragment_bytes =
            sizeof(float) * 2 + sizeof(char) +
            sizeof(std::uint32_t) + sizeof(std::int32_t);
        input.seekg(
            static_cast<std::streamoff>(fragment_count) * fragment_bytes,
            std::ios::cur);
        return static_cast<bool>(input);
    }
    entry.fragments.resize(fragment_count);
    for (auto& fragment : entry.fragments) {
        if (!read_value(input, fragment.mz) ||
            !read_value(input, fragment.intensity) ||
            !read_value(input, fragment.ion_kind) ||
            !read_value(input, fragment.ion_position) ||
            !read_value(input, fragment.charge)) {
            return false;
        }
    }
    return true;
}

void write_entry(
    std::ostream& output, const std::string& key,
    const PredictionCacheEntry& entry) {
    write_value(output, static_cast<std::uint32_t>(key.size()));
    output.write(key.data(), key.size());
    const auto flags = static_cast<std::uint8_t>(
        (entry.has_spectrum ? kHasSpectrum : 0) |
        (entry.has_rt ? kHasRt : 0));
    write_value(output, flags);
    write_value(output, entry.rt);
    write_value(
        output, static_cast<std::uint32_t>(entry.fragments.size()));
    for (const auto& fragment : entry.fragments) {
        write_value(output, fragment.mz);
        write_value(output, fragment.intensity);
        write_value(output, fragment.ion_kind);
        write_value(output, fragment.ion_position);
        write_value(output, fragment.charge);
    }
}

} // namespace

std::filesystem::path prediction_cache_file_path(const Config& config) {
    return config.prediction_cache_path.empty()
        ? std::filesystem::path{}
        : std::filesystem::path(config.prediction_cache_path + ".bin");
}

std::string prediction_cache_peptide_body(const std::string& peptide) {
    const std::string_view view(peptide);
    const auto first_bracket = view.find('[');
    const auto last_bracket = view.rfind(']');
    if (first_bracket != std::string_view::npos &&
        last_bracket > first_bracket &&
        view.find('[', first_bracket + 1) == std::string_view::npos &&
        first_bracket <= 2 && view.size() - last_bracket <= 3) {
        return std::string(view.substr(
            first_bracket + 1, last_bracket - first_bracket - 1));
    }
    const auto first_dot = view.find('.');
    const auto last_dot = view.rfind('.');
    if (first_dot != std::string_view::npos && last_dot > first_dot) {
        return std::string(view.substr(
            first_dot + 1, last_dot - first_dot - 1));
    }
    return peptide;
}

std::string prediction_cache_key(const std::string& peptide, int charge) {
    return prediction_cache_peptide_body(peptide) + '\x1f' +
        std::to_string(charge);
}

PredictionCacheData read_prediction_cache(
    const Config& config, bool load_fragments) {
    PredictionCacheData result;
    const auto path = prediction_cache_file_path(config);
    if (path.empty() || !std::filesystem::is_regular_file(path)) {
        return result;
    }
    std::ifstream input(path, std::ios::binary);
    std::uint64_t magic = 0;
    std::uint64_t fingerprint = 0;
    std::uint64_t count = 0;
    if (!read_value(input, magic) ||
        !read_value(input, fingerprint) ||
        !read_value(input, count) ||
        magic != kPredictionCacheMagic ||
        fingerprint != prediction_cache_fingerprint(config) ||
        count > kMaximumEntries) {
        return result;
    }
    result.entries.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t index = 0; index < count; ++index) {
        std::string key;
        PredictionCacheEntry entry;
        if (!read_entry(input, load_fragments, key, entry)) return {};
        result.entries.insert_or_assign(std::move(key), std::move(entry));
    }
    if (!input) return {};
    result.compatible = true;
    return result;
}

void update_prediction_cache(
    const Config& config,
    const std::unordered_map<std::string, PredictionCacheEntry>& updates) {
    const auto path = prediction_cache_file_path(config);
    if (path.empty() || updates.empty()) return;
    auto cache = read_prediction_cache(config, true);
    if (!cache.compatible) cache.entries.clear();
    for (const auto& [key, update] : updates) {
        auto& entry = cache.entries[key];
        if (update.has_spectrum) {
            entry.has_spectrum = true;
            entry.fragments = update.fragments;
        }
        if (update.has_rt) {
            entry.has_rt = true;
            entry.rt = update.rt;
        }
    }
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const auto temporary = path.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) return;
    write_value(output, kPredictionCacheMagic);
    write_value(output, prediction_cache_fingerprint(config));
    write_value(
        output, static_cast<std::uint64_t>(cache.entries.size()));
    std::vector<std::string> keys;
    keys.reserve(cache.entries.size());
    for (const auto& [key, entry] : cache.entries) keys.push_back(key);
    std::sort(keys.begin(), keys.end());
    for (const auto& key : keys) {
        write_entry(output, key, cache.entries.at(key));
    }
    output.close();
    if (!output) {
        std::error_code error;
        std::filesystem::remove(temporary, error);
        return;
    }
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) std::filesystem::remove(temporary, error);
}

} // namespace aerith
