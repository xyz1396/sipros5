#include "pipeline.hpp"
#include "isotope.hpp"
#include "prediction_cache.hpp"
#include "torch_device.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <omp.h>

#ifdef AERITH_WITH_TORCH
#include <H5Cpp.h>
#include <torch/script.h>
#endif

namespace aerith {

#ifndef AERITH_WITH_TORCH

struct SpectrumPredictionLibrary::Impl {};

SpectrumPredictionLibrary::SpectrumPredictionLibrary()
    : impl_(std::make_unique<Impl>()) {}
SpectrumPredictionLibrary::~SpectrumPredictionLibrary() = default;
SpectrumPredictionLibrary::SpectrumPredictionLibrary(
    SpectrumPredictionLibrary&&) noexcept = default;
SpectrumPredictionLibrary& SpectrumPredictionLibrary::operator=(
    SpectrumPredictionLibrary&&) noexcept = default;
std::string SpectrumPredictionLibrary::device() const { return {}; }
StageTiming SpectrumPredictionLibrary::timing() const { return {}; }
StageTiming SpectrumPredictionLibrary::cache_read_timing() const { return {}; }
StageTiming SpectrumPredictionLibrary::cache_write_timing() const { return {}; }

SpectrumPredictionLibrary SpectralEntropyFeature::predict(
    const Config& config, const Dataset&) {
    if (!config.spectrum_paths.empty()) {
        throw std::runtime_error(
            "Aerith was built without LibTorch; rebuild with Torch_DIR to use --spectra");
    }
    return {};
}

void SpectralEntropyFeature::add(
    const Config& config, Dataset& data,
    const SpectrumPredictionLibrary&) {
    add(config, data);
}

void SpectralEntropyFeature::add(const Config& config, Dataset&) {
    if (!config.spectrum_paths.empty()) {
        throw std::runtime_error(
            "Aerith was built without LibTorch; rebuild with Torch_DIR to use --spectra");
    }
}

#else

constexpr double kProton = 1.007276466621;
constexpr double kWater = 18.0105646837;
constexpr double kMinFragmentMz = 200.0;
constexpr double kMaxFragmentMz = 1800.0;
constexpr double kMinRawRelativeIntensity = 0.001;
constexpr double kPredictedBasePeakCutoff = 0.01;
constexpr std::size_t kMaximumDiannFragments = 100;
constexpr std::size_t kEntropyFragments = 20;
constexpr std::size_t kCpuPredictionBatch = 4096;
constexpr std::size_t kCudaPredictionBatch = 1024;

const std::unordered_map<char, double> kAminoAcidMass{
    {'A', 71.037113805}, {'R', 156.10111105}, {'N', 114.04292747},
    {'D', 115.026943065}, {'C', 103.009184505}, {'E', 129.042593135},
    {'Q', 128.05857754}, {'G', 57.021463735}, {'H', 137.058911875},
    {'I', 113.084063975}, {'L', 113.084063975}, {'K', 128.094963015},
    {'M', 131.040484645}, {'F', 147.068413945}, {'P', 97.052763875},
    {'S', 87.032028435}, {'T', 101.047678505}, {'W', 186.07931298},
    {'Y', 163.063328575}, {'V', 99.068413945}};

const std::unordered_map<char, double> kModificationShift{
    {'~', 15.994915}, {'!', 0.984016}, {'@', 79.966332},
    {'>', 79.966332}, {'<', 79.966332}, {'%', 42.010565},
    {'^', 14.015650}, {'&', 28.031300}, {'*', 42.046950},
    {'(', 28.990164}, {')', 44.985079}, {'/', 57.021464},
    {'$', 45.987721}};

struct ParsedPeptide {
    std::string key;
    std::vector<std::int64_t> tokens;
    std::vector<double> residue_masses;
    double nterm_shift = 0.0;
    int charge = 0;
};

struct Fragment {
    float mz = 0.0f;
    float intensity = 0.0f;
    char ion_kind = '\0';
    std::size_t ion_position = 0;
    int charge = 1;
};

struct SpectrumPredictionLibrary::Impl {
    std::unordered_map<std::string, std::vector<Fragment>> predictions;
    std::string device;
    StageTiming timing;
    StageTiming cache_read_timing;
    StageTiming cache_write_timing;
};

struct SpectrumCacheLoad {
    bool compatible = false;
    std::size_t hits = 0;
};

SpectrumCacheLoad load_spectrum_library_cache(
    const Config& config, const Dataset& unique_peptides,
    SpectrumPredictionLibrary::Impl& library) {
    const auto timing_begin = std::chrono::steady_clock::now();
    const auto timing_cpu_begin = std::clock();
    SpectrumCacheLoad result;
    auto cache = read_prediction_cache(config, true);
    result.compatible = cache.compatible;
    if (cache.compatible) {
        for (const auto& psm : unique_peptides.rows) {
            const auto key = prediction_cache_key(psm.peptide, psm.charge);
            const auto found = cache.entries.find(key);
            if (found == cache.entries.end() || !found->second.has_spectrum ||
                library.predictions.find(key) != library.predictions.end()) {
                continue;
            }
            std::vector<Fragment> fragments;
            fragments.reserve(found->second.fragments.size());
            for (const auto& cached : found->second.fragments) {
                fragments.push_back({
                    cached.mz, cached.intensity, cached.ion_kind,
                    cached.ion_position, cached.charge});
            }
            library.predictions.emplace(key, std::move(fragments));
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

StageTiming append_spectrum_library_cache(
    const Config& config,
    const std::unordered_map<std::string, std::vector<Fragment>>& predictions,
    const std::unordered_set<std::string>* allowed) {
    const auto timing_begin = std::chrono::steady_clock::now();
    const auto timing_cpu_begin = std::clock();
    std::unordered_map<std::string, PredictionCacheEntry> updates;
    updates.reserve(predictions.size());
    for (const auto& [key, fragments] : predictions) {
        if (allowed != nullptr && allowed->find(key) == allowed->end()) {
            continue;
        }
        auto& update = updates[key];
        update.has_spectrum = true;
        update.fragments.reserve(fragments.size());
        for (const auto& fragment : fragments) {
            update.fragments.push_back({
                fragment.mz, fragment.intensity, fragment.ion_kind,
                static_cast<std::uint32_t>(fragment.ion_position),
                static_cast<std::int32_t>(fragment.charge)});
        }
    }
    update_prediction_cache(config, updates);
    return {
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - timing_begin).count(),
        static_cast<double>(std::clock() - timing_cpu_begin) /
            CLOCKS_PER_SEC};
}

SpectrumPredictionLibrary::SpectrumPredictionLibrary()
    : impl_(std::make_unique<Impl>()) {}
SpectrumPredictionLibrary::~SpectrumPredictionLibrary() = default;
SpectrumPredictionLibrary::SpectrumPredictionLibrary(
    SpectrumPredictionLibrary&&) noexcept = default;
SpectrumPredictionLibrary& SpectrumPredictionLibrary::operator=(
    SpectrumPredictionLibrary&&) noexcept = default;
std::string SpectrumPredictionLibrary::device() const {
    return impl_->device;
}
StageTiming SpectrumPredictionLibrary::timing() const {
    return impl_->timing;
}
StageTiming SpectrumPredictionLibrary::cache_read_timing() const {
    return impl_->cache_read_timing;
}
StageTiming SpectrumPredictionLibrary::cache_write_timing() const {
    return impl_->cache_write_timing;
}

struct Spectrum {
    std::vector<float> mz;
    std::vector<float> intensity;
    double lower_mz = -std::numeric_limits<double>::infinity();
    double upper_mz = std::numeric_limits<double>::infinity();
};

template <typename T>
std::vector<T> read_1d(H5::H5File& file, const std::string& name,
                       const H5::PredType& type) {
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("HDF5 dataset is not one-dimensional: " + name);
    }
    hsize_t size = 0;
    space.getSimpleExtentDims(&size, nullptr);
    std::vector<T> values(static_cast<std::size_t>(size));
    if (size != 0) dataset.read(values.data(), type);
    return values;
}

std::vector<std::string> read_fixed_strings(H5::H5File& file,
                                             const std::string& name) {
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("HDF5 dataset is not one-dimensional: " + name);
    }
    const H5::StrType type = dataset.getStrType();
    if (type.isVariableStr()) {
        throw std::runtime_error("Expected a fixed-width HDF5 string table: " + name);
    }
    hsize_t count = 0;
    space.getSimpleExtentDims(&count, nullptr);
    const std::size_t width = type.getSize();
    std::vector<char> buffer(static_cast<std::size_t>(count) * width, '\0');
    if (count != 0) dataset.read(buffer.data(), type);
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(count));
    for (hsize_t i = 0; i < count; ++i) {
        const char* value = buffer.data() + static_cast<std::size_t>(i) * width;
        result.emplace_back(value, std::find(value, value + width, '\0'));
    }
    return result;
}

std::pair<double, double> scan_window(const std::string& scan_filter) {
    // Thermo filters written by Raxport end in e.g. "[87.0000-814.0000]".
    const auto open = scan_filter.rfind('[');
    const auto dash = open == std::string::npos
        ? std::string::npos : scan_filter.find('-', open + 1);
    const auto close = dash == std::string::npos
        ? std::string::npos : scan_filter.find(']', dash + 1);
    if (open == std::string::npos || dash == std::string::npos ||
        close == std::string::npos) {
        return {-std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    }
    try {
        return {std::stod(scan_filter.substr(open + 1, dash - open - 1)),
                std::stod(scan_filter.substr(dash + 1, close - dash - 1))};
    } catch (const std::exception&) {
        return {-std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    }
}

template <typename T>
std::vector<T> read_slice(H5::DataSet& dataset, const H5::PredType& type,
                          hsize_t start, hsize_t count) {
    std::vector<T> values(static_cast<std::size_t>(count));
    if (count == 0) return values;
    H5::DataSpace file_space = dataset.getSpace();
    hsize_t size = 0;
    file_space.getSimpleExtentDims(&size, nullptr);
    if (start > size || count > size - start) {
        throw std::runtime_error("HDF5 peak slice is out of bounds");
    }
    H5::DataSpace memory_space(1, &count);
    file_space.selectHyperslab(H5S_SELECT_SET, &count, &start);
    dataset.read(values.data(), type, memory_space, file_space);
    return values;
}

const std::unordered_map<std::string, std::int64_t>& diann_dictionary() {
    // Exact token table distributed as dict.txt with DIA-NN 2.6.1. Keeping it
    // in the executable prevents the model vocabulary from drifting at runtime.
    static const std::unordered_map<std::string, std::int64_t> dictionary{
        {"0", 0}, {"<", 1}, {">", 2}, {"G", 3}, {"A", 4}, {"V", 5},
        {"I", 6}, {"L", 7}, {"P", 8}, {"F", 9}, {"W", 10}, {"M", 11},
        {"X", 11}, {"S", 13}, {"T", 14}, {"Y", 15}, {"Q", 16}, {"E", 17},
        {"N", 18}, {"D", 19}, {"K", 20}, {"O", 20}, {"R", 21}, {"H", 22},
        {"C", 24}, {"U", 24},
        {"C(UniMod:4)", 25}, {"C(unimod:4)", 25}, {"C(160)", 25},
        {"C(Carbamidomethyl (C))", 25}, {"C(Carbamidomethyl)", 25},
        {"M(UniMod:35)", 26}, {"M(unimod:35)", 26},
        {"M(Oxidation (M))", 26}, {"M(Oxidation)", 26}, {"M(147)", 26},
        {"K(GlyGly (K))", 27}, {"K(UniMod:121)", 27},
        {"K(mTRAQ-Lys0)", 28}, {"K(mTRAQ)", 28}, {"K(UniMod:888)", 28},
        {"<(UniMod:1)", 29}, {"<(unimod:1)", 29},
        {"<(mTRAQ-Nter0)", 30}, {"<(mTRAQ)", 30}, {"<(UniMod:888)", 30},
        {"S(Phospho (STY))", 31}, {"S(UniMod:21)", 31}, {"S(ph)", 31},
        {"T(Phospho (STY))", 32}, {"T(UniMod:21)", 32}, {"T(ph)", 32},
        {"Y(Phospho (STY))", 33}, {"Y(UniMod:21)", 33}, {"Y(ph)", 33},
        {"N(UniMod:7)", 19}, {"N[Deamidation (NQ)]", 19}, {"N(de)", 19},
        {"Q(UniMod:7)", 17}, {"Q[Deamidation (NQ)]", 17}, {"Q(de)", 17},
        {"<(TMT)", 36}, {"K(TMT)", 37},
        {"<(UniMod:255)", 38}, {"<(Dimethyl)", 38},
        {"K(UniMod:255)", 39}, {"K(Dimethyl)", 39},
    };
    return dictionary;
}

std::string_view spectrum_peptide_body_view(std::string_view peptide) {
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

std::string spectrum_peptide_body(const std::string& peptide) {
    return prediction_cache_peptide_body(peptide);
}

std::string token_for(char residue, char modification, bool fixed_cam) {
    if (residue == 'C' && (fixed_cam || modification == '/')) return "C(UniMod:4)";
    if (residue == 'M' && modification == '~') return "M(UniMod:35)";
    if ((residue == 'S' || residue == 'T' || residue == 'Y') &&
        (modification == '@' || modification == '>' || modification == '<')) {
        return std::string(1, residue) + "(UniMod:21)";
    }
    // Match MSBooster's DIA-NN preparation: deamidation is stripped from the
    // neural-network input, then restored below in theoretical fragment m/z.
    if (residue == 'K' && modification == '&') return "K(UniMod:255)";
    return std::string(1, residue);
}

ParsedPeptide parse_peptide(
    const std::string& peptide, int charge,
    const std::unordered_map<std::string, std::int64_t>& dictionary) {
    ParsedPeptide parsed;
    parsed.charge = charge;
    const std::string body = spectrum_peptide_body(peptide);
    parsed.key = body + '\x1f' + std::to_string(charge);

    auto required_token = [&](const std::string& token) {
        const auto found = dictionary.find(token);
        if (found == dictionary.end()) {
            throw std::runtime_error("DIA-NN dictionary is missing token: " + token);
        }
        return found->second;
    };

    std::size_t position = 0;
    std::string nterm = "<";
    if (!body.empty() && body.front() == '%') {
        nterm = "<(UniMod:1)";
        parsed.nterm_shift = kModificationShift.at('%');
        ++position;
    }
    parsed.tokens.push_back(required_token(nterm));

    while (position < body.size()) {
        const char residue = body[position];
        const auto mass = kAminoAcidMass.find(residue);
        if (mass == kAminoAcidMass.end()) {
            // Ignore terminal/search punctuation, but reject an unknown letter.
            if (std::isalpha(static_cast<unsigned char>(residue))) {
                throw std::runtime_error("Unsupported residue in PIN peptide: " + peptide);
            }
            ++position;
            continue;
        }
        ++position;
        char modification = '\0';
        if (position < body.size() && kModificationShift.count(body[position]) != 0) {
            modification = body[position++];
        }

        // Sipros default chemistry has fixed CAM. S-nitrosylation replaces CAM;
        // all other explicit C modifications retain it unless '/' is the CAM token.
        const bool fixed_cam = residue == 'C' && modification != '(';
        double residue_mass = mass->second + (fixed_cam ? 57.021463735 : 0.0);
        if (modification != '\0') {
            if (modification == '(') residue_mass = mass->second + kModificationShift.at('(');
            else if (!(residue == 'C' && modification == '/' && fixed_cam)) {
                residue_mass += kModificationShift.at(modification);
            }
        }
        const std::string token = token_for(residue, modification, fixed_cam);
        const auto encoded = dictionary.find(token);
        parsed.tokens.push_back(encoded == dictionary.end()
                                    ? required_token(std::string(1, residue))
                                    : encoded->second);
        parsed.residue_masses.push_back(residue_mass);
    }
    parsed.tokens.push_back(required_token(">"));
    if (parsed.residue_masses.size() < 4) {
        throw std::runtime_error("DIA-NN fragmentation requires at least four residues: " + peptide);
    }
    return parsed;
}

void merge_close_fragments(std::vector<Fragment>& fragments) {
    if (fragments.size() < 2) return;
    std::sort(fragments.begin(), fragments.end(), [](const Fragment& a, const Fragment& b) {
        return a.mz < b.mz;
    });
    std::vector<Fragment> merged;
    merged.reserve(fragments.size());
    for (const auto& fragment : fragments) {
        if (!merged.empty() && fragment.mz - merged.back().mz < 0.00001f) {
            merged.back().mz = (merged.back().mz + fragment.mz) / 2.0f;
            merged.back().intensity += fragment.intensity;
        } else {
            merged.push_back(fragment);
        }
    }
    fragments.swap(merged);
}

void preprocess_for_entropy(
    std::vector<Fragment>& fragments, bool select_top_fragments = true) {
    merge_close_fragments(fragments);
    if (fragments.empty()) return;
    float maximum = 0.0f;
    for (const auto& fragment : fragments) maximum = std::max(maximum, fragment.intensity);
    if (select_top_fragments) {
        const float cutoff =
            static_cast<float>(kPredictedBasePeakCutoff) * maximum;
        fragments.erase(
            std::remove_if(
                fragments.begin(), fragments.end(),
                [&](const Fragment& fragment) {
                    return fragment.intensity < cutoff;
                }),
            fragments.end());
        if (fragments.size() > kEntropyFragments) {
            std::nth_element(
                fragments.begin(),
                fragments.begin() + kEntropyFragments,
                fragments.end(),
                [](const Fragment& left, const Fragment& right) {
                    return left.intensity > right.intensity;
                });
            fragments.resize(kEntropyFragments);
        }
    }
    maximum = 0.0f;
    for (const auto& fragment : fragments) maximum = std::max(maximum, fragment.intensity);
    if (maximum > 0.0f) {
        for (auto& fragment : fragments) fragment.intensity /= maximum;
    }
    std::sort(fragments.begin(), fragments.end(), [](const Fragment& a, const Fragment& b) {
        return a.mz < b.mz;
    });
}

std::vector<Fragment> predicted_fragments(
    const ParsedPeptide& peptide, const float* values,
    std::size_t output_columns, bool sip_prediction) {
    const std::size_t residues = peptide.residue_masses.size();
    const std::size_t cleavages = residues - 3;
    if (output_columns != 4 * cleavages) {
        throw std::runtime_error(
            "Unexpected DIA-NN fragmentation output shape");
    }
    float raw_maximum = 0.0f;
    for (std::size_t i = 0; i < output_columns; ++i) {
        raw_maximum = std::max(raw_maximum, values[i]);
    }
    std::vector<double> prefix(residues);
    double running = peptide.nterm_shift;
    for (std::size_t i = 0; i < residues; ++i) {
        running += peptide.residue_masses[i];
        prefix[i] = running;
    }
    std::vector<double> suffix(residues + 1, 0.0);
    for (std::size_t i = 1; i <= residues; ++i) {
        suffix[i] =
            suffix[i - 1] + peptide.residue_masses[residues - i];
    }

    struct Candidate {
        float mz;
        float raw;
        char ion_kind;
        std::size_t ion_position;
        int charge;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(output_columns);
    for (std::size_t channel = 0; channel < 4; ++channel) {
        const bool y_ion = channel < 2;
        const int fragment_charge = channel % 2 == 0 ? 1 : 2;
        for (std::size_t position = 0; position < cleavages; ++position) {
            const float raw = values[channel * cleavages + position];
            if (!(raw_maximum > 0.0f) ||
                raw / raw_maximum < kMinRawRelativeIntensity) {
                continue;
            }
            const std::size_t number =
                y_ion ? residues - 1 - position : position + 3;
            const double neutral =
                y_ion ? suffix[number] + kWater : prefix[number - 1];
            const double mz =
                (neutral + fragment_charge * kProton) / fragment_charge;
            if (mz >= kMinFragmentMz && mz <= kMaxFragmentMz) {
                candidates.push_back({
                    static_cast<float>(mz), raw,
                    y_ion ? 'y' : 'b', number, fragment_charge});
            }
        }
    }
    std::sort(
        candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b) {
            return a.raw > b.raw;
        });
    if (candidates.size() > kMaximumDiannFragments) {
        candidates.resize(kMaximumDiannFragments);
    }
    float retained_maximum = 0.0f;
    for (const auto& candidate : candidates) {
        retained_maximum = std::max(retained_maximum, candidate.raw);
    }
    std::vector<Fragment> fragments;
    fragments.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        const auto packed = static_cast<float>(std::lround(
            candidate.raw / retained_maximum * 60000.0f));
        fragments.push_back({
            candidate.mz, packed, candidate.ion_kind,
            candidate.ion_position, candidate.charge});
    }
    preprocess_for_entropy(fragments, !sip_prediction);
    return fragments;
}

std::vector<Fragment> shift_predicted_fragments(
    const Psm& psm, const std::vector<Fragment>& predicted,
    std::size_t top_isotopes) {
    const auto envelopes = product_isotope_envelopes(
        psm.peptide, psm.ms2_isotopic_abundance);
    std::vector<Fragment> shifted;
    shifted.reserve(predicted.size() * 8);
    for (const auto& fragment : predicted) {
        const auto& masses = fragment.ion_kind == 'y'
            ? envelopes.y_mass : envelopes.b_mass;
        const auto& probabilities = fragment.ion_kind == 'y'
            ? envelopes.y_probability : envelopes.b_probability;
        if (fragment.ion_position == 0 ||
            fragment.ion_position > masses.size() ||
            fragment.ion_position > probabilities.size()) {
            continue;
        }
        const auto& envelope_mass = masses[fragment.ion_position - 1];
        const auto& envelope_probability =
            probabilities[fragment.ion_position - 1];
        const std::size_t count = std::min(
            envelope_mass.size(), envelope_probability.size());
        std::vector<std::size_t> isotope_order(count);
        std::iota(isotope_order.begin(), isotope_order.end(), 0);
        std::stable_sort(
            isotope_order.begin(), isotope_order.end(),
            [&](std::size_t left, std::size_t right) {
                return envelope_probability[left] >
                       envelope_probability[right];
            });
        if (isotope_order.size() > top_isotopes) {
            isotope_order.resize(top_isotopes);
        }
        // DIA-NN predicts one intensity per product ion. Renormalize the
        // selected theoretical probabilities so isotope expansion preserves
        // that product ion's total predicted intensity.
        double selected_probability = 0.0;
        for (const auto isotope : isotope_order) {
            selected_probability += envelope_probability[isotope];
        }
        if (!(selected_probability > 0.0)) continue;
        for (const auto isotope : isotope_order) {
            const double relative =
                envelope_probability[isotope] / selected_probability;
            const double mz =
                (envelope_mass[isotope] +
                 static_cast<double>(fragment.charge) * kProton) /
                static_cast<double>(fragment.charge);
            if (mz < kMinFragmentMz || mz > kMaxFragmentMz) continue;
            shifted.push_back({
                static_cast<float>(mz),
                static_cast<float>(fragment.intensity * relative),
                fragment.ion_kind, fragment.ion_position,
                fragment.charge});
        }
    }
    // SIP selection is per product ion. Do not apply the conventional
    // spectrum-wide 1% cutoff or 20-fragment limit after isotope expansion.
    preprocess_for_entropy(shifted, false);
    return shifted;
}

std::unordered_map<std::string, std::vector<Fragment>> predict_fragments(
    const std::filesystem::path& model_path, const Dataset& data,
    std::string& selected_device) {
    const bool sip_prediction = sip_isotope_model_enabled();
    const auto& dictionary = diann_dictionary();
    std::unordered_map<std::string, ParsedPeptide> unique;
    unique.reserve(data.rows.size() / 2);
    for (const auto& row : data.rows) {
        const std::string key = prediction_cache_key(
            row.peptide, row.charge);
        if (unique.find(key) != unique.end()) continue;
        auto parsed = parse_peptide(row.peptide, row.charge, dictionary);
        unique.emplace(std::move(key), std::move(parsed));
    }

    // unordered_map iteration depends on its bucket count, which used to make
    // Torch batch composition depend on how many duplicate PSM rows happened
    // to accompany the same unique peptide set. Canonical key order makes
    // unique-catalog and all-row prediction bit-reproducible.
    std::vector<ParsedPeptide*> ordered;
    ordered.reserve(unique.size());
    for (auto& entry : unique) ordered.push_back(&entry.second);
    std::sort(ordered.begin(), ordered.end(), [](const auto* left,
                                                 const auto* right) {
        return left->key < right->key;
    });
    std::map<std::size_t, std::vector<ParsedPeptide*>> groups;
    for (auto* peptide : ordered) {
        groups[peptide->tokens.size()].push_back(peptide);
    }

    // This checkpoint is DIA-NN 2.6.1 Academia's models/fr.d0.pt, renamed to
    // diann-2.6.1-fragmentation.pt. Source release:
    // https://github.com/vdemichev/DiaNN/releases/tag/2.0
    // It is loaded dynamically through LibTorch.
    return run_torch_prefer_cuda(
        "DIA-NN spectrum prediction", selected_device,
        [&](const torch::Device& device) {
        auto model = load_torch_model_on_device(model_path, device);
        c10::InferenceMode inference_mode;
        std::unordered_map<std::string, std::vector<Fragment>> predictions;
        predictions.reserve(unique.size());
        const std::size_t prediction_batch =
            device.is_cpu() ? kCpuPredictionBatch : kCudaPredictionBatch;

        for (auto& group_entry : groups) {
            auto& group = group_entry.second;
            const std::size_t token_count = group_entry.first;
            for (std::size_t begin = 0; begin < group.size();
                 begin += prediction_batch) {
                const std::size_t count =
                    std::min(prediction_batch, group.size() - begin);
                auto cpu_input = torch::empty(
                {static_cast<long>(count), static_cast<long>(token_count + 1)},
                torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU));
                auto accessor = cpu_input.accessor<std::int64_t, 2>();
                for (std::size_t row = 0; row < count; ++row) {
                    const auto& peptide = *group[begin + row];
                    accessor[row][0] = std::min(peptide.charge, 4) - 2;
                    for (std::size_t column = 0; column < token_count; ++column) {
                        accessor[row][column + 1] = peptide.tokens[column];
                    }
                }
                const auto input = move_torch_input(cpu_input, device);
                const auto output =
                    model.forward({input}).toTensor().to(torch::kCPU).contiguous();
                const std::size_t output_columns =
                    static_cast<std::size_t>(output.size(1));
                const float* values = output.data_ptr<float>();
                std::vector<std::vector<Fragment>> batch_predictions(count);
                std::exception_ptr failure;
                #pragma omp parallel for schedule(static) if(count >= 256)
                for (std::ptrdiff_t row = 0;
                     row < static_cast<std::ptrdiff_t>(count); ++row) {
                    try {
                        const auto local = static_cast<std::size_t>(row);
                        batch_predictions[local] = predicted_fragments(
                            *group[begin + local],
                            values + local * output_columns,
                            output_columns, sip_prediction);
                    } catch (...) {
                        #pragma omp critical(aerith_spectrum_prediction_failure)
                        if (!failure) failure = std::current_exception();
                    }
                }
                if (failure) std::rethrow_exception(failure);
                for (std::size_t row = 0; row < count; ++row) {
                    predictions.emplace(
                        group[begin + row]->key,
                        std::move(batch_predictions[row]));
                }
            }
        }
        return predictions;
    });
}

std::unordered_map<std::uint64_t, Spectrum> load_spectra(
    const std::string& path, const std::unordered_set<std::uint64_t>& requested) {
    H5::Exception::dontPrint();
    H5::H5File file(path, H5F_ACC_RDONLY);
    const auto scan_numbers = read_1d<int>(
        file, "/scans/scan_number", H5::PredType::NATIVE_INT);
    const auto ms_order = read_1d<int>(file, "/scans/ms_order", H5::PredType::NATIVE_INT);
    const auto peak_start = read_1d<long long>(
        file, "/scans/peak_start", H5::PredType::NATIVE_LLONG);
    const auto peak_count = read_1d<int>(
        file, "/scans/peak_count", H5::PredType::NATIVE_INT);
    const auto scan_filter_id = read_1d<int>(
        file, "/scans/scan_filter_id", H5::PredType::NATIVE_INT);
    const auto scan_filters = read_fixed_strings(file, "/string_tables/scan_filter");
    if (ms_order.size() != scan_numbers.size() || peak_start.size() != scan_numbers.size() ||
        peak_count.size() != scan_numbers.size() ||
        scan_filter_id.size() != scan_numbers.size()) {
        throw std::runtime_error("Raxport HDF5 scan datasets have inconsistent lengths");
    }
    H5::DataSet mz_dataset = file.openDataSet("/peaks/mz");
    H5::DataSet intensity_dataset = file.openDataSet("/peaks/intensity");
    std::unordered_map<std::uint64_t, Spectrum> spectra;
    spectra.reserve(requested.size());
    for (std::size_t i = 0; i < scan_numbers.size(); ++i) {
        const auto scan = static_cast<std::uint64_t>(scan_numbers[i]);
        if (ms_order[i] != 2 || requested.count(scan) == 0 ||
            peak_start[i] < 0 || peak_count[i] <= 0) continue;
        const auto mz_double = read_slice<double>(
            mz_dataset, H5::PredType::NATIVE_DOUBLE,
            static_cast<hsize_t>(peak_start[i]), static_cast<hsize_t>(peak_count[i]));
        const auto intensity_double = read_slice<double>(
            intensity_dataset, H5::PredType::NATIVE_DOUBLE,
            static_cast<hsize_t>(peak_start[i]), static_cast<hsize_t>(peak_count[i]));
        Spectrum spectrum;
        if (scan_filter_id[i] >= 0 &&
            static_cast<std::size_t>(scan_filter_id[i]) < scan_filters.size()) {
            const auto window = scan_window(scan_filters[scan_filter_id[i]]);
            spectrum.lower_mz = window.first;
            spectrum.upper_mz = window.second;
        }
        spectrum.mz.reserve(mz_double.size());
        spectrum.intensity.reserve(intensity_double.size());
        for (const auto value : mz_double) spectrum.mz.push_back(static_cast<float>(value));
        for (const auto value : intensity_double) {
            spectrum.intensity.push_back(static_cast<float>(value));
        }
        if (std::is_sorted(spectrum.mz.begin(), spectrum.mz.end())) {
            spectra.emplace(scan, std::move(spectrum));
            continue;
        }
        std::vector<std::size_t> order(spectrum.mz.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return spectrum.mz[a] < spectrum.mz[b];
        });
        Spectrum sorted;
        sorted.lower_mz = spectrum.lower_mz;
        sorted.upper_mz = spectrum.upper_mz;
        sorted.mz.reserve(order.size());
        sorted.intensity.reserve(order.size());
        for (const auto index : order) {
            sorted.mz.push_back(spectrum.mz[index]);
            sorted.intensity.push_back(spectrum.intensity[index]);
        }
        spectra.emplace(scan, std::move(sorted));
    }
    return spectra;
}

double entropy(const std::vector<float>& values) {
    double result = 0.0;
    for (const float value : values) {
        if (value != 0.0f) result += static_cast<double>(value) * std::log(value);
    }
    return -result;
}

std::vector<float> normalize_sum(const std::vector<float>& values) {
    double total = 0.0;
    for (const float value : values) total += value;
    std::vector<float> result(values.size(), 0.0f);
    const float denominator = static_cast<float>(total);
    if (denominator != 0.0f) {
        for (std::size_t i = 0; i < values.size(); ++i) {
            result[i] = values[i] / denominator;
        }
    }
    return result;
}

float entropy_similarity(const std::vector<Fragment>& predicted,
                         const Spectrum& experimental, double ppm) {
    std::vector<const Fragment*> in_window;
    in_window.reserve(predicted.size());
    for (const auto& fragment : predicted) {
        if (fragment.mz >= experimental.lower_mz &&
            fragment.mz <= experimental.upper_mz) {
            in_window.push_back(&fragment);
        }
    }
    if (in_window.size() < 2) return 0.0f;
    std::vector<float> predicted_intensity;
    std::vector<float> matched(in_window.size(), 0.0f);
    predicted_intensity.reserve(in_window.size());
    std::size_t start = 0;
    const double fraction = ppm * 1e-6;
    for (std::size_t i = 0; i < in_window.size(); ++i) {
        predicted_intensity.push_back(in_window[i]->intensity);
        const double minimum = in_window[i]->mz * (1.0 - fraction);
        const double maximum = in_window[i]->mz * (1.0 + fraction);
        while (start < experimental.mz.size() && experimental.mz[start] < minimum) ++start;
        for (std::size_t peak = start;
             peak < experimental.mz.size() && experimental.mz[peak] <= maximum; ++peak) {
            matched[i] = std::max(matched[i], experimental.intensity[peak]);
        }
    }
    if (std::count_if(matched.begin(), matched.end(), [](float value) {
            return value != 0.0f;
        }) < 2) return 0.0f;
    const auto normalized_predicted = normalize_sum(predicted_intensity);
    const auto normalized_matched = normalize_sum(matched);
    std::vector<float> mixture(in_window.size());
    for (std::size_t i = 0; i < mixture.size(); ++i) {
        mixture[i] = (normalized_predicted[i] + normalized_matched[i]) / 2.0f;
    }
    const double score = 1.0 -
        (2.0 * entropy(mixture) - entropy(normalized_matched) -
         entropy(normalized_predicted)) / std::log(4.0);
    // MSBooster writes this feature with four decimal places.
    return static_cast<float>(std::round(score * 10000.0) / 10000.0);
}

std::array<float, 4> sip_entropy_abundance_scores_for_testing() {
    Config config;
    config.sip_isotope = "C13";
    initialize_sip_isotope_model(config);
    Psm psm;
    psm.peptide = "[PEPTIDEK]";
    psm.charge = 2;
    psm.ms2_isotopic_abundance = 50.0;
    std::vector<Fragment> base;
    for (std::size_t position = 2; position <= 7; ++position) {
        base.push_back({
            static_cast<float>(250.0 + position * 50.0),
            static_cast<float>(1.0 / static_cast<double>(position)),
            'b', position, 1});
        base.push_back({
            static_cast<float>(275.0 + position * 50.0),
            static_cast<float>(1.0 / static_cast<double>(position + 1)),
            'y', position, 1});
    }
    const auto matching = shift_predicted_fragments(psm, base, 5);
    Spectrum experimental;
    experimental.lower_mz = kMinFragmentMz;
    experimental.upper_mz = kMaxFragmentMz;
    for (const auto& fragment : matching) {
        experimental.mz.push_back(fragment.mz);
        experimental.intensity.push_back(fragment.intensity);
    }
    psm.ms2_isotopic_abundance = 0.0;
    const auto abundance_zero = shift_predicted_fragments(psm, base, 5);
    Spectrum experimental_zero;
    experimental_zero.lower_mz = kMinFragmentMz;
    experimental_zero.upper_mz = kMaxFragmentMz;
    for (const auto& fragment : abundance_zero) {
        experimental_zero.mz.push_back(fragment.mz);
        experimental_zero.intensity.push_back(fragment.intensity);
    }
    return {
        entropy_similarity(matching, experimental, 20.0),
        entropy_similarity(abundance_zero, experimental, 20.0),
        entropy_similarity(abundance_zero, experimental_zero, 20.0),
        entropy_similarity(matching, experimental_zero, 20.0)};
}

SpectrumPredictionLibrary SpectralEntropyFeature::predict(
    const Config& config, const Dataset& unique_peptides) {
    SpectrumPredictionLibrary library;
    if (config.spectrum_paths.empty()) return library;
    if (config.spectrum_model_path.empty()) {
        throw std::runtime_error(
            "Cannot locate the DIA-NN fragmentation model; use --spectrum-model");
    }
    const std::filesystem::path model_path(config.spectrum_model_path);
    if (!std::filesystem::is_regular_file(model_path)) {
        throw std::runtime_error("DIA-NN fragmentation model does not exist: " +
                                 model_path.string());
    }
    const auto loaded = load_spectrum_library_cache(
        config, unique_peptides, *library.impl_);
    Dataset missing;
    std::unordered_map<std::string, std::size_t> missing_keys;
    missing_keys.reserve(unique_peptides.rows.size());
    for (const auto& psm : unique_peptides.rows) {
        const auto key = prediction_cache_key(psm.peptide, psm.charge);
        if (library.impl_->predictions.find(key) ==
                library.impl_->predictions.end()) {
            upsert_prediction_exemplar(missing, missing_keys, key, psm);
        }
    }
    if (!config.prediction_cache_path.empty()) {
        std::unordered_map<std::string, int> requested_labels;
        requested_labels.reserve(unique_peptides.rows.size());
        for (const auto& psm : unique_peptides.rows) {
            const auto key = prediction_cache_key(psm.peptide, psm.charge);
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
                library.impl_->predictions.find(key) !=
                library.impl_->predictions.end();
            if (label == 1) {
                reused ? ++target_reused : ++target_predict;
            } else {
                reused ? ++decoy_reused : ++decoy_predict;
            }
        }
        std::cerr << "DIA-NN spectrum cache: target reused="
                  << target_reused << ", predict=" << target_predict
                  << "; decoy reused=" << decoy_reused
                  << ", predict=" << decoy_predict << '\n';
    }
    if (missing.rows.empty()) {
        return library;
    }
    const auto prediction_begin = std::chrono::steady_clock::now();
    const std::clock_t prediction_cpu_begin = std::clock();
    std::string prediction_device;
    auto inferred = predict_fragments(
        model_path, missing, prediction_device);
    const std::clock_t prediction_cpu_end = std::clock();
    const auto prediction_end = std::chrono::steady_clock::now();
    library.impl_->timing = {
        std::chrono::duration<double>(prediction_end - prediction_begin).count(),
        static_cast<double>(prediction_cpu_end - prediction_cpu_begin) /
            CLOCKS_PER_SEC};
    std::unordered_set<std::string> target_keys;
    target_keys.reserve(missing.rows.size());
    for (const auto& psm : missing.rows) {
        if (psm.label == 1) {
            target_keys.insert(prediction_cache_key(
                psm.peptide, psm.charge));
        }
    }
    library.impl_->cache_write_timing =
        append_spectrum_library_cache(config, inferred, &target_keys);
    library.impl_->predictions.merge(inferred);
    library.impl_->device = loaded.hits == 0
        ? prediction_device : "cache+" + prediction_device;
    return library;
}

void SpectralEntropyFeature::add(
    const Config& config, Dataset& data,
    const SpectrumPredictionLibrary& library) {
    if (config.spectrum_paths.empty()) return;
    if (std::find(data.feature_names.begin(), data.feature_names.end(),
                  "unweighted_spectral_entropy") != data.feature_names.end()) {
        throw std::runtime_error(
            "PIN already contains unweighted_spectral_entropy; omit --spectra or remove the column");
    }
    const auto& predictions = library.impl_->predictions;
    data.spectrum_prediction_device = library.impl_->device;
    data.spectrum_prediction_timing = library.impl_->timing;
    data.spectrum_cache_read_timing = library.impl_->cache_read_timing;
    data.spectrum_cache_write_timing = library.impl_->cache_write_timing;
    data.feature_names.push_back("unweighted_spectral_entropy");
    data.generated_feature_names.push_back("unweighted_spectral_entropy");

    std::vector<std::vector<std::size_t>> rows_by_file(
        config.spectrum_paths.size());
    for (std::size_t i = 0; i < data.rows.size(); ++i) {
        rows_by_file[data.rows[i].file_id].push_back(i);
    }
    for (std::size_t file = 0; file < config.spectrum_paths.size(); ++file) {
        auto& file_rows = rows_by_file[file];
        std::unordered_set<std::uint64_t> requested;
        requested.reserve(file_rows.size());
        for (const auto index : file_rows) {
            requested.insert(data.rows[index].scan);
        }
        const auto spectra = load_spectra(config.spectrum_paths[file], requested);
        std::exception_ptr failure;
        if (sip_isotope_model_enabled()) {
            // A SIP prediction can contain hundreds of isotope-expanded
            // fragments. Retaining one for every peptide/charge/label key made
            // this workload hold roughly two million large vectors at once.
            // Group equal keys using only row indices, construct one expanded
            // prediction per group, score every PSM in that group, and release
            // the fragments immediately.
            std::sort(
                file_rows.begin(), file_rows.end(),
                [&](std::size_t left, std::size_t right) {
                    const auto& a = data.rows[left];
                    const auto& b = data.rows[right];
                    const auto a_body =
                        spectrum_peptide_body_view(a.peptide);
                    const auto b_body =
                        spectrum_peptide_body_view(b.peptide);
                    if (a_body != b_body) return a_body < b_body;
                    if (a.charge != b.charge) return a.charge < b.charge;
                    return a.ms2_isotopic_abundance <
                           b.ms2_isotopic_abundance;
                });
            const auto same_prediction = [&](std::size_t left,
                                             std::size_t right) {
                const auto& a = data.rows[left];
                const auto& b = data.rows[right];
                return a.charge == b.charge &&
                    a.ms2_isotopic_abundance ==
                        b.ms2_isotopic_abundance &&
                    spectrum_peptide_body_view(a.peptide) ==
                        spectrum_peptide_body_view(b.peptide);
            };
            std::vector<std::size_t> group_starts;
            group_starts.reserve(file_rows.size() / 2 + 2);
            group_starts.push_back(0);
            for (std::size_t position = 1;
                 position < file_rows.size(); ++position) {
                if (!same_prediction(
                        file_rows[position - 1], file_rows[position])) {
                    group_starts.push_back(position);
                }
            }
            group_starts.push_back(file_rows.size());

            #pragma omp parallel for schedule(dynamic, 64)
            for (std::ptrdiff_t task = 0;
                 task < static_cast<std::ptrdiff_t>(
                     group_starts.size() - 1); ++task) {
                try {
                    const std::size_t begin =
                        group_starts[static_cast<std::size_t>(task)];
                    const std::size_t end =
                        group_starts[static_cast<std::size_t>(task) + 1];
                    const auto& exemplar = data.rows[file_rows[begin]];
                    const std::string key = prediction_cache_key(
                        exemplar.peptide, exemplar.charge);
                    const auto base = predictions.find(key);
                    if (base == predictions.end()) {
                        throw std::runtime_error(
                            "Missing DIA-NN prediction for " +
                            exemplar.peptide);
                    }
                    const auto shifted = shift_predicted_fragments(
                        exemplar, base->second,
                        config.product_top_isotopes);
                    for (std::size_t position = begin;
                         position < end; ++position) {
                        auto& row = data.rows[file_rows[position]];
                        const auto spectrum = spectra.find(row.scan);
                        const float score = spectrum == spectra.end()
                            ? 0.0f
                            : entropy_similarity(
                                  shifted, spectrum->second,
                                  config.fragment_ppm);
                        row.features.push_back(score);
                    }
                } catch (...) {
                    #pragma omp critical(aerith_entropy_failure)
                    if (!failure) failure = std::current_exception();
                }
            }
        } else {
            #pragma omp parallel for schedule(dynamic, 512)
            for (std::ptrdiff_t index = 0;
                 index < static_cast<std::ptrdiff_t>(
                     file_rows.size()); ++index) {
                try {
                    auto& row = data.rows[
                        file_rows[static_cast<std::size_t>(index)]];
                    const auto spectrum = spectra.find(row.scan);
                    const std::string key = prediction_cache_key(
                        row.peptide, row.charge);
                    const auto prediction = predictions.find(key);
                    if (prediction == predictions.end()) {
                        throw std::runtime_error(
                            "Missing DIA-NN prediction for " +
                            row.peptide);
                    }
                    const float score = spectrum == spectra.end()
                        ? 0.0f
                        : entropy_similarity(
                              prediction->second, spectrum->second,
                              config.fragment_ppm);
                    row.features.push_back(score);
                } catch (...) {
                    #pragma omp critical(aerith_entropy_failure)
                    if (!failure) failure = std::current_exception();
                }
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
}

void SpectralEntropyFeature::add(const Config& config, Dataset& data) {
    if (config.spectrum_paths.empty()) return;
    auto predictions = predict(config, data);
    add(config, data, predictions);
}

#endif

} // namespace aerith
