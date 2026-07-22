#include "pipeline.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
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

void SpectralEntropyFeature::add(const Config& config, Dataset&) {
    if (!config.spectrum_paths.empty()) {
        throw std::runtime_error(
            "Aerith was built without LibTorch; rebuild with Torch_DIR to use --spectra");
    }
}

#else

namespace {

constexpr double kProton = 1.007276466621;
constexpr double kWater = 18.0105646837;
constexpr double kMinFragmentMz = 200.0;
constexpr double kMaxFragmentMz = 1800.0;
constexpr double kMinRawRelativeIntensity = 0.001;
constexpr double kPredictedBasePeakCutoff = 0.01;
constexpr std::size_t kMaximumDiannFragments = 100;
constexpr std::size_t kEntropyFragments = 20;
constexpr std::size_t kPredictionBatch = 1024;

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
};

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

std::string peptide_body(const std::string& peptide) {
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
    const std::string body = peptide_body(peptide);
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

void preprocess_for_entropy(std::vector<Fragment>& fragments) {
    merge_close_fragments(fragments);
    if (fragments.empty()) return;
    float maximum = 0.0f;
    for (const auto& fragment : fragments) maximum = std::max(maximum, fragment.intensity);
    const float cutoff = static_cast<float>(kPredictedBasePeakCutoff) * maximum;
    fragments.erase(std::remove_if(fragments.begin(), fragments.end(), [&](const Fragment& f) {
        return f.intensity < cutoff;
    }), fragments.end());
    if (fragments.size() > kEntropyFragments) {
        std::nth_element(fragments.begin(), fragments.begin() + kEntropyFragments,
                         fragments.end(), [](const Fragment& a, const Fragment& b) {
                             return a.intensity > b.intensity;
                         });
        fragments.resize(kEntropyFragments);
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

std::unordered_map<std::string, std::vector<Fragment>> predict_fragments(
    const std::filesystem::path& model_path, const Dataset& data) {
    const auto& dictionary = diann_dictionary();
    std::unordered_map<std::string, ParsedPeptide> unique;
    unique.reserve(data.rows.size() / 2);
    for (const auto& row : data.rows) {
        auto parsed = parse_peptide(row.peptide, row.charge, dictionary);
        unique.try_emplace(parsed.key, std::move(parsed));
    }

    std::unordered_map<std::size_t, std::vector<ParsedPeptide*>> groups;
    for (auto& entry : unique) groups[entry.second.tokens.size()].push_back(&entry.second);

    // This checkpoint is DIA-NN 2.6.1 Academia's models/fr.d0.pt, renamed to
    // diann-2.6.1-fragmentation.pt. Source release:
    // https://github.com/vdemichev/DiaNN/releases/tag/2.0
    // It is loaded dynamically through LibTorch.
    torch::jit::script::Module model = torch::jit::load(model_path.string());
    model.eval();
    torch::NoGradGuard no_grad;
    std::unordered_map<std::string, std::vector<Fragment>> predictions;
    predictions.reserve(unique.size());

    for (auto& group_entry : groups) {
        auto& group = group_entry.second;
        const std::size_t token_count = group_entry.first;
        for (std::size_t begin = 0; begin < group.size(); begin += kPredictionBatch) {
            const std::size_t count = std::min(kPredictionBatch, group.size() - begin);
            auto input = torch::empty(
                {static_cast<long>(count), static_cast<long>(token_count + 1)},
                torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU));
            auto accessor = input.accessor<std::int64_t, 2>();
            for (std::size_t row = 0; row < count; ++row) {
                const auto& peptide = *group[begin + row];
                accessor[row][0] = std::min(peptide.charge, 4) - 2;
                for (std::size_t column = 0; column < token_count; ++column) {
                    accessor[row][column + 1] = peptide.tokens[column];
                }
            }
            const auto output = model.forward({input}).toTensor().to(torch::kCPU).contiguous();
            const auto values = output.accessor<float, 2>();
            for (std::size_t row = 0; row < count; ++row) {
                const auto& peptide = *group[begin + row];
                const std::size_t residues = peptide.residue_masses.size();
                const std::size_t cleavages = residues - 3;
                if (static_cast<std::size_t>(output.size(1)) != 4 * cleavages) {
                    throw std::runtime_error("Unexpected DIA-NN fragmentation output shape");
                }
                float raw_maximum = 0.0f;
                for (std::size_t i = 0; i < 4 * cleavages; ++i) {
                    raw_maximum = std::max(raw_maximum, values[row][i]);
                }
                std::vector<double> prefix(residues);
                double running = peptide.nterm_shift;
                for (std::size_t i = 0; i < residues; ++i) {
                    running += peptide.residue_masses[i];
                    prefix[i] = running;
                }
                std::vector<double> suffix(residues + 1, 0.0);
                for (std::size_t i = 1; i <= residues; ++i) {
                    suffix[i] = suffix[i - 1] + peptide.residue_masses[residues - i];
                }

                struct Candidate { float mz; float raw; };
                std::vector<Candidate> candidates;
                candidates.reserve(4 * cleavages);
                for (std::size_t channel = 0; channel < 4; ++channel) {
                    const bool y_ion = channel < 2;
                    const int fragment_charge = channel % 2 == 0 ? 1 : 2;
                    for (std::size_t position = 0; position < cleavages; ++position) {
                        const float raw = values[row][channel * cleavages + position];
                        if (!(raw_maximum > 0.0f) ||
                            raw / raw_maximum < kMinRawRelativeIntensity) continue;
                        const std::size_t number = y_ion
                            ? residues - 1 - position : position + 3;
                        const double neutral = y_ion
                            ? suffix[number] + kWater : prefix[number - 1];
                        const double mz =
                            (neutral + fragment_charge * kProton) / fragment_charge;
                        if (mz >= kMinFragmentMz && mz <= kMaxFragmentMz) {
                            candidates.push_back(
                                {static_cast<float>(mz), raw});
                        }
                    }
                }
                std::sort(candidates.begin(), candidates.end(), [](const Candidate& a,
                                                                   const Candidate& b) {
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
                    fragments.push_back({candidate.mz, packed});
                }
                preprocess_for_entropy(fragments);
                predictions.emplace(peptide.key, std::move(fragments));
            }
        }
    }
    return predictions;
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

} // namespace

void SpectralEntropyFeature::add(const Config& config, Dataset& data) {
    if (config.spectrum_paths.empty()) return;
    if (config.spectrum_model_path.empty()) {
        throw std::runtime_error(
            "Cannot locate the DIA-NN fragmentation model; use --spectrum-model");
    }
    const std::filesystem::path model_path(config.spectrum_model_path);
    if (!std::filesystem::is_regular_file(model_path)) {
        throw std::runtime_error("DIA-NN fragmentation model does not exist: " +
                                 model_path.string());
    }
    if (std::find(data.feature_names.begin(), data.feature_names.end(),
                  "unweighted_spectral_entropy") != data.feature_names.end()) {
        throw std::runtime_error(
            "PIN already contains unweighted_spectral_entropy; omit --spectra or remove the column");
    }

    const auto predictions = predict_fragments(model_path, data);
    data.feature_names.push_back("unweighted_spectral_entropy");
    data.generated_feature_names.push_back("unweighted_spectral_entropy");

    for (std::size_t file = 0; file < config.spectrum_paths.size(); ++file) {
        std::unordered_set<std::uint64_t> requested;
        for (const auto& row : data.rows) {
            if (row.file_id == file) requested.insert(row.scan);
        }
        const auto spectra = load_spectra(config.spectrum_paths[file], requested);
        std::exception_ptr failure;
        #pragma omp parallel for schedule(dynamic, 512)
        for (std::ptrdiff_t index = 0;
             index < static_cast<std::ptrdiff_t>(data.rows.size()); ++index) {
            try {
                auto& row = data.rows[static_cast<std::size_t>(index)];
                if (row.file_id != file) continue;
                const std::string key = peptide_body(row.peptide) + '\x1f' +
                                        std::to_string(row.charge);
                const auto prediction = predictions.find(key);
                const auto spectrum = spectra.find(row.scan);
                if (prediction == predictions.end()) {
                    throw std::runtime_error("Missing DIA-NN prediction for " + row.peptide);
                }
                const float score = spectrum == spectra.end()
                    ? 0.0f
                    : entropy_similarity(prediction->second, spectrum->second,
                                         config.fragment_ppm);
                row.features.push_back(score);
            } catch (...) {
                #pragma omp critical(aerith_entropy_failure)
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
}

#endif

} // namespace aerith
