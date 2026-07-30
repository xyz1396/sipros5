#include "pipeline.hpp"
#include "isotope.hpp"
#include "quantification.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>

#ifdef AERITH_WITH_TORCH
#include <H5Cpp.h>
#endif

namespace aerith {

double psm_intensity(const Psm& psm) {
    if (psm.quantification_attempted) {
        return psm.has_chromatographic_feature &&
                psm.quantified_intensity > 0.0
            ? psm.quantified_intensity : 0.0;
    }
    if (!(psm.log10_precursor_intensity > 0.0) ||
        psm.log10_precursor_intensity > 300.0) {
        return 0.0;
    }
    return std::pow(10.0, psm.log10_precursor_intensity);
}

std::string ion_form(const Psm& psm) {
    std::ostringstream stream;
    const double mass = psm.calculated_mass > 0.0
        ? psm.calculated_mass : psm.exp_mass;
    stream << psm.peptide << '#' << psm.charge << '#'
           << std::fixed << std::setprecision(4) << mass;
    if (psm.sip_abundance_bin >= 0) {
        stream << "#sipbin" << psm.sip_abundance_bin;
    }
    return stream.str();
}

void update_ion_intensity(
    IonIntensityMap& intensities, const std::string& ion, double intensity) {
    auto& current = intensities[ion];
    current = std::max(current, intensity);
}

double top_three_intensity(const IonIntensityMap& intensities) {
    double first = 0.0;
    double second = 0.0;
    double third = 0.0;
    for (const auto& [ion, intensity] : intensities) {
        (void)ion;
        if (intensity >= first) {
            third = second;
            second = first;
            first = intensity;
        } else if (intensity >= second) {
            third = second;
            second = intensity;
        } else if (intensity > third) {
            third = intensity;
        }
    }
    return first + second + third;
}

double summed_intensity(const IonIntensityMap& intensities) {
    return std::accumulate(
        intensities.begin(), intensities.end(), 0.0,
        [](double total, const auto& item) { return total + item.second; });
}

double summed_isotope_apex_intensity(
    const std::array<double, 3>& isotope_apices) {
    return std::accumulate(
        isotope_apices.begin(), isotope_apices.end(), 0.0);
}

double summed_isotope_apex_intensity(
    const std::vector<double>& isotope_apices) {
    // Automatic conventional-LC intensity sums the background-corrected
    // apex-1/apex/apex+1 signal over the traced isotope envelope. Integrated
    // chromatographic area remains available explicitly as intensity mode 1.
    return std::accumulate(
        isotope_apices.begin(), isotope_apices.end(), 0.0);
}

double quantification_median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    const auto middle = values.begin() +
        static_cast<std::ptrdiff_t>(values.size() / 2);
    std::nth_element(values.begin(), middle, values.end());
    const double upper = *middle;
    if (values.size() % 2 != 0) return upper;
    const auto lower = std::max_element(values.begin(), middle);
    return 0.5 * (upper + *lower);
}

double quantification_mz_from_ion_key(const std::string& key) {
    const auto sip_separator = key.rfind("#sipbin");
    const std::string_view base(
        key.data(), sip_separator == std::string::npos
            ? key.size() : sip_separator);
    const auto mass_separator = base.rfind('#');
    if (mass_separator == std::string::npos) return 0.0;
    const auto charge_separator = base.rfind('#', mass_separator - 1);
    if (charge_separator == std::string::npos) return 0.0;
    try {
        const double mass = std::stod(std::string(
            base.substr(mass_separator + 1)));
        const int charge = std::stoi(std::string(base.substr(
            charge_separator + 1,
            mass_separator - charge_separator - 1)));
        constexpr double proton = 1.007276466621;
        return charge <= 0 ? 0.0
            : (mass + static_cast<double>(charge) * proton) /
                  static_cast<double>(charge);
    } catch (const std::exception&) {
        return 0.0;
    }
}

double IntensityNormalizer::factor(
    std::size_t sample, const std::string& ion) const {
    if (sample >= factors.size()) return 1.0;
    const auto found = ion_bins.find(ion);
    const std::size_t bin = found == ion_bins.end()
        ? static_cast<std::size_t>(std::distance(
              upper_mz.begin(),
              std::lower_bound(
                  upper_mz.begin(), upper_mz.end(),
                  quantification_mz_from_ion_key(ion))))
        : found->second;
    return factors[sample][std::min(bins - 1, bin)];
}

IntensityNormalizer build_intensity_normalizer(
    const std::vector<IonIntensityMap>& sample_ions,
    std::vector<std::pair<double, std::string>> ion_mz,
    bool enabled) {
    IntensityNormalizer result;
    result.upper_mz.fill(std::numeric_limits<double>::infinity());
    result.factors.resize(sample_ions.size());
    for (auto& values : result.factors) values.fill(1.0);
    if (!enabled || sample_ions.size() < 2) return result;

    std::sort(ion_mz.begin(), ion_mz.end());
    for (std::size_t index = 0; index < ion_mz.size(); ++index) {
        result.ion_bins[ion_mz[index].second] = std::min(
            IntensityNormalizer::bins - 1,
            index * IntensityNormalizer::bins /
                std::max<std::size_t>(1, ion_mz.size()));
    }
    for (std::size_t bin = 0;
         bin + 1 < IntensityNormalizer::bins && !ion_mz.empty(); ++bin) {
        const auto end = std::min(
            ion_mz.size() - 1,
            (bin + 1) * ion_mz.size() / IntensityNormalizer::bins);
        result.upper_mz[bin] = ion_mz[end].first;
    }

    std::vector<std::pair<double, std::size_t>> totals;
    for (std::size_t sample = 0; sample < sample_ions.size(); ++sample) {
        totals.emplace_back(summed_intensity(sample_ions[sample]), sample);
    }
    std::sort(totals.begin(), totals.end());
    const std::size_t reference = totals[totals.size() / 2].second;

    for (std::size_t sample = 0; sample < sample_ions.size(); ++sample) {
        if (sample == reference) continue;
        std::array<std::vector<double>, IntensityNormalizer::bins> ratios;
        std::vector<double> overall;
        for (const auto& [key, intensity] : sample_ions[sample]) {
            const auto reference_ion = sample_ions[reference].find(key);
            if (reference_ion == sample_ions[reference].end() ||
                !(intensity > 0.0) || !(reference_ion->second > 0.0)) {
                continue;
            }
            const double ratio =
                std::log(intensity / reference_ion->second);
            ratios[result.ion_bins[key]].push_back(ratio);
            overall.push_back(ratio);
        }
        const double fallback = overall.empty()
            ? 0.0 : quantification_median(overall);
        for (std::size_t bin = 0; bin < IntensityNormalizer::bins; ++bin) {
            const double log_ratio = ratios[bin].empty()
                ? fallback : quantification_median(ratios[bin]);
            result.factors[sample][bin] = std::exp(-log_ratio);
        }
    }
    return result;
}

double normalized_intensity(
    const IonIntensityMap& ions, std::size_t sample,
    const IntensityNormalizer& normalizer) {
    double result = 0.0;
    for (const auto& [key, intensity] : ions) {
        result += intensity * normalizer.factor(sample, key);
    }
    return result;
}

double normalized_top_three(
    const IonIntensityMap& ions, std::size_t sample,
    const IntensityNormalizer& normalizer) {
    std::vector<double> values;
    values.reserve(ions.size());
    for (const auto& [key, intensity] : ions) {
        values.push_back(intensity * normalizer.factor(sample, key));
    }
    std::sort(values.begin(), values.end(), std::greater<double>());
    return std::accumulate(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(
            std::min<std::size_t>(3, values.size())),
        0.0);
}

std::vector<double> solve_quantification_linear_system(
    std::vector<std::vector<double>> matrix,
    std::vector<double> values) {
    const std::size_t size = values.size();
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < size; ++row) {
            if (std::abs(matrix[row][column]) >
                std::abs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][column]) < 1e-12) continue;
        std::swap(matrix[pivot], matrix[column]);
        std::swap(values[pivot], values[column]);
        const double divisor = matrix[column][column];
        for (std::size_t entry = column; entry < size; ++entry) {
            matrix[column][entry] /= divisor;
        }
        values[column] /= divisor;
        for (std::size_t row = 0; row < size; ++row) {
            if (row == column) continue;
            const double multiplier = matrix[row][column];
            if (multiplier == 0.0) continue;
            for (std::size_t entry = column; entry < size; ++entry) {
                matrix[row][entry] -=
                    multiplier * matrix[column][entry];
            }
            values[row] -= multiplier * values[column];
        }
    }
    return values;
}

std::vector<double> maxlfq(
    const std::vector<IonIntensityMap>& sample_ions,
    const IntensityNormalizer& normalizer) {
    const std::size_t sample_count = sample_ions.size();
    std::vector<double> result(sample_count, 0.0);
    std::set<std::string> ion_keys;
    for (const auto& ions : sample_ions) {
        for (const auto& [key, intensity] : ions) {
            if (intensity > 0.0) ion_keys.insert(key);
        }
    }
    std::vector<std::vector<double>> intensities(
        ion_keys.size(), std::vector<double>(sample_count, 0.0));
    std::size_t ion_index = 0;
    for (const auto& key : ion_keys) {
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            const auto found = sample_ions[sample].find(key);
            if (found != sample_ions[sample].end()) {
                intensities[ion_index][sample] =
                    found->second * normalizer.factor(sample, key);
            }
        }
        ++ion_index;
    }

    std::vector<std::vector<std::size_t>> graph(sample_count);
    struct Ratio {
        std::size_t first;
        std::size_t second;
        double value;
    };
    std::vector<Ratio> ratios;
    for (std::size_t first = 0; first < sample_count; ++first) {
        for (std::size_t second = first + 1;
             second < sample_count; ++second) {
            std::vector<double> values;
            for (const auto& ion : intensities) {
                if (ion[first] > 0.0 && ion[second] > 0.0) {
                    values.push_back(
                        std::log(ion[first]) - std::log(ion[second]));
                }
            }
            if (values.empty()) continue;
            ratios.push_back(
                {first, second, quantification_median(std::move(values))});
            graph[first].push_back(second);
            graph[second].push_back(first);
        }
    }

    std::vector<bool> visited(sample_count, false);
    for (std::size_t seed = 0; seed < sample_count; ++seed) {
        if (visited[seed]) continue;
        bool has_intensity = false;
        for (const auto& ion : intensities) {
            has_intensity |= ion[seed] > 0.0;
        }
        if (!has_intensity) {
            visited[seed] = true;
            continue;
        }
        std::vector<std::size_t> component{seed};
        visited[seed] = true;
        for (std::size_t position = 0; position < component.size();
             ++position) {
            for (const auto neighbor : graph[component[position]]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    component.push_back(neighbor);
                }
            }
        }
        if (component.size() == 1) {
            for (const auto& ion : intensities) {
                result[seed] = std::max(result[seed], ion[seed]);
            }
            continue;
        }
        std::unordered_map<std::size_t, std::size_t> local;
        for (std::size_t index = 0; index < component.size(); ++index) {
            local[component[index]] = index;
        }
        std::vector<std::vector<double>> matrix(
            component.size(), std::vector<double>(component.size(), 0.0));
        std::vector<double> values(component.size(), 0.0);
        for (const auto& ratio : ratios) {
            if (local.count(ratio.first) == 0 ||
                local.count(ratio.second) == 0) {
                continue;
            }
            const auto first = local.at(ratio.first);
            const auto second = local.at(ratio.second);
            matrix[first][first] += 1.0;
            matrix[second][second] += 1.0;
            matrix[first][second] -= 1.0;
            matrix[second][first] -= 1.0;
            values[first] += ratio.value;
            values[second] -= ratio.value;
        }
        std::fill(matrix.back().begin(), matrix.back().end(), 1.0);
        values.back() = 0.0;
        auto profile =
            solve_quantification_linear_system(matrix, values);

        double anchor = 0.0;
        for (const auto& ion : intensities) {
            double log_sum = 0.0;
            std::size_t observed = 0;
            for (const auto sample : component) {
                if (ion[sample] > 0.0) {
                    log_sum += std::log(ion[sample]);
                    ++observed;
                }
            }
            if (observed == component.size()) {
                anchor = std::max(
                    anchor, std::exp(log_sum / observed));
            }
        }
        if (!(anchor > 0.0)) {
            for (const auto& ion : intensities) {
                for (const auto sample : component) {
                    anchor = std::max(anchor, ion[sample]);
                }
            }
        }
        const double profile_mean = std::accumulate(
            profile.begin(), profile.end(), 0.0) /
            static_cast<double>(profile.size());
        const double shift = std::log(anchor) - profile_mean;
        for (std::size_t index = 0; index < component.size(); ++index) {
            result[component[index]] =
                std::exp(profile[index] + shift);
        }
    }
    return result;
}

#ifndef AERITH_WITH_TORCH

QuantificationResult ChromatographicQuantifier::add(
    const Config& config, Dataset&, const std::vector<double>&) {
    if (!config.spectrum_paths.empty()) {
        throw std::runtime_error(
            "Aerith was built without HDF5 support; chromatographic "
            "quantification requires a Conda/LibTorch build");
    }
    return {};
}

#else

constexpr double kQuantProton = 1.007276466621;
constexpr double kC13Spacing = 1.0033548378;
constexpr double kQuantWater = 18.0105646837;

struct QuantScan {
    std::uint64_t scan = 0;
    double retention = 0.0;
    std::size_t peak_start = 0;
    std::size_t peak_count = 0;
};

struct QuantMs1Data {
    std::vector<QuantScan> scans;
    std::vector<double> mz;
    std::vector<double> intensity;
    std::unordered_map<std::uint64_t, std::uint64_t> parent_scans;
};

double quant_residue_mass(char residue) {
    switch (residue) {
    case 'A': return 71.037113805;
    case 'R': return 156.10111105;
    case 'N': return 114.04292747;
    case 'D': return 115.026943065;
    case 'C': return 103.009184505;
    case 'E': return 129.042593135;
    case 'Q': return 128.05857754;
    case 'G': return 57.021463735;
    case 'H': return 137.058911875;
    case 'I':
    case 'L': return 113.084063975;
    case 'K': return 128.094963015;
    case 'M': return 131.040484645;
    case 'F': return 147.068413945;
    case 'P': return 97.052763875;
    case 'S': return 87.032028435;
    case 'T': return 101.047678505;
    case 'W': return 186.07931298;
    case 'Y': return 163.063328575;
    case 'V': return 99.068413945;
    default: return 0.0;
    }
}

double quant_modification_mass(char token) {
    switch (token) {
    case '~': return 15.994915;
    case '!': return 0.984016;
    case '@':
    case '>':
    case '<': return 79.966332;
    case '%': return 42.010565;
    case '^': return 14.015650;
    case '&': return 28.031300;
    case '*': return 42.046950;
    case '(': return 28.990164;
    case ')': return 44.985079;
    case '/': return 57.021464;
    case '$': return 45.987721;
    default: return 0.0;
    }
}

double quant_theoretical_mass(
    const std::string& decorated, bool fixed_cam) {
    const auto open = decorated.find('[');
    const auto close = decorated.rfind(']');
    const std::string body =
        open != std::string::npos && close != std::string::npos &&
                close > open
            ? decorated.substr(open + 1, close - open - 1)
            : decorated;
    struct Residue {
        char amino_acid = '\0';
        std::vector<char> modifications;
    };
    std::vector<Residue> residues;
    std::vector<char> terminal_modifications;
    for (const char value : body) {
        const char residue = static_cast<char>(
            std::toupper(static_cast<unsigned char>(value)));
        if (quant_residue_mass(residue) > 0.0) {
            residues.push_back({residue, {}});
        } else if (quant_modification_mass(value) != 0.0) {
            if (residues.empty()) terminal_modifications.push_back(value);
            else residues.back().modifications.push_back(value);
        }
    }

    double mass = kQuantWater;
    for (const char token : terminal_modifications) {
        mass += quant_modification_mass(token);
    }
    for (const auto& residue : residues) {
        mass += quant_residue_mass(residue.amino_acid);
        const bool replaces_cam =
            std::find(
                residue.modifications.begin(),
                residue.modifications.end(), '(') !=
                residue.modifications.end() ||
            std::find(
                residue.modifications.begin(),
                residue.modifications.end(), '/') !=
                residue.modifications.end();
        if (fixed_cam && residue.amino_acid == 'C' && !replaces_cam) {
            mass += 57.021464;
        }
        for (const char token : residue.modifications) {
            mass += quant_modification_mass(token);
        }
    }
    return mass;
}

template <typename T>
std::vector<T> quant_read_1d(
    H5::H5File& file, const std::string& name,
    const H5::PredType& type) {
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1) {
        throw std::runtime_error(
            "HDF5 dataset is not one-dimensional: " + name);
    }
    hsize_t size = 0;
    space.getSimpleExtentDims(&size, nullptr);
    std::vector<T> values(static_cast<std::size_t>(size));
    if (size != 0) dataset.read(values.data(), type);
    return values;
}

template <typename T>
std::vector<T> quant_read_slice(
    H5::DataSet& dataset, const H5::PredType& type,
    hsize_t start, hsize_t count) {
    std::vector<T> result(static_cast<std::size_t>(count));
    if (count == 0) return result;
    H5::DataSpace file_space = dataset.getSpace();
    hsize_t extent = 0;
    file_space.getSimpleExtentDims(&extent, nullptr);
    if (start > extent || count > extent - start) {
        throw std::runtime_error("HDF5 MS1 peak slice is out of bounds");
    }
    H5::DataSpace memory_space(1, &count);
    file_space.selectHyperslab(H5S_SELECT_SET, &count, &start);
    dataset.read(result.data(), type, memory_space, file_space);
    return result;
}

QuantMs1Data load_quant_ms1(const std::string& path) {
    H5::Exception::dontPrint();
    H5::H5File file(path, H5F_ACC_RDONLY);
    const auto scan_numbers = quant_read_1d<int>(
        file, "/scans/scan_number", H5::PredType::NATIVE_INT);
    const auto ms_order = quant_read_1d<int>(
        file, "/scans/ms_order", H5::PredType::NATIVE_INT);
    const auto retention = quant_read_1d<double>(
        file, "/scans/retention_time", H5::PredType::NATIVE_DOUBLE);
    const auto peak_start = quant_read_1d<long long>(
        file, "/scans/peak_start", H5::PredType::NATIVE_LLONG);
    const auto peak_count = quant_read_1d<int>(
        file, "/scans/peak_count", H5::PredType::NATIVE_INT);
    const auto parent_scan = quant_read_1d<int>(
        file, "/scans/parent_scan_number", H5::PredType::NATIVE_INT);
    if (ms_order.size() != scan_numbers.size() ||
        retention.size() != scan_numbers.size() ||
        peak_start.size() != scan_numbers.size() ||
        peak_count.size() != scan_numbers.size() ||
        parent_scan.size() != scan_numbers.size()) {
        throw std::runtime_error(
            "Raxport HDF5 scan datasets have inconsistent lengths");
    }

    QuantMs1Data result;
    std::size_t total_peaks = 0;
    for (std::size_t index = 0; index < scan_numbers.size(); ++index) {
        if (ms_order[index] == 1 && peak_start[index] >= 0 &&
            peak_count[index] > 0) {
            total_peaks += static_cast<std::size_t>(peak_count[index]);
        }
        if (parent_scan[index] > 0) {
            result.parent_scans.emplace(
                static_cast<std::uint64_t>(scan_numbers[index]),
                static_cast<std::uint64_t>(parent_scan[index]));
        }
    }
    result.scans.reserve(
        static_cast<std::size_t>(
            std::count(ms_order.begin(), ms_order.end(), 1)));
    result.mz.reserve(total_peaks);
    result.intensity.reserve(total_peaks);
    H5::DataSet mz_dataset = file.openDataSet("/peaks/mz");
    H5::DataSet intensity_dataset = file.openDataSet("/peaks/intensity");
    for (std::size_t index = 0; index < scan_numbers.size(); ++index) {
        if (ms_order[index] != 1 || peak_start[index] < 0 ||
            peak_count[index] <= 0) {
            continue;
        }
        const auto mz = quant_read_slice<double>(
            mz_dataset, H5::PredType::NATIVE_DOUBLE,
            static_cast<hsize_t>(peak_start[index]),
            static_cast<hsize_t>(peak_count[index]));
        const auto intensity = quant_read_slice<double>(
            intensity_dataset, H5::PredType::NATIVE_DOUBLE,
            static_cast<hsize_t>(peak_start[index]),
            static_cast<hsize_t>(peak_count[index]));
        const auto output_start = result.mz.size();
        result.mz.insert(result.mz.end(), mz.begin(), mz.end());
        result.intensity.insert(
            result.intensity.end(), intensity.begin(), intensity.end());
        result.scans.push_back({
            static_cast<std::uint64_t>(scan_numbers[index]),
            retention[index], output_start, mz.size()});
    }
    return result;
}

struct QuantPeak {
    double intensity = 0.0;
    double mz = 0.0;
};

QuantPeak trace_peak(
    const QuantMs1Data& data, const QuantScan& scan,
    double target_mz, double ppm) {
    const auto mz_begin =
        data.mz.begin() + static_cast<std::ptrdiff_t>(scan.peak_start);
    const auto mz_end =
        mz_begin + static_cast<std::ptrdiff_t>(scan.peak_count);
    const double tolerance = target_mz * ppm * 1e-6;
    const auto begin =
        std::lower_bound(mz_begin, mz_end, target_mz - tolerance);
    QuantPeak result;
    for (auto peak = begin; peak != mz_end && *peak <= target_mz + tolerance;
         ++peak) {
        const auto index = static_cast<std::size_t>(
            std::distance(data.mz.begin(), peak));
        if (data.intensity[index] > result.intensity) {
            result.intensity = data.intensity[index];
            result.mz = *peak;
        }
    }
    return result;
}

struct ResampledChromatogram {
    std::vector<double> retention;
    std::vector<std::vector<double>> isotope;
};

ResampledChromatogram resample_chromatogram(
    const std::vector<QuantScan>& scans, std::size_t offset,
    const std::vector<std::vector<double>>& isotope) {
    ResampledChromatogram result;
    if (isotope.empty() || isotope.front().empty()) return result;
    result.isotope.resize(isotope.size());
    const std::size_t count = isotope.front().size();
    if (count == 1) {
        result.retention.push_back(scans[offset].retention);
        for (std::size_t mass = 0; mass < result.isotope.size(); ++mass) {
            result.isotope[mass].push_back(isotope[mass].front());
        }
        return result;
    }
    const double first = scans[offset].retention;
    const double last = scans[offset + count - 1].retention;
    // Keep the original point count and round the regular-grid spacing upward
    // to a whole millisecond for deterministic interpolation.
    const double step = std::ceil(
        std::max(0.0, last - first) * 60000.0 /
        static_cast<double>(count - 1)) / 60000.0;
    if (!(step > 0.0)) {
        result.retention.push_back(first);
        for (std::size_t mass = 0; mass < result.isotope.size(); ++mass) {
            result.isotope[mass].push_back(isotope[mass].front());
        }
        return result;
    }
    result.retention.resize(count);
    for (auto& trace : result.isotope) {
        trace.resize(count, 0.0);
    }
    std::size_t source = 0;
    for (std::size_t index = 0; index < count; ++index) {
        const double rt = std::min(last, first + index * step);
        result.retention[index] = rt;
        while (source + 1 < count &&
               scans[offset + source + 1].retention < rt) {
            ++source;
        }
        const std::size_t next = std::min(count - 1, source + 1);
        const double left_rt = scans[offset + source].retention;
        const double right_rt = scans[offset + next].retention;
        const double fraction = right_rt > left_rt
            ? std::clamp((rt - left_rt) / (right_rt - left_rt), 0.0, 1.0)
            : 0.0;
        for (std::size_t mass = 0; mass < result.isotope.size(); ++mass) {
            result.isotope[mass][index] =
                isotope[mass][source] +
                fraction * (isotope[mass][next] - isotope[mass][source]);
        }
    }
    return result;
}

std::vector<double> savitzky_golay_smooth(
    const std::vector<double>& values) {
    // Five-point, second-order Savitzky-Golay kernel.
    static constexpr double kernel[] = {
        -3.0 / 35.0, 12.0 / 35.0, 17.0 / 35.0,
        12.0 / 35.0, -3.0 / 35.0};
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t index = 0; index < values.size(); ++index) {
        for (int shift = -2; shift <= 2; ++shift) {
            const auto position =
                static_cast<std::ptrdiff_t>(index) + shift;
            if (position < 0 ||
                position >= static_cast<std::ptrdiff_t>(values.size())) {
                continue;
            }
            result[index] += values[static_cast<std::size_t>(position)] *
                kernel[shift + 2];
        }
        result[index] = std::max(0.0, result[index]);
    }
    return result;
}

struct ChromatographicPeak {
    std::size_t apex = 0;
    std::size_t left = 0;
    std::size_t right = 0; // exclusive
};

bool derivative_valley(
    const std::vector<double>& derivative, std::size_t index) {
    return (derivative[index - 1] < 0.0 && derivative[index] >= 0.0) ||
        (derivative[index - 1] <= 0.0 && derivative[index] > 0.0);
}

bool derivative_peak(
    const std::vector<double>& derivative, std::size_t index) {
    return (derivative[index - 1] > 0.0 && derivative[index] <= 0.0) ||
        (derivative[index - 1] >= 0.0 && derivative[index] < 0.0);
}

bool detect_chromatographic_peak(
    const std::vector<double>& raw, const std::vector<double>& retention,
    std::size_t anchor, double search_start, double search_end,
    ChromatographicPeak& result) {
    if (raw.size() < 3 || raw.size() != retention.size() ||
        anchor >= raw.size()) {
        return false;
    }
    const auto smooth = savitzky_golay_smooth(raw);
    std::vector<double> derivative(smooth.size() - 1, 0.0);
    for (std::size_t index = 0; index < derivative.size(); ++index) {
        derivative[index] = smooth[index + 1] - smooth[index];
    }

    std::vector<std::size_t> maxima;
    if (smooth[anchor] > 0.0) maxima.push_back(anchor);
    double neighboring_peak = smooth[anchor];
    std::size_t right_valley = smooth.size() - 1;
    if (derivative.size() >= 4 && anchor <= derivative.size() - 4) {
        for (std::size_t index = anchor + 1;
             index + 1 < derivative.size(); ++index) {
            if (derivative_valley(derivative, index)) {
                for (std::size_t next = index + 1;
                     next < derivative.size(); ++next) {
                    if (!derivative_peak(derivative, next)) continue;
                    if (smooth[index] <
                        0.5 * std::min(neighboring_peak, smooth[next])) {
                        right_valley = index;
                    }
                    break;
                }
            } else if (derivative_peak(derivative, index)) {
                maxima.push_back(index);
                neighboring_peak = smooth[index];
            }
            if (right_valley + 1 < smooth.size()) break;
        }
    }

    std::size_t left_valley = 0;
    if (anchor >= 4) {
        for (std::size_t index = anchor - 1; index > 1; --index) {
            if (derivative_valley(derivative, index)) {
                for (std::size_t previous = index - 1;
                     previous > 0; --previous) {
                    if (!derivative_peak(derivative, previous)) continue;
                    if (smooth[index] <
                        0.5 * std::min(
                            neighboring_peak, smooth[previous])) {
                        left_valley = index;
                    }
                    break;
                }
            } else if (derivative_peak(derivative, index)) {
                maxima.push_back(index);
                neighboring_peak = smooth[index];
            }
            if (left_valley > 0) break;
        }
    }

    std::size_t apex = smooth.size();
    double apex_value = -1.0;
    for (const auto candidate : maxima) {
        if (retention[candidate] <= search_start ||
            retention[candidate] >= search_end ||
            smooth[candidate] <= apex_value) {
            continue;
        }
        apex = candidate;
        apex_value = smooth[candidate];
    }
    if (apex == smooth.size()) return false;

    // Preserve low-baseline Orbitrap features by using the smaller of 1% of
    // the apex and one intensity unit.
    const double floor = std::min(apex_value * 0.01, 1.0);
    std::size_t right = apex;
    while (right < right_valley && smooth[right] >= floor) ++right;
    right = std::min(right, right_valley);
    std::size_t left = apex;
    while (left > left_valley && smooth[left] >= floor) --left;
    left = std::max(left, left_valley);
    if (right <= left + 1 || anchor < left || anchor > right) return false;

    // Adjust the smoothed apex against its two immediate raw neighbors while
    // preserving deterministic right-neighbor tie/order behavior.
    std::size_t raw_apex = apex;
    const double original_apex = raw[apex];
    if (apex > left && raw[apex - 1] > original_apex) {
        raw_apex = apex - 1;
    }
    if (apex + 1 <= right && apex + 1 < raw.size() &&
        raw[apex + 1] > original_apex) {
        raw_apex = apex + 1;
    }
    result = {raw_apex, left, std::min(raw.size(), right + 1)};
    return true;
}

double background_corrected_chromatographic_area(
    const std::vector<double>& retention, const std::vector<double>& trace,
    std::size_t left, std::size_t right) {
    if (trace.empty() || trace.size() != retention.size() ||
        right <= left + 1 || right > trace.size()) {
        return 0.0;
    }
    const double background = std::min(trace[left], trace[right - 1]);
    double raw_area = 0.0;
    double background_area = 0.0;
    for (std::size_t index = left; index + 1 < right; ++index) {
        const double elapsed_seconds =
            std::max(0.0, retention[index + 1] - retention[index]) * 60.0;
        raw_area +=
            0.5 * (trace[index] + trace[index + 1]) * elapsed_seconds;
        background_area += 0.5 *
            (std::min(background, trace[index]) +
             std::min(background, trace[index + 1])) * elapsed_seconds;
    }
    return std::max(0.0, raw_area - background_area);
}

double background_corrected_chromatographic_apex(
    const std::vector<double>& trace, std::size_t apex,
    std::size_t left, std::size_t right) {
    if (trace.empty() || left >= right || right > trace.size() ||
        apex < left || apex >= right) {
        return 0.0;
    }
    const double background = std::min(trace[left], trace[right - 1]);
    double intensity = 0.0;
    for (int shift = -1; shift <= 1; ++shift) {
        const auto index = static_cast<std::ptrdiff_t>(apex) + shift;
        if (index < static_cast<std::ptrdiff_t>(left) ||
            index >= static_cast<std::ptrdiff_t>(right)) {
            continue;
        }
        intensity += std::max(
            0.0, trace[static_cast<std::size_t>(index)] - background);
    }
    return intensity;
}

double interpolated_half_width(
    const std::vector<double>& retention,
    const std::vector<double>& trace, std::size_t apex,
    std::size_t left, std::size_t right) {
    if (trace.empty() || trace[apex] <= 0.0) return 0.0;
    const double half = trace[apex] * 0.5;
    double left_rt = retention[apex];
    for (std::size_t index = apex; index > left; --index) {
        const double current = trace[index];
        const double previous = trace[index - 1];
        if (previous <= half && current >= half) {
            const double fraction = current == previous
                ? 0.0 : (half - previous) / (current - previous);
            left_rt = retention[index - 1] +
                fraction *
                    (retention[index] - retention[index - 1]);
            break;
        }
        left_rt = retention[index - 1];
    }
    double right_rt = retention[apex];
    for (std::size_t index = apex; index < right; ++index) {
        const double current = trace[index];
        const double next = trace[index + 1];
        if (current >= half && next <= half) {
            const double fraction = current == next
                ? 0.0 : (half - current) / (next - current);
            right_rt = retention[index] +
                fraction *
                    (retention[index + 1] - retention[index]);
            break;
        }
        right_rt = retention[index + 1];
    }
    return std::max(0.0, (right_rt - left_rt) * 60.0);
}

void quantify_psm(
    const Config& config, const QuantMs1Data& ms1, Psm& psm,
    double rt_window = -1.0, double mz_shift = 0.0) {
    psm.quantification_attempted = true;
    psm.has_chromatographic_feature = false;
    psm.quantified_intensity = 0.0;
    psm.apex_retention = 0.0;
    psm.retention_start = 0.0;
    psm.retention_end = 0.0;
    psm.retention_fwhm = 0.0;
    psm.apex_scan = 0;
    psm.traced_scans = 0;
    psm.quant_mass_error_ppm = 0.0;
    psm.quant_isotope_kl = 0.0;
    psm.quant_isotope_correlation = 0.0;
    psm.quant_isotope_fraction = 0.0;
    psm.quant_isotope_apex_spread = 0.0;
    const auto parent = ms1.parent_scans.find(psm.scan);
    if (parent != ms1.parent_scans.end()) psm.parent_scan = parent->second;
    if (psm.charge <= 0 || ms1.scans.empty()) return;
    psm.calculated_mass =
        quant_theoretical_mass(psm.peptide, config.fixed_cam);
    const double unshifted_mz =
        (psm.calculated_mass +
         static_cast<double>(psm.charge) * kQuantProton) /
        static_cast<double>(psm.charge);
    psm.calculated_mz = unshifted_mz;
    struct QuantTheoreticalPeak {
        double mz = 0.0;
        double probability = 0.0;
    };
    std::vector<QuantTheoreticalPeak> theoretical_peaks;
    if (sip_isotope_model_enabled()) {
        const auto peaks = precursor_isotope_peaks(
            psm.peptide, psm.ms2_isotopic_abundance,
            config.quant_top_isotopes);
        theoretical_peaks.reserve(peaks.size());
        for (const auto& peak : peaks) {
            theoretical_peaks.push_back({
                (peak.neutral_mass +
                 static_cast<double>(psm.charge) * kQuantProton) /
                    static_cast<double>(psm.charge) +
                    mz_shift,
                peak.probability});
        }
    } else {
        const double precursor_mz = unshifted_mz + mz_shift;
        const double lambda =
            std::max(0.05, psm.calculated_mass * 0.00049);
        theoretical_peaks = {
            {precursor_mz, 1.0},
            {precursor_mz + kC13Spacing /
                static_cast<double>(psm.charge), lambda},
            {precursor_mz + 2.0 * kC13Spacing /
                static_cast<double>(psm.charge),
             0.5 * lambda * lambda}};
    }
    const double active_window =
        rt_window > 0.0 ? rt_window : config.quant_rt_window;
    // MBR uses a narrow aligned target region for apex selection, but peak
    // boundaries may extend beyond that region.
    // Permit detected peak boundaries to extend beyond the target RT region.
    // Keep one normal tracing window on either side of the apex search region
    // so broad features are not truncated at the alignment tolerance.
    const double extraction_window =
        std::max(active_window, config.quant_rt_window) +
        config.quant_rt_window;
    const double first_rt = psm.retention - extraction_window;
    const double last_rt = psm.retention + extraction_window;
    const auto first = std::lower_bound(
        ms1.scans.begin(), ms1.scans.end(), first_rt,
        [](const QuantScan& scan, double rt) {
            return scan.retention < rt;
        });
    const auto last = std::upper_bound(
        first, ms1.scans.end(), last_rt,
        [](double rt, const QuantScan& scan) {
            return rt < scan.retention;
        });
    if (first == last) return;
    const auto offset =
        static_cast<std::size_t>(std::distance(ms1.scans.begin(), first));
    const auto count =
        static_cast<std::size_t>(std::distance(first, last));
    std::vector<std::vector<double>> isotope(
        theoretical_peaks.size(), std::vector<double>(count, 0.0));
    std::vector<double> anchor_mz(count, 0.0);
    for (std::size_t scan = 0; scan < count; ++scan) {
        for (std::size_t mass = 0; mass < isotope.size(); ++mass) {
            const auto peak = trace_peak(
                ms1, ms1.scans[offset + scan],
                theoretical_peaks[mass].mz,
                config.quant_mz_ppm);
            isotope[mass][scan] = peak.intensity;
            if (mass == 0) anchor_mz[scan] = peak.mz;
        }
    }
    const auto chromatogram =
        resample_chromatogram(ms1.scans, offset, isotope);
    if (chromatogram.retention.empty()) return;

    std::vector<double> theoretical;
    theoretical.reserve(theoretical_peaks.size());
    for (const auto& peak : theoretical_peaks) {
        theoretical.push_back(peak.probability);
    }
    // Regular quantification anchors to M+0. SIP quantification orders the
    // exact source-aware precursor distribution by probability, so index zero
    // is the theoretical envelope apex. Other isotope traces validate the
    // envelope and contribute intensity without moving the reported apex.
    constexpr std::size_t base_isotope = 0;
    const auto anchor_position = std::lower_bound(
        chromatogram.retention.begin(), chromatogram.retention.end(),
        psm.retention);
    std::size_t anchor = anchor_position == chromatogram.retention.end()
        ? chromatogram.retention.size() - 1
        : static_cast<std::size_t>(
              std::distance(
                  chromatogram.retention.begin(), anchor_position));
    if (anchor > 0 &&
        std::abs(chromatogram.retention[anchor - 1] - psm.retention) <
            std::abs(chromatogram.retention[anchor] - psm.retention)) {
        --anchor;
    }
    ChromatographicPeak base_peak;
    if (!detect_chromatographic_peak(
            chromatogram.isotope[base_isotope], chromatogram.retention,
            anchor, psm.retention - active_window,
            psm.retention + active_window, base_peak)) {
        return;
    }
    const std::size_t left = base_peak.left;
    const std::size_t right = base_peak.right;
    const std::size_t apex = base_peak.apex;
    const auto base_smooth =
        savitzky_golay_smooth(chromatogram.isotope[base_isotope]);

    std::size_t traced = 0;
    std::size_t isotope_count = 0;
    for (std::size_t scan = left; scan < right; ++scan) {
        if (chromatogram.isotope[base_isotope][scan] > 0.0) ++traced;
    }

    std::vector<double> isotope_totals(
        chromatogram.isotope.size(), 0.0);
    std::vector<double> isotope_apices(
        chromatogram.isotope.size(), 0.0);
    std::vector<std::size_t> isotope_peak_apex(
        chromatogram.isotope.size(), apex);
    std::vector<bool> isotope_detected(
        chromatogram.isotope.size(), false);
    for (std::size_t mass = 0;
         mass < chromatogram.isotope.size(); ++mass) {
        ChromatographicPeak peak = base_peak;
        if (mass != base_isotope &&
            !detect_chromatographic_peak(
                chromatogram.isotope[mass], chromatogram.retention, apex,
                chromatogram.retention[left],
                chromatogram.retention[right - 1], peak)) {
            continue;
        }
        std::size_t isotope_scans = 0;
        for (std::size_t scan = peak.left; scan < peak.right; ++scan) {
            if (chromatogram.isotope[mass][scan] > 0.0) ++isotope_scans;
        }
        if (isotope_scans < config.quant_min_scans) continue;
        isotope_totals[mass] = background_corrected_chromatographic_area(
            chromatogram.retention, chromatogram.isotope[mass],
            peak.left, peak.right);
        isotope_apices[mass] = background_corrected_chromatographic_apex(
            chromatogram.isotope[mass], peak.apex, peak.left, peak.right);
        if (isotope_totals[mass] > 0.0 ||
            isotope_apices[mass] > 0.0) {
            ++isotope_count;
            isotope_detected[mass] = true;
            isotope_peak_apex[mass] = peak.apex;
        }
    }
    if (traced < config.quant_min_scans ||
        isotope_count < config.quant_min_isotopes ||
        (isotope_totals[base_isotope] <= 0.0 &&
         isotope_apices[base_isotope] <= 0.0)) {
        return;
    }
    const double area_intensity = std::accumulate(
        isotope_totals.begin(), isotope_totals.end(), 0.0);
    const double apex_intensity = std::accumulate(
        isotope_apices.begin(), isotope_apices.end(), 0.0);
    double intensity = area_intensity;
    if (config.quant_intensity_mode == 0) {
        intensity = apex_intensity;
    } else if (config.quant_intensity_mode == 2) {
        intensity = summed_isotope_apex_intensity(isotope_apices);
    }
    if (!(intensity > 0.0) || !std::isfinite(intensity)) return;

    const double apex_rt = chromatogram.retention[apex];
    const auto nearest = std::lower_bound(
        ms1.scans.begin() + static_cast<std::ptrdiff_t>(offset),
        ms1.scans.begin() + static_cast<std::ptrdiff_t>(offset + count),
        apex_rt,
        [](const QuantScan& scan, double rt) {
            return scan.retention < rt;
        });
    std::size_t nearest_index = nearest == ms1.scans.begin() +
            static_cast<std::ptrdiff_t>(offset + count)
        ? count - 1
        : static_cast<std::size_t>(
              std::distance(
                  ms1.scans.begin() + static_cast<std::ptrdiff_t>(offset),
                  nearest));
    if (nearest_index > 0 &&
        std::abs(ms1.scans[offset + nearest_index - 1].retention - apex_rt) <
            std::abs(ms1.scans[offset + nearest_index].retention - apex_rt)) {
        --nearest_index;
    }
    const auto& apex_scan = ms1.scans[offset + nearest_index];
    psm.quantified_intensity = intensity;
    psm.apex_retention = apex_rt * 60.0;
    psm.apex_scan = apex_scan.scan;
    psm.retention_start = chromatogram.retention[left] * 60.0;
    psm.retention_end = chromatogram.retention[right - 1] * 60.0;
    psm.retention_fwhm = interpolated_half_width(
        chromatogram.retention, base_smooth, apex, left, right - 1);
    psm.traced_scans = traced;
    if (anchor_mz[nearest_index] > 0.0) {
        psm.quant_mass_error_ppm =
            (anchor_mz[nearest_index] -
             theoretical_peaks[base_isotope].mz) /
                theoretical_peaks[base_isotope].mz * 1e6;
    }
    const double isotope_sum = std::accumulate(
        isotope_totals.begin(), isotope_totals.end(), 0.0);
    const double theoretical_sum = std::accumulate(
        theoretical.begin(), theoretical.end(), 0.0);
    psm.quant_isotope_kl = 0.0;
    for (std::size_t mass = 0; mass < isotope_totals.size(); ++mass) {
        const double observed =
            std::max(1e-12, isotope_totals[mass] / isotope_sum);
        const double expected =
            std::max(1e-12, theoretical[mass] / theoretical_sum);
        psm.quant_isotope_kl += observed * std::log(observed / expected);
    }
    psm.quant_isotope_fraction = static_cast<double>(isotope_count) /
        static_cast<double>(chromatogram.isotope.size());
    double correlation_sum = 0.0;
    double apex_spread_sum = 0.0;
    std::size_t correlation_count = 0;
    std::size_t apex_count = 0;
    for (std::size_t mass = 1; mass < chromatogram.isotope.size(); ++mass) {
        if (!isotope_detected[mass]) continue;
        double base_mean = 0.0;
        double isotope_mean = 0.0;
        const double point_count = static_cast<double>(right - left);
        for (std::size_t scan = left; scan < right; ++scan) {
            base_mean += std::log1p(
                chromatogram.isotope[base_isotope][scan]);
            isotope_mean += std::log1p(
                chromatogram.isotope[mass][scan]);
        }
        base_mean /= point_count;
        isotope_mean /= point_count;
        double covariance = 0.0;
        double base_square = 0.0;
        double isotope_square = 0.0;
        for (std::size_t scan = left; scan < right; ++scan) {
            const double base_value = std::log1p(
                chromatogram.isotope[base_isotope][scan]) - base_mean;
            const double isotope_value = std::log1p(
                chromatogram.isotope[mass][scan]) - isotope_mean;
            covariance += base_value * isotope_value;
            base_square += base_value * base_value;
            isotope_square += isotope_value * isotope_value;
        }
        const double denominator = std::sqrt(base_square * isotope_square);
        if (denominator > 0.0) {
            correlation_sum += covariance / denominator;
            ++correlation_count;
        }
        apex_spread_sum += std::abs(
            chromatogram.retention[isotope_peak_apex[mass]] -
            chromatogram.retention[apex]) * 60.0;
        ++apex_count;
    }
    if (correlation_count > 0) {
        psm.quant_isotope_correlation = correlation_sum /
            static_cast<double>(correlation_count);
    }
    if (apex_count > 0) {
        psm.quant_isotope_apex_spread = apex_spread_sum /
            static_cast<double>(apex_count);
    }
    psm.has_chromatographic_feature = true;
}

std::string quant_peptide_body(const std::string& peptide) {
    const auto open = peptide.find('[');
    const auto close = peptide.rfind(']');
    if (open != std::string::npos && close != std::string::npos &&
        close > open) {
        return peptide.substr(open + 1, close - open - 1);
    }
    return peptide;
}

int quant_sip_abundance_bin(const Config& config, double abundance) {
    if (config.sip_isotope.empty()) return -1;
    const double width = config.mbr_sip_bin_width;
    const int bin_count = std::max(
        1, static_cast<int>(std::ceil(100.0 / width)));
    const double bounded = std::clamp(abundance, 0.0, 100.0);
    return std::min(
        bin_count - 1, static_cast<int>(std::floor(bounded / width)));
}

double quant_sip_bin_center(const Config& config, int bin) {
    if (config.sip_isotope.empty() || bin < 0) return 0.0;
    return std::min(
        100.0, (static_cast<double>(bin) + 0.5) *
            config.mbr_sip_bin_width);
}

std::string quant_base_ion_key(const Psm& psm) {
    return quant_peptide_body(psm.peptide) + "#" +
        std::to_string(psm.charge);
}

std::string quant_ion_key(const Psm& psm) {
    std::string key = quant_base_ion_key(psm);
    if (psm.sip_abundance_bin >= 0) {
        key += "#sipbin" + std::to_string(psm.sip_abundance_bin);
    }
    return key;
}

double quant_median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    const auto middle = values.begin() +
        static_cast<std::ptrdiff_t>(values.size() / 2);
    std::nth_element(values.begin(), middle, values.end());
    const double upper = *middle;
    if (values.size() % 2 != 0) return upper;
    const auto lower = std::max_element(values.begin(), middle);
    return 0.5 * (upper + *lower);
}

std::vector<double> quant_ranks(const std::vector<double>& values) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t left,
                                                     std::size_t right) {
        return values[left] < values[right];
    });
    std::vector<double> ranks(values.size(), 0.0);
    for (std::size_t begin = 0; begin < order.size();) {
        std::size_t end = begin + 1;
        while (end < order.size() &&
               values[order[end]] == values[order[begin]]) {
            ++end;
        }
        const double rank =
            0.5 * static_cast<double>(begin + end - 1);
        for (std::size_t index = begin; index < end; ++index) {
            ranks[order[index]] = rank;
        }
        begin = end;
    }
    return ranks;
}

double quant_spearman(
    const std::vector<double>& first,
    const std::vector<double>& second) {
    if (first.size() != second.size() || first.size() < 3) return 0.0;
    const auto first_rank = quant_ranks(first);
    const auto second_rank = quant_ranks(second);
    const double first_mean = std::accumulate(
        first_rank.begin(), first_rank.end(), 0.0) / first_rank.size();
    const double second_mean = std::accumulate(
        second_rank.begin(), second_rank.end(), 0.0) / second_rank.size();
    double numerator = 0.0;
    double first_square = 0.0;
    double second_square = 0.0;
    for (std::size_t index = 0; index < first_rank.size(); ++index) {
        const double first_delta = first_rank[index] - first_mean;
        const double second_delta = second_rank[index] - second_mean;
        numerator += first_delta * second_delta;
        first_square += first_delta * first_delta;
        second_square += second_delta * second_delta;
    }
    const double denominator = std::sqrt(first_square * second_square);
    return denominator > 0.0 ? numerator / denominator : 0.0;
}

using QuantIonMap = std::unordered_map<std::string, std::size_t>;

struct QuantAlignment {
    std::size_t donor = 0;
    double correlation = 0.0;
    std::vector<std::pair<double, double>> donor_rt_delta;
    double median_delta = 0.0;
    double mad = 0.0;
};

QuantAlignment build_quant_alignment(
    std::size_t donor, const QuantIonMap& acceptor_ions,
    const QuantIonMap& donor_ions, const Dataset& data) {
    QuantAlignment result;
    result.donor = donor;
    std::vector<double> acceptor_rt;
    std::vector<double> donor_rt;
    std::vector<double> acceptor_intensity;
    std::vector<double> donor_intensity;
    std::vector<double> deltas;
    for (const auto& [key, donor_row] : donor_ions) {
        const auto found = acceptor_ions.find(key);
        if (found == acceptor_ions.end()) continue;
        const auto& donor_psm = data.rows[donor_row];
        const auto& acceptor_psm = data.rows[found->second];
        const double donor_time = donor_psm.apex_retention / 60.0;
        const double acceptor_time = acceptor_psm.apex_retention / 60.0;
        donor_rt.push_back(donor_time);
        acceptor_rt.push_back(acceptor_time);
        donor_intensity.push_back(donor_psm.quantified_intensity);
        acceptor_intensity.push_back(acceptor_psm.quantified_intensity);
        deltas.push_back(acceptor_time - donor_time);
        result.donor_rt_delta.emplace_back(
            donor_time, acceptor_time - donor_time);
    }
    const std::size_t ion_union =
        donor_ions.size() + acceptor_ions.size() - deltas.size();
    const double overlap = ion_union == 0
        ? 0.0
        : static_cast<double>(deltas.size()) /
            static_cast<double>(ion_union);
    result.correlation = overlap * 0.5 *
        (quant_spearman(donor_rt, acceptor_rt) +
         quant_spearman(donor_intensity, acceptor_intensity));
    result.median_delta = quant_median(deltas);
    for (auto& delta : deltas) {
        delta = std::abs(delta - result.median_delta);
    }
    result.mad = quant_median(std::move(deltas));
    return result;
}

std::pair<double, double> quant_transfer_region(
    const Config& config, const QuantAlignment& alignment,
    double donor_rt) {
    std::vector<double> local;
    for (const auto& [anchor_rt, delta] : alignment.donor_rt_delta) {
        if (std::abs(anchor_rt - donor_rt) <= config.mbr_rt_window) {
            local.push_back(delta);
        }
    }
    const double shift = local.size() >= 5
        ? quant_median(local) : alignment.median_delta;
    double mad = alignment.mad;
    if (local.size() >= 5) {
        for (auto& value : local) value = std::abs(value - shift);
        mad = quant_median(std::move(local));
    }
    // Center the transfer window at donor RT plus the robust median shift and
    // span two median absolute deviations. The floor keeps a zero-MAD
    // alignment wide enough to contain one LC peak.
    const double half_window = std::clamp(
        2.0 * mad, 0.05, config.mbr_rt_window);
    return {donor_rt + shift, half_window};
}

struct QuantTransferJob {
    std::string key;
    std::size_t donor_row = 0;
    std::size_t acceptor_row = std::numeric_limits<std::size_t>::max();
    double predicted_rt = 0.0;
    double trace_rt = 0.0;
    double rt_window = 0.0;
};

using QuantTransferFeatures = std::array<double, 4>;

struct QuantTransferCandidate {
    std::string key;
    Psm psm;
    QuantTransferFeatures features{};
    double score = 0.0;
    double probability = 0.0;
    double qvalue = 1.0;
    double false_prior = 1.0;
    double decoy_score = std::numeric_limits<double>::quiet_NaN();
    std::size_t donor_row = 0;
    bool identified = false;
    bool blocked_by_direct = false;
    bool selected_bin = false;
    bool has_decoy = false;
    double predicted_rt = 0.0;
    QuantTransferFeatures decoy_features{};
};

QuantTransferFeatures quant_transfer_features(
    const Psm& psm, const Psm& donor, double predicted_rt) {
    (void)donor;
    // IonQuant's natural-abundance LDA uses exactly these four transformed
    // measurements. Aerith stores KL with natural logarithms, so convert it
    // to IonQuant's base-10 divergence and apply the same 0.01 floor before
    // the signed square-root transform.
    const double log_isotope_kl = std::log10(std::max(
        0.01, psm.quant_isotope_kl / std::log(10.0)));
    const double isotope_fit = std::copysign(
        std::sqrt(std::abs(log_isotope_kl)), log_isotope_kl);
    return {
        std::log10(std::max(1.0, psm.quantified_intensity)),
        isotope_fit,
        std::sqrt(std::abs(psm.quant_mass_error_ppm)),
        std::sqrt(std::abs(
            psm.apex_retention / 60.0 - predicted_rt))};
}

struct QuantFeatureLocation {
    double mz = 0.0;
    double apex_retention = 0.0;
    int charge = 0;
};

bool trace_quant_decoy(
    const Config& config, const QuantMs1Data& ms1,
    const Psm& donor, double predicted_rt, double rt_window,
    const std::vector<QuantFeatureLocation>& acceptor_features,
    Psm& result) {
    // IonQuant 1.11.18 tries the first traceable shifted envelope from
    // +11 through +7 Da. Keeping the same null search space is important for
    // comparable posterior calibration.
    for (int shift = 11; shift >= 7; --shift) {
        result = donor;
        result.id.clear();
        result.file_id = donor.file_id;
        result.scan = 0;
        result.parent_scan = 0;
        result.retention = predicted_rt;
        quantify_psm(
            config, ms1, result, rt_window,
            static_cast<double>(shift) * 1.0005 /
                static_cast<double>(std::max(1, donor.charge)));
        if (!result.has_chromatographic_feature) continue;
        // IonQuant does not treat a shifted envelope as null evidence when it
        // coincides with a real acceptor feature (same charge, within 0.01
        // m/z and two seconds). Its feature index includes both identified
        // and newly detected targets, so apply the same exclusion here.
        const double shifted_mz = result.calculated_mz +
            static_cast<double>(shift) * 1.0005 /
                static_cast<double>(std::max(1, donor.charge));
        const auto first = std::lower_bound(
            acceptor_features.begin(), acceptor_features.end(),
            shifted_mz - 0.01,
            [](const QuantFeatureLocation& feature, double mz) {
                return feature.mz < mz;
            });
        bool collision = false;
        for (auto feature = first;
             feature != acceptor_features.end() &&
                 feature->mz < shifted_mz + 0.01;
             ++feature) {
            if (feature->charge == donor.charge &&
                std::abs(feature->apex_retention -
                    result.apex_retention) < 2.0) {
                collision = true;
                break;
            }
        }
        if (!collision) return true;
    }
    return false;
}

struct QuantLdaModel {
    QuantTransferFeatures mean{};
    QuantTransferFeatures scale{};
    QuantTransferFeatures weight{};
    bool valid = false;

    double predict(const QuantTransferFeatures& features) const {
        double result = 0.0;
        for (std::size_t index = 0; index < features.size(); ++index) {
            result += weight[index] *
                (features[index] - mean[index]) / scale[index];
        }
        return result;
    }
};

using QuantTransferMatrix = std::array<
    std::array<double, std::tuple_size_v<QuantTransferFeatures>>,
    std::tuple_size_v<QuantTransferFeatures>>;

bool solve_quant_linear_system(
    QuantTransferMatrix matrix, QuantTransferFeatures right,
    QuantTransferFeatures& solution) {
    constexpr std::size_t count =
        std::tuple_size_v<QuantTransferFeatures>;
    for (std::size_t column = 0; column < count; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < count; ++row) {
            if (std::abs(matrix[row][column]) >
                std::abs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][column]) < 1e-12) return false;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(right[pivot], right[column]);
        }
        const double diagonal = matrix[column][column];
        for (std::size_t index = column; index < count; ++index) {
            matrix[column][index] /= diagonal;
        }
        right[column] /= diagonal;
        for (std::size_t row = 0; row < count; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            if (factor == 0.0) continue;
            for (std::size_t index = column; index < count; ++index) {
                matrix[row][index] -= factor * matrix[column][index];
            }
            right[row] -= factor * right[column];
        }
    }
    solution = right;
    return true;
}

QuantLdaModel fit_quant_lda(
    const std::vector<QuantTransferFeatures>& targets,
    const std::vector<QuantTransferFeatures>& decoys) {
    QuantLdaModel result;
    // IonQuant trains with unpaired type +2/-2 observations and only needs
    // ten positive and four negative examples. A small ridge term keeps the
    // covariance solve stable for sparse SIP-bin feature sets.
    if (targets.size() < 10 || decoys.size() < 4) return result;
    const double total =
        static_cast<double>(targets.size() + decoys.size());
    for (std::size_t feature = 0; feature < result.mean.size(); ++feature) {
        for (const auto& row : targets) result.mean[feature] += row[feature];
        for (const auto& row : decoys) result.mean[feature] += row[feature];
        result.mean[feature] /= total;
        double square = 0.0;
        for (const auto& row : targets) {
            square += std::pow(row[feature] - result.mean[feature], 2);
        }
        for (const auto& row : decoys) {
            square += std::pow(row[feature] - result.mean[feature], 2);
        }
        result.scale[feature] = std::sqrt(
            square / std::max(1.0, total - 1.0));
        result.scale[feature] = std::max(1e-6, result.scale[feature]);
    }
    QuantTransferFeatures target_mean{};
    QuantTransferFeatures decoy_mean{};
    for (const auto& row : targets) {
        for (std::size_t feature = 0; feature < row.size(); ++feature) {
            target_mean[feature] +=
                (row[feature] - result.mean[feature]) /
                result.scale[feature];
        }
    }
    for (const auto& row : decoys) {
        for (std::size_t feature = 0; feature < row.size(); ++feature) {
            decoy_mean[feature] +=
                (row[feature] - result.mean[feature]) /
                result.scale[feature];
        }
    }
    for (std::size_t feature = 0; feature < result.mean.size(); ++feature) {
        target_mean[feature] /= static_cast<double>(targets.size());
        decoy_mean[feature] /= static_cast<double>(decoys.size());
    }
    QuantTransferMatrix covariance{};
    const auto add_covariance = [&](const auto& rows,
                                    const auto& class_mean) {
        if (rows.size() < 2) return;
        const double denominator = static_cast<double>(rows.size() - 1);
        for (const auto& row : rows) {
            QuantTransferFeatures delta{};
            for (std::size_t feature = 0; feature < row.size(); ++feature) {
                delta[feature] =
                    (row[feature] - result.mean[feature]) /
                        result.scale[feature] -
                    class_mean[feature];
            }
            for (std::size_t left = 0; left < delta.size(); ++left) {
                for (std::size_t right = 0; right < delta.size(); ++right) {
                    covariance[left][right] +=
                        delta[left] * delta[right] / denominator;
                }
            }
        }
    };
    add_covariance(targets, target_mean);
    add_covariance(decoys, decoy_mean);
    double trace = 0.0;
    for (std::size_t feature = 0; feature < result.mean.size(); ++feature) {
        trace += covariance[feature][feature];
    }
    const double ridge = std::max(
        1e-4, 1e-3 * trace /
            static_cast<double>(result.mean.size()));
    QuantTransferFeatures difference{};
    for (std::size_t feature = 0; feature < result.mean.size(); ++feature) {
        covariance[feature][feature] += ridge;
        difference[feature] = target_mean[feature] - decoy_mean[feature];
    }
    if (!solve_quant_linear_system(
            covariance, difference, result.weight)) {
        return result;
    }
    double length = 0.0;
    for (const double weight : result.weight) length += weight * weight;
    length = std::sqrt(length);
    if (!(length > 0.0) || !std::isfinite(length)) return result;
    for (auto& weight : result.weight) weight /= length;
    result.valid = true;
    return result;
}

std::pair<double, double> quant_mean_sd(
    const std::vector<double>& values) {
    if (values.empty()) return {0.0, 0.0};
    const double mean = std::accumulate(
        values.begin(), values.end(), 0.0) /
        static_cast<double>(values.size());
    double square = 0.0;
    for (const double value : values) {
        const double delta = value - mean;
        square += delta * delta;
    }
    return {
        mean,
        std::sqrt(square /
            static_cast<double>(std::max<std::size_t>(1, values.size() - 1)))};
}

std::vector<double> mbr_posterior_probabilities(
    const std::vector<double>& target_scores,
    const std::vector<double>& decoy_scores,
    const std::vector<double>& identified_target_scores,
    const std::vector<double>& identified_decoy_scores,
    double* false_prior_result) {
    (void)identified_decoy_scores;
    // IonQuant fits this model only when +1, -1, and +2 each contain more
    // than ten observations. Its -2 scores are passed to the model but are
    // not used by the one-dimensional mixture implementation.
    if (target_scores.size() <= 10 || decoy_scores.size() <= 10 ||
        identified_target_scores.size() <= 10) {
        return {};
    }
    constexpr std::size_t grid_size = 1000;
    constexpr double sqrt_two_pi = 2.5066282746310005024;
    const auto [target_mean, target_sd] = quant_mean_sd(target_scores);
    const auto [decoy_mean, decoy_sd] = quant_mean_sd(decoy_scores);
    std::vector<double> sorted_targets = target_scores;
    std::vector<double> sorted_decoys = decoy_scores;
    std::vector<double> sorted_identified = identified_target_scores;
    std::sort(sorted_targets.begin(), sorted_targets.end());
    std::sort(sorted_decoys.begin(), sorted_decoys.end());
    std::sort(sorted_identified.begin(), sorted_identified.end());
    const auto median = [](const std::vector<double>& values) {
        const auto middle = values.size() / 2;
        return values.size() % 2 == 0
            ? 0.5 * (values[middle - 1] + values[middle])
            : values[middle];
    };
    // SSJ's EmpiricalDist defines the quartiles as the medians of the two
    // half samples (excluding the overall median for odd-sized samples).
    const std::size_t half = sorted_targets.size() / 2;
    const auto half_median = [&](std::size_t begin) {
        const auto middle = half / 2;
        return half % 2 == 0
            ? 0.5 * (sorted_targets[begin + middle - 1] +
                     sorted_targets[begin + middle])
            : sorted_targets[begin + middle];
    };
    const double first_quartile = half_median(0);
    const double third_quartile =
        half_median(sorted_targets.size() - half);
    double bandwidth = 0.99 * std::min(
        target_sd, (third_quartile - first_quartile) / 1.34) /
        std::pow(static_cast<double>(target_scores.size()), 0.2);
    // This asymmetric grid is intentional: IonQuant uses the lower extent of
    // +1/-1 and the upper extent of +1/+2.
    const float minimum_float = static_cast<float>(std::min(
        sorted_targets.front(), sorted_decoys.front()));
    const float maximum_float = static_cast<float>(std::max(
        sorted_targets.back(), sorted_identified.back()));
    const double minimum = minimum_float;
    const double maximum = maximum_float;
    const double score_range = maximum - minimum;
    if (!(score_range > 1e-9) || !(bandwidth > 0.0) ||
        !std::isfinite(bandwidth)) return {};
    const float step_float = static_cast<float>(
        score_range / static_cast<double>(grid_size));
    std::vector<double> grid(grid_size);
    for (std::size_t index = 0; index < grid_size; ++index) {
        grid[index] = static_cast<double>(static_cast<float>(
            minimum_float + static_cast<float>(index) * step_float));
    }
    const auto floor_bin = [&](double score) {
        const auto found = std::upper_bound(grid.begin(), grid.end(), score);
        return found == grid.begin() ? std::size_t{0} :
            static_cast<std::size_t>(found - grid.begin() - 1);
    };
    std::vector<std::size_t> target_bins(target_scores.size());
    for (std::size_t index = 0; index < target_scores.size(); ++index) {
        target_bins[index] = floor_bin(target_scores[index]);
    }
    const double kernel_scale = 1.0 / (bandwidth * sqrt_two_pi);
    std::vector<double> false_density(grid_size, 0.0);
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t grid_index = 0;
         grid_index < static_cast<std::ptrdiff_t>(grid_size); ++grid_index) {
        double sum = 0.0;
        for (const double score : decoy_scores) {
            const double z = (grid[static_cast<std::size_t>(grid_index)] -
                score) / bandwidth;
            sum += std::exp(-0.5 * z * z);
        }
        false_density[static_cast<std::size_t>(grid_index)] =
            sum * kernel_scale / static_cast<double>(decoy_scores.size());
    }
    // Preserve each target's exact location in the weighted KDE. Binning
    // these observations materially changes IonQuant's fitted prior.
    std::vector<double> target_kernel(
        grid_size * target_scores.size());
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t grid_index = 0;
         grid_index < static_cast<std::ptrdiff_t>(grid_size); ++grid_index) {
        const auto offset = static_cast<std::size_t>(grid_index) *
            target_scores.size();
        for (std::size_t target = 0; target < target_scores.size(); ++target) {
            const double z = (grid[static_cast<std::size_t>(grid_index)] -
                target_scores[target]) / bandwidth;
            target_kernel[offset + target] = std::exp(-0.5 * z * z);
        }
    }
    const double lower_threshold = std::max(
        target_mean - 2.5 * target_sd,
        decoy_mean - decoy_sd);
    const auto target_lower = static_cast<double>(std::count_if(
        target_scores.begin(), target_scores.end(),
        [&](double value) { return value <= lower_threshold; })) /
        static_cast<double>(target_scores.size());
    const auto decoy_lower = static_cast<double>(std::count_if(
        decoy_scores.begin(), decoy_scores.end(),
        [&](double value) { return value <= lower_threshold; })) /
        static_cast<double>(decoy_scores.size());
    const double initial_false_prior = decoy_lower > 0.0
        ? target_lower / decoy_lower : 0.0;

    // IonQuant seeds a two-cluster fit of the +1 scores at the medians of
    // -1 and +2. The target-like cluster supplies the initial positive
    // density before the weighted-KDE refinement.
    double null_center = median(sorted_decoys);
    const double identified_center = median(sorted_identified);
    double true_center = identified_center;
    std::vector<unsigned char> cluster(target_scores.size(), 2);
    for (unsigned int iteration = 0; iteration < 30; ++iteration) {
        double sums[2]{};
        std::size_t counts[2]{};
        bool changed = false;
        for (std::size_t index = 0; index < target_scores.size(); ++index) {
            const double score = target_scores[index];
            const unsigned char label =
                std::abs(score - true_center) <
                    std::abs(score - null_center) ? 1 : 0;
            changed |= cluster[index] != label;
            cluster[index] = label;
            sums[label] += score;
            ++counts[label];
        }
        const double next_null = counts[0] > 0
            ? sums[0] / static_cast<double>(counts[0]) : null_center;
        const double next_true = counts[1] > 0
            ? sums[1] / static_cast<double>(counts[1]) : true_center;
        null_center = next_null;
        true_center = next_true;
        if (!changed) break;
    }
    double component_mean[2]{null_center, true_center};
    double component_variance[2]{};
    double component_weight[2]{};
    std::size_t component_count[2]{};
    for (std::size_t index = 0; index < target_scores.size(); ++index) {
        const auto label = cluster[index];
        const double delta = target_scores[index] - component_mean[label];
        component_variance[label] += delta * delta;
        ++component_count[label];
    }
    for (std::size_t component = 0; component < 2; ++component) {
        if (component_count[component] == 0) return {};
        component_variance[component] /=
            static_cast<double>(component_count[component]);
        component_weight[component] =
            static_cast<double>(component_count[component]) /
            static_cast<double>(target_scores.size());
    }
    // IonQuant follows k-means with a two-normal EM fit and deliberately uses
    // component 1 (the component seeded at the +2 median) as the initial true
    // distribution without multiplying by its fitted mixture weight.
    std::vector<std::array<double, 2>> responsibilities(
        target_scores.size());
    double initial_likelihood = 0.0;
    for (const double score : target_scores) {
        double mixture = 0.0;
        for (std::size_t component = 0; component < 2; ++component) {
            const double z = (score - component_mean[component]) /
                std::sqrt(component_variance[component]);
            mixture += component_weight[component] *
                std::exp(-0.5 * z * z) /
                (std::sqrt(component_variance[component]) * sqrt_two_pi);
        }
        initial_likelihood += std::log(std::max(1e-300, mixture));
    }
    double previous_gaussian_likelihood = initial_likelihood;
    const double gaussian_tolerance =
        std::abs(initial_likelihood) * 1e-4;
    for (unsigned int iteration = 0; iteration < 30; ++iteration) {
        for (std::size_t index = 0; index < target_scores.size(); ++index) {
            double total_density = 0.0;
            for (std::size_t component = 0; component < 2; ++component) {
                const double z =
                    (target_scores[index] - component_mean[component]) /
                    std::sqrt(component_variance[component]);
                responsibilities[index][component] =
                    component_weight[component] *
                    std::exp(-0.5 * z * z) /
                    (std::sqrt(component_variance[component]) * sqrt_two_pi);
                total_density += responsibilities[index][component];
            }
            for (std::size_t component = 0; component < 2; ++component) {
                responsibilities[index][component] /=
                    std::max(1e-300, total_density);
            }
        }
        for (std::size_t component = 0; component < 2; ++component) {
            double responsibility_sum = 0.0;
            double weighted_sum = 0.0;
            for (std::size_t index = 0; index < target_scores.size(); ++index) {
                responsibility_sum += responsibilities[index][component];
                weighted_sum += responsibilities[index][component] *
                    target_scores[index];
            }
            if (!(responsibility_sum > 0.0)) return {};
            component_mean[component] = weighted_sum / responsibility_sum;
            double weighted_square = 0.0;
            for (std::size_t index = 0; index < target_scores.size(); ++index) {
                const double delta =
                    target_scores[index] - component_mean[component];
                weighted_square += responsibilities[index][component] *
                    delta * delta;
            }
            component_variance[component] =
                weighted_square / responsibility_sum;
            component_weight[component] = responsibility_sum /
                static_cast<double>(target_scores.size());
        }
        double likelihood = 0.0;
        for (const double score : target_scores) {
            double mixture = 0.0;
            for (std::size_t component = 0; component < 2; ++component) {
                const double z =
                    (score - component_mean[component]) /
                    std::sqrt(component_variance[component]);
                mixture += component_weight[component] *
                    std::exp(-0.5 * z * z) /
                    (std::sqrt(component_variance[component]) * sqrt_two_pi);
            }
            likelihood += std::log(std::max(1e-300, mixture));
        }
        if (std::abs(likelihood - previous_gaussian_likelihood) <=
            gaussian_tolerance) {
            break;
        }
        previous_gaussian_likelihood = likelihood;
    }
    const double initial_true_mean = component_mean[1];
    const double initial_true_sd = std::sqrt(component_variance[1]);
    std::vector<double> weighted_false(grid_size, 0.0);
    std::vector<double> weighted_true(grid_size, 0.0);
    for (std::size_t index = 0; index < grid_size; ++index) {
        weighted_false[index] =
            initial_false_prior * false_density[index];
        const double score = grid[index];
        const double z = (score - initial_true_mean) / initial_true_sd;
        weighted_true[index] = (1.0 - initial_false_prior) *
            std::exp(-0.5 * z * z) /
            (initial_true_sd * sqrt_two_pi);
    }
    std::vector<double> posterior_true(target_scores.size(), 0.0);
    const auto update_posteriors = [&]() {
        double true_sum = 0.0;
        for (std::size_t index = 0; index < target_scores.size(); ++index) {
            const auto bin = target_bins[index];
            const double null_value = weighted_false[bin];
            const double true_value = weighted_true[bin];
            const double denominator = null_value + true_value;
            posterior_true[index] = denominator > 0.0
                ? std::clamp(true_value / denominator, 0.0, 1.0)
                : 0.0;
            true_sum += posterior_true[index];
        }
        return true_sum;
    };
    double true_sum = update_posteriors();
    const auto rebuild_true_density = [&](double true_prior) {
        std::vector<double> result(grid_size, 0.0);
        const double weight_sum = std::accumulate(
            posterior_true.begin(), posterior_true.end(), 0.0);
        if (!(weight_sum > 0.0)) return result;
        #pragma omp parallel for schedule(static)
        for (std::ptrdiff_t grid_index = 0;
             grid_index < static_cast<std::ptrdiff_t>(grid_size);
             ++grid_index) {
            const auto offset = static_cast<std::size_t>(grid_index) *
                target_scores.size();
            double sum = 0.0;
            for (std::size_t target = 0;
                 target < target_scores.size(); ++target) {
                sum += target_kernel[offset + target] *
                    posterior_true[target];
            }
            result[static_cast<std::size_t>(grid_index)] =
                true_prior * sum * kernel_scale / weight_sum;
        }
        return result;
    };
    weighted_true = rebuild_true_density(1.0 - initial_false_prior);
    double previous_likelihood = 0.0;
    for (const auto bin : target_bins) {
        const double mixture = weighted_false[bin] + weighted_true[bin];
        if (mixture > 0.0) previous_likelihood += std::log(mixture);
    }
    const double likelihood_tolerance =
        std::abs(previous_likelihood) * 1e-5;
    for (unsigned int iteration = 0; iteration < 50; ++iteration) {
        true_sum = update_posteriors();
        weighted_true = rebuild_true_density(1.0 - initial_false_prior);
        double likelihood = 0.0;
        for (const auto bin : target_bins) {
            const double mixture = weighted_false[bin] + weighted_true[bin];
            if (mixture > 0.0) likelihood += std::log(mixture);
        }
        if (std::abs(likelihood - previous_likelihood) <=
            likelihood_tolerance) {
            break;
        }
        previous_likelihood = likelihood;
    }
    const double final_true_prior =
        true_sum / static_cast<double>(target_scores.size());
    const double final_false_prior = 1.0 - final_true_prior;
    // IonQuant recomputes the two final weighted densities using the last
    // posterior weights, then constructs a descending lookup table.
    update_posteriors();
    weighted_true = rebuild_true_density(final_true_prior);
    for (std::size_t index = 0; index < grid_size; ++index) {
        weighted_false[index] = final_false_prior * false_density[index];
    }
    if (false_prior_result != nullptr) {
        *false_prior_result = final_false_prior;
    }
    std::vector<double> grid_probability(grid_size, 0.0);
    for (std::size_t index = 0; index < grid_size; ++index) {
        const double denominator = weighted_true[index] +
            weighted_false[index];
        grid_probability[index] = denominator == 0.0 ? 0.0 :
            0.999999 * weighted_true[index] / denominator;
    }
    const auto above_decoy_mean = std::upper_bound(
        grid.begin(), grid.end(), decoy_mean);
    if (above_decoy_mean != grid.end()) {
        std::size_t index = static_cast<std::size_t>(
            above_decoy_mean - grid.begin());
        while (index > 0) {
            grid_probability[index - 1] = std::min(
                grid_probability[index], grid_probability[index - 1]);
            --index;
        }
    }
    std::vector<double> result(target_scores.size());
    for (std::size_t index = 0; index < target_scores.size(); ++index) {
        const float score = static_cast<float>(target_scores[index]);
        const auto found = std::lower_bound(
            grid.begin(), grid.end(), static_cast<double>(score));
        result[index] = found == grid.end() ? 1.0 :
            grid_probability[static_cast<std::size_t>(
                found - grid.begin())];
    }
    return result;
}

bool calibrate_quant_transfer_probabilities(
    std::vector<QuantTransferCandidate*>& targets,
    const std::vector<double>& decoy_scores,
    const std::vector<double>& identified_target_scores,
    const std::vector<double>& identified_decoy_scores,
    double* false_prior) {
    std::vector<double> target_scores;
    target_scores.reserve(targets.size());
    for (const auto* target : targets) target_scores.push_back(target->score);
    const auto probabilities = mbr_posterior_probabilities(
        target_scores, decoy_scores, identified_target_scores,
        identified_decoy_scores, false_prior);
    if (probabilities.size() != targets.size()) return false;
    for (std::size_t index = 0; index < targets.size(); ++index) {
        targets[index]->probability = probabilities[index];
    }
    return true;
}

void assign_quant_transfer_qvalues(
    std::vector<QuantTransferCandidate*>& targets) {
    std::sort(targets.begin(), targets.end(), [](const auto* left,
                                                 const auto* right) {
        return left->probability > right->probability;
    });
    double false_sum = 0.0;
    std::vector<double> raw_fdr(targets.size(), 1.0);
    for (std::size_t index = 0; index < targets.size(); ++index) {
        false_sum += 1.0 - targets[index]->probability;
        raw_fdr[index] = false_sum / static_cast<double>(index + 1);
    }
    double minimum = 1.0;
    for (std::size_t index = targets.size(); index-- > 0;) {
        minimum = std::min(minimum, raw_fdr[index]);
        targets[index]->qvalue = minimum;
    }
}

QuantificationResult ChromatographicQuantifier::add(
    const Config& config, Dataset& data, const std::vector<double>& q) {
    using Clock = std::chrono::steady_clock;
    QuantificationResult result;
    if (config.spectrum_paths.empty()) return result;
    if (config.spectrum_paths.size() != data.input_paths.size() ||
        q.size() != data.rows.size()) {
        throw std::runtime_error(
            "Chromatographic quantification received inconsistent inputs");
    }
    const auto elapsed_seconds = [](Clock::time_point begin,
                                    Clock::time_point end) {
        return std::chrono::duration<double>(end - begin).count();
    };
    const auto cpu_seconds = [](std::clock_t begin, std::clock_t end) {
        return static_cast<double>(end - begin) / CLOCKS_PER_SEC;
    };
    const auto add_elapsed = [&](StageTiming& timing,
                                 Clock::time_point wall_begin,
                                 std::clock_t cpu_begin) {
        timing.wall_seconds += elapsed_seconds(wall_begin, Clock::now());
        timing.cpu_seconds += cpu_seconds(cpu_begin, std::clock());
    };
    const auto record_stage = [&](std::string name, const StageTiming& timing,
                                  bool uses_omp = false,
                                  bool uses_simd = false) {
        result.stages.push_back({
            std::move(name), timing, uses_omp, uses_simd});
    };
    const std::size_t file_count = config.spectrum_paths.size();
    std::vector<QuantIonMap> ions(file_count);
    std::vector<std::unordered_set<std::string>> identified(file_count);
    std::vector<std::unordered_set<std::string>> identified_base(file_count);
    StageTiming load_identified_timing;
    StageTiming select_timing;
    StageTiming trace_identified_timing;
    std::size_t accepted_psms = 0;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        auto& psm = data.rows[row];
        psm.sip_abundance_bin = quant_sip_abundance_bin(
            config, psm.ms2_isotopic_abundance);
    }
    for (std::size_t file = 0; file < file_count; ++file) {
        auto wall_begin = Clock::now();
        auto cpu_begin = std::clock();
        const auto ms1 = load_quant_ms1(config.spectrum_paths[file]);
        add_elapsed(load_identified_timing, wall_begin, cpu_begin);

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        std::vector<std::size_t> rows;
        for (std::size_t row = 0; row < data.rows.size(); ++row) {
            if (data.rows[row].file_id == file &&
                data.rows[row].label == 1 &&
                q[row] <= config.q_threshold) {
                rows.push_back(row);
            }
        }
        accepted_psms += rows.size();
        add_elapsed(select_timing, wall_begin, cpu_begin);

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        #pragma omp parallel for schedule(dynamic, 32)
        for (std::ptrdiff_t index = 0;
             index < static_cast<std::ptrdiff_t>(rows.size()); ++index) {
            quantify_psm(config, ms1, data.rows[rows[
                static_cast<std::size_t>(index)]]);
        }
        add_elapsed(trace_identified_timing, wall_begin, cpu_begin);
    }

    record_stage(
        "Load MS1 data for identified runs (" +
            std::to_string(file_count) + " files)",
        load_identified_timing);
    record_stage(
        "Select accepted PSMs (" +
            std::to_string(accepted_psms) + ")",
        select_timing);
    record_stage(
        "Trace identified XICs + detect peaks/intensity",
        trace_identified_timing, true);

    const auto index_wall_begin = Clock::now();
    const auto index_cpu_begin = std::clock();
    std::size_t chromatographic_features = 0;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto& psm = data.rows[row];
        if (psm.label != 1 || q[row] > config.q_threshold ||
            psm.file_id >= file_count) {
            continue;
        }
        const auto key = quant_ion_key(psm);
        identified[psm.file_id].insert(key);
        identified_base[psm.file_id].insert(
            quant_base_ion_key(psm));
        if (!psm.has_chromatographic_feature) continue;
        ++chromatographic_features;
        auto [found, inserted] = ions[psm.file_id].emplace(key, row);
        if (!inserted &&
            psm.quantified_intensity >
                data.rows[found->second].quantified_intensity) {
            found->second = row;
        }
    }
    StageTiming index_timing;
    add_elapsed(index_timing, index_wall_begin, index_cpu_begin);
    std::size_t indexed_ions = 0;
    for (const auto& sample : ions) indexed_ions += sample.size();
    record_stage(
        "Index quantified evidence (" +
            std::to_string(chromatographic_features) + " PSMs, " +
            std::to_string(indexed_ions) + " ions)",
        index_timing);
    if (!config.mbr || file_count < 2) return result;

    data.transferred_ions.clear();
    StageTiming alignment_timing;
    StageTiming scheduling_timing;
    StageTiming load_mbr_timing;
    StageTiming trace_mbr_timing;
    StageTiming score_mbr_timing;
    std::size_t alignment_count = 0;
    std::size_t transfer_jobs = 0;
    std::size_t scored_transfer_candidates = 0;
    std::size_t calibrated_transfer_candidates = 0;
    std::size_t training_target_count = 0;
    std::size_t training_decoy_count = 0;
    std::size_t transfer_decoy_count = 0;
    std::vector<double> calibration_false_priors(
        file_count, std::numeric_limits<double>::quiet_NaN());
    std::vector<std::size_t> accepted_transfers(file_count, 0);
    std::vector<QuantTransferFeatures> calibration_lda_weights(file_count);
    std::vector<bool> calibration_lda_valid(file_count, false);
    struct QuantPooledEvidence {
        std::vector<std::size_t> target_indices;
    };
    std::vector<QuantTransferCandidate> pooled_candidates;
    pooled_candidates.reserve(indexed_ions);
    std::map<int, QuantPooledEvidence> pooled_evidence;
    for (std::size_t acceptor = 0; acceptor < file_count; ++acceptor) {
        auto wall_begin = Clock::now();
        auto cpu_begin = std::clock();
        std::vector<QuantAlignment> alignments;
        for (std::size_t donor = 0; donor < file_count; ++donor) {
            if (donor == acceptor) continue;
            auto alignment = build_quant_alignment(
                donor, ions[acceptor], ions[donor], data);
            if (alignment.correlation > config.mbr_min_correlation &&
                alignment.donor_rt_delta.size() >= 5) {
                alignments.push_back(std::move(alignment));
            }
        }
        std::sort(
            alignments.begin(), alignments.end(),
            [](const auto& left, const auto& right) {
                return left.correlation > right.correlation;
            });
        if (alignments.size() > config.mbr_top_runs) {
            alignments.resize(config.mbr_top_runs);
        }
        alignment_count += alignments.size();
        add_elapsed(alignment_timing, wall_begin, cpu_begin);
        if (alignments.empty()) continue;

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        std::vector<QuantTransferJob> jobs;
        // Trace every eligible donor. A peptide can align differently from
        // each run, so choose its donor only after scoring acceptor evidence.
        for (const auto& alignment : alignments) {
            for (const auto& [key, donor_row] : ions[alignment.donor]) {
                const auto& donor_psm = data.rows[donor_row];
                const auto region = quant_transfer_region(
                    config, alignment,
                    donor_psm.apex_retention / 60.0);
                const auto acceptor_ion = ions[acceptor].find(key);
                if (acceptor_ion != ions[acceptor].end()) {
                    QuantTransferJob job;
                    job.key = key;
                    job.donor_row = donor_row;
                    job.acceptor_row = acceptor_ion->second;
                    job.predicted_rt = region.first;
                    job.trace_rt = region.first;
                    job.rt_window = region.second;
                    jobs.push_back(std::move(job));
                    continue;
                }
                // IonQuant enumerates chromatographic features throughout
                // the aligned interval before LDA selects the best evidence.
                // Probe five overlapping local basins to reproduce that
                // behavior without retracing every MS1 point independently.
                static constexpr double offsets[] = {
                    -0.8, -0.4, 0.0, 0.4, 0.8};
                for (const double offset : offsets) {
                    QuantTransferJob job;
                    job.key = key;
                    job.donor_row = donor_row;
                    job.predicted_rt = region.first;
                    job.trace_rt = region.first + offset * region.second;
                    job.rt_window = std::max(0.025, region.second * 0.5);
                    jobs.push_back(std::move(job));
                }
            }
        }
        transfer_jobs += jobs.size();
        add_elapsed(scheduling_timing, wall_begin, cpu_begin);

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        const auto ms1 = load_quant_ms1(config.spectrum_paths[acceptor]);
        add_elapsed(load_mbr_timing, wall_begin, cpu_begin);

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        std::vector<QuantTransferCandidate> candidates(jobs.size());
        #pragma omp parallel for schedule(dynamic, 8)
        for (std::ptrdiff_t index = 0;
             index < static_cast<std::ptrdiff_t>(jobs.size()); ++index) {
            const auto& job = jobs[static_cast<std::size_t>(index)];
            auto& candidate = candidates[static_cast<std::size_t>(index)];
            candidate.key = job.key;
            candidate.donor_row = job.donor_row;
            candidate.predicted_rt = job.predicted_rt;
            candidate.identified =
                identified[acceptor].count(job.key) != 0;
            const auto& donor_psm = data.rows[job.donor_row];
            candidate.blocked_by_direct =
                identified_base[acceptor].count(
                    quant_base_ion_key(donor_psm)) != 0;
            if (job.acceptor_row !=
                std::numeric_limits<std::size_t>::max()) {
                candidate.psm = data.rows[job.acceptor_row];
            } else {
                candidate.psm = donor_psm;
                candidate.psm.file_id = acceptor;
                candidate.psm.id.clear();
                candidate.psm.scan = 0;
                candidate.psm.parent_scan = 0;
                candidate.psm.retention = job.trace_rt;
                quantify_psm(
                    config, ms1, candidate.psm, job.rt_window);
                candidate.psm.ms1_isotopic_abundance =
                    quant_sip_bin_center(
                        config, candidate.psm.sip_abundance_bin);
            }
            if (candidate.psm.has_chromatographic_feature) {
                candidate.features = quant_transfer_features(
                    candidate.psm, donor_psm, job.predicted_rt);
            }
        }
        std::vector<QuantFeatureLocation> acceptor_features;
        acceptor_features.reserve(
            ions[acceptor].size() + candidates.size());
        for (const auto& [key, row] : ions[acceptor]) {
            (void)key;
            const auto& feature = data.rows[row];
            acceptor_features.push_back({
                feature.calculated_mz, feature.apex_retention,
                feature.charge});
        }
        for (const auto& candidate : candidates) {
            if (!candidate.psm.has_chromatographic_feature) continue;
            acceptor_features.push_back({
                candidate.psm.calculated_mz,
                candidate.psm.apex_retention,
                candidate.psm.charge});
        }
        std::sort(
            acceptor_features.begin(), acceptor_features.end(),
            [](const auto& left, const auto& right) {
                return left.mz < right.mz;
            });
        #pragma omp parallel for schedule(dynamic, 8)
        for (std::ptrdiff_t index = 0;
             index < static_cast<std::ptrdiff_t>(jobs.size()); ++index) {
            const auto& job = jobs[static_cast<std::size_t>(index)];
            auto& candidate = candidates[static_cast<std::size_t>(index)];
            if (!candidate.psm.has_chromatographic_feature) continue;
            const auto& donor_psm = data.rows[job.donor_row];
            Psm decoy;
            const double decoy_rt =
                candidate.psm.apex_retention / 60.0;
            const double decoy_window = std::clamp(
                candidate.psm.retention_fwhm / 120.0,
                0.025, job.rt_window);
            candidate.has_decoy = trace_quant_decoy(
                config, ms1, candidate.psm, decoy_rt,
                decoy_window, acceptor_features, decoy);
            if (candidate.has_decoy) {
                candidate.decoy_features = quant_transfer_features(
                    decoy, donor_psm, job.predicted_rt);
            }
        }
        add_elapsed(trace_mbr_timing, wall_begin, cpu_begin);

        wall_begin = Clock::now();
        cpu_begin = std::clock();
        struct QuantBinEvidence {
            std::vector<QuantTransferFeatures> training_targets;
            std::vector<QuantTransferFeatures> training_decoys;
        };
        std::map<int, QuantBinEvidence> bin_evidence;
        for (const auto& candidate : candidates) {
            if (!candidate.identified) continue;
            auto& bin = bin_evidence[
                candidate.psm.sip_abundance_bin];
            if (candidate.psm.has_chromatographic_feature) {
                bin.training_targets.push_back(candidate.features);
            }
            if (candidate.has_decoy) {
                bin.training_decoys.push_back(candidate.decoy_features);
            }
        }
        std::map<int, QuantLdaModel> initial_models;
        for (const auto& [sip_bin, evidence] : bin_evidence) {
            initial_models.emplace(
                sip_bin,
                fit_quant_lda(
                    evidence.training_targets, evidence.training_decoys));
        }
        // IonQuant scores all donor observations once, retains the best
        // positive and negative evidence at an acceptor ion, and refits the
        // final LDA. This prevents ions present in several donors from being
        // over-weighted in training.
        std::map<int, std::unordered_map<std::string, std::size_t>>
            best_training_targets;
        std::map<int, std::unordered_map<std::string, std::size_t>>
            best_training_decoys;
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            const auto& candidate = candidates[index];
            if (!candidate.identified) continue;
            const int sip_bin = candidate.psm.sip_abundance_bin;
            const auto model = initial_models.find(sip_bin);
            if (model == initial_models.end() || !model->second.valid) {
                continue;
            }
            if (candidate.psm.has_chromatographic_feature) {
                const double score = model->second.predict(
                    candidate.features);
                auto& best = best_training_targets[sip_bin];
                const auto found = best.find(candidate.key);
                if (found == best.end() || score > model->second.predict(
                        candidates[found->second].features)) {
                    best[candidate.key] = index;
                }
            }
            if (candidate.has_decoy) {
                const double score = model->second.predict(
                    candidate.decoy_features);
                auto& best = best_training_decoys[sip_bin];
                const auto found = best.find(candidate.key);
                if (found == best.end() || score > model->second.predict(
                        candidates[found->second].decoy_features)) {
                    best[candidate.key] = index;
                }
            }
        }
        std::map<int, QuantLdaModel> models = initial_models;
        for (const auto& [sip_bin, target_rows] :
             best_training_targets) {
            QuantBinEvidence selected;
            selected.training_targets.reserve(target_rows.size());
            for (const auto& [key, index] : target_rows) {
                (void)key;
                selected.training_targets.push_back(
                    candidates[index].features);
            }
            const auto decoy_rows = best_training_decoys.find(sip_bin);
            if (decoy_rows != best_training_decoys.end()) {
                selected.training_decoys.reserve(decoy_rows->second.size());
                for (const auto& [key, index] : decoy_rows->second) {
                    (void)key;
                    selected.training_decoys.push_back(
                        candidates[index].decoy_features);
                }
            }
            training_target_count += selected.training_targets.size();
            training_decoy_count += selected.training_decoys.size();
            auto final_model = fit_quant_lda(
                selected.training_targets, selected.training_decoys);
            if (final_model.valid) models[sip_bin] = final_model;
        }
        const auto natural_model = models.find(
            config.sip_isotope.empty() ? -1 : 0);
        if (natural_model != models.end() && natural_model->second.valid) {
            calibration_lda_weights[acceptor] = natural_model->second.weight;
            calibration_lda_valid[acceptor] = true;
        }
        std::map<int, std::vector<double>> identified_target_scores;
        std::map<int, std::vector<double>> identified_decoy_scores;
        for (const auto& [sip_bin, target_rows] : best_training_targets) {
            const auto model = models.find(sip_bin);
            if (model == models.end() || !model->second.valid) continue;
            auto& scores = identified_target_scores[sip_bin];
            scores.reserve(target_rows.size());
            for (const auto& [key, index] : target_rows) {
                (void)key;
                scores.push_back(model->second.predict(
                    candidates[index].features));
            }
        }
        for (const auto& [sip_bin, decoy_rows] : best_training_decoys) {
            const auto model = models.find(sip_bin);
            if (model == models.end() || !model->second.valid) continue;
            auto& scores = identified_decoy_scores[sip_bin];
            scores.reserve(decoy_rows.size());
            for (const auto& [key, index] : decoy_rows) {
                (void)key;
                scores.push_back(model->second.predict(
                    candidates[index].decoy_features));
            }
        }
        std::unordered_map<std::string, std::size_t> best_donors;
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            auto& candidate = candidates[index];
            if (candidate.identified || candidate.blocked_by_direct ||
                !candidate.psm.has_chromatographic_feature) {
                continue;
            }
            const auto model = models.find(candidate.psm.sip_abundance_bin);
            if (model == models.end() || !model->second.valid) continue;
            candidate.score = model->second.predict(candidate.features);
            const auto ion = quant_ion_key(candidate.psm);
            const auto found = best_donors.find(ion);
            if (found == best_donors.end()) {
                best_donors.emplace(ion, index);
                continue;
            }
            const auto& current = candidates[found->second];
            if (candidate.score > current.score) found->second = index;
        }
        std::unordered_map<std::string, std::size_t> selected_bins;
        for (const auto& [ion, index] : best_donors) {
            (void)ion;
            const auto& candidate = candidates[index];
            const auto base = quant_base_ion_key(candidate.psm);
            const auto found = selected_bins.find(base);
            if (found == selected_bins.end()) {
                selected_bins.emplace(base, index);
                continue;
            }
            const auto& current = candidates[found->second];
            const auto quality = std::tuple{
                candidate.psm.quant_isotope_kl,
                std::abs(candidate.psm.apex_retention / 60.0 -
                    candidate.predicted_rt),
                -candidate.psm.quantified_intensity};
            const auto current_quality = std::tuple{
                current.psm.quant_isotope_kl,
                std::abs(current.psm.apex_retention / 60.0 -
                    current.predicted_rt),
                -current.psm.quantified_intensity};
            if (quality < current_quality) found->second = index;
        }
        for (const auto& [base, index] : selected_bins) {
            (void)base;
            candidates[index].selected_bin = true;
        }
        std::map<int, std::vector<QuantTransferCandidate*>>
            transfer_targets;
        std::map<int, std::vector<double>> transfer_decoys;
        for (auto& candidate : candidates) {
            if (!candidate.selected_bin) continue;
            const auto model = models.find(candidate.psm.sip_abundance_bin);
            if (model == models.end() || !model->second.valid) continue;
            transfer_targets[candidate.psm.sip_abundance_bin]
                .push_back(&candidate);
            ++scored_transfer_candidates;
            if (candidate.has_decoy) {
                candidate.decoy_score = model->second.predict(
                    candidate.decoy_features);
                transfer_decoys[candidate.psm.sip_abundance_bin].push_back(
                    candidate.decoy_score);
                ++transfer_decoy_count;
            }
        }
        double false_prior_sum = 0.0;
        std::size_t false_prior_weight = 0;
        for (auto& [sip_bin, targets] : transfer_targets) {
            const auto found = transfer_decoys.find(sip_bin);
            const std::vector<double> empty;
            const auto& decoys = found == transfer_decoys.end()
                ? empty : found->second;
            const auto found_identified_targets =
                identified_target_scores.find(sip_bin);
            const auto& training_targets =
                found_identified_targets == identified_target_scores.end()
                ? empty : found_identified_targets->second;
            const auto found_identified_decoys =
                identified_decoy_scores.find(sip_bin);
            const auto& training_decoys =
                found_identified_decoys == identified_decoy_scores.end()
                ? empty : found_identified_decoys->second;
            double false_prior = 1.0;
            if (!calibrate_quant_transfer_probabilities(
                    targets, decoys, training_targets, training_decoys,
                    &false_prior)) {
                continue;
            }
            false_prior_sum += false_prior *
                static_cast<double>(targets.size());
            false_prior_weight += targets.size();
            calibrated_transfer_candidates += targets.size();
            for (auto* candidate : targets) {
                candidate->false_prior = false_prior;
                const auto pooled_index = pooled_candidates.size();
                pooled_candidates.push_back(*candidate);
                pooled_evidence[sip_bin].target_indices.push_back(
                    pooled_index);
            }
        }
        if (false_prior_weight > 0) {
            calibration_false_priors[acceptor] = false_prior_sum /
                static_cast<double>(false_prior_weight);
        }
        add_elapsed(score_mbr_timing, wall_begin, cpu_begin);
    }
    double accepted_expected_false = 0.0;
    double null_tail_expected_false = 0.0;
    std::map<int, std::size_t> accepted_by_sip_bin;
    std::map<int, double> posterior_false_by_sip_bin;
    std::map<int, double> null_tail_false_by_sip_bin;
    {
        const auto wall_begin = Clock::now();
        const auto cpu_begin = std::clock();
        // Models remain acceptor/bin specific, but transferred-ion FDR is
        // estimated over all acceptors within the same SIP-abundance bin.
        struct NullTailGroup {
            double false_prior_sum = 0.0;
            double accepted_score_threshold =
                std::numeric_limits<double>::infinity();
            std::size_t accepted = 0;
            std::size_t decoys = 0;
            std::size_t decoys_above = 0;
        };
        std::map<std::pair<std::size_t, int>, NullTailGroup>
            null_tail_groups;
        for (auto& [sip_bin, evidence] : pooled_evidence) {
            std::vector<QuantTransferCandidate*> targets;
            targets.reserve(evidence.target_indices.size());
            for (const auto index : evidence.target_indices) {
                targets.push_back(&pooled_candidates[index]);
            }
            assign_quant_transfer_qvalues(targets);
            for (const auto* candidate : targets) {
                auto& group = null_tail_groups[{
                    candidate->psm.file_id, sip_bin}];
                group.false_prior_sum += candidate->false_prior;
                if (candidate->qvalue <= config.mbr_ion_fdr) {
                    ++group.accepted;
                    group.accepted_score_threshold = std::min(
                        group.accepted_score_threshold, candidate->score);
                    accepted_expected_false +=
                        1.0 - candidate->probability;
                    ++accepted_by_sip_bin[sip_bin];
                    posterior_false_by_sip_bin[sip_bin] +=
                        1.0 - candidate->probability;
                }
            }
            for (const auto* candidate : targets) {
                if (!std::isfinite(candidate->decoy_score)) continue;
                auto& group = null_tail_groups[{
                    candidate->psm.file_id, sip_bin}];
                ++group.decoys;
                if (group.accepted > 0 && candidate->decoy_score >=
                        group.accepted_score_threshold) {
                    ++group.decoys_above;
                }
            }
            for (const auto* candidate : targets) {
                if (candidate->qvalue > config.mbr_ion_fdr) continue;
                const auto& donor = data.rows[candidate->donor_row];
                data.transferred_ions.push_back({
                    candidate->psm, candidate->score,
                    candidate->qvalue, candidate->donor_row,
                    donor.id});
                if (candidate->psm.file_id < accepted_transfers.size()) {
                    ++accepted_transfers[candidate->psm.file_id];
                }
            }
        }
        for (const auto& [key, group] : null_tail_groups) {
            (void)key;
            if (group.accepted == 0 || group.decoys == 0) continue;
            const double expected_false = group.false_prior_sum *
                static_cast<double>(group.decoys_above) /
                static_cast<double>(group.decoys);
            null_tail_expected_false += expected_false;
            null_tail_false_by_sip_bin[key.second] += expected_false;
        }
        add_elapsed(score_mbr_timing, wall_begin, cpu_begin);
    }
    record_stage(
        "Build MBR run alignments (" +
            std::to_string(alignment_count) + ")",
        alignment_timing);
    record_stage(
        "Schedule MBR transfers (" +
            std::to_string(transfer_jobs) + " jobs)",
        scheduling_timing);
    record_stage(
        "Reload acceptor MS1 data for MBR", load_mbr_timing);
    record_stage(
        "Trace MBR target/decoy XICs", trace_mbr_timing, true);
    std::string calibration_summary;
    for (std::size_t file = 0; file < file_count; ++file) {
        calibration_summary += file == 0 ? "; " : ", ";
        calibration_summary += "run" + std::to_string(file + 1) +
            " pi0=";
        if (std::isfinite(calibration_false_priors[file])) {
            calibration_summary += std::to_string(
                calibration_false_priors[file]);
        } else {
            calibration_summary += "NA";
        }
        calibration_summary += " accepted=" +
            std::to_string(accepted_transfers[file]);
        if (calibration_lda_valid[file]) {
            calibration_summary += " lda=";
            for (std::size_t feature = 0;
                 feature < calibration_lda_weights[file].size(); ++feature) {
                if (feature != 0) calibration_summary += "/";
                calibration_summary += std::to_string(
                    calibration_lda_weights[file][feature]);
            }
        }
    }
    if (!config.sip_isotope.empty()) {
        calibration_summary += "; SIP-bin FDR audit";
        for (const auto& [sip_bin, accepted] : accepted_by_sip_bin) {
            const double denominator = static_cast<double>(accepted);
            calibration_summary += " b" + std::to_string(sip_bin) +
                "=" + std::to_string(accepted) +
                "/post:" + std::to_string(
                    posterior_false_by_sip_bin[sip_bin] / denominator) +
                "/null:" + std::to_string(
                    null_tail_false_by_sip_bin[sip_bin] / denominator);
        }
    }
    record_stage(
        "Fit covariance MBR LDA + four-population per-run probability/"
        "global ion FDR (" +
            std::to_string(training_target_count) + " +2, " +
            std::to_string(training_decoy_count) + " -2; " +
            std::to_string(scored_transfer_candidates) + " +1, " +
            std::to_string(transfer_decoy_count) + " -1; " +
            std::to_string(calibrated_transfer_candidates) +
            " calibrated, " +
            std::to_string(data.transferred_ions.size()) + " accepted" +
            "; posterior expected false=" +
            std::to_string(accepted_expected_false) + " (" +
            std::to_string(data.transferred_ions.empty() ? 0.0 :
                accepted_expected_false /
                    static_cast<double>(data.transferred_ions.size())) +
            "), shifted-decoy null-tail expected false=" +
            std::to_string(null_tail_expected_false) + " (" +
            std::to_string(data.transferred_ions.empty() ? 0.0 :
                null_tail_expected_false /
                    static_cast<double>(data.transferred_ions.size())) +
            ")" +
            calibration_summary + ")",
        score_mbr_timing);
    return result;
}

#endif

} // namespace aerith
