#include "pipeline.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <omp.h>

namespace aerith {

struct FastaEntry {
    std::string protein;
    std::string header;
    std::string sequence;
    std::string protein_id;
    std::string entry_name;
    std::string gene;
    std::string organism;
    std::string description;
    std::string existence;
};

struct PeptideEvidence {
    std::set<std::string> proteins;
    double probability = 0.0;
};

struct ProteinScore {
    std::string protein;
    double score = 0.0;
    double probability = 0.0;
    double qvalue = 1.0;
    bool decoy = false;
    bool picked = true;
};

struct ReportStats {
    std::set<std::string> total_peptides;
    std::set<std::string> unique_peptides;
    std::set<std::string> razor_peptides;
    std::size_t total_spectra = 0;
    std::size_t unique_spectra = 0;
    std::size_t razor_spectra = 0;
    std::unordered_map<std::string, double> total_ion_intensities;
    std::unordered_map<std::string, double> unique_ion_intensities;
    std::unordered_map<std::string, double> razor_ion_intensities;
    double total_intensity = 0.0;
    double unique_intensity = 0.0;
    double razor_intensity = 0.0;
    std::set<std::string> razor_assigned_modifications;
    std::set<std::string> razor_observed_modifications;
};

struct ModificationInfo {
    std::string sequence;
    std::string modified_peptide;
    std::vector<std::string> assigned;
    std::vector<std::string> observed;
};

std::string trim(std::string value) {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

bool starts_with(const std::string& value, const std::string& prefix) {
    return !prefix.empty() && value.rfind(prefix, 0) == 0;
}

bool ptm_mass(char token, double& mass) {
    switch (token) {
    case '~': mass = 15.994915; return true;
    case '!': mass = 0.984016; return true;
    case '@':
    case '>':
    case '<': mass = 79.966332; return true;
    case '%': mass = 42.010565; return true;
    case '^': mass = 14.015650; return true;
    case '&': mass = 28.031300; return true;
    case '*': mass = 42.046950; return true;
    case '(': mass = 28.990164; return true;
    case ')': mass = 44.985079; return true;
    case '/': mass = 57.021464; return true;
    case '$': mass = 45.987721; return true;
    default: return false;
    }
}

std::string mass_text(double mass) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << mass;
    return stream.str();
}

ModificationInfo modification_info(
    const std::string& decorated, bool fixed_cam) {
    const auto open = decorated.find('[');
    const auto close = decorated.rfind(']');
    const std::string body =
        open != std::string::npos && close != std::string::npos && close > open
        ? decorated.substr(open + 1, close - open - 1)
        : decorated;
    struct Residue {
        char amino_acid = '\0';
        std::vector<char> tokens;
    };
    std::vector<Residue> residues;
    std::vector<char> nterm_tokens;
    for (const char value : body) {
        if (std::isalpha(static_cast<unsigned char>(value)) != 0) {
            residues.push_back({static_cast<char>(
                std::toupper(static_cast<unsigned char>(value))), {}});
            continue;
        }
        double mass = 0.0;
        if (!ptm_mass(value, mass)) continue;
        if (residues.empty()) nterm_tokens.push_back(value);
        else residues.back().tokens.push_back(value);
    }

    ModificationInfo result;
    for (const auto token : nterm_tokens) {
        double mass = 0.0;
        if (!ptm_mass(token, mass)) continue;
        result.modified_peptide += "n[+" + mass_text(mass) + "]";
        result.assigned.push_back("N-term(" + mass_text(mass) + ")");
    }
    for (std::size_t index = 0; index < residues.size(); ++index) {
        const auto& residue = residues[index];
        result.sequence.push_back(residue.amino_acid);
        result.modified_peptide.push_back(residue.amino_acid);
        const bool replaces_cam =
            std::find(residue.tokens.begin(), residue.tokens.end(), '(') !=
                residue.tokens.end() ||
            std::find(residue.tokens.begin(), residue.tokens.end(), '/') !=
                residue.tokens.end();
        if (fixed_cam && residue.amino_acid == 'C' && !replaces_cam) {
            constexpr double cam = 57.021464;
            result.modified_peptide += "[+" + mass_text(cam) + "]";
            result.assigned.push_back(
                std::to_string(index + 1) + "C(" + mass_text(cam) + ")");
        }
        for (const auto token : residue.tokens) {
            double mass = 0.0;
            if (!ptm_mass(token, mass)) continue;
            result.modified_peptide += "[+" + mass_text(mass) + "]";
            result.assigned.push_back(
                std::to_string(index + 1) + residue.amino_acid +
                "(" + mass_text(mass) + ")");
        }
    }
    if (result.assigned.empty()) result.modified_peptide.clear();
    return result;
}

double psm_intensity(const Psm& psm) {
    if (!(psm.log10_precursor_intensity > 0.0) ||
        psm.log10_precursor_intensity > 300.0) {
        return 0.0;
    }
    return std::pow(10.0, psm.log10_precursor_intensity);
}

std::string ion_form(const Psm& psm) {
    std::ostringstream stream;
    stream << psm.peptide << '#' << psm.charge << '#'
           << std::fixed << std::setprecision(4) << psm.exp_mass;
    return stream.str();
}

void update_ion_intensity(
    std::unordered_map<std::string, double>& intensities,
    const std::string& ion, double intensity) {
    auto& current = intensities[ion];
    current = std::max(current, intensity);
}

double top_three_intensity(
    const std::unordered_map<std::string, double>& intensities) {
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

std::string join_values(const std::set<std::string>& values) {
    std::ostringstream stream;
    bool first = true;
    for (const auto& value : values) {
        if (!first) stream << ", ";
        first = false;
        stream << value;
    }
    return stream.str();
}

std::string sanitize_protein(std::string protein, const std::string& decoy_prefix) {
    protein = trim(std::move(protein));
    std::string active_prefix;
    for (const auto& prefix : {decoy_prefix, std::string("DECOY_"),
                               std::string("Decoy_")}) {
        if (starts_with(protein, prefix)) {
            active_prefix = decoy_prefix;
            protein.erase(0, prefix.size());
            break;
        }
    }
    if (protein.find('|') == std::string::npos) {
        std::replace_if(protein.begin(), protein.end(), [](char ch) {
            return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n';
        }, '_');
        protein = "sp|" + protein + "|" + protein;
    }
    return active_prefix + protein;
}

std::vector<std::string> protein_ids(
    const Psm& psm, const std::string& decoy_prefix) {
    std::string value = psm.proteins;
    if (!value.empty() && value.front() == '{') value.erase(value.begin());
    if (!value.empty() && value.back() == '}') value.pop_back();
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto comma = value.find(',', begin);
        const auto end = comma == std::string::npos ? value.size() : comma;
        auto protein = trim(value.substr(begin, end - begin));
        if (!protein.empty()) {
            result.push_back(sanitize_protein(std::move(protein), decoy_prefix));
        }
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::string annotation_value(
    const std::string& annotation, const std::string& key) {
    const auto begin = annotation.find(key);
    if (begin == std::string::npos) return {};
    const auto value_begin = begin + key.size();
    auto end = annotation.size();
    for (const auto& next : {" OS=", " OX=", " GN=", " PE=", " SV="}) {
        const auto found = annotation.find(next, value_begin);
        if (found != std::string::npos) end = std::min(end, found);
    }
    return trim(annotation.substr(value_begin, end - value_begin));
}

std::string protein_existence(const std::string& code) {
    if (code == "1") return "1:Experimental evidence at protein level";
    if (code == "2") return "2:Experimental evidence at transcript level";
    if (code == "3") return "3:Protein inferred from homology";
    if (code == "4") return "4:Protein predicted";
    if (code == "5") return "5:Protein uncertain";
    return {};
}

FastaEntry make_fasta_entry(
    const std::string& raw_header, const std::string& sequence,
    const std::string& decoy_prefix) {
    const auto space = raw_header.find_first_of(" \t");
    const auto raw_id = raw_header.substr(0, space);
    const auto annotation = space == std::string::npos
        ? std::string{} : trim(raw_header.substr(space + 1));
    FastaEntry entry;
    entry.protein = sanitize_protein(raw_id, decoy_prefix);
    entry.header = entry.protein + (annotation.empty() ? "" : " " + annotation);
    entry.sequence = sequence;
    auto core = entry.protein;
    if (starts_with(core, decoy_prefix)) core.erase(0, decoy_prefix.size());
    const auto first = core.find('|');
    const auto second = first == std::string::npos
        ? std::string::npos : core.find('|', first + 1);
    if (first != std::string::npos && second != std::string::npos) {
        entry.protein_id = core.substr(first + 1, second - first - 1);
        entry.entry_name = core.substr(second + 1);
    } else {
        entry.protein_id = core;
        entry.entry_name = core;
    }
    entry.gene = annotation_value(" " + annotation, " GN=");
    entry.organism = annotation_value(" " + annotation, " OS=");
    const auto os = annotation.find(" OS=");
    if (os != std::string::npos) entry.description = trim(annotation.substr(0, os));
    else if (!starts_with(annotation, "OS=")) entry.description = trim(annotation);
    if (entry.description.empty()) {
        entry.description = entry.organism.empty()
            ? entry.entry_name : "OS=" + entry.organism;
    }
    entry.existence = protein_existence(annotation_value(" " + annotation, " PE="));
    return entry;
}

std::unordered_map<std::string, FastaEntry> read_fasta(
    const std::string& path, const std::string& decoy_prefix,
    bool force_decoy) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("Cannot read protein database: " + path);
    std::unordered_map<std::string, FastaEntry> result;
    std::string line;
    std::string header;
    std::string sequence;
    const auto flush = [&]() {
        if (header.empty()) return;
        auto normalized_header = header;
        const auto first_space = normalized_header.find_first_of(" \t");
        const auto raw_id = normalized_header.substr(0, first_space);
        if (force_decoy &&
            !starts_with(sanitize_protein(raw_id, decoy_prefix), decoy_prefix)) {
            normalized_header = decoy_prefix + normalized_header;
        }
        auto entry = make_fasta_entry(
            normalized_header, sequence, decoy_prefix);
        result[entry.protein] = std::move(entry);
    };
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            flush();
            header = trim(line.substr(1));
            sequence.clear();
        } else {
            for (const char ch : line) {
                if (ch != ' ' && ch != '\t' && ch != '\r') sequence += ch;
            }
        }
    }
    flush();
    return result;
}

std::filesystem::path default_output_dir(const Config& config) {
    if (!config.protein_output_dir.empty()) return config.protein_output_dir;
    if (config.output_prefixes.empty()) return std::filesystem::current_path();
    auto parent = std::filesystem::path(config.output_prefixes.front()).parent_path();
    if (config.output_prefixes.size() > 1 && parent.has_parent_path()) {
        return parent.parent_path();
    }
    return parent.empty() ? std::filesystem::current_path() : parent;
}

double protein_probability(
    const std::set<std::string>& peptides,
    const std::unordered_map<std::string, PeptideEvidence>& evidence) {
    long double failure = 1.0L;
    for (const auto& peptide : peptides) {
        failure *= 1.0L - std::clamp(evidence.at(peptide).probability, 0.0, 1.0);
    }
    return static_cast<double>(1.0L - failure);
}

double proteinprophet_peptide_probability(double probability) {
    // ProteinProphet represents its peptide-confidence histogram at the
    // midpoint of 0.002-wide bins (0.991, 0.993, 0.995, ...).
    const double midpoint = std::round(std::clamp(probability, 0.0, 1.0) * 500.0)
        * 0.002 - 0.001;
    return std::clamp(midpoint, 0.0, 0.999);
}

double coverage(const FastaEntry& entry, const std::set<std::string>& peptides) {
    if (entry.sequence.empty()) return 0.0;
    std::vector<bool> covered(entry.sequence.size(), false);
    auto normalized_protein = entry.sequence;
    std::replace(normalized_protein.begin(), normalized_protein.end(), 'I', 'L');
    for (const auto& peptide : peptides) {
        auto normalized_peptide = peptide;
        std::replace(normalized_peptide.begin(), normalized_peptide.end(), 'I', 'L');
        std::size_t position = 0;
        while (!normalized_peptide.empty() &&
               (position = normalized_protein.find(normalized_peptide, position)) !=
                   std::string::npos) {
            const auto end = std::min(covered.size(), position + normalized_peptide.size());
            std::fill(covered.begin() + static_cast<std::ptrdiff_t>(position),
                      covered.begin() + static_cast<std::ptrdiff_t>(end), true);
            ++position;
        }
    }
    const auto count = std::count(covered.begin(), covered.end(), true);
    return 100.0 * static_cast<double>(count) /
           static_cast<double>(entry.sequence.size());
}

std::string signature(const std::set<std::string>& peptides) {
    std::string value;
    for (const auto& peptide : peptides) {
        value += peptide;
        value.push_back('\x1f');
    }
    return value;
}

std::string join_values(const std::vector<std::string>& values) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) stream << ", ";
        stream << values[index];
    }
    return stream.str();
}

std::pair<std::size_t, std::size_t> protein_coordinates(
    const FastaEntry& entry, std::string peptide) {
    auto protein = entry.sequence;
    std::replace(protein.begin(), protein.end(), 'I', 'L');
    std::replace(peptide.begin(), peptide.end(), 'I', 'L');
    const auto position = protein.find(peptide);
    if (position == std::string::npos) return {0, 0};
    return {position + 1, position + peptide.size()};
}

std::pair<std::string, std::string> peptide_flanks(
    const std::string& decorated) {
    const auto open = decorated.find('[');
    const auto close = decorated.rfind(']');
    std::string previous = open == std::string::npos
        ? "-" : decorated.substr(0, open);
    std::string next = close == std::string::npos
        ? "-" : decorated.substr(close + 1);
    if (previous.empty()) previous = "-";
    if (next.empty()) next = "-";
    return {previous, next};
}

void write_psm_report(
    const std::filesystem::path& directory, std::size_t file,
    double psm_threshold, const Config& config, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep,
    const std::unordered_map<std::string, FastaEntry>& fasta,
    const std::unordered_map<std::string, PeptideEvidence>& peptides,
    const std::unordered_map<std::string, std::string>& assignment,
    const std::set<std::string>& accepted) {
    const auto all_files = std::numeric_limits<std::size_t>::max();
    const auto report_path = directory / "psm.tsv";
    std::ofstream report(report_path);
    if (!report) {
        throw std::runtime_error("Cannot create output: " + report_path.string());
    }
    report
        << "Spectrum\tSpectrum File\tPeptide\tModified Peptide\tExtended Peptide"
        << "\tPrev AA\tNext AA\tPeptide Length\tCharge\tRetention"
        << "\tObserved Mass\tCalibrated Observed Mass\tObserved M/Z"
        << "\tCalibrated Observed M/Z\tCalculated Peptide Mass\tCalculated M/Z"
        << "\tDelta Mass\tExpectation\tHyperscore\tNextscore\tProbability"
        << "\tQvalue\tNumber of Enzymatic Termini\tNumber of Missed Cleavages"
        << "\tProtein Start\tProtein End\tIntensity\tAssigned Modifications"
        << "\tObserved Modifications\tClass\tPurity\tIs Decoy\tIs Contaminant"
        << "\tIs Unique\tProtein\tProtein ID\tEntry Name\tGene"
        << "\tProtein Description\tMapped Genes\tMapped Proteins\n";
    report << std::fixed;
    std::vector<std::size_t> rows;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto& psm = data.rows[row];
        if ((file == all_files || psm.file_id == file) && psm.label == 1 &&
            1.0 - pep[row] >= psm_threshold) {
            rows.push_back(row);
        }
    }
    std::stable_sort(rows.begin(), rows.end(), [&](std::size_t left,
                                                   std::size_t right) {
        const auto& lhs = data.rows[left];
        const auto& rhs = data.rows[right];
        if (lhs.file_id != rhs.file_id) return lhs.file_id < rhs.file_id;
        if (lhs.scan != rhs.scan) return lhs.scan < rhs.scan;
        return lhs.rank < rhs.rank;
    });
    constexpr double proton = 1.007276466621;
    for (const auto row : rows) {
        const auto& psm = data.rows[row];
        const auto peptide = stripped_peptide(psm.peptide);
        const auto evidence = peptides.find(peptide);
        const auto razor = assignment.find(peptide);
        if (evidence == peptides.end() || razor == assignment.end() ||
            accepted.count(razor->second) == 0) {
            continue;
        }
        const auto fasta_it = fasta.find(razor->second);
        FastaEntry fallback;
        if (fasta_it == fasta.end()) {
            fallback = make_fasta_entry(
                razor->second, {}, config.decoy_prefix);
        }
        const auto& entry =
            fasta_it == fasta.end() ? fallback : fasta_it->second;
        const auto modifications =
            modification_info(psm.peptide, config.fixed_cam);
        const auto [previous, next] = peptide_flanks(psm.peptide);
        const auto [protein_start, protein_end] =
            protein_coordinates(entry, peptide);
        const double mz =
            (psm.exp_mass + static_cast<double>(psm.charge) * proton) /
            static_cast<double>(psm.charge);
        std::vector<std::string> mapped_proteins;
        std::set<std::string> mapped_genes;
        for (const auto& protein : protein_ids(psm, config.decoy_prefix)) {
            if (protein == razor->second || starts_with(
                    protein, config.decoy_prefix)) {
                continue;
            }
            mapped_proteins.push_back(protein);
            const auto mapped = fasta.find(protein);
            if (mapped != fasta.end() && !mapped->second.gene.empty()) {
                mapped_genes.insert(mapped->second.gene);
            }
        }
        const std::string spectrum_file =
            psm.file_id < config.output_prefixes.size()
            ? std::filesystem::path(
                  config.output_prefixes[psm.file_id]).filename().string()
            : std::filesystem::path(
                  data.input_paths[psm.file_id]).filename().string();
        report << psm.id << '\t' << spectrum_file << '\t'
               << modifications.sequence << '\t'
               << modifications.modified_peptide << '\t'
               << previous << '.'
               << (modifications.modified_peptide.empty()
                       ? modifications.sequence
                       : modifications.modified_peptide)
               << '.' << next << '\t' << previous << '\t' << next << '\t'
               << modifications.sequence.size() << '\t' << psm.charge << '\t'
               << std::setprecision(4) << psm.retention << '\t'
               << psm.exp_mass << '\t' << psm.exp_mass << '\t'
               << mz << '\t' << mz << '\t' << psm.exp_mass << '\t' << mz
               << '\t' << psm.mass_error << "\t0\t"
               << scores[row] << '\t' << psm.diff_score << '\t'
               << 1.0 - pep[row] << '\t' << std::setprecision(14) << q[row]
               << "\t2\t" << psm.missed_cleavages << '\t'
               << protein_start << '\t' << protein_end << '\t'
               << std::setprecision(2) << psm_intensity(psm) << '\t'
               << join_values(modifications.assigned) << '\t'
               << join_values(modifications.observed)
               << "\tTarget\t0.00\tfalse\tfalse\t"
               << (evidence->second.proteins.size() == 1 ? "true" : "false")
               << '\t' << razor->second << '\t' << entry.protein_id << '\t'
               << entry.entry_name << '\t' << entry.gene << '\t'
               << entry.description << '\t' << join_values(mapped_genes)
               << '\t';
        for (std::size_t index = 0; index < mapped_proteins.size(); ++index) {
            if (index) report << ", ";
            report << mapped_proteins[index];
        }
        report << '\n';
    }
}

void write_report(
    const std::filesystem::path& directory, std::size_t file,
    const Config& config, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep,
    const std::unordered_map<std::string, FastaEntry>& fasta,
    const std::unordered_map<std::string, PeptideEvidence>& peptides,
    const std::unordered_map<std::string, std::string>& assignment,
    const std::unordered_map<std::string, ProteinScore>& protein_scores,
    const std::set<std::string>& accepted) {
    const auto all_files = std::numeric_limits<std::size_t>::max();
    struct RankedPsm {
        double probability;
        bool decoy;
    };
    std::vector<RankedPsm> ranked;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        if ((file == all_files || data.rows[row].file_id == file) && pep[row] < 0.5) {
            ranked.push_back({1.0 - pep[row], data.rows[row].label == -1});
        }
    }
    std::sort(ranked.begin(), ranked.end(), [](const auto& left, const auto& right) {
        return left.probability > right.probability;
    });
    std::size_t target_count = 0;
    std::size_t decoy_count = 0;
    double psm_threshold = 1.0;
    for (std::size_t begin = 0; begin < ranked.size();) {
        std::size_t end = begin;
        while (end < ranked.size() && ranked[end].probability == ranked[begin].probability) {
            ranked[end].decoy ? ++decoy_count : ++target_count;
            ++end;
        }
        if (target_count > 0 && static_cast<double>(decoy_count) / target_count <= 0.01) {
            psm_threshold = ranked[begin].probability;
        }
        begin = end;
    }
    std::unordered_map<std::string, ReportStats> stats;
    std::unordered_map<std::string, std::set<std::string>>
        observed_protein_peptides;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        const auto& psm = data.rows[row];
        if ((file != all_files && psm.file_id != file) ||
            psm.label != 1 || 1.0 - pep[row] < psm_threshold) {
            continue;
        }
        const auto peptide = stripped_peptide(psm.peptide);
        const auto evidence = peptides.find(peptide);
        const auto razor = assignment.find(peptide);
        if (evidence == peptides.end() || razor == assignment.end() ||
            accepted.count(razor->second) == 0) {
            continue;
        }
        const auto mappings = protein_ids(psm, config.decoy_prefix);
        const bool unique = evidence->second.proteins.size() == 1;
        const double intensity = psm_intensity(psm);
        const auto ion = ion_form(psm);
        const auto modifications =
            modification_info(psm.peptide, config.fixed_cam);
        for (const auto& protein : mappings) {
            if (!starts_with(protein, config.decoy_prefix)) {
                observed_protein_peptides[protein].insert(peptide);
            }
            if (accepted.count(protein) == 0) continue;
            auto& value = stats[protein];
            value.total_peptides.insert(peptide);
            ++value.total_spectra;
            update_ion_intensity(
                value.total_ion_intensities, ion, intensity);
            if (unique) {
                value.unique_peptides.insert(peptide);
                ++value.unique_spectra;
                update_ion_intensity(
                    value.unique_ion_intensities, ion, intensity);
            }
            if (razor->second == protein) {
                value.razor_peptides.insert(peptide);
                ++value.razor_spectra;
                update_ion_intensity(
                    value.razor_ion_intensities, ion, intensity);
                value.razor_assigned_modifications.insert(
                    modifications.assigned.begin(),
                    modifications.assigned.end());
                value.razor_observed_modifications.insert(
                    modifications.observed.begin(),
                    modifications.observed.end());
            }
        }
    }
    for (auto& [protein, value] : stats) {
        (void)protein;
        value.total_intensity =
            top_three_intensity(value.total_ion_intensities);
        value.unique_intensity =
            top_three_intensity(value.unique_ion_intensities);
        value.razor_intensity =
            top_three_intensity(value.razor_ion_intensities);
    }

    std::vector<std::string> proteins;
    for (const auto& [protein, value] : stats) {
        if (!value.razor_peptides.empty()) proteins.push_back(protein);
    }
    std::sort(proteins.begin(), proteins.end(), [&](const auto& left, const auto& right) {
        const auto& left_score = protein_scores.at(left);
        const auto& right_score = protein_scores.at(right);
        if (left_score.probability != right_score.probability) {
            return left_score.probability > right_score.probability;
        }
        if (left_score.score != right_score.score) return left_score.score > right_score.score;
        return left < right;
    });

    std::unordered_map<std::string, std::vector<std::string>> equivalent;
    for (const auto& [protein, observed_peptides] :
         observed_protein_peptides) {
        equivalent[signature(observed_peptides)].push_back(protein);
    }
    for (auto& [peptide_signature, members] : equivalent) {
        (void)peptide_signature;
        std::sort(members.begin(), members.end());
    }
    std::filesystem::create_directories(directory);
    const auto table_path = directory / "protein.tsv";
    std::ofstream table(table_path);
    if (!table) throw std::runtime_error("Cannot create output: " + table_path.string());
    table << "Protein\tProtein ID\tEntry Name\tGene\tLength\tIs Decoy\tIs Contaminant"
          << "\tOrganism\tProtein Description\tProtein Existence\tCoverage"
          << "\tProtein Probability\tTop Peptide Probability\tProtein Qvalue"
          << "\tTotal Peptides\tUnique Peptides\tRazor Peptides\tTotal Spectral Count"
          << "\tUnique Spectral Count\tRazor Spectral Count\tTotal Intensity"
          << "\tUnique Intensity\tRazor Intensity\tRazor Assigned Modifications"
          << "\tRazor Observed Modifications\tIndistinguishable Proteins\n";
    table << std::fixed;
    for (const auto& protein : proteins) {
        const auto fasta_it = fasta.find(protein);
        FastaEntry fallback;
        if (fasta_it == fasta.end()) {
            fallback = make_fasta_entry(protein, {}, config.decoy_prefix);
        }
        const auto& entry = fasta_it == fasta.end() ? fallback : fasta_it->second;
        const auto& value = stats.at(protein);
        const auto& score = protein_scores.at(protein);
        std::vector<std::string> indistinguishable;
        for (const auto& other : equivalent.at(signature(value.total_peptides))) {
            if (other != protein) indistinguishable.push_back(other);
        }
        table << protein << '\t' << entry.protein_id << '\t' << entry.entry_name
              << '\t' << entry.gene << '\t' << entry.sequence.size()
              << "\tfalse\tfalse\t" << entry.organism << '\t' << entry.description
              << '\t' << entry.existence << '\t' << std::setprecision(2)
              << coverage(entry, value.total_peptides) << '\t' << std::setprecision(4)
              << score.probability << '\t' << score.score << '\t'
              << std::setprecision(14) << score.qvalue << '\t'
              << value.total_peptides.size() << '\t' << value.unique_peptides.size()
              << '\t' << value.razor_peptides.size() << '\t' << value.total_spectra
              << '\t' << value.unique_spectra << '\t' << value.razor_spectra
              << '\t' << std::setprecision(0) << value.total_intensity
              << '\t' << value.unique_intensity << '\t' << value.razor_intensity
              << '\t' << join_values(value.razor_assigned_modifications)
              << '\t' << join_values(value.razor_observed_modifications)
              << '\t';
        for (std::size_t i = 0; i < indistinguishable.size(); ++i) {
            if (i) table << ", ";
            table << indistinguishable[i];
        }
        table << '\n';
    }
    write_psm_report(
        directory, file, psm_threshold, config, data, scores, q, pep, fasta,
        peptides, assignment, accepted);

    const auto fasta_path = directory / "protein.fas";
    std::ofstream selected_fasta(fasta_path);
    if (!selected_fasta) {
        throw std::runtime_error("Cannot create output: " + fasta_path.string());
    }
    for (const auto& protein : proteins) {
        const auto found = fasta.find(protein);
        if (found == fasta.end()) continue;
        selected_fasta << '>' << found->second.header << '\n'
                       << found->second.sequence << '\n';
    }
}

ProteinAssemblyResult ProteinAssembler::write(
    const Config& config, const Dataset& data,
    const std::vector<double>& scores, const std::vector<double>& q,
    const std::vector<double>& pep) {
    using Clock = std::chrono::steady_clock;
    if (data.rows.size() != scores.size() || scores.size() != q.size() ||
        q.size() != pep.size()) {
        throw std::runtime_error("Protein assembly received inconsistent PSM arrays");
    }
    ProteinAssemblyResult result;
    const auto record_stage = [&](std::string name, Clock::time_point wall_begin,
                                  std::clock_t cpu_begin, bool uses_omp = false,
                                  bool uses_simd = false) {
        const auto wall_end = Clock::now();
        const auto cpu_end = std::clock();
        result.stages.push_back({
            std::move(name),
            {std::chrono::duration<double>(wall_end - wall_begin).count(),
             static_cast<double>(cpu_end - cpu_begin) / CLOCKS_PER_SEC},
            uses_omp,
            uses_simd});
    };

    auto wall_begin = Clock::now();
    auto cpu_begin = std::clock();
    auto fasta = read_fasta(
        config.database_path, config.decoy_prefix, false);
    auto decoy_fasta = read_fasta(
        config.decoy_database_path, config.decoy_prefix, true);
    for (auto& [protein, entry] : decoy_fasta) {
        fasta.insert_or_assign(protein, std::move(entry));
    }
    record_stage("Read and annotate target/decoy FASTAs", wall_begin, cpu_begin);

    wall_begin = Clock::now();
    cpu_begin = std::clock();
    std::unordered_map<std::string, PeptideEvidence> peptides;
    std::unordered_map<std::string, std::set<std::string>> protein_peptides;
    std::unordered_map<std::string, std::set<std::string>> unique_peptides;
    for (std::size_t row = 0; row < data.rows.size(); ++row) {
        if (!(pep[row] < 0.5)) continue;
        const auto sequence = stripped_peptide(data.rows[row].peptide);
        if (sequence.empty()) continue;
        auto& evidence = peptides[sequence];
        evidence.probability = std::max(
            evidence.probability,
            proteinprophet_peptide_probability(1.0 - pep[row]));
        for (const auto& protein : protein_ids(data.rows[row], config.decoy_prefix)) {
            evidence.proteins.insert(protein);
            protein_peptides[protein].insert(sequence);
        }
    }
    for (const auto& [peptide, evidence] : peptides) {
        if (evidence.proteins.size() == 1) {
            unique_peptides[*evidence.proteins.begin()].insert(peptide);
        }
    }
    record_stage("Build peptide-protein evidence maps", wall_begin, cpu_begin);

    wall_begin = Clock::now();
    cpu_begin = std::clock();
    std::unordered_map<std::string, std::string> assignment;
    std::unordered_map<std::string, std::set<std::string>> assigned_peptides;
    // ProteinProphet's exact EM/NSP model is tied to its XML parser. This
    // deterministic in-memory analogue preserves its observable decisions:
    // favor the protein with the largest peptide set, use unique evidence and
    // stable accession ordering as ties, then score the resulting razor set.
    for (const auto& [peptide, evidence] : peptides) {
        if (evidence.proteins.empty()) continue;
        auto best = evidence.proteins.begin();
        for (auto candidate = std::next(evidence.proteins.begin());
             candidate != evidence.proteins.end(); ++candidate) {
            const auto best_total = protein_peptides[*best].size();
            const auto candidate_total = protein_peptides[*candidate].size();
            const auto best_unique = unique_peptides[*best].size();
            const auto candidate_unique = unique_peptides[*candidate].size();
            if (candidate_total > best_total ||
                (candidate_total == best_total && candidate_unique > best_unique) ||
                (candidate_total == best_total && candidate_unique == best_unique &&
                 (candidate->size() < best->size() ||
                  (candidate->size() == best->size() && *candidate < *best)))) {
                best = candidate;
            }
        }
        assignment[peptide] = *best;
        assigned_peptides[*best].insert(peptide);
    }
    record_stage("Assign unique and razor peptides", wall_begin, cpu_begin);

    wall_begin = Clock::now();
    cpu_begin = std::clock();
    std::unordered_map<std::string, ProteinScore> protein_scores;
    for (const auto& [protein, assigned] : assigned_peptides) {
        if (assigned.empty()) continue;
        ProteinScore value;
        value.protein = protein;
        value.decoy = starts_with(protein, config.decoy_prefix);
        for (const auto& peptide : assigned) {
            value.score = std::max(value.score, peptides.at(peptide).probability);
        }
        value.probability = protein_probability(assigned, peptides);
        protein_scores[protein] = std::move(value);
    }

    for (auto& [protein, value] : protein_scores) {
        if (value.decoy) continue;
        const auto decoy = config.decoy_prefix + protein;
        const auto found = protein_scores.find(decoy);
        if (found == protein_scores.end()) continue;
        if (value.score > found->second.score) found->second.picked = false;
        else if (found->second.score > value.score) value.picked = false;
    }
    std::vector<ProteinScore*> ordered;
    for (auto& [protein, value] : protein_scores) {
        (void)protein;
        if (value.picked) ordered.push_back(&value);
    }
    std::sort(ordered.begin(), ordered.end(), [](const auto* left, const auto* right) {
        if (left->score != right->score) return left->score > right->score;
        return left->protein < right->protein;
    });
    std::size_t targets = 0;
    std::size_t decoys = 0;
    std::size_t accepted_end = 0;
    std::vector<double> raw_fdr(ordered.size(), 1.0);
    for (std::size_t begin = 0; begin < ordered.size();) {
        std::size_t end = begin;
        while (end < ordered.size() && ordered[end]->score == ordered[begin]->score) {
            ordered[end]->decoy ? ++decoys : ++targets;
            ++end;
        }
        const double fdr = targets == 0 ? 1.0
            : static_cast<double>(decoys) / static_cast<double>(targets);
        for (std::size_t i = begin; i < end; ++i) raw_fdr[i] = fdr;
        if (fdr <= 0.01) accepted_end = end;
        begin = end;
    }
    double minimum_q = 1.0;
    for (std::size_t i = ordered.size(); i-- > 0;) {
        minimum_q = std::min(minimum_q, raw_fdr[i]);
        ordered[i]->qvalue = minimum_q;
    }
    std::set<std::string> accepted;
    for (std::size_t i = 0; i < accepted_end; ++i) {
        if (!ordered[i]->decoy) accepted.insert(ordered[i]->protein);
    }
    record_stage("Picked competition and block protein FDR", wall_begin, cpu_begin);

    wall_begin = Clock::now();
    cpu_begin = std::clock();
    const auto output_dir = default_output_dir(config);
    const auto all_files = std::numeric_limits<std::size_t>::max();
    std::vector<std::pair<std::filesystem::path, std::size_t>> reports;
    reports.emplace_back(output_dir, all_files);
    if (config.output_prefixes.size() == data.input_paths.size()) {
        for (std::size_t file = 0; file < data.input_paths.size(); ++file) {
            const auto directory =
                std::filesystem::path(config.output_prefixes[file]).parent_path();
            reports.emplace_back(directory, file);
        }
    }
    std::set<std::filesystem::path> report_directories;
    bool distinct_report_directories = true;
    for (const auto& [directory, file] : reports) {
        (void)file;
        distinct_report_directories &= report_directories.insert(directory).second;
        std::filesystem::create_directories(directory);
    }
    const bool parallel_reports = reports.size() > 1 &&
        distinct_report_directories && omp_get_max_threads() > 1;
    std::exception_ptr report_failure;
    #pragma omp parallel for schedule(dynamic) if(parallel_reports)
    for (std::ptrdiff_t report = 0;
         report < static_cast<std::ptrdiff_t>(reports.size()); ++report) {
        try {
            const auto& [directory, file] = reports[static_cast<std::size_t>(report)];
            write_report(
                directory, file, config, data, scores, q, pep, fasta, peptides,
                assignment, protein_scores, accepted);
        } catch (...) {
            #pragma omp critical(aerith_protein_report_failure)
            if (!report_failure) report_failure = std::current_exception();
        }
    }
    if (report_failure) std::rethrow_exception(report_failure);
    record_stage(
        "Build and write combined/sample reports", wall_begin, cpu_begin,
        parallel_reports);

    result.proteins = accepted.size();
    result.output_dir = output_dir.string();
    return result;
}

} // namespace aerith
