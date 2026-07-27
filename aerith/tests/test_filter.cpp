#include "filter.hpp"
#include "pipeline.hpp"
#include "quantification.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
    const auto split_tabs = [](const std::string& line) {
        std::vector<std::string> fields;
        std::size_t begin = 0;
        while (true) {
            const auto tab = line.find('\t', begin);
            fields.push_back(line.substr(begin, tab - begin));
            if (tab == std::string::npos) break;
            begin = tab + 1;
        }
        return fields;
    };
    assert(aerith::stripped_peptide("K[PEPTM~IDE]R") == "PEPTMIDE");
    assert(aerith::stripped_peptide("-.ACD[+57.0]EF.-") == "ACDEF");
    assert(aerith::stripped_peptide("[%M~PEPTIDE]R") == "MPEPTIDE");
    assert(aerith::stripped_peptide("K[PEPN!IDE]") == "PEPNIDE");
    assert(aerith::stripped_peptide("[PEPTIDE]R") == "PEPTIDE");
    aerith::Psm intensity_psm;
    intensity_psm.log10_precursor_intensity = 3.0;
    assert(std::abs(aerith::psm_intensity(intensity_psm) - 1000.0) < 1e-12);
    intensity_psm.quantification_attempted = true;
    assert(aerith::psm_intensity(intensity_psm) == 0.0);
    intensity_psm.has_chromatographic_feature = true;
    intensity_psm.quantified_intensity = 500.0;
    assert(aerith::psm_intensity(intensity_psm) == 500.0);
    const std::array<double, 3> isotope_heights{25.0, 4.0, 10.0};
    assert(std::abs(
        aerith::summed_isotope_apex_intensity(isotope_heights) - 39.0) <
        1e-12);

    const std::vector<double> scores{10.0, 9.0, 8.0, 7.0, 6.0};
    const std::vector<int> labels{1, 1, -1, 1, -1};
    const auto q = aerith::target_decoy_qvalues(scores, labels);
    assert(q.size() == scores.size());
    assert(std::abs(q[0] - 0.5) < 1e-12);
    assert(std::abs(q[1] - 0.5) < 1e-12);
    assert(q[2] <= q[4]);
    assert(q[3] <= q[4]);

    double pi0 = 0.0;
    const auto mixmax = aerith::mixmax_qvalues(scores, labels, &pi0);
    assert(mixmax.size() == scores.size());
    assert(pi0 >= 0.0 && pi0 <= 1.0);
    assert(mixmax[0] <= mixmax[1]);
    assert(mixmax[1] <= mixmax[3]);

    const auto collision_root = std::filesystem::temp_directory_path() /
        "aerith_target_decoy_peptide_collision_unit";
    std::filesystem::remove_all(collision_root);
    std::filesystem::create_directories(collision_root);
    const auto write_pin = [](
        const std::filesystem::path& path,
        const std::vector<std::array<std::string, 6>>& rows) {
        std::ofstream pin(path);
        pin << "SpecId\tLabel\tScanNr\tretentiontime\tExpMass\tWDPscores"
               "\tranks\tparentCharges\tPeptide\tProteins\n";
        for (const auto& row : rows) {
            pin << row[0] << '\t' << row[1] << '\t' << row[2]
                << "\t1.0\t1000.0\t" << row[3] << "\t1\t2\t"
                << row[4] << '\t' << row[5] << '\n';
        }
    };
    const auto target0 = collision_root / "target0.pin";
    const auto decoy0 = collision_root / "decoy0.pin";
    const auto target1 = collision_root / "target1.pin";
    const auto decoy1 = collision_root / "decoy1.pin";
    write_pin(target0, {{
        "sample0.1.1", "1", "1", "10.0", "K[PEPTM~IDE]R", "{sp|P1|ONE}"}});
    write_pin(decoy0, {
        {"sample0.1.1", "-1", "1", "10.0", "R[PEPTMIDE]K", "{Decoy_P1}"},
        {"sample0.2.1", "-1", "2", "9.0", "R[UNIQUE]K", "{Decoy_P2}"}});
    write_pin(target1, {{
        "sample1.1.1", "1", "1", "10.0", "K[ANOTHER]R", "{sp|P2|TWO}"}});
    write_pin(decoy1, {
        {"sample1.2.1", "-1", "2", "9.0", "R[PEPTMIDE]K", "{Decoy_P1}"},
        {"sample1.3.1", "-1", "3", "8.0", "R[ANOTHER]K", "{Decoy_P2}"}});
    aerith::Config collision_config;
    collision_config.target_pins = {target0.string(), target1.string()};
    collision_config.decoy_pins = {decoy0.string(), decoy1.string()};
    const auto collision_data = aerith::PinReader::read(collision_config);
    assert(collision_data.rows.size() == 3);
    assert(collision_data.removed_decoy_peptide_collisions == 3);
    assert(std::count_if(
        collision_data.rows.begin(), collision_data.rows.end(),
        [](const auto& row) { return row.label == 1; }) == 2);
    assert(std::count_if(
        collision_data.rows.begin(), collision_data.rows.end(),
        [](const auto& row) { return row.label == -1; }) == 1);
    assert(std::any_of(
        collision_data.rows.begin(), collision_data.rows.end(),
        [](const auto& row) {
            return row.label == -1 &&
                aerith::stripped_peptide(row.peptide) == "UNIQUE";
        }));
    std::filesystem::remove_all(collision_root);

    const auto root = std::filesystem::temp_directory_path() /
        "aerith_native_protein_unit";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root / "sample");
    const auto database = root / "target.faa";
    const auto decoy_database = root / "decoy.faa";
    {
        std::ofstream fasta(database);
        fasta << ">sp|P1|ONE Protein one OS=Test GN=one PE=1\n"
              << "MAAAAAAABBBBBBBK\n"
              << ">sp|P2|TWO Protein two OS=Test GN=two PE=4\n"
              << "MBBBBBBBCCCCCCCK\n";
    }
    {
        std::ofstream fasta(decoy_database);
        fasta << ">sp|D1|DECOY Decoy protein OS=Test GN=decoy PE=4\n"
              << "MKKKKKKKLLLLLLLK\n";
    }
    aerith::Dataset data;
    data.input_paths.push_back("sample.pin");
    const auto add_psm = [&](const std::string& id, const std::string& peptide,
                             const std::string& proteins) {
        aerith::Psm psm;
        psm.id = id;
        psm.peptide = peptide;
        psm.proteins = proteins;
        psm.label = 1;
        psm.charge = 2;
        psm.exp_mass = 1000.0;
        psm.log10_precursor_intensity = 3.0;
        data.rows.push_back(std::move(psm));
    };
    add_psm("sample.1.1", "K[AAAAAAA]R", "{sp|P1|ONE}");
    add_psm("sample.2.1", "K[BBBBBBB]R", "{sp|P1|ONE,sp|P2|TWO}");
    add_psm("sample.3.1", "K[CCCCCCC]R", "{sp|P2|TWO}");
    aerith::TransferredIon transfer;
    transfer.psm = data.rows.front();
    transfer.psm.charge = 3;
    transfer.psm.quantified_intensity = 500.0;
    transfer.psm.has_chromatographic_feature = true;
    transfer.psm.apex_retention = 60.0;
    transfer.psm.retention_start = 55.0;
    transfer.psm.retention_end = 65.0;
    transfer.psm.retention_fwhm = 5.0;
    transfer.psm.apex_scan = 10;
    transfer.psm.traced_scans = 3;
    transfer.qvalue = 0.005;
    data.transferred_ions.push_back(std::move(transfer));
    aerith::Config config;
    config.database_path = database.string();
    config.decoy_database_path = decoy_database.string();
    config.protein_output_dir = root.string();
    config.output_prefixes.push_back((root / "sample" / "sample").string());
    const std::vector<double> protein_q(data.rows.size(), 0.001);
    const std::vector<double> protein_pep(data.rows.size(), 0.001);
    const std::vector<double> protein_scores(data.rows.size(), 1.0);
    const auto assembly =
        aerith::ProteinAssembler::write(
            config, data, protein_scores, protein_q, protein_pep);
    assert(assembly.proteins == 2);
    assert(assembly.stages.size() == 5);
    assert(assembly.stages.front().name ==
           "Read and annotate target/decoy FASTAs");
    assert(assembly.stages.back().name ==
           "Build and write combined/sample reports");
    assert(std::none_of(
        assembly.stages.begin(), assembly.stages.end(),
        [](const auto& stage) { return stage.uses_simd; }));
    assert(std::filesystem::exists(root / "combined_protein.tsv"));
    assert(std::filesystem::exists(root / "combined_protein.fas"));
    assert(std::filesystem::exists(root / "combined_psm.tsv"));
    assert(std::filesystem::exists(root / "combined_ion.tsv"));
    assert(std::filesystem::exists(root / "combined_modified_peptide.tsv"));
    assert(std::filesystem::exists(root / "combined_peptide.tsv"));
    assert(std::filesystem::exists(root / "sample" / "protein.tsv"));
    assert(std::filesystem::exists(root / "sample" / "psm.tsv"));
    assert(std::filesystem::exists(root / "sample" / "ion.tsv"));
    assert(std::filesystem::exists(root / "sample" / "peptide.tsv"));
    assert(std::filesystem::exists(root / "sample" / "modified_peptide.tsv"));
    assert(!std::filesystem::exists(root / "sample" / "sample.pep.xml"));
    assert(!std::filesystem::exists(root / "aerith.log"));
    std::ifstream protein_table(root / "sample" / "protein.tsv");
    const std::string report(
        (std::istreambuf_iterator<char>(protein_table)),
        std::istreambuf_iterator<char>());
    assert(report.find("sp|P1|ONE") != std::string::npos);
    assert(report.find("sp|P2|TWO") != std::string::npos);
    assert(report.find("\t2500\t1500\t2500\t") != std::string::npos);
    std::ifstream combined_protein_table(root / "combined_protein.tsv");
    const std::string combined_protein_report(
        (std::istreambuf_iterator<char>(combined_protein_table)),
        std::istreambuf_iterator<char>());
    assert(combined_protein_report.find("\tProtein Coverage\t") !=
           std::string::npos);
    std::ifstream psm_table(root / "combined_psm.tsv");
    const std::string psm_report(
        (std::istreambuf_iterator<char>(psm_table)),
        std::istreambuf_iterator<char>());
    assert(psm_report.find("\tIntensity\tAssigned Modifications") !=
           std::string::npos);
    assert(psm_report.find("1C(57.0215)") != std::string::npos);
    assert(psm_report.find("\tTarget\t0.00\tfalse\tfalse\t") !=
           std::string::npos);
    const auto header_end = psm_report.find('\n');
    const auto row_end = psm_report.find('\n', header_end + 1);
    assert(header_end != std::string::npos);
    assert(row_end != std::string::npos);
    const auto psm_header = split_tabs(psm_report.substr(0, header_end));
    const auto psm_row = split_tabs(
        psm_report.substr(header_end + 1, row_end - header_end - 1));
    const auto calibrated_mass = std::find(
        psm_header.begin(), psm_header.end(), "Calibrated Observed Mass");
    const auto calibrated_mz = std::find(
        psm_header.begin(), psm_header.end(), "Calibrated Observed M/Z");
    assert(calibrated_mass != psm_header.end());
    assert(calibrated_mz != psm_header.end());
    assert(psm_row[calibrated_mass - psm_header.begin()].empty());
    assert(psm_row[calibrated_mz - psm_header.begin()].empty());
    std::ifstream ion_table(root / "sample" / "ion.tsv");
    const std::string ion_report(
        (std::istreambuf_iterator<char>(ion_table)),
        std::istreambuf_iterator<char>());
    const auto ion_header_end = ion_report.find('\n');
    const auto ion_header = split_tabs(
        ion_report.substr(0, ion_header_end));
    const auto probability = std::find(
        ion_header.begin(), ion_header.end(), "Probability");
    const auto expectation = std::find(
        ion_header.begin(), ion_header.end(), "Expectation");
    const auto match_type = std::find(
        ion_header.begin(), ion_header.end(), "Match Type");
    assert(probability != ion_header.end());
    assert(expectation != ion_header.end());
    assert(match_type != ion_header.end());
    bool found_mbr = false;
    std::size_t line_begin = ion_header_end + 1;
    while (line_begin < ion_report.size()) {
        const auto line_end = ion_report.find('\n', line_begin);
        const auto fields = split_tabs(ion_report.substr(
            line_begin, line_end - line_begin));
        if (fields[match_type - ion_header.begin()] == "MBR") {
            assert(fields[probability - ion_header.begin()].empty());
            assert(fields[expectation - ion_header.begin()].empty());
            found_mbr = true;
        }
        if (line_end == std::string::npos) break;
        line_begin = line_end + 1;
    }
    assert(found_mbr);
    std::ifstream combined_peptide_table(root / "combined_peptide.tsv");
    const std::string combined_peptide_report(
        (std::istreambuf_iterator<char>(combined_peptide_table)),
        std::istreambuf_iterator<char>());
    assert(combined_peptide_report.find("sample Spectral Count") !=
           std::string::npos);
    assert(combined_peptide_report.find("sample MaxLFQ Intensity") !=
           std::string::npos);
    aerith::Summary summary;
    summary.reporting_fdr = 0.01;
    summary.target_ids = 3;
    summary.distinct_target_peptides = 3;
    summary.distinct_target_peptide_forms = 3;
    summary.distinct_target_ptm_peptides = 1;
    summary.target_ptm_psms = 1;
    summary.removed_decoy_peptide_collisions = 154;
    summary.quantification_stages.push_back({
        "Trace identified XICs + detect peaks/intensity",
        {2.0, 8.0}, true, false});
    summary.protein_assembly_stages = assembly.stages;
    std::ostringstream log;
    aerith::print_summary(log, summary);
    const auto log_text = log.str();
    assert(log_text.find("Distinct naked peptides") != std::string::npos);
    assert(log_text.find("Distinct PTM peptide forms") != std::string::npos);
    assert(log_text.find("Removed colliding decoy PSMs  154") !=
           std::string::npos);
    assert(log_text.find("Timing by stage (seconds)") != std::string::npos);
    assert(log_text.find("Quantification total") !=
           std::string::npos);
    assert(log_text.find("Trace identified XICs + detect peaks/intensity") !=
           std::string::npos);
    assert(log_text.find("Protein assembly: Read and annotate") !=
           std::string::npos);
    assert(log_text.find("Protein assembly optimization detail") ==
           std::string::npos);
    assert(log_text.find("aerith.log") == std::string::npos);
    std::filesystem::remove_all(root);

    std::cout << "aerith unit tests passed\n";
    return 0;
}
