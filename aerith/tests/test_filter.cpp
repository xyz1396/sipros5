#include "filter.hpp"
#include "isotope.hpp"
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
#ifdef AERITH_TEST_WITH_TORCH
    aerith::Psm rt_psm;
    rt_psm.peptide = "R[QMDVVEQMMPGLK]D";
    rt_psm.charge = 2;
    const auto unmodified_rt_tokens =
        aerith::diann_rt_tokens_for_testing(rt_psm);
    assert(unmodified_rt_tokens ==
           (std::vector<std::int64_t>{
               1, 16, 11, 19, 5, 5, 17, 16, 11, 11, 8, 3, 7, 20, 2}));
    rt_psm.peptide = "R[QM~DVVEQMMPGLK]D";
    const auto oxidized_rt_tokens =
        aerith::diann_rt_tokens_for_testing(rt_psm);
    assert(oxidized_rt_tokens[2] == 26);
    assert(oxidized_rt_tokens[8] == 11);
    assert(oxidized_rt_tokens[9] == 11);
#endif
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
        pin << "SpecId\tLabel\tScanNr\tretentiontime\tExpMass\tObservedMass\tWDPscores\tdiffScores"
               "\tranks\tparentCharges\tPeptide\tProteins\n";
        for (const auto& row : rows) {
            pin << row[0] << '\t' << row[1] << '\t' << row[2]
                << "\t1.0\t1000.0\t999.5\t" << row[3] << "\t0\t1\t2\t"
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
    aerith::Config sip_result_config;
    sip_result_config.output_prefixes = {
        (collision_root / "SIP").string()};
    sip_result_config.filtered_only = false;
    sip_result_config.sip_isotope = "C13";
    const std::vector<double> collision_scores(
        collision_data.rows.size(), 1.0);
    std::vector<double> collision_q(collision_data.rows.size(), 0.0);
    const auto filtered_target = std::find_if(
        collision_data.rows.begin(), collision_data.rows.end(),
        [](const auto& row) { return row.label == 1; });
    assert(filtered_target != collision_data.rows.end());
    collision_q[static_cast<std::size_t>(
        filtered_target - collision_data.rows.begin())] = 0.02;
    const std::vector<double> collision_pep(
        collision_data.rows.size(), 0.0);
    aerith::ResultWriter::write(
        sip_result_config, collision_data, collision_scores,
        collision_q, collision_pep, {}, {});
    std::ifstream sip_result(collision_root / "SIP_filtered_psms.tsv");
    std::string sip_header_line;
    assert(static_cast<bool>(std::getline(
        sip_result, sip_header_line)));
    const auto sip_header = split_tabs(sip_header_line);
    assert(std::find(sip_header.begin(), sip_header.end(), "SVMscore") !=
           sip_header.end());
    for (const std::string& removed : {"score", "Label", "diffScores"}) {
        assert(std::find(sip_header.begin(), sip_header.end(), removed) ==
               sip_header.end());
    }
    std::size_t sip_rows = 0;
    std::string sip_line;
    while (std::getline(sip_result, sip_line)) {
        if (!sip_line.empty()) ++sip_rows;
    }
    assert(sip_rows == 1);
    sip_result.close();
    const auto check_score_table = [](
        const std::filesystem::path& path, std::size_t expected_rows) {
        std::ifstream table(path);
        std::string header_line;
        assert(static_cast<bool>(std::getline(table, header_line)));
        assert(split_tabs(header_line) == std::vector<std::string>({
            "PSMId", "SVMscore", "q-value", "posterior_error_prob",
            "peptide", "proteinIds"}));
        std::size_t rows = 0;
        std::string line;
        while (std::getline(table, line)) {
            if (!line.empty()) ++rows;
        }
        assert(rows == expected_rows);
    };
    check_score_table(collision_root / "SIP_target_psms.tsv", 2);
    check_score_table(collision_root / "SIP_decoy_psms.tsv", 1);
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
    data.input_paths.push_back("negative_control.pin");
    const auto add_psm = [&](const std::string& id, const std::string& peptide,
                             const std::string& proteins) {
        aerith::Psm psm;
        psm.id = id;
        psm.peptide = peptide;
        psm.proteins = proteins;
        psm.label = 1;
        psm.charge = 2;
        psm.exp_mass = 1000.0;
        psm.observed_mass = 999.5;
        psm.log10_precursor_intensity = 3.0;
        psm.ms1_isotopic_abundance = 49.5;
        psm.ms2_isotopic_abundance = 50.0;
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
    config.output_prefixes.push_back(
        (root / "negative_control" / "negative_control").string());
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
    assert(std::filesystem::exists(
        root / "combined_protein_with_PSM.tsv"));
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
    const auto protein_header_end = report.find('\n');
    const auto protein_row_end = report.find('\n', protein_header_end + 1);
    assert(protein_header_end != std::string::npos);
    assert(protein_row_end != std::string::npos);
    const auto protein_header = split_tabs(
        report.substr(0, protein_header_end));
    const auto protein_row = split_tabs(report.substr(
        protein_header_end + 1,
        protein_row_end - protein_header_end - 1));
    assert(protein_header.size() == protein_row.size());
    assert(std::find(
               protein_header.begin(), protein_header.end(),
               "Razor Observed Modifications") == protein_header.end());
    std::ifstream combined_protein_table(root / "combined_protein.tsv");
    const std::string combined_protein_report(
        (std::istreambuf_iterator<char>(combined_protein_table)),
        std::istreambuf_iterator<char>());
    assert(combined_protein_report.find("\tProtein Coverage\t") !=
           std::string::npos);
    const auto report_row = [&](const std::string& report,
                                const std::string& protein) {
        std::istringstream rows(report);
        std::string row;
        std::getline(rows, row);
        while (std::getline(rows, row)) {
            const auto fields = split_tabs(row);
            if (!fields.empty() && fields.front() == protein) return fields;
        }
        return std::vector<std::string>{};
    };
    const auto combined_header_end = combined_protein_report.find('\n');
    const auto combined_header = split_tabs(
        combined_protein_report.substr(0, combined_header_end));
    const auto combined_column = [&](const std::string& name) {
        const auto found =
            std::find(combined_header.begin(), combined_header.end(), name);
        assert(found != combined_header.end());
        return static_cast<std::size_t>(found - combined_header.begin());
    };
    const auto p1_combined =
        report_row(combined_protein_report, "sp|P1|ONE");
    const auto p2_combined =
        report_row(combined_protein_report, "sp|P2|TWO");
    assert(!p1_combined.empty());
    assert(!p2_combined.empty());
    assert(p1_combined[combined_column("sample Spectral Count")] == "2");
    assert(p1_combined[combined_column("sample Unique Spectral Count")] == "1");
    assert(p1_combined[combined_column("sample Total Spectral Count")] == "2");
    assert(p1_combined[combined_column("sample Intensity")] == "2500");
    assert(p2_combined[combined_column("sample Spectral Count")] == "1");
    assert(p2_combined[combined_column("sample Unique Spectral Count")] == "1");
    assert(p2_combined[combined_column("sample Total Spectral Count")] == "2");
    assert(p2_combined[combined_column("sample Intensity")] == "1000");
    std::ifstream combined_psm_protein_table(
        root / "combined_protein_with_PSM.tsv");
    const std::string combined_psm_protein_report(
        (std::istreambuf_iterator<char>(combined_psm_protein_table)),
        std::istreambuf_iterator<char>());
    assert(combined_psm_protein_report.find(
        "sample_PeptideSequences") != std::string::npos);
    assert(combined_psm_protein_report.find(
        "sample_PSMids") != std::string::npos);
    assert(combined_psm_protein_report.find(
        "sample_MS1IsotopicAbundances") != std::string::npos);
    assert(combined_psm_protein_report.find(
        "sample_MS2IsotopicAbundances") != std::string::npos);
    assert(combined_psm_protein_report.find(
        "sample_SIP_intensity") != std::string::npos);
    assert(combined_psm_protein_report.find(
        "_log10_precursorIntensities") == std::string::npos);
    assert(combined_psm_protein_report.find(
        "_ProteinAbundance") == std::string::npos);
    {
        std::istringstream rows(combined_psm_protein_report);
        std::string row;
        std::getline(rows, row);
        const auto expected_tabs =
            static_cast<std::size_t>(std::count(row.begin(), row.end(), '\t'));
        while (std::getline(rows, row)) {
            assert(static_cast<std::size_t>(
                std::count(row.begin(), row.end(), '\t')) == expected_tabs);
        }
    }
    const auto annotated_header_end =
        combined_psm_protein_report.find('\n');
    const auto annotated_header = split_tabs(
        combined_psm_protein_report.substr(0, annotated_header_end));
    const auto peptide_column = std::find(
        annotated_header.begin(), annotated_header.end(),
        "sample_PeptideSequences");
    const auto psm_id_column = std::find(
        annotated_header.begin(), annotated_header.end(),
        "sample_PSMids");
    assert(peptide_column != annotated_header.end());
    assert(psm_id_column != annotated_header.end());
    const auto peptide_index =
        static_cast<std::size_t>(peptide_column - annotated_header.begin());
    const auto psm_id_index =
        static_cast<std::size_t>(psm_id_column - annotated_header.begin());
    const auto p1_annotated =
        report_row(combined_psm_protein_report, "sp|P1|ONE");
    const auto p2_annotated =
        report_row(combined_psm_protein_report, "sp|P2|TWO");
    assert(p1_annotated[peptide_index] == "AAAAAAA,BBBBBBB,AAAAAAA");
    assert(p2_annotated[peptide_index] == "CCCCCCC");
    assert(p1_annotated[psm_id_index] ==
           "sample.1.1,sample.2.1,AAAAAAA_MBR");
    assert(p2_annotated[psm_id_index] == "sample.3.1");
    const auto ms1_column = std::find(
        annotated_header.begin(), annotated_header.end(),
        "sample_MS1IsotopicAbundances");
    const auto ms2_column = std::find(
        annotated_header.begin(), annotated_header.end(),
        "sample_MS2IsotopicAbundances");
    const auto sip_intensity_column = std::find(
        annotated_header.begin(), annotated_header.end(),
        "sample_SIP_intensity");
    assert(ms1_column != annotated_header.end());
    assert(ms2_column != annotated_header.end());
    assert(sip_intensity_column != annotated_header.end());
    const auto ms1_index =
        static_cast<std::size_t>(ms1_column - annotated_header.begin());
    const auto ms2_index =
        static_cast<std::size_t>(ms2_column - annotated_header.begin());
    const auto sip_intensity_index = static_cast<std::size_t>(
        sip_intensity_column - annotated_header.begin());
    assert(p1_annotated[ms1_index] == "49.5,49.5,49.5");
    assert(p1_annotated[ms2_index] == "50,50,50");
    assert(p1_annotated[sip_intensity_index] == "0,0,500");
    std::ifstream psm_table(root / "combined_psm.tsv");
    const std::string psm_report(
        (std::istreambuf_iterator<char>(psm_table)),
        std::istreambuf_iterator<char>());
    assert(psm_report.find("\tIntensity\tAssigned Modifications") !=
           std::string::npos);
    assert(psm_report.find("1C(57.0215)") != std::string::npos);
    assert(psm_report.find("\tTarget\t") !=
           std::string::npos);
    const auto header_end = psm_report.find('\n');
    const auto row_end = psm_report.find('\n', header_end + 1);
    assert(header_end != std::string::npos);
    assert(row_end != std::string::npos);
    const auto psm_header = split_tabs(psm_report.substr(0, header_end));
    const auto psm_row = split_tabs(
        psm_report.substr(header_end + 1, row_end - header_end - 1));
    assert(psm_header.size() == psm_row.size());
    const auto observed_mass = std::find(
        psm_header.begin(), psm_header.end(), "Observed Mass");
    const auto observed_mz = std::find(
        psm_header.begin(), psm_header.end(), "Observed M/Z");
    assert(observed_mass != psm_header.end());
    assert(observed_mz != psm_header.end());
    assert(psm_row[observed_mass - psm_header.begin()] == "999.5000");
    assert(psm_row[observed_mz - psm_header.begin()] == "500.7573");
    for (const std::string& removed : {
             "Calibrated Observed Mass", "Calibrated Observed M/Z",
             "SpectralSim", "RTScore", "Expectation",
             "Hyperscore", "Nextscore", "Observed Modifications",
             "Is Decoy", "Is Contaminant", "Purity"}) {
        assert(std::find(psm_header.begin(), psm_header.end(), removed) ==
               psm_header.end());
    }
    assert(std::find(psm_header.begin(), psm_header.end(), "SVMscore") !=
           psm_header.end());
    std::ifstream peptide_table(root / "sample" / "peptide.tsv");
    std::string peptide_header_line;
    std::string peptide_row_line;
    assert(static_cast<bool>(std::getline(
        peptide_table, peptide_header_line)));
    assert(static_cast<bool>(std::getline(
        peptide_table, peptide_row_line)));
    const auto peptide_header = split_tabs(peptide_header_line);
    const auto peptide_row = split_tabs(peptide_row_line);
    assert(peptide_header.size() == peptide_row.size());
    assert(std::find(
               peptide_header.begin(), peptide_header.end(),
               "Observed Modifications") == peptide_header.end());
    std::ifstream ion_table(root / "sample" / "ion.tsv");
    const std::string ion_report(
        (std::istreambuf_iterator<char>(ion_table)),
        std::istreambuf_iterator<char>());
    const auto ion_header_end = ion_report.find('\n');
    const auto ion_header = split_tabs(
        ion_report.substr(0, ion_header_end));
    const auto probability = std::find(
        ion_header.begin(), ion_header.end(), "Probability");
    const auto match_type = std::find(
        ion_header.begin(), ion_header.end(), "Match Type");
    assert(probability != ion_header.end());
    assert(std::find(
               ion_header.begin(), ion_header.end(), "Expectation") ==
           ion_header.end());
    assert(std::find(
               ion_header.begin(), ion_header.end(),
               "Observed Modifications") == ion_header.end());
    assert(match_type != ion_header.end());
    bool found_mbr = false;
    std::size_t line_begin = ion_header_end + 1;
    while (line_begin < ion_report.size()) {
        const auto line_end = ion_report.find('\n', line_begin);
        const auto fields = split_tabs(ion_report.substr(
            line_begin, line_end - line_begin));
        assert(fields.size() == ion_header.size());
        if (fields[match_type - ion_header.begin()] == "MBR") {
            assert(fields[probability - ion_header.begin()].empty());
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
    summary.files = 2;
    summary.psms = 13;
    summary.targets = 7;
    summary.decoys = 6;
    summary.threads = 4;
    summary.score_model = "test_model";
    summary.reporting_fdr = 0.01;
    summary.target_ids = 3;
    summary.distinct_target_peptides = 3;
    summary.distinct_target_peptide_forms = 3;
    summary.distinct_target_ptm_peptides = 1;
    summary.target_ptm_psms = 1;
    summary.mbr_ions = 1;
    summary.protein_ids = 2;
    summary.removed_decoy_peptide_collisions = 154;
    summary.negative_control_input_psms = 12;
    summary.negative_control_threshold_filtered_psms = 2;
    summary.negative_control_candidates = 10;
    summary.negative_control_targets = 6;
    summary.negative_control_decoys = 4;
    summary.negative_control_target_ids = 5;
    summary.negative_control_label_threshold = 2.0;
    summary.negative_control_target_output_path =
        "SIP_target_psms.tsv";
    summary.negative_control_decoy_output_path =
        "SIP_decoy_psms.tsv";
    summary.negative_control_output_path =
        "SIP_filtered_psms.tsv";
    summary.negative_control_timing = {1.0, 2.0};
    summary.negative_control_stages.push_back({
        "Fit and score SIP-Negative-control SVM folds",
        {0.5, 1.5}, true, false});
    summary.negative_control_feature_names = {"WDPscores"};
    summary.negative_control_model.name = "SIP-Negative-control";
    summary.negative_control_model.feature_weights = {{
        {0.5, -0.1}, {0.4, -0.2}, {0.6, -0.3}}};
    summary.quantification_stages.push_back({
        "Trace identified XICs + detect peaks/intensity",
        {2.0, 8.0}, true, false});
    summary.protein_assembly_stages = assembly.stages;
    std::ostringstream log;
    aerith::print_summary(log, summary);
    const auto log_text = log.str();
    assert(log_text.find("Distinct naked peptides") != std::string::npos);
    assert(log_text.find("Distinct PTM peptide forms") != std::string::npos);
    assert(log_text.find("Removed colliding decoy PSMs") !=
           std::string::npos);
    assert(log_text.find("SIP-Negative-control filtering total") !=
           std::string::npos);
    assert(log_text.find("SIP_target_psms.tsv") != std::string::npos);
    assert(log_text.find("SIP_decoy_psms.tsv") != std::string::npos);
    assert(log_text.find("SIP_filtered_psms.tsv") != std::string::npos);
    assert(log_text.find("SIP-Negative-control threshold") !=
           std::string::npos);
    assert(log_text.find("Primary-filtered target PSMs") !=
           std::string::npos);
    assert(log_text.find("PSMs below SIP threshold") !=
           std::string::npos);
    const auto log_line = [&](const std::string& label) {
        const auto begin = log_text.find(label);
        assert(begin != std::string::npos);
        const auto line_begin = log_text.rfind('\n', begin);
        const auto line_end = log_text.find('\n', begin);
        return log_text.substr(
            line_begin == std::string::npos ? 0 : line_begin + 1,
            line_end - (line_begin == std::string::npos
                ? 0 : line_begin + 1));
    };
    const auto assert_value_column = [&](const std::string& label,
                                         const std::string& value) {
        const auto line = log_line(label);
        const auto column = line.rfind(value);
        assert(column == 34);
    };
    assert_value_column("Files", "2");
    assert_value_column("Input PSMs", "13");
    assert_value_column("Targets", "7");
    assert_value_column("Decoys", "6");
    assert_value_column("Removed colliding decoy PSMs", "154");
    assert_value_column("OpenMP threads", "4");
    assert_value_column("Score model", "test_model");
    assert_value_column("Reporting FDR", "1%");
    assert_value_column("Target PSMs", "3");
    assert_value_column("Distinct naked peptides", "3");
    assert_value_column("Distinct peptide forms", "3");
    assert_value_column("Distinct PTM peptide forms", "1");
    assert_value_column("PSMs carrying PTMs", "1");
    assert_value_column("MBR transferred ions", "1");
    assert_value_column("Protein IDs at 1% FDR", "2");
    assert_value_column("SIP-Negative-control threshold", "2%");
    assert_value_column("Primary-filtered target PSMs", "12");
    assert_value_column("PSMs below SIP threshold", "2");
    assert_value_column("SIP-Negative-control candidates", "10");
    assert_value_column("SIP-Negative-control targets", "6");
    assert_value_column("SIP-Negative-control decoys", "4");
    assert_value_column("SIP-Negative-control target IDs", "5");
    assert_value_column("RT feature source", "DIA-NN prediction");
    assert_value_column("Aerith internal RT model", "skipped");
    assert(log_text.find("Fit and score SIP-Negative-control SVM folds") !=
           std::string::npos);
    assert(log_text.find("Recompute negative-control RT features") ==
           std::string::npos);
    assert(log_text.find("SIP-Negative-control SVM feature weights") !=
           std::string::npos);
    assert(log_text.find("Sample: SIP-Negative-control") ==
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
    aerith::Dataset sip_data;
    sip_data.rows = data.rows;
    for (auto& psm : sip_data.rows) psm.sample_name = "sample";
    sip_data.rows[0].quantification_attempted = true;
    sip_data.rows[0].has_chromatographic_feature = true;
    sip_data.rows[0].quantified_intensity = 123.0;
    sip_data.rows[1].quantification_attempted = true;
    sip_data.rows[1].has_chromatographic_feature = true;
    sip_data.rows[1].quantified_intensity = 999.0;
    sip_data.rows[2].quantification_attempted = true;
    sip_data.rows[2].has_chromatographic_feature = true;
    sip_data.rows[2].quantified_intensity = 456.0;
    auto sip_mapping_q = protein_q;
    sip_mapping_q[1] = 0.02;
    aerith::Config sip_mapping_config;
    sip_mapping_config.protein_reference_path =
        (root / "combined_protein.tsv").string();
    sip_mapping_config.sip_protein_output_path =
        (root / "combined_protein_with_SIP_filtered_PSM.tsv").string();
    sip_mapping_config.negative_control_samples = {"negative_control"};
    aerith::ProteinAssembler::write_sip_psm_mapping(
        sip_mapping_config, sip_data, sip_mapping_q);
    std::ifstream sip_mapping_table(
        root / "combined_protein_with_SIP_filtered_PSM.tsv");
    const std::string sip_mapping_report(
        (std::istreambuf_iterator<char>(sip_mapping_table)),
        std::istreambuf_iterator<char>());
    assert(sip_mapping_report.find(
        "SIP_sample_PSMids") != std::string::npos);
    assert(sip_mapping_report.find(
        "SIP_sample_PeptideSequences") != std::string::npos);
    assert(sip_mapping_report.find(
        "SIP_sample_MS1IsotopicAbundances") != std::string::npos);
    assert(sip_mapping_report.find(
        "SIP_sample_MS2IsotopicAbundances") != std::string::npos);
    assert(sip_mapping_report.find(
        "SIP_sample_SIP_intensity") != std::string::npos);
    assert(sip_mapping_report.find(
        "_log10_precursorIntensities") == std::string::npos);
    assert(sip_mapping_report.find(
        "negative_control Spectral Count") == std::string::npos);
    assert(sip_mapping_report.find(
        "SIP_negative_control_") == std::string::npos);
    {
        std::istringstream rows(sip_mapping_report);
        std::string row;
        std::getline(rows, row);
        const auto expected_tabs =
            static_cast<std::size_t>(std::count(row.begin(), row.end(), '\t'));
        while (std::getline(rows, row)) {
            assert(static_cast<std::size_t>(
                std::count(row.begin(), row.end(), '\t')) == expected_tabs);
        }
    }
    const auto sip_header_end = sip_mapping_report.find('\n');
    const auto sip_mapping_header = split_tabs(
        sip_mapping_report.substr(0, sip_header_end));
    const auto sip_column = [&](const std::string& name) {
        const auto found =
            std::find(
                sip_mapping_header.begin(), sip_mapping_header.end(), name);
        assert(found != sip_mapping_header.end());
        return static_cast<std::size_t>(
            found - sip_mapping_header.begin());
    };
    const auto sip_p1 =
        report_row(sip_mapping_report, "sp|P1|ONE");
    const auto sip_p2 =
        report_row(sip_mapping_report, "sp|P2|TWO");
    assert(sip_p1[sip_column("Combined Total Peptides")] == "1");
    assert(sip_p1[sip_column("Combined Spectral Count")] == "1");
    assert(sip_p1[sip_column("Combined Unique Spectral Count")] == "1");
    assert(sip_p1[sip_column("Combined Total Spectral Count")] == "1");
    assert(sip_p1[sip_column("sample Spectral Count")] == "1");
    assert(sip_p1[sip_column("sample Unique Spectral Count")] == "1");
    assert(sip_p1[sip_column("sample Total Spectral Count")] == "1");
    assert(sip_p1[sip_column("sample Intensity")] == "123");
    assert(sip_p1[sip_column("sample MaxLFQ Intensity")] == "123");
    assert(sip_p1[sip_column("SIP_sample_PSMids")] == "sample.1.1");
    assert(sip_p1[sip_column("SIP_sample_PeptideSequences")] == "AAAAAAA");
    assert(sip_p1[sip_column("SIP_sample_MS1IsotopicAbundances")] == "49.5");
    assert(sip_p1[sip_column("SIP_sample_MS2IsotopicAbundances")] == "50");
    assert(sip_p1[sip_column("SIP_sample_SIP_intensity")] == "123");
    assert(sip_p2[sip_column("Combined Total Peptides")] == "1");
    assert(sip_p2[sip_column("Combined Spectral Count")] == "1");
    assert(sip_p2[sip_column("Combined Unique Spectral Count")] == "1");
    assert(sip_p2[sip_column("Combined Total Spectral Count")] == "1");
    assert(sip_p2[sip_column("sample Intensity")] == "456");
    assert(sip_p2[sip_column("sample MaxLFQ Intensity")] == "456");
    assert(sip_p2[sip_column("SIP_sample_PSMids")] == "sample.3.1");
    assert(sip_p2[sip_column("SIP_sample_PeptideSequences")] == "CCCCCCC");
    assert(sip_mapping_report.find("sample.2.1") == std::string::npos);

    aerith::Config sip_isotope_config;
    sip_isotope_config.sip_isotope = "C13";
    aerith::initialize_sip_isotope_model(sip_isotope_config);
    const auto precursor_peaks =
        aerith::precursor_isotope_peaks("K[AAAAAAA]I", 50.0, 6);
    assert(precursor_peaks.size() == 6);
    assert(std::is_sorted(
        precursor_peaks.begin(), precursor_peaks.end(),
        [](const auto& left, const auto& right) {
            return left.probability > right.probability;
        }));
    const auto product_envelopes =
        aerith::product_isotope_envelopes("K[AAAAAAA]I", 50.0);
    assert(product_envelopes.b_mass.size() == 6);
    assert(product_envelopes.y_mass.size() == 6);

    protein_table.close();
    combined_protein_table.close();
    combined_psm_protein_table.close();
    psm_table.close();
    peptide_table.close();
    ion_table.close();
    combined_peptide_table.close();
    sip_mapping_table.close();
    std::filesystem::remove_all(root);

    std::cout << "aerith unit tests passed\n";
    return 0;
}
