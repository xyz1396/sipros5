#include "filter.hpp"
#include "pipeline.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
    assert(aerith::stripped_peptide("K[PEPTM~IDE]R") == "PEPTMIDE");
    assert(aerith::stripped_peptide("-.ACD[+57.0]EF.-") == "ACDEF");
    assert(aerith::stripped_peptide("[%M~PEPTIDE]R") == "MPEPTIDE");
    assert(aerith::stripped_peptide("K[PEPN!IDE]") == "PEPNIDE");
    assert(aerith::stripped_peptide("[PEPTIDE]R") == "PEPTIDE");

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
    assert(std::filesystem::exists(root / "protein.tsv"));
    assert(std::filesystem::exists(root / "protein.fas"));
    assert(std::filesystem::exists(root / "psm.tsv"));
    assert(std::filesystem::exists(root / "sample" / "protein.tsv"));
    assert(std::filesystem::exists(root / "sample" / "psm.tsv"));
    assert(!std::filesystem::exists(root / "sample" / "sample.pep.xml"));
    std::ifstream protein_table(root / "protein.tsv");
    const std::string report(
        (std::istreambuf_iterator<char>(protein_table)),
        std::istreambuf_iterator<char>());
    assert(report.find("sp|P1|ONE") != std::string::npos);
    assert(report.find("sp|P2|TWO") != std::string::npos);
    assert(report.find("\t2000\t1000\t2000\t") != std::string::npos);
    std::ifstream psm_table(root / "psm.tsv");
    const std::string psm_report(
        (std::istreambuf_iterator<char>(psm_table)),
        std::istreambuf_iterator<char>());
    assert(psm_report.find("\tIntensity\tAssigned Modifications") !=
           std::string::npos);
    assert(psm_report.find("1C(57.0215)") != std::string::npos);
    assert(psm_report.find("\tTarget\t0.00\tfalse\tfalse\t") !=
           std::string::npos);
    aerith::Summary summary;
    summary.reporting_fdr = 0.01;
    summary.target_ids = 3;
    summary.distinct_target_peptides = 3;
    summary.distinct_target_peptide_forms = 3;
    summary.distinct_target_ptm_peptides = 1;
    summary.target_ptm_psms = 1;
    summary.protein_assembly_stages = assembly.stages;
    std::ostringstream log;
    aerith::print_summary(log, summary);
    const auto log_text = log.str();
    assert(log_text.find("Distinct naked peptides") != std::string::npos);
    assert(log_text.find("Distinct PTM peptide forms") != std::string::npos);
    assert(log_text.find("Timing by stage (seconds)") != std::string::npos);
    assert(log_text.find("Protein assembly: Read and annotate") !=
           std::string::npos);
    assert(log_text.find("Protein assembly optimization detail") ==
           std::string::npos);
    std::filesystem::remove_all(root);

    std::cout << "aerith unit tests passed\n";
    return 0;
}
