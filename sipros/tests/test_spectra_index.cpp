#include "spectraindex.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#ifdef _WIN32
#include "windows_posix_compat.h"
#else
#include <unistd.h>
#endif
#include <vector>



bool expect(bool condition, const std::string &message)
{
	if (!condition)
		std::cerr << message << '\n';
	return condition;
}


int main()
{
	const std::filesystem::path path = std::filesystem::temp_directory_path() /
		("sipros_spectra_index_test_" + std::to_string(getpid()) + ".sfi");
	std::error_code ec;
	std::filesystem::remove(path, ec);

	sipros::SpectraIndexMetadata metadata;
	metadata.chemistryProfileId = "builtin-sip-fixed-cam-v2";
	metadata.recordKind = "target";
	metadata.targetSipAbundancePct =
		sipros::SpectraIndexMetadata::MixedSipAbundancePct;
	metadata.sipAtom = 'C';
	metadata.sipIsotopeMassNumber = 13;
	metadata.probabilityCutoff = 0.01;
	metadata.envelopeTopN = 3;
	metadata.label = 1;

	std::vector<sipros::SpectraIndexRecordInput> records(2);
	records[0].psmId = "scan.10.10.2";
	records[0].peptide = "[PEPTIDE]";
	records[0].proteins = "{P1}";
	records[0].retentionMinutes = 12.5;
	records[0].sipAbundancePct = 37.0;
	records[0].charge = 2;
	records[0].precursors = {
		{500.2, 0.4}, {500.7, 1.0}, {500.9, 0.8}, {501.1, 0.1}};
	records[0].fragments = {
		{300.1234, 1.0F, 120.0F, 2, static_cast<uint8_t>('b'), {0, 0, 0}},
		{301.1234, 0.8F, 96.0F, 2, static_cast<uint8_t>('b'), {0, 0, 0}},
		{302.1234, 0.6F, 72.0F, 2, static_cast<uint8_t>('b'), {0, 0, 0}},
		{303.1234, 0.1F, 12.0F, 2, static_cast<uint8_t>('b'), {0, 0, 0}},
		{700.4321, 0.5F, 60.0F, 5, static_cast<uint8_t>('y'), {0, 0, 0}}};
	records[1].psmId = "scan.20.20.3";
	records[1].peptide = "[ANOTHER]";
	records[1].proteins = "{P2}";
	records[1].retentionMinutes = 22.0;
	records[1].sipAbundancePct = 38.0;
	records[1].charge = 3;
	records[1].precursors = {{600.3, 1.0}};
	records[1].fragments = {
		{400.2222, 1.0F, 80.0F, 3, static_cast<uint8_t>('b'), {0, 0, 0}}};
	for (size_t i = 2; i < 300; ++i)
	{
		records.emplace_back();
		auto &record = records.back();
		record.psmId = "scan." + std::to_string(i) + "." + std::to_string(i) + ".2";
		record.peptide = "[TESTPEPTIDE]";
		record.proteins = "{P_TEST}";
		record.retentionMinutes = 30.0;
		record.sipAbundancePct = 40.0;
		record.charge = 2;
		record.precursors = {{800.0 + static_cast<double>(i), 1.0}};
		record.fragments = {{1000.0 + static_cast<double>(i), 1.0F, 1.0F,
			1, static_cast<uint8_t>('b'), {0, 0, 0}}};
	}

	std::string error;
	sipros::SpectraIndexBuildStats buildStats;
	bool ok = expect(sipros::SpectraIndex::write(
		path.string(), metadata, records, error, 4, &buildStats), error);
	ok = expect(records[0].psmId.empty() && records[0].precursors.empty() &&
		records[0].fragments.empty() && records[1].psmId.empty() &&
		records[1].precursors.empty() && records[1].fragments.empty(),
		"SFI writer retained consumed peak buffers") && ok;
	sipros::SpectraIndex index;
	ok = expect(index.load(path.string(), error), error) && ok;
	ok = expect(index.recordCount() == 300, "wrong SFI record count") && ok;
	ok = expect(buildStats.threadsUsed == 2 && buildStats.blockCount == 2,
		"SFI packed product index was not built in parallel") && ok;
	ok = expect(buildStats.fileBytes == std::filesystem::file_size(path),
		"SFI build statistics contain the wrong file size") && ok;
	ok = expect(buildStats.rtBinCount == index.rtBinCount() && index.rtBinCount() >= 3,
		"SFI RT-bin build statistics or mapped count are wrong") && ok;
	ok = expect(index.metadata().targetSipAbundancePct ==
		sipros::SpectraIndexMetadata::MixedSipAbundancePct,
		"mixed-abundance SFI metadata was not preserved") && ok;
	ok = expect(index.metadata().envelopeTopN == 3,
		"compact-envelope top-N metadata was not preserved") && ok;
	const auto firstRange = index.precursorMzRange(500.69, 500.71);
	ok = expect(firstRange.second - firstRange.first == 1,
		"precursor apex lookup did not find exactly one record") && ok;
	const uint32_t firstId = firstRange.first;
	ok = expect(index.psmId(firstId) == "scan.10.10.2", "wrong mapped PSM id") && ok;
	ok = expect(index.peptide(firstId) == "[PEPTIDE]", "wrong mapped peptide") && ok;
	ok = expect(std::fabs(index.record(firstId).sipAbundancePct - 37.0) < 1e-12,
		"wrong per-record SIP abundance") && ok;
	ok = expect(std::fabs(index.record(firstId).retentionMinutes - 12.5) < 1e-12,
		"wrong indexed retention time") && ok;
	const auto fragments = index.fragments(firstId);
	ok = expect(fragments.second - fragments.first == 4 &&
		std::fabs(fragments.first->mz() - 300.123) < 1e-12,
		"compact fragment envelope was not top-3 pruned/quantized") && ok;
	ok = expect(fragments.first->ionPosition() == 2 &&
		fragments.first->ionKind() == static_cast<uint8_t>('b') &&
		std::fabs(fragments.first->theoreticalIntensity() - 1.0F) < 1e-7,
		"packed fragment metadata or theoretical intensity changed") && ok;
	auto experimental = index.experimentalIntensities(firstId);
	ok = expect(std::fabs(experimental.next() - 120.0F) < 1e-6 &&
		std::fabs(experimental.next() - 96.0F) < 1e-6 &&
		std::fabs(experimental.next() - 72.0F) < 1e-6 &&
		std::fabs(experimental.next() - 60.0F) < 1e-6,
		"sparse experimental intensities did not round-trip") && ok;
	const auto precursorPeaks = index.precursors(firstId);
	ok = expect(precursorPeaks.second - precursorPeaks.first == 3,
		"compact precursor envelope was not top-3 pruned") && ok;
	const uint32_t block = index.precursorBlockForRecord(firstId);
	const auto rtBins = index.rtBins(block, 12.0, 13.0);
	ok = expect(rtBins.first != nullptr && rtBins.second - rtBins.first == 1,
		"RT-aware lookup did not select exactly one product segment") && ok;
	const auto postings = rtBins.first == nullptr
		? std::make_pair(static_cast<const sipros::SpectraIndexFragmentPosting *>(nullptr),
			static_cast<const sipros::SpectraIndexFragmentPosting *>(nullptr))
		: index.fragmentRange(block, *rtBins.first, 300.12, 300.13);
	ok = expect(postings.first != nullptr && postings.second > postings.first,
		"product-ion sparse posting lookup returned no hit") && ok;
	if (postings.first != nullptr && postings.second > postings.first)
	{
		ok = expect(index.postingRecordId(block, *postings.first) == firstId,
			"product-ion posting points to the wrong record") && ok;
	}
	const auto wrongRtBins = index.rtBins(block, 20.0, 23.0);
	if (wrongRtBins.first != nullptr)
	{
		const auto wrongRtPostings = index.fragmentRange(
			block, *wrongRtBins.first, 300.12, 300.13);
		ok = expect(wrongRtPostings.first == nullptr,
			"RT-aware product lookup leaked a fragment from another RT segment") && ok;
	}
	index.close();
	std::filesystem::remove(path, ec);
	return ok ? 0 : 1;
}
