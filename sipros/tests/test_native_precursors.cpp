#include "RaxportHdf5Reader.h"
#include "ms2scan.h"
#include "proNovoConfig.h"

#include <H5Cpp.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace
{

void check(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

template <typename T>
void write1D(H5::H5File &file, const std::string &path,
			 const H5::PredType &type, const std::vector<T> &values)
{
	const hsize_t size = values.size();
	H5::DataSpace space(1, &size);
	H5::DataSet data = file.createDataSet(path, type, space);
	if (!values.empty())
	{
		data.write(values.data(), type);
	}
}

void writeFixture(const fs::path &path, int reactionCharge)
{
	H5::H5File file(path.string(), H5F_ACC_TRUNC);
	const int schemaVersion = 6;
	H5::DataSpace scalar(H5S_SCALAR);
	H5::Attribute attribute = file.createAttribute(
		"schema_version", H5::PredType::NATIVE_INT, scalar);
	attribute.write(H5::PredType::NATIVE_INT, &schemaVersion);
	file.createGroup("/scans");
	file.createGroup("/peaks");
	file.createGroup("/reactions");
	file.createGroup("/precursor_candidates");

	write1D(file, "/scans/scan_number", H5::PredType::NATIVE_INT,
		std::vector<int>{100, 101});
	write1D(file, "/scans/ms_order", H5::PredType::NATIVE_INT,
		std::vector<int>{1, 2});
	write1D(file, "/scans/retention_time", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{1.0, 1.1});
	write1D(file, "/scans/tic", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{60.0, 30.0});
	write1D(file, "/scans/parent_scan_number", H5::PredType::NATIVE_INT,
		std::vector<int>{0, 100});
	write1D(file, "/scans/reaction_start", H5::PredType::NATIVE_LLONG,
		std::vector<long long>{-1, 0});
	write1D(file, "/scans/reaction_count", H5::PredType::NATIVE_INT,
		std::vector<int>{0, 1});
	write1D(file, "/scans/peak_start", H5::PredType::NATIVE_LLONG,
		std::vector<long long>{0, 3});
	write1D(file, "/scans/peak_count", H5::PredType::NATIVE_INT,
		std::vector<int>{3, 2});

	write1D(file, "/peaks/mz", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{499.5, 500.0, 500.01, 100.0, 200.0});
	write1D(file, "/peaks/intensity", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{10.0, 20.0, 30.0, 10.0, 20.0});
	write1D(file, "/peaks/charge", H5::PredType::NATIVE_INT,
		std::vector<int>{2, 0, 25, 0, 0});

	write1D(file, "/reactions/precursor_mass", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{500.0});
	write1D(file, "/reactions/charge_state", H5::PredType::NATIVE_INT,
		std::vector<int>{reactionCharge});
	write1D(file, "/reactions/isolation_width", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{2.0});
	write1D(file, "/reactions/isolation_width_offset", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{0.0});
	write1D(file, "/reactions/candidate_start", H5::PredType::NATIVE_LLONG,
		std::vector<long long>{0});
	write1D(file, "/reactions/candidate_count", H5::PredType::NATIVE_INT,
		std::vector<int>{0});
	write1D(file, "/precursor_candidates/charge", H5::PredType::NATIVE_INT,
		std::vector<int>{});
	write1D(file, "/precursor_candidates/mz", H5::PredType::NATIVE_DOUBLE,
		std::vector<double>{});
	write1D(file, "/precursor_candidates/charge_source", H5::PredType::NATIVE_INT,
		std::vector<int>{});
	write1D(file, "/precursor_candidates/isotope_match_count", H5::PredType::NATIVE_INT,
		std::vector<int>{});
}

std::vector<MS2Scan *> readNative(const fs::path &path)
{
	std::vector<MS2Scan *> scans;
	std::string error;
	sipros::RaxportReadOptions options;
	options.precursorSource = sipros::PrecursorSource::Ms1Neighborhood;
	options.ms1NeighborhoodRadius = 2;
	check(sipros::readRaxportHdf5Scans(
		path.string(), scans, nullptr, error, nullptr, options), error);
	return scans;
}

std::vector<MS2Scan *> readCandidates(const fs::path &path)
{
	std::vector<MS2Scan *> scans;
	std::string error;
	sipros::RaxportReadOptions options;
	options.precursorSource = sipros::PrecursorSource::RaxportCandidates;
	check(sipros::readRaxportHdf5Scans(
		path.string(), scans, nullptr, error, nullptr, options), error);
	return scans;
}

void destroy(std::vector<MS2Scan *> &scans)
{
	for (MS2Scan *scan : scans)
	{
		delete scan;
	}
	scans.clear();
}

} // namespace

int main()
{
	const unsigned long processId =
#if defined(__unix__) || defined(__APPLE__)
		static_cast<unsigned long>(getpid());
#else
		0;
#endif
	const fs::path firstPath = fs::temp_directory_path() /
		("sipros_native_precursor_low_" + std::to_string(processId) + ".h5");
	const fs::path secondPath = fs::temp_directory_path() /
		("sipros_native_precursor_high_" + std::to_string(processId) + ".h5");
	std::vector<MS2Scan *> first;
	std::vector<MS2Scan *> second;
	try
	{
		check(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"could not load Regular profile");
		writeFixture(firstPath, 1);
		writeFixture(secondPath, 7);
		first = readNative(firstPath);
		second = readNative(secondPath);
		check(first.size() == 1 && second.size() == 1,
			"native fixture did not produce one MS2 scan");
		const MS2Scan &lowReaction = *first.front();
		const MS2Scan &highReaction = *second.front();
		check(!lowReaction.bUseReactionChargeForScoring &&
			!highReaction.bUseReactionChargeForScoring,
			"native mode left reaction-charge scoring enabled");
		check(lowReaction.iParentChargeState == 1 &&
			highReaction.iParentChargeState == 7,
			"fixture did not preserve reaction charge as metadata");
		check(lowReaction.iParentChargeStates ==
			highReaction.iParentChargeStates &&
			lowReaction.dParentMZs == highReaction.dParentMZs,
			"reaction charge changed native MS1 hypotheses");
		check(lowReaction.iParentChargeStates ==
			std::vector<int>({2, 2, 3, 4, 25}),
			"native charges did not preserve positive peak charge and zero -> 2/3/4");
		check(lowReaction.iMaxCandidateCharge == 25 &&
			highReaction.iMaxCandidateCharge == 25,
			"reaction charge contaminated the candidate charge limit");
		check(lowReaction.iParentChargeStates.back() == 25,
			"a positive MS1 peak charge above seven was discarded");
		check(std::abs(lowReaction.dParentMass - highReaction.dParentMass) < 1e-12 &&
			std::abs(lowReaction.dParentNeutralMass -
				highReaction.dParentNeutralMass) < 1e-12,
			"reaction charge contaminated native precursor mass limits");
		check(std::abs(lowReaction.dParentMass - 500.01 * 25.0) < 1e-9,
			"native precursor mass limit did not use the charge-25 MS1 peak");

		sipros::RaxportMs1Data featureMs1;
		std::string featureError;
		check(sipros::readRaxportHdf5Ms1Data(
			firstPath.string(), featureMs1, featureError), featureError);
		check(featureMs1.scans.size() == 1 &&
			featureMs1.scans.front().scanNumber == 100,
			"MS1-only feature reader still depended on legacy precursor candidates");
		auto candidateOnly = readCandidates(firstPath);
		check(candidateOnly.empty(),
			"Raxport candidate mode fell back to the reaction precursor");
		destroy(candidateOnly);

		destroy(first);
		destroy(second);
		std::error_code ignored;
		fs::remove(firstPath, ignored);
		fs::remove(secondPath, ignored);
		std::cout << "ok: native MS1 hypotheses are reaction-charge invariant\n";
		return 0;
	}
	catch (const std::exception &ex)
	{
		destroy(first);
		destroy(second);
		std::error_code ignored;
		fs::remove(firstPath, ignored);
		fs::remove(secondPath, ignored);
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
