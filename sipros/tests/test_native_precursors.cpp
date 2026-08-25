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

void writeFixture(const std::filesystem::path &path, int reactionCharge,
				  bool hasNearbyPrecursor = true)
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
		hasNearbyPrecursor
			? std::vector<double>{499.5, 500.0, 500.01, 100.0, 200.0}
			: std::vector<double>{400.0, 410.0, 420.0, 100.0, 200.0});
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

std::vector<MS2Scan *> readNative(
	const std::filesystem::path &path,
	sipros::RaxportMs1Data *ms1Data = nullptr)
{
	std::vector<MS2Scan *> scans;
	std::string error;
	sipros::RaxportReadOptions options;
	options.precursorSource = sipros::PrecursorSource::IsolationWindow;
	check(options.precursorMatchScanRadius == 5,
		"Regular post-score MS1 radius is not five scans");
	check(sipros::readRaxportHdf5Scans(
		path.string(), scans, ms1Data, error, nullptr, options), error);
	return scans;
}

std::vector<MS2Scan *> readCandidates(const std::filesystem::path &path)
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

void checkPostScoreFiveScanMatch()
{
	sipros::RaxportMs1Data data;
	data.scans.resize(12);
	for (size_t i = 0; i < data.scans.size(); ++i)
	{
		data.scans[i].scanNumber = 100 + static_cast<int>(i);
		data.scans[i].retentionTime = 10.0 + 0.1 * static_cast<double>(i);
		data.scanNumberToIndex.emplace(data.scans[i].scanNumber, i);
	}
	const int charge = 2;
	const double calculatedMass = 1000.0;
	const double targetMz = calculatedMass / charge +
		ProNovoConfig::getProtonMass();
	for (size_t i : {size_t{0}, size_t{10}, size_t{11}})
	{
		data.scans[i].mz = {targetMz};
		data.scans[i].intensity = {100.0 + static_cast<double>(i)};
		data.scans[i].charge = {charge};
	}

	const auto boundary = sipros::findRaxportPrecursorMatch(
		data, 105, 11.09, charge, calculatedMass, 1.0033548,
		{0, 1, 2}, 0.01, 5);
	check(boundary.found && boundary.ms1ScanNumber == 110,
		"post-score precursor search did not include the +5 MS1 boundary");
	check(std::abs(boundary.rtDiffSeconds - 5.4) < 1e-9,
		"post-score precursor RT difference is not in seconds");
	const auto tooNarrow = sipros::findRaxportPrecursorMatch(
		data, 105, 11.09, charge, calculatedMass, 1.0033548,
		{0, 1, 2}, 0.01, 4);
	check(!tooNarrow.found && tooNarrow.rtDiffSeconds == -1.0,
		"post-score precursor search escaped its configured scan radius");
	data.scans[10].charge = {0};
	const auto unknownCharge = sipros::findRaxportPrecursorMatch(
		data, 105, 11.09, charge, calculatedMass, 1.0033548,
		{0}, 0.01, 5);
	check(unknownCharge.found,
		"unassigned MS1 peak did not support a post-score charge-2 PSM");
}


int main()
{
	const unsigned long processId =
#if defined(__unix__) || defined(__APPLE__)
		static_cast<unsigned long>(getpid());
#else
		0;
#endif
	const std::filesystem::path firstPath = std::filesystem::temp_directory_path() /
		("sipros_native_precursor_low_" + std::to_string(processId) + ".h5");
	const std::filesystem::path secondPath = std::filesystem::temp_directory_path() /
		("sipros_native_precursor_high_" + std::to_string(processId) + ".h5");
	const std::filesystem::path noPrecursorPath = std::filesystem::temp_directory_path() /
		("sipros_native_precursor_absent_" + std::to_string(processId) + ".h5");
	std::vector<MS2Scan *> first;
	std::vector<MS2Scan *> second;
	std::vector<MS2Scan *> noPrecursor;
	sipros::RaxportMs1Data firstMs1;
	try
	{
		check(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"could not load Regular profile");
		checkPostScoreFiveScanMatch();
		writeFixture(firstPath, 1);
		writeFixture(secondPath, 7);
		writeFixture(noPrecursorPath, 2, false);
		first = readNative(firstPath, &firstMs1);
		second = readNative(secondPath);
		noPrecursor = readNative(noPrecursorPath);
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
		check(lowReaction.vIsolationWindowsMz ==
			std::vector<std::pair<double, double>>({{500.0, 2.0}}) &&
			highReaction.vIsolationWindowsMz ==
			std::vector<std::pair<double, double>>({{500.0, 2.0}}),
			"native mode did not preserve the acquisition isolation window");
		check(lowReaction.iParentChargeStates.empty() &&
			highReaction.iParentChargeStates.empty() &&
			lowReaction.dParentMZs.empty() && highReaction.dParentMZs.empty(),
			"window-only mode created candidates from MS1 peaks");
		check(lowReaction.iMaxCandidateCharge == 0 &&
			highReaction.iMaxCandidateCharge == 0,
			"reader assigned a precursor charge before window search");
		check(std::abs(lowReaction.dParentMass - highReaction.dParentMass) < 1e-12 &&
			std::abs(lowReaction.dParentNeutralMass -
				highReaction.dParentNeutralMass) < 1e-12,
			"reaction charge contaminated window-only precursor limits");
		check(lowReaction.dParentMass == 0.0 &&
			lowReaction.dParentNeutralMass == 0.0,
			"reader assigned a searchable mass before window expansion");
		check(noPrecursor.size() == 1 &&
			noPrecursor.front()->iParentChargeStates.empty() &&
			noPrecursor.front()->vIsolationWindowsMz ==
				std::vector<std::pair<double, double>>({{500.0, 2.0}}),
			"isolation-window scan without an MS1 hypothesis was discarded");
		check(firstMs1.scans.size() == 1 &&
			firstMs1.scans.front().scanNumber == 100 &&
			firstMs1.scans.front().mz.size() == 3,
			"window-only scan read did not retain MS1 data for post-score matching");

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
		destroy(noPrecursor);
		std::error_code ignored;
		std::filesystem::remove(firstPath, ignored);
		std::filesystem::remove(secondPath, ignored);
		std::filesystem::remove(noPrecursorPath, ignored);
		std::cout << "ok: Regular candidates are isolation-window-only and retain MS1 data\n";
		return 0;
	}
	catch (const std::exception &ex)
	{
		destroy(first);
		destroy(second);
		destroy(noPrecursor);
		std::error_code ignored;
		std::filesystem::remove(firstPath, ignored);
		std::filesystem::remove(secondPath, ignored);
		std::filesystem::remove(noPrecursorPath, ignored);
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
