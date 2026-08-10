#include "PinWriter.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void check(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<std::string> splitTabs(const std::string &line)
{
    std::vector<std::string> fields;
    size_t begin = 0;
    while (begin <= line.size())
    {
        const size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string::npos
                ? std::string::npos
                : end - begin));
        if (end == std::string::npos)
            break;
        begin = end + 1;
    }
    return fields;
}

void checkPin(const std::filesystem::path &path,
              int expectedRawCount,
              double expectedFitScore,
              double expectedObservedMass,
              double expectedRtDiffSeconds,
              double expectedAbundanceDiff)
{
    std::ifstream input(path);
    check(static_cast<bool>(input),
          "cannot read test PIN: " + path.string());
    std::string headerLine;
    std::string rowLine;
    check(static_cast<bool>(std::getline(input, headerLine)),
          "test PIN has no header");
    check(static_cast<bool>(std::getline(input, rowLine)),
          "test PIN has no data row");

    const std::vector<std::string> header =
        splitTabs(headerLine);
    const std::vector<std::string> row =
        splitTabs(rowLine);
    check(header.size() == row.size(),
          "PIN header/data field counts differ");

    const auto rawIt = std::find(
        header.begin(), header.end(),
        "isotopicPeakNumbers");
    const auto scoreIt = std::find(
        header.begin(), header.end(),
        "MS1IsotopeFitScore");
    const auto abundanceIt = std::find(
        header.begin(), header.end(),
        "MS1IsotopicAbundances");
    const auto observedMassIt = std::find(
        header.begin(), header.end(), "ObservedMass");
    const auto rtDiffIt = std::find(
        header.begin(), header.end(), "absPrecursorRtDiffSeconds");
    const auto ms2AbundanceIt = std::find(
        header.begin(), header.end(), "MS2IsotopicAbundances");
    const auto abundanceDiffIt = std::find(
        header.begin(), header.end(), "isotopicAbundanceDiffs");
    check(rawIt != header.end() &&
              rtDiffIt == rawIt + 1 &&
              scoreIt == rtDiffIt + 1 &&
              abundanceIt == scoreIt + 1,
          "MS1 raw-count/RT-distance/fit-score columns are missing or misordered");
    check(observedMassIt != header.end(),
          "PIN is missing ObservedMass metadata");
    check(rtDiffIt != header.end(),
          "PIN is missing absPrecursorRtDiffSeconds");
    check(ms2AbundanceIt != header.end() &&
              abundanceDiffIt == ms2AbundanceIt + 1,
          "isotopicAbundanceDiffs is missing or misordered");
    check(std::find(header.begin(), header.end(), "hasMs1Precursor") ==
              header.end(),
          "obsolete hasMs1Precursor feature is still present");

    const size_t rawIndex = static_cast<size_t>(
        rawIt - header.begin());
    const size_t scoreIndex = static_cast<size_t>(
        scoreIt - header.begin());
    check(std::stoi(row[rawIndex]) == expectedRawCount,
          "PIN changed the raw isotopic peak count");
    check(std::fabs(
              std::stod(row[scoreIndex]) -
              expectedFitScore) < 1e-6,
          "PIN changed the MS1 isotope fit score");
    check(std::fabs(std::stod(row[static_cast<size_t>(
                            observedMassIt - header.begin())]) -
                    expectedObservedMass) < 1e-6,
          "PIN changed the observed precursor mass");
    check(std::fabs(std::stod(row[static_cast<size_t>(
                            rtDiffIt - header.begin())]) -
                    expectedRtDiffSeconds) < 1e-6,
          "PIN changed the precursor RT-distance feature");
    check(std::fabs(std::stod(row[static_cast<size_t>(
                            abundanceDiffIt - header.begin())]) -
                    expectedAbundanceDiff) < 1e-6,
          "PIN changed the signed MS1-minus-MS2 abundance feature");
    const double printedMs1 = std::stod(row[static_cast<size_t>(
        abundanceIt - header.begin())]);
    const double printedMs2 = std::stod(row[static_cast<size_t>(
        ms2AbundanceIt - header.begin())]);
    const double printedDifference = std::stod(row[static_cast<size_t>(
        abundanceDiffIt - header.begin())]);
    check(std::fabs(printedDifference - (printedMs1 - printedMs2)) < 2e-6,
          "printed abundance difference is inconsistent with printed MS1-MS2");
	if (expectedRawCount == 1)
	{
		const auto intensityIt = std::find(
			header.begin(), header.end(), "log10_precursorIntensities");
		check(intensityIt != header.end() &&
			intensityIt + 2 < header.end() &&
			*(intensityIt + 1) == "ddaResidualRank" &&
			*(intensityIt + 2) == "ddaResidualScore",
			"DDA+ PIN columns are missing or misordered");
		const size_t rtDiffIndex = static_cast<size_t>(
			rtDiffIt - header.begin());
		const size_t intensityIndex = static_cast<size_t>(
			intensityIt - header.begin());
		check(std::fabs(std::stod(row[rtDiffIndex]) - 3.25) < 1e-6 &&
			row[intensityIndex + 1] == "2" &&
			std::fabs(std::stod(row[intensityIndex + 2]) - 12.5) < 1e-6,
			"DDA+ PIN feature values changed");
	}

    const auto matchedBIt = std::find(
        header.begin(), header.end(), "matchedBIons");
    if (matchedBIt != header.end())
    {
        const auto matchedYIt = matchedBIt + 1;
        const auto consecutiveBIt = matchedBIt + 2;
        const auto consecutiveYIt = matchedBIt + 3;
        check(consecutiveYIt < header.end() &&
                  *matchedYIt == "matchedYIons" &&
                  *consecutiveBIt == "maxConsecutiveBIons" &&
                  *consecutiveYIt == "maxConsecutiveYIons",
              "classic PIN B/Y ion columns are missing or misordered");
        const size_t featureIndex = static_cast<size_t>(
            matchedBIt - header.begin());
        check(row[featureIndex] == "3" &&
                  row[featureIndex + 1] == "4" &&
                  row[featureIndex + 2] == "2" &&
                  row[featureIndex + 3] == "3",
              "classic PIN B/Y ion feature values changed");
    }
}

sipPSM makeClassicPsm()
{
    sipPSM psm;
    psm.fileNames = {"sample"};
    psm.scanNumbers = {42};
    psm.calculatedParentMasses = {1000.0};
    psm.measuredParentMasses = {999.5};
    psm.retentionTimes = {12.5f};
    psm.ranks = {1};
    psm.parentCharges = {2};
    psm.massErrors = {0.5};
    psm.isotopicMassWindowShifts = {0};
    psm.mzShiftFromisolationWindowCenters = {0.0};
    psm.peptideLengths = {7};
    psm.missCleavageSiteNumbers = {0};
    psm.PTMnumbers = {0};
    psm.isotopicPeakNumbers = {1};
    psm.MS1IsotopeFitScores = {0.0};
    psm.MS1IsotopicAbundances = {1.2};
    psm.MS2IsotopicAbundances = {1.1};
    psm.isotopicAbundanceDiffs = {0.1};
    psm.WDPscores = {10.0f};
    psm.XcorrScores = {2.0f};
    psm.MVHscores = {20.0f};
    psm.diffScores = {0.0f};
    psm.matchedBIons = {3};
    psm.matchedYIons = {4};
    psm.maxConsecutiveBIons = {2};
    psm.maxConsecutiveYIons = {3};
    psm.precursorIntensities = {1000.0};
	psm.precursorRtDiffSeconds = {3.25};
	psm.ddaResidualRanks = {2};
	psm.ddaResidualScores = {12.5f};
    psm.identifiedPeptides = {"K[PEPTIDE]R"};
    psm.originalPeptides = {"K[PEPTIDE]R"};
    psm.proteinNames = {"{protein}"};
    psm.isDecoys = {false};
    return psm;
}
} // namespace

int main()
{
    try
    {
        const auto nonce =
            std::chrono::steady_clock::now()
                .time_since_epoch().count();
        const std::filesystem::path directory =
            std::filesystem::temp_directory_path();
        const std::filesystem::path classicPath =
            directory /
            ("sipros_classic_pin_schema_" +
             std::to_string(nonce) + ".pin");
        const std::filesystem::path spectraPath =
            directory /
            ("sipros_spectra_pin_schema_" +
             std::to_string(nonce) + ".pin");

        PinWriter::writePecorlatorPin(
            classicPath.string(), {makeClassicPsm()}, false);
        checkPin(classicPath, 1, 0.0, 999.5, 3.25, 0.1);

        PinWriter::SearchSpectraPinRow spectraRow;
        spectraRow.label = 1;
        spectraRow.scanNumber = 43;
        spectraRow.rank = 1;
        spectraRow.isotopicPeakNumber = 4;
        spectraRow.ms1IsotopeFitScore = 0.75;
        spectraRow.ms1IsotopicAbundance = 2.3;
        spectraRow.ms2IsotopicAbundance = 1.1;
        spectraRow.observedMass = 998.5;
        spectraRow.precursorRtDiffSeconds = 7.5;
        spectraRow.peptide = "K[PEPTIDE]R";
        spectraRow.proteins = "{protein}";
        check(PinWriter::writeSearchSpectraPin(
                  spectraPath.string(), "sample",
                  {spectraRow}) == 1,
              "search-spectra PIN writer returned no rows");
        checkPin(spectraPath, 4, 0.75, 998.5, 7.5, 1.2);

        std::filesystem::remove(classicPath);
        std::filesystem::remove(spectraPath);
        std::cout
            << "ok: both PIN schemas preserve raw MS1 counts "
               "and 0-1 fit scores"
            << std::endl;
        return EXIT_SUCCESS;
    }
    catch (const std::exception &error)
    {
        std::cerr << "FAIL: " << error.what()
                  << std::endl;
        return EXIT_FAILURE;
    }
}
