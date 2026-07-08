#include "SiprosSearchRunner.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "averagine.h"
#include "CometSearchMod.h"
#include "MVH.h"
#include "ms2scanvector.h"
#include "proNovoConfig.h"
#include "PSMfeatureExtractor.h"

namespace fs = std::filesystem;

namespace sipros
{

std::string TextUtils::toLower(std::string value)
{
	for (char &ch : value)
	{
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	return value;
}

std::string TextUtils::trim(const std::string &value)
{
	size_t start = 0;
	while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start])) != 0)
	{
		++start;
	}
	size_t end = value.size();
	while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1])) != 0)
	{
		--end;
	}
	return value.substr(start, end - start);
}

std::vector<std::string> TextUtils::splitTab(const std::string &line)
{
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, static_cast<char>(9)))
	{
		fields.push_back(field);
	}
	if (!line.empty() && line.back() == static_cast<char>(9))
	{
		fields.emplace_back();
	}
	return fields;
}

std::vector<std::string> SiprosSearchRunner::listFilesWithExtensions(const std::string &directory,
												 const std::vector<std::string> &extensions)
{
	std::vector<std::string> files;
	if (!fs::is_directory(directory))
	{
		return files;
	}
	for (const auto &entry : fs::directory_iterator(directory))
	{
		if (!entry.is_regular_file())
		{
			continue;
		}
		const std::string ext = entry.path().extension().string();
		for (const std::string &expected : extensions)
		{
			if (ext == expected)
			{
				files.push_back(entry.path().string());
				break;
			}
		}
	}
	std::sort(files.begin(), files.end());
	return files;
}

static int atomIndex(char atom)
{
	const std::string atoms = "CHONPS";
	const size_t pos = atoms.find(static_cast<char>(std::toupper(static_cast<unsigned char>(atom))));
	return pos == std::string::npos ? -1 : static_cast<int>(pos);
}

bool TextUtils::parseSipAtomSpec(const std::string &spec, char &sipAtom, int &sipIsotopeMassNumber)
{
	const std::string value = TextUtils::trim(spec);
	if (value.empty())
	{
		return false;
	}
	const char atom = static_cast<char>(std::toupper(static_cast<unsigned char>(value[0])));
	if (atomIndex(atom) < 0)
	{
		return false;
	}
	sipAtom = atom;
	if (value.size() == 1)
	{
		sipIsotopeMassNumber = atom == 'O' ? 18 : atom == 'S' ? 34 : -1;
		return true;
	}
	for (size_t i = 1; i < value.size(); ++i)
	{
		if (!std::isdigit(static_cast<unsigned char>(value[i])))
		{
			return false;
		}
	}
	try
	{
		sipIsotopeMassNumber = std::stoi(value.substr(1));
	}
	catch (const std::exception &)
	{
		return false;
	}
	return (atom == 'C' && sipIsotopeMassNumber == 13) ||
		   (atom == 'H' && sipIsotopeMassNumber == 2) ||
		   (atom == 'N' && sipIsotopeMassNumber == 15) ||
		   (atom == 'O' && sipIsotopeMassNumber == 18) ||
		   (atom == 'S' && (sipIsotopeMassNumber == 33 || sipIsotopeMassNumber == 34));
}

static bool parseDouble(const std::string &text, double &value)
{
	char *end = nullptr;
	value = std::strtod(text.c_str(), &end);
	return end != text.c_str() && end != nullptr && *end == 0;
}

static bool parseInt(const std::string &text, int &value)
{
	char *end = nullptr;
	long parsed = std::strtol(text.c_str(), &end, 10);
	if (end == text.c_str() || end == nullptr || *end != 0)
	{
		return false;
	}
	value = static_cast<int>(parsed);
	return true;
}

static bool isHdf5Path(const std::string &path)
{
	const std::string lower = TextUtils::toLower(fs::path(path).extension().string());
	return lower == ".h5" || lower == ".hdf5";
}

static std::string normalizeSipSpec(const std::string &spec, char sipAtom, int sipIsotopeMassNumber)
{
	std::ostringstream out;
	out << static_cast<char>(std::toupper(static_cast<unsigned char>(sipAtom)));
	if (sipIsotopeMassNumber > 0)
	{
		out << sipIsotopeMassNumber;
	}
	else
	{
		std::string trimmed = TextUtils::trim(spec);
		if (trimmed.size() > 1)
		{
			out << trimmed.substr(1);
		}
	}
	return out.str();
}

static std::string formatSipSearchName(const std::string &sipSpec, double pct)
{
	std::ostringstream out;
	out << "SIP_" << sipSpec << "_"
		<< std::setw(7) << std::setfill('0') << std::fixed << std::setprecision(3) << pct
		<< "Pct";
	return out.str();
}

static std::vector<double> parsePctRange(const std::string &rangeSpec, double stepPct)
{
	if (rangeSpec.empty())
	{
		throw std::runtime_error("SIP range is required; use -b <pct|lower-upper>");
	}
	if (stepPct <= 0.0)
	{
		throw std::runtime_error("SIP step must be positive");
	}
	const size_t dash = rangeSpec.find('-');
	double lower = 0.0;
	double upper = 0.0;
	if (dash == std::string::npos)
	{
		if (!parseDouble(rangeSpec, lower))
		{
			throw std::runtime_error("Invalid SIP pct value: " + rangeSpec);
		}
		return {lower};
	}
	if (!parseDouble(rangeSpec.substr(0, dash), lower) || !parseDouble(rangeSpec.substr(dash + 1), upper))
	{
		throw std::runtime_error("Invalid SIP range: " + rangeSpec);
	}
	if (upper < lower)
	{
		throw std::runtime_error("SIP range upper bound is lower than the lower bound: " + rangeSpec);
	}
	std::vector<double> values;
	for (int step = 0;; ++step)
	{
		const double value = lower + static_cast<double>(step) * stepPct;
		if (value > upper + 1e-9)
		{
			break;
		}
		values.push_back(value);
	}
	if (values.empty() || std::abs(values.back() - upper) > 1e-6)
	{
		values.push_back(upper);
	}
	return values;
}

static bool applySipAbundance(char sipAtom, double fraction)
{
	const int index = atomIndex(sipAtom);
	if (index < 0 || index >= static_cast<int>(ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution.size()))
	{
		return false;
	}
	ProNovoConfig::getSetSIPelement() = std::string(1, sipAtom);
	const bool updated = averagine::changeAtomProbability(
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[index].vProb, sipAtom, fraction);
	if (updated)
	{
		ProNovoConfig::setDeductionCoefficient();
		ProNovoConfig::getSetSIPelement() = std::string(1, sipAtom);
	}
	return updated;
}

static std::string fileNameForSpecId(const std::string &sampleName, const std::string &searchName, bool isDecoy)
{
	std::string base = sampleName;
	if (searchName.find("Pct") == std::string::npos)
	{
		return isDecoy ? base + ".decoy" : base;
	}
	std::string sanitized = searchName;
	std::replace(sanitized.begin(), sanitized.end(), '.', '_');
	base += "." + sanitized;
	return isDecoy ? base + ".decoy" : base;
}

static bool scoredPsmBetter(const ScoredPsmRow &left, const ScoredPsmRow &right)
{
	if (left.wdpScore != right.wdpScore)
	{
		return left.wdpScore > right.wdpScore;
	}
	if (left.xcorrScore != right.xcorrScore)
	{
		return left.xcorrScore > right.xcorrScore;
	}
	if (left.mvhScore != right.mvhScore)
	{
		return left.mvhScore > right.mvhScore;
	}
	if (left.searchName != right.searchName)
	{
		return left.searchName < right.searchName;
	}
	return left.identifiedPeptide < right.identifiedPeptide;
}

static void pruneScoredRowsByScan(std::vector<ScoredPsmRow> &rows, int topPsmsPerScan)
{
	const int topN = topPsmsPerScan > 0 ? topPsmsPerScan : ProNovoConfig::INTTOPKEEP;
	if (topN <= 0 || rows.empty())
	{
		return;
	}

	std::map<int, std::vector<ScoredPsmRow>> groupedRows;
	for (ScoredPsmRow &row : rows)
	{
		groupedRows[row.scanNumber].push_back(std::move(row));
	}

	std::vector<ScoredPsmRow> prunedRows;
	prunedRows.reserve(std::min(rows.size(), groupedRows.size() * static_cast<size_t>(topN)));
	for (auto &entry : groupedRows)
	{
		std::vector<ScoredPsmRow> &scanRows = entry.second;
		std::sort(scanRows.begin(), scanRows.end(), scoredPsmBetter);
		std::set<std::string> selectedPeptides;
		int rank = 1;
		for (ScoredPsmRow &row : scanRows)
		{
			if (!selectedPeptides.insert(row.identifiedPeptide).second)
			{
				continue;
			}
			row.rank = rank;
			prunedRows.push_back(std::move(row));
			++rank;
			if (rank > topN)
			{
				break;
			}
		}
	}
	rows.swap(prunedRows);
}

static sipPSM buildSipPsm(const std::string &sampleName,
						  const std::vector<ScoredPsmRow> &rows,
						  int topPsmsPerScan)
{
	sipPSM psm;
	psm.fileName = sampleName;
	std::vector<ScoredPsmRow> prunedRows = rows;
	pruneScoredRowsByScan(prunedRows, topPsmsPerScan);
	for (const ScoredPsmRow &row : prunedRows)
	{
		psm.fileNames.push_back(fileNameForSpecId(sampleName, row.searchName, row.isDecoy));
		psm.scanNumbers.push_back(row.scanNumber);
		psm.precursorScanNumbers.push_back(row.precursorScanNumber);
		psm.parentCharges.push_back(row.parentCharge);
		psm.isolationWindowCenterMZs.push_back(row.isolationWindowCenterMZ);
		psm.measuredParentMasses.push_back(row.measuredParentMass);
		psm.calculatedParentMasses.push_back(row.calculatedParentMass);
		psm.searchNames.push_back(row.searchName);
		psm.ranks.push_back(row.rank);
		psm.scores.push_back(row.wdpScore);
		psm.identifiedPeptides.push_back(row.identifiedPeptide);
		psm.originalPeptides.push_back(row.originalPeptide);
		psm.nakePeptides.push_back(row.nakedPeptide);
		psm.proteinNames.push_back(row.proteinNames);
		psm.retentionTimes.push_back(row.retentionTime);
		psm.MVHscores.push_back(row.mvhScore);
		psm.XcorrScores.push_back(row.xcorrScore);
		psm.WDPscores.push_back(row.wdpScore);
		psm.isDecoys.push_back(row.isDecoy);
	}
	return psm;
}

void SiprosSearchRunner::printUsage(std::ostream &out)
{
	out << "Usage:\n";
	out << "  sipros -f sample.h5 -c config.cfg -fasta proteins.fasta -o out "
		<< "[-a C13 -b 1-5 -s 1] [--pin-label 1|-1]\n";
	out << "  sipros -w hdf5_directory -c config.cfg -fasta proteins.fasta -o out "
		<< "[-a C13 -b 1-5 -s 1] [--pin-label 1|-1]\n\n";
	out << "Parameters:\n";
	out << "  -c <config.cfg>             one base configuration file\n";
	out << "  -f <sample.h5>              one Raxport HDF5 scan file\n";
	out << "  -w <directory>              directory of Raxport HDF5 scan files\n";
	out << "  -fasta <proteins.fasta>     one FASTA database; target/decoy orchestration is external\n";
	out << "  -o <directory>              output directory, default: out\n";
	out << "  -a <SIP atom/isotope>       SIP element such as C13, H2, N15, O18, S33, S34\n";
	out << "  -b <pct|lower-upper>        SIP percentage or inclusive range\n";
	out << "  -s, --step <pct>            SIP percentage step\n";
	out << "  --pin-label <1|-1>          label written to PIN rows; default: 1\n";
	out << "  --top-psms-per-scan <N>     final WDP-ranked PSMs retained per scan across SIP pct\n";
	out << "  --silent                    silence most standard output\n";
	out << "Unsupported legacy scan inputs (.FT1/.FT2/.mzML) and -g config directories are rejected.\n";
}

bool SiprosSearchRunner::initializeArguments(int argc, char **argv,
									 DatabaseSearchArguments &args,
									 std::ostream &out,
									 std::ostream &err) const
{
	std::vector<std::string> values;
	while (argc--)
	{
		values.push_back(*argv++);
	}
	bool valid = true;
	auto requireValue = [&](int &index, const std::string &option) -> std::string
	{
		if (index + 1 >= static_cast<int>(values.size()))
		{
			err << option << " requires a value\n";
			valid = false;
			return {};
		}
		return values[++index];
	};

	for (int i = 1; i <= static_cast<int>(values.size()) - 1; ++i)
	{
		const std::string &option = values[i];
		if (option == "-w")
		{
			args.workingDirectory = requireValue(i, option);
		}
		else if (option == "-c")
		{
			args.configFile = requireValue(i, option);
		}
		else if (option == "-f")
		{
			args.singleWorkingFile = requireValue(i, option);
		}
		else if (option == "-fasta")
		{
			args.fastaFile = requireValue(i, option);
		}
		else if (option == "-o")
		{
			args.outputDirectory = requireValue(i, option);
		}
		else if (option == "-a")
		{
			args.sipElementSpec = requireValue(i, option);
		}
		else if (option == "-b")
		{
			args.sipRangeSpec = requireValue(i, option);
		}
		else if (option == "-s" || option == "--step")
		{
			const std::string value = requireValue(i, option);
			if (valid && !parseDouble(value, args.sipStepPct))
			{
				err << "Invalid SIP step: " << value << "\n";
				return false;
			}
		}
		else if (option == "--pin-label")
		{
			const std::string value = requireValue(i, option);
			if (valid && (!parseInt(value, args.pinLabel) || (args.pinLabel != 1 && args.pinLabel != -1)))
			{
				err << "--pin-label must be 1 or -1\n";
				return false;
			}
		}
		else if (option == "--top-psms-per-scan")
		{
			const std::string value = requireValue(i, option);
			if (valid && (!parseInt(value, args.topPsmsPerScan) || args.topPsmsPerScan <= 0))
			{
				err << "--top-psms-per-scan must be a positive integer\n";
				return false;
			}
		}
		else if (option == "--silent" || option == "-q")
		{
			args.screenOutput = false;
		}
		else if (option == "-g")
		{
			err << "-g config-directory input was removed; use -c plus -a/-b/-s for SIP percentage search\n";
			return false;
		}
		else if ((option == "-h") || (option == "--help"))
		{
			args.showHelp = true;
			printUsage(out);
			return true;
		}
		else if (option == "-p")
		{
			const std::string value = requireValue(i, option);
			if (valid)
			{
				MVH::ProbabilityCutOff = atof(value.c_str());
				CometSearchMod::ProbabilityCutOff = MVH::ProbabilityCutOff;
			}
		}
		else
		{
			err << "Unknown option " << option << "\n\n";
			return false;
		}
	}
	if (!valid)
	{
		return false;
	}
	if (args.configFile.empty())
	{
		err << "One base configuration file is required with -c\n";
		return false;
	}
	if (args.workingDirectory.empty() && args.singleWorkingFile.empty())
	{
		args.workingDirectory = ".";
	}
	if (!args.workingDirectory.empty() && !args.singleWorkingFile.empty())
	{
		err << "Either one input scan file or a directory of input scan files must be specified\n";
		return false;
	}
	if (!args.singleWorkingFile.empty())
	{
		if (!isHdf5Path(args.singleWorkingFile))
		{
			err << "Raxport HDF5 scan input required (.h5 or .hdf5); unsupported file: "
				<< args.singleWorkingFile << "\n";
			return false;
		}
		args.scanFiles.push_back(args.singleWorkingFile);
	}
	else
	{
		args.scanFiles = listFilesWithExtensions(args.workingDirectory, {".h5", ".H5", ".hdf5", ".HDF5"});
		if (args.scanFiles.empty())
		{
			err << "no Raxport HDF5 scan file in the working directory\n";
			return false;
		}
	}
	if (args.outputDirectory.empty())
	{
		args.outputDirectory = "out";
	}
	return true;
}

int SiprosSearchRunner::runScan(const std::string &scanFile,
							  const DatabaseSearchArguments &args) const
{
	if (!ProNovoConfig::setFilename(args.configFile))
	{
		std::cerr << "Could not load config file " << args.configFile << std::endl;
		return 1;
	}
	if (!args.fastaFile.empty())
	{
		ProNovoConfig::setFASTAfilename(args.fastaFile);
	}
	const int topKeep = args.topPsmsPerScan > 0 ? args.topPsmsPerScan : ProNovoConfig::INTTOPKEEP;
	fs::create_directories(args.outputDirectory);

	MS2ScanVector scanVector(scanFile, args.outputDirectory, args.configFile, args.screenOutput);
	if (args.screenOutput)
	{
		std::cout << "Reading Raxport HDF5 scan file: " << scanFile << "\n";
		std::cout << "Using fasta file: " << ProNovoConfig::getFASTAfilename() << "\n";
		std::cout << "Using Configuration file: " << args.configFile << "\n";
	}
	if (!scanVector.loadMassData())
	{
		std::cerr << "Error: Failed to load file: " << scanFile << "\n";
		return 1;
	}

	std::vector<ScoredPsmRow> scoredRows;
	const bool configIsSip = ProNovoConfig::getSearchType() == "SIP";
	const bool hasSipControls = !args.sipElementSpec.empty();
	const bool isDecoyLabel = args.pinLabel < 0;
	if (hasSipControls || configIsSip)
	{
		if (!hasSipControls)
		{
			std::cerr << "SIP search requires -a <SIP atom/isotope>, -b <pct|lower-upper>, and -s <step>\n";
			return 1;
		}
		if (!configIsSip)
		{
			std::cerr << "SIP controls require a SIP Search_Type config file\n";
			return 1;
		}
		char sipAtom = 0;
		int sipIsotopeMassNumber = -1;
		if (!TextUtils::parseSipAtomSpec(args.sipElementSpec, sipAtom, sipIsotopeMassNumber))
		{
			std::cerr << "Unsupported SIP element/isotope: " << args.sipElementSpec << "\n";
			return 1;
		}
		std::vector<double> pctValues;
		try
		{
			pctValues = parsePctRange(args.sipRangeSpec, args.sipStepPct);
		}
		catch (const std::exception &ex)
		{
			std::cerr << ex.what() << "\n";
			return 1;
		}
		const std::string sipSpec = normalizeSipSpec(args.sipElementSpec, sipAtom, sipIsotopeMassNumber);
		for (double pct : pctValues)
		{
			scanVector.clearSearchResults();
			if (!applySipAbundance(sipAtom, pct / 100.0))
			{
				std::cerr << "Could not apply SIP abundance for " << args.sipElementSpec << " at " << pct << " percent\n";
				return 1;
			}
			ProNovoConfig::setSearchName(formatSipSearchName(sipSpec, pct));
			if (args.screenOutput)
			{
				std::cout << "Searching " << scanFile << " at " << ProNovoConfig::getSearchName() << "\n";
			}
			scanVector.startProcessingWdpSip();
			scanVector.appendScoredPsmRows(scoredRows, isDecoyLabel, topKeep);
			pruneScoredRowsByScan(scoredRows, topKeep);
		}
	}
	else
	{
		scanVector.clearSearchResults();
		if (ProNovoConfig::getSearchName().empty() || ProNovoConfig::getSearchName() == "Null" ||
			ProNovoConfig::getSearchName() == "NULL" || ProNovoConfig::getSearchName() == "null")
		{
			ProNovoConfig::setSearchName("SE");
		}
		scanVector.startProcessingMvh();
		scanVector.appendScoredPsmRows(scoredRows, isDecoyLabel, topKeep);
		pruneScoredRowsByScan(scoredRows, topKeep);
	}

	const std::string sampleName = fs::path(scanFile).stem().string();
	sipPSM psm = buildSipPsm(sampleName, scoredRows, topKeep);
	PSMfeatureExtractor extractor;
	extractor.sipPSMs.push_back(std::move(psm));
	if (!extractor.sipPSMs.front().scanNumbers.empty())
	{
		extractor.extractFeaturesForPsm(scanFile, extractor.sipPSMs.front());
	}
	const fs::path pinPath = fs::path(args.outputDirectory) / (sampleName + ".pin");
	extractor.writePecorlatorPin(pinPath.string(), false);
	if (args.screenOutput)
	{
		std::cout << "Wrote PIN: " << pinPath.string() << " (" << scoredRows.size() << " scored rows retained after top-N pruning)\n";
	}
	return 0;
}

} // namespace sipros
