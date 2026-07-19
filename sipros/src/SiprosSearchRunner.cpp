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

#include "CometSearchMod.h"
#include "MVH.h"
#include "fragmentindex.h"
#include "ms2scanvector.h"
#include "performancelog.h"
#include "proNovoConfig.h"
#include "PSMfeatureExtractor.h"
#include "PinWriter.h"

namespace fs = std::filesystem;

namespace sipros
{

static std::string formatElapsedSeconds(double seconds)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(1) << seconds << "s";
	return out.str();
}

static PerformanceTiming performanceTiming(double wallSeconds,
										double cpuSeconds)
{
	return {wallSeconds, cpuSeconds};
}

static std::string formatMebibytes(uint64_t bytes)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(1)
		<< static_cast<double>(bytes) / (1024.0 * 1024.0) << " MiB";
	return out.str();
}

static const char *precursorSourceName(PrecursorSource source)
{
	return source == PrecursorSource::Ms1Neighborhood
		? "ms1-neighborhood"
		: "raxport-candidates";
}

static std::string enabledPtmSummary()
{
	std::map<std::string, std::string> enabledPtms;
	if (!ProNovoConfig::getPTMinfo(enabledPtms) || enabledPtms.empty())
	{
		return "none";
	}

	std::ostringstream out;
	bool first = true;
	std::set<std::string> reportedKeys;
	const auto append = [&](const std::string &text)
	{
		if (!first)
		{
			out << ", ";
		}
		out << text;
		first = false;
	};
	for (const ProNovoConfig::PtmDefinition &definition :
		 ProNovoConfig::getPtmCatalog())
	{
		for (const auto &ptm : enabledPtms)
		{
			if (ptm.first.substr(0, 1) != definition.token)
			{
				continue;
			}
			std::ostringstream item;
			item << definition.name << " [token " << std::quoted(definition.token)
				 << ", sites " << ptm.second << "]";
			append(item.str());
			reportedKeys.insert(ptm.first);
			break;
		}
	}
	for (const auto &ptm : enabledPtms)
	{
		if (reportedKeys.find(ptm.first) == reportedKeys.end())
		{
			append(ptm.first + " [sites " + ptm.second + "]");
		}
	}
	return out.str();
}

static std::string enabledFixedPtmSummary()
{
	const std::vector<std::string> enabled =
		ProNovoConfig::getEnabledFixedPtmNames();
	if (enabled.empty())
	{
		return "none";
	}

	std::ostringstream out;
	for (size_t index = 0; index < enabled.size(); ++index)
	{
		if (index > 0)
		{
			out << ", ";
		}
		out << enabled[index];
	}
	return out.str();
}

static void printPtmCatalog(std::ostream &out)
{
	out << "Built-in fixed PTMs:\n";
	out << "  " << std::left << std::setw(24) << "Name"
		<< std::setw(10) << "Sites" << std::right << std::setw(14)
		<< "Mass shift" << "  Notes\n";
	for (const ProNovoConfig::FixedPtmDefinition &ptm :
		 ProNovoConfig::getFixedPtmCatalog())
	{
		std::ostringstream notes;
		notes << ptm.description;
		if (ptm.profileDefault)
		{
			notes << " [profile default]";
		}
		out << "  " << std::left << std::setw(24) << ptm.name
			<< std::setw(10) << ptm.sites << std::right << std::setw(14)
			<< std::fixed << std::setprecision(6)
			<< ptm.externalMonoisotopicShift << "  " << notes.str() << "\n";
	}
	out << "\n";

	out << "Built-in variable PTMs:\n";
	out << "  " << std::left << std::setw(36) << "Name"
		<< std::setw(8) << "Token" << std::setw(10) << "Sites"
		<< std::right << std::setw(14) << "Mass shift" << "  Notes\n";
	for (const ProNovoConfig::PtmDefinition &ptm : ProNovoConfig::getPtmCatalog())
	{
		std::ostringstream notes;
		notes << ptm.description;
		if (ptm.regularDefault)
		{
			notes << " [Regular default]";
		}
		if (!ptm.selectable)
		{
			notes << " [not selectable]";
		}
		out << "  " << std::left << std::setw(36) << ptm.name
			<< std::setw(8) << ptm.token << std::setw(10) << ptm.sites
			<< std::right << std::setw(14) << std::fixed << std::setprecision(6)
			<< ptm.externalMonoisotopicShift << "  " << notes.str() << "\n";
	}
	out << "\nUse names with --fixed-ptm and --ptm. Variable-PTM tokens are also\n"
		<< "accepted when quoted for the shell.\n";
}


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

bool TextUtils::parseSipAtomSpec(const std::string &spec,
                                     char &sipAtom,
                                     int &sipIsotopeMassNumber)
{
    std::string value = TextUtils::trim(spec);
    for (char &c : value)
        c = static_cast<char>(
            std::toupper(static_cast<unsigned char>(c)));

    if (value == "C13")
    {
        sipAtom = 'C';
        sipIsotopeMassNumber = 13;
    }
    else if (value == "H2")
    {
        sipAtom = 'H';
        sipIsotopeMassNumber = 2;
    }
    else if (value == "N15")
    {
        sipAtom = 'N';
        sipIsotopeMassNumber = 15;
    }
    else if (value == "O18")
    {
        sipAtom = 'O';
        sipIsotopeMassNumber = 18;
    }
    else if (value == "S34")
    {
        sipAtom = 'S';
        sipIsotopeMassNumber = 34;
    }
    else
    {
        return false;
    }
    return true;
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
	if (!std::isfinite(stepPct) || stepPct <= 0.0)
	{
		throw std::runtime_error("SIP step must be a finite positive number");
	}
	const auto validatePct = [&](double value)
	{
		if (!std::isfinite(value) || value < 0.0 || value > 100.0)
		{
			throw std::runtime_error(
				"SIP enrichment must be within [0,100]: " + rangeSpec);
		}
	};
	const size_t dash = rangeSpec.find('-');
	double lower = 0.0;
	double upper = 0.0;
	if (dash == std::string::npos)
	{
		if (!parseDouble(rangeSpec, lower))
		{
			throw std::runtime_error("Invalid SIP pct value: " + rangeSpec);
		}
		validatePct(lower);
		return {lower};
	}
	if (!parseDouble(rangeSpec.substr(0, dash), lower) || !parseDouble(rangeSpec.substr(dash + 1), upper))
	{
		throw std::runtime_error("Invalid SIP range: " + rangeSpec);
	}
	validatePct(lower);
	validatePct(upper);
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

static std::string uniquePeptideKey(const ScoredPsmRow &row)
{
	return row.nakedPeptide.empty() ? row.identifiedPeptide : row.nakedPeptide;
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
			if (!selectedPeptides.insert(uniquePeptideKey(row)).second)
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
						  const std::vector<ScoredPsmRow> &rows)
{
	sipPSM psm;
	psm.fileName = sampleName;
	const size_t count = rows.size();
	psm.fileNames.reserve(count);
	psm.scanNumbers.reserve(count);
	psm.precursorScanNumbers.reserve(count);
	psm.parentCharges.reserve(count);
	psm.isolationWindowCenterMZs.reserve(count);
	psm.measuredParentMasses.reserve(count);
	psm.calculatedParentMasses.reserve(count);
	psm.precursorNeutronMasses.reserve(count);
	psm.MS2IsotopicAbundances.reserve(count);
	psm.ranks.reserve(count);
	psm.scores.reserve(count);
	psm.identifiedPeptides.reserve(count);
	psm.originalPeptides.reserve(count);
	psm.nakePeptides.reserve(count);
	psm.proteinNames.reserve(count);
	psm.retentionTimes.reserve(count);
	psm.MVHscores.reserve(count);
	psm.XcorrScores.reserve(count);
	psm.WDPscores.reserve(count);
	psm.isDecoys.reserve(count);
	for (const ScoredPsmRow &row : rows)
	{
		psm.fileNames.push_back(fileNameForSpecId(sampleName, row.searchName, row.isDecoy));
		psm.scanNumbers.push_back(row.scanNumber);
		psm.precursorScanNumbers.push_back(row.precursorScanNumber);
		psm.parentCharges.push_back(row.parentCharge);
		psm.isolationWindowCenterMZs.push_back(row.isolationWindowCenterMZ);
		psm.measuredParentMasses.push_back(row.measuredParentMass);
		psm.calculatedParentMasses.push_back(row.calculatedParentMass);
		psm.precursorNeutronMasses.push_back(row.precursorNeutronMass);
		psm.MS2IsotopicAbundances.push_back(row.ms2IsotopicAbundancePct);
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

void SiprosSearchRunner::printUsage(std::ostream &out, const std::string &prog)
{
	out << "Usage:\n";
	out << "  " << prog << " -f sample.h5 -fasta proteins.fasta -o out "
		<< "[-a C13 -b 1-5 -s 1] [--pin-label 1|-1] [--pin-output name.pin]\n";
	out << "  " << prog << " -w hdf5_directory -fasta proteins.fasta -o out "
		<< "[-a C13 -b 1-5 -s 1] [--pin-label 1|-1]\n";
	out << "  " << prog << " -fasta proteins.fasta --fragment-index-cache cache.sfi "
		<< "--prepare-only\n";
	out << "  " << prog << " --list-ptms\n\n";
	out << "Parameters:\n";
	out << "  -f <sample.h5>              Raxport HDF5 scan file; repeat for multiple files\n";
	out << "  -w <directory>              directory of Raxport HDF5 scan files\n";
	out << "  -fasta <proteins.fasta>     one FASTA database; target/decoy orchestration is external\n";
	out << "  -o <directory>              output directory, default: out\n";
	out << "  --pin-output <name.pin>      PIN filename for a single -f input; default: <sample>.pin\n";
	out << "  --fragment-index-cache <p>  load/build a reusable v5 .sfi peptide cache\n";
	out << "  --rebuild-fragment-index    ignore an existing cache and rebuild it\n";
	out << "  --prepare-only              prepare a Regular peptide cache without searching H5 files\n";
	out << "  --precursor-source <mode>    profile policy: Regular=ms1-neighborhood, SIP=raxport-candidates\n";
	out << "                              Regular MS1 mode uses linked parent +/-2 scans and peak charge\n";
	out << "  -a <SIP atom/isotope>       SIP isotope: C13, H2, N15, O18, or S34\n";
	out << "  -b <pct|lower-upper>        SIP percentage or inclusive range\n";
	out << "  -s, --step <pct>            SIP percentage step\n";
	out << "  --tolerance-ms1 <Da>        parent mass tolerance (default: 0.01 Da)\n";
	out << "  --tolerance-ms2 <Da>        fragment mass tolerance (default: 0.01 Da)\n";
	out << "  --fixed-ptm <name>          exact fixed-PTM selection; repeat to select several\n";
	out << "                              special values: default, none, all\n";
	out << "  --ptm <name>                exact variable-PTM selection; repeat to select several\n";
	out << "                              use names when possible; special values: default, none, all\n";
	out << "  --max-ptm-count <N>         maximum variable PTMs per peptide\n";
	out << "  --list-ptms                 list built-in PTM names, tokens, and sites, then exit\n";
	out << "  --pin-label <1|-1>          label written to PIN rows; default: 1\n";
	out << "  --top-psms-per-scan <N>     final WDP-ranked PSMs retained per scan across SIP pct\n";
	out << "\nWithout --fixed-ptm, compiled fixed CAM is used. Without --ptm, the compiled\n"
		<< "profile defaults are used (Regular: oxidation and deamidation; SIP: none).\n";
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
		else if (option == "-f")
		{
			args.scanFiles.push_back(requireValue(i, option));
		}
		else if (option == "-fasta")
		{
			args.fastaFile = requireValue(i, option);
		}
		else if (option == "-o")
		{
			args.outputDirectory = requireValue(i, option);
		}
		else if (option == "--pin-output")
		{
			args.pinOutputFile = requireValue(i, option);
		}
		else if (option == "--fragment-index-cache")
		{
			args.fragmentIndexCache = requireValue(i, option);
		}
		else if (option == "--rebuild-fragment-index")
		{
			args.rebuildFragmentIndex = true;
		}
		else if (option == "--prepare-only")
		{
			args.prepareOnly = true;
		}
		else if (option == "--precursor-source")
		{
			args.precursorSourceProvided = true;
			const std::string value = TextUtils::toLower(
				TextUtils::trim(requireValue(i, option)));
			if (!valid)
			{
				return false;
			}
			if (value == "raxport-candidates")
			{
				args.raxportReadOptions.precursorSource =
					PrecursorSource::RaxportCandidates;
			}
			else if (value == "ms1-neighborhood")
			{
				args.raxportReadOptions.precursorSource =
					PrecursorSource::Ms1Neighborhood;
			}
			else
			{
				err << "--precursor-source must be raxport-candidates or ms1-neighborhood\n";
				return false;
			}
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
			args.sipStepProvided = true;
		}
		else if (option == "--tolerance-ms1" || option == "--tolerance-ms2")
		{
			double tolerance = 0.0;
			const std::string value = requireValue(i, option);
			if (valid && (!parseDouble(value, tolerance) || tolerance <= 0.0))
			{
				err << option << " must be a positive number in Da\n";
				return false;
			}
			if (option == "--tolerance-ms1")
			{
				args.toleranceMs1Da = tolerance;
				args.toleranceMs1Provided = true;
			}
			else
			{
				args.toleranceMs2Da = tolerance;
				args.toleranceMs2Provided = true;
			}
		}
		else if (option == "--fixed-ptm" || option == "--ptm")
		{
			const std::string value = requireValue(i, option);
			if (!valid)
			{
				return false;
			}
			const std::string selector = TextUtils::trim(value);
			if (selector.empty())
			{
				err << option << " requires a non-empty name";
				if (option == "--ptm")
				{
					err << ", token";
				}
				err << ", or special value\n";
				return false;
			}
			if (option == "--fixed-ptm")
			{
				args.fixedPtmSelectors.push_back(selector);
			}
			else
			{
				args.ptmSelectors.push_back(selector);
			}
		}
		else if (option == "--max-ptm-count")
		{
			const std::string value = requireValue(i, option);
			if (valid && (!parseInt(value, args.maxPtmCountOverride) ||
				args.maxPtmCountOverride < 0))
			{
				err << "--max-ptm-count must be a non-negative integer\n";
				return false;
			}
		}
		else if (option == "--list-ptms")
		{
			args.listPtms = true;
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
		else if ((option == "-h") || (option == "--help"))
		{
			args.showHelp = true;
			printUsage(out, values.empty() ? "sipros search-fasta" : values[0]);
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
	if (args.listPtms)
	{
		printPtmCatalog(out);
		args.showHelp = true;
		return true;
	}
	if (args.fastaFile.empty())
	{
		err << "One FASTA database is required with -fasta\n";
		return false;
	}
	const bool anySipControl = !args.sipElementSpec.empty() ||
		!args.sipRangeSpec.empty() || args.sipStepProvided;
	const bool allSipControls = !args.sipElementSpec.empty() &&
		!args.sipRangeSpec.empty() && args.sipStepProvided;
	if (anySipControl && !allSipControls)
	{
		err << "SIP search requires -a, -b, and -s together; no values are inferred\n";
		return false;
	}
	if (allSipControls)
	{
		if (!args.fragmentIndexCache.empty() ||
			args.rebuildFragmentIndex)
		{
			err << "SIP FASTA search uses only legacy H5 precursor candidates; "
				<< "fragment-index options are not allowed\n";
			return false;
		}
		if (args.precursorSourceProvided &&
			args.raxportReadOptions.precursorSource !=
				PrecursorSource::RaxportCandidates)
		{
			err << "SIP FASTA search requires --precursor-source raxport-candidates\n";
			return false;
		}
		args.raxportReadOptions.precursorSource =
			PrecursorSource::RaxportCandidates;
	}
	else
	{
		if (args.precursorSourceProvided &&
			args.raxportReadOptions.precursorSource !=
				PrecursorSource::Ms1Neighborhood)
		{
			err << "Regular FASTA search requires --precursor-source ms1-neighborhood\n";
			return false;
		}
		args.raxportReadOptions.precursorSource =
			PrecursorSource::Ms1Neighborhood;
	}
	if (args.rebuildFragmentIndex && args.fragmentIndexCache.empty())
	{
		err << "--rebuild-fragment-index requires --fragment-index-cache\n";
		return false;
	}
	if (args.prepareOnly)
	{
		if (allSipControls)
		{
			err << "--prepare-only is available only for Regular FASTA peptide caches\n";
			return false;
		}
		if (args.fragmentIndexCache.empty())
		{
			err << "--prepare-only requires --fragment-index-cache\n";
			return false;
		}
		if (!args.workingDirectory.empty() || !args.scanFiles.empty())
		{
			err << "--prepare-only does not accept -f or -w scan inputs\n";
			return false;
		}
		if (!args.pinOutputFile.empty())
		{
			err << "--prepare-only does not write --pin-output\n";
			return false;
		}
	}
	else if (args.workingDirectory.empty() && args.scanFiles.empty())
	{
		args.workingDirectory = ".";
	}
	if (!args.workingDirectory.empty() && !args.scanFiles.empty())
	{
		err << "Use either repeatable -f inputs or one -w directory, not both\n";
		return false;
	}
	if (args.prepareOnly)
	{
		// Cache preparation intentionally has no scan input.
	}
	else if (!args.scanFiles.empty())
	{
		for (const std::string &scanFile : args.scanFiles)
		{
			if (!isHdf5Path(scanFile))
			{
				err << "Raxport HDF5 scan input required (.h5 or .hdf5); unsupported file: "
					<< scanFile << "\n";
				return false;
			}
		}
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
	if (!args.pinOutputFile.empty())
	{
		const fs::path pinOutput(args.pinOutputFile);
		if (pinOutput.has_parent_path() ||
			pinOutput.filename().string() != args.pinOutputFile ||
			pinOutput.extension() != ".pin")
		{
			err << "--pin-output must be a .pin filename within the output directory\n";
			return false;
		}
		if (args.scanFiles.size() != 1 || !args.workingDirectory.empty())
		{
			err << "--pin-output requires exactly one -f input file\n";
			return false;
		}
	}
	if (!args.fragmentIndexCache.empty())
	{
		const fs::path cachePath(args.fragmentIndexCache);
		if (cachePath.extension() != ".sfi")
		{
			err << "--fragment-index-cache must use the .sfi extension\n";
			return false;
		}
	}
	return true;
}

int SiprosSearchRunner::prepare(const DatabaseSearchArguments &args)
{
	if (prepared_)
	{
		return preparedSipMode_ == !args.sipElementSpec.empty() ? 0 : 1;
	}

	preparedSipMode_ = !args.sipElementSpec.empty();
	if (!ProNovoConfig::load(preparedSipMode_
			? ProNovoConfig::Profile::Sip
			: ProNovoConfig::Profile::Regular))
	{
		std::cerr << "Could not initialize the built-in Sipros profile.\n";
		return 1;
	}
	std::string ptmError;
	if (!ProNovoConfig::configureFixedPtms(args.fixedPtmSelectors, ptmError))
	{
		std::cerr << "Invalid fixed PTM selection: " << ptmError << "\n";
		return 1;
	}
	if (!ProNovoConfig::configureVariablePtms(
			args.ptmSelectors, args.maxPtmCountOverride, ptmError))
	{
		std::cerr << "Invalid PTM selection: " << ptmError << "\n";
		return 1;
	}
	ProNovoConfig::setFASTAfilename(args.fastaFile);
	if (args.toleranceMs1Provided || args.toleranceMs2Provided)
	{
		ProNovoConfig::setMassAccuracy(
			args.toleranceMs1Provided
				? args.toleranceMs1Da
				: ProNovoConfig::getMassAccuracyParentIon(),
			args.toleranceMs2Provided
				? args.toleranceMs2Da
				: ProNovoConfig::getMassAccuracyFragmentIon());
	}

	if (preparedSipMode_)
	{
		std::string chemistryError;
		if (!ProNovoConfig::validatePreparationChemistry(
				ProNovoConfig::configIsotopologue, chemistryError))
		{
			std::cerr << chemistryError << "\n";
			return 1;
		}
	}
	else
	{
		fragmentIndex_ = std::make_shared<FragmentIndex>();
		std::string indexError;
		if (!fragmentIndex_->loadOrBuild(
				args.fragmentIndexCache,
				args.rebuildFragmentIndex,
				indexError))
		{
			std::cerr << "Fragment-index preparation failed: "
					  << indexError << "\n";
			fragmentIndex_.reset();
			return 1;
		}
		const FragmentIndexStats &stats = fragmentIndex_->stats();
		std::cout << "\nRegular-search peptide cache\n"
				  << "  File    : "
				  << (args.fragmentIndexCache.empty()
						  ? "memory only"
						  : args.fragmentIndexCache)
				  << "\n"
				  << "  Status  : "
				  << (stats.loadedFromCache ? "loaded existing cache" : "generated cache")
				  << "\n"
				  << "  Contents: "
				  << formatPerformanceCount(fragmentIndex_->peptideCount())
				  << " peptides, "
				  << formatPerformanceCount(fragmentIndex_->fragmentCount())
				  << " fragment postings, "
				  << formatPerformanceCount(fragmentIndex_->fragmentBinCount())
				  << " occupied fragment bins, "
				  << formatPerformanceCount(fragmentIndex_->precursorBlockCount())
				  << " precursor blocks";
		if (stats.cacheBytes > 0)
		{
			std::cout << ", " << formatMebibytes(stats.cacheBytes);
		}
		std::cout << "\n";
		printPerformanceHeader(
			std::cout, "Peptide-cache timing", omp_get_max_threads());
		if (stats.loadedFromCache)
		{
			printPerformanceStage(
				std::cout, "Load peptide cache",
				performanceTiming(stats.loadSeconds, stats.loadCpuSeconds));
		}
		else
		{
			printPerformanceStage(
				std::cout, "Read FASTA",
				performanceTiming(stats.fastaParseSeconds,
					stats.fastaParseCpuSeconds));
			printPerformanceStage(
				std::cout, "Count digested peptides",
				performanceTiming(stats.digestCountSeconds,
					stats.digestCountCpuSeconds));
			printPerformanceStage(
				std::cout, "Generate peptide records",
				performanceTiming(stats.digestFillSeconds,
					stats.digestFillCpuSeconds));
			printPerformanceStage(
				std::cout, "Calculate precursor masses",
				performanceTiming(stats.precursorSeconds,
					stats.precursorCpuSeconds));
			printPerformanceStage(
				std::cout, "Generate fragment postings",
				performanceTiming(stats.fragmentSeconds,
					stats.fragmentCpuSeconds));
			printPerformanceStage(
				std::cout, "Sort fragment postings",
				performanceTiming(stats.sortSeconds,
					stats.sortCpuSeconds));
			if (stats.saveSeconds > 0.0)
			{
				printPerformanceStage(
					std::cout, "Save peptide cache",
					performanceTiming(stats.saveSeconds,
						stats.saveCpuSeconds));
			}
			if (stats.mapSeconds > 0.0)
			{
				printPerformanceStage(
					std::cout, "Map generated peptide cache",
					performanceTiming(stats.mapSeconds,
						stats.mapCpuSeconds));
			}
			printPerformanceStage(
				std::cout, "Generate peptide cache (total)",
				performanceTiming(stats.generateSeconds,
					stats.generateCpuSeconds));
		}
		printPerformanceFooter(std::cout);
	}

	prepared_ = true;
	return 0;
}

int SiprosSearchRunner::runScan(const std::string &scanFile,
								  const DatabaseSearchArguments &args)
{
	const bool directSipMode = !args.sipElementSpec.empty();
	if (!prepared_ && prepare(args) != 0)
	{
		return 1;
	}
	if (preparedSipMode_ != directSipMode)
	{
		std::cerr << "Prepared search profile does not match this scan request.\n";
		return 1;
	}
	const int topKeep = args.topPsmsPerScan > 0 ? args.topPsmsPerScan : ProNovoConfig::INTTOPKEEP;
	const bool isDecoyLabel = args.pinLabel < 0;
	char sipAtom = 0;
	int sipIsotopeMassNumber = -1;
	if (directSipMode)
	{
		if (!TextUtils::parseSipAtomSpec(
				args.sipElementSpec, sipAtom, sipIsotopeMassNumber))
		{
			std::cerr << "Unsupported SIP element/isotope: "
					  << args.sipElementSpec << "\n";
			return 1;
		}
	}
	fs::create_directories(args.outputDirectory);

	MS2ScanVector scanVector(scanFile, args.outputDirectory, args.raxportReadOptions);
	scanVector.configureFragmentIndex(fragmentIndex_);
	const PerformanceTimer fileTimer;
	if (directSipMode)
	{
		std::cout << "\nSipros FASTA SIP search\n";
		std::cout << "  Scan file : " << scanFile << "\n";
		std::cout << "  FASTA     : " << ProNovoConfig::getFASTAfilename() << (isDecoyLabel ? " (decoy)" : " (target)") << "\n";
		std::cout << "  Profile   : built-in SIP\n";
		std::cout << "  SIP       : " << args.sipElementSpec << " " << args.sipRangeSpec << " step " << args.sipStepPct << "\n";
		std::cout << "  Fixed PTMs: " << enabledFixedPtmSummary() << "\n";
		std::cout << "  Var PTMs  : " << enabledPtmSummary()
				  << " (max " << ProNovoConfig::getMaxPTMcount() << " per peptide)\n";
		std::cout << "  Precursors: "
				  << precursorSourceName(args.raxportReadOptions.precursorSource) << "\n";
		std::cout << "  TopN      : " << topKeep << " unique peptides per scan across SIP pct\n";
	}
	else
	{
		std::cout << "\nRegular FASTA search\n"
				  << "  Scan file : " << scanFile << "\n"
				  << "  FASTA     : " << ProNovoConfig::getFASTAfilename()
				  << (isDecoyLabel ? " (decoy)" : " (target)") << "\n"
				  << "  Profile   : built-in Regular\n"
				  << "  Fixed PTMs: " << enabledFixedPtmSummary() << "\n"
				  << "  Var PTMs  : " << enabledPtmSummary()
				  << " (max " << ProNovoConfig::getMaxPTMcount() << " per peptide)\n"
				  << "  Precursors: "
				  << precursorSourceName(args.raxportReadOptions.precursorSource) << "\n"
				  << "  Search    : lossless peptide/fragment cache\n";
		if (!args.fragmentIndexCache.empty())
		{
			std::cout << "  Cache     : " << args.fragmentIndexCache << "\n";
		}
		printPerformanceHeader(
			std::cout,
			"Search timing: " + fs::path(scanFile).filename().string(),
			omp_get_max_threads());
	}

	const PerformanceTimer loadTimer;
	bool loadedScans = scanVector.loadMassData();
	if (!loadedScans)
	{
		std::cerr << "Error: Failed to load file: " << scanFile << "\n";
		return 1;
	}
	if (!directSipMode)
	{
		std::ostringstream loadDetail;
		loadDetail << formatPerformanceCount(scanVector.scanCount())
				   << " MS2 scans, "
				   << formatPerformanceCount(
						  scanVector.precursorHypothesisCount())
				   << " precursor hypotheses";
		printPerformanceStage(
			std::cout, "Load HDF5 scans", loadTimer.elapsed(), loadDetail.str());
	}

	std::vector<ScoredPsmRow> scoredRows;
	if (directSipMode)
	{
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
		for (size_t pctIndex = 0; pctIndex < pctValues.size(); ++pctIndex)
		{
			const double pct = pctValues[pctIndex];
			scanVector.clearSearchResults();
			if (!ProNovoConfig::applySipAbundance(sipAtom, pct / 100.0))
			{
				std::cerr << "Could not apply SIP abundance for " << args.sipElementSpec << " at " << pct << " percent\n";
				return 1;
			}
			ProNovoConfig::setSearchName(formatSipSearchName(sipSpec, pct));
			std::cout << "  [" << (pctIndex + 1) << "/" << pctValues.size() << "] "
					  << ProNovoConfig::getSearchName() << "\n";
			const double pctBegin = omp_get_wtime();
			scanVector.startProcessingWdpSip();
			scanVector.appendScoredPsmRows(scoredRows, isDecoyLabel, topKeep, pct);
			pruneScoredRowsByScan(scoredRows, topKeep);
			std::cout << "  Result: done in " << formatElapsedSeconds(omp_get_wtime() - pctBegin)
					  << ", cumulative retained rows: " << scoredRows.size() << "\n";
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
		const PerformanceTimer collectTimer;
		scanVector.appendScoredPsmRows(scoredRows, isDecoyLabel, topKeep, 1.07);
		pruneScoredRowsByScan(scoredRows, topKeep);
		printPerformanceStage(
			std::cout, "Collect top PSMs", collectTimer.elapsed(),
			formatPerformanceCount(scoredRows.size()) + " retained rows");
	}

	const PerformanceTimer pinTotalTimer;
	const std::string sampleName = fs::path(scanFile).stem().string();
	const PerformanceTimer pinBuildTimer;
	sipPSM psm = buildSipPsm(sampleName, scoredRows);
	const PerformanceTiming pinBuildTiming = pinBuildTimer.elapsed();
	PSMfeatureExtractor extractor;
	extractor.sipPSMs.push_back(std::move(psm));
	PerformanceTiming pinMs1LoadTiming;
	PerformanceTiming pinFeatureTiming;
	if (!extractor.sipPSMs.front().scanNumbers.empty())
	{
		const PerformanceTimer pinMs1LoadTimer;
		extractor.ms1Data = scanVector.releaseMs1Data();
		pinMs1LoadTiming = pinMs1LoadTimer.elapsed();
		const PerformanceTimer pinFeatureTimer;
		extractor.mSipPSM = &extractor.sipPSMs.front();
		extractor.initializeFeatureVectors(extractor.sipPSMs.front());
		extractor.extractFeaturesOfEachPSM();
		pinFeatureTiming = pinFeatureTimer.elapsed();
	}
	std::string pinFileName = args.pinOutputFile;
	if (pinFileName.empty())
	{
		pinFileName = sampleName + ".pin";
	}
	const fs::path pinPath = fs::path(args.outputDirectory) / pinFileName;
	const PerformanceTimer pinWriteTimer;
	PinWriter::writePecorlatorPin(pinPath.string(), extractor.sipPSMs, false);
	const PerformanceTiming pinWriteTiming = pinWriteTimer.elapsed();
	if (!directSipMode)
	{
		printPerformanceStage(
			std::cout, "Build PIN columns", pinBuildTiming,
			formatPerformanceCount(scoredRows.size()) + " retained rows");
		printPerformanceStage(
			std::cout, "Reuse loaded MS1 for PIN", pinMs1LoadTiming,
			formatPerformanceCount(scoredRows.size()) + " retained rows");
		printPerformanceStage(
			std::cout, "Calculate PIN features", pinFeatureTiming,
			formatPerformanceCount(scoredRows.size()) + " retained rows");
		printPerformanceStage(
			std::cout, "Format and write PIN", pinWriteTiming,
			pinPath.string());
		printPerformanceStage(
			std::cout, "PIN output (total)", pinTotalTimer.elapsed(),
			pinPath.string());
		printPerformanceStage(
			std::cout, "Search file (total)", fileTimer.elapsed());
		printPerformanceFooter(std::cout);
	}
	else
	{
		std::cout << "Wrote PIN: " << pinPath.string() << " ("
				  << scoredRows.size()
				  << " scored rows retained after top-N pruning)\n";
	}
	return 0;
}

} // namespace sipros
