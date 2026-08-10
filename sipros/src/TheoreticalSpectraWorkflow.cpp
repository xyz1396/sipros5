#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <omp.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <H5Cpp.h>
#include "SiprosWorkflows.h"
#include "SiprosSearchRunner.h"
#include "isotopologue.h"
#include "proNovoConfig.h"
#include "spectraindex.h"

namespace fs = std::filesystem;

namespace
{
constexpr int kSpectraHdf5FormatVersion = 2;

enum class SipAbundanceMode
{
	InputRow,
	FixedUser
};

struct Args
{
	std::string inputPath;
	std::string outputPath;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	SipAbundanceMode sipAbundanceMode = SipAbundanceMode::InputRow;
	double fixedSipAbundancePct = 0.0;
	double probCutoff = 0.01;
	int threads = 0;
	std::vector<std::string> fixedPtmSelectors;
};

struct PsmRow
{
	std::string psmId;
	double sipPct = 0.0;
	std::string peptide;
	int precursorCharge = 1;
	std::string retention = "0";
	std::string proteins = "{UNKNOWN}";
	double probability = 0.0;
	double svmScore = -std::numeric_limits<double>::infinity();
	size_t order = 0;
};

void printUsage(const char *prog)
{
	std::cerr << "Usage: " << prog
			  << " -i <input.tsv|report_dir> -o <output.sfi> -a <SIP atom/isotope, e.g. C13,H2,O18,N15,S34> [-b <pct>] [-p <prob cutoff>] [--fixed-ptm <name|default|none|all>] [-t <threads>]\n";
	std::cerr << "Required columns in input: (PSMId, SpecId, or Spectrum), Peptide";
	std::cerr << " [MS2IsotopicAbundances required unless -b/--sip-abundance is set]\n";
	std::cerr << "Directory input recursively combines files named psm.tsv. Output is a memory-mapped SIP fragment index (.sfi).\n";
	std::cerr << "--fixed-ptm is repeatable; omit it to use the compiled default (carbamidomethyl C).\n";
}

bool parseArgs(int argc, char **argv, Args &args)
{
	for (int i = 1; i < argc; ++i)
	{
		const std::string opt = argv[i];
		if (opt == "-h" || opt == "--help")
		{
			printUsage(argv[0]);
			return false;
		}
		if (opt == "-i")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			args.inputPath = argv[++i];
		}
		else if (opt == "-o")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			args.outputPath = argv[++i];
		}
		else if (opt == "-a" || opt == "--sip-atom")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			const std::string sipSpec = sipros::TextUtils::trim(argv[++i]);
			if (!sipros::TextUtils::parseSipAtomSpec(sipSpec, args.sipAtom, args.sipIsotopeMassNumber))
			{
				std::cerr << "Invalid SIP atom/isotope: " << sipSpec
						  << ". Use C13,H2,N15,O18,S34.\n";
				return false;
			}
		}
		else if (opt == "-b" || opt == "--sip-abundance")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing numeric value for option: " << opt << "\n";
				return false;
			}
			const std::string abundanceSpec = sipros::TextUtils::trim(argv[++i]);
			try
			{
				args.fixedSipAbundancePct = std::stod(abundanceSpec);
				args.sipAbundanceMode = SipAbundanceMode::FixedUser;
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid SIP abundance value: " << abundanceSpec
						  << ". Provide a percentage in [0, 100].\n";
				return false;
			}
			if (args.fixedSipAbundancePct < 0.0 || args.fixedSipAbundancePct > 100.0)
			{
				std::cerr << "SIP abundance percentage must be in [0, 100].\n";
				return false;
			}
		}
		else if (opt == "-p" || opt == "--prob-cutoff")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			try
			{
				args.probCutoff = std::stod(argv[++i]);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid probability cutoff value.\n";
				return false;
			}
		}
		else if (opt == "--fixed-ptm")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			const std::string selector = sipros::TextUtils::trim(argv[++i]);
			if (selector.empty())
			{
				std::cerr << "Fixed PTM selector must not be empty.\n";
				return false;
			}
			args.fixedPtmSelectors.push_back(selector);
		}
		else if (opt == "-t" || opt == "--threads")
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << opt << "\n";
				return false;
			}
			try
			{
				args.threads = std::stoi(argv[++i]);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid thread count value.\n";
				return false;
			}
		}
		else
		{
			std::cerr << "Unknown option: " << opt << "\n";
			return false;
		}
	}

	if (args.inputPath.empty() || args.outputPath.empty() || args.sipAtom == '\0')
	{
		printUsage(argv[0]);
		return false;
	}
	if (args.probCutoff < 0.0 || args.probCutoff > 1.0)
	{
		std::cerr << "Probability cutoff must be in [0, 1].\n";
		return false;
	}
	if (args.threads < 0)
	{
		std::cerr << "Thread count must be >= 0.\n";
		return false;
	}
	return true;
}

double parseFirstDouble(const std::string &s)
{
	const std::string t = sipros::TextUtils::trim(s);
	if (t.empty())
	{
		throw std::runtime_error("Empty MS2IsotopicAbundances value.");
	}
	size_t consumed = 0;
	const double v = std::stod(t, &consumed);
	if (consumed == 0)
	{
		throw std::runtime_error("Invalid MS2IsotopicAbundances value: " + s);
	}
	return v;
}

bool parseDoubleField(const std::string &text, double &value)
{
	try
	{
		const std::string trimmed = sipros::TextUtils::trim(text);
		if (trimmed.empty())
		{
			return false;
		}
		size_t consumed = 0;
		value = std::stod(trimmed, &consumed);
		return consumed > 0;
	}
	catch (const std::exception &)
	{
		return false;
	}
}

bool parseIntField(const std::string &text, int &value)
{
	try
	{
		const std::string trimmed = sipros::TextUtils::trim(text);
		if (trimmed.empty())
		{
			return false;
		}
		size_t consumed = 0;
		value = std::stoi(trimmed, &consumed);
		return consumed > 0;
	}
	catch (const std::exception &)
	{
		return false;
	}
}

bool convertModifiedPeptide(const std::string &plainPeptide,
							const std::string &modifiedPeptide,
							const std::string &assignedModifications,
							std::string &normalizedPeptide,
							std::string &reason)
{
	return ProNovoConfig::translatePsmPeptide(
		plainPeptide, modifiedPeptide, assignedModifications,
		normalizedPeptide, reason);
}

std::string normalizePeptide(const std::string &raw)
{
	const std::string p = sipros::TextUtils::trim(raw);
	if (p.empty())
	{
		return p;
	}

	const size_t lb = p.find('[');
	const size_t rb = p.rfind(']');
	if (lb != std::string::npos && rb != std::string::npos && lb < rb)
	{
		return p.substr(lb, rb - lb + 1);
	}
	if (!p.empty() && p.front() == '[' && p.back() == ']')
	{
		return p;
	}

	const size_t dot1 = p.find('.');
	if (dot1 != std::string::npos)
	{
		const size_t dot2 = p.find('.', dot1 + 1);
		if (dot2 != std::string::npos && dot2 > dot1 + 1)
		{
			return "[" + p.substr(dot1 + 1, dot2 - dot1 - 1) + "]";
		}
	}

	return "[" + p + "]";
}

size_t optionalColumn(const std::unordered_map<std::string, size_t> &columns,
					  const std::vector<std::string> &names)
{
	for (const std::string &name : names)
	{
		const auto found = columns.find(name);
		if (found != columns.end())
		{
			return found->second;
		}
	}
	return std::string::npos;
}

std::vector<fs::path> collectInputFiles(const std::string &inputPath)
{
	const fs::path path(inputPath);
	if (!fs::exists(path))
	{
		throw std::runtime_error("Input path does not exist: " + inputPath);
	}

	std::vector<fs::path> files;
	if (fs::is_regular_file(path))
	{
		files.push_back(path);
	}
	else if (fs::is_directory(path))
	{
		for (const auto &entry : fs::recursive_directory_iterator(path))
		{
			if (entry.is_regular_file() && entry.path().filename() == "psm.tsv")
			{
				files.push_back(entry.path());
			}
		}
	}
	std::sort(files.begin(), files.end());
	if (files.empty())
	{
		throw std::runtime_error("No input TSV files found under: " + inputPath);
	}
	return files;
}

std::vector<std::string> splitProteinList(const std::string &proteins)
{
	std::string inner = sipros::TextUtils::trim(proteins);
	if (inner.size() >= 2 && inner.front() == '{' && inner.back() == '}')
	{
		inner = inner.substr(1, inner.size() - 2);
	}
	std::vector<std::string> output;
	std::stringstream stream(inner);
	std::string protein;
	while (std::getline(stream, protein, ','))
	{
		protein = sipros::TextUtils::trim(protein);
		if (!protein.empty())
		{
			output.push_back(protein);
		}
	}
	return output;
}

std::string formatProteinList(const std::vector<std::string> &proteins)
{
	std::string output = "{";
	for (const std::string &protein : proteins)
	{
		if (protein.empty())
		{
			continue;
		}
		if (output.size() > 1)
		{
			output += ',';
		}
		output += protein;
	}
	output += '}';
	return output;
}

std::string mergeProteinLists(const std::string &left, const std::string &right)
{
	std::set<std::string> seen;
	std::vector<std::string> merged;
	for (const std::string &protein : splitProteinList(left))
	{
		if (seen.insert(protein).second)
		{
			merged.push_back(protein);
		}
	}
	for (const std::string &protein : splitProteinList(right))
	{
		if (seen.insert(protein).second)
		{
			merged.push_back(protein);
		}
	}
	return merged.empty() ? "{UNKNOWN}" : formatProteinList(merged);
}

std::string combineProteins(const std::string &primary, const std::string &mapped)
{
	return mergeProteinLists(primary, mapped);
}

void readInputFile(const fs::path &path,
				   bool requireSipPct,
				   std::vector<PsmRow> &rows,
				   size_t &skippedUnsupportedMods)
{
	std::ifstream in(path);
	if (!in)
	{
		throw std::runtime_error("Cannot open input file: " + path.string());
	}

	std::string headerLine;
	if (!std::getline(in, headerLine))
	{
		throw std::runtime_error("Input file is empty: " + path.string());
	}
	const std::vector<std::string> headers = sipros::TextUtils::splitTab(headerLine);
	std::unordered_map<std::string, size_t> col;
	for (size_t i = 0; i < headers.size(); ++i)
	{
		col[sipros::TextUtils::trim(headers[i])] = i;
	}

	const size_t idxPSMId = optionalColumn(col, {"PSMId", "SpecId", "Spectrum"});
	if (idxPSMId == std::string::npos)
	{
		throw std::runtime_error(
			"Missing required ID column PSMId, SpecId, or Spectrum in " + path.string());
	}
	const size_t idxPeptide = optionalColumn(col, {"Peptide"});
	if (idxPeptide == std::string::npos)
	{
		throw std::runtime_error("Missing required column Peptide in " + path.string());
	}
	size_t idxMS2Pct = std::string::npos;
	if (requireSipPct)
	{
		idxMS2Pct = optionalColumn(col, {"MS2IsotopicAbundances"});
		if (idxMS2Pct == std::string::npos)
		{
			throw std::runtime_error(
				"Missing required column MS2IsotopicAbundances in " + path.string());
		}
	}
	const size_t idxCharge = optionalColumn(col, {"parentCharges", "Charge"});
	const size_t idxModifiedPeptide = optionalColumn(col, {"Modified Peptide"});
	const size_t idxAssignedModifications =
		optionalColumn(col, {"Assigned Modifications"});
	const size_t idxRetention = optionalColumn(
		col, {"Retention", "RetentionTime", "retentionTime"});
	const size_t idxProteins = optionalColumn(
		col, {"Proteins", "ProteinNames", "proteinNames", "ProteinName",
			  "proteinName", "Protein", "protein"});
	const size_t idxMappedProteins = optionalColumn(col, {"Mapped Proteins"});
	const size_t idxProbability = optionalColumn(col, {"Probability"});
	const size_t idxSvmScore = optionalColumn(col, {"SVMscore"});
	if (idxSvmScore == std::string::npos)
	{
		throw std::runtime_error(
			"Missing required column SVMscore in " + path.string());
	}
	const bool reportStyleInput = col.find("Spectrum") != col.end();

	std::string line;
	size_t lineNo = 1;
	while (std::getline(in, line))
	{
		++lineNo;
		if (line.empty())
		{
			continue;
		}
		const std::vector<std::string> f = sipros::TextUtils::splitTab(line);
		size_t need = std::max(idxPSMId, idxPeptide);
		if (requireSipPct)
		{
			need = std::max(need, idxMS2Pct);
		}
		for (const size_t optional : {idxCharge, idxModifiedPeptide,
									  idxAssignedModifications, idxRetention,
									  idxProteins, idxMappedProteins, idxProbability,
									  idxSvmScore})
		{
			if (optional != std::string::npos)
			{
				need = std::max(need, optional);
			}
		}
		if (f.size() <= need)
		{
			continue;
		}
		PsmRow row;
		row.psmId = sipros::TextUtils::trim(f[idxPSMId]);
		if (reportStyleInput)
		{
			const std::string modified = idxModifiedPeptide == std::string::npos
										 ? std::string()
										 : f[idxModifiedPeptide];
			const std::string assignedModifications =
				idxAssignedModifications == std::string::npos
					? std::string()
					: f[idxAssignedModifications];
			std::string reason;
			if (!convertModifiedPeptide(
					f[idxPeptide], modified, assignedModifications,
					row.peptide, reason))
			{
				++skippedUnsupportedMods;
				continue;
			}
		}
		else
		{
			row.peptide = normalizePeptide(f[idxPeptide]);
		}
		if (row.psmId.empty() || row.peptide.empty())
		{
			continue;
		}
		try
		{
			if (requireSipPct)
			{
				row.sipPct = parseFirstDouble(f[idxMS2Pct]);
			}
			if (idxCharge != std::string::npos)
			{
				int charge = 0;
				if (parseIntField(f[idxCharge], charge) && charge > 0)
				{
					row.precursorCharge = charge;
				}
			}
			if (idxRetention != std::string::npos)
			{
				const std::string retention = sipros::TextUtils::trim(f[idxRetention]);
				if (!retention.empty())
				{
					row.retention = retention;
				}
			}
			const std::string primary = idxProteins == std::string::npos
										? std::string()
										: f[idxProteins];
			const std::string mapped = idxMappedProteins == std::string::npos
									   ? std::string()
									   : f[idxMappedProteins];
			row.proteins = combineProteins(primary, mapped);
			if (idxProbability != std::string::npos)
			{
				(void)parseDoubleField(f[idxProbability], row.probability);
			}
			if (!parseDoubleField(f[idxSvmScore], row.svmScore))
			{
				throw std::runtime_error("invalid SVMscore");
			}
		}
		catch (const std::exception &)
		{
			if (requireSipPct)
			{
				std::cerr << "Skipping " << path << ':' << lineNo
						  << ": invalid MS2IsotopicAbundances or charge.\n";
			}
			else
			{
				std::cerr << "Skipping " << path << ':' << lineNo << ": invalid charge.\n";
			}
			continue;
		}
		row.order = rows.size();
		rows.push_back(std::move(row));
	}
}

std::vector<PsmRow> readInputRows(const std::string &inputPath,
								  bool requireSipPct,
								  size_t &inputFiles,
								  size_t &skippedUnsupportedMods)
{
	const std::vector<fs::path> paths = collectInputFiles(inputPath);
	inputFiles = paths.size();
	std::vector<PsmRow> rows;
	for (const fs::path &path : paths)
	{
		readInputFile(path, requireSipPct, rows, skippedUnsupportedMods);
	}
	return rows;
}

std::string peptideMassClassKey(const std::string &peptide)
{
	std::string key = peptide;
	for (char &residue : key)
	{
		if (residue == 'I' || residue == 'J')
		{
			residue = 'L';
		}
	}
	return key;
}

bool isBetterPsm(const PsmRow &candidate, const PsmRow &current)
{
	constexpr double epsilon = 1e-15;
	if (candidate.probability > current.probability + epsilon)
	{
		return true;
	}
	if (std::abs(candidate.probability - current.probability) <= epsilon)
	{
		if (candidate.svmScore > current.svmScore + epsilon)
		{
			return true;
		}
		if (std::abs(candidate.svmScore - current.svmScore) <= epsilon)
		{
			return candidate.order < current.order;
		}
	}
	return false;
}

std::vector<PsmRow> selectBestRowsByPeptideCharge(const std::vector<PsmRow> &rows)
{
	std::unordered_map<std::string, size_t> bestIndex;
	std::vector<PsmRow> selected;
	selected.reserve(rows.size());
	for (const PsmRow &row : rows)
	{
		const std::string key = peptideMassClassKey(row.peptide) + '\t' +
								std::to_string(row.precursorCharge);
		const auto found = bestIndex.find(key);
		if (found == bestIndex.end())
		{
			bestIndex[key] = selected.size();
			selected.push_back(row);
			continue;
		}

		PsmRow &current = selected[found->second];
		const std::string mergedProteins =
			mergeProteinLists(current.proteins, row.proteins);
		if (isBetterPsm(row, current))
		{
			current = row;
		}
		current.proteins = mergedProteins;
	}

	std::sort(selected.begin(), selected.end(),
			  [](const PsmRow &left, const PsmRow &right)
			  { return left.order < right.order; });
	return selected;
}

struct FragmentEntry
{
	double mz = 0.0;
	double intensity = 0.0;
	char ionKind = 'b';
	size_t position = 0;
};

struct SpectrumRecord
{
	std::string psmId;
	std::string retention;
	int charge = 1;
	std::string peptide;
	std::string proteins;
	double sipAbundancePct = 0.0;
	std::vector<double> precursorMz;
	std::vector<double> precursorIntensity;
	std::vector<double> fragmentMz;
	std::vector<double> theoreticalIntensity;
	std::vector<char> ionKinds;
	std::vector<uint64_t> ionPositions;
};

struct Hdf5OutputData
{
	std::vector<std::string> psmIds;
	std::vector<std::string> retentions;
	std::vector<int> charges;
	std::vector<std::string> peptides;
	std::vector<std::string> proteins;
	std::vector<double> sipAbundancePct;
	std::vector<double> precursorMz;
	std::vector<double> precursorIntensity;
	std::vector<uint64_t> precursorOffset;
	std::vector<uint64_t> precursorCount;
	std::vector<double> fragmentMz;
	std::vector<double> theoreticalIntensity;
	std::vector<char> ionKind;
	std::vector<uint64_t> ionPosition;
	std::vector<uint64_t> fragmentOffset;
	std::vector<uint64_t> fragmentCount;
};

void buildPrecursorChargePeaks(const IsotopeDistribution &dist,
							   int charge,
							   double probCutoff,
							   std::vector<double> &mzValues,
							   std::vector<double> &intensities)
{
	const double proton = ProNovoConfig::getProtonMass();
	std::vector<std::pair<double, double>> peaks;
	peaks.reserve(dist.vMass.size());
	for (size_t i = 0; i < dist.vMass.size() && i < dist.vProb.size(); ++i)
	{
		if (dist.vProb[i] < probCutoff)
		{
			continue;
		}
		const double mz = (dist.vMass[i] + static_cast<double>(charge) * proton) / static_cast<double>(charge);
		peaks.emplace_back(mz, dist.vProb[i]);
	}
	std::sort(peaks.begin(), peaks.end(),
			  [](const std::pair<double, double> &a, const std::pair<double, double> &b)
			  { return a.first < b.first; });
	double maximumIntensity = 0.0;
	for (const auto &peak : peaks)
	{
		maximumIntensity = std::max(maximumIntensity, peak.second);
	}
	mzValues.reserve(peaks.size());
	intensities.reserve(peaks.size());
	for (const auto &peak : peaks)
	{
		mzValues.push_back(peak.first);
		intensities.push_back(
			maximumIntensity > 0.0 ? peak.second / maximumIntensity : 0.0);
	}
}

void buildChargeOneFragments(const std::vector<std::vector<double>> &bMass,
							 const std::vector<std::vector<double>> &bProb,
							 const std::vector<std::vector<double>> &yMass,
							 const std::vector<std::vector<double>> &yProb,
							 double probCutoff,
							 SpectrumRecord &record)
{
	const double proton = ProNovoConfig::getProtonMass();
	std::vector<FragmentEntry> entries;
	entries.reserve(512);

	const auto collect = [&](const std::vector<std::vector<double>> &masses,
							 const std::vector<std::vector<double>> &probs,
							 const char ionKind)
	{
		for (size_t i = 0; i < masses.size() && i < probs.size(); ++i)
		{
			if (masses[i].empty() || probs[i].empty())
			{
				continue;
			}
			const double maximumProbability = *std::max_element(
				probs[i].begin(), probs[i].end());
			if (!(maximumProbability > 0.0))
			{
				continue;
			}
			for (size_t j = 0; j < masses[i].size() && j < probs[i].size(); ++j)
			{
				if (probs[i][j] < probCutoff)
				{
					continue;
				}
					FragmentEntry entry;
					entry.mz = masses[i][j] + proton; // charge 1 only
				entry.intensity = probs[i][j] / maximumProbability;
				entry.ionKind = ionKind;
				entry.position = i + 1;
				entries.push_back(entry);
				}
		}
	};

	collect(bMass, bProb, 'b');
	collect(yMass, yProb, 'y');

	std::sort(entries.begin(), entries.end(),
			  [](const FragmentEntry &a, const FragmentEntry &b)
			  { return a.mz < b.mz; });
	record.fragmentMz.reserve(entries.size());
	record.theoreticalIntensity.reserve(entries.size());
	record.ionKinds.reserve(entries.size());
	record.ionPositions.reserve(entries.size());
	for (const FragmentEntry &entry : entries)
	{
		record.fragmentMz.push_back(entry.mz);
		record.theoreticalIntensity.push_back(entry.intensity);
		record.ionKinds.push_back(entry.ionKind);
		record.ionPositions.push_back(static_cast<uint64_t>(entry.position));
	}
}

H5::DataSpace createDataspace(size_t count)
{
	const hsize_t dimension = static_cast<hsize_t>(count);
	return H5::DataSpace(1, &dimension);
}

H5::DataSpace createScalarDataspace()
{
	return H5::DataSpace(H5S_SCALAR);
}

H5::DSetCreatPropList createCompressedDatasetProperties(size_t count,
												 size_t chunkSize = 262144)
{
	H5::DSetCreatPropList properties;
	if (count > 0)
	{
		const hsize_t chunk = static_cast<hsize_t>(std::min(count, chunkSize));
		properties.setChunk(1, &chunk);
		properties.setShuffle();
		properties.setDeflate(6);
	}
	return properties;
}

template <typename T>
void writeVectorDataset(const H5::Group &group,
						const char *name,
						const H5::DataType &type,
						const std::vector<T> &values)
{
	H5::DataSpace space = createDataspace(values.size());
	H5::DataSet dataset = values.empty()
							 ? group.createDataSet(name, type, space)
							 : group.createDataSet(
								   name, type, space,
								   createCompressedDatasetProperties(values.size()));
	if (!values.empty())
	{
		dataset.write(values.data(), type);
	}
}

void writeZeroDoubleDataset(const H5::Group &group, const char *name, size_t count)
{
	H5::DataSpace space = createDataspace(count);
	if (count == 0)
	{
		group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
		return;
	}
	H5::DSetCreatPropList properties = createCompressedDatasetProperties(count);
	const double fillValue = 0.0;
	properties.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);
	group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, properties);
}

void writeStringDataset(const H5::Group &group,
						const char *name,
						const std::vector<std::string> &values)
{
	size_t width = 1;
	for (const std::string &value : values)
	{
		width = std::max(width, value.size() + 1);
	}
	H5::StrType type(H5::PredType::C_S1, width);
	type.setStrpad(H5T_STR_NULLTERM);
	type.setCset(H5T_CSET_UTF8);

	std::vector<char> flat(values.size() * width, '\0');
	for (size_t i = 0; i < values.size(); ++i)
	{
		std::memcpy(flat.data() + i * width,
					values[i].c_str(),
					std::min(values[i].size(), width - 1));
	}
	H5::DataSpace space = createDataspace(values.size());
	H5::DataSet dataset = values.empty()
							 ? group.createDataSet(name, type, space)
							 : group.createDataSet(
								   name, type, space,
								   createCompressedDatasetProperties(values.size()));
	if (!values.empty())
	{
		dataset.write(flat.data(), type);
	}
}

void writeStringAttribute(const H5::H5Object &object,
						  const char *name,
						  const std::string &value)
{
	H5::StrType type(H5::PredType::C_S1, std::max<size_t>(1, value.size() + 1));
	type.setStrpad(H5T_STR_NULLTERM);
	type.setCset(H5T_CSET_UTF8);
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attribute = object.createAttribute(name, type, space);
	attribute.write(type, value.c_str());
}

template <typename T>
void writeScalarAttribute(const H5::H5Object &object,
						  const char *name,
						  const H5::DataType &type,
						  const T &value)
{
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attribute = object.createAttribute(name, type, space);
	attribute.write(type, &value);
}

Hdf5OutputData flattenRecords(const std::vector<SpectrumRecord> &records,
							  const std::vector<char> &ok)
{
	Hdf5OutputData output;
	size_t recordCount = 0;
	size_t precursorValues = 0;
	size_t fragmentValues = 0;
	for (size_t i = 0; i < records.size() && i < ok.size(); ++i)
	{
		if (ok[i])
		{
			++recordCount;
			precursorValues += records[i].precursorMz.size();
			fragmentValues += records[i].fragmentMz.size();
		}
	}

	output.psmIds.reserve(recordCount);
	output.retentions.reserve(recordCount);
	output.charges.reserve(recordCount);
	output.peptides.reserve(recordCount);
	output.proteins.reserve(recordCount);
	output.sipAbundancePct.reserve(recordCount);
	output.precursorOffset.reserve(recordCount);
	output.precursorCount.reserve(recordCount);
	output.fragmentOffset.reserve(recordCount);
	output.fragmentCount.reserve(recordCount);
	output.precursorMz.reserve(precursorValues);
	output.precursorIntensity.reserve(precursorValues);
	output.fragmentMz.reserve(fragmentValues);
	output.theoreticalIntensity.reserve(fragmentValues);
	output.ionKind.reserve(fragmentValues);
	output.ionPosition.reserve(fragmentValues);

	for (size_t i = 0; i < records.size() && i < ok.size(); ++i)
	{
		if (!ok[i])
		{
			continue;
		}
		const SpectrumRecord &record = records[i];
		output.psmIds.push_back(record.psmId);
		output.retentions.push_back(record.retention);
		output.charges.push_back(record.charge);
		output.peptides.push_back(record.peptide);
		output.proteins.push_back(record.proteins);
		output.sipAbundancePct.push_back(record.sipAbundancePct);
		output.precursorOffset.push_back(static_cast<uint64_t>(output.precursorMz.size()));
		output.precursorCount.push_back(static_cast<uint64_t>(record.precursorMz.size()));
		output.fragmentOffset.push_back(static_cast<uint64_t>(output.fragmentMz.size()));
		output.fragmentCount.push_back(static_cast<uint64_t>(record.fragmentMz.size()));
		output.precursorMz.insert(output.precursorMz.end(),
							  record.precursorMz.begin(), record.precursorMz.end());
		output.precursorIntensity.insert(output.precursorIntensity.end(),
									 record.precursorIntensity.begin(),
									 record.precursorIntensity.end());
		output.fragmentMz.insert(output.fragmentMz.end(),
							 record.fragmentMz.begin(), record.fragmentMz.end());
		output.theoreticalIntensity.insert(output.theoreticalIntensity.end(),
									   record.theoreticalIntensity.begin(),
									   record.theoreticalIntensity.end());
		output.ionKind.insert(output.ionKind.end(),
						  record.ionKinds.begin(), record.ionKinds.end());
		output.ionPosition.insert(output.ionPosition.end(),
							  record.ionPositions.begin(), record.ionPositions.end());
	}
	return output;
}

fs::path resolvedOutputPath(const std::string &outputPath)
{
	fs::path path(outputPath);
	const std::string extension = sipros::TextUtils::toLower(path.extension().string());
	if (extension != ".sfi")
	{
		path += ".sfi";
	}
	return path;
}

bool writeSpectraSfi(const fs::path &path,
					 const Hdf5OutputData &data,
					 char sipAtom,
					 int sipIsotopeMassNumber,
					 double targetSipAbundancePct,
					 double probCutoff,
					 int threads,
					 sipros::SpectraIndexBuildStats &buildStats)
{
	std::vector<sipros::SpectraIndexRecordInput> records;
	records.reserve(data.psmIds.size());
	std::vector<double> retentionValues(data.psmIds.size(), 0.0);
	for (size_t i = 0; i < data.psmIds.size(); ++i)
	{
		retentionValues[i] = i < data.retentions.size()
			? parseFirstDouble(data.retentions[i]) : 0.0;
	}
	for (size_t i = 0; i < data.psmIds.size(); ++i)
	{
		sipros::SpectraIndexRecordInput record;
		record.psmId = data.psmIds[i];
		record.peptide = i < data.peptides.size() ? data.peptides[i] : std::string();
		record.proteins = i < data.proteins.size() ? data.proteins[i] : std::string();
		record.charge = i < data.charges.size() ? data.charges[i] : 1;
		// Accepted PSM reports define Retention in seconds. Store the index RT
		// in minutes so search does not need a unit heuristic.
		record.retentionMinutes = retentionValues[i] / 60.0;
		record.sipAbundancePct = i < data.sipAbundancePct.size()
			? data.sipAbundancePct[i]
			: (std::isfinite(targetSipAbundancePct) ? targetSipAbundancePct : 0.0);
		const size_t precursorOffset = static_cast<size_t>(data.precursorOffset[i]);
		const size_t precursorCount = static_cast<size_t>(data.precursorCount[i]);
		record.precursors.reserve(precursorCount);
		for (size_t k = 0; k < precursorCount; ++k)
		{
			record.precursors.push_back({data.precursorMz[precursorOffset + k],
				data.precursorIntensity[precursorOffset + k]});
		}
		const size_t fragmentOffset = static_cast<size_t>(data.fragmentOffset[i]);
		const size_t fragmentCount = static_cast<size_t>(data.fragmentCount[i]);
		record.fragments.reserve(fragmentCount);
		for (size_t k = 0; k < fragmentCount; ++k)
		{
			const size_t index = fragmentOffset + k;
			sipros::SpectraIndexFragmentPeakInput fragment;
			fragment.mz = data.fragmentMz[index];
			fragment.theoreticalIntensity = static_cast<float>(data.theoreticalIntensity[index]);
			fragment.experimentalIntensity = 0.0F;
			fragment.ionKind = static_cast<uint8_t>(data.ionKind[index]);
			fragment.ionPosition = static_cast<uint32_t>(data.ionPosition[index]);
			record.fragments.push_back(fragment);
		}
		records.push_back(std::move(record));
	}
	sipros::SpectraIndexMetadata metadata;
	metadata.chemistryProfileId = ProNovoConfig::getChemistryProfileId();
	metadata.recordKind = "theoretical";
	metadata.targetSipAbundancePct = std::isfinite(targetSipAbundancePct)
		? targetSipAbundancePct : 0.0;
	metadata.sipAtom = sipAtom;
	metadata.sipIsotopeMassNumber = sipIsotopeMassNumber;
	metadata.probabilityCutoff = probCutoff;
	metadata.envelopeTopN = 3;
	metadata.label = 1;
	std::string error;
	if (!sipros::SpectraIndex::write(
			path.string(), metadata, records, error, threads, &buildStats))
	{
		std::cerr << error << '\n';
		return false;
	}
	return true;
}

bool writeSpectraHdf5(const fs::path &path,
					  const Hdf5OutputData &data,
					  char sipAtom,
					  int sipIsotopeMassNumber,
					  double targetSipAbundancePct,
					  const std::string &abundanceSource,
					  double probCutoff)
{
	try
	{
		if (!path.parent_path().empty())
		{
			fs::create_directories(path.parent_path());
		}
		H5::H5File file(path.string(), H5F_ACC_TRUNC);
		const int formatVersion = kSpectraHdf5FormatVersion;
		writeScalarAttribute(file, "format_version", H5::PredType::NATIVE_INT,
							 formatVersion);
		writeStringAttribute(file, "chemistry_profile_id",
						 ProNovoConfig::getChemistryProfileId());
		writeStringAttribute(file, "record_kind", "theoretical");
		writeScalarAttribute(file, "target_sip_abundance_pct",
							 H5::PredType::NATIVE_DOUBLE, targetSipAbundancePct);
		writeStringAttribute(file, "sip_abundance_source", abundanceSource);
		writeStringAttribute(file, "sip_atom", std::string(1, sipAtom));
		writeScalarAttribute(file, "sip_isotope_mass_number",
							 H5::PredType::NATIVE_INT, sipIsotopeMassNumber);
		writeScalarAttribute(file, "prob_cutoff", H5::PredType::NATIVE_DOUBLE,
							 probCutoff);
		const double ppmTolerance = 0.0;
		const uint64_t minimumMatchedEnvelopes = 0;
		writeScalarAttribute(file, "ppm_tolerance", H5::PredType::NATIVE_DOUBLE,
							 ppmTolerance);
		writeScalarAttribute(file, "min_matched_envelopes",
							 H5::PredType::NATIVE_UINT64, minimumMatchedEnvelopes);

		H5::Group records = file.createGroup("records");
		H5::Group precursor = file.createGroup("precursor");
		H5::Group fragments = file.createGroup("fragments");
		writeStringDataset(records, "psm_id", data.psmIds);
		writeStringDataset(records, "retention", data.retentions);
		writeVectorDataset(records, "charge", H5::PredType::NATIVE_INT, data.charges);
		writeStringDataset(records, "peptide", data.peptides);
		writeStringDataset(records, "proteins", data.proteins);
		writeVectorDataset(records, "sip_abundance_pct", H5::PredType::NATIVE_DOUBLE,
						   data.sipAbundancePct);

		writeVectorDataset(precursor, "mz", H5::PredType::NATIVE_DOUBLE,
						   data.precursorMz);
		writeVectorDataset(precursor, "intensity", H5::PredType::NATIVE_DOUBLE,
						   data.precursorIntensity);
		writeVectorDataset(precursor, "offset", H5::PredType::NATIVE_UINT64,
						   data.precursorOffset);
		writeVectorDataset(precursor, "count", H5::PredType::NATIVE_UINT64,
						   data.precursorCount);

		writeVectorDataset(fragments, "mz", H5::PredType::NATIVE_DOUBLE,
						   data.fragmentMz);
		writeVectorDataset(fragments, "theoretical_intensity",
						   H5::PredType::NATIVE_DOUBLE, data.theoreticalIntensity);
		writeZeroDoubleDataset(fragments, "experimental_intensity",
						   data.fragmentMz.size());
		writeVectorDataset(fragments, "ion_kind", H5::PredType::NATIVE_CHAR,
						   data.ionKind);
		writeVectorDataset(fragments, "ion_position", H5::PredType::NATIVE_UINT64,
						   data.ionPosition);
		writeVectorDataset(fragments, "offset", H5::PredType::NATIVE_UINT64,
						   data.fragmentOffset);
		writeVectorDataset(fragments, "count", H5::PredType::NATIVE_UINT64,
						   data.fragmentCount);
		file.close();
		return true;
	}
	catch (const H5::Exception &exception)
	{
		std::cerr << exception.getDetailMsg() << '\n';
		return false;
	}
	catch (const std::exception &exception)
	{
		std::cerr << exception.what() << '\n';
		return false;
	}
}

bool buildPrecursorDistributionFromProductIons(
	const Isotopologue &iso,
	const std::vector<std::vector<double>> &yMass,
	const std::vector<std::vector<double>> &yProb,
	const std::vector<std::vector<double>> &bMass,
	const std::vector<std::vector<double>> &bProb,
	IsotopeDistribution &precursorDist)
{
	if (bMass.empty() || bProb.empty() || yMass.empty() || yProb.empty())
		return false;
	precursorDist = iso.sum(
		IsotopeDistribution(bMass.back(), bProb.back()),
		IsotopeDistribution(yMass.front(), yProb.front()));
	return !precursorDist.vMass.empty() &&
		precursorDist.vMass.size() == precursorDist.vProb.size();
}
} // namespace

int TheoreticalSpectraWorkflow::run(int argc, char **argv)
{
	for (int i = 1; i < argc; ++i)
	{
		const std::string opt = argv[i];
		if (opt == "-h" || opt == "--help")
		{
			printUsage(argv[0]);
			return 0;
		}
	}
	Args args;
	if (!parseArgs(argc, argv, args))
	{
		return 1;
	}

	if (!ProNovoConfig::load(ProNovoConfig::Profile::Sip))
	{
		std::cerr << "Could not initialize the built-in SIP profile.\n";
		return 1;
	}
	std::string fixedPtmError;
	if (!ProNovoConfig::configureFixedPtms(
			args.fixedPtmSelectors, fixedPtmError))
	{
		std::cerr << "Invalid fixed PTM selection: " << fixedPtmError << "\n";
		return 1;
	}
	char sipAtom = args.sipAtom;
	int sipIsotopeMassNumber = args.sipIsotopeMassNumber;
	std::string chemistryError;
	if (!ProNovoConfig::validatePreparationChemistry(
			ProNovoConfig::configIsotopologue, chemistryError))
	{
		std::cerr << chemistryError << "\n";
		return 1;
	}

	if (ProNovoConfig::atomIndex(sipAtom) < 0)
	{
		std::cerr << "Invalid SIP atom '" << sipAtom << "'. Valid targets: C13,H2,N15,O18,S34\n";
		return 1;
	}

	int sipIsotopeIndex = 1;
	try
	{
		sipIsotopeIndex = ProNovoConfig::resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, sipAtom, sipIsotopeMassNumber);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	const double fixedSipAbundancePct = args.fixedSipAbundancePct;

	std::cout << "SIP atom: " << sipAtom;
	if (sipIsotopeMassNumber > 0)
	{
		std::cout << sipIsotopeMassNumber;
	}
	std::cout << "\n";
	if (args.sipAbundanceMode == SipAbundanceMode::InputRow)
	{
		std::cout << "SIP abundance source: input MS2IsotopicAbundances\n";
	}
	else
	{
		std::cout << "SIP abundance source: user-defined (" << fixedSipAbundancePct << "%)\n";
	}

	if (args.threads > 0)
	{
		omp_set_num_threads(args.threads);
	}

	std::vector<PsmRow> rows;
	size_t inputFileCount = 0;
	size_t skippedUnsupportedMods = 0;
	try
	{
		rows = readInputRows(args.inputPath,
						 args.sipAbundanceMode == SipAbundanceMode::InputRow,
						 inputFileCount,
						 skippedUnsupportedMods);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}
	const size_t inputRowCount = rows.size();
	rows = selectBestRowsByPeptideCharge(rows);

	std::cout << "Input files: " << inputFileCount
			  << "; PSM rows read: " << inputRowCount
			  << "; retained peptide/charge rows: " << rows.size();
	if (skippedUnsupportedMods > 0)
	{
		std::cout << "; unsupported modified peptides skipped: "
				  << skippedUnsupportedMods;
	}
	std::cout << '\n';

	Isotopologue baseIso = ProNovoConfig::configIsotopologue;
	if (args.sipAbundanceMode == SipAbundanceMode::FixedUser)
	{
		try
		{
			ProNovoConfig::setSipAbundance(baseIso, sipAtom, sipIsotopeIndex, fixedSipAbundancePct);
		}
		catch (const std::exception &ex)
		{
			std::cerr << "Failed to apply user-defined SIP abundance: " << ex.what() << "\n";
			return 1;
		}
	}
	std::vector<SpectrumRecord> records(rows.size());
	std::vector<char> ok(rows.size(), 0);

	const size_t targetAtomIndex = static_cast<size_t>(
		ProNovoConfig::atomIndex(sipAtom));
#pragma omp parallel
	{
		Isotopologue localIso = baseIso;
#pragma omp for schedule(dynamic)
		for (int i = 0; i < static_cast<int>(rows.size()); ++i)
		{
			const auto &row = rows[static_cast<size_t>(i)];
			try
			{
				if (args.sipAbundanceMode == SipAbundanceMode::InputRow)
				{
					// Restore the pristine target categories before every row so
					// endpoint O18/S34 rows cannot affect later enrichments.
					localIso.vAtomIsotopicDistribution[
						targetAtomIndex].vProb =
						baseIso.vAtomIsotopicDistribution[
							targetAtomIndex].vProb;
					ProNovoConfig::setSipAbundance(localIso, sipAtom, sipIsotopeIndex, row.sipPct);
				}
				std::vector<std::vector<double>> yMass, yProb, bMass, bProb;
				if (!localIso.computeProductIon(row.peptide, yMass, yProb, bMass, bProb))
				{
#pragma omp critical
					{
						std::cerr << "Skipping PSM " << row.psmId << ": computeProductIon failed.\n";
					}
					continue;
				}
				IsotopeDistribution precursorDist;
				if (!buildPrecursorDistributionFromProductIons(
						localIso, yMass, yProb, bMass, bProb,
						precursorDist))
				{
#pragma omp critical
					{
						std::cerr << "Skipping PSM " << row.psmId
								  << ": precursor reconstruction from product ions failed.\n";
					}
					continue;
				}

				SpectrumRecord record;
				record.psmId = row.psmId;
				record.retention = row.retention;
				record.charge = row.precursorCharge;
				record.peptide = row.peptide;
				record.proteins = row.proteins;
				record.sipAbundancePct =
					args.sipAbundanceMode == SipAbundanceMode::InputRow
						? row.sipPct
						: fixedSipAbundancePct;
				buildPrecursorChargePeaks(precursorDist,
								  row.precursorCharge,
								  args.probCutoff,
								  record.precursorMz,
								  record.precursorIntensity);
				buildChargeOneFragments(
					bMass, bProb, yMass, yProb, args.probCutoff, record);
				records[static_cast<size_t>(i)] = std::move(record);
				ok[static_cast<size_t>(i)] = 1;
			}
			catch (const std::exception &ex)
			{
#pragma omp critical
				{
					std::cerr << "Skipping PSM " << row.psmId << ": " << ex.what() << "\n";
				}
			}
		}
	}

	const Hdf5OutputData output = flattenRecords(records, ok);
	const fs::path outputPath = resolvedOutputPath(args.outputPath);
	const double fileAbundance =
		args.sipAbundanceMode == SipAbundanceMode::InputRow
			? std::numeric_limits<double>::quiet_NaN()
			: fixedSipAbundancePct;
	sipros::SpectraIndexBuildStats buildStats;
	if (!writeSpectraSfi(outputPath,
						  output,
						  sipAtom,
						  sipIsotopeMassNumber,
						  fileAbundance,
						  args.probCutoff,
						  args.threads > 0 ? args.threads : omp_get_max_threads(),
						  buildStats))
	{
		return 1;
	}

	std::cerr << "Wrote theoretical spectra for " << output.psmIds.size()
			  << " PSMs to " << outputPath << "\n"
			  << "SFI v5 compact RT-aware index (top 3 peaks/envelope): fragments=" << buildStats.fragmentCount
			  << " x 16 bytes, packed_product_postings="
			  << buildStats.productPostingCount << " x 4 bytes, sparse_rt_bins="
			  << buildStats.rtBinCount << ", blocks="
			  << buildStats.blockCount << "\n"
			  << "  parallel product-index build: " << std::fixed
			  << std::setprecision(3) << buildStats.productIndexSeconds << "s ("
			  << buildStats.threadsUsed << " threads); total="
			  << buildStats.totalSeconds << "s; file="
			  << static_cast<double>(buildStats.fileBytes) /
				 (1024.0 * 1024.0 * 1024.0) << " GiB\n";
	return 0;
}
