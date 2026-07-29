// sipros search-spectra
//
// Re-score Raxport HDF5 MS2 scans against a pre-generated SIP spectra
// library (HDF5). HDF5 records are loaded into bounded in-memory batches.
// Output: one Percolator PIN per HDF5 sample file, 30 columns.

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cerrno>
#include <ctime>
#include <filesystem>
#include <chrono>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <numeric>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/resource.h>

#ifndef MFD_CLOEXEC
#define MFD_CLOEXEC 0x0001U
#endif

#include <omp.h>
#include "H5Cpp.h"

#include "proNovoConfig.h"
#include "ms2scanvector.h"
#include "ms2scan.h"
#include "peptide.h"
#include "MVH.h"
#include "CometSearchMod.h"
#include "SiprosWorkflows.h"
#include "PeptideIsotopeCalculator.h"
#include "RaxportHdf5Reader.h"
#include "PinWriter.h"
#include "PSMfeatureExtractor.h"

namespace fs = std::filesystem;

namespace
{
constexpr int kSpectraHdf5FormatVersion = 2;

// -------------------- Args --------------------

struct Args
{
	std::string workingDir;
	std::string singleHdf5;
	std::string hdf5Dir;
	std::string outputDir;
	int threads = 0;
	double rtToleranceMin = 5.0;
	// Separate MS1 (precursor) and MS2 (fragment) tolerances. Default 10 ppm
	// each. --tolerance / --tolerance-unit are shortcuts that set both.
	double toleranceMs1 = 10.0;
	bool toleranceMs1Ppm = true;
	double toleranceMs2 = 10.0;
	bool toleranceMs2Ppm = true;
	int isotopeShiftWindow = 3;
	int scoreEnvelopeTopN = 2;
	int topPsmsPerScan = 20;
};

// One tolerance, two units. Resolve to a Da window for a specific m/z.
struct MassTolerance
{
	bool ppm = true;
	double value = 10.0;
	double daAt(double mz) const
	{
		return ppm ? (value * mz / 1.0e6) : value;
	}
};

double timevalSeconds(const timeval &tv)
{
	return static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec) / 1.0e6;
}

double processTreeCpuSeconds()
{
	rusage self{};
	rusage children{};
	getrusage(RUSAGE_SELF, &self);
	getrusage(RUSAGE_CHILDREN, &children);
	return timevalSeconds(self.ru_utime) + timevalSeconds(self.ru_stime) +
		   timevalSeconds(children.ru_utime) + timevalSeconds(children.ru_stime);
}

struct TimingEntry
{
	std::string group;
	std::string label;
	std::string workName;
	double wallSeconds = 0.0;
	double cpuSeconds = 0.0;
	long memoryMbDelta = 0;
	size_t workItems = 0;
};

struct TimingLogger
{
	std::vector<TimingEntry> entries;

	void printHeader() const
	{
		std::cout << "\nTiming log\n"
				  << "  " << std::left << std::setw(34) << "Region"
				  << std::right << std::setw(10) << "Wall(s)"
				  << std::setw(10) << "CPU(s)"
				  << std::setw(9) << "CPU/Wall"
				  << std::setw(11) << "Mem(MB)"
				  << "  Work\n";
	}

	template <typename Fn>
	void run(const std::string &group,
			 const std::string &label,
			 size_t workItems,
			 const std::string &workName,
			 Fn &&fn)
	{
		const long memStart = static_cast<long>(checkMemoryUsage());
		const double cpuStart = processTreeCpuSeconds();
		const double wallStart = omp_get_wtime();
		fn();
		const double wallSeconds = omp_get_wtime() - wallStart;
		const double cpuSeconds = processTreeCpuSeconds() - cpuStart;
		const long memEnd = static_cast<long>(checkMemoryUsage());

		TimingEntry entry;
		entry.group = group;
		entry.label = label;
		entry.workName = workName;
		entry.wallSeconds = wallSeconds;
		entry.cpuSeconds = cpuSeconds;
		entry.memoryMbDelta = memEnd - memStart;
		entry.workItems = workItems;
		entries.push_back(entry);
		printEntry(entry);
	}

	static void printEntry(const TimingEntry &entry)
	{
		const double speedup = entry.wallSeconds > 0.0 ? entry.cpuSeconds / entry.wallSeconds : 0.0;
		std::cout << "  " << std::left << std::setw(34) << entry.label
				  << std::right << std::setw(10) << std::fixed << std::setprecision(3) << entry.wallSeconds
				  << std::setw(10) << std::fixed << std::setprecision(3) << entry.cpuSeconds
				  << std::setw(9) << std::fixed << std::setprecision(2) << speedup
				  << std::setw(11) << std::showpos << entry.memoryMbDelta << std::noshowpos;
		if (entry.workItems > 0 && !entry.workName.empty())
		{
			const double rate = entry.wallSeconds > 0.0
									? static_cast<double>(entry.workItems) / entry.wallSeconds
									: 0.0;
			std::cout << "  " << entry.workItems << ' ' << entry.workName
					  << " (" << std::fixed << std::setprecision(0) << rate << "/s)";
		}
		std::cout << '\n';
	}

	void printSummary(double totalWallSeconds,
					  double totalCpuSeconds,
					  size_t totalRecords,
					  size_t totalCandidates,
					  size_t retainedPsms,
					  size_t writtenRows) const
	{
		std::map<std::string, TimingEntry> totals;
		std::vector<std::string> groupOrder;
		for (const TimingEntry &entry : entries)
		{
			if (totals.find(entry.group) == totals.end())
				groupOrder.push_back(entry.group);
			TimingEntry &sum = totals[entry.group];
			sum.group = entry.group;
			sum.label = entry.group;
			sum.wallSeconds += entry.wallSeconds;
			sum.cpuSeconds += entry.cpuSeconds;
			sum.memoryMbDelta += entry.memoryMbDelta;
			sum.workItems += entry.workItems;
			if (sum.workName.empty())
				sum.workName = entry.workName;
		}

		std::cout << "\nTiming summary\n"
				  << "  " << std::left << std::setw(28) << "Stage"
				  << std::right << std::setw(10) << "Wall(s)"
				  << std::setw(9) << "Share"
				  << std::setw(10) << "CPU(s)"
				  << std::setw(9) << "CPU/Wall"
				  << std::setw(11) << "Mem(MB)"
				  << '\n';
		for (const std::string &group : groupOrder)
		{
			const TimingEntry &entry = totals[group];
			const double share = totalWallSeconds > 0.0 ? entry.wallSeconds / totalWallSeconds * 100.0 : 0.0;
			const double speedup = entry.wallSeconds > 0.0 ? entry.cpuSeconds / entry.wallSeconds : 0.0;
			std::cout << "  " << std::left << std::setw(28) << entry.label
					  << std::right << std::setw(10) << std::fixed << std::setprecision(3) << entry.wallSeconds
					  << std::setw(8) << std::fixed << std::setprecision(1) << share << "%"
					  << std::setw(10) << std::fixed << std::setprecision(3) << entry.cpuSeconds
					  << std::setw(9) << std::fixed << std::setprecision(2) << speedup
					  << std::setw(11) << std::showpos << entry.memoryMbDelta << std::noshowpos
					  << '\n';
		}

		const double totalSpeedup = totalWallSeconds > 0.0 ? totalCpuSeconds / totalWallSeconds : 0.0;
		std::cout << "  " << std::left << std::setw(28) << "End-to-end"
				  << std::right << std::setw(10) << std::fixed << std::setprecision(3) << totalWallSeconds
				  << std::setw(8) << "100.0%"
				  << std::setw(10) << std::fixed << std::setprecision(3) << totalCpuSeconds
				  << std::setw(9) << std::fixed << std::setprecision(2) << totalSpeedup
				  << std::setw(11) << "-" << '\n';

		std::cout << "\nThroughput\n"
				  << "  HDF5 records loaded:       " << totalRecords << '\n'
				  << "  Assigned candidates:       " << totalCandidates << '\n'
				  << "  Retained top PSMs:         " << retainedPsms << '\n'
				  << "  PIN rows written:          " << writtenRows << '\n';
		if (totalWallSeconds > 0.0)
		{
			std::cout << "  Records/sec end-to-end:    "
					  << std::fixed << std::setprecision(0)
					  << static_cast<double>(totalRecords) / totalWallSeconds << '\n'
					  << "  Candidates/sec end-to-end: "
					  << std::fixed << std::setprecision(0)
					  << static_cast<double>(totalCandidates) / totalWallSeconds << '\n';
		}
		std::cout << '\n';
	}
};

void printUsage(const char *prog)
{
	std::cerr
		<< "Usage:\n  " << prog
		<< " -w <Raxport HDF5 dir> [-f <single.h5>]\n"
		<< "    -h5 <SIP spectra dir> -o <PIN output dir> [-t <N>] [--rt-tolerance <min>]\n"
		<< "    [--tolerance-ms1 <N>] [--tolerance-ms1-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance-ms2 <N>] [--tolerance-ms2-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance <N>] [--tolerance-unit ppm|da]            shortcut: set BOTH MS1 and MS2\n"
		<< "    [--isotope-shift-window <N>]                           precursor isotope shifts +/-N (default: 3)\n"
		<< "    [--score-envelope-top-n <N>]                           Xcorr/MVH envelope peaks (default: 2)\n"
		<< "    [--top-psms-per-scan <N>]                              WDP winners per scan/label (default: 20)\n";
}

Args parseArgs(int argc, char **argv)
{
	Args a;
	std::vector<std::string> v;
	for (int i = 1; i < argc; ++i)
		v.emplace_back(argv[i]);
	for (size_t i = 0; i < v.size(); ++i)
	{
		const std::string &k = v[i];
		auto next = [&]() -> std::string
		{
			if (i + 1 >= v.size())
			{
				std::cerr << "Missing value after " << k << "\n";
				std::exit(1);
			}
			return v[++i];
		};
		if (k == "-w")
			a.workingDir = next();
		else if (k == "-f")
			a.singleHdf5 = next();
		else if (k == "-h5")
			a.hdf5Dir = next();
		else if (k == "-o")
			a.outputDir = next();
		else if (k == "-t")
			a.threads = std::atoi(next().c_str());
		else if (k == "--rt-tolerance")
			a.rtToleranceMin = std::atof(next().c_str());
		else if (k == "--tolerance")
		{
			double v = std::atof(next().c_str());
			a.toleranceMs1 = v;
			a.toleranceMs2 = v;
		}
		else if (k == "--tolerance-unit")
		{
			std::string u = next();
			for (char &c : u)
				c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
			bool ppm;
			if (u == "ppm")
				ppm = true;
			else if (u == "da" || u == "dalton")
				ppm = false;
			else
			{
				std::cerr << "--tolerance-unit must be 'ppm' or 'da' (got '" << u << "')\n";
				std::exit(1);
			}
			a.toleranceMs1Ppm = ppm;
			a.toleranceMs2Ppm = ppm;
		}
		else if (k == "--tolerance-ms1")
			a.toleranceMs1 = std::atof(next().c_str());
		else if (k == "--tolerance-ms2")
			a.toleranceMs2 = std::atof(next().c_str());
		else if (k == "--tolerance-ms1-unit")
		{
			std::string u = next();
			for (char &c : u)
				c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
			if (u == "ppm")
				a.toleranceMs1Ppm = true;
			else if (u == "da" || u == "dalton")
				a.toleranceMs1Ppm = false;
			else
			{
				std::cerr << "--tolerance-ms1-unit must be 'ppm' or 'da' (got '" << u << "')\n";
				std::exit(1);
			}
		}
		else if (k == "--tolerance-ms2-unit")
		{
			std::string u = next();
			for (char &c : u)
				c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
			if (u == "ppm")
				a.toleranceMs2Ppm = true;
			else if (u == "da" || u == "dalton")
				a.toleranceMs2Ppm = false;
			else
			{
				std::cerr << "--tolerance-ms2-unit must be 'ppm' or 'da' (got '" << u << "')\n";
				std::exit(1);
			}
		}
		else if (k == "--isotope-shift-window")
			a.isotopeShiftWindow = std::atoi(next().c_str());
		else if (k == "--score-envelope-top-n")
			a.scoreEnvelopeTopN = std::atoi(next().c_str());
		else if (k == "--top-psms-per-scan")
			a.topPsmsPerScan = std::atoi(next().c_str());
		else if (k == "-h" || k == "--help")
		{
			printUsage(argv[0]);
			std::exit(0);
		}
		else
		{
			std::cerr << "Unknown option: " << k << "\n";
			printUsage(argv[0]);
			std::exit(1);
		}
	}
	if (a.singleHdf5.empty() && a.workingDir.empty())
		a.workingDir = ".";
	if (a.outputDir.empty())
		a.outputDir = a.workingDir;
	if (a.hdf5Dir.empty())
	{
		std::cerr << "-h5 <SIP spectra dir> is required\n";
		std::exit(1);
	}
	if (a.threads <= 0)
	{
		long n = sysconf(_SC_NPROCESSORS_ONLN);
		a.threads = (n > 0 ? static_cast<int>(n) : 1);
	}
	if (a.scoreEnvelopeTopN <= 0)
	{
		std::cerr << "--score-envelope-top-n must be > 0\n";
		std::exit(1);
	}
	if (a.topPsmsPerScan <= 0)
	{
		std::cerr << "--top-psms-per-scan must be > 0\n";
		std::exit(1);
	}
	if (a.isotopeShiftWindow < 0)
	{
		std::cerr << "--isotope-shift-window must be >= 0\n";
		std::exit(1);
	}
	return a;
}

std::vector<std::string> listFiles(const std::string &dir,
								   const std::vector<std::string> &exts)
{
	std::vector<std::string> out;
	if (!fs::is_directory(dir))
		return out;
	for (const auto &e : fs::directory_iterator(dir))
	{
		if (!e.is_regular_file())
			continue;
		const std::string ext = e.path().extension().string();
		for (const auto &p : exts)
		{
			if (ext == p)
			{
				out.push_back(e.path().string());
				break;
			}
		}
	}
	std::sort(out.begin(), out.end());
	return out;
}

bool isDecoyHdf5(const std::string &path)
{
	std::string stem = fs::path(path).filename().string();
	std::string lower;
	lower.reserve(stem.size());
	for (char c : stem)
		lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
	// match either trailing _decoy.h5 or any filename containing "decoy"
	return lower.find("decoy") != std::string::npos;
}

// -------------------- SipRecord --------------------

struct SipRecord
{
	std::string psmId;
	std::string retention;
	std::string peptide;     // bracketed/decorated form from HDF5
	std::string nakedPeptide; // letters-only used for scoring
	std::string proteins;
	int charge = 1;
	double rtMinutes = 0.0;
	double topPrecursorMz = 0.0;
	double topPrecursorIntensity = 0.0;
	double sumPrecursorIntensity = 0.0;
	std::vector<double> precursorMz;
	std::vector<double> precursorIntensity;
	// indexed by ion position (1-based ⇒ pos-1 row); each row holds an envelope
	std::vector<std::vector<double>> vvdYionMass;
	std::vector<std::vector<double>> vvdYionProb;
	std::vector<std::vector<double>> vvdBionMass;
	std::vector<std::vector<double>> vvdBionProb;
	// flattened fragment view for entropy/cosine
	std::vector<double> fragMz;
	std::vector<double> fragExpInt;
	size_t nYpositions = 0;
	size_t nBpositions = 0;
};

struct Hdf5FileMeta
{
	int formatVersion = 0;
	std::string chemistryProfileId;
	double targetSipAbundancePct = 0.0;
	std::string sipAtom;
	int sipIsotopeMassNumber = -1;
};

std::string readRequiredStringAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path);
int readRequiredIntAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path);
double readRequiredDoubleAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path);
void readSpectraFormatMetadata(
	const H5::H5Object &obj,
	const std::string &path,
	int &formatVersion,
	std::string &chemistryProfileId);
void readAndValidateSpectraFormatMetadata(
	const H5::H5Object &obj,
	const std::string &path,
	int &formatVersion,
	std::string &chemistryProfileId);

bool validateSpectraLibraryMetadata(
	const std::vector<std::string> &paths,
	std::string &canonicalSipIsotope,
	std::string &error)
{
	canonicalSipIsotope.clear();
	if (paths.empty())
	{
		error = "No spectra-library metadata was provided.";
		return false;
	}

	std::string canonicalChemistryProfileId;
	try
	{
		H5::H5File firstFile(paths.front(), H5F_ACC_RDONLY);
		int firstFormatVersion = 0;
		readSpectraFormatMetadata(
			firstFile, paths.front(), firstFormatVersion,
			canonicalChemistryProfileId);
	}
	catch (const H5::Exception &e)
	{
		error = "Cannot read spectra library metadata from " + paths.front() +
			": " + (e.getCDetailMsg() ? e.getCDetailMsg() : "HDF5 error");
		return false;
	}
	catch (const std::exception &e)
	{
		error = e.what();
		return false;
	}

	std::string chemistryError;
	if (!ProNovoConfig::configureChemistryProfileId(
			canonicalChemistryProfileId, chemistryError))
	{
		error = "Cannot use spectra library " + paths.front() + ": " +
			chemistryError;
		return false;
	}

	std::string firstIsotopePath;
	for (const std::string &path : paths)
	{
		try
		{
			H5::H5File file(path, H5F_ACC_RDONLY);
			int formatVersion = 0;
			std::string chemistryProfileId;
			readSpectraFormatMetadata(
				file, path, formatVersion, chemistryProfileId);
			if (chemistryProfileId != canonicalChemistryProfileId)
			{
				error = "Mixed chemistry_profile_id values in spectra libraries: " +
					paths.front() + " uses '" + canonicalChemistryProfileId +
					"', while " + path + " uses '" + chemistryProfileId +
					"'. Search one preparation chemistry at a time.";
				return false;
			}
			if (chemistryProfileId != ProNovoConfig::getChemistryProfileId())
			{
				error = "Spectra-library chemistry_profile_id '" +
					chemistryProfileId +
					"' was not reproduced by the compiled chemistry.";
				return false;
			}
			const double targetSipAbundancePct =
				readRequiredDoubleAttribute(
					file, "target_sip_abundance_pct", path);
			if (!std::isfinite(targetSipAbundancePct) ||
				targetSipAbundancePct < 0.0 || targetSipAbundancePct > 100.0)
			{
				error = "Invalid target_sip_abundance_pct in " + path +
					": " + std::to_string(targetSipAbundancePct) +
					"; expected a percentage in [0, 100].";
				return false;
			}
			const std::string sipAtom =
				readRequiredStringAttribute(file, "sip_atom", path);
			const int isotopeMassNumber =
				readRequiredIntAttribute(file, "sip_isotope_mass_number", path);
			const std::string currentSipIsotope =
				PSMfeatureExtractor::canonicalSipIsotope(
					sipAtom, isotopeMassNumber);
			if (currentSipIsotope.empty())
			{
				error = "Unsupported SIP isotope metadata in " + path +
					": atom=" + sipAtom + ", mass_number=" +
					std::to_string(isotopeMassNumber) +
					". Supported labels: C13,H2,N15,O18,S34.";
				return false;
			}
			if (canonicalSipIsotope.empty())
			{
				canonicalSipIsotope = currentSipIsotope;
				firstIsotopePath = path;
			}
			else if (currentSipIsotope != canonicalSipIsotope)
			{
				error = "Mixed SIP isotope targets in spectra libraries: " +
					firstIsotopePath + " uses " + canonicalSipIsotope +
					", while " + path + " uses " + currentSipIsotope +
					". Search one isotope target at a time.";
				return false;
			}
			if (!ProNovoConfig::validatePreparationChemistry(
					ProNovoConfig::configIsotopologue, error))
			{
				error = "Cannot use spectra library " + path + ": " + error;
				return false;
			}
		}
		catch (const H5::Exception &e)
		{
			error = "Cannot read spectra library metadata from " + path + ": " +
				(e.getCDetailMsg() ? e.getCDetailMsg() : "HDF5 error");
			return false;
		}
		catch (const std::exception &e)
		{
			error = e.what();
			return false;
		}
	}
	return true;
}

// Strip [..]-style decorations, leaving naked letters for scoring.
// SIP peptides use the form "[XXX]" for the N/C terminus and modifications.
std::string nakedPeptideOf(const std::string &decorated)
{
	std::string out;
	out.reserve(decorated.size());
	for (char c : decorated)
	{
		if (std::isalpha(static_cast<unsigned char>(c)))
			out.push_back(c);
	}
	return out;
}

std::string trimCopy(const std::string &s)
{
	size_t begin = 0;
	while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])))
		++begin;
	size_t end = s.size();
	while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1])))
		--end;
	return s.substr(begin, end - begin);
}

std::vector<std::string> splitProteinList(const std::string &proteins)
{
	std::string inner = trimCopy(proteins);
	if (inner.size() >= 2 && inner.front() == '{' && inner.back() == '}')
		inner = inner.substr(1, inner.size() - 2);

	std::vector<std::string> out;
	std::stringstream ss(inner);
	std::string token;
	while (std::getline(ss, token, ','))
	{
		token = trimCopy(token);
		if (!token.empty())
			out.push_back(token);
	}
	return out;
}

std::string formatProteinList(const std::vector<std::string> &proteins)
{
	std::string out = "{";
	for (const std::string &protein : proteins)
	{
		if (protein.empty())
			continue;
		if (out.size() > 1)
			out += ",";
		out += protein;
	}
	out += "}";
	return out;
}

std::string normalizeProteinList(const std::string &proteins)
{
	std::vector<std::string> parsed = splitProteinList(proteins);
	std::set<std::string> seen;
	std::vector<std::string> unique;
	unique.reserve(parsed.size());
	for (const std::string &protein : parsed)
	{
		if (seen.insert(protein).second)
			unique.push_back(protein);
	}
	return formatProteinList(unique);
}

std::string mergeProteinLists(const std::string &a, const std::string &b)
{
	std::set<std::string> seen;
	std::vector<std::string> merged;
	for (const std::string &protein : splitProteinList(a))
	{
		if (seen.insert(protein).second)
			merged.push_back(protein);
	}
	for (const std::string &protein : splitProteinList(b))
	{
		if (seen.insert(protein).second)
			merged.push_back(protein);
	}
	return formatProteinList(merged);
}

// Read an entire 1-D string dataset (fixed-length) into a vector<string>.
std::vector<std::string> readStringDataset(H5::Group &g, const char *name)
{
	H5::DataSet ds = g.openDataSet(name);
	H5::DataSpace sp = ds.getSpace();
	hsize_t dims[1] = {0};
	sp.getSimpleExtentDims(dims);
	size_t n = static_cast<size_t>(dims[0]);
	H5::StrType type = ds.getStrType();
	if (type.isVariableStr())
	{
		std::vector<char *> buf(n, nullptr);
		ds.read(buf.data(), type);
		std::vector<std::string> out;
		out.reserve(n);
		for (size_t i = 0; i < n; ++i)
		{
			out.emplace_back(buf[i] ? buf[i] : "");
			if (buf[i])
				std::free(buf[i]);
		}
		return out;
	}
	const size_t width = type.getSize();
	std::vector<char> flat(n * width, '\0');
	ds.read(flat.data(), type);
	std::vector<std::string> out;
	out.reserve(n);
	for (size_t i = 0; i < n; ++i)
	{
		const char *p = flat.data() + i * width;
		size_t len = ::strnlen(p, width);
		out.emplace_back(p, len);
	}
	return out;
}

template <typename T>
std::vector<T> readVecDataset(H5::Group &g, const char *name, const H5::PredType &type)
{
	H5::DataSet ds = g.openDataSet(name);
	H5::DataSpace sp = ds.getSpace();
	hsize_t dims[1] = {0};
	sp.getSimpleExtentDims(dims);
	std::vector<T> v(static_cast<size_t>(dims[0]));
	if (!v.empty())
		ds.read(v.data(), type);
	return v;
}

std::string readRequiredStringAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path)
{
	if (!obj.attrExists(name))
	{
		throw std::runtime_error(
			"Missing required HDF5 attribute '" + std::string(name) +
			"' in " + path);
	}
	H5::Attribute a = obj.openAttribute(name);
	H5::StrType type = a.getStrType();
	std::string value;
	if (type.isVariableStr())
	{
		char *p = nullptr;
		a.read(type, &p);
		value = p ? p : "";
		if (p)
			std::free(p);
	}
	else
	{
		std::vector<char> buf(type.getSize() + 1, '\0');
		a.read(type, buf.data());
		value = std::string(buf.data());
	}
	if (value.empty())
	{
		throw std::runtime_error(
			"Required HDF5 attribute '" + std::string(name) +
			"' is empty in " + path);
	}
	return value;
}

double readRequiredDoubleAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path)
{
	if (!obj.attrExists(name))
	{
		throw std::runtime_error(
			"Missing required HDF5 attribute '" + std::string(name) +
			"' in " + path);
	}
	double x = 0.0;
	H5::Attribute a = obj.openAttribute(name);
	a.read(H5::PredType::NATIVE_DOUBLE, &x);
	return x;
}

int readRequiredIntAttribute(
	const H5::H5Object &obj, const char *name, const std::string &path)
{
	if (!obj.attrExists(name))
	{
		throw std::runtime_error(
			"Missing required HDF5 attribute '" + std::string(name) +
			"' in " + path);
	}
	int x = 0;
	H5::Attribute a = obj.openAttribute(name);
	a.read(H5::PredType::NATIVE_INT, &x);
	return x;
}

void readSpectraFormatMetadata(
	const H5::H5Object &obj,
	const std::string &path,
	int &formatVersion,
	std::string &chemistryProfileId)
{
	formatVersion = readRequiredIntAttribute(obj, "format_version", path);
	if (formatVersion != kSpectraHdf5FormatVersion)
	{
		throw std::runtime_error(
			"Unsupported spectra HDF5 format_version in " + path + ": " +
			std::to_string(formatVersion) + "; expected " +
			std::to_string(kSpectraHdf5FormatVersion) +
			" (source-aware chemistry schema). Regenerate this spectra "
			"library with the current Sipros build.");
	}
	if (!obj.attrExists("chemistry_profile_id"))
	{
		throw std::runtime_error(
			"Missing required HDF5 attribute 'chemistry_profile_id' in " +
			path + ". This library does not identify its preparation "
			"chemistry; regenerate it with the current Sipros build.");
	}
	chemistryProfileId =
		readRequiredStringAttribute(obj, "chemistry_profile_id", path);
}

void readAndValidateSpectraFormatMetadata(
	const H5::H5Object &obj,
	const std::string &path,
	int &formatVersion,
	std::string &chemistryProfileId)
{
	readSpectraFormatMetadata(
		obj, path, formatVersion, chemistryProfileId);
	const std::string expected = ProNovoConfig::getChemistryProfileId();
	if (chemistryProfileId != expected)
	{
		throw std::runtime_error(
			"Incompatible chemistry_profile_id in " + path + ": '" +
			chemistryProfileId + "'; expected '" + expected +
			"'. Regenerate the spectra library with the same Sipros "
			"chemistry profile used for search.");
	}
}

// Best-effort parse of "1.23" minutes (or "PT74S"-like) into minutes.
double parseRetentionMinutes(const std::string &s)
{
	try
	{
		size_t pos = 0;
		double v = std::stod(s, &pos);
		// allow trailing units silently
		return v;
	}
	catch (...)
	{
		return 0.0;
	}
}

// Load all records from one HDF5 file.
bool loadHdf5File(const std::string &path,
				  std::vector<SipRecord> &out,
				  Hdf5FileMeta &meta,
				  std::string &errOut)
{
	try
	{
		H5::H5File f(path, H5F_ACC_RDONLY);
		readAndValidateSpectraFormatMetadata(
			f, path, meta.formatVersion, meta.chemistryProfileId);
		meta.targetSipAbundancePct = readRequiredDoubleAttribute(
			f, "target_sip_abundance_pct", path);
		meta.sipAtom = readRequiredStringAttribute(f, "sip_atom", path);
		meta.sipIsotopeMassNumber =
			readRequiredIntAttribute(f, "sip_isotope_mass_number", path);
		const std::string canonicalSipIsotope =
			PSMfeatureExtractor::canonicalSipIsotope(
				meta.sipAtom, meta.sipIsotopeMassNumber);
		if (canonicalSipIsotope.empty())
		{
			throw std::runtime_error(
				"Unsupported SIP isotope metadata in " + path +
				": atom=" + meta.sipAtom +
				", mass_number=" +
				std::to_string(meta.sipIsotopeMassNumber) +
				". Supported labels: C13,H2,N15,O18,S34.");
		}
		const std::string expectedSipIsotope =
			PSMfeatureExtractor::canonicalSipIsotope(
				ProNovoConfig::getSetSIPelement());
		if (canonicalSipIsotope != expectedSipIsotope)
		{
			throw std::runtime_error(
				"Incompatible SIP isotope metadata in " + path + ": " +
				canonicalSipIsotope + "; expected " + expectedSipIsotope +
				" from the validated spectra-library set.");
		}
		meta.sipAtom = canonicalSipIsotope;
		std::string chemistryError;
		if (!ProNovoConfig::validatePreparationChemistry(
				ProNovoConfig::configIsotopologue, chemistryError))
		{
			throw std::runtime_error(chemistryError);
		}

		H5::Group records = f.openGroup("records");
		H5::Group precursor = f.openGroup("precursor");
		H5::Group fragments = f.openGroup("fragments");

		auto psmIds = readStringDataset(records, "psm_id");
		auto retentions = readStringDataset(records, "retention");
		auto peptides = readStringDataset(records, "peptide");
		if (!records.nameExists("proteins"))
			throw std::runtime_error("Missing required HDF5 dataset records/proteins in " + path);
		auto proteins = readStringDataset(records, "proteins");
		auto charges = readVecDataset<int>(records, "charge", H5::PredType::NATIVE_INT);

		auto preMz = readVecDataset<double>(precursor, "mz", H5::PredType::NATIVE_DOUBLE);
		auto preInt = readVecDataset<double>(precursor, "intensity", H5::PredType::NATIVE_DOUBLE);
		auto preOff = readVecDataset<uint64_t>(precursor, "offset", H5::PredType::NATIVE_UINT64);
		auto preCnt = readVecDataset<uint64_t>(precursor, "count", H5::PredType::NATIVE_UINT64);

		auto fragMz = readVecDataset<double>(fragments, "mz", H5::PredType::NATIVE_DOUBLE);
		auto fragTh = readVecDataset<double>(fragments, "theoretical_intensity", H5::PredType::NATIVE_DOUBLE);
		auto fragEx = readVecDataset<double>(fragments, "experimental_intensity", H5::PredType::NATIVE_DOUBLE);
		auto fragOff = readVecDataset<uint64_t>(fragments, "offset", H5::PredType::NATIVE_UINT64);
		auto fragCnt = readVecDataset<uint64_t>(fragments, "count", H5::PredType::NATIVE_UINT64);
		auto fragPos = readVecDataset<uint64_t>(fragments, "ion_position", H5::PredType::NATIVE_UINT64);
		// ion_kind is a char dataset; H5 stores it as int8 or string
		std::vector<char> fragKind;
		{
			H5::DataSet ds = fragments.openDataSet("ion_kind");
			H5::DataSpace sp = ds.getSpace();
			hsize_t dims[1] = {0};
			sp.getSimpleExtentDims(dims);
			fragKind.assign(static_cast<size_t>(dims[0]), 0);
			H5::DataType dt = ds.getDataType();
			ds.read(fragKind.data(), dt);
		}

		size_t n = psmIds.size();
		// Auto-detect retention units: typical LC runs are < 200 minutes; if
		// the max parsed value is much larger, the file is using seconds.
		bool retentionIsSeconds = false;
		{
			double maxRt = 0.0;
			for (const auto &s : retentions)
			{
				double v = parseRetentionMinutes(s);
				if (v > maxRt)
					maxRt = v;
			}
			if (maxRt > 200.0)
				retentionIsSeconds = true;
		}
		out.clear();
		out.reserve(n);
		if (proteins.size() != n)
		{
			throw std::runtime_error("HDF5 dataset records/proteins has " +
									 std::to_string(proteins.size()) + " entries for " +
									 std::to_string(n) + " records in " + path);
		}
		for (size_t i = 0; i < n; ++i)
		{
			SipRecord r;
			r.psmId = psmIds[i];
			r.retention = i < retentions.size() ? retentions[i] : "";
			r.rtMinutes = parseRetentionMinutes(r.retention);
			if (retentionIsSeconds)
				r.rtMinutes /= 60.0;
			r.peptide = i < peptides.size() ? peptides[i] : "";
			r.nakedPeptide = nakedPeptideOf(r.peptide);
			if (splitProteinList(proteins[i]).empty())
			{
				throw std::runtime_error("Missing protein names for HDF5 record " +
										 r.psmId + " in " + path);
			}
			r.proteins = normalizeProteinList(proteins[i]);
			r.charge = i < charges.size() ? charges[i] : 1;

			// precursor envelope
			size_t pOff = i < preOff.size() ? static_cast<size_t>(preOff[i]) : 0;
			size_t pCnt = i < preCnt.size() ? static_cast<size_t>(preCnt[i]) : 0;
			r.precursorMz.assign(preMz.begin() + pOff, preMz.begin() + pOff + pCnt);
			r.precursorIntensity.assign(preInt.begin() + pOff, preInt.begin() + pOff + pCnt);
			r.topPrecursorMz = 0.0;
			r.topPrecursorIntensity = 0.0;
			r.sumPrecursorIntensity = 0.0;
			for (size_t k = 0; k < pCnt; ++k)
			{
				double intensity = r.precursorIntensity[k];
				r.sumPrecursorIntensity += intensity;
				if (intensity > r.topPrecursorIntensity)
				{
					r.topPrecursorIntensity = intensity;
					r.topPrecursorMz = r.precursorMz[k];
				}
			}

			// fragments: fragMz keeps the HDF5 charge-1 m/z for direct
			// spectrum-shape features. The old Sipros WDP/MVH/Xcorr helpers
			// expect ion masses before their final +proton m/z conversion, so
			// vvd?ionMass stores (HDF5 m/z - proton).
			size_t fOff = i < fragOff.size() ? static_cast<size_t>(fragOff[i]) : 0;
			size_t fCnt = i < fragCnt.size() ? static_cast<size_t>(fragCnt[i]) : 0;
			r.fragMz.assign(fragMz.begin() + fOff, fragMz.begin() + fOff + fCnt);
			r.fragExpInt.assign(fragEx.begin() + fOff, fragEx.begin() + fOff + fCnt);
			const double protonMass = ProNovoConfig::getProtonMass();

			// MVH/Xcorr scoring functions require vvd?ionMass.size() ==
			// peptideLength - 1. Size to that and leave inner vectors empty
			// where HDF5 has no fragments for that position.
			const size_t pepLen = r.nakedPeptide.size();
			const size_t need = pepLen > 0 ? pepLen - 1 : 0;
			r.nYpositions = need;
			r.nBpositions = need;
			r.vvdYionMass.assign(need, {});
			r.vvdYionProb.assign(need, {});
			r.vvdBionMass.assign(need, {});
			r.vvdBionProb.assign(need, {});

			// per-position sums for probability normalization
			std::vector<double> sumY(need, 0.0), sumB(need, 0.0);
			for (size_t k = 0; k < fCnt; ++k)
			{
				char kind = fragKind[fOff + k];
				size_t pos = static_cast<size_t>(fragPos[fOff + k]);
				double m = fragMz[fOff + k];
				double t = fragTh[fOff + k];
				if (pos == 0 || pos > need)
					continue;
				if (kind == 'y' || kind == 'Y')
				{
					r.vvdYionMass[pos - 1].push_back(m - protonMass);
					r.vvdYionProb[pos - 1].push_back(t);
					sumY[pos - 1] += t;
				}
				else if (kind == 'b' || kind == 'B')
				{
					r.vvdBionMass[pos - 1].push_back(m - protonMass);
					r.vvdBionProb[pos - 1].push_back(t);
					sumB[pos - 1] += t;
				}
			}
			for (size_t p = 0; p < need; ++p)
				if (sumY[p] > 0.0)
					for (double &x : r.vvdYionProb[p])
						x /= sumY[p];
			for (size_t p = 0; p < need; ++p)
				if (sumB[p] > 0.0)
					for (double &x : r.vvdBionProb[p])
						x /= sumB[p];

			// Sipros scoring functions index vvdYionMass[i] etc. and call
			// findProductIonSIP which dereferences vdIonMass[0] — empty
			// envelopes (e.g. for decoy records missing some positions)
			// would segfault. Fill any empty envelope with one zero-mass
			// placeholder so the search misses cleanly.
			for (size_t p = 0; p < need; ++p)
			{
				if (r.vvdYionMass[p].empty())
				{
					r.vvdYionMass[p].push_back(0.0);
					r.vvdYionProb[p].push_back(0.0);
				}
				if (r.vvdBionMass[p].empty())
				{
					r.vvdBionMass[p].push_back(0.0);
					r.vvdBionProb[p].push_back(0.0);
				}
			}

			out.push_back(std::move(r));
		}
		return true;
	}
	catch (const H5::Exception &e)
	{
		errOut = e.getCDetailMsg() ? e.getCDetailMsg() : "HDF5 error";
		return false;
	}
	catch (const std::exception &e)
	{
		errOut = e.what();
		return false;
	}
}

// -------------------- Matching --------------------

struct PrecursorMatch
{
	bool ok = false;
	int isotopicShift = 0;
	double mzShift = 0.0;        // (real_mz - sip_mz) Da
	double mzAbsErrPpm = 0.0;
	double deltaRT = 0.0;        // abs(scan.RT - record.rtMinutes) in minutes
	int matchedCharge = 0;
	double matchedMz = 0.0;
};

struct ScanPrecursor
{
	double mz = 0.0;
	int charge = 0;
};

double scanRtMinutes(const MS2Scan *scan)
{
	return parseRetentionMinutes(scan->sRTime);
}

std::vector<ScanPrecursor> scanPrecursors(const MS2Scan *scan)
{
	std::vector<ScanPrecursor> precursors;
	if (!scan->dParentMZs.empty())
	{
		size_t k = std::min(scan->dParentMZs.size(), scan->iParentChargeStates.size());
		precursors.reserve(k);
		for (size_t i = 0; i < k; ++i)
			precursors.push_back({scan->dParentMZs[i], scan->iParentChargeStates[i]});
	}
	if (precursors.empty())
		precursors.push_back({scan->dParentMZ, scan->iParentChargeState});
	return precursors;
}

PrecursorMatch matchPrecursor(const MS2Scan *scan,
							  const SipRecord &rec,
							  const MassTolerance &mzTol,
							  double rtTolMin,
							  int isotopeShiftWindow)
{
	PrecursorMatch best;
	double bestErr = std::numeric_limits<double>::infinity();
	const double neutron = ProNovoConfig::getNeutronMass();
	const double scanRt = scanRtMinutes(scan);
	const double dRt = scanRt - rec.rtMinutes;
	if (std::fabs(dRt) > rtTolMin)
		return best;

	for (const ScanPrecursor &precursor : scanPrecursors(scan))
	{
		double realMz = precursor.mz;
		int realCharge = precursor.charge;
		if (realCharge <= 0)
			continue;
		if (rec.charge > 0 && realCharge != rec.charge)
			continue;
		// Tolerance is on the observed m/z (Da window scaled by ppm when
		// ppm mode). The matched precursor sets the scale.
		const double tolMzDa = mzTol.daAt(realMz);
		// Scale to neutral mass space.
		const double tolNeutralDa = tolMzDa * realCharge;
		double dm = realMz - rec.topPrecursorMz;
		double dmNeutral = dm * realCharge;
		int shift = static_cast<int>(std::round(dmNeutral / neutron));
		if (std::abs(shift) > isotopeShiftWindow)
			continue;
		double resid = dmNeutral - shift * neutron;
		if (std::fabs(resid) > tolNeutralDa)
			continue;
		double err = std::fabs(resid);
		if (err < bestErr)
		{
			bestErr = err;
			best.ok = true;
			best.isotopicShift = shift;
			best.mzShift = dm;
			best.mzAbsErrPpm = (rec.topPrecursorMz > 0) ? std::fabs(resid) / rec.topPrecursorMz * 1e6 : 0.0;
			best.deltaRT = std::fabs(dRt);
			best.matchedCharge = realCharge;
			best.matchedMz = realMz;
		}
	}
	return best;
}

// -------------------- Auxiliary scores --------------------

struct EnvCounts
{
	int matchedY = 0;
	int matchedB = 0;
};

EnvCounts countMatchedEnvelopes(const MS2Scan *scan,
								const SipRecord &rec,
								const MassTolerance &tol)
{
	EnvCounts c;
	const auto &mz = scan->vdMZ;
	if (mz.empty())
		return c;
	const double protonMass = ProNovoConfig::getProtonMass();
	auto envelopeMatches = [&](const std::vector<double> &envMz) -> bool
	{
		for (double em : envMz)
		{
			if (em <= 0.0)
				continue;
			const double expectedMz = em + protonMass;
			const double tolDa = tol.daAt(expectedMz);
			auto it = std::lower_bound(mz.begin(), mz.end(), expectedMz - tolDa);
			while (it != mz.end() && *it <= expectedMz + tolDa)
			{
				if (std::fabs(*it - expectedMz) <= tolDa)
					return true;
				++it;
			}
		}
		return false;
	};
	for (const auto &env : rec.vvdYionMass)
		if (envelopeMatches(env))
			++c.matchedY;
	for (const auto &env : rec.vvdBionMass)
		if (envelopeMatches(env))
			++c.matchedB;
	return c;
}

// Build aligned (p, q) pairs by matching each HDF5 fragment to its nearest
// scan peak within tolDa, then normalize to unit sum.
void alignSpectra(const std::vector<double> &fragMz,
				  const std::vector<double> &fragInt,
				  const std::vector<double> &scanMz,
				  const std::vector<double> &scanIntensity,
				  const MassTolerance &tol,
				  std::vector<double> &p,
				  std::vector<double> &q)
{
	p.clear();
	q.clear();
	if (fragMz.empty() || scanMz.empty())
		return;
	for (size_t i = 0; i < fragMz.size(); ++i)
	{
		double m = fragMz[i];
		const double tolDa = tol.daAt(m);
		auto it = std::lower_bound(scanMz.begin(), scanMz.end(), m - tolDa);
		double bestErr = tolDa + 1.0;
		double bestI = 0.0;
		bool found = false;
		while (it != scanMz.end() && *it <= m + tolDa)
		{
			double err = std::fabs(*it - m);
			if (err < bestErr)
			{
				bestErr = err;
				bestI = scanIntensity[static_cast<size_t>(it - scanMz.begin())];
				found = true;
			}
			++it;
		}
		if (found)
		{
			p.push_back(fragInt[i]);
			q.push_back(bestI);
		}
	}
	auto norm = [](std::vector<double> &v)
	{
		double s = 0.0;
		for (double x : v)
			s += x;
		if (s > 0.0)
			for (double &x : v)
				x /= s;
	};
	norm(p);
	norm(q);
}

// HDF5-driven sparse alignment for cosine: keep every positive HDF5 fragment
// intensity and pair it with the nearest observed scan peak within tolerance,
// or zero when no observed peak is present.
void alignSpectraSparseCosine(const std::vector<double> &fragMz,
							  const std::vector<double> &fragInt,
							  const std::vector<double> &scanMz,
							  const std::vector<double> &scanIntensity,
							  const MassTolerance &tol,
							  std::vector<double> &p,
							  std::vector<double> &q)
{
	p.clear();
	q.clear();
	const size_t nFrag = std::min(fragMz.size(), fragInt.size());
	if (nFrag == 0)
		return;

	const size_t nScan = std::min(scanMz.size(), scanIntensity.size());
	for (size_t i = 0; i < nFrag; ++i)
	{
		const double m = fragMz[i];
		const double intensity = fragInt[i];
		if (intensity <= 0.0 || !std::isfinite(intensity) || !std::isfinite(m))
			continue;

		double bestI = 0.0;
		if (nScan > 0)
		{
			const double tolDa = tol.daAt(m);
			auto it = std::lower_bound(scanMz.begin(), scanMz.begin() + static_cast<std::ptrdiff_t>(nScan),
									   m - tolDa);
			double bestErr = tolDa + 1.0;
			while (it != scanMz.begin() + static_cast<std::ptrdiff_t>(nScan) && *it <= m + tolDa)
			{
				const double err = std::fabs(*it - m);
				if (err < bestErr)
				{
					bestErr = err;
					bestI = scanIntensity[static_cast<size_t>(it - scanMz.begin())];
				}
				++it;
			}
		}

		p.push_back(intensity);
		q.push_back(bestI > 0.0 && std::isfinite(bestI) ? bestI : 0.0);
	}

	auto norm = [](std::vector<double> &v)
	{
		double s = 0.0;
		for (double x : v)
			s += x;
		if (s > 0.0)
			for (double &x : v)
				x /= s;
	};
	norm(p);
	norm(q);
}

double computeCosine(const std::vector<double> &p, const std::vector<double> &q)
{
	double dot = 0.0, np = 0.0, nq = 0.0;
	for (size_t i = 0; i < p.size(); ++i)
	{
		dot += p[i] * q[i];
		np += p[i] * p[i];
		nq += q[i] * q[i];
	}
	if (np <= 0.0 || nq <= 0.0)
		return 0.0;
	return dot / std::sqrt(np * nq);
}

double computeEntropy(const std::vector<double> &p, const std::vector<double> &q)
{
	auto H = [](const std::vector<double> &v)
	{
		double h = 0.0;
		for (double x : v)
			if (x > 0.0)
				h += -x * std::log(x);
		return h;
	};
	if (p.empty() || q.empty())
		return 0.0;
	double Hp = H(p);
	double Hq = H(q);
	std::vector<double> m(p.size());
	for (size_t i = 0; i < p.size(); ++i)
		m[i] = 0.5 * (p[i] + q[i]);
	double Hm = H(m);
	// Jensen-Shannon similarity score, mapped to higher = more similar.
	// JS divergence = Hm - 0.5*Hp - 0.5*Hq. Similarity = 1 - JS/ln(2).
	double js = Hm - 0.5 * Hp - 0.5 * Hq;
	double sim = 1.0 - js / std::log(2.0);
	if (sim < 0.0)
		sim = 0.0;
	if (sim > 1.0)
		sim = 1.0;
	return sim;
}

// -------------------- Per-peptide feature helpers --------------------

int nearestPrecursorPeakIndex(const SipRecord &rec, double mz)
{
	if (rec.precursorMz.empty())
		return -1;
	int bestIdx = -1;
	double bestErr = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < rec.precursorMz.size(); ++i)
	{
		const double err = std::fabs(rec.precursorMz[i] - mz);
		if (err < bestErr)
		{
			bestErr = err;
			bestIdx = static_cast<int>(i);
		}
	}
	return bestIdx;
}

int topPrecursorPeakIndex(const SipRecord &rec)
{
	if (rec.precursorIntensity.empty())
		return nearestPrecursorPeakIndex(rec, rec.topPrecursorMz);
	return static_cast<int>(std::distance(
		rec.precursorIntensity.begin(),
		std::max_element(rec.precursorIntensity.begin(), rec.precursorIntensity.end())));
}

struct ShardPsmRow
{
	int32_t scanIdx = 0;
	int32_t parentCharge = 0;
	int32_t isotopicShift = 0;
	int32_t matchedY = 0;
	int32_t matchedB = 0;
	int32_t peptideLength = 0;
	int32_t missCleavage = 0;
	int32_t ptmCount = 0;
	int32_t isotopicPeakNumbers = 0;
	double ms1IsotopeFitScore = 0.0;
	double wdp = 0.0;
	double xcorr = 0.0;
	double mvh = 0.0;
	double entropy = 0.0;
	double cosine = 0.0;
	double deltaRT = 0.0;
	double mzShiftDa = 0.0;
	double mzAbsErrPpm = 0.0;
	double expMassNeutral = 0.0;
	double calcMassNeutral = 0.0;
	double rtScan = 0.0;
	double ms1IsotopicAbundancePct = 0.0;
	double log10PrecursorIntensity = 0.0;
	double matchedMz = 0.0;
	// strings: psmId, peptide (decorated), proteins, sipAtom — pascal-style
	std::string psmId;
	std::string peptide;
	std::string proteins;
	std::string sipAtom;
};

// -------------------- Parent: PIN rows --------------------

constexpr double kMinPinWdpScore = 0.5;

struct MergedRow
{
	ShardPsmRow row;
	int32_t label = 0;
	double ms2Pct = 0.0;
	int scanNumber = 0; // real MS2 scan ID
	int rank = 0;
};

std::vector<PinWriter::SearchSpectraPinRow> makeSearchSpectraPinRows(const std::vector<MergedRow> &rows)
{
	std::vector<PinWriter::SearchSpectraPinRow> pinRows;
	pinRows.reserve(rows.size());
	for (const MergedRow &m : rows)
	{
		const ShardPsmRow &r = m.row;
		PinWriter::SearchSpectraPinRow row;
		row.label = m.label;
		row.scanNumber = m.scanNumber;
		row.rank = m.rank;
		row.parentCharge = r.parentCharge;
		row.isotopicMassWindowShift = r.isotopicShift;
		row.peptideLength = r.peptideLength;
		row.missCleavageSiteNumber = r.missCleavage;
		row.ptmNumber = r.ptmCount;
		row.isotopicPeakNumber = r.isotopicPeakNumbers;
		row.ms1IsotopeFitScore = r.ms1IsotopeFitScore;
		row.matchedYEnvelope = r.matchedY;
		row.matchedBEnvelope = r.matchedB;
		row.expMass = r.calcMassNeutral;
		row.observedMass = r.expMassNeutral;
		row.retentionTime = r.rtScan;
		row.massError = r.mzAbsErrPpm;
		row.mzShiftFromIsolationWindowCenter = r.mzShiftDa;
		row.ms1IsotopicAbundance = r.ms1IsotopicAbundancePct;
		row.ms2IsotopicAbundance = m.ms2Pct;
		row.wdpScore = r.wdp;
		row.xcorrScore = r.xcorr;
		row.mvhScore = r.mvh;
		row.entropyScore = r.entropy;
		row.cosineScore = r.cosine;
		row.deltaRT = r.deltaRT;
		row.log10PrecursorIntensity = r.log10PrecursorIntensity;
		row.peptide = r.peptide;
		row.proteins = r.proteins;
		pinRows.push_back(std::move(row));
	}
	return pinRows;
}

// -------------------- Parent: per-HDF5 sample driver --------------------

// Holds one in-memory record paired with its source HDF5's label/metadata.
struct LabeledRecord
{
	SipRecord rec;
	int label = 0; // +1 / -1
	double ms2Pct = 0.0;
	std::string sipAtom;
};

bool writeFull(int fd, const void *data, size_t len)
{
	const char *p = static_cast<const char *>(data);
	while (len > 0)
	{
		ssize_t n = ::write(fd, p, len);
		if (n < 0)
		{
			if (errno == EINTR)
				continue;
			return false;
		}
		if (n == 0)
			return false;
		p += n;
		len -= static_cast<size_t>(n);
	}
	return true;
}

bool readFull(int fd, void *data, size_t len)
{
	char *p = static_cast<char *>(data);
	while (len > 0)
	{
		ssize_t n = ::read(fd, p, len);
		if (n < 0)
		{
			if (errno == EINTR)
				continue;
			return false;
		}
		if (n == 0)
			return false;
		p += n;
		len -= static_cast<size_t>(n);
	}
	return true;
}

template <typename T>
bool writePod(int fd, const T &x)
{
	return writeFull(fd, &x, sizeof(T));
}

template <typename T>
bool readPod(int fd, T &x)
{
	return readFull(fd, &x, sizeof(T));
}

bool writeString(int fd, const std::string &s)
{
	const uint64_t n = static_cast<uint64_t>(s.size());
	return writePod(fd, n) && (n == 0 || writeFull(fd, s.data(), static_cast<size_t>(n)));
}

bool readString(int fd, std::string &s)
{
	uint64_t n = 0;
	if (!readPod(fd, n))
		return false;
	s.assign(static_cast<size_t>(n), '\0');
	return n == 0 || readFull(fd, s.data(), static_cast<size_t>(n));
}

int createAnonymousMemfd(const char *name)
{
#ifdef SYS_memfd_create
	return static_cast<int>(::syscall(SYS_memfd_create, name, MFD_CLOEXEC));
#else
	errno = ENOSYS;
	return -1;
#endif
}

struct FlatMappedWriter
{
	unsigned char *data = nullptr;
	size_t size = 0;
	size_t pos = 0;
	bool ok = true;

	template <typename T>
	void pod(const T &x)
	{
		bytes(&x, sizeof(T));
	}

	void bytes(const void *p, size_t n)
	{
		if (!ok || n > size - pos)
		{
			ok = false;
			return;
		}
		if (n > 0)
			std::memcpy(data + pos, p, n);
		pos += n;
	}

	void str(const std::string &s)
	{
		const uint64_t n = static_cast<uint64_t>(s.size());
		pod(n);
		if (n > 0)
			bytes(s.data(), static_cast<size_t>(n));
	}

	void doubles(const std::vector<double> &v)
	{
		const uint64_t n = static_cast<uint64_t>(v.size());
		pod(n);
		if (n > 0)
			bytes(v.data(), static_cast<size_t>(n) * sizeof(double));
	}

	void rows(const std::vector<std::vector<double>> &v)
	{
		const uint64_t n = static_cast<uint64_t>(v.size());
		pod(n);
		for (const auto &row : v)
			doubles(row);
	}
};

struct FlatReader
{
	const unsigned char *data = nullptr;
	size_t size = 0;
	size_t pos = 0;

	bool bytes(void *out, size_t n)
	{
		if (n > size - pos)
			return false;
		std::memcpy(out, data + pos, n);
		pos += n;
		return true;
	}

	template <typename T>
	bool pod(T &x)
	{
		return bytes(&x, sizeof(T));
	}

	bool str(std::string &s)
	{
		uint64_t n = 0;
		if (!pod(n) || n > size - pos)
			return false;
		s.assign(reinterpret_cast<const char *>(data + pos), static_cast<size_t>(n));
		pos += static_cast<size_t>(n);
		return true;
	}

	bool doubles(std::vector<double> &v)
	{
		uint64_t n = 0;
		if (!pod(n))
			return false;
		const size_t bytesNeeded = static_cast<size_t>(n) * sizeof(double);
		if (bytesNeeded > size - pos)
			return false;
		v.resize(static_cast<size_t>(n));
		if (bytesNeeded > 0)
			std::memcpy(v.data(), data + pos, bytesNeeded);
		pos += bytesNeeded;
		return true;
	}

	bool rows(std::vector<std::vector<double>> &v)
	{
		uint64_t n = 0;
		if (!pod(n))
			return false;
		v.resize(static_cast<size_t>(n));
		for (auto &row : v)
			if (!doubles(row))
				return false;
		return true;
	}
};

size_t flatStringSize(const std::string &s)
{
	return sizeof(uint64_t) + s.size();
}

size_t flatDoublesSize(const std::vector<double> &v)
{
	return sizeof(uint64_t) + v.size() * sizeof(double);
}

size_t flatRowsSize(const std::vector<std::vector<double>> &v)
{
	size_t n = sizeof(uint64_t);
	for (const auto &row : v)
		n += flatDoublesSize(row);
	return n;
}

size_t flatSipRecordSize(const SipRecord &r)
{
	size_t n = 0;
	n += flatStringSize(r.psmId);
	n += flatStringSize(r.retention);
	n += flatStringSize(r.peptide);
	n += flatStringSize(r.nakedPeptide);
	n += flatStringSize(r.proteins);
	n += sizeof(r.charge);
	n += sizeof(r.rtMinutes);
	n += sizeof(r.topPrecursorMz);
	n += sizeof(r.topPrecursorIntensity);
	n += sizeof(r.sumPrecursorIntensity);
	n += flatDoublesSize(r.precursorMz);
	n += flatDoublesSize(r.precursorIntensity);
	n += flatRowsSize(r.vvdYionMass);
	n += flatRowsSize(r.vvdYionProb);
	n += flatRowsSize(r.vvdBionMass);
	n += flatRowsSize(r.vvdBionProb);
	n += flatDoublesSize(r.fragMz);
	n += flatDoublesSize(r.fragExpInt);
	n += sizeof(uint64_t);
	n += sizeof(uint64_t);
	return n;
}

size_t flatHdf5BlobSize(const std::vector<SipRecord> &records,
						const Hdf5FileMeta &meta)
{
	size_t n = 0;
	n += sizeof(meta.formatVersion);
	n += flatStringSize(meta.chemistryProfileId);
	n += sizeof(meta.targetSipAbundancePct);
	n += flatStringSize(meta.sipAtom);
	n += sizeof(meta.sipIsotopeMassNumber);
	n += sizeof(uint64_t);
	for (const auto &record : records)
		n += flatSipRecordSize(record);
	return n;
}

template <typename Writer>
void writeFlatSipRecord(Writer &w, const SipRecord &r)
{
	const uint64_t nY = static_cast<uint64_t>(r.nYpositions);
	const uint64_t nB = static_cast<uint64_t>(r.nBpositions);
	w.str(r.psmId);
	w.str(r.retention);
	w.str(r.peptide);
	w.str(r.nakedPeptide);
	w.str(r.proteins);
	w.pod(r.charge);
	w.pod(r.rtMinutes);
	w.pod(r.topPrecursorMz);
	w.pod(r.topPrecursorIntensity);
	w.pod(r.sumPrecursorIntensity);
	w.doubles(r.precursorMz);
	w.doubles(r.precursorIntensity);
	w.rows(r.vvdYionMass);
	w.rows(r.vvdYionProb);
	w.rows(r.vvdBionMass);
	w.rows(r.vvdBionProb);
	w.doubles(r.fragMz);
	w.doubles(r.fragExpInt);
	w.pod(nY);
	w.pod(nB);
}

template <typename Writer>
void writeFlatHdf5Blob(Writer &w,
					   const std::vector<SipRecord> &records,
					   const Hdf5FileMeta &meta)
{
	w.pod(meta.formatVersion);
	w.str(meta.chemistryProfileId);
	w.pod(meta.targetSipAbundancePct);
	w.str(meta.sipAtom);
	w.pod(meta.sipIsotopeMassNumber);
	const uint64_t n = static_cast<uint64_t>(records.size());
	w.pod(n);
	for (const auto &record : records)
		writeFlatSipRecord(w, record);
}

bool readFlatSipRecord(FlatReader &r, SipRecord &out)
{
	uint64_t nY = 0;
	uint64_t nB = 0;
	if (!r.str(out.psmId) ||
		!r.str(out.retention) ||
		!r.str(out.peptide) ||
		!r.str(out.nakedPeptide) ||
		!r.str(out.proteins) ||
		!r.pod(out.charge) ||
		!r.pod(out.rtMinutes) ||
		!r.pod(out.topPrecursorMz) ||
		!r.pod(out.topPrecursorIntensity) ||
		!r.pod(out.sumPrecursorIntensity) ||
		!r.doubles(out.precursorMz) ||
		!r.doubles(out.precursorIntensity) ||
		!r.rows(out.vvdYionMass) ||
		!r.rows(out.vvdYionProb) ||
		!r.rows(out.vvdBionMass) ||
		!r.rows(out.vvdBionProb) ||
		!r.doubles(out.fragMz) ||
		!r.doubles(out.fragExpInt) ||
		!r.pod(nY) ||
		!r.pod(nB))
		return false;
	out.nYpositions = static_cast<size_t>(nY);
	out.nBpositions = static_cast<size_t>(nB);
	return true;
}

bool readFlatHdf5Blob(const void *data,
					  size_t size,
					  std::vector<SipRecord> &records,
					  Hdf5FileMeta &meta)
{
	FlatReader r{static_cast<const unsigned char *>(data), size, 0};
	uint64_t n = 0;
	if (!r.pod(meta.formatVersion) ||
		!r.str(meta.chemistryProfileId) ||
		!r.pod(meta.targetSipAbundancePct) ||
		!r.str(meta.sipAtom) ||
		!r.pod(meta.sipIsotopeMassNumber) ||
		!r.pod(n))
		return false;
	records.clear();
	records.reserve(static_cast<size_t>(n));
	for (uint64_t i = 0; i < n; ++i)
	{
		SipRecord record;
		if (!readFlatSipRecord(r, record))
			return false;
		records.push_back(std::move(record));
	}
	return r.pos == r.size;
}

bool writeFlatHdf5BlobToSharedMemory(int fd,
									 const std::vector<SipRecord> &records,
									 const Hdf5FileMeta &meta,
									 uint64_t &blobSize,
									 std::string &err)
{
	const size_t n = flatHdf5BlobSize(records, meta);
	blobSize = static_cast<uint64_t>(n);
	if (::ftruncate(fd, static_cast<off_t>(n)) != 0)
	{
		err = std::string("ftruncate failed: ") + std::strerror(errno);
		return false;
	}
	if (n == 0)
		return true;
	void *p = ::mmap(nullptr, n, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (p == MAP_FAILED)
	{
		err = std::string("mmap write failed: ") + std::strerror(errno);
		return false;
	}
	FlatMappedWriter w{static_cast<unsigned char *>(p), n, 0, true};
	writeFlatHdf5Blob(w, records, meta);
	const bool ok = w.ok && w.pos == n;
	::munmap(p, n);
	if (!ok)
	{
		err = "flat HDF5 shared-memory write overflow";
		return false;
	}
	return true;
}

bool readBlobFromSharedMemory(int fd,
							  size_t blobSize,
							  std::vector<SipRecord> &records,
							  Hdf5FileMeta &meta,
							  std::string &err)
{
	if (blobSize == 0)
	{
		err = "empty HDF5 shared-memory blob";
		return false;
	}
	void *p = ::mmap(nullptr, blobSize, PROT_READ, MAP_SHARED, fd, 0);
	if (p == MAP_FAILED)
	{
		err = std::string("mmap read failed: ") + std::strerror(errno);
		return false;
	}
	const bool ok = readFlatHdf5Blob(p, blobSize, records, meta);
	::munmap(p, blobSize);
	if (!ok)
	{
		err = "cannot parse flat HDF5 shared-memory blob";
		return false;
	}
	if (meta.formatVersion != kSpectraHdf5FormatVersion)
	{
		err = "child metadata carries unsupported spectra format_version " +
			std::to_string(meta.formatVersion);
		return false;
	}
	const std::string expectedChemistryProfileId =
		ProNovoConfig::getChemistryProfileId();
	if (meta.chemistryProfileId != expectedChemistryProfileId)
	{
		err = "child metadata carries incompatible chemistry_profile_id '" +
			meta.chemistryProfileId + "'; expected '" +
			expectedChemistryProfileId + "'";
		return false;
	}
	const std::string actualSipIsotope =
		PSMfeatureExtractor::canonicalSipIsotope(
			meta.sipAtom, meta.sipIsotopeMassNumber);
	const std::string expectedSipIsotope =
		PSMfeatureExtractor::canonicalSipIsotope(
			ProNovoConfig::getSetSIPelement());
	if (actualSipIsotope.empty() || actualSipIsotope != expectedSipIsotope)
	{
		err = "child metadata carries incompatible SIP isotope atom='" +
			meta.sipAtom + "', mass_number=" +
			std::to_string(meta.sipIsotopeMassNumber) +
			"; expected '" + expectedSipIsotope + "'";
		return false;
	}
	return true;
}

bool writeChildStatus(int fd,
					  bool ok,
					  uint64_t blobSize,
					  const std::string &err)
{
	const int32_t status = ok ? 0 : 1;
	return writePod(fd, status) &&
		   writePod(fd, blobSize) &&
		   (!ok ? writeString(fd, err) : true);
}

bool readChildStatus(int fd,
					 bool &ok,
					 uint64_t &blobSize,
					 std::string &err)
{
	int32_t status = 1;
	if (!readPod(fd, status) ||
		!readPod(fd, blobSize))
		return false;
	ok = status == 0;
	if (!ok && !readString(fd, err))
		return false;
	return true;
}

void appendLabeledRecords(const std::string &path,
						  std::vector<SipRecord> &records,
						  const Hdf5FileMeta &meta,
						  std::vector<LabeledRecord> &out)
{
	out.reserve(out.size() + records.size());
	const int label = isDecoyHdf5(path) ? -1 : +1;
	for (auto &record : records)
	{
		LabeledRecord lr;
		lr.rec = std::move(record);
		lr.label = label;
		lr.ms2Pct = meta.targetSipAbundancePct;
		lr.sipAtom = meta.sipAtom;
		out.push_back(std::move(lr));
	}
}

struct Hdf5LoadChild
{
	size_t fileIdx = 0;
	pid_t pid = -1;
	int readFd = -1;
	int shmFd = -1;
};

struct Hdf5ReadyBlob
{
	std::string path;
	int shmFd = -1;
	uint64_t blobSize = 0;
};

struct Hdf5DeserializeResult
{
	bool ok = false;
	std::string error;
	std::vector<LabeledRecord> records;
};

// Load a bounded set of HDF5 files directly into RAM. Nothing is written to
// the output directory until the final PIN file is materialized. Each file in
// this batch is loaded by a separate child process so HDF5's process-local
// global lock does not serialize the batch. Large HDF5 payloads come back as
// one flat blob in anonymous shared memory; the pipe carries only status.
int loadHdf5Batch(const std::vector<std::string> &hdf5Files,
				  size_t begin,
				  size_t end,
				  std::vector<LabeledRecord> &out)
{
	int failures = 0;
	std::vector<Hdf5LoadChild> children;
	for (size_t i = begin; i < end && i < hdf5Files.size(); ++i)
	{
		int fds[2] = {-1, -1};
		if (::pipe(fds) != 0)
		{
			std::cerr << "Cannot create pipe for " << hdf5Files[i]
					  << ": " << std::strerror(errno) << "\n";
			++failures;
			continue;
		}
		int shmFd = createAnonymousMemfd("sipros_hdf5_batch");
		if (shmFd < 0)
		{
			std::cerr << "Cannot create shared-memory blob for " << hdf5Files[i]
					  << ": " << std::strerror(errno) << "\n";
			::close(fds[0]);
			::close(fds[1]);
			++failures;
			continue;
		}

		pid_t pid = ::fork();
		if (pid < 0)
		{
			std::cerr << "Cannot fork HDF5 reader for " << hdf5Files[i]
					  << ": " << std::strerror(errno) << "\n";
			::close(fds[0]);
			::close(fds[1]);
			::close(shmFd);
			++failures;
			continue;
		}

		if (pid == 0)
		{
			::close(fds[0]);
			std::vector<SipRecord> records;
			Hdf5FileMeta meta;
			std::string err;
			const bool ok = loadHdf5File(hdf5Files[i], records, meta, err);
			uint64_t blobSize = 0;
			bool copied = false;
			if (ok)
			{
				copied = writeFlatHdf5BlobToSharedMemory(shmFd, records, meta,
														 blobSize, err);
			}
			const bool wrote = writeChildStatus(fds[1], ok && copied, blobSize, err);
			::close(fds[1]);
			::close(shmFd);
			::_exit(wrote ? 0 : 2);
		}

		::close(fds[1]);
		Hdf5LoadChild child;
		child.fileIdx = i;
		child.pid = pid;
		child.readFd = fds[0];
		child.shmFd = shmFd;
		children.push_back(child);
	}

	std::vector<Hdf5ReadyBlob> readyBlobs;
	readyBlobs.reserve(children.size());
	for (const Hdf5LoadChild &child : children)
	{
		std::string err;
		bool childOk = false;
		uint64_t blobSize = 0;
		const bool statusOk = readChildStatus(child.readFd, childOk, blobSize, err);
		::close(child.readFd);

		int status = 0;
		while (::waitpid(child.pid, &status, 0) < 0 && errno == EINTR)
		{
		}

		const std::string &path = hdf5Files[child.fileIdx];
		if (!statusOk)
		{
			std::cerr << "Cannot load " << path
					  << ": failed to read child process status\n";
			::close(child.shmFd);
			++failures;
			continue;
		}
		if (!childOk)
		{
			std::cerr << "Cannot load " << path << ": " << err << "\n";
			::close(child.shmFd);
			++failures;
			continue;
		}
		if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
		{
			std::cerr << "Cannot load " << path
					  << ": HDF5 reader process exited abnormally\n";
			::close(child.shmFd);
			++failures;
			continue;
		}

		Hdf5ReadyBlob ready;
		ready.path = path;
		ready.shmFd = child.shmFd;
		ready.blobSize = blobSize;
		readyBlobs.push_back(std::move(ready));
	}

	if (!readyBlobs.empty())
	{
		std::vector<Hdf5DeserializeResult> results(readyBlobs.size());
#pragma omp parallel for schedule(dynamic, 1)
		for (size_t i = 0; i < readyBlobs.size(); ++i)
		{
			const Hdf5ReadyBlob &ready = readyBlobs[i];
			Hdf5DeserializeResult result;

			std::vector<SipRecord> records;
			Hdf5FileMeta meta;
			std::string err;
			if (readBlobFromSharedMemory(ready.shmFd, static_cast<size_t>(ready.blobSize),
										 records, meta, err))
			{
				appendLabeledRecords(ready.path, records, meta, result.records);
				result.ok = true;
			}
			else
			{
				result.error = "Cannot load " + ready.path + ": " + err;
			}
			::close(ready.shmFd);
			results[i] = std::move(result);
		}

		size_t recordsToAppend = 0;
		for (const Hdf5DeserializeResult &result : results)
		{
			if (!result.ok)
			{
				std::cerr << result.error << "\n";
				++failures;
				continue;
			}
			recordsToAppend += result.records.size();
		}

		out.reserve(out.size() + recordsToAppend);
		for (Hdf5DeserializeResult &result : results)
		{
			if (!result.ok)
				continue;
			out.insert(out.end(),
					   std::make_move_iterator(result.records.begin()),
					   std::make_move_iterator(result.records.end()));
		}
	}
	return failures;
}

// Per-scan scoring state.
struct ScoringPsm
{
	size_t scanIdx = 0;
	int scanNumber = 0;
	int label = 0;
	double wdp = -std::numeric_limits<double>::infinity();
	ShardPsmRow row;
	double ms2Pct = 0.0;
	std::string nakedPeptide;
	std::string sipAtom;
	PrecursorMatch match;
	SipRecord rec;
};

struct CandidateMatch
{
	size_t recordIdx = 0;
	PrecursorMatch match;
};

struct MzIndexedRecord
{
	double mz = 0.0;
	size_t recordIdx = 0;
};

bool wouldKeepTopUniquePeptide(const std::vector<ScoringPsm> &top,
							   const std::string &nakedPeptide,
							   double wdp,
							   int limit)
{
	auto samePeptide = std::find_if(top.begin(), top.end(),
									[&](const ScoringPsm &x)
									{ return x.nakedPeptide == nakedPeptide; });
	if (samePeptide != top.end())
		return true;
	return top.size() < static_cast<size_t>(limit) || wdp > top.back().wdp;
}

void pushTopUniquePeptide(std::vector<ScoringPsm> &top, ScoringPsm &&psm, int limit)
{
	auto byWdpDesc = [](const ScoringPsm &a, const ScoringPsm &b)
	{ return a.wdp > b.wdp; };

	auto samePeptide = std::find_if(top.begin(), top.end(),
									[&](const ScoringPsm &x)
									{ return x.nakedPeptide == psm.nakedPeptide; });
	if (samePeptide != top.end())
	{
		psm.rec.proteins = mergeProteinLists(samePeptide->rec.proteins, psm.rec.proteins);
		if (psm.wdp > samePeptide->wdp)
		{
			*samePeptide = std::move(psm);
			std::sort(top.begin(), top.end(), byWdpDesc);
		}
		else
		{
			samePeptide->rec.proteins = std::move(psm.rec.proteins);
		}
		return;
	}

	if (top.size() < static_cast<size_t>(limit))
	{
		top.push_back(std::move(psm));
		std::sort(top.begin(), top.end(), byWdpDesc);
	}
	else if (psm.wdp > top.back().wdp)
	{
		top.back() = std::move(psm);
		std::sort(top.begin(), top.end(), byWdpDesc);
	}
}

std::vector<MzIndexedRecord> buildMzIndex(const std::vector<LabeledRecord> &records)
{
	std::vector<MzIndexedRecord> index;
	index.reserve(records.size());
	for (size_t r = 0; r < records.size(); ++r)
	{
		const SipRecord &rec = records[r].rec;
		if (rec.charge <= 0 || rec.nakedPeptide.empty() || !std::isfinite(rec.topPrecursorMz))
			continue;
		index.push_back({rec.topPrecursorMz, r});
	}
	std::sort(index.begin(), index.end(),
			  [](const MzIndexedRecord &a, const MzIndexedRecord &b)
			  { return a.mz < b.mz; });
	return index;
}

size_t assignCandidatesToScans(const std::vector<MS2Scan *> &scans,
							   const std::vector<LabeledRecord> &records,
							   const std::vector<MzIndexedRecord> &mzIndex,
							   const MassTolerance &mzTol,
							   double rtTol,
							   int isotopeShiftWindow,
							   std::vector<std::vector<CandidateMatch>> &out)
{
	out.clear();
	out.resize(scans.size());

	size_t assignedCount = 0;
	const double neutron = ProNovoConfig::getNeutronMass();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : assignedCount)
	for (size_t s = 0; s < scans.size(); ++s)
	{
		MS2Scan *scan = scans[s];
		const double scanRt = scanRtMinutes(scan);
		std::vector<size_t> candidateIdx;

		for (const ScanPrecursor &precursor : scanPrecursors(scan))
		{
			if (precursor.charge <= 0)
				continue;
			const double tolMzDa = mzTol.daAt(precursor.mz);
			for (int shift = -isotopeShiftWindow; shift <= isotopeShiftWindow; ++shift)
			{
				const double center = precursor.mz - shift * neutron / precursor.charge;
				const double lo = center - tolMzDa;
				const double hi = center + tolMzDa;
				auto first = std::lower_bound(
					mzIndex.begin(), mzIndex.end(), lo,
					[](const MzIndexedRecord &entry, double value)
					{ return entry.mz < value; });
				for (auto it = first; it != mzIndex.end() && it->mz <= hi; ++it)
				{
					const SipRecord &rec = records[it->recordIdx].rec;
					if (rec.charge == precursor.charge &&
						std::fabs(scanRt - rec.rtMinutes) <= rtTol)
						candidateIdx.push_back(it->recordIdx);
				}
			}
		}

		std::sort(candidateIdx.begin(), candidateIdx.end());
		candidateIdx.erase(std::unique(candidateIdx.begin(), candidateIdx.end()),
						   candidateIdx.end());

		auto &assigned = out[s];
		assigned.reserve(candidateIdx.size());
		for (size_t r : candidateIdx)
		{
			PrecursorMatch match = matchPrecursor(scan, records[r].rec, mzTol, rtTol,
												  isotopeShiftWindow);
			if (match.ok)
				assigned.push_back({r, match});
		}
		assignedCount += assigned.size();
	}
	return assignedCount;
}

void keepTopNEnvelopePeaks(std::vector<std::vector<double>> &masses,
						   std::vector<std::vector<double>> &probs,
						   int topN)
{
	for (size_t i = 0; i < masses.size(); ++i)
	{
		const size_t n = std::min(masses[i].size(), probs[i].size());
		if (n == 0)
		{
			masses[i].assign(1, 0.0);
			probs[i].assign(1, 0.0);
			continue;
		}

		std::vector<size_t> order(n);
		for (size_t j = 0; j < n; ++j)
			order[j] = j;
		std::stable_sort(order.begin(), order.end(),
						 [&](size_t a, size_t b)
						 { return probs[i][a] > probs[i][b]; });

		const size_t keep = std::min(static_cast<size_t>(topN), n);
		std::vector<double> keptMasses;
		std::vector<double> keptProbs;
		keptMasses.reserve(keep);
		keptProbs.reserve(keep);
		double totalProb = 0.0;
		for (size_t j = 0; j < keep; ++j)
		{
			const size_t idx = order[j];
			keptMasses.push_back(masses[i][idx]);
			keptProbs.push_back(probs[i][idx]);
			totalProb += probs[i][idx];
		}
		if (totalProb > 0.0)
			for (double &p : keptProbs)
				p /= totalProb;

		masses[i] = std::move(keptMasses);
		probs[i] = std::move(keptProbs);
	}
}

void trimIonVectorsForXcorr(const MS2Scan *scan,
							std::vector<std::vector<double>> &masses,
							std::vector<std::vector<double>> &probs)
{
	if (!scan->pQuery)
		return;

	const int scanArraySize = scan->pQuery->_spectrumInfoInternal.iArraySize;
	const double invBinWidth = scan->isMS2HighRes
								   ? ProNovoConfig::dHighResInverseBinWidth
								   : ProNovoConfig::dLowResInverseBinWidth;
	const double oneMinusBinOffset = scan->isMS2HighRes
										 ? ProNovoConfig::dHighResOneMinusBinOffset
										 : ProNovoConfig::dLowResOneMinusBinOffset;
	const int maxFragCharge = scan->pQuery->_spectrumInfoInternal.iMaxFragCharge > 0
								  ? scan->pQuery->_spectrumInfoInternal.iMaxFragCharge
								  : 1;

	for (size_t i = 0; i < masses.size(); ++i)
	{
		std::vector<double> keptMasses;
		std::vector<double> keptProbs;
		keptMasses.reserve(masses[i].size());
		keptProbs.reserve(probs[i].size());

		for (size_t j = 0; j < masses[i].size(); ++j)
		{
			bool keep = true;
			for (int z = 1; z <= maxFragCharge; ++z)
			{
				double mz = (masses[i][j] + z * 1.00727646688) / z;
				int bin = static_cast<int>(mz * invBinWidth + oneMinusBinOffset);
				if (bin < 0 || bin >= scanArraySize)
				{
					keep = false;
					break;
				}
			}
			if (keep)
			{
				keptMasses.push_back(masses[i][j]);
				keptProbs.push_back(probs[i][j]);
			}
		}

		if (keptMasses.empty())
		{
			keptMasses.push_back(0.0);
			keptProbs.push_back(0.0);
		}
		masses[i] = std::move(keptMasses);
		probs[i] = std::move(keptProbs);
	}
}

ShardPsmRow makeScoringRow(size_t scanIdx,
                           const MS2Scan *scan,
                           const SipRecord &rec,
                           const LabeledRecord &labeledRecord,
                           const PrecursorMatch &match,
                           const EnvCounts &envCounts,
                           const sipros::RaxportMs1Data *ms1Data,
                           const MassTolerance &mzTol,
                           double wdp,
                           double xcorr,
                           double mvh,
                           double entropy,
                           double cosine)
{
    double baseMass = 0.0;
    const int precursorCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
    std::vector<isotopicPeak> ms1Peaks;
    if (ms1Data && scan && scan->iParentScanID > 0 && precursorCharge > 0)
    {
        PeptideIsotopeCalculator calculator;
        const std::string peptideForComposition = PSMfeatureExtractor::peptideBodyWithPtms(rec.peptide);
        baseMass = calculator.calPrecursorBaseMass(peptideForComposition);
        const double monoPrecursorMz =
            baseMass / precursorCharge + ProNovoConfig::getProtonMass();
        int ms1ScanNumber = scan->iParentScanID;
        const auto mzToleranceDaAt = [&](double mz)
        { return mzTol.daAt(mz); };
        ms1Peaks = PSMfeatureExtractor::findMs1IsotopicPeaks(
            ms1Data,
            ms1ScanNumber,
            precursorCharge,
            monoPrecursorMz,
            match.matchedMz,
            calculator.pepComposition,
            labeledRecord.sipAtom,
            labeledRecord.ms2Pct,
            mzToleranceDaAt);
    }
    const PSMfeatureExtractor::Ms1AbundanceResult ms1Abundance =
        PSMfeatureExtractor::getSIPelementAbundanceFromMS1Peaks(
            ms1Peaks, baseMass, rec.peptide, precursorCharge, labeledRecord.sipAtom,
            labeledRecord.ms2Pct);

    ShardPsmRow row;
    row.scanIdx = static_cast<int32_t>(scanIdx);
    row.parentCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
    row.isotopicShift = match.isotopicShift;
    row.matchedY = envCounts.matchedY;
    row.matchedB = envCounts.matchedB;
    row.peptideLength = static_cast<int>(rec.nakedPeptide.size());
    row.missCleavage = PSMfeatureExtractor::countMissCleavage(rec.nakedPeptide);
    row.ptmCount = PSMfeatureExtractor::countPTM(rec.peptide);
    row.isotopicPeakNumbers =
        ms1Abundance.rawIsotopicPeakCount;
    row.ms1IsotopeFitScore = ms1Abundance.fitScore;
    row.wdp = wdp;
    row.xcorr = xcorr;
    row.mvh = mvh;
    row.entropy = entropy;
    row.cosine = cosine;
    row.deltaRT = match.deltaRT;
    row.mzShiftDa = match.mzShift;
    row.mzAbsErrPpm = match.mzAbsErrPpm;
    row.expMassNeutral = match.matchedMz * match.matchedCharge - match.matchedCharge * ProNovoConfig::getProtonMass();
    row.calcMassNeutral = rec.topPrecursorMz * rec.charge - rec.charge * ProNovoConfig::getProtonMass();
    row.rtScan = scanRtMinutes(scan);
    row.ms1IsotopicAbundancePct = ms1Abundance.abundancePct;
    row.log10PrecursorIntensity = (rec.sumPrecursorIntensity > 0.0)
                                      ? std::log10(rec.sumPrecursorIntensity)
                                      : 0.0;
    row.matchedMz = match.matchedMz;
    row.psmId = rec.psmId;
    row.peptide = rec.peptide;
    row.proteins = rec.proteins;
    row.sipAtom = labeledRecord.sipAtom;
    return row;
}

int processOneHdf5(const Args &args, const std::string &scanPath,
                   const std::vector<std::string> &hdf5Files)
{
	const double processWallStart = omp_get_wtime();
	const double processCpuStart = processTreeCpuSeconds();
	TimingLogger timing;

	const std::string sampleBasename = fs::path(scanPath).stem().string();

	// ppm mode is approximated to a fixed Da window at representative m/z
	// (1000 m/z for both MS1 and MS2). My matchPrecursor / alignSpectra /
	// countMatchedEnvelopes use the ppm value DIRECTLY (per-peak scaling),
	// so they remain exact ppm; only the upstream scoring funcs see the
	// approximated Da window.
	const double parentTolRefMz = 1000.0;
	const double fragTolRefMz = 1000.0;
	const double parentTolDa = args.toleranceMs1Ppm
								   ? args.toleranceMs1 * parentTolRefMz / 1.0e6
								   : args.toleranceMs1;
	const double fragmentToleranceDa = args.toleranceMs2Ppm
									? args.toleranceMs2 * fragTolRefMz / 1.0e6
									: args.toleranceMs2;

	if (!ProNovoConfig::load(ProNovoConfig::Profile::Sip))
	{
		std::cerr << "Cannot initialize the built-in SIP profile\n";
		return 2;
	}
	std::string libraryMetadataError;
	std::string canonicalSipIsotope;
	if (!validateSpectraLibraryMetadata(
			hdf5Files, canonicalSipIsotope, libraryMetadataError))
	{
		std::cerr << libraryMetadataError << "\n";
		return 2;
	}
	if (!ProNovoConfig::selectSipTarget(
			canonicalSipIsotope[0],
			std::stoi(canonicalSipIsotope.substr(1)),
			libraryMetadataError))
	{
		std::cerr << "Cannot activate spectra-library SIP isotope "
				  << canonicalSipIsotope << ": " << libraryMetadataError << "\n";
		return 2;
	}
	ProNovoConfig::setMassAccuracy(parentTolDa, fragmentToleranceDa);

	std::cout << "  Chemistry profile: "
			  << ProNovoConfig::getChemistryProfileId() << "\n";
	std::cout << "  Tolerance: MS1=" << args.toleranceMs1
			  << (args.toleranceMs1Ppm ? " ppm" : " Da")
			  << ", MS2=" << args.toleranceMs2
			  << (args.toleranceMs2Ppm ? " ppm" : " Da")
				  << "  (scoring tolerances: parent=" << parentTolDa
			  << " Da, fragment=" << fragmentToleranceDa << " Da)\n";
	std::cout << "  Xcorr/MVH envelope top-N peaks: "
			  << args.scoreEnvelopeTopN << "\n";
	std::cout << "  WDP top PSMs per scan/label: "
			  << args.topPsmsPerScan << "\n";
	std::cout << "  Precursor isotope shift window: +/-"
			  << args.isotopeShiftWindow << "\n";

	// -------- Read HDF5 scans and score assigned SIP spectra --------
	const int nThreads = std::max(1, args.threads);
	omp_set_num_threads(nThreads);
	timing.printHeader();

	std::vector<MS2Scan *> scans;
	sipros::RaxportMs1Data ms1Data;
	std::string readError;
	bool readHdf5Ok = false;
	timing.run("Read HDF5", "read HDF5 MS1/MS2 scans", 0, "", [&]()
			   { readHdf5Ok = sipros::readRaxportHdf5Scans(scanPath, scans, &ms1Data, readError); });
	if (!readHdf5Ok)
	{
		std::cerr << readError << "\n";
		return 3;
	}
	if (scans.empty())
	{
		std::cerr << "No ms_order == 2 scans found in Raxport HDF5 file: " << scanPath << "\n";
		return 3;
	}
	for (MS2Scan *scan : scans)
	{
		if (scan->dParentMass > ProNovoConfig::dMaxMS2ScanMass)
		{
			ProNovoConfig::dMaxMS2ScanMass = scan->dParentMass;
		}
		if (scan->iParentChargeState > ProNovoConfig::iMaxPercusorCharge)
		{
			ProNovoConfig::iMaxPercusorCharge = scan->iParentChargeState;
		}
		const size_t nPrecursors = std::min(scan->iParentChargeStates.size(), scan->dParentMZs.size());
		for (size_t i = 0; i < nPrecursors; ++i)
		{
			const int charge = scan->iParentChargeStates[i];
			const double chargedMass = scan->dParentMZs[i] * charge;
			if (chargedMass > ProNovoConfig::dMaxMS2ScanMass)
			{
				ProNovoConfig::dMaxMS2ScanMass = chargedMass;
			}
			if (charge > ProNovoConfig::iMaxPercusorCharge)
			{
				ProNovoConfig::iMaxPercusorCharge = charge;
			}
		}
	}

	const sipros::RaxportMs1Data *ms1DataPtr = &ms1Data;
	std::cout << "HDF5 sample: " << scanPath << "  (" << scans.size()
			  << " MS2 scans, " << ms1Data.scans.size() << " MS1 scans)\n";
	std::cout << "HDF5 spectra library: " << (fs::path(args.hdf5Dir) / "*.h5").string()
			  << " (" << hdf5Files.size() << " percentages)\n";
	auto preprocessScans = [&]()
	{
#pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < scans.size(); ++i)
			scans[i]->preprocess();
	};
	timing.run("Preprocess HDF5", "preprocess HDF5 MS2 scans",
			   scans.size(), "scans", preprocessScans);

	MassTolerance mzTol{args.toleranceMs1Ppm, args.toleranceMs1};
	MassTolerance fragTol{args.toleranceMs2Ppm, args.toleranceMs2};
	const double fragTolDa = ProNovoConfig::getMassAccuracyFragmentIon();
	const double rtTol = args.rtToleranceMin;

	// MVH ln-table — must be called once before parallel region.
	double totalPeakSpace = ProNovoConfig::maxObservedMz - ProNovoConfig::minObservedMz;
	int totalPeakBins = static_cast<int>(std::round(totalPeakSpace / (fragTolDa * 2.0)));
	if (totalPeakBins < 1)
		totalPeakBins = 1;
	MVH::initialLnTable(totalPeakBins);

	// Xcorr buffer sizing
	int iArraySizePreprocess = static_cast<int>((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0) * ProNovoConfig::dHighResInverseBinWidth);
	int iArraySizeScoreSip = static_cast<int>((ProNovoConfig::dMaxMS2ScanMass * 2 + 100) * ProNovoConfig::dHighResInverseBinWidth);
	int iArraySizeScore = static_cast<int>((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	// Per-thread scratch
	std::vector<std::unique_ptr<double[]>> tRaw(nThreads), tFXc(nThreads), tCorr(nThreads), tSmooth(nThreads), tPeak(nThreads);
	std::vector<std::vector<bool>> tDuplSip(nThreads);
	std::vector<std::vector<double>> tBinIonSip(nThreads);
	std::vector<std::vector<int>> tBinSip(nThreads);
	std::vector<std::vector<double>> tMvhIonMasses(nThreads);
	std::vector<std::vector<char>> tMvhSeqs(nThreads);
	std::vector<std::unique_ptr<multimap<double, double>>> tMvhSorted(nThreads);
	for (int i = 0; i < nThreads; ++i)
	{
		tRaw[i].reset(new double[iArraySizePreprocess]());
		tFXc[i].reset(new double[iArraySizePreprocess]());
		tCorr[i].reset(new double[iArraySizePreprocess]());
		tSmooth[i].reset(new double[iArraySizePreprocess]());
		tPeak[i].reset(new double[iArraySizePreprocess]());
		tDuplSip[i].assign(iArraySizeScoreSip, false);
		tBinIonSip[i].assign(iArraySizeScoreSip, 0.0);
		tMvhSorted[i] = std::make_unique<multimap<double, double>>();
	}

	auto prepareScanScoring = [&]()
	{
#pragma omp parallel for schedule(dynamic, 1)
		for (size_t s = 0; s < scans.size(); ++s)
		{
			const int tid = omp_get_thread_num();
			MS2Scan *scan = scans[s];

			// Per-scan Xcorr preprocess (Query allocated, owned by scan).
			scan->preprocessMvh(tMvhSorted[tid].get());

			Query *pQuery = new Query();
			if (CometSearchMod::Preprocess(pQuery, scan,
										   tRaw[tid].get(), tFXc[tid].get(),
										   tCorr[tid].get(), tSmooth[tid].get(),
										   tPeak[tid].get()))
			{
				scan->pQuery = pQuery;
			}
			else
			{
				delete pQuery;
				scan->pQuery = nullptr;
			}
		}
	};
	timing.run("Prepare scan scoring", "prepare scan scoring",
			   scans.size(), "scans", prepareScanScoring);

	// Per-scan WDP winners split by label. Entries are self-contained so the
	// current HDF5 cache batch can be freed after each scoring pass.
	std::vector<std::vector<ScoringPsm>> topTarget(scans.size());
	std::vector<std::vector<ScoringPsm>> topDecoy(scans.size());

	const size_t cacheBatchSize = static_cast<size_t>(std::max(1, args.threads));
	size_t totalRecordsScored = 0;
	size_t totalAssignedCandidates = 0;
	size_t batchesScored = 0;
	int hdf5LoadFailures = 0;

	for (size_t batchBegin = 0; batchBegin < hdf5Files.size(); batchBegin += cacheBatchSize)
	{
		const size_t batchEnd = std::min(hdf5Files.size(), batchBegin + cacheBatchSize);
		std::vector<LabeledRecord> batchRecords;
		++batchesScored;
		auto loadBatch = [&]()
		{
			hdf5LoadFailures += loadHdf5Batch(hdf5Files, batchBegin, batchEnd,
											  batchRecords);
		};
		const size_t totalBatches = (hdf5Files.size() + cacheBatchSize - 1) / cacheBatchSize;
		std::ostringstream loadLabel;
		loadLabel << "load HDF5 batch " << batchesScored << '/' << totalBatches;
		timing.run("Load HDF5", loadLabel.str(),
				   batchEnd - batchBegin, "files", loadBatch);

		const size_t nRec = batchRecords.size();
		totalRecordsScored += nRec;
		if (nRec == 0)
			continue;

		std::vector<MzIndexedRecord> mzIndex = buildMzIndex(batchRecords);
		std::vector<std::vector<CandidateMatch>> scanCandidates(scans.size());
		size_t assignedCandidates = 0;
		auto assignBatch = [&]()
		{
			assignedCandidates = assignCandidatesToScans(scans, batchRecords,
														 mzIndex, mzTol, rtTol,
														 args.isotopeShiftWindow,
														 scanCandidates);
		};
		std::ostringstream assignLabel;
		assignLabel << "assign candidates batch " << batchesScored << '/' << totalBatches;
		timing.run("Assign candidates", assignLabel.str(),
				   nRec, "records", assignBatch);
		totalAssignedCandidates += assignedCandidates;

		auto scoreBatch = [&]()
		{
#pragma omp parallel for schedule(dynamic, 1)
			for (size_t s = 0; s < scans.size(); ++s)
			{
				MS2Scan *scan = scans[s];

				for (const CandidateMatch &candidate : scanCandidates[s])
				{
					const LabeledRecord &lr = batchRecords[candidate.recordIdx];
					const SipRecord &rec = lr.rec;
					const PrecursorMatch &m = candidate.match;

					std::vector<std::vector<double>> vvdYM = rec.vvdYionMass;
					std::vector<std::vector<double>> vvdYP = rec.vvdYionProb;
					std::vector<std::vector<double>> vvdBM = rec.vvdBionMass;
					std::vector<std::vector<double>> vvdBP = rec.vvdBionProb;

					std::string pepForScoring = rec.nakedPeptide;
					int charge = m.matchedCharge;

					double wdp = scan->scoreWeightSumHighMS2(
						&pepForScoring, charge, &vvdYM, &vvdYP, &vvdBM, &vvdBP);

					std::vector<ScoringPsm> &topList = lr.label == -1 ? topDecoy[s] : topTarget[s];
					if (!wouldKeepTopUniquePeptide(topList, rec.nakedPeptide, wdp, args.topPsmsPerScan))
						continue;

					ScoringPsm psm;
					psm.scanIdx = s;
					psm.scanNumber = scan->iScanId;
					psm.label = lr.label;
					psm.wdp = wdp;
					psm.ms2Pct = lr.ms2Pct;
					psm.nakedPeptide = rec.nakedPeptide;
					psm.sipAtom = lr.sipAtom;
					psm.match = m;
					psm.rec = rec;

					pushTopUniquePeptide(topList, std::move(psm), args.topPsmsPerScan);
				}
			}
		};
		std::ostringstream rankLabel;
		rankLabel << "rank WDP batch " << batchesScored << '/' << totalBatches;
		timing.run("Rank WDP candidates", rankLabel.str(),
				   assignedCandidates, "candidates", scoreBatch);
	}
	std::cout << "  Scoring (" << nThreads << " threads, "
			  << batchesScored << " HDF5 batches, "
			  << totalRecordsScored << " records, "
			  << totalAssignedCandidates << " assigned candidates)\n";

	size_t retainedPsms = 0;
	for (const auto &v : topTarget)
		retainedPsms += v.size();
	for (const auto &v : topDecoy)
		retainedPsms += v.size();

	auto scoreRetained = [&](std::vector<std::vector<ScoringPsm>> &perScan)
	{
#pragma omp parallel for schedule(dynamic, 1)
		for (size_t s = 0; s < perScan.size(); ++s)
		{
			const int tid = omp_get_thread_num();
			MS2Scan *scan = scans[s];
			for (ScoringPsm &psm : perScan[s])
			{
				std::vector<std::vector<double>> vvdYM = psm.rec.vvdYionMass;
				std::vector<std::vector<double>> vvdYP = psm.rec.vvdYionProb;
				std::vector<std::vector<double>> vvdBM = psm.rec.vvdBionMass;
				std::vector<std::vector<double>> vvdBP = psm.rec.vvdBionProb;

				keepTopNEnvelopePeaks(vvdYM, vvdYP, args.scoreEnvelopeTopN);
				keepTopNEnvelopePeaks(vvdBM, vvdBP, args.scoreEnvelopeTopN);

				double xcorr = 0.0;
				if (scan->pQuery)
				{
					auto xYM = vvdYM, xYP = vvdYP, xBM = vvdBM, xBP = vvdBP;
					trimIonVectorsForXcorr(scan, xYM, xYP);
					trimIonVectorsForXcorr(scan, xBM, xBP);
					CometSearchMod::ScorePeptidesSIPNoCancelOut(
						xYM, xYP, xBM, xBP, scan,
						tDuplSip[tid], tBinIonSip[tid], tBinSip[tid], xcorr);
				}

				double mvhScore = 0.0;
				if (!scan->bSkip)
				{
					std::string pepForScoring = psm.rec.nakedPeptide;
					MVH::ScoreSequenceVsSpectrumSIP(
						pepForScoring, psm.match.matchedCharge, scan,
						&tMvhIonMasses[tid], vvdYM, vvdYP, vvdBM, vvdBP, mvhScore,
						&tMvhSeqs[tid]);
				}

				std::vector<double> p, q;
				alignSpectra(psm.rec.fragMz, psm.rec.fragExpInt, scan->vdMZ, scan->vdIntensity,
							 fragTol, p, q);
				const double entropy = computeEntropy(p, q);
				alignSpectraSparseCosine(psm.rec.fragMz, psm.rec.fragExpInt,
										 scan->vdMZ, scan->vdIntensity,
										 fragTol, p, q);
				const double cosine = computeCosine(p, q);
				const EnvCounts ec = countMatchedEnvelopes(scan, psm.rec, fragTol);

				LabeledRecord labeledRecord;
				labeledRecord.label = psm.label;
				labeledRecord.ms2Pct = psm.ms2Pct;
				labeledRecord.sipAtom = psm.sipAtom;
				psm.row = makeScoringRow(psm.scanIdx, scan, psm.rec, labeledRecord, psm.match, ec,
										  ms1DataPtr, mzTol, psm.wdp, xcorr, mvhScore, entropy, cosine);
			}
		}
	};
	timing.run("Score retained PSMs", "score retained PSMs",
			   retainedPsms, "PSMs", [&]()
			   {
				   scoreRetained(topTarget);
				   scoreRetained(topDecoy);
			   });

	// -------- Materialize PIN rows --------
	std::vector<MergedRow> kept;
	auto appendPassing = [](std::vector<const ScoringPsm *> &out,
							const std::vector<ScoringPsm> &psms)
	{
		for (const ScoringPsm &psm : psms)
		{
			if (psm.wdp >= kMinPinWdpScore)
				out.push_back(&psm);
		}
	};

	timing.run("Materialize PIN rows", "materialize PIN rows",
			   retainedPsms, "PSMs", [&]()
			   {
				   for (size_t s = 0; s < scans.size(); ++s)
				   {
					   std::vector<const ScoringPsm *> scanPsms;
					   scanPsms.reserve(topTarget[s].size() + topDecoy[s].size());
					   appendPassing(scanPsms, topTarget[s]);
					   appendPassing(scanPsms, topDecoy[s]);
					   if (scanPsms.empty())
						   continue;

					   std::sort(scanPsms.begin(), scanPsms.end(),
								 [](const ScoringPsm *a, const ScoringPsm *b)
								 {
									 if (a->wdp != b->wdp)
										 return a->wdp > b->wdp;
									 if (a->label != b->label)
										 return a->label > b->label;
									 return a->nakedPeptide < b->nakedPeptide;
								 });

					   for (size_t i = 0; i < scanPsms.size(); ++i)
					   {
						   const ScoringPsm &e = *scanPsms[i];
						   MergedRow m;
						   m.label = e.label;
						   m.ms2Pct = e.ms2Pct;
						   m.scanNumber = e.scanNumber;
						   m.rank = static_cast<int>(i + 1);
						   m.row = e.row;
						   kept.push_back(std::move(m));
					   }
				   }
			   });

	fs::path pinPath = fs::path(args.outputDir) / (sampleBasename + ".pin");
	size_t writtenRows = 0;
	timing.run("Write PIN", "write PIN file",
			   kept.size(), "rows", [&]()
			   {
				   const std::vector<PinWriter::SearchSpectraPinRow> pinRows = makeSearchSpectraPinRows(kept);
				   writtenRows = PinWriter::writeSearchSpectraPin(pinPath.string(), sampleBasename, pinRows);
			   });
	std::cout << "Wrote " << pinPath << " (" << writtenRows << " rows, WDP >= "
			  << std::fixed << std::setprecision(1) << kMinPinWdpScore << ")\n";

	// Cleanup
	timing.run("Cleanup", "cleanup scans",
			   scans.size(), "scans", [&]()
			   {
				   for (auto *s : scans)
				   {
					   if (s->pQuery)
					   {
						   delete s->pQuery;
						   s->pQuery = nullptr;
					   }
					   delete s;
				   }
			   });

	const double processWallSeconds = omp_get_wtime() - processWallStart;
	const double processCpuSeconds = processTreeCpuSeconds() - processCpuStart;
	timing.printSummary(processWallSeconds, processCpuSeconds,
						totalRecordsScored, totalAssignedCandidates,
						retainedPsms, writtenRows);
	return hdf5LoadFailures == 0 ? 0 : 1;
}

} // namespace

int SearchSpectraWorkflow::run(int argc, char **argv)
{
	Args args = parseArgs(argc, argv);

	std::vector<std::string> scanFiles;
	if (!args.singleHdf5.empty())
		scanFiles.push_back(args.singleHdf5);
	else
		scanFiles = listFiles(args.workingDir, {".h5", ".H5", ".hdf5", ".HDF5"});
	if (scanFiles.empty())
	{
		std::cerr << "No Raxport HDF5 scan files found.\n";
		return 1;
	}

	std::vector<std::string> hdf5Files = listFiles(args.hdf5Dir, {".h5", ".H5", ".hdf5", ".HDF5"});
	if (hdf5Files.empty())
	{
		std::cerr << "No HDF5 spectra library files in " << args.hdf5Dir << "\n";
		return 1;
	}
	// target first, decoy second (no semantic effect — every shard knows its
	// own label — but predictable ordering is convenient for debugging).
	std::stable_partition(hdf5Files.begin(), hdf5Files.end(),
						  [](const std::string &p)
						  { return !isDecoyHdf5(p); });

	if (!fs::is_directory(args.outputDir))
		fs::create_directories(args.outputDir);

	const double t0 = omp_get_wtime();
	int rc = 0;
	for (const auto &scanPath : scanFiles)
	{
		int r = processOneHdf5(args, scanPath, hdf5Files);
		if (r != 0)
			rc = r;
	}
	std::cout << "sipros search-spectra finished in " << (omp_get_wtime() - t0) << "s\n";
	return rc;
}
