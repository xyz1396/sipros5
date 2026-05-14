// sipros_search_spectra
//
// Re-score real experimental FT2 scans against a pre-generated SIP spectra
// library (HDF5). HDF5 records are loaded into bounded in-memory batches.
// Output: one Percolator PIN per FT2 file, 28 columns.

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
#include "averagine.h"

namespace fs = std::filesystem;

namespace
{

// -------------------- Args --------------------

struct Args
{
	std::string workingDir;
	std::string singleFt2;
	std::string configFile;
	std::string configDir;
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
	int topPsmsPerScan = 10;
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
		<< " -w <FT2 dir> [-f <single.FT2>] -c <config.cfg> | -g <config dir>\n"
		<< "    -h5 <SIP spectra dir> -o <PIN output dir> [-t <N>] [--rt-tolerance <min>]\n"
		<< "    [--tolerance-ms1 <N>] [--tolerance-ms1-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance-ms2 <N>] [--tolerance-ms2-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance <N>] [--tolerance-unit ppm|da]            shortcut: set BOTH MS1 and MS2\n"
		<< "    [--isotope-shift-window <N>]                           precursor isotope shifts +/-N (default: 3)\n"
		<< "    [--score-envelope-top-n <N>]                           Xcorr/MVH envelope peaks (default: 2)\n"
		<< "    [--top-psms-per-scan <N>]                              WDP winners per scan/label (default: 10)\n";
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
			a.singleFt2 = next();
		else if (k == "-c")
			a.configFile = next();
		else if (k == "-g")
			a.configDir = next();
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
	if (a.singleFt2.empty() && a.workingDir.empty())
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

// -------------------- FT2 parser (self-contained) --------------------
//
// Mirrors the subset of MS2ScanVector::ReadFT2File / saveFT2Scan needed for
// SIP search: parse S/Z/I/D/peak lines and build MS2Scan objects with the
// public fields populated. Avoids touching MS2ScanVector's private state.

const char *skipWs(const char *p)
{
	while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
		++p;
	return p;
}

const char *skipToken(const char *p)
{
	p = skipWs(p);
	while (*p && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n')
		++p;
	return p;
}

bool nextToken(const char *&p, const char *&begin, const char *&end)
{
	begin = skipWs(p);
	if (!*begin)
	{
		p = begin;
		return false;
	}
	end = begin;
	while (*end && *end != ' ' && *end != '\t' && *end != '\r' && *end != '\n')
		++end;
	p = end;
	return true;
}

bool tokenEquals(const char *begin, const char *end, const char *literal)
{
	const size_t len = static_cast<size_t>(end - begin);
	return std::strlen(literal) == len && std::strncmp(begin, literal, len) == 0;
}

bool parseDoubleToken(const char *&p, double &value)
{
	p = skipWs(p);
	if (!*p)
		return false;
	char *end = nullptr;
	value = std::strtod(p, &end);
	if (end == p)
		return false;
	p = end;
	return true;
}

bool parseIntToken(const char *&p, int &value)
{
	p = skipWs(p);
	if (!*p)
		return false;
	char *end = nullptr;
	long parsed = std::strtol(p, &end, 10);
	if (end == p)
		return false;
	value = static_cast<int>(parsed);
	p = end;
	return true;
}

void reserveScanPeakStorage(MS2Scan *scan)
{
	constexpr size_t kInitialPeakReserve = 512;
	scan->vdMZ.reserve(kInitialPeakReserve);
	scan->vdIntensity.reserve(kInitialPeakReserve);
	scan->viCharge.reserve(kInitialPeakReserve);
}

bool parsePeakLine(const std::string &line,
				   double &mz,
				   double &intensity,
				   int &charge,
				   bool &highRes)
{
	const char *p = line.c_str();
	if (!parseDoubleToken(p, mz))
		return false;
	if (!parseDoubleToken(p, intensity))
		return false;

	highRes = false;
	charge = 0;
	double ignored = 0.0;
	for (int token = 3; token <= 5; ++token)
	{
		if (!parseDoubleToken(p, ignored))
			return true;
	}
	if (parseIntToken(p, charge))
		highRes = true;
	return true;
}

void finalizeFt2Scan(MS2Scan *scan, std::vector<MS2Scan *> &all)
{
	if (!scan)
		return;
	if (scan->vdIntensity.empty())
	{
		delete scan;
		return;
	}
	scan->isMS2HighRes = (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1);
	if (scan->iParentChargeState == 0)
	{
		if (!scan->iParentChargeStates.empty())
		{
			// DIA: use max charged mass across the list
			double maxMass = 0.0;
			int maxCharge = 0;
			size_t k = std::min(scan->dParentMZs.size(), scan->iParentChargeStates.size());
			for (size_t i = 0; i < k; ++i)
			{
				int c = scan->iParentChargeStates[i];
				double mz = scan->dParentMZs[i];
				double m = mz * c;
				if (m > maxMass)
					maxMass = m;
				if (c > maxCharge)
					maxCharge = c;
			}
			scan->dParentMass = maxMass;
			scan->dParentNeutralMass = maxMass - maxCharge * ProNovoConfig::getProtonMass();
			if (maxMass > ProNovoConfig::dMaxMS2ScanMass)
				ProNovoConfig::dMaxMS2ScanMass = maxMass;
			if (maxCharge > ProNovoConfig::iMaxPercusorCharge)
				ProNovoConfig::iMaxPercusorCharge = maxCharge;
			all.push_back(scan);
		}
		else
			delete scan;
	}
	else
	{
		double charged = scan->dParentMZ * scan->iParentChargeState;
		scan->dParentMass = charged;
		scan->dParentNeutralMass = charged - scan->iParentChargeState * ProNovoConfig::getProtonMass();
		if (charged > ProNovoConfig::dMaxMS2ScanMass)
			ProNovoConfig::dMaxMS2ScanMass = charged;
		if (scan->iParentChargeState > ProNovoConfig::iMaxPercusorCharge)
			ProNovoConfig::iMaxPercusorCharge = scan->iParentChargeState;
		all.push_back(scan);
	}
}

bool readFt2File(const std::string &path, std::vector<MS2Scan *> &out)
{
	std::ifstream in(path);
	if (!in)
		return false;
	MS2Scan *cur = nullptr;
	bool first = true;
	std::string line;
	line.reserve(256);
	while (std::getline(in, line))
	{
		if (line.empty())
			continue;
		const char *lineStart = skipWs(line.c_str());
		if (!*lineStart)
			continue;
		char c = *lineStart;
		if (c >= '0' && c <= '9')
		{
			if (!cur)
				continue;
			double mz = 0.0;
			double inten = 0.0;
			int charge = 0;
			bool highRes = false;
			if (!parsePeakLine(line, mz, inten, charge, highRes))
				continue;
			cur->isMS2HighRes = highRes;
			if (mz > ProNovoConfig::maxObservedMz)
				ProNovoConfig::maxObservedMz = mz;
			if (mz < ProNovoConfig::minObservedMz)
				ProNovoConfig::minObservedMz = mz;
			cur->vdMZ.push_back(mz);
			cur->vdIntensity.push_back(inten);
			cur->viCharge.push_back(charge);
		}
		else if (c == 'S')
		{
			if (!first)
				finalizeFt2Scan(cur, out);
			first = false;
			cur = new MS2Scan;
			cur->sFT2Filename = path;
			cur->iParentChargeState = 0;
			cur->dParentNeutralMass = 0;
			reserveScanPeakStorage(cur);
			const char *p = lineStart;
			p = skipToken(p); // S
			int scanId = 0;
			if (parseIntToken(p, scanId))
			{
				const char *mzBegin = nullptr;
				const char *mzEnd = nullptr;
				const char *mzToken = p;
				double parentMz = 0.0;
				if (nextToken(mzToken, mzBegin, mzEnd))
				{
					if (parseDoubleToken(p, parentMz))
					{
						cur->iScanId = scanId;
						cur->dParentMZ = parentMz;
						const size_t mzLen = static_cast<size_t>(mzEnd - mzBegin);
						cur->isMS1HighRes = (mzLen >= 2 && mzBegin[mzLen - 2] != '.');
					}
				}
			}
		}
		else if (c == 'Z' && cur)
		{
			const char *p = lineStart;
			p = skipToken(p); // Z
			int charge = 0;
			if (parseIntToken(p, charge))
				cur->iParentChargeState = charge;

			p = skipToken(p); // legacy mass token before DIA charge/mz pairs
			while (*skipWs(p))
			{
				int diaCharge = 0;
				double diaMz = 0.0;
				if (!parseIntToken(p, diaCharge) || !parseDoubleToken(p, diaMz))
					break;
				cur->iParentChargeStates.push_back(diaCharge);
				cur->dParentMZs.push_back(diaMz);
			}
		}
		else if (c == 'I' && cur)
		{
			const char *p = lineStart;
			const char *keyBegin = nullptr;
			const char *keyEnd = nullptr;
			const char *valueBegin = nullptr;
			const char *valueEnd = nullptr;
			p = skipToken(p); // I
			if (nextToken(p, keyBegin, keyEnd) && nextToken(p, valueBegin, valueEnd))
			{
				std::string value(valueBegin, static_cast<size_t>(valueEnd - valueBegin));
				if (tokenEquals(keyBegin, keyEnd, "ScanType"))
					cur->setScanType(value);
				else if (tokenEquals(keyBegin, keyEnd, "RetentionTime") || tokenEquals(keyBegin, keyEnd, "RTime"))
					cur->setRTime(value);
			}
		}
		else if (c == 'D' && cur)
		{
			const char *p = lineStart;
			const char *keyBegin = nullptr;
			const char *keyEnd = nullptr;
			p = skipToken(p); // D
			if (nextToken(p, keyBegin, keyEnd) && tokenEquals(keyBegin, keyEnd, "ParentScanNumber"))
			{
				int parentScanId = 0;
				if (parseIntToken(p, parentScanId))
					cur->iParentScanID = parentScanId;
			}
		}
	}
	if (!first)
		finalizeFt2Scan(cur, out);
	return true;
}

// -------------------- FT1 parser and MS1 isotope helpers --------------------

struct Ft1Scan
{
	int scanNumber = 0;
	double retentionTime = 0.0;
	std::vector<double> mz;
	std::vector<double> intensity;
	std::vector<int> charge;
};

struct Ft1Data
{
	std::vector<Ft1Scan> scans;
	std::unordered_map<int, size_t> scanNumberToIndex;
};

void finalizeFt1Scan(Ft1Scan &scan, Ft1Data &data)
{
	if (scan.scanNumber <= 0 || scan.mz.empty())
		return;
	std::vector<size_t> order(scan.mz.size());
	for (size_t i = 0; i < order.size(); ++i)
		order[i] = i;
	std::stable_sort(order.begin(), order.end(),
					 [&](size_t a, size_t b)
					 { return scan.mz[a] < scan.mz[b]; });

	Ft1Scan sorted;
	sorted.scanNumber = scan.scanNumber;
	sorted.retentionTime = scan.retentionTime;
	sorted.mz.reserve(order.size());
	sorted.intensity.reserve(order.size());
	sorted.charge.reserve(order.size());
	for (size_t idx : order)
	{
		sorted.mz.push_back(scan.mz[idx]);
		sorted.intensity.push_back(scan.intensity[idx]);
		sorted.charge.push_back(idx < scan.charge.size() ? scan.charge[idx] : 0);
	}

	data.scanNumberToIndex[sorted.scanNumber] = data.scans.size();
	data.scans.push_back(std::move(sorted));
}

bool readFt1File(const std::string &path, Ft1Data &data)
{
	std::ifstream in(path);
	if (!in)
		return false;
	data.scans.clear();
	data.scanNumberToIndex.clear();

	Ft1Scan cur;
	bool haveScan = false;
	std::string line;
	line.reserve(256);
	while (std::getline(in, line))
	{
		if (line.empty())
			continue;
		const char *lineStart = skipWs(line.c_str());
		if (!*lineStart)
			continue;
		const char c = *lineStart;
		if (c >= '0' && c <= '9')
		{
			if (!haveScan)
				continue;
			double mz = 0.0;
			double inten = 0.0;
			int charge = 0;
			bool highRes = false;
			if (!parsePeakLine(line, mz, inten, charge, highRes))
				continue;
			cur.mz.push_back(mz);
			cur.intensity.push_back(inten);
			cur.charge.push_back(charge);
		}
		else if (c == 'S')
		{
			if (haveScan)
				finalizeFt1Scan(cur, data);
			cur = Ft1Scan();
			haveScan = true;
			const char *p = lineStart;
			p = skipToken(p); // S
			int scanId = 0;
			if (parseIntToken(p, scanId))
				cur.scanNumber = scanId;
		}
		else if (c == 'I' && haveScan)
		{
			const char *p = lineStart;
			const char *keyBegin = nullptr;
			const char *keyEnd = nullptr;
			p = skipToken(p); // I
			if (nextToken(p, keyBegin, keyEnd) &&
				(tokenEquals(keyBegin, keyEnd, "RetentionTime") || tokenEquals(keyBegin, keyEnd, "RTime")))
			{
				double rt = 0.0;
				if (parseDoubleToken(p, rt))
					cur.retentionTime = rt;
			}
		}
	}
	if (haveScan)
		finalizeFt1Scan(cur, data);
	return true;
}

std::string ft1PathForFt2(const std::string &ft2Path)
{
	fs::path p(ft2Path);
	p.replace_extension(".FT1");
	if (fs::exists(p))
		return p.string();
	p.replace_extension(".ft1");
	if (fs::exists(p))
		return p.string();
	return "";
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
	double targetSipAbundancePct = 0.0;
	std::string sipAtom = "C";
	int sipIsotopeMassNumber = 13;
};

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

// Remove only terminus brackets so averagine can still see PTM symbols.
std::string peptideBodyWithPtms(const std::string &decorated)
{
	std::string out;
	out.reserve(decorated.size());
	for (char c : decorated)
	{
		if (c != '[' && c != ']')
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

void ensureDefaultNTermAcetylation()
{
	Isotopologue &iso = ProNovoConfig::configIsotopologue;
	if (iso.mResidueAtomicComposition.find("%") == iso.mResidueAtomicComposition.end())
	{
		// Element_List is CHONPS in Sipros configs; % is acetylation C2H2O.
		iso.mResidueAtomicComposition["%"] = {2, 2, 1, 0, 0, 0};
	}

	IsotopeDistribution dist;
	iso.computeIsotopicDistribution(iso.mResidueAtomicComposition["%"], dist);
	iso.vResidueIsotopicDistribution["%"] = dist;

	if (std::find(ProNovoConfig::vsSingleResidueNames.begin(),
				  ProNovoConfig::vsSingleResidueNames.end(),
				  "%") == ProNovoConfig::vsSingleResidueNames.end())
	{
		ProNovoConfig::vsSingleResidueNames.push_back("%");
		ProNovoConfig::vdSingleResidueMasses.push_back(dist.getMostAbundantMass());
	}
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

std::string readStringAttribute(const H5::H5Object &obj, const char *name)
{
	if (!obj.attrExists(name))
		return "";
	H5::Attribute a = obj.openAttribute(name);
	H5::StrType type = a.getStrType();
	if (type.isVariableStr())
	{
		char *p = nullptr;
		a.read(type, &p);
		std::string s = p ? p : "";
		if (p)
			std::free(p);
		return s;
	}
	std::vector<char> buf(type.getSize() + 1, '\0');
	a.read(type, buf.data());
	return std::string(buf.data());
}

double readDoubleAttribute(const H5::H5Object &obj, const char *name, double dflt = 0.0)
{
	if (!obj.attrExists(name))
		return dflt;
	double x = dflt;
	H5::Attribute a = obj.openAttribute(name);
	a.read(H5::PredType::NATIVE_DOUBLE, &x);
	return x;
}

int readIntAttribute(const H5::H5Object &obj, const char *name, int dflt = 0)
{
	if (!obj.attrExists(name))
		return dflt;
	int x = dflt;
	H5::Attribute a = obj.openAttribute(name);
	a.read(H5::PredType::NATIVE_INT, &x);
	return x;
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
		meta.targetSipAbundancePct = readDoubleAttribute(f, "target_sip_abundance_pct", 0.0);
		meta.sipAtom = readStringAttribute(f, "sip_atom");
		if (meta.sipAtom.empty())
			meta.sipAtom = "C";
		meta.sipIsotopeMassNumber = readIntAttribute(f, "sip_isotope_mass_number", 13);

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

int sipAtomIndex(const std::string &sipAtom)
{
	if (sipAtom.empty())
		return 0;
	const char atom = static_cast<char>(std::toupper(static_cast<unsigned char>(sipAtom.front())));
	if (atom == 'H')
		return 1;
	if (atom == 'O')
		return 2;
	if (atom == 'N')
		return 3;
	if (atom == 'S')
		return 5;
	return 0;
}

int sipNominalShiftPerAtom(const std::string &sipAtom)
{
	if (sipAtom.empty())
		return 1;
	const char atom = static_cast<char>(std::toupper(static_cast<unsigned char>(sipAtom.front())));
	return (atom == 'O' || atom == 'S') ? 2 : 1;
}

double expectedNaturalNominalShiftExceptTarget(const std::array<int, 6> &atomCounts,
											   int targetAtomIndex,
											   int targetIsotopeIndex)
{
	double expectedShift = 0.0;
	const auto &atomDistributions = ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution;
	const size_t atomCount = std::min(atomCounts.size(), atomDistributions.size());
	for (size_t atomIdx = 0; atomIdx < atomCount; ++atomIdx)
	{
		if (atomCounts[atomIdx] <= 0)
			continue;
		const auto &probs = atomDistributions[atomIdx].vProb;
		for (size_t isotopeIdx = 1; isotopeIdx < probs.size(); ++isotopeIdx)
		{
			if (static_cast<int>(atomIdx) == targetAtomIndex &&
				static_cast<int>(isotopeIdx) == targetIsotopeIndex)
			{
				continue;
			}
			expectedShift += static_cast<double>(atomCounts[atomIdx]) *
							 static_cast<double>(isotopeIdx) * probs[isotopeIdx];
		}
	}
	return expectedShift;
}

int countMissCleavage(const std::string &naked)
{
	if (naked.size() <= 1)
		return 0;
	int n = 0;
	for (size_t i = 0; i + 1 < naked.size(); ++i)
		if (naked[i] == 'K' || naked[i] == 'R')
			++n;
	return n;
}

int countPTM(const std::string &decorated)
{
	int n = 0;
	for (char c : decorated)
		if (!std::isalpha(static_cast<unsigned char>(c)) && c != '[' && c != ']')
			++n;
	return n;
}

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

struct IsotopicPeak
{
	double mz = 0.0;
	int charge = 0;
	double intensity = 0.0;
	int isotopeIndex = -1;
};

int ft1PeakCharge(const Ft1Scan &scan, size_t idx)
{
	if (idx < scan.charge.size())
		return scan.charge[idx];
	return 0;
}

size_t findFt1Peak(const Ft1Scan &scan,
				   double targetMz,
				   const MassTolerance &mzTol,
				   int requiredCharge = -1)
{
	if (scan.mz.empty())
		return std::numeric_limits<size_t>::max();

	size_t best = std::numeric_limits<size_t>::max();
	double bestIntensity = 0.0;
	const double mzTolerance = mzTol.daAt(targetMz);
	auto first = std::lower_bound(scan.mz.begin(), scan.mz.end(), targetMz - mzTolerance);
	for (auto it = first; it != scan.mz.end() && *it <= targetMz + mzTolerance; ++it)
	{
		const size_t idx = static_cast<size_t>(it - scan.mz.begin());
		if (requiredCharge >= 0 && ft1PeakCharge(scan, idx) != requiredCharge)
		{
			continue;
		}
		if (std::fabs(scan.mz[idx] - targetMz) <= mzTolerance &&
			scan.intensity[idx] > bestIntensity)
		{
			bestIntensity = scan.intensity[idx];
			best = idx;
		}
	}
	return best;
}

std::vector<IsotopicPeak> findFt1IsotopicPeaks(const Ft1Data *ft1Data,
											   int &ms1ScanNumber,
											   int precursorCharge,
											   double monoPrecursorMz,
											   double matchedPrecursorMz,
											   int maxNominalShift,
											   const MassTolerance &mzTol)
{
	std::vector<IsotopicPeak> peaks;
	if (!ft1Data || precursorCharge <= 0)
		return peaks;
	auto scanIt = ft1Data->scanNumberToIndex.find(ms1ScanNumber);
	if (scanIt == ft1Data->scanNumberToIndex.end())
		return peaks;

	const double neutronMz = ProNovoConfig::getNeutronMass() / precursorCharge;
	constexpr int kMs1QuantWindowNeutrons = 10;
	constexpr int kMs1AssignmentIndexTolerance = 1;
	const int assignedIndex = static_cast<int>(std::round(
		(matchedPrecursorMz - monoPrecursorMz) / neutronMz));
	if (assignedIndex < 0 || assignedIndex > maxNominalShift)
		return peaks;

	size_t scanIdx = scanIt->second;
	for (int attempt = 0; attempt < 3; ++attempt)
	{
		const Ft1Scan &scan = ft1Data->scans[scanIdx];
		peaks.clear();
		for (int isotopeIndex = 0; isotopeIndex <= maxNominalShift; ++isotopeIndex)
		{
			const double expectedMz = monoPrecursorMz + isotopeIndex * neutronMz;
			if (std::fabs(expectedMz - matchedPrecursorMz) >
				kMs1QuantWindowNeutrons * neutronMz)
			{
				continue;
			}

			size_t idx = findFt1Peak(scan, expectedMz, mzTol, precursorCharge);
			if (idx == std::numeric_limits<size_t>::max())
			{
				idx = findFt1Peak(scan, expectedMz, mzTol);
			}
			if (idx != std::numeric_limits<size_t>::max())
			{
				peaks.push_back({scan.mz[idx],
								 ft1PeakCharge(scan, idx),
								 scan.intensity[idx],
								 isotopeIndex});
			}
		}
		bool hasAssignmentAnchor = false;
		bool hasExactAssignmentAnchor = false;
		for (const IsotopicPeak &peak : peaks)
		{
			const int indexDelta = std::abs(peak.isotopeIndex - assignedIndex);
			if (indexDelta <= kMs1AssignmentIndexTolerance)
				hasAssignmentAnchor = true;
			if (indexDelta == 0)
				hasExactAssignmentAnchor = true;
		}
		if (!peaks.empty() && hasAssignmentAnchor &&
			(peaks.size() > 1 || hasExactAssignmentAnchor))
		{
			ms1ScanNumber = scan.scanNumber;
			break;
		}
		peaks.clear();
		if (scanIdx == 0)
			break;
		--scanIdx;
	}
	if (peaks.empty())
		return peaks;

	std::sort(peaks.begin(), peaks.end(),
			  [](const IsotopicPeak &a, const IsotopicPeak &b)
			  { return a.mz < b.mz; });
	peaks.erase(std::unique(peaks.begin(), peaks.end(),
							[](const IsotopicPeak &a, const IsotopicPeak &b)
							{
								return std::fabs(a.mz - b.mz) <= 1e-12;
							}),
				peaks.end());
	return peaks;
}

std::vector<IsotopicPeak> ft1IsotopicPeaksForMatch(const Ft1Data *ft1Data,
												   const MS2Scan *scan,
												   const SipRecord &rec,
												   const PrecursorMatch &match,
												   const MassTolerance &mzTol,
												   const std::string &sipAtom,
												   double &baseMass)
{
	baseMass = 0.0;
	if (!ft1Data || !scan || scan->iParentScanID <= 0)
		return {};

	averagine avg;
	baseMass = avg.calPrecursorBaseMass(peptideBodyWithPtms(rec.peptide));
	int ms1ScanNumber = scan->iParentScanID;
	const int precursorCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
	avg.calPepAtomCounts(peptideBodyWithPtms(rec.peptide));
	const int atomIndex = sipAtomIndex(sipAtom);
	const int targetNominalShift = sipNominalShiftPerAtom(sipAtom);
	const double maxTargetShift = atomIndex < static_cast<int>(avg.pepAtomCounts.size())
									  ? avg.pepAtomCounts[atomIndex] * targetNominalShift
									  : 0.0;
	const int maxNominalShift = std::min(512, std::max(20, static_cast<int>(std::ceil(maxTargetShift)) + 20));
	const double monoPrecursorMz = baseMass / precursorCharge + ProNovoConfig::getProtonMass();
	return findFt1IsotopicPeaks(
		ft1Data, ms1ScanNumber, precursorCharge, monoPrecursorMz,
		match.matchedMz, maxNominalShift, mzTol);
}

struct Ms1AbundanceResult
{
	double abundancePct = 0.0;
	int isotopicPeakCount = 0;
};

Ms1AbundanceResult ms1AbundanceFromFt1Peaks(const std::vector<IsotopicPeak> &peaks,
											double baseMass,
											const SipRecord &rec,
											const PrecursorMatch &match,
											const std::string &sipAtom)
{
	if (peaks.empty())
	{
		return {};
	}
	const int precursorCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
	if (precursorCharge <= 0)
	{
		return {};
	}

	averagine avg;
	const char atom = sipAtom.empty()
						  ? 'C'
						  : static_cast<char>(std::toupper(static_cast<unsigned char>(sipAtom.front())));
	const int targetIsotopeIndex = sipNominalShiftPerAtom(std::string(1, atom));
	avg.calPepAtomCounts(peptideBodyWithPtms(rec.peptide));
	const int atomIndex = sipAtomIndex(sipAtom);
	const double atomNumber = avg.pepAtomCounts[atomIndex];
	if (atomNumber <= 0.0)
	{
		return {};
	}

	const double maxIsotopeIndex = atomNumber * targetIsotopeIndex;
	const double baseMz = baseMass / precursorCharge + ProNovoConfig::getProtonMass();

	double sumIntensity = 0.0;
	double weightedIsotopeIndex = 0.0;
	int validPeakCount = 0;
	for (size_t i = 0; i < peaks.size(); ++i)
	{
		const IsotopicPeak &peak = peaks[i];
		if (peak.intensity <= 0.0)
		{
			continue;
		}
		const int isotopeIndex = peak.isotopeIndex >= 0
									 ? peak.isotopeIndex
									 : static_cast<int>(std::round(
										   (peak.mz - baseMz) /
										   ProNovoConfig::getNeutronMass() * precursorCharge));
		if (isotopeIndex < 0)
		{
			continue;
		}
		sumIntensity += peak.intensity;
		weightedIsotopeIndex += peak.intensity * std::min(static_cast<double>(isotopeIndex), maxIsotopeIndex);
		++validPeakCount;
	}
	if (validPeakCount == 0 || sumIntensity <= 0.0)
	{
		return {};
	}

	const double meanIsotopeIndex = weightedIsotopeIndex / sumIntensity;
	const double naturalOtherShift = expectedNaturalNominalShiftExceptTarget(
		avg.pepAtomCounts, atomIndex, targetIsotopeIndex);
	double pct = (meanIsotopeIndex - naturalOtherShift) / maxIsotopeIndex * 100.0;
	if (!std::isfinite(pct))
	{
		return {};
	}
	pct = std::min(100.0, std::max(0.0, pct));
	return {pct, validPeakCount};
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

// -------------------- Parent: PIN writer --------------------

const char *kPinHeader =
	"SpecId\tLabel\tScanNr\tExpMass\tretentiontime\tranks\tparentCharges\t"
	"massErrors\tisotopicMassWindowShifts\tmzShiftFromisolationWindowCenters\t"
	"peptideLengths\tmissCleavageSiteNumbers\tPTMnumbers\tisotopicPeakNumbers\t"
	"MS1IsotopicAbundances\tMS2IsotopicAbundances\tisotopicAbundanceDiffs\t"
	"WDPscores\tXcorrScores\tMVHscores\tentropyScores\tcosineScores\t"
	"matchedYenvelopes\tmatchedBenvelopes\tdeltaRT\t"
	"diffScores\tlog10_precursorIntensities\tPeptide\tProteins\n";

constexpr double kMinPinWdpScore = 0.5;

struct MergedRow
{
	ShardPsmRow row;
	int32_t label = 0;
	double ms2Pct = 0.0;
	int scanNumber = 0; // real MS2 scan ID
	int rank = 0;
};

size_t writePinFileWithDiff(const std::string &path,
							const std::string &ft2Basename,
							const std::vector<MergedRow> &rows)
{
	std::vector<const MergedRow *> outputRows;
	outputRows.reserve(rows.size());
	for (const auto &m : rows)
		outputRows.push_back(&m);

	// Target and decoy rows are ranked together per scan, so diffScores use
	// the best WDP among both labels for that scan.
	std::map<int, double> topByScan;
	for (const MergedRow *m : outputRows)
	{
		auto it = topByScan.find(m->scanNumber);
		if (it == topByScan.end() || m->row.wdp > it->second)
			topByScan[m->scanNumber] = m->row.wdp;
	}

	std::ostringstream pin;
	pin.imbue(std::locale::classic());
	pin << kPinHeader;
	pin << std::fixed << std::setprecision(6);
	for (const MergedRow *m : outputRows)
	{
		const auto &r = m->row;
		double top = topByScan[m->scanNumber];
		double diff = top - r.wdp;
		pin << ft2Basename << "." << m->scanNumber << "." << m->rank << '\t';
		pin << m->label << '\t';
		pin << m->scanNumber << '\t';
		pin << r.calcMassNeutral << '\t';
		pin << r.rtScan << '\t';
		pin << m->rank << '\t';
		pin << r.parentCharge << '\t';
		pin << r.mzAbsErrPpm << '\t';
		pin << r.isotopicShift << '\t';
		pin << r.mzShiftDa << '\t';
		pin << r.peptideLength << '\t';
		pin << r.missCleavage << '\t';
		pin << r.ptmCount << '\t';
		pin << r.isotopicPeakNumbers << '\t';
		pin << r.ms1IsotopicAbundancePct << '\t';
		pin << m->ms2Pct << '\t';
		pin << (r.ms1IsotopicAbundancePct - m->ms2Pct) << '\t';
		pin << r.wdp << '\t';
		pin << r.xcorr << '\t';
		pin << r.mvh << '\t';
		pin << r.entropy << '\t';
		pin << r.cosine << '\t';
		pin << r.matchedY << '\t';
		pin << r.matchedB << '\t';
		pin << r.deltaRT << '\t';
		pin << diff << '\t';
		pin << r.log10PrecursorIntensity << '\t';
		pin << r.peptide << '\t';
		pin << r.proteins << '\n';
	}

	const std::string pinText = pin.str();
	std::ofstream os(path);
	if (!os)
	{
		std::cerr << "Cannot write " << path << "\n";
		return 0;
	}
	os.write(pinText.data(), static_cast<std::streamsize>(pinText.size()));
	if (!os)
	{
		std::cerr << "Failed while writing " << path << "\n";
		return 0;
	}
	return outputRows.size();
}

// -------------------- Parent: per-FT2 driver --------------------

// Holds one in-memory record paired with its source HDF5's label/metadata.
struct LabeledRecord
{
	SipRecord rec;
	int label = 0; // +1 / -1
	double ms2Pct = 0.0;
	std::string sipAtom = "C";
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
	if (!r.pod(meta.targetSipAbundancePct) ||
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
		err = "cannot parse flat HDF5 shared-memory blob";
	return ok;
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
	const std::string sipAtom = meta.sipAtom.empty() ? "C" : std::string(1, meta.sipAtom[0]);
	for (auto &record : records)
	{
		LabeledRecord lr;
		lr.rec = std::move(record);
		lr.label = label;
		lr.ms2Pct = meta.targetSipAbundancePct;
		lr.sipAtom = sipAtom;
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
						   const Ft1Data *ft1Data,
						   const MassTolerance &mzTol,
						   double wdp,
						   double xcorr,
						   double mvh,
						   double entropy,
						   double cosine)
{
	double baseMass = 0.0;
	const std::vector<IsotopicPeak> ft1Peaks = ft1IsotopicPeaksForMatch(
		ft1Data, scan, rec, match, mzTol, labeledRecord.sipAtom, baseMass);
	const Ms1AbundanceResult ms1Abundance = ms1AbundanceFromFt1Peaks(
		ft1Peaks, baseMass, rec, match, labeledRecord.sipAtom);

	ShardPsmRow row;
	row.scanIdx = static_cast<int32_t>(scanIdx);
	row.parentCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
	row.isotopicShift = match.isotopicShift;
	row.matchedY = envCounts.matchedY;
	row.matchedB = envCounts.matchedB;
	row.peptideLength = static_cast<int>(rec.nakedPeptide.size());
	row.missCleavage = countMissCleavage(rec.nakedPeptide);
	row.ptmCount = countPTM(rec.peptide);
	row.isotopicPeakNumbers = ms1Abundance.isotopicPeakCount;
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

int processOneFt2(const Args &args, const std::string &ft2Path,
				  const std::vector<std::string> &hdf5Files,
				  const Ft1Data *ft1Data)
{
	const double processWallStart = omp_get_wtime();
	const double processCpuStart = processTreeCpuSeconds();
	TimingLogger timing;

	const std::string ft2Basename = fs::path(ft2Path).stem().string();

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
	const double fragTolDaCfg = args.toleranceMs2Ppm
									? args.toleranceMs2 * fragTolRefMz / 1.0e6
									: args.toleranceMs2;

	if (!ProNovoConfig::setFilename(args.configFile))
	{
		std::cerr << "Cannot load config " << args.configFile << "\n";
		return 2;
	}
	ensureDefaultNTermAcetylation();
	ProNovoConfig::setMassAccuracy(parentTolDa, fragTolDaCfg);

	std::cout << "  Tolerance: MS1=" << args.toleranceMs1
			  << (args.toleranceMs1Ppm ? " ppm" : " Da")
			  << ", MS2=" << args.toleranceMs2
			  << (args.toleranceMs2Ppm ? " ppm" : " Da")
			  << "  (scoring config: parent=" << parentTolDa
			  << " Da, fragment=" << fragTolDaCfg << " Da)\n";
	std::cout << "  Xcorr/MVH envelope top-N peaks: "
			  << args.scoreEnvelopeTopN << "\n";
	std::cout << "  WDP top PSMs per scan/label: "
			  << args.topPsmsPerScan << "\n";
	std::cout << "  Precursor isotope shift window: +/-"
			  << args.isotopeShiftWindow << "\n";

	// -------- Read FT2 and score assigned SIP spectra --------
	const int nThreads = std::max(1, args.threads);
	omp_set_num_threads(nThreads);
	timing.printHeader();

	std::vector<MS2Scan *> scans;
	bool readFt2Ok = false;
	timing.run("Read FT2", "read FT2 scans", 0, "", [&]()
			   { readFt2Ok = readFt2File(ft2Path, scans); });
	if (!readFt2Ok)
	{
		std::cerr << "Parent: cannot read " << ft2Path << "\n";
		return 3;
	}
	std::cout << "FT2: " << ft2Path << "  (" << scans.size() << " scans)\n";
	std::cout << "HDF5: " << (fs::path(args.hdf5Dir) / "*.h5").string()
			  << " (" << hdf5Files.size() << " percentages)\n";
	auto preprocessScans = [&]()
	{
#pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < scans.size(); ++i)
			scans[i]->preprocess();
	};
	timing.run("Preprocess FT2", "preprocess FT2 scans",
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
										  ft1Data, mzTol, psm.wdp, xcorr, mvhScore, entropy, cosine);
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

	fs::path pinPath = fs::path(args.outputDir) / (ft2Basename + ".pin");
	size_t writtenRows = 0;
	timing.run("Write PIN", "write PIN file",
			   kept.size(), "rows", [&]()
			   {
				   writtenRows = writePinFileWithDiff(pinPath.string(), ft2Basename, kept);
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

int main(int argc, char **argv)
{
	Args args = parseArgs(argc, argv);

	if (args.configFile.empty() && !args.configDir.empty())
	{
		auto cfgs = listFiles(args.configDir, {".cfg", ".CFG"});
		if (cfgs.empty())
		{
			std::cerr << "No .cfg in " << args.configDir << "\n";
			return 1;
		}
		args.configFile = cfgs.front();
	}
	if (args.configFile.empty())
	{
		std::cerr << "-c <config> or -g <config dir> required\n";
		return 1;
	}

	std::vector<std::string> ft2Files;
	if (!args.singleFt2.empty())
		ft2Files.push_back(args.singleFt2);
	else
		ft2Files = listFiles(args.workingDir, {".FT2", ".ft2"});
	if (ft2Files.empty())
	{
		std::cerr << "No FT2 files found.\n";
		return 1;
	}

	std::vector<std::string> hdf5Files = listFiles(args.hdf5Dir, {".h5"});
	if (hdf5Files.empty())
	{
		std::cerr << "No .h5 files in " << args.hdf5Dir << "\n";
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
	for (const auto &ft2 : ft2Files)
	{
		Ft1Data ft1Data;
		Ft1Data *ft1Ptr = nullptr;
		const std::string ft1Path = ft1PathForFt2(ft2);
		if (!ft1Path.empty())
		{
			if (readFt1File(ft1Path, ft1Data))
			{
				ft1Ptr = &ft1Data;
				std::cout << "FT1: " << ft1Path << "  (" << ft1Data.scans.size() << " scans)\n";
			}
			else
			{
				std::cerr << "Cannot read FT1 " << ft1Path
						  << "; MS1IsotopicAbundances will be 0 when FT1 peaks cannot be found.\n";
			}
		}
		else
		{
			std::cerr << "No matching FT1 file for " << ft2
					  << "; MS1IsotopicAbundances will be 0 when FT1 peaks cannot be found.\n";
		}
		int r = processOneFt2(args, ft2, hdf5Files, ft1Ptr);
		if (r != 0)
			rc = r;
	}
	std::cout << "sipros_search_spectra finished in " << (omp_get_wtime() - t0) << "s\n";
	return rc;
}
