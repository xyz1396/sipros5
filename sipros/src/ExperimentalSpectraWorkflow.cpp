#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <clocale>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <omp.h>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <H5Cpp.h>
#include "SiprosWorkflows.h"
#include "SiprosSearchRunner.h"
#include "isotopologue.h"
#include "proNovoConfig.h"
#include "RaxportHdf5Reader.h"
#include "spectraindex.h"
#include "ms2scan.h"
#include "performancelog.h"

#if !defined(_WIN32)
#include <cerrno>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#if !defined(_WIN32)
#include <unistd.h>
#endif
#endif

namespace fs = std::filesystem;

namespace
{
constexpr int kSpectraHdf5FormatVersion = 2;

#if !defined(_WIN32)
double timevalSeconds(const timeval &tv)
{
	return static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec) / 1.0e6;
}
#endif

double processTreeCpuSeconds()
{
#if !defined(_WIN32)
	rusage self{};
	rusage children{};
	getrusage(RUSAGE_SELF, &self);
	getrusage(RUSAGE_CHILDREN, &children);
	return timevalSeconds(self.ru_utime) + timevalSeconds(self.ru_stime) +
		   timevalSeconds(children.ru_utime) + timevalSeconds(children.ru_stime);
#else
        return sipros::processCpuSeconds();
#endif
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
					  size_t inputRows,
					  size_t retainedRows,
					  size_t baselineMatchedRows,
					  size_t outputFilesWritten) const
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
				  << "  PSM rows read:              " << inputRows << '\n'
				  << "  Retained peptide/charge:    " << retainedRows << '\n'
				  << "  Baseline matched spectra:   " << baselineMatchedRows << '\n'
				  << "  SFI output files written:   " << outputFilesWritten << '\n';
		if (totalWallSeconds > 0.0)
		{
			std::cout << "  Retained rows/sec:          "
					  << std::fixed << std::setprecision(0)
					  << static_cast<double>(retainedRows) / totalWallSeconds << '\n';
		}
		std::cout << '\n';
	}
};

struct Args
{
	std::string inputPath;
	std::string hdf5Path;
	std::string outputPath;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double fixedSipAbundancePct = 0.0;
	bool sipAbundanceRange = true;
	double sipAbundanceStartPct = 0.0;
	double sipAbundanceEndPct = 100.0;
	double sipAbundanceStepPct = 1.0;
	double probCutoff = 0.01;
	double ppmTolerance = 10.0;
	size_t minMatchedEnvelopes = 3;
	size_t envelopeTopN = 3;
	int threads = 0;
	bool writeDecoy = false;
	unsigned int decoySeed = 1;
	std::vector<std::string> fixedPtmSelectors;
};

struct PsmRow
{
	std::string psmId;
	std::string sample;
	int scanNumber = 0;
	std::string peptide;
	std::string proteins;
	int precursorCharge = 1;
	std::string retentionText = "0";
	double probability = 0.0;
	double svmScore = -std::numeric_limits<double>::infinity();
	double wdpScore = -std::numeric_limits<double>::infinity();
	size_t order = 0;
};

struct ReadStats
{
	size_t totalRows = 0;
	size_t shortRows = 0;
	size_t invalidRows = 0;
	size_t unsupportedMods = 0;
	size_t decoyRows = 0;
};

struct ObservedPeak
{
	double mz = 0.0;
	double intensity = 0.0;
	int charge = 1;
};

struct ObservedScan
{
	std::vector<ObservedPeak> peaks;
};

struct TheoreticalPeak
{
	double mz = 0.0;
	double probability = 0.0;
};

struct FragmentEntry
{
	double mz = 0.0;
	double theoreticalIntensity = 0.0;
	double experimentalIntensity = 0.0;
	char ionKind = 'b';
	size_t position = 0;
};

struct SpectrumOutputRecord
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
	std::vector<double> experimentalIntensity;
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
	std::vector<double> precursorMz;
	std::vector<double> precursorIntensity;
	std::vector<uint64_t> precursorOffset;
	std::vector<uint64_t> precursorCount;
	std::vector<double> fragmentMz;
	std::vector<double> theoreticalIntensity;
	std::vector<double> experimentalIntensity;
	std::vector<char> ionKind;
	std::vector<uint64_t> ionPosition;
	std::vector<uint64_t> fragmentOffset;
	std::vector<uint64_t> fragmentCount;

	void reserveRecords(size_t count)
	{
		psmIds.reserve(count);
		retentions.reserve(count);
		charges.reserve(count);
		peptides.reserve(count);
		proteins.reserve(count);
		precursorOffset.reserve(count);
		precursorCount.reserve(count);
		fragmentOffset.reserve(count);
		fragmentCount.reserve(count);
	}

	size_t recordCount() const
	{
		return psmIds.size();
	}
};

struct Hdf5OutputMetadata
{
	std::string recordKind;
	std::string chemistryProfileId;
	double targetSipAbundancePct = 0.0;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double probCutoff = 0.0;
	double ppmTolerance = 0.0;
	uint64_t minMatchedEnvelopes = 0;
	uint32_t envelopeTopN = 3;
};

struct Hdf5SampleTask
{
	std::string sample;
	fs::path path;
	std::unordered_set<int> requestedScans;
	std::vector<size_t> rowIndices;
};

struct OutputFileJob
{
	fs::path path;
	Hdf5OutputMetadata metadata;
	bool decoy = false;
	std::vector<double> targetAbundancesPct;
	size_t written = 0;
	sipros::SpectraIndexBuildStats buildStats;
	double generationSeconds = 0.0;
	double combineSeconds = 0.0;
	bool success = false;
	std::string error;
};

struct OutputJobStats
{
	uint64_t targetFailed = 0;
	uint64_t decoyComputeFailed = 0;
};

struct MatchedEnvelopeSet
{
	std::vector<std::pair<char, size_t>> envelopeKeys;
	std::map<std::pair<char, size_t>, double> apexIntensityByEnvelope;
	size_t retainedEnvelopes = 0;
	size_t matchedEnvelopes = 0;
};

struct PeptideTokens
{
	std::string nTermPrefix;
	std::vector<std::string> residues;
	std::string cTermSuffix;
};

struct ProcessingStats
{
	std::atomic<size_t> missingHdf5{0};
	std::atomic<size_t> missingScan{0};
	std::atomic<size_t> computeFailed{0};
	std::atomic<size_t> precursorFailed{0};
	std::atomic<size_t> unmatched{0};
	std::atomic<size_t> targetFailed{0};
	std::atomic<size_t> decoyCollision{0};
	std::atomic<size_t> decoyComputeFailed{0};
};

std::string peptideMassClassKey(const std::string &peptide);

void printUsage(const char *prog)
{
	std::cerr << "Usage: " << prog
			  << " -i <*_filtered_psms.tsv|filtered_dir> -f <h5_file|h5_dir> -o <output.sfi|output_dir/>"
			  << " -a <SIP atom/isotope, e.g. C13,H2,O18,N15,S34>"
			  << " [-b <fixed SIP pct|lower-upper, default 0-100>] [-s|--step <pct, default 1.0>]"
			  << " [-p <prob cutoff, default 0.01>]"
			  << " [--ppm <match tolerance, default 10>] [--min-matched-envelopes <N, default 3>]"
			  << " [--envelope-top-n <N, default 3>]"
			  << " [--decoy] [--decoy-seed <N, default 1>]"
			  << " [--fixed-ptm <name|default|none|all>] [-t <threads>]\n";
	std::cerr << "HDF5 MS2 matching is always performed at the natural-abundance C13 baseline; -b controls shifted output abundance(s).\n";
	std::cerr << "The output is one target SFI plus, with --decoy, one generated-decoy SFI; each contains every requested SIP percentage.\n";
	std::cerr << "--fixed-ptm is repeatable; omit it to use the compiled default (carbamidomethyl C).\n";
}

bool parseDoubleStrict(const std::string &text, double &value)
{
	try
	{
		const std::string t = sipros::TextUtils::trim(text);
		if (t.empty())
		{
			return false;
		}
		size_t consumed = 0;
		value = std::stod(t, &consumed);
		return consumed == t.size();
	}
	catch (const std::exception &)
	{
		return false;
	}
}

bool parseSipAbundanceSpec(const std::string &spec, Args &args)
{
	const std::string t = sipros::TextUtils::trim(spec);
	const size_t dash = t.find('-', 1);
	if (dash != std::string::npos)
	{
		double startPct = 0.0;
		double endPct = 0.0;
		if (!parseDoubleStrict(t.substr(0, dash), startPct) ||
			!parseDoubleStrict(t.substr(dash + 1), endPct))
		{
			return false;
		}
		args.sipAbundanceRange = true;
		args.sipAbundanceStartPct = startPct;
		args.sipAbundanceEndPct = endPct;
		args.fixedSipAbundancePct = startPct;
		return true;
	}

	double pct = 0.0;
	if (!parseDoubleStrict(t, pct))
	{
		return false;
	}
	args.sipAbundanceRange = false;
	args.fixedSipAbundancePct = pct;
	args.sipAbundanceStartPct = pct;
	args.sipAbundanceEndPct = pct;
	return true;
}

bool parseArgs(int argc, char **argv, Args &args)
{
	for (int i = 1; i < argc; ++i)
	{
		const std::string opt = argv[i];
		const auto requireValue = [&](const std::string &name) -> bool
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Missing value for option: " << name << "\n";
				return false;
			}
			return true;
		};

		if (opt == "-h" || opt == "--help")
		{
			printUsage(argv[0]);
			return false;
		}
		if (opt == "-i")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			args.inputPath = argv[++i];
		}
		else if (opt == "-f")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			args.hdf5Path = argv[++i];
		}
		else if (opt == "-o")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			args.outputPath = argv[++i];
		}
		else if (opt == "-a" || opt == "--sip-atom")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			const std::string sipSpec = argv[++i];
			if (!sipros::TextUtils::parseSipAtomSpec(sipSpec, args.sipAtom, args.sipIsotopeMassNumber))
			{
				std::cerr << "Invalid SIP atom/isotope: " << sipSpec << "\n";
				return false;
			}
		}
		else if (opt == "-b" || opt == "--sip-abundance")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			const std::string abundanceSpec = argv[++i];
			if (!parseSipAbundanceSpec(abundanceSpec, args))
			{
				std::cerr << "Invalid SIP abundance value: " << abundanceSpec
						  << ". Use a percentage or a range like 0-10.\n";
				return false;
			}
		}
		else if (opt == "-s" || opt == "--step")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			if (!parseDoubleStrict(argv[++i], args.sipAbundanceStepPct))
			{
				std::cerr << "Invalid SIP abundance step percentage.\n";
				return false;
			}
		}
		else if (opt == "-p" || opt == "--prob-cutoff")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				args.probCutoff = std::stod(argv[++i]);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid probability cutoff.\n";
				return false;
			}
		}
		else if (opt == "--ppm")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				args.ppmTolerance = std::stod(argv[++i]);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid ppm tolerance.\n";
				return false;
			}
		}
		else if (opt == "--min-matched-envelopes")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				args.minMatchedEnvelopes = static_cast<size_t>(std::stoul(argv[++i]));
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid minimum matched envelope count.\n";
				return false;
			}
		}
		else if (opt == "--envelope-top-n")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				args.envelopeTopN = static_cast<size_t>(std::stoul(argv[++i]));
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid SFI envelope top-N peak count.\n";
				return false;
			}
		}
		else if (opt == "--fixed-ptm")
		{
			if (!requireValue(opt))
			{
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
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				args.threads = std::stoi(argv[++i]);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid thread count.\n";
				return false;
			}
		}
		else if (opt == "--decoy")
		{
			args.writeDecoy = true;
		}
		else if (opt == "--decoy-seed")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			try
			{
				const long long seed = std::stoll(argv[++i]);
				if (seed < 0)
				{
					std::cerr << "Decoy seed must be >= 0.\n";
					return false;
				}
				args.decoySeed = static_cast<unsigned int>(seed);
			}
			catch (const std::exception &)
			{
				std::cerr << "Invalid decoy seed.\n";
				return false;
			}
		}
		else
		{
			std::cerr << "Unknown option: " << opt << "\n";
			return false;
		}
	}

	if (args.inputPath.empty() || args.hdf5Path.empty() || args.outputPath.empty() ||
		args.sipAtom == '\0')
	{
		printUsage(argv[0]);
		return false;
	}
	if (args.fixedSipAbundancePct < 0.0 || args.fixedSipAbundancePct > 100.0)
	{
		std::cerr << "SIP abundance percentage must be in [0, 100].\n";
		return false;
	}
	if (args.sipAbundanceRange)
	{
		if (args.sipAbundanceStartPct < 0.0 || args.sipAbundanceStartPct > 100.0 ||
			args.sipAbundanceEndPct < 0.0 || args.sipAbundanceEndPct > 100.0)
		{
			std::cerr << "SIP abundance range bounds must be in [0, 100].\n";
			return false;
		}
		if (args.sipAbundanceStartPct > args.sipAbundanceEndPct)
		{
			std::cerr << "SIP abundance range lower bound must be <= upper bound.\n";
			return false;
		}
	}
	if (args.sipAbundanceStepPct <= 0.0)
	{
		std::cerr << "SIP abundance step must be > 0.\n";
		return false;
	}
	if (args.probCutoff < 0.0 || args.probCutoff > 1.0)
	{
		std::cerr << "Probability cutoff must be in [0, 1].\n";
		return false;
	}
	if (args.ppmTolerance <= 0.0)
	{
		std::cerr << "PPM tolerance must be > 0.\n";
		return false;
	}
	if (args.minMatchedEnvelopes == 0)
	{
		std::cerr << "Minimum matched envelope count must be > 0.\n";
		return false;
	}
	if (args.envelopeTopN == 0)
	{
		std::cerr << "SFI envelope top-N peak count must be > 0.\n";
		return false;
	}
	if (args.threads < 0)
	{
		std::cerr << "Thread count must be >= 0.\n";
		return false;
	}
	return true;
}

bool parseDoubleField(const std::string &s, double &value)
{
	try
	{
		size_t consumed = 0;
		const std::string t = sipros::TextUtils::trim(s);
		if (t.empty())
		{
			return false;
		}
		value = std::stod(t, &consumed);
		return consumed > 0;
	}
	catch (const std::exception &)
	{
		return false;
	}
}

bool parseIntField(const std::string &s, int &value)
{
	try
	{
		size_t consumed = 0;
		const std::string t = sipros::TextUtils::trim(s);
		if (t.empty())
		{
			return false;
		}
		value = std::stoi(t, &consumed);
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

std::vector<std::string> splitDot(const std::string &s)
{
	std::vector<std::string> parts;
	size_t start = 0;
	while (true)
	{
		const size_t dot = s.find('.', start);
		if (dot == std::string::npos)
		{
			parts.push_back(s.substr(start));
			break;
		}
		parts.push_back(s.substr(start, dot - start));
		start = dot + 1;
	}
	return parts;
}

bool parseSpectrumId(const std::string &spectrum, std::string &sample, int &scanNumber, int &charge)
{
	const std::vector<std::string> parts = splitDot(sipros::TextUtils::trim(spectrum));
	if (parts.size() < 3)
	{
		return false;
	}
	sample = parts[0];
	if (sample.empty())
	{
		return false;
	}

	// Standard IDs end in .scan.charge (or .scan.scan.charge), while SIP
	// FASTA-search IDs may insert an abundance token before those numeric
	// fields: sample.SIP_C13_050_000Pct.scan.charge. Parse from the end so both
	// forms feed spectra generation without rewriting Aerith output.
	int parsedCharge = 0;
	if (!parseIntField(parts.back(), parsedCharge) || parsedCharge <= 0)
		return false;
	charge = parsedCharge;
	bool foundScan = false;
	for (size_t reverseIndex = parts.size() - 1;
		 reverseIndex > 0; --reverseIndex)
	{
		int parsedScan = 0;
		if (parseIntField(parts[reverseIndex - 1], parsedScan) &&
			parsedScan > 0)
		{
			scanNumber = parsedScan;
			foundScan = true;
			break;
		}
	}
	return foundScan;
}

std::vector<fs::path> collectPsmFiles(const std::string &inputPath)
{
	const fs::path path(inputPath);
	if (!fs::exists(path))
	{
		throw std::runtime_error("Input path does not exist: " + inputPath);
	}

	std::vector<fs::path> files;
	const std::string filteredSuffix = "_filtered_psms.tsv";
	const auto isFilteredPsm = [&](const fs::path &candidate)
	{
		const std::string filename = candidate.filename().string();
		return filename.size() >= filteredSuffix.size() &&
			filename.compare(filename.size() - filteredSuffix.size(),
				filteredSuffix.size(), filteredSuffix) == 0 &&
			filename != "SIP_filtered_psms.tsv";
	};
	if (fs::is_regular_file(path))
	{
		if (!isFilteredPsm(path))
		{
			throw std::runtime_error(
				"Experimental SFI input must be an Aerith "
				"*_filtered_psms.tsv file: " + inputPath);
		}
		files.push_back(path);
	}
	else if (fs::is_directory(path))
	{
		std::vector<fs::path> filteredFiles;
		for (const auto &entry : fs::recursive_directory_iterator(path))
		{
			if (!entry.is_regular_file())
			{
				continue;
			}
			if (isFilteredPsm(entry.path()))
			{
				filteredFiles.push_back(entry.path());
			}
		}
		files = std::move(filteredFiles);
	}
	std::sort(files.begin(), files.end());
	if (files.empty())
	{
		throw std::runtime_error(
			"No *_filtered_psms.tsv files found under: " + inputPath);
	}
	return files;
}

size_t getRequiredColumn(const std::unordered_map<std::string, size_t> &columns, const std::string &name)
{
	const auto it = columns.find(name);
	if (it == columns.end())
	{
		throw std::runtime_error("Missing required column: " + name);
	}
	return it->second;
}

size_t getOptionalColumn(const std::unordered_map<std::string, size_t> &columns, const std::string &name)
{
	const auto it = columns.find(name);
	return it == columns.end() ? std::string::npos : it->second;
}

size_t getRequiredColumnAny(const std::unordered_map<std::string, size_t> &columns,
							const std::vector<std::string> &names)
{
	for (const std::string &name : names)
	{
		const size_t idx = getOptionalColumn(columns, name);
		if (idx != std::string::npos)
		{
			return idx;
		}
	}
	std::string msg = "Missing required PSM TSV column; expected one of:";
	for (const std::string &name : names)
	{
		msg += " ";
		msg += name;
	}
	throw std::runtime_error(msg);
}

std::string requireProteinNames(const std::string &value, const std::string &context)
{
	std::string proteins = sipros::TextUtils::trim(value);
	if (proteins.empty())
	{
		throw std::runtime_error("Missing protein names for " + context);
	}
	return proteins;
}

std::vector<std::string> splitProteinList(const std::string &proteins)
{
	std::string inner = sipros::TextUtils::trim(proteins);
	if (inner.size() >= 2 && inner.front() == '{' && inner.back() == '}')
	{
		inner = inner.substr(1, inner.size() - 2);
	}

	std::vector<std::string> out;
	std::stringstream ss(inner);
	std::string token;
	while (std::getline(ss, token, ','))
	{
		token = sipros::TextUtils::trim(token);
		if (!token.empty())
		{
			out.push_back(token);
		}
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

std::string combineProteinColumns(const std::string &protein,
								  const std::string &mappedProteins,
								  const std::string &context)
{
	const std::string primary = requireProteinNames(protein, context);
	return mergeProteinLists(primary, mappedProteins);
}

std::string requireRetentionSeconds(const std::string &value, const std::string &context)
{
	std::string seconds = sipros::TextUtils::trim(value);
	if (seconds.empty())
	{
		throw std::runtime_error("Missing retention time in seconds for " + context);
	}

	size_t parsedChars = 0;
	try
	{
		(void)std::stod(seconds, &parsedChars);
	}
	catch (const std::exception &)
	{
		throw std::runtime_error("Retention time must be numeric seconds for " + context +
								 ": " + value);
	}

	if (!sipros::TextUtils::trim(seconds.substr(parsedChars)).empty())
	{
		throw std::runtime_error("Retention time must be numeric seconds without units for " +
								 context + ": " + value);
	}
	return seconds;
}

std::string decoyProteinNames(const std::string &value)
{
	std::string proteins = requireProteinNames(value, "decoy spectrum");
	const bool braced = proteins.size() >= 2 && proteins.front() == '{' && proteins.back() == '}';
	std::string inner = braced ? proteins.substr(1, proteins.size() - 2) : proteins;
	std::stringstream in(inner);
	std::string token;
	std::string out;
	while (std::getline(in, token, ','))
	{
		token = sipros::TextUtils::trim(token);
		if (token.empty())
		{
			continue;
		}
		if (token.rfind("DECOY_", 0) != 0)
		{
			token = "DECOY_" + token;
		}
		if (!out.empty())
		{
			out += ",";
		}
		out += token;
	}
	if (out.empty())
	{
		throw std::runtime_error("Missing protein names for decoy spectrum");
	}
	return braced ? "{" + out + "}" : out;
}

void readPsmFile(const fs::path &path, std::vector<PsmRow> &rows, ReadStats &stats, size_t &nextOrder)
{
	std::ifstream in(path);
	if (!in)
	{
		throw std::runtime_error("Cannot open PSM TSV: " + path.string());
	}

	std::string headerLine;
	if (!std::getline(in, headerLine))
	{
		throw std::runtime_error("PSM TSV is empty: " + path.string());
	}

	const std::vector<std::string> headers = sipros::TextUtils::splitTab(headerLine);
	std::unordered_map<std::string, size_t> columns;
	for (size_t i = 0; i < headers.size(); ++i)
	{
		columns[sipros::TextUtils::trim(headers[i])] = i;
	}

	const size_t idxSpectrum = getRequiredColumnAny(
		columns, {"Spectrum", "PSMId", "SpecId"});
	const size_t idxPeptide = getRequiredColumn(columns, "Peptide");
	const size_t idxCharge = getRequiredColumnAny(
		columns, {"Charge", "parentCharges"});
	const size_t idxRetention = getRequiredColumnAny(
		columns, {"Retention", "retentiontime"});
	const size_t idxPosteriorError =
		getRequiredColumn(columns, "posterior_error_prob");
	const size_t idxModifiedPeptide = getOptionalColumn(columns, "Modified Peptide") != std::string::npos
		? getOptionalColumn(columns, "Modified Peptide")
		: getOptionalColumn(columns, "ModifiedPeptide");
	const size_t idxAssignedModifications = getOptionalColumn(columns, "Assigned Modifications") != std::string::npos
		? getOptionalColumn(columns, "Assigned Modifications")
		: getOptionalColumn(columns, "AssignedModifications");
	const size_t idxSvmScore = getRequiredColumn(columns, "SVMscore");
	const size_t idxWdpScore = getRequiredColumn(columns, "WDPscores");
	const size_t idxLabel = getOptionalColumn(columns, "Label");
	const size_t idxMappedProteins = getOptionalColumn(columns, "Mapped Proteins");
	const size_t idxProteins = getRequiredColumnAny(columns, {"Proteins", "ProteinNames", "proteinNames",
															  "ProteinName", "proteinName", "Protein", "protein"});

	size_t requiredMax = std::max({
		idxSpectrum, idxPeptide, idxCharge, idxRetention,
		idxPosteriorError, idxSvmScore, idxWdpScore, idxProteins});
	if (idxLabel != std::string::npos)
		requiredMax = std::max(requiredMax, idxLabel);
	if (idxModifiedPeptide != std::string::npos)
	{
		requiredMax = std::max(requiredMax, idxModifiedPeptide);
	}
	if (idxAssignedModifications != std::string::npos)
	{
		requiredMax = std::max(requiredMax, idxAssignedModifications);
	}

	std::string line;
	while (std::getline(in, line))
	{
		if (line.empty())
		{
			continue;
		}
		++stats.totalRows;
		const std::vector<std::string> fields = sipros::TextUtils::splitTab(line);
		if (fields.size() <= requiredMax)
		{
			++stats.shortRows;
			continue;
		}

		PsmRow row;
		if (idxLabel != std::string::npos)
		{
			int label = 0;
			if (!parseIntField(fields[idxLabel], label))
			{
				++stats.invalidRows;
				continue;
			}
			if (label != 1)
			{
				++stats.decoyRows;
				continue;
			}
		}
		row.psmId = sipros::TextUtils::trim(fields[idxSpectrum]);
		if (row.psmId.empty())
		{
			++stats.invalidRows;
			continue;
		}

		if (!parseSpectrumId(row.psmId, row.sample, row.scanNumber, row.precursorCharge))
		{
			++stats.invalidRows;
			continue;
		}

		int charge = 0;
		if (parseIntField(fields[idxCharge], charge) && charge > 0)
		{
			row.precursorCharge = charge;
		}
		row.precursorCharge = std::max(1, row.precursorCharge);
		if (sipros::TextUtils::trim(headers[idxRetention]) == "retentiontime")
		{
			double retentionMinutes = 0.0;
			if (!parseDoubleField(fields[idxRetention], retentionMinutes) ||
				!std::isfinite(retentionMinutes))
			{
				++stats.invalidRows;
				continue;
			}
			std::ostringstream retentionSeconds;
			retentionSeconds << std::setprecision(17)
				<< retentionMinutes * 60.0;
			row.retentionText = retentionSeconds.str();
		}
		else
		{
			row.retentionText = requireRetentionSeconds(fields[idxRetention],
				path.string() + " PSM " + row.psmId);
		}
		const std::string mappedProteins =
			(idxMappedProteins != std::string::npos && fields.size() > idxMappedProteins)
				? fields[idxMappedProteins]
				: std::string();
		row.proteins = combineProteinColumns(fields[idxProteins], mappedProteins,
											 path.string() + " PSM " + row.psmId);
		double posteriorError = 0.0;
		if (!parseDoubleField(fields[idxPosteriorError], posteriorError) ||
			!std::isfinite(posteriorError) || posteriorError < 0.0 ||
			posteriorError > 1.0)
		{
			throw std::runtime_error(
				"Invalid posterior_error_prob in " + path.string() +
				" for PSM " + row.psmId);
		}
		row.probability = 1.0 - posteriorError;
		if (!parseDoubleField(fields[idxSvmScore], row.svmScore) ||
			!std::isfinite(row.svmScore))
		{
			throw std::runtime_error(
				"Invalid SVMscore in " + path.string() +
				" for PSM " + row.psmId);
		}
		if (!parseDoubleField(fields[idxWdpScore], row.wdpScore) ||
			!std::isfinite(row.wdpScore))
		{
			throw std::runtime_error(
				"Invalid WDPscores in " + path.string() +
				" for PSM " + row.psmId);
		}

		const std::string modified = (idxModifiedPeptide != std::string::npos && fields.size() > idxModifiedPeptide)
										 ? fields[idxModifiedPeptide]
										 : std::string();
		const std::string assignedModifications =
			(idxAssignedModifications != std::string::npos &&
			 fields.size() > idxAssignedModifications)
				? fields[idxAssignedModifications]
				: std::string();
		std::string reason;
		if (!convertModifiedPeptide(
				fields[idxPeptide], modified, assignedModifications,
				row.peptide, reason))
		{
			++stats.unsupportedMods;
			continue;
		}
		row.order = nextOrder++;
		rows.push_back(std::move(row));
	}
}

std::vector<PsmRow> readInputRows(const std::string &inputPath, ReadStats &stats)
{
	const std::vector<fs::path> files = collectPsmFiles(inputPath);
	std::vector<PsmRow> rows;
	size_t nextOrder = 0;
	for (const fs::path &path : files)
	{
		readPsmFile(path, rows, stats, nextOrder);
	}
	return rows;
}

std::map<std::string, size_t> requestedScanCountsBySample(const std::vector<PsmRow> &rows)
{
	std::map<std::string, std::set<int>> scansBySample;
	for (const PsmRow &row : rows)
	{
		scansBySample[row.sample].insert(row.scanNumber);
	}

	std::map<std::string, size_t> counts;
	for (const auto &kv : scansBySample)
	{
		counts[kv.first] = kv.second.size();
	}
	return counts;
}

bool isBetterPsm(const PsmRow &candidate, const PsmRow &current)
{
	constexpr double eps = 1e-15;
	if (candidate.probability > current.probability + eps)
	{
		return true;
	}
	if (std::abs(candidate.probability - current.probability) <= eps)
	{
		if (candidate.svmScore > current.svmScore + eps)
		{
			return true;
		}
		if (std::abs(candidate.svmScore - current.svmScore) <= eps)
		{
			if (candidate.wdpScore > current.wdpScore + eps)
			{
				return true;
			}
			if (std::abs(candidate.wdpScore - current.wdpScore) <= eps)
			{
				return candidate.order < current.order;
			}
		}
	}
	return false;
}

std::vector<PsmRow> selectBestRowsByPeptideCharge(const std::vector<PsmRow> &rows)
{
	std::unordered_map<std::string, size_t> bestIndex;
	std::vector<PsmRow> selected;
	for (const PsmRow &row : rows)
	{
		const std::string key = peptideMassClassKey(row.peptide) + "\t" + std::to_string(row.precursorCharge);
		const auto it = bestIndex.find(key);
		if (it == bestIndex.end())
		{
			bestIndex[key] = selected.size();
			selected.push_back(row);
			continue;
		}
		PsmRow &current = selected[it->second];
		const std::string mergedProteins = mergeProteinLists(current.proteins, row.proteins);
		if (isBetterPsm(row, current))
		{
			current = row;
			current.proteins = mergedProteins;
		}
		else
		{
			current.proteins = mergedProteins;
		}
	}

	std::sort(selected.begin(), selected.end(),
			  [](const PsmRow &a, const PsmRow &b)
			  { return a.order < b.order; });
	return selected;
}

std::unordered_map<std::string, fs::path> collectHdf5Files(const std::string &hdf5Path)
{
	const fs::path path(hdf5Path);
	if (!fs::exists(path))
	{
		throw std::runtime_error("Raxport HDF5 path does not exist: " + hdf5Path);
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
			if (!entry.is_regular_file())
			{
				continue;
			}
			const std::string ext = sipros::TextUtils::toLower(entry.path().extension().string());
			if (ext == ".h5" || ext == ".hdf5")
			{
				files.push_back(entry.path());
			}
		}
	}
	std::sort(files.begin(), files.end());
	if (files.empty())
	{
		throw std::runtime_error("No Raxport HDF5 files found under: " + hdf5Path);
	}

	std::unordered_map<std::string, fs::path> byStem;
	for (const fs::path &file : files)
	{
		const std::string stem = file.stem().string();
		if (byStem.find(stem) == byStem.end())
		{
			byStem[stem] = file;
		}
	}
	return byStem;
}

std::unordered_map<int, ObservedScan> readRequestedScans(const fs::path &hdf5File,
                                                         const std::unordered_set<int> &requestedScans)
{
	std::unordered_map<int, ObservedScan> scans;
	if (requestedScans.empty())
	{
		return scans;
	}

	std::vector<MS2Scan *> ms2Scans;
	std::string error;
	if (!sipros::readRaxportHdf5Scans(hdf5File.string(), ms2Scans, nullptr, error, &requestedScans))
	{
		throw std::runtime_error(error);
	}

	for (MS2Scan *scan : ms2Scans)
	{
		if (scan == nullptr)
		{
			continue;
		}
		ObservedScan observed;
		observed.peaks.reserve(scan->vdMZ.size());
		for (size_t i = 0; i < scan->vdMZ.size(); ++i)
		{
			ObservedPeak peak;
			peak.mz = scan->vdMZ[i];
			peak.intensity = i < scan->vdIntensity.size() ? scan->vdIntensity[i] : 0.0;
			const int charge = i < scan->viCharge.size() ? scan->viCharge[i] : 0;
			peak.charge = charge;
			observed.peaks.push_back(peak);
		}
		scans[scan->iScanId] = std::move(observed);
		delete scan;
	}
	ms2Scans.clear();
	return scans;
}

std::vector<Hdf5SampleTask> buildHdf5SampleTasks(const std::vector<PsmRow> &rows,
											   const std::unordered_map<std::string, fs::path> &hdf5FilesBySample,
											   ProcessingStats &processingStats)
{
	struct SampleRows
	{
		std::unordered_set<int> requestedScans;
		std::vector<size_t> rowIndices;
	};

	std::map<std::string, SampleRows> groupedBySample;
	for (size_t rowIndex = 0; rowIndex < rows.size(); ++rowIndex)
	{
		const PsmRow &row = rows[rowIndex];
		SampleRows &group = groupedBySample[row.sample];
		group.requestedScans.insert(row.scanNumber);
		group.rowIndices.push_back(rowIndex);
	}

	std::vector<Hdf5SampleTask> tasks;
	tasks.reserve(groupedBySample.size());
	for (auto &kv : groupedBySample)
	{
		const auto hdf5It = hdf5FilesBySample.find(kv.first);
		if (hdf5It == hdf5FilesBySample.end())
		{
			std::cerr << "No Raxport HDF5 file found for sample: " << kv.first << "\n";
			processingStats.missingHdf5.fetch_add(kv.second.rowIndices.size(), std::memory_order_relaxed);
			continue;
		}

		Hdf5SampleTask task;
		task.sample = kv.first;
		task.path = hdf5It->second;
		task.requestedScans = std::move(kv.second.requestedScans);
		task.rowIndices = std::move(kv.second.rowIndices);
		tasks.push_back(std::move(task));
	}
	return tasks;
}

void buildPrecursorChargePeaks(const IsotopeDistribution &dist,
							   int charge,
							   double probCutoff,
							   size_t envelopeTopN,
							   std::vector<double> &mzs,
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
	if (peaks.size() > envelopeTopN)
	{
		std::stable_sort(peaks.begin(), peaks.end(),
			[](const std::pair<double, double> &a,
			   const std::pair<double, double> &b)
			{
				if (a.second != b.second)
					return a.second > b.second;
				return a.first < b.first;
			});
		peaks.resize(envelopeTopN);
	}
	std::sort(peaks.begin(), peaks.end(),
			  [](const std::pair<double, double> &a, const std::pair<double, double> &b)
			  { return a.first < b.first; });

	double apexProbability = 0.0;
	for (const auto &peak : peaks)
	{
		apexProbability = std::max(apexProbability, peak.second);
	}

	mzs.clear();
	intensities.clear();
	mzs.reserve(peaks.size());
	intensities.reserve(peaks.size());
	for (const auto &peak : peaks)
	{
		mzs.push_back(peak.first);
		intensities.push_back(apexProbability > 0.0 ? peak.second / apexProbability : peak.second);
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

bool findMatchedPeak(const std::vector<ObservedPeak> &peaks,
					 double targetMz,
					 double ppmTolerance,
					 int theoreticalCharge,
					 ObservedPeak &matchedPeak)
{
	const double tolerance = targetMz * ppmTolerance / 1000000.0;
	const double lowerMz = targetMz - tolerance;
	const double upperMz = targetMz + tolerance;
	auto it = std::lower_bound(peaks.begin(), peaks.end(), lowerMz,
							   [](const ObservedPeak &peak, double value)
							   { return peak.mz < value; });

	bool found = false;
	double bestError = std::numeric_limits<double>::infinity();
	for (; it != peaks.end() && it->mz <= upperMz; ++it)
	{
		if (!sipros::observedPeakChargeMatches(
				it->charge, theoreticalCharge))
		{
			continue;
		}
		const double error = std::abs(it->mz - targetMz);
		if (!found || error < bestError ||
			(std::abs(error - bestError) <= 1e-12 && it->intensity > matchedPeak.intensity))
		{
			found = true;
			bestError = error;
			matchedPeak = *it;
		}
	}
	return found;
}

bool collectMatchedEnvelopeApexIntensities(const std::vector<std::vector<double>> &bMass,
										   const std::vector<std::vector<double>> &bProb,
										   const std::vector<std::vector<double>> &yMass,
										   const std::vector<std::vector<double>> &yProb,
										   const ObservedScan &scan,
										   double probCutoff,
										   double ppmTolerance,
										   size_t minMatchedEnvelopes,
										   MatchedEnvelopeSet &matchedSet)
{
	const double proton = ProNovoConfig::getProtonMass();
	matchedSet = MatchedEnvelopeSet();

	const auto collect = [&](const std::vector<std::vector<double>> &masses,
							 const std::vector<std::vector<double>> &probs,
							 char ionKind)
	{
		for (size_t i = 0; i < masses.size() && i < probs.size(); ++i)
		{
			std::vector<TheoreticalPeak> envelope;
			envelope.reserve(std::min(masses[i].size(), probs[i].size()));
			for (size_t j = 0; j < masses[i].size() && j < probs[i].size(); ++j)
			{
				if (probs[i][j] < probCutoff)
				{
					continue;
				}
				TheoreticalPeak peak;
				peak.mz = masses[i][j] + proton;
				peak.probability = probs[i][j];
				envelope.push_back(peak);
			}
			if (envelope.empty())
			{
				continue;
			}
			++matchedSet.retainedEnvelopes;
			const std::pair<char, size_t> envelopeKey{ionKind, i + 1};
			matchedSet.envelopeKeys.push_back(envelopeKey);
			matchedSet.apexIntensityByEnvelope[envelopeKey] = 0.0;

			size_t apexIndex = 0;
			for (size_t j = 1; j < envelope.size(); ++j)
			{
				if (envelope[j].probability > envelope[apexIndex].probability)
				{
					apexIndex = j;
				}
			}

			ObservedPeak matched;
			const double apexProbability = envelope[apexIndex].probability;
			if (apexProbability <= 0.0)
			{
				continue;
			}

			const bool matchedEnvelope = findMatchedPeak(
				scan.peaks, envelope[apexIndex].mz, ppmTolerance, 1, matched);
			if (matchedEnvelope)
			{
				++matchedSet.matchedEnvelopes;
				matchedSet.apexIntensityByEnvelope[envelopeKey] = matched.intensity;
			}
		}
	};

	collect(bMass, bProb, 'b');
	collect(yMass, yProb, 'y');
	if (matchedSet.retainedEnvelopes == 0)
	{
		return false;
	}
	if (matchedSet.matchedEnvelopes < minMatchedEnvelopes)
	{
		return false;
	}
	return true;
}

void matchBaselineInHdf5Batches(const std::vector<PsmRow> &rows,
							   const std::unordered_map<std::string, fs::path> &hdf5FilesBySample,
							   const Isotopologue &baselineIso,
							   const Args &args,
							   int effectiveThreads,
							   std::vector<MatchedEnvelopeSet> &matchedEnvelopeSets,
							   std::vector<char> &baselineOk,
							   ProcessingStats &processingStats,
							   TimingLogger &timing)
{
	const std::vector<Hdf5SampleTask> tasks = buildHdf5SampleTasks(rows, hdf5FilesBySample, processingStats);
	const size_t batchSize = static_cast<size_t>(std::max(1, effectiveThreads));

	for (size_t batchStart = 0; batchStart < tasks.size(); batchStart += batchSize)
	{
		const size_t batchEnd = std::min(batchStart + batchSize, tasks.size());
		const size_t batchCount = batchEnd - batchStart;
		std::vector<std::unordered_map<int, ObservedScan>> loadedScans(batchCount);
		std::vector<std::string> errors(batchCount);

		{
			std::ostringstream label;
			label << "load HDF5 batch " << (batchStart / batchSize + 1)
				  << '/' << ((tasks.size() + batchSize - 1) / batchSize);
			timing.run("Load HDF5", label.str(), batchCount, "files", [&]()
			{
#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < static_cast<int>(batchCount); ++i)
				{
					const size_t localTaskIndex = static_cast<size_t>(i);
					const Hdf5SampleTask &task = tasks[batchStart + localTaskIndex];
#pragma omp critical(cerr_output)
					{
						std::cerr << "Reading " << task.requestedScans.size()
								  << " requested HDF5 MS2 scans from " << task.path << "\n";
					}
					try
					{
						loadedScans[localTaskIndex] = readRequestedScans(task.path, task.requestedScans);
					}
					catch (const std::exception &ex)
					{
						errors[localTaskIndex] = ex.what();
					}
				}
			});
		}

		for (const std::string &error : errors)
		{
			if (!error.empty())
			{
				throw std::runtime_error(error);
			}
		}

		struct BatchRowRef
		{
			size_t localTaskIndex = 0;
			size_t rowIndex = 0;
		};

		std::vector<BatchRowRef> batchRows;
		for (size_t localTaskIndex = 0; localTaskIndex < batchCount; ++localTaskIndex)
		{
			const Hdf5SampleTask &task = tasks[batchStart + localTaskIndex];
			for (size_t rowIndex : task.rowIndices)
			{
				batchRows.push_back({localTaskIndex, rowIndex});
			}
		}

		{
			std::ostringstream label;
			label << "baseline match batch " << (batchStart / batchSize + 1)
				  << '/' << ((tasks.size() + batchSize - 1) / batchSize);
			timing.run("Baseline Match", label.str(), batchRows.size(), "rows", [&]()
			{
#pragma omp parallel
				{
					Isotopologue localIso = baselineIso;
#pragma omp for schedule(dynamic)
					for (int i = 0; i < static_cast<int>(batchRows.size()); ++i)
					{
						const BatchRowRef &ref = batchRows[static_cast<size_t>(i)];
						const PsmRow &row = rows[ref.rowIndex];
						const auto scanIt = loadedScans[ref.localTaskIndex].find(row.scanNumber);
						if (scanIt == loadedScans[ref.localTaskIndex].end())
						{
							++processingStats.missingScan;
							continue;
						}

						std::vector<std::vector<double>> yMass, yProb, bMass, bProb;
						if (!localIso.computeProductIon(row.peptide, yMass, yProb, bMass, bProb))
						{
							++processingStats.computeFailed;
							continue;
						}

						if (!collectMatchedEnvelopeApexIntensities(bMass, bProb, yMass, yProb,
																   scanIt->second, args.probCutoff, args.ppmTolerance,
																   args.minMatchedEnvelopes,
																   matchedEnvelopeSets[ref.rowIndex]))
						{
							++processingStats.unmatched;
							continue;
						}

						baselineOk[ref.rowIndex] = 1;
					}
				}
			});
		}
	}
}

bool buildShiftedChargeOneFragmentEntries(const std::vector<std::vector<double>> &bMass,
										  const std::vector<std::vector<double>> &bProb,
										  const std::vector<std::vector<double>> &yMass,
										  const std::vector<std::vector<double>> &yProb,
										  const MatchedEnvelopeSet &matchedSet,
										  double probCutoff,
										  size_t envelopeTopN,
										  std::vector<FragmentEntry> &entries)
{
	const double proton = ProNovoConfig::getProtonMass();
	entries.clear();
	entries.reserve(512);

	const auto collect = [&](const std::vector<std::vector<double>> &masses,
							 const std::vector<std::vector<double>> &probs,
							 char ionKind)
	{
		for (size_t i = 0; i < masses.size() && i < probs.size(); ++i)
		{
			std::vector<TheoreticalPeak> envelope;
			envelope.reserve(std::min(masses[i].size(), probs[i].size()));
			for (size_t j = 0; j < masses[i].size() && j < probs[i].size(); ++j)
			{
				if (probs[i][j] < probCutoff)
				{
					continue;
				}
				TheoreticalPeak peak;
				peak.mz = masses[i][j] + proton;
				peak.probability = probs[i][j];
				envelope.push_back(peak);
			}
			if (envelope.empty())
			{
				continue;
			}
			if (envelope.size() > envelopeTopN)
			{
				std::stable_sort(envelope.begin(), envelope.end(),
					[](const TheoreticalPeak &a, const TheoreticalPeak &b)
					{
						if (a.probability != b.probability)
							return a.probability > b.probability;
						return a.mz < b.mz;
					});
				envelope.resize(envelopeTopN);
				std::sort(envelope.begin(), envelope.end(),
					[](const TheoreticalPeak &a, const TheoreticalPeak &b)
					{ return a.mz < b.mz; });
			}

			size_t apexIndex = 0;
			for (size_t j = 1; j < envelope.size(); ++j)
			{
				if (envelope[j].probability > envelope[apexIndex].probability)
				{
					apexIndex = j;
				}
			}

			const double apexProbability = envelope[apexIndex].probability;
			if (apexProbability <= 0.0)
			{
				continue;
			}

			const auto matchIt = matchedSet.apexIntensityByEnvelope.find({ionKind, i + 1});
			const double matchedApexIntensity =
				matchIt == matchedSet.apexIntensityByEnvelope.end() ? 0.0 : matchIt->second;

			for (const TheoreticalPeak &peak : envelope)
			{
				const double relativeIntensity = peak.probability / apexProbability;
				FragmentEntry entry;
				entry.mz = peak.mz;
				entry.theoreticalIntensity = relativeIntensity;
				entry.experimentalIntensity = matchedApexIntensity * relativeIntensity;
				entry.ionKind = ionKind;
				entry.position = i + 1;
				entries.push_back(entry);
			}
		}
	};

	collect(bMass, bProb, 'b');
	collect(yMass, yProb, 'y');
	if (entries.empty())
	{
		return false;
	}

	double maxIntensity = 0.0;
	for (const FragmentEntry &entry : entries)
	{
		maxIntensity = std::max(maxIntensity, entry.experimentalIntensity);
	}
	if (maxIntensity > 0.0)
	{
		for (FragmentEntry &entry : entries)
		{
			entry.experimentalIntensity /= maxIntensity;
		}
	}

	std::sort(entries.begin(), entries.end(),
			  [](const FragmentEntry &a, const FragmentEntry &b)
			  { return a.mz < b.mz; });
	return true;
}

bool parsePeptideTokens(const std::string &peptide, PeptideTokens &tokens)
{
	tokens = PeptideTokens();
	if (peptide.size() < 2 || peptide.front() != '[')
	{
		return false;
	}

	const size_t close = peptide.find(']');
	if (close == std::string::npos || close <= 1)
	{
		return false;
	}

	size_t i = 1;
	while (i < close && !std::isalpha(static_cast<unsigned char>(peptide[i])))
	{
		tokens.nTermPrefix.push_back(peptide[i]);
		++i;
	}

	while (i < close)
	{
		if (!std::isalpha(static_cast<unsigned char>(peptide[i])))
		{
			return false;
		}
		std::string token;
		token.push_back(peptide[i++]);
		while (i < close && !std::isalpha(static_cast<unsigned char>(peptide[i])))
		{
			token.push_back(peptide[i++]);
		}
		tokens.residues.push_back(token);
	}

	if (close + 1 < peptide.size())
	{
		tokens.cTermSuffix = peptide.substr(close + 1);
	}
	return !tokens.residues.empty();
}

std::string buildPeptideFromTokens(const PeptideTokens &tokens)
{
	std::string peptide = "[";
	peptide += tokens.nTermPrefix;
	for (const std::string &token : tokens.residues)
	{
		peptide += token;
	}
	peptide += "]";
	peptide += tokens.cTermSuffix;
	return peptide;
}

std::string peptideMassClassKey(const std::string &peptide)
{
	std::string key = peptide;
	for (char &c : key)
	{
		if (c == 'I' || c == 'J')
		{
			c = 'L';
		}
	}
	return key;
}

std::string adjustedDecoyRetention(const std::string &retentionText,
								   unsigned int seed,
								   size_t rowIndex,
								   bool decoyAddedResidue)
{
	const std::string secondsText = requireRetentionSeconds(retentionText, "decoy spectrum");
	const double retentionSeconds = std::stod(secondsText);

	std::mt19937 rng(seed ^
					 static_cast<unsigned int>((rowIndex + 1) * 2654435761u) ^
					 0x9E3779B9u);
	std::uniform_real_distribution<double> jitterSecondsDist(15.0, 60.0);
	const double signedJitterSeconds = (decoyAddedResidue ? 1.0 : -1.0) * jitterSecondsDist(rng);
	const double adjusted = retentionSeconds + signedJitterSeconds;

	std::ostringstream out;
	out.imbue(std::locale::classic());
	out << std::setprecision(10) << adjusted;
	return out.str();
}

std::vector<double> shuffledApexIntensityPool(const MatchedEnvelopeSet &matchedSet,
											  unsigned int seed,
											  size_t rowIndex)
{
	std::vector<double> intensities;
	intensities.reserve(matchedSet.envelopeKeys.size());
	for (const auto &key : matchedSet.envelopeKeys)
	{
		const auto it = matchedSet.apexIntensityByEnvelope.find(key);
		intensities.push_back(it == matchedSet.apexIntensityByEnvelope.end() ? 0.0 : it->second);
	}

	std::mt19937 rng(seed ^ static_cast<unsigned int>((rowIndex + 1) * 2654435761u));
	std::shuffle(intensities.begin(), intensities.end(), rng);
	return intensities;
}

double meanIntensity(const std::vector<double> &intensities)
{
	if (intensities.empty())
	{
		return 0.0;
	}
	double sum = 0.0;
	for (double intensity : intensities)
	{
		sum += intensity;
	}
	return sum / static_cast<double>(intensities.size());
}

bool generateDecoyPeptide(const std::string &targetPeptide,
						  const std::unordered_set<std::string> &experimentalPeptideMassClasses,
						  unsigned int seed,
						  size_t rowIndex,
						  std::string &decoyPeptide,
						  bool &decoyAddedResidue)
{
	PeptideTokens baseTokens;
	if (!parsePeptideTokens(targetPeptide, baseTokens))
	{
		return false;
	}

	constexpr size_t maxAttempts = 100;
	const std::string standardResidues = "ACDEFGHIKLMNPQRSTVWY";
	const int minPeptideLength = ProNovoConfig::getMinPeptideLength();
	const std::string targetMassClass = peptideMassClassKey(targetPeptide);
	for (size_t attempt = 0; attempt < maxAttempts; ++attempt)
	{
		PeptideTokens tokens = baseTokens;
		std::mt19937 rng(seed ^
						 static_cast<unsigned int>((rowIndex + 1) * 2654435761u) ^
						 static_cast<unsigned int>((attempt + 1) * 2246822519u));

		if (tokens.residues.size() > 3)
		{
			std::shuffle(tokens.residues.begin() + 1, tokens.residues.end() - 1, rng);
		}

		bool deleteResidue = (std::uniform_int_distribution<int>(0, 1)(rng) == 0);
		if (deleteResidue &&
			(tokens.residues.size() <= 2 ||
			 static_cast<int>(tokens.residues.size() - 1) < minPeptideLength))
		{
			deleteResidue = false;
		}

		if (deleteResidue)
		{
			std::uniform_int_distribution<size_t> deleteDist(1, tokens.residues.size() - 2);
			tokens.residues.erase(tokens.residues.begin() + static_cast<std::ptrdiff_t>(deleteDist(rng)));
		}
		else
		{
			std::uniform_int_distribution<size_t> aaDist(0, standardResidues.size() - 1);
			std::string inserted(1, standardResidues[aaDist(rng)]);
			const size_t lastInsertPosition = tokens.residues.size() > 1 ? tokens.residues.size() - 1 : tokens.residues.size();
			std::uniform_int_distribution<size_t> insertDist(1, lastInsertPosition);
			tokens.residues.insert(tokens.residues.begin() + static_cast<std::ptrdiff_t>(insertDist(rng)), inserted);
		}

		const std::string candidate = buildPeptideFromTokens(tokens);
		const std::string candidateMassClass = peptideMassClassKey(candidate);
		if (candidateMassClass != targetMassClass &&
			experimentalPeptideMassClasses.find(candidateMassClass) == experimentalPeptideMassClasses.end())
		{
			decoyPeptide = candidate;
			decoyAddedResidue = !deleteResidue;
			return true;
		}
	}
	return false;
}

bool buildDecoyChargeOneFragmentEntries(const std::vector<std::vector<double>> &bMass,
										const std::vector<std::vector<double>> &bProb,
										const std::vector<std::vector<double>> &yMass,
										const std::vector<std::vector<double>> &yProb,
										const std::vector<double> &shuffledApexIntensities,
										double meanApexIntensity,
										bool useOneMeanApex,
										double probCutoff,
										size_t envelopeTopN,
										std::vector<FragmentEntry> &entries)
{
	if (shuffledApexIntensities.empty() && (!useOneMeanApex || meanApexIntensity <= 0.0))
	{
		return false;
	}

	const double proton = ProNovoConfig::getProtonMass();
	entries.clear();
	entries.reserve(512);
	size_t nextEnvelopeIndex = 0;
	bool meanApexUsed = false;

	const auto collect = [&](const std::vector<std::vector<double>> &masses,
							 const std::vector<std::vector<double>> &probs,
							 char ionKind)
	{
		for (size_t i = 0; i < masses.size() && i < probs.size(); ++i)
		{
			std::vector<TheoreticalPeak> envelope;
			envelope.reserve(std::min(masses[i].size(), probs[i].size()));
			for (size_t j = 0; j < masses[i].size() && j < probs[i].size(); ++j)
			{
				if (probs[i][j] < probCutoff)
				{
					continue;
				}
				TheoreticalPeak peak;
				peak.mz = masses[i][j] + proton;
				peak.probability = probs[i][j];
				envelope.push_back(peak);
			}
			if (envelope.empty())
			{
				continue;
			}
			if (envelope.size() > envelopeTopN)
			{
				std::stable_sort(envelope.begin(), envelope.end(),
					[](const TheoreticalPeak &a, const TheoreticalPeak &b)
					{
						if (a.probability != b.probability)
							return a.probability > b.probability;
						return a.mz < b.mz;
					});
				envelope.resize(envelopeTopN);
				std::sort(envelope.begin(), envelope.end(),
					[](const TheoreticalPeak &a, const TheoreticalPeak &b)
					{ return a.mz < b.mz; });
			}

			size_t apexIndex = 0;
			for (size_t j = 1; j < envelope.size(); ++j)
			{
				if (envelope[j].probability > envelope[apexIndex].probability)
				{
					apexIndex = j;
				}
			}
			const double apexProbability = envelope[apexIndex].probability;
			if (apexProbability <= 0.0)
			{
				continue;
			}

			double assignedApexIntensity = 0.0;
			if (useOneMeanApex && !meanApexUsed)
			{
				assignedApexIntensity = meanApexIntensity;
				meanApexUsed = true;
			}
			else if (nextEnvelopeIndex < shuffledApexIntensities.size())
			{
				assignedApexIntensity = shuffledApexIntensities[nextEnvelopeIndex];
				++nextEnvelopeIndex;
			}
			else
			{
				return;
			}

			for (const TheoreticalPeak &peak : envelope)
			{
				const double relativeIntensity = peak.probability / apexProbability;
				FragmentEntry entry;
				entry.mz = peak.mz;
				entry.theoreticalIntensity = relativeIntensity;
				entry.experimentalIntensity = assignedApexIntensity * relativeIntensity;
				entry.ionKind = ionKind;
				entry.position = i + 1;
				entries.push_back(entry);
			}
		}
	};

	collect(bMass, bProb, 'b');
	collect(yMass, yProb, 'y');
	if (entries.empty())
	{
		return false;
	}

	double maxIntensity = 0.0;
	for (const FragmentEntry &entry : entries)
	{
		maxIntensity = std::max(maxIntensity, entry.experimentalIntensity);
	}
	if (maxIntensity > 0.0)
	{
		for (FragmentEntry &entry : entries)
		{
			entry.experimentalIntensity /= maxIntensity;
		}
	}

	std::sort(entries.begin(), entries.end(),
			  [](const FragmentEntry &a, const FragmentEntry &b)
			  { return a.mz < b.mz; });
	return true;
}

std::vector<double> makeTargetAbundances(const Args &args)
{
	if (!args.sipAbundanceRange)
	{
		return {args.fixedSipAbundancePct};
	}

	std::vector<double> abundances;
	const double eps = 1e-9;
	for (double pct = args.sipAbundanceStartPct;
		 pct <= args.sipAbundanceEndPct + eps;
		 pct += args.sipAbundanceStepPct)
	{
		abundances.push_back(std::min(pct, args.sipAbundanceEndPct));
	}
	if (abundances.empty() ||
		std::abs(abundances.back() - args.sipAbundanceEndPct) > eps)
	{
		abundances.push_back(args.sipAbundanceEndPct);
	}
	return abundances;
}

std::string formatAbundancePctForPath(double pct)
{
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(3) << std::setw(7) << std::setfill('0') << pct;
	return ss.str();
}

bool hasTrailingPathSeparator(const std::string &path)
{
	return !path.empty() && (path.back() == '/' || path.back() == '\\');
}

fs::path appendSfiExtensionIfNeeded(const fs::path &path)
{
	if (sipros::TextUtils::toLower(path.extension().string()) == ".sfi")
	{
		return path;
	}
	fs::path sfiPath = path;
	sfiPath += ".sfi";
	return sfiPath;
}

fs::path resolveOutputBasePath(const std::string &outputPath, bool multipleOutputFiles)
{
	const fs::path path(outputPath);
	if (hasTrailingPathSeparator(outputPath) ||
		(fs::exists(path) && fs::is_directory(path)))
	{
		return path / "spectra.sfi";
	}

	const fs::path sfiPath = appendSfiExtensionIfNeeded(path);
	if (!multipleOutputFiles)
	{
		return sfiPath;
	}

	if (sipros::TextUtils::toLower(path.extension().string()) == ".sfi")
	{
		throw std::runtime_error("Multiple output files requested, but -o names a single .sfi file. Use an output directory path, for example -o out or -o out/.");
	}

	return path / "spectra.sfi";
}

std::string sipLabelForPath(char sipAtom, int sipIsotopeMassNumber)
{
	std::string label(1, sipAtom);
	if (sipIsotopeMassNumber > 0)
	{
		label += std::to_string(sipIsotopeMassNumber);
	}
	return label;
}

fs::path outputPathForAbundance(const fs::path &basePath, const std::string &sipLabel, double pct)
{
	const std::string filename = basePath.stem().string() + "_" +
								 sipLabel + "_" +
								 formatAbundancePctForPath(pct) + "Pct" +
								 basePath.extension().string();
	const fs::path parent = basePath.parent_path();
	return parent.empty() ? fs::path(filename) : parent / filename;
}

fs::path decoyOutputPathForAbundance(const fs::path &basePath, const std::string &sipLabel, double pct)
{
	const std::string filename = basePath.stem().string() + "_Decoy_" +
								 sipLabel + "_" +
								 formatAbundancePctForPath(pct) + "Pct" +
								 basePath.extension().string();
	const fs::path parent = basePath.parent_path();
	return parent.empty() ? fs::path(filename) : parent / filename;
}

fs::path decoyOutputPath(const fs::path &basePath)
{
	const std::string filename = basePath.stem().string() + "_Decoy" +
		basePath.extension().string();
	const fs::path parent = basePath.parent_path();
	return parent.empty() ? fs::path(filename) : parent / filename;
}

std::string spectraHeaderId(const std::string &psmId, double targetSipAbundancePct, bool decoy)
{
	std::string id;
	if (decoy)
	{
		id += "DECOY_";
	}
	id += psmId;
	id += "_";
	id += formatAbundancePctForPath(targetSipAbundancePct);
	id += "Pct";
	return id;
}

void copyFragmentEntriesToRecord(const std::vector<FragmentEntry> &entries, SpectrumOutputRecord &record)
{
	record.fragmentMz.clear();
	record.theoreticalIntensity.clear();
	record.experimentalIntensity.clear();
	record.ionKinds.clear();
	record.ionPositions.clear();
	record.fragmentMz.reserve(entries.size());
	record.theoreticalIntensity.reserve(entries.size());
	record.experimentalIntensity.reserve(entries.size());
	record.ionKinds.reserve(entries.size());
	record.ionPositions.reserve(entries.size());
	for (const FragmentEntry &entry : entries)
	{
		record.fragmentMz.push_back(entry.mz);
		record.theoreticalIntensity.push_back(entry.theoreticalIntensity);
		record.experimentalIntensity.push_back(entry.experimentalIntensity);
		record.ionKinds.push_back(entry.ionKind);
		record.ionPositions.push_back(static_cast<uint64_t>(entry.position));
	}
}

H5::DataSpace createDataspace(size_t count)
{
	const hsize_t dim = static_cast<hsize_t>(count);
	return H5::DataSpace(1, &dim);
}

H5::DataSpace createScalarDataspace()
{
	return H5::DataSpace(H5S_SCALAR);
}

H5::DSetCreatPropList createCompressedDatasetProperties(size_t count, size_t chunkSize)
{
	H5::DSetCreatPropList plist;
	if (count > 0)
	{
		const hsize_t chunk = static_cast<hsize_t>(std::min(count, chunkSize));
		plist.setChunk(1, &chunk);
		plist.setShuffle();
		plist.setDeflate(6);
	}
	return plist;
}

template <typename T>
void writeVectorDataset(const H5::Group &group,
						const char *name,
						const H5::DataType &type,
						const std::vector<T> &values,
						bool compress = false,
						size_t chunkSize = 262144)
{
	H5::DataSpace space = createDataspace(values.size());
	H5::DataSet dataset = compress && !values.empty()
							  ? group.createDataSet(name, type, space, createCompressedDatasetProperties(values.size(), chunkSize))
							  : group.createDataSet(name, type, space);
	if (!values.empty())
	{
		dataset.write(values.data(), type);
	}
}

void writeStringDataset(const H5::Group &group,
						const char *name,
						const std::vector<std::string> &values,
						bool compress = false,
						size_t chunkSize = 262144)
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
		std::memcpy(flat.data() + i * width, values[i].c_str(), std::min(values[i].size(), width - 1));
	}

	H5::DataSpace space = createDataspace(values.size());
	H5::DataSet dataset = compress && !values.empty()
							  ? group.createDataSet(name, type, space, createCompressedDatasetProperties(values.size(), chunkSize))
							  : group.createDataSet(name, type, space);
	if (!values.empty())
	{
		dataset.write(flat.data(), type);
	}
}

void writeStringAttribute(const H5::H5Object &object, const char *name, const std::string &value)
{
	H5::StrType type(H5::PredType::C_S1, std::max<size_t>(1, value.size() + 1));
	type.setStrpad(H5T_STR_NULLTERM);
	type.setCset(H5T_CSET_UTF8);
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attr = object.createAttribute(name, type, space);
	attr.write(type, value.c_str());
}

void writeIntAttribute(const H5::H5Object &object, const char *name, int value)
{
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attr = object.createAttribute(name, H5::PredType::NATIVE_INT, space);
	attr.write(H5::PredType::NATIVE_INT, &value);
}

void writeUInt64Attribute(const H5::H5Object &object, const char *name, uint64_t value)
{
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attr = object.createAttribute(name, H5::PredType::NATIVE_UINT64, space);
	attr.write(H5::PredType::NATIVE_UINT64, &value);
}

void writeDoubleAttribute(const H5::H5Object &object, const char *name, double value)
{
	H5::DataSpace space = createScalarDataspace();
	H5::Attribute attr = object.createAttribute(name, H5::PredType::NATIVE_DOUBLE, space);
	attr.write(H5::PredType::NATIVE_DOUBLE, &value);
}

void writeSpectraAttributes(const H5::H5Object &file, const Hdf5OutputMetadata &metadata)
{
	writeIntAttribute(file, "format_version", kSpectraHdf5FormatVersion);
	writeStringAttribute(file, "chemistry_profile_id", metadata.chemistryProfileId);
	writeStringAttribute(file, "record_kind", metadata.recordKind);
	writeDoubleAttribute(file, "target_sip_abundance_pct", metadata.targetSipAbundancePct);
	writeStringAttribute(file, "sip_atom", std::string(1, metadata.sipAtom));
	writeIntAttribute(file, "sip_isotope_mass_number", metadata.sipIsotopeMassNumber);
	writeDoubleAttribute(file, "prob_cutoff", metadata.probCutoff);
	writeDoubleAttribute(file, "ppm_tolerance", metadata.ppmTolerance);
	writeUInt64Attribute(file, "min_matched_envelopes", metadata.minMatchedEnvelopes);
}

Hdf5OutputData buildOutputDataFromRecords(const std::vector<SpectrumOutputRecord> &records,
										  const std::vector<char> &ok)
{
	Hdf5OutputData output;
	const size_t n = std::min(records.size(), ok.size());
	constexpr size_t invalidIndex = std::numeric_limits<size_t>::max();
	std::vector<size_t> recordIndex(n, invalidIndex);
	std::vector<size_t> precursorOffsetByRow(n, 0);
	std::vector<size_t> fragmentOffsetByRow(n, 0);

	size_t recordCount = 0;
	size_t precursorValueCount = 0;
	size_t fragmentValueCount = 0;
	for (size_t i = 0; i < n; ++i)
	{
		if (!ok[i])
		{
			continue;
		}
		recordIndex[i] = recordCount++;
		precursorOffsetByRow[i] = precursorValueCount;
		fragmentOffsetByRow[i] = fragmentValueCount;
		precursorValueCount += records[i].precursorMz.size();
		fragmentValueCount += records[i].fragmentMz.size();
	}

	output.psmIds.resize(recordCount);
	output.retentions.resize(recordCount);
	output.charges.resize(recordCount);
	output.peptides.resize(recordCount);
	output.proteins.resize(recordCount);
	output.precursorOffset.resize(recordCount);
	output.precursorCount.resize(recordCount);
	output.fragmentOffset.resize(recordCount);
	output.fragmentCount.resize(recordCount);
	output.precursorMz.resize(precursorValueCount);
	output.precursorIntensity.resize(precursorValueCount);
	output.fragmentMz.resize(fragmentValueCount);
	output.theoreticalIntensity.resize(fragmentValueCount);
	output.experimentalIntensity.resize(fragmentValueCount);
	output.ionKind.resize(fragmentValueCount);
	output.ionPosition.resize(fragmentValueCount);

	for (int i = 0; i < static_cast<int>(n); ++i)
	{
		const size_t rowIndex = static_cast<size_t>(i);
		const size_t outIndex = recordIndex[rowIndex];
		if (outIndex == invalidIndex)
		{
			continue;
		}

		const SpectrumOutputRecord &record = records[rowIndex];
		const size_t precursorOffset = precursorOffsetByRow[rowIndex];
		const size_t fragmentOffset = fragmentOffsetByRow[rowIndex];

		output.psmIds[outIndex] = record.psmId;
		output.retentions[outIndex] = record.retention;
		output.charges[outIndex] = record.charge;
		output.peptides[outIndex] = record.peptide;
		output.proteins[outIndex] = requireProteinNames(record.proteins,
														"HDF5 record " + record.psmId);
		output.precursorOffset[outIndex] = static_cast<uint64_t>(precursorOffset);
		output.precursorCount[outIndex] = static_cast<uint64_t>(record.precursorMz.size());
		output.fragmentOffset[outIndex] = static_cast<uint64_t>(fragmentOffset);
		output.fragmentCount[outIndex] = static_cast<uint64_t>(record.fragmentMz.size());

		std::copy(record.precursorMz.begin(), record.precursorMz.end(),
				  output.precursorMz.begin() + static_cast<std::ptrdiff_t>(precursorOffset));
		std::copy(record.precursorIntensity.begin(), record.precursorIntensity.end(),
				  output.precursorIntensity.begin() + static_cast<std::ptrdiff_t>(precursorOffset));
		std::copy(record.fragmentMz.begin(), record.fragmentMz.end(),
				  output.fragmentMz.begin() + static_cast<std::ptrdiff_t>(fragmentOffset));
		std::copy(record.theoreticalIntensity.begin(), record.theoreticalIntensity.end(),
				  output.theoreticalIntensity.begin() + static_cast<std::ptrdiff_t>(fragmentOffset));
		std::copy(record.experimentalIntensity.begin(), record.experimentalIntensity.end(),
				  output.experimentalIntensity.begin() + static_cast<std::ptrdiff_t>(fragmentOffset));
		std::copy(record.ionKinds.begin(), record.ionKinds.end(),
				  output.ionKind.begin() + static_cast<std::ptrdiff_t>(fragmentOffset));
		std::copy(record.ionPositions.begin(), record.ionPositions.end(),
				  output.ionPosition.begin() + static_cast<std::ptrdiff_t>(fragmentOffset));
	}

	return output;
}

bool writeSpectraHdf5File(const fs::path &path,
						  const Hdf5OutputData &data,
						  const Hdf5OutputMetadata &metadata)
{
	try
	{
		const fs::path parent = path.parent_path();
		if (!parent.empty())
		{
			fs::create_directories(parent);
		}

		H5::H5File file(path.string(), H5F_ACC_TRUNC);
		writeSpectraAttributes(file, metadata);

		H5::Group recordsGroup = file.createGroup("records");
		H5::Group precursorGroup = file.createGroup("precursor");
		H5::Group fragmentsGroup = file.createGroup("fragments");

		writeStringDataset(recordsGroup, "psm_id", data.psmIds, true);
		writeStringDataset(recordsGroup, "retention", data.retentions, true);
		writeVectorDataset(recordsGroup, "charge", H5::PredType::NATIVE_INT, data.charges, true);
		writeStringDataset(recordsGroup, "peptide", data.peptides, true);
		writeStringDataset(recordsGroup, "proteins", data.proteins, true);

		writeVectorDataset(precursorGroup, "mz", H5::PredType::NATIVE_DOUBLE, data.precursorMz, true);
		writeVectorDataset(precursorGroup, "intensity", H5::PredType::NATIVE_DOUBLE, data.precursorIntensity, true);
		writeVectorDataset(precursorGroup, "offset", H5::PredType::NATIVE_UINT64, data.precursorOffset, true);
		writeVectorDataset(precursorGroup, "count", H5::PredType::NATIVE_UINT64, data.precursorCount, true);

		writeVectorDataset(fragmentsGroup, "mz", H5::PredType::NATIVE_DOUBLE, data.fragmentMz, true);
		writeVectorDataset(fragmentsGroup, "theoretical_intensity", H5::PredType::NATIVE_DOUBLE, data.theoreticalIntensity, true);
		writeVectorDataset(fragmentsGroup, "experimental_intensity", H5::PredType::NATIVE_DOUBLE, data.experimentalIntensity, true);
		writeVectorDataset(fragmentsGroup, "ion_kind", H5::PredType::NATIVE_CHAR, data.ionKind, true);
		writeVectorDataset(fragmentsGroup, "ion_position", H5::PredType::NATIVE_UINT64, data.ionPosition, true);
		writeVectorDataset(fragmentsGroup, "offset", H5::PredType::NATIVE_UINT64, data.fragmentOffset, true);
		writeVectorDataset(fragmentsGroup, "count", H5::PredType::NATIVE_UINT64, data.fragmentCount, true);

		file.close();
		return true;
	}
	catch (const H5::Exception &ex)
	{
		std::cerr << ex.getDetailMsg() << "\n";
		return false;
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return false;
	}
}

void appendSpectraIndexRecord(
	SpectrumOutputRecord &&source,
	std::vector<sipros::SpectraIndexRecordInput> &records)
{
	double retentionSeconds = 0.0;
	if (!parseDoubleField(source.retention, retentionSeconds))
		throw std::runtime_error(
			"Invalid retention time in SFI record " + source.psmId);
	sipros::SpectraIndexRecordInput record;
	record.psmId = std::move(source.psmId);
	record.peptide = std::move(source.peptide);
	record.proteins = requireProteinNames(
		source.proteins, "SFI record " + record.psmId);
	record.charge = source.charge;
	// PSM input retention is validated as numeric seconds by
	// requireRetentionSeconds(); SFI stores one canonical unit (minutes).
	record.retentionMinutes = retentionSeconds / 60.0;
	record.sipAbundancePct = source.sipAbundancePct;
	const size_t precursorCount = std::min(
		source.precursorMz.size(), source.precursorIntensity.size());
	record.precursors.reserve(precursorCount);
	for (size_t k = 0; k < precursorCount; ++k)
		record.precursors.push_back(
			{source.precursorMz[k], source.precursorIntensity[k]});
	const size_t fragmentCount = std::min({source.fragmentMz.size(),
		source.theoreticalIntensity.size(), source.experimentalIntensity.size(),
		source.ionKinds.size(), source.ionPositions.size()});
	record.fragments.reserve(fragmentCount);
	for (size_t k = 0; k < fragmentCount; ++k)
	{
		sipros::SpectraIndexFragmentPeakInput fragment;
		fragment.mz = source.fragmentMz[k];
		fragment.theoreticalIntensity =
			static_cast<float>(source.theoreticalIntensity[k]);
		fragment.experimentalIntensity =
			static_cast<float>(source.experimentalIntensity[k]);
		fragment.ionKind = static_cast<uint8_t>(source.ionKinds[k]);
		fragment.ionPosition = static_cast<uint32_t>(source.ionPositions[k]);
		record.fragments.push_back(fragment);
	}
	records.push_back(std::move(record));
}

bool writeGeneratedSpectraFile(
	const fs::path &path,
	std::vector<sipros::SpectraIndexRecordInput> &indexRecords,
	const Hdf5OutputMetadata &metadata,
	size_t &written,
	int threads,
	sipros::SpectraIndexBuildStats &buildStats)
{
	sipros::SpectraIndexMetadata indexMetadata;
	indexMetadata.chemistryProfileId = metadata.chemistryProfileId;
	indexMetadata.recordKind = metadata.recordKind;
	indexMetadata.targetSipAbundancePct = metadata.targetSipAbundancePct;
	indexMetadata.sipAtom = metadata.sipAtom;
	indexMetadata.sipIsotopeMassNumber = metadata.sipIsotopeMassNumber;
	indexMetadata.probabilityCutoff = metadata.probCutoff;
	indexMetadata.generationPpmTolerance = metadata.ppmTolerance;
	indexMetadata.minimumMatchedEnvelopes = metadata.minMatchedEnvelopes;
	indexMetadata.envelopeTopN = metadata.envelopeTopN;
	indexMetadata.label = metadata.recordKind == "decoy" ? -1 : 1;
	std::string error;
	written = indexRecords.size();
	if (!sipros::SpectraIndex::write(
			path.string(), indexMetadata, indexRecords, error, threads, &buildStats,
			true))
	{
		std::cerr << error << "\n";
		return false;
	}
	return true;
}

void addOutputJobStats(ProcessingStats &processingStats, const OutputJobStats &jobStats)
{
	processingStats.targetFailed.fetch_add(static_cast<size_t>(jobStats.targetFailed), std::memory_order_relaxed);
	processingStats.decoyComputeFailed.fetch_add(static_cast<size_t>(jobStats.decoyComputeFailed), std::memory_order_relaxed);
}

bool generateAndWriteOutputFileJob(OutputFileJob &job,
								   const std::vector<PsmRow> &rows,
								   const std::vector<char> &baselineOk,
								   const std::vector<MatchedEnvelopeSet> &matchedEnvelopeSets,
								   const std::vector<std::string> &decoyPeptides,
								   const std::vector<char> &decoyAddedResidues,
								   const Isotopologue &pristineIso,
								   const Args &args,
								   char sipAtom,
								   int targetSipIsotopeIndex,
								   int workerThreads,
								   OutputJobStats &jobStats)
{
	try
	{
		const size_t abundanceCount = job.targetAbundancesPct.size();
		const std::vector<double> &targetAbundances =
			job.targetAbundancesPct;

		// Abundance-only scheduling leaves nearly every core idle for a fixed
		// abundance. Split each abundance into deterministic row chunks so both
		// fixed and ranged generation use all threads without sharing mutable
		// Isotopologue state or output vectors.
		constexpr size_t GenerationRowsPerTask = 256;
		const size_t chunksPerAbundance =
			(rows.size() + GenerationRowsPerTask - 1) /
			GenerationRowsPerTask;
		const size_t taskCount = abundanceCount * chunksPerAbundance;
		std::vector<std::vector<sipros::SpectraIndexRecordInput>>
			recordsByTask(taskCount);
		std::vector<OutputJobStats> statsByTask(taskCount);
		std::vector<std::string> errorsByTask(taskCount);
		std::atomic<size_t> nextTask{0};
		const size_t parallelism = std::min(
			taskCount,
			static_cast<size_t>(std::max(1, workerThreads)));
		const double generationStart = omp_get_wtime();
		std::vector<std::thread> workers;
		workers.reserve(parallelism);
		for (size_t worker = 0; worker < parallelism; ++worker)
		{
			workers.emplace_back([&]()
			{
				Isotopologue localIso = pristineIso;
				size_t activeAbundance = abundanceCount;
				while (true)
				{
					const size_t taskIndex =
						nextTask.fetch_add(1, std::memory_order_relaxed);
					if (taskIndex >= taskCount)
						break;
					const size_t abundanceIndex =
						taskIndex / chunksPerAbundance;
					const size_t chunkIndex = taskIndex % chunksPerAbundance;
					const size_t rowBegin = chunkIndex * GenerationRowsPerTask;
					const size_t rowEnd = std::min(
						rows.size(), rowBegin + GenerationRowsPerTask);
					auto &records = recordsByTask[taskIndex];
					auto &localStats = statsByTask[taskIndex];
					records.reserve(rowEnd - rowBegin);
					try
					{
						if (activeAbundance != abundanceIndex)
						{
							localIso = pristineIso;
							ProNovoConfig::setSipAbundance(
								localIso, sipAtom, targetSipIsotopeIndex,
								targetAbundances[abundanceIndex]);
							activeAbundance = abundanceIndex;
						}
						const double abundancePct =
							targetAbundances[abundanceIndex];
						for (size_t rowIndex = rowBegin;
							 rowIndex < rowEnd; ++rowIndex)
						{
							if (!baselineOk[rowIndex])
								continue;
							const PsmRow &row = rows[rowIndex];
							if (!job.decoy)
							{
								std::vector<std::vector<double>> yMass, yProb, bMass, bProb;
								if (!localIso.computeProductIon(
										row.peptide, yMass, yProb, bMass, bProb))
								{
									++localStats.targetFailed;
									continue;
								}
								IsotopeDistribution precursorDist;
								if (!buildPrecursorDistributionFromProductIons(
										localIso, yMass, yProb, bMass, bProb,
										precursorDist))
								{
									++localStats.targetFailed;
									continue;
								}
								SpectrumOutputRecord record;
								record.psmId = spectraHeaderId(
									row.psmId, abundancePct, false);
								record.retention = row.retentionText;
								record.charge = row.precursorCharge;
								record.peptide = row.peptide;
								record.proteins = row.proteins;
								record.sipAbundancePct = abundancePct;
								buildPrecursorChargePeaks(
									precursorDist, row.precursorCharge,
									args.probCutoff, args.envelopeTopN,
									record.precursorMz, record.precursorIntensity);
								std::vector<FragmentEntry> entries;
								if (!buildShiftedChargeOneFragmentEntries(
										bMass, bProb, yMass, yProb,
										matchedEnvelopeSets[rowIndex], args.probCutoff,
										args.envelopeTopN, entries))
								{
									++localStats.targetFailed;
									continue;
								}
								copyFragmentEntriesToRecord(entries, record);
								appendSpectraIndexRecord(std::move(record), records);
								continue;
							}

							if (decoyPeptides[rowIndex].empty())
								continue;
							const std::vector<double> apexPool =
								shuffledApexIntensityPool(
									matchedEnvelopeSets[rowIndex],
									args.decoySeed, rowIndex);
							const double meanApex = meanIntensity(apexPool);
							if (apexPool.empty() || meanApex <= 0.0)
							{
								++localStats.decoyComputeFailed;
								continue;
							}
							std::vector<std::vector<double>> decoyYMass, decoyYProb;
							std::vector<std::vector<double>> decoyBMass, decoyBProb;
							if (!localIso.computeProductIon(
									decoyPeptides[rowIndex], decoyYMass, decoyYProb,
									decoyBMass, decoyBProb))
							{
								++localStats.decoyComputeFailed;
								continue;
							}
							IsotopeDistribution decoyPrecursorDist;
							if (!buildPrecursorDistributionFromProductIons(
									localIso, decoyYMass, decoyYProb,
									decoyBMass, decoyBProb, decoyPrecursorDist))
							{
								++localStats.decoyComputeFailed;
								continue;
							}
							SpectrumOutputRecord decoyRecord;
							decoyRecord.psmId = spectraHeaderId(
								row.psmId, abundancePct, true);
							decoyRecord.retention = adjustedDecoyRetention(
								row.retentionText, args.decoySeed, rowIndex,
								decoyAddedResidues[rowIndex] != 0);
							decoyRecord.charge = row.precursorCharge;
							decoyRecord.peptide = decoyPeptides[rowIndex];
							decoyRecord.proteins = decoyProteinNames(row.proteins);
							decoyRecord.sipAbundancePct = abundancePct;
							buildPrecursorChargePeaks(
								decoyPrecursorDist, row.precursorCharge,
								args.probCutoff, args.envelopeTopN,
								decoyRecord.precursorMz,
								decoyRecord.precursorIntensity);
							std::vector<FragmentEntry> decoyEntries;
							if (!buildDecoyChargeOneFragmentEntries(
									decoyBMass, decoyBProb, decoyYMass, decoyYProb,
									apexPool, meanApex,
									decoyAddedResidues[rowIndex] != 0,
									args.probCutoff, args.envelopeTopN, decoyEntries))
							{
								++localStats.decoyComputeFailed;
								continue;
							}
							copyFragmentEntriesToRecord(decoyEntries, decoyRecord);
							appendSpectraIndexRecord(std::move(decoyRecord), records);
						}
					}
					catch (const std::exception &ex)
					{
						errorsByTask[taskIndex] = ex.what();
					}
					catch (...)
					{
						errorsByTask[taskIndex] = "unknown generation error";
					}
				}
			});
		}
		for (std::thread &worker : workers)
			worker.join();
		job.generationSeconds = omp_get_wtime() - generationStart;

		const double combineStart = omp_get_wtime();
		size_t totalRecordCount = 0;
		for (size_t taskIndex = 0; taskIndex < taskCount; ++taskIndex)
		{
			if (!errorsByTask[taskIndex].empty())
			{
				const size_t abundanceIndex = taskIndex / chunksPerAbundance;
				throw std::runtime_error(
					"SIP abundance " +
					std::to_string(job.targetAbundancesPct[abundanceIndex]) +
					" failed: " + errorsByTask[taskIndex]);
			}
			jobStats.targetFailed += statsByTask[taskIndex].targetFailed;
			jobStats.decoyComputeFailed +=
				statsByTask[taskIndex].decoyComputeFailed;
			totalRecordCount += recordsByTask[taskIndex].size();
		}

		std::vector<sipros::SpectraIndexRecordInput> combinedRecords;
		combinedRecords.reserve(totalRecordCount);
		for (size_t taskIndex = 0; taskIndex < taskCount; ++taskIndex)
		{
			auto &records = recordsByTask[taskIndex];
			combinedRecords.insert(combinedRecords.end(),
				std::make_move_iterator(records.begin()),
				std::make_move_iterator(records.end()));
			std::vector<sipros::SpectraIndexRecordInput>().swap(records);
		}
		job.combineSeconds = omp_get_wtime() - combineStart;

		job.metadata.targetSipAbundancePct =
			job.targetAbundancesPct.size() == 1
				? job.targetAbundancesPct.front()
				: sipros::SpectraIndexMetadata::MixedSipAbundancePct;
		job.success = writeGeneratedSpectraFile(
			job.path, combinedRecords, job.metadata, job.written,
			workerThreads, job.buildStats);
		return job.success;
	}
	catch (const std::exception &ex)
	{
		job.error = ex.what();
		job.success = false;
		return false;
	}
}

#if !defined(_WIN32)
struct OutputJobChildMessage
{
	uint64_t written = 0;
	uint64_t targetFailed = 0;
	uint64_t decoyComputeFailed = 0;
	sipros::SpectraIndexBuildStats buildStats;
	double generationSeconds = 0.0;
	double combineSeconds = 0.0;
	uint8_t success = 0;
	char error[512] = {0};
};

struct ActiveOutputChild
{
	pid_t pid = -1;
	int readFd = -1;
	size_t jobIndex = 0;
};

bool writeAllToFd(int fd, const void *buffer, size_t size)
{
	const char *cursor = static_cast<const char *>(buffer);
	size_t remaining = size;
	while (remaining > 0)
	{
		const ssize_t written = ::write(fd, cursor, remaining);
		if (written < 0)
		{
			if (errno == EINTR)
			{
				continue;
			}
			return false;
		}
		if (written == 0)
		{
			return false;
		}
		cursor += written;
		remaining -= static_cast<size_t>(written);
	}
	return true;
}

bool readAllFromFd(int fd, void *buffer, size_t size)
{
	char *cursor = static_cast<char *>(buffer);
	size_t remaining = size;
	while (remaining > 0)
	{
		const ssize_t bytesRead = ::read(fd, cursor, remaining);
		if (bytesRead < 0)
		{
			if (errno == EINTR)
			{
				continue;
			}
			return false;
		}
		if (bytesRead == 0)
		{
			return false;
		}
		cursor += bytesRead;
		remaining -= static_cast<size_t>(bytesRead);
	}
	return true;
}

void closeFdIfOpen(int &fd)
{
	if (fd >= 0)
	{
		::close(fd);
		fd = -1;
	}
}

OutputJobChildMessage makeOutputJobChildMessage(const OutputFileJob &job,
												const OutputJobStats &jobStats)
{
	OutputJobChildMessage message;
	message.written = static_cast<uint64_t>(job.written);
	message.targetFailed = jobStats.targetFailed;
	message.decoyComputeFailed = jobStats.decoyComputeFailed;
	message.buildStats = job.buildStats;
	message.generationSeconds = job.generationSeconds;
	message.combineSeconds = job.combineSeconds;
	message.success = job.success ? 1 : 0;
	const std::string error = job.error.empty() && !job.success ? "child writer failed" : job.error;
	if (!error.empty())
	{
		std::strncpy(message.error, error.c_str(), sizeof(message.error) - 1);
		message.error[sizeof(message.error) - 1] = '\0';
	}
	return message;
}

void applyOutputJobChildMessage(OutputFileJob &job,
								const OutputJobChildMessage &message,
								int childStatus,
								bool messageRead,
								ProcessingStats &processingStats)
{
	if (messageRead)
	{
		job.written = static_cast<size_t>(message.written);
		job.success = message.success != 0 && WIFEXITED(childStatus) && WEXITSTATUS(childStatus) == 0;
		job.error = message.error;
		job.buildStats = message.buildStats;
		job.generationSeconds = message.generationSeconds;
		job.combineSeconds = message.combineSeconds;
		OutputJobStats jobStats;
		jobStats.targetFailed = message.targetFailed;
		jobStats.decoyComputeFailed = message.decoyComputeFailed;
		addOutputJobStats(processingStats, jobStats);
	}
	else
	{
		job.success = false;
		job.error = "failed to read child writer result";
	}

	if (!WIFEXITED(childStatus))
	{
		job.success = false;
		job.error = "child writer process did not exit normally";
	}
	else if (WEXITSTATUS(childStatus) != 0 && job.error.empty())
	{
		job.success = false;
		job.error = "child writer process exited with status " + std::to_string(WEXITSTATUS(childStatus));
	}
}

bool waitForOneOutputChild(std::vector<ActiveOutputChild> &activeChildren,
						   std::vector<OutputFileJob> &outputJobs,
						   ProcessingStats &processingStats)
{
	while (true)
	{
		int status = 0;
		const pid_t pid = ::waitpid(-1, &status, 0);
		if (pid < 0)
		{
			if (errno == EINTR)
			{
				continue;
			}
			std::cerr << "waitpid failed while waiting for output writer process: "
					  << std::strerror(errno) << "\n";
			return false;
		}

		const auto childIt = std::find_if(activeChildren.begin(), activeChildren.end(),
										  [pid](const ActiveOutputChild &child)
										  { return child.pid == pid; });
		if (childIt == activeChildren.end())
		{
			continue;
		}

		OutputJobChildMessage message;
		const bool messageRead = readAllFromFd(childIt->readFd, &message, sizeof(message));
		closeFdIfOpen(childIt->readFd);
		applyOutputJobChildMessage(outputJobs[childIt->jobIndex], message, status,
								   messageRead, processingStats);
		activeChildren.erase(childIt);
		return true;
	}
}

bool runOutputJobsWithFork(std::vector<OutputFileJob> &outputJobs,
						   const std::vector<PsmRow> &rows,
						   const std::vector<char> &baselineOk,
						   const std::vector<MatchedEnvelopeSet> &matchedEnvelopeSets,
						   const std::vector<std::string> &decoyPeptides,
						   const std::vector<char> &decoyAddedResidues,
						   const Isotopologue &pristineIso,
						   const Args &args,
						   char sipAtom,
						   int targetSipIsotopeIndex,
						   int effectiveThreads,
						   ProcessingStats &processingStats)
{
	const size_t processLimit = static_cast<size_t>(std::max(1, effectiveThreads));
	std::vector<ActiveOutputChild> activeChildren;
	activeChildren.reserve(processLimit);

	for (size_t jobIndex = 0; jobIndex < outputJobs.size(); ++jobIndex)
	{
		while (activeChildren.size() >= processLimit)
		{
			if (!waitForOneOutputChild(activeChildren, outputJobs, processingStats))
			{
				return false;
			}
		}

		int pipeFd[2] = {-1, -1};
		if (::pipe(pipeFd) != 0)
		{
			std::cerr << "Cannot create output writer pipe: " << std::strerror(errno) << "\n";
			return false;
		}

		std::cout.flush();
		std::cerr.flush();
		const pid_t pid = ::fork();
		if (pid < 0)
		{
			closeFdIfOpen(pipeFd[0]);
			closeFdIfOpen(pipeFd[1]);
			std::cerr << "Cannot fork output writer process: " << std::strerror(errno) << "\n";
			return false;
		}

		if (pid == 0)
		{
			closeFdIfOpen(pipeFd[0]);
			omp_set_num_threads(1);
			OutputJobStats childStats;
			const int threadsPerJob = std::max(
				1, effectiveThreads /
					static_cast<int>(std::min(
						outputJobs.size(), static_cast<size_t>(effectiveThreads))));
			generateAndWriteOutputFileJob(outputJobs[jobIndex], rows, baselineOk, matchedEnvelopeSets,
										  decoyPeptides, decoyAddedResidues, pristineIso,
										  args, sipAtom, targetSipIsotopeIndex,
										  threadsPerJob, childStats);
			const OutputJobChildMessage message =
				makeOutputJobChildMessage(outputJobs[jobIndex], childStats);
			const bool wroteMessage = writeAllToFd(pipeFd[1], &message, sizeof(message));
			closeFdIfOpen(pipeFd[1]);
			std::_Exit(outputJobs[jobIndex].success && wroteMessage ? 0 : 1);
		}

		closeFdIfOpen(pipeFd[1]);
		activeChildren.push_back({pid, pipeFd[0], jobIndex});
	}

	while (!activeChildren.empty())
	{
		if (!waitForOneOutputChild(activeChildren, outputJobs, processingStats))
		{
			return false;
		}
	}
	return true;
}
#endif

bool runOutputJobs(std::vector<OutputFileJob> &outputJobs,
				   const std::vector<PsmRow> &rows,
				   const std::vector<char> &baselineOk,
				   const std::vector<MatchedEnvelopeSet> &matchedEnvelopeSets,
				   const std::vector<std::string> &decoyPeptides,
				   const std::vector<char> &decoyAddedResidues,
				   const Isotopologue &pristineIso,
				   const Args &args,
				   char sipAtom,
				   int targetSipIsotopeIndex,
				   int effectiveThreads,
				   ProcessingStats &processingStats)
{
#if !defined(_WIN32)
	return runOutputJobsWithFork(outputJobs, rows, baselineOk, matchedEnvelopeSets,
								 decoyPeptides, decoyAddedResidues, pristineIso, args,
								 sipAtom, targetSipIsotopeIndex, effectiveThreads, processingStats);
#else
	std::vector<OutputJobStats> jobStats(outputJobs.size());
	for (size_t i = 0; i < outputJobs.size(); ++i)
	{
		generateAndWriteOutputFileJob(outputJobs[i], rows, baselineOk,
									  matchedEnvelopeSets, decoyPeptides, decoyAddedResidues,
									  pristineIso, args, sipAtom, targetSipIsotopeIndex,
									  effectiveThreads, jobStats[i]);
	}
	for (const OutputJobStats &stats : jobStats)
	{
		addOutputJobStats(processingStats, stats);
	}
	return true;
#endif
}

} // namespace

int ExperimentalSpectraWorkflow::run(int argc, char **argv)
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

	std::setlocale(LC_ALL, "C");
	std::ios_base::sync_with_stdio(false);

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

	int targetSipIsotopeIndex = 1;
	try
	{
		targetSipIsotopeIndex = ProNovoConfig::resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, sipAtom, sipIsotopeMassNumber);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	int baselineSipIsotopeIndex = 1;
	try
	{
		baselineSipIsotopeIndex = ProNovoConfig::resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, 'C', 13);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Failed to resolve baseline C13 isotope: " << ex.what() << "\n";
		return 1;
	}

	Isotopologue pristineIso = ProNovoConfig::configIsotopologue;
	std::string chemistryError;
	if (!ProNovoConfig::validatePreparationChemistry(
			pristineIso, chemistryError))
	{
		std::cerr << chemistryError << "\n";
		return 1;
	}
	Isotopologue baselineIso = pristineIso;
	try
	{
		const double baselineC13Pct = ProNovoConfig::getIsotopeAbundancePct(
			pristineIso, 'C', baselineSipIsotopeIndex);
		ProNovoConfig::setSipAbundance(baselineIso, 'C', baselineSipIsotopeIndex, baselineC13Pct);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Failed to apply baseline natural C13 abundance: " << ex.what() << "\n";
		return 1;
	}

	const std::vector<double> targetAbundances = makeTargetAbundances(args);

	if (args.threads > 0)
	{
		omp_set_num_threads(args.threads);
	}
	const int effectiveThreads = std::max(1, omp_get_max_threads());
	const double processWallStart = omp_get_wtime();
	const double processCpuStart = processTreeCpuSeconds();
	TimingLogger timing;
	timing.printHeader();

	ReadStats readStats;
	std::vector<PsmRow> allRows;
	try
	{
		timing.run("Read PSM", "read PSM rows", 0, "", [&]()
		{
			allRows = readInputRows(args.inputPath, readStats);
		});
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	std::unordered_set<std::string> experimentalPeptideMassClasses;
	experimentalPeptideMassClasses.reserve(allRows.size());
	for (const PsmRow &row : allRows)
	{
		experimentalPeptideMassClasses.insert(peptideMassClassKey(row.peptide));
	}

	std::vector<PsmRow> rows;
	timing.run("Select PSM", "PEP/SVM/WDP peptide/charge", allRows.size(), "rows", [&]()
	{
		rows = selectBestRowsByPeptideCharge(allRows);
	});
	std::cout << "PSM: " << args.inputPath << "  ("
			  << readStats.totalRows << " rows; " << rows.size()
			  << " retained peptide/charge rows)\n";
	if (readStats.shortRows > 0 || readStats.invalidRows > 0 ||
		readStats.unsupportedMods > 0 || readStats.decoyRows > 0)
	{
		std::cerr << "PSM skipped rows: short=" << readStats.shortRows
				  << ", invalid=" << readStats.invalidRows
				  << ", unsupported_mods=" << readStats.unsupportedMods
				  << ", decoy=" << readStats.decoyRows << "\n";
	}

	std::unordered_map<std::string, fs::path> hdf5FilesBySample;
	try
	{
		timing.run("Collect HDF5", "collect HDF5 files", 0, "", [&]()
		{
			hdf5FilesBySample = collectHdf5Files(args.hdf5Path);
		});
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}
	const std::map<std::string, size_t> requestedScanCounts = requestedScanCountsBySample(rows);
	std::vector<std::pair<std::string, fs::path>> sortedHdf5Files(hdf5FilesBySample.begin(), hdf5FilesBySample.end());
	std::sort(sortedHdf5Files.begin(), sortedHdf5Files.end(),
			  [](const auto &a, const auto &b)
			  { return a.second.string() < b.second.string(); });
	if (sortedHdf5Files.size() == 1)
	{
		const auto &entry = sortedHdf5Files.front();
		const auto countIt = requestedScanCounts.find(entry.first);
		const size_t requestedScans = countIt == requestedScanCounts.end() ? 0 : countIt->second;
		std::cout << "HDF5: " << entry.second << "  ("
				  << requestedScans << " requested scans)\n";
	}
	else
	{
		std::cout << "HDF5: " << args.hdf5Path << "  ("
				  << sortedHdf5Files.size() << " files)\n";
		for (const auto &entry : sortedHdf5Files)
		{
			const auto countIt = requestedScanCounts.find(entry.first);
			const size_t requestedScans = countIt == requestedScanCounts.end() ? 0 : countIt->second;
			std::cout << "  HDF5: " << entry.second << "  ("
					  << requestedScans << " requested scans)\n";
		}
	}
	std::cout << (args.writeDecoy ? 2 : 1)
			  << " SFI output files planned\n";

	std::vector<MatchedEnvelopeSet> matchedEnvelopeSets(rows.size());
	std::vector<char> baselineOk(rows.size(), 0);
	ProcessingStats processingStats;

	try
	{
		matchBaselineInHdf5Batches(rows, hdf5FilesBySample, baselineIso, args, effectiveThreads,
								  matchedEnvelopeSets, baselineOk, processingStats, timing);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	size_t baselineMatchedRows = 0;
	for (char value : baselineOk)
	{
		if (value)
		{
			++baselineMatchedRows;
		}
	}

	std::vector<std::string> decoyPeptides(rows.size());
	std::vector<char> decoyAddedResidues(rows.size(), 0);
	size_t decoyPeptidesGenerated = 0;
	if (args.writeDecoy)
	{
		size_t generatedCount = 0;
		{
			timing.run("Generate Decoy", "generate decoy peptides", rows.size(), "rows", [&]()
			{
				size_t generatedInRun = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : generatedInRun)
				for (int i = 0; i < static_cast<int>(rows.size()); ++i)
				{
					const size_t rowIndex = static_cast<size_t>(i);
					if (!baselineOk[rowIndex])
					{
						continue;
					}
					bool addedResidue = false;
					if (generateDecoyPeptide(rows[rowIndex].peptide, experimentalPeptideMassClasses,
										 args.decoySeed, rowIndex, decoyPeptides[rowIndex],
											 addedResidue))
					{
						decoyAddedResidues[rowIndex] = addedResidue ? 1 : 0;
						++generatedInRun;
					}
					else
					{
						++processingStats.decoyCollision;
					}
				}
				generatedCount += generatedInRun;
			});
		}
		decoyPeptidesGenerated = generatedCount;
	}

	const size_t plannedOutputFileCount = args.writeDecoy ? 2 : 1;
	const bool multipleOutputFiles = plannedOutputFileCount > 1;
	fs::path outputBasePath;
	try
	{
		outputBasePath = resolveOutputBasePath(args.outputPath, multipleOutputFiles);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}
	size_t outputFilesWritten = 0;
	size_t decoyFilesWritten = 0;
	std::vector<OutputFileJob> outputJobs;
	outputJobs.reserve(plannedOutputFileCount);
	Hdf5OutputMetadata outputMetadata;
	outputMetadata.recordKind = "target";
	outputMetadata.chemistryProfileId = ProNovoConfig::getChemistryProfileId();
	outputMetadata.targetSipAbundancePct = targetAbundances.size() == 1
		? targetAbundances.front()
		: sipros::SpectraIndexMetadata::MixedSipAbundancePct;
	outputMetadata.sipAtom = sipAtom;
	outputMetadata.sipIsotopeMassNumber = sipIsotopeMassNumber;
	outputMetadata.probCutoff = args.probCutoff;
	outputMetadata.ppmTolerance = args.ppmTolerance;
	outputMetadata.minMatchedEnvelopes = static_cast<uint64_t>(args.minMatchedEnvelopes);
	outputMetadata.envelopeTopN = static_cast<uint32_t>(args.envelopeTopN);

	OutputFileJob targetJob;
	targetJob.path = outputBasePath;
	targetJob.metadata = outputMetadata;
	targetJob.decoy = false;
	targetJob.targetAbundancesPct = targetAbundances;
	outputJobs.push_back(std::move(targetJob));

	if (args.writeDecoy)
	{
		Hdf5OutputMetadata decoyMetadata = outputMetadata;
		decoyMetadata.recordKind = "decoy";
		OutputFileJob decoyJob;
		decoyJob.path = decoyOutputPath(outputBasePath);
		decoyJob.metadata = decoyMetadata;
		decoyJob.decoy = true;
		decoyJob.targetAbundancesPct = targetAbundances;
		outputJobs.push_back(std::move(decoyJob));
	}

	{
		bool outputOk = false;
		timing.run("Generate/write SFI", "generate spectra and write SFI files",
			outputJobs.size(), "files", [&]()
		{
			outputOk = runOutputJobs(outputJobs, rows, baselineOk, matchedEnvelopeSets,
									 decoyPeptides, decoyAddedResidues, pristineIso, args,
									 sipAtom, targetSipIsotopeIndex, effectiveThreads, processingStats);
		});
		if (!outputOk)
		{
			return 1;
		}
	}

	bool allOutputSucceeded = true;
	for (const OutputFileJob &job : outputJobs)
	{
		if (!job.success)
		{
			if (!job.error.empty())
			{
				std::cerr << job.error << "\n";
			}
			std::cerr << "Cannot write output file: " << job.path << "\n";
			allOutputSucceeded = false;
			continue;
		}
		const auto &build = job.buildStats;
		std::ostringstream buildLog;
		buildLog << "SFI v5 compact RT-aware index (top "
				 << job.metadata.envelopeTopN << " peaks/envelope): "
				 << job.path << '\n'
				 << "  records=" << build.recordCount
				 << ", precursors=" << build.precursorCount
				 << ", fragments=" << build.fragmentCount << " x 16 bytes"
				 << ", packed_product_postings=" << build.productPostingCount
				 << " x 4 bytes, sparse_rt_bins=" << build.rtBinCount
				 << ", blocks=" << build.blockCount << '\n'
					 << "  parallel spectra generation: " << std::fixed
				 << std::setprecision(3) << job.generationSeconds << "s (up to "
				 << std::max(1, effectiveThreads /
					static_cast<int>(outputJobs.size())) << " threads/job)\n"
					 << "  combine generated records: " << job.combineSeconds << "s\n"
					 << "  parallel compact/validate: "
					 << build.compactValidateSeconds
					 << "s (generation supplied compact envelopes; "
					 << build.threadsUsed << " threads)\n"
					 << "  precursor ordering: " << build.orderSeconds << "s; reserve="
					 << build.reserveSeconds << "s\n"
				 << "  parallel product-index build: " << std::fixed
				 << std::setprecision(3) << build.productIndexSeconds << "s ("
				 << build.threadsUsed << " threads)\n"
				 << "  parallel flatten=" << build.flattenSeconds << "s, layout/checksum="
				 << build.layoutChecksumSeconds << "s, parallel write="
				 << build.writeSeconds << "s (" << build.threadsUsed
				 << " threads), total=" << build.totalSeconds << "s, file="
				 << std::setprecision(3)
				 << static_cast<double>(build.fileBytes) / (1024.0 * 1024.0 * 1024.0)
				 << " GiB\n";
		std::cerr << buildLog.str();

		if (job.decoy)
		{
			++decoyFilesWritten;
			std::cerr << "Wrote decoy shifted spectra for " << job.written
					  << " abundance-specific PSM records across "
					  << job.targetAbundancesPct.size() << " SIP percentages to "
					  << job.path << "\n";
		}
		else
		{
			++outputFilesWritten;
			std::cerr << "Wrote shifted matched experimental spectra for " << job.written
					  << " abundance-specific PSM records across "
					  << job.targetAbundancesPct.size() << " SIP percentages to "
					  << job.path << "\n";
		}
	}

	if (!allOutputSucceeded)
	{
		return 1;
	}

	std::cerr << "Input peptide spectra: " << rows.size()
			  << "; baseline matched peptide spectra extracted: " << baselineMatchedRows
			  << "; output files written: " << outputFilesWritten
			  << "; minimum matched envelopes: " << args.minMatchedEnvelopes
			  << "; SFI envelope top-N: " << args.envelopeTopN;
	if (args.writeDecoy)
	{
		std::cerr << "; decoy peptides generated: " << decoyPeptidesGenerated
				  << "; decoy files written: " << decoyFilesWritten;
	}
	std::cerr << "\n";
	std::cerr << "Skipped during processing: missing_hdf5=" << processingStats.missingHdf5.load()
			  << ", missing_scan=" << processingStats.missingScan.load()
			  << ", compute_failed=" << processingStats.computeFailed.load()
			  << ", precursor_failed=" << processingStats.precursorFailed.load()
			  << ", insufficient_matched_envelopes=" << processingStats.unmatched.load()
			  << ", target_failed=" << processingStats.targetFailed.load();
	if (args.writeDecoy)
	{
		std::cerr << ", decoy_collision=" << processingStats.decoyCollision.load()
				  << ", decoy_compute_failed=" << processingStats.decoyComputeFailed.load();
	}
	std::cerr << "\n";
	const double processWallSeconds = omp_get_wtime() - processWallStart;
	const double processCpuSeconds = processTreeCpuSeconds() - processCpuStart;
	timing.printSummary(processWallSeconds, processCpuSeconds,
						readStats.totalRows, rows.size(), baselineMatchedRows,
						outputFilesWritten + decoyFilesWritten);
	return 0;
}
