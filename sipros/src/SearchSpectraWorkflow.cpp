// sipros search-spectra
//
// Re-score Raxport HDF5 MS2 scans against memory-mapped SIP spectra indexes
// (.sfi). Candidate discovery uses isolation-window/RT ranges and sparse
// product-ion postings before the existing MVH, Xcorr, and WDP scoring stages.
// MS1 evidence is attached only after scoring retained PSMs.

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cctype>
#include <filesystem>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <tuple>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#include <omp.h>

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
#include "spectraindex.h"
#include "performancelog.h"

namespace sipros
{
#ifdef _MSC_VER
unsigned populationCount64(uint64_t value)
{
	return static_cast<unsigned>(__popcnt64(value));
}

unsigned trailingZeroCount64(uint64_t value)
{
	unsigned long index = 0;
	_BitScanForward64(&index, value);
	return static_cast<unsigned>(index);
}
#else
unsigned populationCount64(uint64_t value)
{
	return static_cast<unsigned>(__builtin_popcountll(value));
}

unsigned trailingZeroCount64(uint64_t value)
{
	return static_cast<unsigned>(__builtin_ctzll(value));
}
#endif
// First-gate settings validated for the regular FASTA DDA search.
constexpr uint32_t kDdaWindowMinMatchedFragments = 4;
constexpr size_t kDdaWindowTopPeaks = 200;
constexpr int kPrecursorMatchScanRadius = 5;

// -------------------- Args --------------------

struct Args
{
	std::string workingDir;
	std::string singleHdf5;
	std::string sfiDir;
	std::string outputDir;
	int sfiLabel = 0;
	int threads = 0;
	double rtToleranceMin = 5.0;
	// Separate MS1 (precursor) and MS2 (fragment) tolerances. Default 10 ppm
	// each. --tolerance / --tolerance-unit are shortcuts that set both.
	double toleranceMs1 = 10.0;
	bool toleranceMs1Ppm = true;
	double toleranceMs2 = 10.0;
	bool toleranceMs2Ppm = true;
	int mvhCascadeTopN = 150;
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

struct TimingEntry
{
	std::string group;
	std::string detailLabel;
	std::string workName;
	sipros::PerformanceTiming timing;
	size_t workItems = 0;
};

struct TimingLogger
{
	std::vector<TimingEntry> entries;

	void printHeader(const std::string &title, int threadCount) const
	{
		sipros::printPerformanceHeader(std::cout, title, threadCount);
	}

	template <typename Fn>
	void run(const std::string &group,
			 const std::string &label,
			 size_t workItems,
			 const std::string &workName,
			 Fn &&fn)
	{
		const sipros::PerformanceTimer timer;
		fn();

		TimingEntry entry;
		entry.group = group;
		entry.detailLabel = label;
		entry.workName = workName;
		entry.timing = timer.elapsed();
		entry.workItems = workItems;
		entries.push_back(entry);
		printEntry(entry);
	}

	static void printEntry(const TimingEntry &entry)
	{
		std::string detail = entry.detailLabel;
		if (entry.workItems > 0 && !entry.workName.empty())
		{
			if (!detail.empty())
				detail += "; ";
			detail += sipros::formatPerformanceCount(entry.workItems) +
				" " + entry.workName;
		}
		else if (!entry.workName.empty())
		{
			if (!detail.empty())
				detail += "; ";
			detail += entry.workName;
		}
		sipros::printPerformanceStage(
			std::cout, entry.group, entry.timing, detail);
	}

	void printFooter() const
	{
		sipros::printPerformanceFooter(std::cout);
	}
};

void printUsage(const char *prog)
{
	std::cerr
		<< "Usage:\n  " << prog
		<< " -w <Raxport HDF5 dir> [-f <single.h5>]\n"
		<< "    --sfi <SIP spectra-index dir> -o <PIN output dir> [-t <N>] [--rt-tolerance <min>]\n"
		<< "    --sfi-label target|decoy                             search one label set and write _target/_decoy.pin\n"
		<< "    [--tolerance-ms1 <N>] [--tolerance-ms1-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance-ms2 <N>] [--tolerance-ms2-unit ppm|da]   (default: 10 ppm)\n"
		<< "    [--tolerance <N>] [--tolerance-unit ppm|da]            shortcut: set BOTH MS1 and MS2\n"
		<< "    [--mvh-cascade-top-n <N>]                              MVH candidates per scan sent to Xcorr/WDP (default: 150)\n"
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
		else if (k == "--sfi")
			a.sfiDir = next();
		else if (k == "--sfi-label")
		{
			std::string label = next();
			for (char &c : label)
				c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
			if (label == "target" || label == "1" || label == "+1")
				a.sfiLabel = 1;
			else if (label == "decoy" || label == "-1")
				a.sfiLabel = -1;
			else
			{
				std::cerr << "--sfi-label must be target or decoy\n";
				std::exit(1);
			}
		}
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
		else if (k == "--mvh-cascade-top-n")
			a.mvhCascadeTopN = std::atoi(next().c_str());
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
	if (a.sfiDir.empty())
	{
		std::cerr << "--sfi <SIP spectra-index dir> is required; HDF5 spectra libraries are not supported\n";
		std::exit(1);
	}
	if (a.sfiLabel == 0)
	{
		std::cerr << "--sfi-label target|decoy is required; paired workflow processes run concurrently\n";
		std::exit(1);
	}
	if (a.threads <= 0)
	{
		a.threads = std::max(1, omp_get_num_procs());
	}
	if (a.mvhCascadeTopN <= 0)
	{
		std::cerr << "--mvh-cascade-top-n must be > 0\n";
		std::exit(1);
	}
	if (a.topPsmsPerScan <= 0)
	{
		std::cerr << "--top-psms-per-scan must be > 0\n";
		std::exit(1);
	}
	return a;
}

std::vector<std::string> listFiles(const std::string &dir,
								   const std::vector<std::string> &exts)
{
	std::vector<std::string> out;
	if (!std::filesystem::is_directory(dir))
		return out;
	for (const auto &e : std::filesystem::directory_iterator(dir))
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

bool isDecoySfi(const std::string &path)
{
	std::string name = std::filesystem::path(path).filename().string();
	std::transform(name.begin(), name.end(), name.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return name.find("decoy") != std::string::npos;
}

// -------------------- SipRecord --------------------

struct SipRecord
{
	std::string peptide;     // bracketed/decorated form from the spectra index
	std::string nakedPeptide; // letters-only used for scoring
	std::string proteins;
	int charge = 1;
	double topPrecursorMz = 0.0;
	double sumPrecursorIntensity = 0.0;
	// indexed by ion position (1-based ⇒ pos-1 row); each row holds an envelope
	std::vector<std::vector<double>> vvdYionMass;
	std::vector<std::vector<double>> vvdYionProb;
	std::vector<std::vector<double>> vvdYionExpIntensity;
	std::vector<std::vector<double>> vvdBionMass;
	std::vector<std::vector<double>> vvdBionProb;
	std::vector<std::vector<double>> vvdBionExpIntensity;
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

// -------------------- Matching --------------------

struct PrecursorMatch
{
	bool ok = false;
	int isotopicShift = 0;
	double mzShift = 0.0;
	double mzAbsErrPpm = 0.0;
	double deltaRT = 0.0;        // abs(scan RT - SFI record RT) in minutes
	double precursorRtDiffSeconds = -1.0;
	int matchedCharge = 0;
	double matchedMz = 0.0;
	int matchedMs1ScanNumber = 0;
	double isolationCenterMz = 0.0;
};

double scanRtMinutes(const MS2Scan *scan)
{
	return parseRetentionMinutes(scan->sRTime);
}

bool ms1PeakChargeMatches(int nativeCharge, int precursorCharge)
{
	return nativeCharge > 0
		? nativeCharge == precursorCharge
		: nativeCharge == 0 && precursorCharge >= 2 && precursorCharge <= 4;
}

PrecursorMatch matchNearbyEnvelopePrecursor(
	const sipros::RaxportMs1Data &ms1Data,
	int parentScanNumber,
	double ms2RetentionTimeMinutes,
	int precursorCharge,
	const std::vector<double> &precursorNeutralMasses,
	const std::vector<double> &precursorProbabilities,
	const MassTolerance &mzTolerance,
	int scanRadius,
	const PrecursorMatch &candidateMatch)
{
	PrecursorMatch best = candidateMatch;
	best.precursorRtDiffSeconds = -1.0;
	best.matchedMs1ScanNumber = 0;
	const size_t envelopeSize = std::min(
		precursorNeutralMasses.size(), precursorProbabilities.size());
	if (scanRadius < 0 || precursorCharge <= 0 || envelopeSize == 0 ||
		ms1Data.scans.empty())
	{
		return best;
	}
	const auto parent = ms1Data.scanNumberToIndex.find(parentScanNumber);
	if (parent == ms1Data.scanNumberToIndex.end())
		return best;

	size_t modalIndex = envelopeSize;
	double maximumProbability = 0.0;
	for (size_t index = 0; index < envelopeSize; ++index)
	{
		if (std::isfinite(precursorNeutralMasses[index]) &&
			std::isfinite(precursorProbabilities[index]) &&
			precursorProbabilities[index] > maximumProbability)
		{
			modalIndex = index;
			maximumProbability = precursorProbabilities[index];
		}
	}
	if (modalIndex == envelopeSize)
		return best;
	const double modalNeutralMass = precursorNeutralMasses[modalIndex];
	const double probabilityFloor = maximumProbability * 1e-12;
	const double proton = ProNovoConfig::getProtonMass();
	const double neutron = ProNovoConfig::getNeutronMass();
	const int64_t first = std::max<int64_t>(
		0, static_cast<int64_t>(parent->second) - scanRadius);
	const int64_t last = std::min<int64_t>(
		static_cast<int64_t>(ms1Data.scans.size()) - 1,
		static_cast<int64_t>(parent->second) + scanRadius);
	const auto rangeBegin = ms1Data.scans.begin() + first;
	const auto rangeEnd = ms1Data.scans.begin() + last + 1;
	const auto insertion = std::lower_bound(
		rangeBegin, rangeEnd, ms2RetentionTimeMinutes,
		[](const sipros::RaxportMs1Scan &scan, double retentionTime)
		{ return scan.retentionTime < retentionTime; });
	int64_t right = static_cast<int64_t>(
		std::distance(ms1Data.scans.begin(), insertion));
	int64_t left = right - 1;
	bool found = false;
	double bestMassError = std::numeric_limits<double>::infinity();
	double bestIntensity = 0.0;
	while (left >= first || right <= last)
	{
		const double leftRtDiff = left >= first
			? std::abs(ms1Data.scans[static_cast<size_t>(left)].retentionTime -
				ms2RetentionTimeMinutes)
			: std::numeric_limits<double>::infinity();
		const double rightRtDiff = right <= last
			? std::abs(ms1Data.scans[static_cast<size_t>(right)].retentionTime -
				ms2RetentionTimeMinutes)
			: std::numeric_limits<double>::infinity();
		const int64_t scanIndex = leftRtDiff <= rightRtDiff ? left-- : right++;
		const sipros::RaxportMs1Scan &ms1Scan =
			ms1Data.scans[static_cast<size_t>(scanIndex)];
		const double rtDiffSeconds =
			std::abs(ms1Scan.retentionTime - ms2RetentionTimeMinutes) * 60.0;
		if (found && rtDiffSeconds > best.precursorRtDiffSeconds)
			break;

		for (size_t envelopeIndex = 0;
			 envelopeIndex < envelopeSize; ++envelopeIndex)
		{
			const double expectedNeutralMass =
				precursorNeutralMasses[envelopeIndex];
			if (precursorProbabilities[envelopeIndex] < probabilityFloor ||
				!std::isfinite(expectedNeutralMass))
			{
				continue;
			}
			const double expectedMz =
				expectedNeutralMass / precursorCharge + proton;
			const double toleranceMz = mzTolerance.daAt(expectedMz);
			auto peak = std::lower_bound(
				ms1Scan.mz.begin(), ms1Scan.mz.end(), expectedMz - toleranceMz);
			for (; peak != ms1Scan.mz.end() && *peak <= expectedMz + toleranceMz;
				 ++peak)
			{
				const size_t peakIndex = static_cast<size_t>(
					peak - ms1Scan.mz.begin());
				const int nativeCharge = peakIndex < ms1Scan.charge.size()
					? ms1Scan.charge[peakIndex] : 0;
				if (!ms1PeakChargeMatches(nativeCharge, precursorCharge))
					continue;
				const double massError =
					std::abs(*peak - expectedMz) * precursorCharge;
				const double intensity = peakIndex < ms1Scan.intensity.size()
					? ms1Scan.intensity[peakIndex] : 0.0;
				const bool better = !found ||
					rtDiffSeconds < best.precursorRtDiffSeconds ||
					(rtDiffSeconds == best.precursorRtDiffSeconds &&
					 (massError < bestMassError ||
					  (massError == bestMassError && intensity > bestIntensity)));
				if (!better)
					continue;
				found = true;
				bestMassError = massError;
				bestIntensity = intensity;
				best.ok = true;
				best.isotopicShift = static_cast<int>(std::lround(
					(expectedNeutralMass - modalNeutralMass) / neutron));
				best.mzAbsErrPpm = expectedMz > 0.0
					? std::abs(*peak - expectedMz) / expectedMz * 1.0e6 : 0.0;
				best.precursorRtDiffSeconds = rtDiffSeconds;
				best.matchedCharge = precursorCharge;
				best.matchedMz = *peak;
				best.matchedMs1ScanNumber = ms1Scan.scanNumber;
				if (candidateMatch.isolationCenterMz > 0.0)
					best.mzShift = std::abs(
						*peak - candidateMatch.isolationCenterMz);
			}
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
								const std::vector<std::vector<double>> &yIonMass,
								const std::vector<std::vector<double>> &bIonMass,
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
	for (const auto &env : yIonMass)
		if (envelopeMatches(env))
			++c.matchedY;
	for (const auto &env : bIonMass)
		if (envelopeMatches(env))
			++c.matchedB;
	return c;
}

EnvCounts countMatchedEnvelopes(const MS2Scan *scan,
								const SipRecord &rec,
								const MassTolerance &tol)
{
	return countMatchedEnvelopes(
		scan, rec.vvdYionMass, rec.vvdBionMass, tol);
}

// Match each library fragment once, then construct the two established feature
// alignments: entropy uses matched pairs only, while cosine retains unmatched
// positive library peaks with zero observed intensity.
void alignSpectraFeatures(
	const std::vector<double> &fragMz,
	const std::vector<double> &fragInt,
	const std::vector<double> &scanMz,
	const std::vector<double> &scanIntensity,
	const MassTolerance &tol,
	std::vector<double> &entropyP,
	std::vector<double> &entropyQ,
	std::vector<double> &cosineP,
	std::vector<double> &cosineQ)
{
	entropyP.clear();
	entropyQ.clear();
	cosineP.clear();
	cosineQ.clear();
	const size_t nFrag = std::min(fragMz.size(), fragInt.size());
	if (nFrag == 0)
		return;

	const size_t nScan = std::min(scanMz.size(), scanIntensity.size());
	const auto scanEnd = scanMz.begin() + static_cast<std::ptrdiff_t>(nScan);
	for (size_t i = 0; i < nFrag; ++i)
	{
		const double m = fragMz[i];
		const double intensity = fragInt[i];
		double bestI = 0.0;
		bool found = false;
		if (nScan > 0)
		{
			const double tolDa = tol.daAt(m);
			auto it = std::lower_bound(scanMz.begin(), scanEnd, m - tolDa);
			double bestErr = tolDa + 1.0;
			while (it != scanEnd && *it <= m + tolDa)
			{
				const double err = std::fabs(*it - m);
				if (err < bestErr)
				{
					bestErr = err;
					bestI = scanIntensity[static_cast<size_t>(it - scanMz.begin())];
					found = true;
				}
				++it;
			}
		}
		if (found)
		{
			entropyP.push_back(intensity);
			entropyQ.push_back(bestI);
		}
		if (intensity > 0.0 && std::isfinite(intensity) && std::isfinite(m))
		{
			cosineP.push_back(intensity);
			cosineQ.push_back(
				bestI > 0.0 && std::isfinite(bestI) ? bestI : 0.0);
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
	norm(entropyP);
	norm(entropyQ);
	norm(cosineP);
	norm(cosineQ);
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

bool buildPrecursorEnvelopeFromProductIons(
	const Isotopologue &isotopologue,
	const std::vector<std::vector<double>> &yIonMass,
	const std::vector<std::vector<double>> &yIonProb,
	const std::vector<std::vector<double>> &bIonMass,
	const std::vector<std::vector<double>> &bIonProb,
	IsotopeDistribution &precursor)
{
	if (bIonMass.empty() || bIonProb.empty() ||
		yIonMass.empty() || yIonProb.empty() ||
		bIonMass.back().empty() || bIonProb.back().empty() ||
		yIonMass.front().empty() || yIonProb.front().empty())
		return false;
	precursor = isotopologue.sum(
		IsotopeDistribution(bIonMass.back(), bIonProb.back()),
		IsotopeDistribution(yIonMass.front(), yIonProb.front()));
	return !precursor.vMass.empty() &&
		precursor.vMass.size() == precursor.vProb.size();
}

// Expand the compact SFI envelope intensities onto the full-resolution WDP
// isotope masses. SFI retains the empirical apex for every b/y envelope; WDP
// supplies the complete isotope shape used to distribute that apex intensity.
void buildHighResolutionFeatureSpectrum(
	const SipRecord &record,
	const std::vector<std::vector<double>> &yIonMass,
	const std::vector<std::vector<double>> &yIonProb,
	const std::vector<std::vector<double>> &bIonMass,
	const std::vector<std::vector<double>> &bIonProb,
	std::vector<std::pair<double, double>> &peaks,
	std::vector<double> &fragmentMz,
	std::vector<double> &fragmentIntensity)
{
	peaks.clear();
	const double proton = ProNovoConfig::getProtonMass();
	const auto appendSeries = [&peaks, proton](
		const std::vector<std::vector<double>> &masses,
		const std::vector<std::vector<double>> &probabilities,
		const std::vector<std::vector<double>> &compactExperimental)
	{
		const size_t envelopeCount = std::min(
			masses.size(), probabilities.size());
		for (size_t envelopeIndex = 0;
			 envelopeIndex < envelopeCount; ++envelopeIndex)
		{
			const size_t peakCount = std::min(
				masses[envelopeIndex].size(),
				probabilities[envelopeIndex].size());
			if (peakCount == 0 ||
				envelopeIndex >= compactExperimental.size() ||
				compactExperimental[envelopeIndex].empty())
				continue;
			const double empiricalApex = *std::max_element(
				compactExperimental[envelopeIndex].begin(),
				compactExperimental[envelopeIndex].end());
			const double theoreticalApex = *std::max_element(
				probabilities[envelopeIndex].begin(),
				probabilities[envelopeIndex].begin() +
					static_cast<std::ptrdiff_t>(peakCount));
			if (!(empiricalApex > 0.0) || !(theoreticalApex > 0.0))
				continue;
			for (size_t peakIndex = 0; peakIndex < peakCount; ++peakIndex)
			{
				if (!(masses[envelopeIndex][peakIndex] > 0.0) ||
					!(probabilities[envelopeIndex][peakIndex] > 0.0))
					continue;
				peaks.emplace_back(
					masses[envelopeIndex][peakIndex] + proton,
					empiricalApex *
						probabilities[envelopeIndex][peakIndex] /
						theoreticalApex);
			}
		}
	};
	appendSeries(yIonMass, yIonProb, record.vvdYionExpIntensity);
	appendSeries(bIonMass, bIonProb, record.vvdBionExpIntensity);
	std::sort(peaks.begin(), peaks.end(),
		[](const auto &left, const auto &right)
		{ return left.first < right.first; });
	fragmentMz.clear();
	fragmentIntensity.clear();
	fragmentMz.reserve(peaks.size());
	fragmentIntensity.reserve(peaks.size());
	for (const auto &peak : peaks)
	{
		fragmentMz.push_back(peak.first);
		fragmentIntensity.push_back(peak.second);
	}
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
	double precursorRtDiffSeconds = -1.0;
	// strings use the spectra-index decorated peptide/protein representation.
	std::string peptide;
	std::string proteins;
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
		row.precursorRtDiffSeconds = r.precursorRtDiffSeconds;
		row.peptide = r.peptide;
		row.proteins = r.proteins;
		pinRows.push_back(std::move(row));
	}
	return pinRows;
}

// -------------------- Parent: per-HDF5 sample driver --------------------

// Per-scan scoring state.
struct ScoringPsm
{
	size_t scanIdx = 0;
	size_t libraryIndex = 0;
	uint32_t sfiRecordId = 0;
	int scanNumber = 0;
	int label = 0;
	double wdp = -std::numeric_limits<double>::infinity();
	double xcorr = 0.0;
	double mvh = 0.0;
	int xcorrGeometryKey = 0;
	bool featuresReady = false;
	ShardPsmRow row;
	double ms2Pct = 0.0;
	std::string nakedPeptide;
	std::string proteins;
	std::string sipAtom;
};

int scoringGeometryKey(const MS2Scan *scan)
{
	const bool hasQuery = scan->pQuery != nullptr;
	const int maximumFragmentCharge = hasQuery
		? scan->pQuery->_spectrumInfoInternal.iMaxFragCharge : 0;
	return (hasQuery ? 1 << 16 : 0) |
		(scan->isMS2HighRes ? 1 << 8 : 0) |
		(maximumFragmentCharge & 0xff);
}

struct CandidateMatch
{
	uint32_t recordIdx = 0;
};

struct SfiGateCandidate
{
	uint32_t recordIdx = 0;
	uint32_t recordRtBin = 0;
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
		psm.proteins = mergeProteinLists(samePeptide->proteins, psm.proteins);
		if (psm.wdp > samePeptide->wdp)
		{
			*samePeptide = std::move(psm);
			std::sort(top.begin(), top.end(), byWdpDesc);
		}
		else
		{
			samePeptide->proteins = std::move(psm.proteins);
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

// The cascade needs only sequence metadata for WDP. Full SFI fragment and
// precursor envelopes are materialized after WDP ranking for the small set of
// retained PSMs. Keeping this object compact avoids millions of nested-vector
// allocations and deep copies in the hot ranking path.
struct CascadeRecord
{
	std::string peptide;
	std::string nakedPeptide;
	uint32_t sfiRecordId = 0;
	int charge = 1;
	double ms2Pct = 0.0;
};

SipRecord materializeSfiRecord(const sipros::SpectraIndex &index,
							   uint32_t recordId)
{
	const sipros::SpectraIndexRecord &source = index.record(recordId);
	SipRecord record;
	record.peptide = std::string(index.peptide(recordId));
	record.nakedPeptide = nakedPeptideOf(record.peptide);
	record.charge = source.charge;
	record.topPrecursorMz = source.topPrecursorMz;
	record.sumPrecursorIntensity = source.sumPrecursorIntensity;

	const size_t peptideLength = record.nakedPeptide.size();
	const size_t positions = peptideLength > 0 ? peptideLength - 1 : 0;
	record.vvdYionMass.assign(positions, {});
	record.vvdYionProb.assign(positions, {});
	record.vvdYionExpIntensity.assign(positions, {});
	record.vvdBionMass.assign(positions, {});
	record.vvdBionProb.assign(positions, {});
	record.vvdBionExpIntensity.assign(positions, {});
	std::vector<double> sumY(positions, 0.0);
	std::vector<double> sumB(positions, 0.0);
	const double proton = ProNovoConfig::getProtonMass();
	const auto fragmentRange = index.fragments(recordId);
	for (const auto *fragment = fragmentRange.first;
		 fragment != fragmentRange.second; ++fragment)
	{
		const double fragmentMz = fragment->mz();
		const size_t position = fragment->ionPosition;
		if (position == 0 || position > positions)
			continue;
		const double mass = fragmentMz - proton;
		const double probability = fragment->theoreticalIntensity;
		if (fragment->ionKind == static_cast<uint8_t>('y') ||
			fragment->ionKind == static_cast<uint8_t>('Y'))
		{
			record.vvdYionMass[position - 1].push_back(mass);
			record.vvdYionProb[position - 1].push_back(probability);
			record.vvdYionExpIntensity[position - 1].push_back(
				fragment->experimentalIntensity);
			sumY[position - 1] += probability;
		}
		else if (fragment->ionKind == static_cast<uint8_t>('b') ||
				 fragment->ionKind == static_cast<uint8_t>('B'))
		{
			record.vvdBionMass[position - 1].push_back(mass);
			record.vvdBionProb[position - 1].push_back(probability);
			record.vvdBionExpIntensity[position - 1].push_back(
				fragment->experimentalIntensity);
			sumB[position - 1] += probability;
		}
	}
	for (size_t position = 0; position < positions; ++position)
	{
		if (sumY[position] > 0.0)
			for (double &probability : record.vvdYionProb[position])
				probability /= sumY[position];
		if (sumB[position] > 0.0)
			for (double &probability : record.vvdBionProb[position])
				probability /= sumB[position];
		if (record.vvdYionMass[position].empty())
		{
			record.vvdYionMass[position].push_back(0.0);
			record.vvdYionProb[position].push_back(0.0);
			record.vvdYionExpIntensity[position].push_back(0.0);
		}
		if (record.vvdBionMass[position].empty())
		{
			record.vvdBionMass[position].push_back(0.0);
			record.vvdBionProb[position].push_back(0.0);
			record.vvdBionExpIntensity[position].push_back(0.0);
		}
	}
	return record;
}

PrecursorMatch matchSfiIsolationWindow(
	const MS2Scan *scan,
	const sipros::SpectraIndexRecord &record,
	double rtTolerance)
{
	PrecursorMatch best;
	double bestMzShift = std::numeric_limits<double>::infinity();
	const double scanRt = scanRtMinutes(scan);
	const double deltaRt = scanRt - record.retentionMinutes;
	if (std::fabs(deltaRt) > rtTolerance || record.charge < 1 ||
		record.charge > 4)
	{
		return best;
	}
	for (const auto &window : scan->vIsolationWindowsMz)
	{
		const double center = window.first;
		const double width = window.second;
		if (!(center > 0.0) || !(width > 0.0) ||
			!std::isfinite(center) || !std::isfinite(width) ||
			record.topPrecursorMz < center - width / 2.0 ||
			record.topPrecursorMz > center + width / 2.0)
		{
			continue;
		}
		const double mzShift = std::fabs(record.topPrecursorMz - center);
		if (mzShift >= bestMzShift)
			continue;
		bestMzShift = mzShift;
		best.ok = true;
		best.mzShift = mzShift;
		best.deltaRT = std::fabs(deltaRt);
		best.matchedCharge = record.charge;
		best.matchedMz = record.topPrecursorMz;
		best.isolationCenterMz = center;
	}
	return best;
}

struct SfiGateCounters
{
	size_t rtRejected = 0;
	size_t precursorCandidates = 0;
	size_t rtProductBins = 0;
	size_t productMassIntervals = 0;
	size_t productRangeLookups = 0;
	size_t productPostingMatches = 0;
};

struct SfiCandidateMetadata
{
	std::vector<double> retentionMinutes;
	std::vector<int32_t> charge;
	std::vector<uint32_t> rtBin;

	void build(const sipros::SpectraIndex &index)
	{
		const size_t count = static_cast<size_t>(index.recordCount());
		retentionMinutes.resize(count);
		charge.resize(count);
		rtBin.resize(count);
#pragma omp parallel for schedule(static)
		for (size_t recordId = 0; recordId < count; ++recordId)
		{
			const auto &record = index.record(static_cast<uint32_t>(recordId));
			retentionMinutes[recordId] = record.retentionMinutes;
			charge[recordId] = record.charge;
			rtBin[recordId] = static_cast<uint32_t>(std::floor(
				std::max(0.0, record.retentionMinutes) /
				sipros::SpectraIndex::rtBinWidthMinutes()));
		}
	}
};

struct ProductMassBinRange
{
	uint32_t lower = 0;
	uint32_t upper = 0;
};

ProductMassBinRange productMassBinRange(double lowerMz, double upperMz)
{
	constexpr double binWidth = 0.001;
	constexpr double maximumBin = static_cast<double>((1U << 24) - 1U);
	const double lower = std::floor(std::max(0.0, lowerMz) / binWidth);
	const double upper = std::ceil(std::max(0.0, upperMz) / binWidth);
	return {
		static_cast<uint32_t>(std::min(lower, maximumBin)),
		static_cast<uint32_t>(std::min(upper, maximumBin))};
}

void selectSfiDdaWindowGatePeaks(
	const MS2Scan &scan,
	std::vector<double> &gatePeakMzs,
	std::vector<std::pair<double, double>> &intensityMz)
{
	gatePeakMzs.clear();
	intensityMz.clear();
	if (scan.pPeakList == nullptr)
		return;
	const std::vector<double> &filteredMzs = scan.pPeakList->pPeaks;
	if (filteredMzs.size() <= kDdaWindowTopPeaks)
	{
		gatePeakMzs = filteredMzs;
		return;
	}

	// Regular search applies the top-200 gate after MVH preprocessing. Recover
	// those peaks' raw intensities, retain the most
	// intense 200, then restore m/z order for sparse-index range lookups.
	intensityMz.reserve(filteredMzs.size());
	size_t rawIndex = 0;
	for (double mz : filteredMzs)
	{
		while (rawIndex < scan.vdMZ.size() && scan.vdMZ[rawIndex] < mz)
			++rawIndex;
		const double intensity = rawIndex < scan.vdMZ.size() &&
			rawIndex < scan.vdIntensity.size() && scan.vdMZ[rawIndex] == mz
				? scan.vdIntensity[rawIndex] : 0.0;
		intensityMz.push_back({intensity, mz});
	}
	std::partial_sort(
		intensityMz.begin(),
		intensityMz.begin() +
			static_cast<std::ptrdiff_t>(kDdaWindowTopPeaks),
		intensityMz.end(),
		[](const auto &left, const auto &right)
		{
			if (left.first != right.first)
				return left.first > right.first;
			return left.second < right.second;
		});
	gatePeakMzs.reserve(kDdaWindowTopPeaks);
	for (size_t i = 0; i < kDdaWindowTopPeaks; ++i)
		gatePeakMzs.push_back(intensityMz[i].second);
	std::sort(gatePeakMzs.begin(), gatePeakMzs.end());
}

size_t assignSfiCandidatesToScans(
	const std::vector<MS2Scan *> &scans,
	const sipros::SpectraIndex &index,
	const SfiCandidateMetadata &candidateMetadata,
	const MassTolerance &fragmentTolerance,
	double rtTolerance,
	uint32_t minimumProductIons,
	std::vector<uint32_t> &candidateRecordIds,
	std::vector<std::vector<CandidateMatch>> &out,
	SfiGateCounters &counters)
{
	out.assign(scans.size(), {});
	size_t precursorCandidates = 0;
	size_t rtRejected = 0;
	size_t rtProductBins = 0;
	size_t productMassIntervals = 0;
	size_t productRangeLookups = 0;
	size_t productPostingMatches = 0;
	size_t survivors = 0;
	if (candidateMetadata.retentionMinutes.size() < index.recordCount() ||
		candidateMetadata.charge.size() < index.recordCount() ||
		candidateMetadata.rtBin.size() < index.recordCount())
	{
		throw std::runtime_error("incomplete SFI candidate metadata");
	}
	const double proton = ProNovoConfig::getProtonMass();
	const size_t recordWordCount =
		(static_cast<size_t>(index.recordCount()) + 63U) / 64U;
	const int gateThreadCount = std::max(1, omp_get_max_threads());
	// Each scan is owned by one OpenMP worker. Record its surviving SFI IDs in
	// that worker's bitset so the post-gate unique-record reduction does not
	// append and sort tens of millions of IDs on one core.
	std::vector<uint64_t> threadRecordBits(
		static_cast<size_t>(gateThreadCount) * recordWordCount, 0);
	std::vector<std::vector<SfiGateCandidate>> candidateScratch(
		static_cast<size_t>(gateThreadCount));
	std::vector<std::vector<double>> gatePeakScratch(
		static_cast<size_t>(gateThreadCount));
	std::vector<std::vector<std::pair<double, double>>> intensityMzScratch(
		static_cast<size_t>(gateThreadCount));
	std::vector<std::vector<ProductMassBinRange>> massRangeScratch(
		static_cast<size_t>(gateThreadCount));
#pragma omp parallel for schedule(guided, 4) reduction(+ : rtRejected, precursorCandidates, rtProductBins, productMassIntervals, productRangeLookups, productPostingMatches, survivors)
	for (size_t scanIndex = 0; scanIndex < scans.size(); ++scanIndex)
	{
		const int threadIndex = omp_get_thread_num();
		MS2Scan *scan = scans[scanIndex];
		const double scanRt = scanRtMinutes(scan);
		auto &candidates = candidateScratch[
			static_cast<size_t>(threadIndex)];
		candidates.clear();
		for (const auto &window : scan->vIsolationWindowsMz)
		{
			const double center = window.first;
			const double width = window.second;
			if (!(center > 0.0) || !(width > 0.0) ||
				!std::isfinite(center) || !std::isfinite(width))
			{
				continue;
			}
			const auto range = index.precursorMzRange(
				center - width / 2.0, center + width / 2.0);
			for (uint32_t recordId = range.first;
				 recordId < range.second; ++recordId)
			{
				const int recordCharge = candidateMetadata.charge[recordId];
				if (recordCharge < 1 || recordCharge > 4)
					continue;
				const double recordRt =
					candidateMetadata.retentionMinutes[recordId];
				if (std::fabs(scanRt - recordRt) > rtTolerance)
				{
					++rtRejected;
					continue;
				}
				candidates.push_back(
					{recordId, candidateMetadata.rtBin[recordId]});
			}
		}
		const auto candidateOrder =
			[](const SfiGateCandidate &left,
				const SfiGateCandidate &right)
			{
				return left.recordIdx < right.recordIdx;
			};
		// One acquisition window yields record IDs directly in SFI precursor-m/z
		// order, so its hundreds of millions of candidates need no comparison
		// sort. Multiple overlapping reactions require ordering before adjacent
		// deduplication.
		if (scan->vIsolationWindowsMz.size() > 1)
		{
			std::sort(candidates.begin(), candidates.end(), candidateOrder);
			candidates.erase(std::unique(candidates.begin(), candidates.end(),
				[](const SfiGateCandidate &left, const SfiGateCandidate &right)
				{ return left.recordIdx == right.recordIdx; }), candidates.end());
		}
		precursorCandidates += candidates.size();
		auto &gatePeakMzs = gatePeakScratch[static_cast<size_t>(threadIndex)];
		auto &intensityMz = intensityMzScratch[
			static_cast<size_t>(threadIndex)];
		selectSfiDdaWindowGatePeaks(*scan, gatePeakMzs, intensityMz);
		auto &massRanges = massRangeScratch[
			static_cast<size_t>(threadIndex)];
		massRanges.clear();
		if (!candidates.empty())
		{
			massRanges.reserve(gatePeakMzs.size());
			for (double observedMz : gatePeakMzs)
			{
				const double indexedMz = observedMz - proton + proton;
				const double tolerance = fragmentTolerance.daAt(observedMz);
				massRanges.push_back(productMassBinRange(
					indexedMz - tolerance, indexedMz + tolerance));
			}
			std::sort(massRanges.begin(), massRanges.end(),
				[](const ProductMassBinRange &left, const ProductMassBinRange &right)
				{
					if (left.lower != right.lower)
						return left.lower < right.lower;
					return left.upper < right.upper;
				});
			productMassIntervals += massRanges.size();
		}

		std::vector<CandidateMatch> gated;
		gated.reserve(std::min<size_t>(candidates.size(), 256U));
		for (size_t begin = 0; begin < candidates.size();)
		{
			constexpr uint32_t recordsPerBlock =
				sipros::SpectraIndex::recordBlockCapacity();
			const uint32_t block = static_cast<uint32_t>(
				candidates[begin].recordIdx) / recordsPerBlock;
			const uint32_t recordBegin = block * recordsPerBlock;
			size_t end = begin + 1;
			while (end < candidates.size() &&
				static_cast<uint32_t>(candidates[end].recordIdx) /
					recordsPerBlock == block)
				++end;
			std::array<uint32_t, sipros::SpectraIndex::recordBlockCapacity()> hits;
			hits.fill(std::numeric_limits<uint32_t>::max());
			for (size_t i = begin; i < end; ++i)
			{
				const size_t localRecord = candidates[i].recordIdx - recordBegin;
				hits[localRecord] = 0;
			}
			if (minimumProductIons > 0)
			{
				const auto rtBins = index.rtBins(
					block, scanRt - rtTolerance, scanRt + rtTolerance);
				if (rtBins.first != nullptr)
				{
					for (const auto *rtBin = rtBins.first;
						 rtBin != rtBins.second; ++rtBin)
					{
						// A precursor block is mass ordered and can contain records from
						// many runs/RT regions.  Only visit an RT segment when it owns an
						// exact precursor candidate, and stop visiting that segment once
						// all of its candidates reach the product-ion threshold.  This is
						// the same bounded-block/early-exit strategy used by the regular
						// FASTA peptide-cache query.
						size_t candidatesBelowThreshold = 0;
						for (size_t i = begin; i < end; ++i)
						{
							const uint32_t candidateRtBin = candidates[i].recordRtBin;
							if (candidateRtBin == rtBin->rtBin &&
								hits[candidates[i].recordIdx - recordBegin] <
									minimumProductIons)
							{
								++candidatesBelowThreshold;
							}
						}
						if (candidatesBelowThreshold == 0)
							continue;
						++rtProductBins;
						if (massRanges.empty())
							continue;
						const auto segment = index.productPostings(block, *rtBin);
						if (segment.first == nullptr)
							continue;
						constexpr size_t intervalsPerLookup = 1;
						const sipros::SpectraIndexFragmentPosting *postingSearchBegin =
							segment.first;
						for (size_t groupBegin = 0;
							 groupBegin < massRanges.size() && candidatesBelowThreshold != 0;
							 groupBegin += intervalsPerLookup)
						{
							const size_t groupEnd = std::min(
								massRanges.size(), groupBegin + intervalsPerLookup);
							const uint32_t lowerPacked =
								massRanges[groupBegin].lower << 8U;
							const uint32_t upperPacked =
								(massRanges[groupEnd - 1].upper << 8U) | 0xffU;
							const sipros::SpectraIndexFragmentPosting *first =
								postingSearchBegin;
							while (first != segment.second &&
								first->packed < lowerPacked)
							{
								++first;
							}
							auto *last = first;
							while (last != segment.second &&
								last->packed <= upperPacked)
							{
								++last;
							}
							postingSearchBegin = groupEnd < massRanges.size() &&
								massRanges[groupEnd].lower <= massRanges[groupEnd - 1].upper
								? first : last;
							++productRangeLookups;
							productPostingMatches += static_cast<size_t>(last - first);
							for (const auto *posting = first; posting != last; ++posting)
							{
								const uint8_t localRecordId = posting->localRecordId();
								uint32_t &hit = hits[localRecordId];
								if (hit < minimumProductIons &&
									++hit == minimumProductIons)
								{
									--candidatesBelowThreshold;
								}
							}
						}
					}
				}
			}
			for (size_t i = begin; i < end; ++i)
			{
				if (minimumProductIons == 0 ||
					hits[candidates[i].recordIdx - recordBegin] >= minimumProductIons)
				{
					gated.push_back({candidates[i].recordIdx});
					++survivors;
				}
			}
			begin = end;
		}
		uint64_t *recordBits = threadRecordBits.data() +
			static_cast<size_t>(threadIndex) * recordWordCount;
		for (const CandidateMatch &candidate : gated)
		{
			const size_t recordId = candidate.recordIdx;
			recordBits[recordId >> 6U] |= uint64_t{1} << (recordId & 63U);
		}
		out[scanIndex] = std::move(gated);
	}
	std::vector<uint64_t> combinedRecordBits(recordWordCount, 0);
	std::vector<size_t> recordWordOffsets(recordWordCount + 1U, 0);
#pragma omp parallel for schedule(static)
	for (size_t word = 0; word < recordWordCount; ++word)
	{
		uint64_t combined = 0;
		for (int thread = 0; thread < gateThreadCount; ++thread)
			combined |= threadRecordBits[
				static_cast<size_t>(thread) * recordWordCount + word];
		combinedRecordBits[word] = combined;
		recordWordOffsets[word + 1U] = static_cast<size_t>(
			populationCount64(combined));
	}
	for (size_t word = 0; word < recordWordCount; ++word)
		recordWordOffsets[word + 1U] += recordWordOffsets[word];
	candidateRecordIds.resize(recordWordOffsets.back());
#pragma omp parallel for schedule(static)
	for (size_t word = 0; word < recordWordCount; ++word)
	{
		uint64_t bits = combinedRecordBits[word];
		size_t output = recordWordOffsets[word];
		while (bits != 0)
		{
			const unsigned bit = trailingZeroCount64(bits);
			candidateRecordIds[output++] = static_cast<uint32_t>(
				(word << 6U) + bit);
			bits &= bits - 1U;
		}
	}
	counters.rtRejected += rtRejected;
	counters.precursorCandidates += precursorCandidates;
	counters.rtProductBins += rtProductBins;
	counters.productMassIntervals += productMassIntervals;
	counters.productRangeLookups += productRangeLookups;
	counters.productPostingMatches += productPostingMatches;
	return survivors;
}

void materializeSfiCascadeMetadata(
	const sipros::SpectraIndex &index,
	const std::vector<uint32_t> &recordIds,
	std::vector<std::vector<CandidateMatch>> &scanCandidates,
	std::vector<CascadeRecord> &materialized)
{
	materialized.clear();
	materialized.resize(recordIds.size());
	std::vector<uint32_t> materializedPositions(
		static_cast<size_t>(index.recordCount()),
		std::numeric_limits<uint32_t>::max());
#pragma omp parallel for schedule(guided, 32)
	for (size_t i = 0; i < recordIds.size(); ++i)
	{
		const uint32_t recordId = recordIds[i];
		CascadeRecord record;
		record.peptide = std::string(index.peptide(recordId));
		record.nakedPeptide = nakedPeptideOf(record.peptide);
		record.sfiRecordId = recordId;
		record.charge = index.record(recordId).charge;
		record.ms2Pct = index.record(recordId).sipAbundancePct;
		materialized[i] = std::move(record);
		materializedPositions[recordId] = static_cast<uint32_t>(i);
	}
#pragma omp parallel for schedule(guided, 32)
	for (size_t scanIndex = 0; scanIndex < scanCandidates.size(); ++scanIndex)
	{
		auto &candidates = scanCandidates[scanIndex];
		for (CandidateMatch &candidate : candidates)
		{
			const uint32_t recordId = static_cast<uint32_t>(candidate.recordIdx);
			candidate.recordIdx = materializedPositions[recordId];
		}
	}
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

void buildSfiMvhIonCache(
	const sipros::SpectraIndex &index,
	const std::vector<uint32_t> &recordIds,
	int envelopeTopN,
	std::vector<std::vector<double>> &recordIons)
{
	recordIons.clear();
	recordIons.resize(static_cast<size_t>(index.recordCount()));
#pragma omp parallel for schedule(guided, 32)
	for (size_t i = 0; i < recordIds.size(); ++i)
	{
		const uint32_t recordId = recordIds[i];
		const std::string peptide(index.peptide(recordId));
		std::string nakedPeptide = nakedPeptideOf(peptide);
		const size_t positions = nakedPeptide.empty()
			? 0 : nakedPeptide.size() - 1U;
		std::vector<std::vector<double>> yMass(positions);
		std::vector<std::vector<double>> yProb(positions);
		std::vector<std::vector<double>> bMass(positions);
		std::vector<std::vector<double>> bProb(positions);
		std::vector<double> sumY(positions, 0.0);
		std::vector<double> sumB(positions, 0.0);
		const double proton = ProNovoConfig::getProtonMass();
		const auto fragmentRange = index.fragments(recordId);
		for (const auto *fragment = fragmentRange.first;
			 fragment != fragmentRange.second; ++fragment)
		{
			const size_t position = fragment->ionPosition;
			if (position == 0 || position > positions)
				continue;
			const double mass = fragment->mz() - proton;
			const double probability = fragment->theoreticalIntensity;
			if (fragment->ionKind == static_cast<uint8_t>('y') ||
				fragment->ionKind == static_cast<uint8_t>('Y'))
			{
				yMass[position - 1U].push_back(mass);
				yProb[position - 1U].push_back(probability);
				sumY[position - 1U] += probability;
			}
			else if (fragment->ionKind == static_cast<uint8_t>('b') ||
				 fragment->ionKind == static_cast<uint8_t>('B'))
			{
				bMass[position - 1U].push_back(mass);
				bProb[position - 1U].push_back(probability);
				sumB[position - 1U] += probability;
			}
		}
		for (size_t position = 0; position < positions; ++position)
		{
			if (sumY[position] > 0.0)
				for (double &probability : yProb[position])
					probability /= sumY[position];
			if (sumB[position] > 0.0)
				for (double &probability : bProb[position])
					probability /= sumB[position];
			if (yMass[position].empty())
			{
				yMass[position].push_back(0.0);
				yProb[position].push_back(0.0);
			}
			if (bMass[position].empty())
			{
				bMass[position].push_back(0.0);
				bProb[position].push_back(0.0);
			}
		}
		keepTopNEnvelopePeaks(
			yMass, yProb, envelopeTopN);
		keepTopNEnvelopePeaks(
			bMass, bProb, envelopeTopN);
		std::vector<char> residues;
		MVH::CalculateSequenceIonsSIP(
			nakedPeptide, index.record(recordId).charge,
			MVH::bUseSmartPlusThreeModel, &recordIons[recordId],
			yMass, yProb, bMass, bProb, &residues);
	}
}

struct RankedMvhCandidate
{
	CandidateMatch candidate;
	double score = 0.0;
};

bool betterMvhCandidate(
	const RankedMvhCandidate &left,
	const RankedMvhCandidate &right)
{
	if (left.score != right.score)
		return left.score > right.score;
	return left.candidate.recordIdx < right.candidate.recordIdx;
}

size_t scoreAndRetainTopMvhCandidates(
	std::vector<std::vector<CandidateMatch>> &scanCandidates,
	const std::vector<MS2Scan *> &scans,
	const std::vector<std::vector<double>> &recordMvhIons,
	std::vector<std::vector<double>> &candidateMvh,
	size_t perScanLimit,
	size_t &acceptedCount)
{
	size_t retained = 0;
	size_t acceptedTotal = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : retained, acceptedTotal)
	for (size_t scanIndex = 0; scanIndex < scanCandidates.size(); ++scanIndex)
	{
		auto &candidates = scanCandidates[scanIndex];
		auto &scores = candidateMvh[scanIndex];
		std::vector<RankedMvhCandidate> heap;
		heap.reserve(std::min(perScanLimit, candidates.size()));
		if (!scans[scanIndex]->bSkip)
		{
			for (const CandidateMatch &candidate : candidates)
			{
				double score = 0.0;
				if (!MVH::ScoreIonsVsSpectrum(
						recordMvhIons[candidate.recordIdx], scans[scanIndex], score))
				{
					continue;
				}
				++acceptedTotal;
				RankedMvhCandidate ranked{candidate, score};
				if (heap.size() < perScanLimit)
				{
					heap.push_back(std::move(ranked));
					std::push_heap(
						heap.begin(), heap.end(), betterMvhCandidate);
				}
				else if (perScanLimit != 0 &&
					betterMvhCandidate(ranked, heap.front()))
				{
					std::pop_heap(
						heap.begin(), heap.end(), betterMvhCandidate);
					heap.back() = std::move(ranked);
					std::push_heap(
						heap.begin(), heap.end(), betterMvhCandidate);
				}
			}
		}
		std::sort(heap.begin(), heap.end(), betterMvhCandidate);
		const size_t keep = heap.size();
		std::vector<CandidateMatch> keptCandidates;
		std::vector<double> keptScores;
		keptCandidates.reserve(keep);
		keptScores.reserve(keep);
		for (RankedMvhCandidate &ranked : heap)
		{
			keptCandidates.push_back(std::move(ranked.candidate));
			keptScores.push_back(ranked.score);
		}
		candidates = std::move(keptCandidates);
		scores = std::move(keptScores);
		retained += keep;
	}
	acceptedCount = acceptedTotal;
	return retained;
}

std::vector<uint32_t> uniqueSfiCandidateRecordIds(
	const std::vector<std::vector<CandidateMatch>> &scanCandidates)
{
	std::vector<uint32_t> recordIds;
	for (const auto &candidates : scanCandidates)
		for (const CandidateMatch &candidate : candidates)
			recordIds.push_back(static_cast<uint32_t>(candidate.recordIdx));
	std::sort(recordIds.begin(), recordIds.end());
	recordIds.erase(std::unique(recordIds.begin(), recordIds.end()),
		recordIds.end());
	return recordIds;
}

ShardPsmRow makeScoringRow(size_t scanIdx,
                           const MS2Scan *scan,
                           const SipRecord &rec,
                           const std::string &sipAtom,
                           double ms2Pct,
                           const PrecursorMatch &match,
                           const EnvCounts &envCounts,
                           const sipros::RaxportMs1Data *ms1Data,
                           const MassTolerance &mzTol,
						   const std::vector<double> &precursorNeutralMasses,
						   const std::vector<double> &precursorProbabilities,
                           double wdp,
                           double xcorr,
                           double mvh,
                           double entropy,
                           double cosine)
{
    double baseMass = 0.0;
    const int precursorCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
    std::vector<isotopicPeak> ms1Peaks;
    const bool missingMatchedPrecursor =
        match.precursorRtDiffSeconds < 0.0;
    const int ms1AnchorScanNumber = missingMatchedPrecursor && scan
        ? scan->iParentScanID : match.matchedMs1ScanNumber;
    const double ms1AnchorMz = missingMatchedPrecursor
        ? match.isolationCenterMz : match.matchedMz;
    if (ms1Data && scan && ms1AnchorScanNumber > 0 &&
        ms1AnchorMz > 0.0 && precursorCharge > 0)
    {
        PeptideIsotopeCalculator calculator;
        const std::string peptideForComposition = PSMfeatureExtractor::peptideBodyWithPtms(rec.peptide);
        baseMass = calculator.calPrecursorBaseMass(peptideForComposition);
        int ms1ScanNumber = ms1AnchorScanNumber;
        const auto mzToleranceDaAt = [&](double mz)
        { return mzTol.daAt(mz); };
        ms1Peaks = PSMfeatureExtractor::findMs1IsotopicPeaksFromEnvelope(
            ms1Data,
            ms1ScanNumber,
            precursorCharge,
            baseMass,
            ms1AnchorMz,
            precursorNeutralMasses,
            precursorProbabilities,
            mzToleranceDaAt);
    }
    const PSMfeatureExtractor::Ms1AbundanceResult ms1Abundance =
        PSMfeatureExtractor::getSIPelementAbundanceFromMS1PeaksWithEnvelope(
            ms1Peaks, baseMass, rec.peptide, precursorCharge,
            sipAtom, ms2Pct,
            precursorNeutralMasses, precursorProbabilities);
    const double ms1AbundancePct = missingMatchedPrecursor &&
            ms1Abundance.rawIsotopicPeakCount == 0
        ? PSMfeatureExtractor::estimateSIPelementAbundanceFromIsolationCenter(
              rec.peptide, precursorCharge, match.isolationCenterMz,
              sipAtom, ms2Pct)
        : ms1Abundance.abundancePct;

    ShardPsmRow row;
    row.scanIdx = static_cast<int32_t>(scanIdx);
    row.parentCharge = rec.charge > 0 ? rec.charge : match.matchedCharge;
    row.isotopicShift = match.isotopicShift;
    row.matchedY = envCounts.matchedY;
    row.matchedB = envCounts.matchedB;
    row.peptideLength = static_cast<int>(rec.nakedPeptide.size());
    row.missCleavage = PSMfeatureExtractor::countMissCleavage(rec.nakedPeptide);
    row.ptmCount = PSMfeatureExtractor::countPTM(rec.peptide);
    row.isotopicPeakNumbers = missingMatchedPrecursor
        ? 0 : ms1Abundance.rawIsotopicPeakCount;
    row.ms1IsotopeFitScore = missingMatchedPrecursor
        ? 0.0 : ms1Abundance.fitScore;
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
    row.ms1IsotopicAbundancePct = ms1AbundancePct;
    row.log10PrecursorIntensity = (rec.sumPrecursorIntensity > 0.0)
                                      ? std::log10(rec.sumPrecursorIntensity)
                                      : 0.0;
	row.precursorRtDiffSeconds = match.precursorRtDiffSeconds;
    row.peptide = rec.peptide;
    row.proteins = rec.proteins;
    return row;
}

bool loadAndValidateSfiLibraries(
	const std::vector<std::string> &paths,
	std::vector<sipros::SpectraIndex> &indices,
	std::string &canonicalSipIsotope,
	int expectedLabel,
	std::string &error)
{
	indices.clear();
	canonicalSipIsotope.clear();
	if (paths.empty())
	{
		error = "No .sfi SIP spectra libraries were provided.";
		return false;
	}
	indices.reserve(paths.size());
	std::string chemistryProfile;
	uint32_t envelopeTopN = 0;
	for (const std::string &path : paths)
	{
		sipros::SpectraIndex index;
		if (!index.load(path, error))
			return false;
		const auto &metadata = index.metadata();
		if (chemistryProfile.empty())
		{
			chemistryProfile = metadata.chemistryProfileId;
			std::string chemistryError;
			if (!ProNovoConfig::configureChemistryProfileId(
					chemistryProfile, chemistryError))
			{
				error = "Cannot use spectra index " + path + ": " + chemistryError;
				return false;
			}
		}
		else if (metadata.chemistryProfileId != chemistryProfile)
		{
			error = "Mixed chemistry_profile_id values in SFI spectra libraries.";
			return false;
		}
		if (metadata.envelopeTopN == 0)
		{
			error = "SFI library has no compact-envelope top-N metadata: " + path;
			return false;
		}
		if (envelopeTopN == 0)
			envelopeTopN = metadata.envelopeTopN;
		else if (metadata.envelopeTopN != envelopeTopN)
		{
			error = "Mixed envelope_top_n values in SFI spectra libraries.";
			return false;
		}
		if (metadata.label != 1 && metadata.label != -1)
		{
			error = "Invalid target/decoy label in " + path;
			return false;
		}
		if (metadata.label != expectedLabel)
		{
			error = "SFI label metadata does not match --sfi-label in " + path;
			return false;
		}
		const std::string isotope = PSMfeatureExtractor::canonicalSipIsotope(
			std::string(1, metadata.sipAtom), metadata.sipIsotopeMassNumber);
		if (isotope.empty())
		{
			error = "Unsupported SIP isotope metadata in " + path;
			return false;
		}
		if (canonicalSipIsotope.empty())
			canonicalSipIsotope = isotope;
		else if (canonicalSipIsotope != isotope)
		{
			error = "Mixed SIP isotope targets in SFI spectra libraries.";
			return false;
		}
		indices.push_back(std::move(index));
	}
	return true;
}

std::string fixedPtmSummary()
{
	const auto names = ProNovoConfig::getEnabledFixedPtmNames();
	if (names.empty())
		return "none";
	std::ostringstream output;
	for (size_t index = 0; index < names.size(); ++index)
	{
		if (index != 0)
			output << ", ";
		output << names[index];
	}
	output << " (implicit in every applicable SFI peptide)";
	return output.str();
}

int processOneHdf5(const Args &args, const std::string &scanPath,
                   const std::vector<std::string> &sfiFiles)
{
	const sipros::PerformanceTimer processTimer;
	TimingLogger timing;

	const std::string sampleBasename =
		std::filesystem::path(scanPath).stem().string();

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
	std::vector<sipros::SpectraIndex> spectraIndices;
	if (!loadAndValidateSfiLibraries(
			sfiFiles, spectraIndices, canonicalSipIsotope, args.sfiLabel,
			libraryMetadataError))
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
	const uint32_t storedEnvelopeTopN =
		spectraIndices.front().metadata().envelopeTopN;
	int sipIsotopeIndex = -1;
	try
	{
		sipIsotopeIndex = ProNovoConfig::resolveSipIsotopeIndex(
			ProNovoConfig::configIsotopologue,
			canonicalSipIsotope[0],
			std::stoi(canonicalSipIsotope.substr(1)));
	}
	catch (const std::exception &exception)
	{
		std::cerr << "Cannot initialize high-resolution WDP envelopes: "
				  << exception.what() << "\n";
		return 2;
	}
	const Isotopologue pristineWdpIsotopologue =
		ProNovoConfig::configIsotopologue;
	ProNovoConfig::setMassAccuracy(parentTolDa, fragmentToleranceDa);

	const int nThreads = std::max(1, args.threads);
	omp_set_num_threads(nThreads);
	uint64_t totalRtBins = 0;
	for (const auto &index : spectraIndices)
		totalRtBins += index.rtBinCount();
	const char *labelName = args.sfiLabel == -1 ? "decoy" : "target";
	std::ostringstream gateDescription;
	gateDescription
		<< "window-only z1-z4, top " << kDdaWindowTopPeaks
		<< " MS2 peaks, singly charged b/y ions, "
		<< std::max<uint32_t>(
			kDdaWindowMinMatchedFragments,
			static_cast<uint32_t>(
				std::max(0, ProNovoConfig::MinMatchedFragments)))
		<< " matches";
	std::cout << "\nSFI spectra search\n"
			  << "  Scan file : " << scanPath << "\n"
			  << "  SFI       : "
			  << (std::filesystem::path(args.sfiDir) / "*.sfi").string()
			  << " (" << labelName << ")\n"
			  << "  Profile   : " << ProNovoConfig::getChemistryProfileId() << "\n"
			  << "  Fixed PTMs: " << fixedPtmSummary() << "\n"
			  << "  Var PTMs  : encoded in SFI peptide forms\n"
			  << "  SIP       : " << canonicalSipIsotope
			  << " abundance stored per SFI record\n"
			  << "  Tolerance : MS1=" << args.toleranceMs1
			  << (args.toleranceMs1Ppm ? " ppm" : " Da")
			  << ", MS2=" << args.toleranceMs2
			  << (args.toleranceMs2Ppm ? " ppm" : " Da")
			  << " (scoring: parent=" << parentTolDa
			  << " Da, fragment=" << fragmentToleranceDa << " Da)\n";
	std::cout << "  RT window : +/-" << args.rtToleranceMin << " minutes; "
			  << sipros::SpectraIndex::rtBinWidthMinutes()
			  << " minute bins; " << sipros::formatPerformanceCount(totalRtBins)
			  << " sparse segments\n"
			  << "  Gate      : " << gateDescription.str() << "\n"
			  << "  Cascade   : top " << args.mvhCascadeTopN
			  << " MVH candidates per scan to Xcorr/WDP; "
			  << storedEnvelopeTopN << " compact peaks/envelope; top "
			  << args.topPsmsPerScan << " PSMs per label\n"
			  << "  WDP       : regenerated full high-resolution SIP envelopes\n"
			  << "  MS1       : post-score parent +/-"
			  << kPrecursorMatchScanRadius
			  << " scans; absolute RT-distance feature; no MS1 gate\n"
			  << "  Search    : RT-aware isolation-window/product SFI index\n";

	// -------- Read HDF5 scans and score assigned SIP spectra --------
	timing.printHeader(
		"Search timing: " +
			std::filesystem::path(scanPath).filename().string(), nThreads);

	std::vector<MS2Scan *> scans;
	sipros::RaxportMs1Data ms1Data;
	std::string readError;
	sipros::RaxportReadOptions readOptions;
	readOptions.precursorSource = sipros::PrecursorSource::IsolationWindow;
	readOptions.precursorMatchScanRadius = kPrecursorMatchScanRadius;
	bool readHdf5Ok = false;
	timing.run("Load HDF5 scans", "read HDF5 MS1/MS2 scans", 0, "", [&]()
			   { readHdf5Ok = sipros::readRaxportHdf5Scans(
					 scanPath, scans, &ms1Data, readError, nullptr, readOptions); });
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
		const double proton = ProNovoConfig::getProtonMass();
		for (const auto &window : scan->vIsolationWindowsMz)
		{
			const double upperMz = window.first + window.second / 2.0;
			if (!(upperMz > 0.0) || !std::isfinite(upperMz))
				continue;
			for (int charge = 1; charge <= 4; ++charge)
			{
				const double chargedMass = upperMz * charge;
				const double neutralMass = chargedMass -
					static_cast<double>(charge) * proton;
				scan->dParentMass = std::max(
					scan->dParentMass, chargedMass);
				scan->dParentNeutralMass = std::max(
					scan->dParentNeutralMass, neutralMass);
				scan->iMaxCandidateCharge = std::max(
					scan->iMaxCandidateCharge, charge);
				ProNovoConfig::dMaxMS2ScanMass = std::max(
					ProNovoConfig::dMaxMS2ScanMass, chargedMass);
			}
		}
		ProNovoConfig::iMaxPercusorCharge = std::max(
			ProNovoConfig::iMaxPercusorCharge, 4);
	}

	const sipros::RaxportMs1Data *ms1DataPtr = &ms1Data;
	auto preprocessScans = [&]()
	{
#pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < scans.size(); ++i)
			scans[i]->preprocess();
	};
	timing.run("Preprocess spectra", "preprocess HDF5 MS2 scans",
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

	// Window-only XCorr contains no information beyond the highest observed
	// fragment and its background-correlation tail. Size both preprocessing
	// and sparse theoretical-ion scratch to that common bound instead of the
	// much larger hypothetical z=4 neutral precursor mass.
	double maximumObservedMz = 0.0;
	for (const MS2Scan *scan : scans)
	{
		if (!scan->vdMZ.empty())
			maximumObservedMz = std::max(
				maximumObservedMz,
				*std::max_element(scan->vdMZ.begin(), scan->vdMZ.end()));
	}
	const double correlationTailMz =
		static_cast<double>(ProNovoConfig::iXcorrProcessingOffset + 12) /
		ProNovoConfig::dHighResInverseBinWidth;
	const double maximumXcorrMz = maximumObservedMz + correlationTailMz;
	int iArraySizePreprocess = static_cast<int>((maximumXcorrMz + 3 + 2.0) * ProNovoConfig::dHighResInverseBinWidth);
	int iArraySizeScoreSip = iArraySizePreprocess;
	int iArraySizeScore = static_cast<int>((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	// Per-thread scratch
	std::vector<std::unique_ptr<double[]>> tRaw(nThreads), tFXc(nThreads), tCorr(nThreads), tSmooth(nThreads), tPeak(nThreads);
	std::vector<std::vector<unsigned char>> tDuplSip(nThreads);
	std::vector<std::vector<double>> tBinIonSip(nThreads);
	std::vector<std::vector<int>> tBinSip(nThreads);
	std::vector<std::unique_ptr<multimap<double, double>>> tMvhSorted(nThreads);
	// First-touch the large per-worker buffers on the NUMA node of the worker
	// that will use them. Serial initialization put every page on node 0 and
	// made the second socket pay remote-memory latency.
#pragma omp parallel for schedule(static)
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
			const double maximumObservedMz = !scan->vdMZ.empty()
				? *std::max_element(scan->vdMZ.begin(), scan->vdMZ.end())
				: -1.0;

			Query *pQuery = new Query();
			if (CometSearchMod::Preprocess(pQuery, scan,
										   tRaw[tid].get(), tFXc[tid].get(),
										   tCorr[tid].get(), tSmooth[tid].get(),
										   tPeak[tid].get(), maximumObservedMz))
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

	// Each optimized process searches exactly one SFI label.
	std::vector<std::vector<ScoringPsm>> topPsms(scans.size());

	size_t totalLibraryRecords = 0;
	for (const auto &index : spectraIndices)
		totalLibraryRecords += static_cast<size_t>(index.recordCount());
	size_t totalAssignedCandidates = 0;
	size_t totalMvhAccepted = 0;
	size_t totalCascadeCandidates = 0;
	size_t totalWdpEnvelopeFailures = 0;
	SfiGateCounters gateCounters;
	std::unordered_set<double> uniqueWdpAbundances;
	uniqueWdpAbundances.reserve(256);
	auto collectWdpAbundances = [&]()
	{
		std::vector<std::unordered_set<double>> perThread(
			static_cast<size_t>(nThreads));
		for (auto &values : perThread)
			values.reserve(256);
		for (const sipros::SpectraIndex &index : spectraIndices)
		{
#pragma omp parallel for schedule(static)
			for (uint32_t recordId = 0;
				 recordId < index.recordCount(); ++recordId)
			{
				perThread[static_cast<size_t>(omp_get_thread_num())].emplace(
					index.record(recordId).sipAbundancePct);
			}
		}
		for (const auto &values : perThread)
			uniqueWdpAbundances.insert(values.begin(), values.end());
	};
	timing.run("Collect SFI abundances",
		"collect distinct WDP isotope states",
		totalLibraryRecords, "records", collectWdpAbundances);
	std::vector<double> wdpAbundances(
		uniqueWdpAbundances.begin(), uniqueWdpAbundances.end());
	std::sort(wdpAbundances.begin(), wdpAbundances.end());
	std::vector<Isotopologue> wdpIsotopologues(
		wdpAbundances.size(), pristineWdpIsotopologue);
	std::vector<std::string> wdpIsotopologueErrors(wdpAbundances.size());
	auto prepareWdpIsotopologues = [&]()
	{
#pragma omp parallel for schedule(dynamic, 1)
		for (size_t abundanceIndex = 0;
			 abundanceIndex < wdpAbundances.size(); ++abundanceIndex)
		{
			try
			{
				ProNovoConfig::setSipAbundance(
					wdpIsotopologues[abundanceIndex],
					canonicalSipIsotope[0], sipIsotopeIndex,
					wdpAbundances[abundanceIndex]);
			}
			catch (const std::exception &exception)
			{
				wdpIsotopologueErrors[abundanceIndex] = exception.what();
			}
		}
	};
	timing.run("Prepare WDP isotope states",
		"prepare full-resolution WDP isotope states",
		wdpAbundances.size(), "SIP abundances", prepareWdpIsotopologues);
	for (size_t abundanceIndex = 0;
		 abundanceIndex < wdpIsotopologueErrors.size(); ++abundanceIndex)
	{
		if (!wdpIsotopologueErrors[abundanceIndex].empty())
		{
			std::cerr << "Cannot prepare WDP isotope state for "
					  << wdpAbundances[abundanceIndex] << "%: "
					  << wdpIsotopologueErrors[abundanceIndex] << "\n";
			return 2;
		}
	}

	for (size_t libraryIndex = 0; libraryIndex < spectraIndices.size(); ++libraryIndex)
	{
		const sipros::SpectraIndex &spectraIndex = spectraIndices[libraryIndex];
		SfiCandidateMetadata candidateMetadata;
		timing.run("Prepare candidate metadata",
			"compact SFI RT/charge metadata",
			static_cast<size_t>(spectraIndex.recordCount()),
			"records", [&]() { candidateMetadata.build(spectraIndex); });
		std::vector<CascadeRecord> batchRecords;
		std::vector<uint32_t> candidateRecordIds;
		const size_t nRec = static_cast<size_t>(spectraIndex.recordCount());
		if (nRec == 0)
			continue;

		std::vector<std::vector<CandidateMatch>> scanCandidates(scans.size());
		size_t assignedCandidates = 0;
		auto assignBatch = [&]()
		{
			const uint32_t minimumProductIons = std::max<uint32_t>(
				kDdaWindowMinMatchedFragments,
				static_cast<uint32_t>(
					std::max(0, ProNovoConfig::MinMatchedFragments)));
			assignedCandidates = assignSfiCandidatesToScans(
				scans, spectraIndex, candidateMetadata,
				fragTol, rtTol,
				minimumProductIons,
				candidateRecordIds, scanCandidates, gateCounters);
		};
		std::ostringstream assignLabel;
		assignLabel << "precursor/product gate SFI " << (libraryIndex + 1)
					<< '/' << spectraIndices.size();
		timing.run("Query spectra index", assignLabel.str(),
				   nRec, "records", assignBatch);
		totalAssignedCandidates += assignedCandidates;

		std::vector<std::vector<double>> candidateMvh(scans.size());
		std::vector<std::vector<double>> recordMvhIons;
		auto prepareMvhCache = [&]()
		{
			buildSfiMvhIonCache(spectraIndex, candidateRecordIds,
				static_cast<int>(storedEnvelopeTopN), recordMvhIons);
		};
		timing.run("Prepare MVH ions", "prepare cached SFI MVH ions",
			candidateRecordIds.size(), "records", prepareMvhCache);

		const size_t cascadeLimit = std::max<size_t>(
			static_cast<size_t>(args.mvhCascadeTopN),
			static_cast<size_t>(args.topPsmsPerScan));
		size_t batchMvhAccepted = 0;
		size_t cascadeCandidates = 0;
		auto scoreMvhBatch = [&]()
		{
			cascadeCandidates = scoreAndRetainTopMvhCandidates(
				scanCandidates, scans, recordMvhIons,
				candidateMvh,
				cascadeLimit, batchMvhAccepted);
		};
		std::ostringstream mvhLabel;
		mvhLabel << "score MVH SFI " << (libraryIndex + 1)
				 << '/' << spectraIndices.size();
		timing.run("MVH scoring and cascade", mvhLabel.str(),
				   assignedCandidates, "candidates", scoreMvhBatch);
		totalMvhAccepted += batchMvhAccepted;
		totalCascadeCandidates += cascadeCandidates;
		std::vector<std::vector<double>>().swap(recordMvhIons);
		candidateRecordIds = uniqueSfiCandidateRecordIds(scanCandidates);
		auto materializeCascade = [&]()
		{
			materializeSfiCascadeMetadata(spectraIndex, candidateRecordIds,
				scanCandidates, batchRecords);
		};
		timing.run("Prepare cascade metadata", "materialize MVH top-N sequence metadata",
			candidateRecordIds.size(), "records", materializeCascade);

		std::vector<std::vector<std::vector<double>>> wdpYMass(nThreads);
		std::vector<std::vector<std::vector<double>>> wdpYProb(nThreads);
		std::vector<std::vector<std::vector<double>>> wdpBMass(nThreads);
		std::vector<std::vector<std::vector<double>>> wdpBProb(nThreads);
		std::vector<std::vector<double>> candidateWdp(scans.size());
		std::vector<size_t> recordOccurrenceOffsets;
		std::vector<uint64_t> recordOccurrences;
		auto prepareScoreReuse = [&]()
		{
			std::vector<uint32_t> counts(batchRecords.size(), 0);
			for (size_t s = 0; s < scanCandidates.size(); ++s)
			{
				candidateWdp[s].assign(scanCandidates[s].size(),
					std::numeric_limits<double>::quiet_NaN());
				for (const CandidateMatch &candidate : scanCandidates[s])
					++counts[candidate.recordIdx];
			}
			recordOccurrenceOffsets.assign(batchRecords.size() + 1, 0);
			for (size_t recordIndex = 0; recordIndex < counts.size(); ++recordIndex)
				recordOccurrenceOffsets[recordIndex + 1] =
					recordOccurrenceOffsets[recordIndex] + counts[recordIndex];
			recordOccurrences.resize(recordOccurrenceOffsets.back());
			std::vector<size_t> cursor = recordOccurrenceOffsets;
			for (size_t s = 0; s < scanCandidates.size(); ++s)
			{
				for (size_t i = 0; i < scanCandidates[s].size(); ++i)
				{
					const size_t recordIndex = scanCandidates[s][i].recordIdx;
					recordOccurrences[cursor[recordIndex]++] =
						(static_cast<uint64_t>(s) << 32U) |
						static_cast<uint32_t>(i);
				}
			}
		};
		timing.run("Prepare score reuse", "group candidates by SFI record",
			cascadeCandidates, "candidates", prepareScoreReuse);

		size_t batchWdpEnvelopeFailures = 0;

		auto scoreWdpBatch = [&]()
		{
			size_t wdpEnvelopeFailures = 0;
#pragma omp parallel for schedule(guided, 32) reduction(+ : wdpEnvelopeFailures)
			for (size_t recordIndex = 0;
				 recordIndex < batchRecords.size(); ++recordIndex)
			{
				const int tid = omp_get_thread_num();
				const size_t occurrenceBegin =
					recordOccurrenceOffsets[recordIndex];
				const size_t occurrenceEnd =
					recordOccurrenceOffsets[recordIndex + 1];
				if (occurrenceBegin == occurrenceEnd)
					continue;
				const CascadeRecord &rec = batchRecords[recordIndex];
				const auto abundance = std::lower_bound(
					wdpAbundances.begin(), wdpAbundances.end(), rec.ms2Pct);
				if (abundance == wdpAbundances.end() || *abundance != rec.ms2Pct)
				{
					wdpEnvelopeFailures += occurrenceEnd - occurrenceBegin;
					continue;
				}
				const size_t abundanceIndex = static_cast<size_t>(
					abundance - wdpAbundances.begin());
				if (!wdpIsotopologues[abundanceIndex].computeProductIon(
						rec.peptide,
						wdpYMass[tid], wdpYProb[tid],
						wdpBMass[tid], wdpBProb[tid]))
				{
					wdpEnvelopeFailures += occurrenceEnd - occurrenceBegin;
					continue;
				}
				for (size_t occurrenceIndex = occurrenceBegin;
					 occurrenceIndex < occurrenceEnd; ++occurrenceIndex)
				{
					const uint64_t occurrence =
						recordOccurrences[occurrenceIndex];
					const size_t s = static_cast<size_t>(occurrence >> 32U);
					const size_t i = static_cast<uint32_t>(occurrence);
					candidateWdp[s][i] = scans[s]->scoreWeightSumHighMS2(
						&rec.nakedPeptide, rec.charge,
						&wdpYMass[tid], &wdpYProb[tid],
						&wdpBMass[tid], &wdpBProb[tid]);
				}
			}
			batchWdpEnvelopeFailures += wdpEnvelopeFailures;
		};
		std::ostringstream wdpLabel;
		wdpLabel << "score WDP with record reuse SFI " << (libraryIndex + 1)
				 << '/' << spectraIndices.size();
		timing.run("WDP scoring", wdpLabel.str(),
				   cascadeCandidates, "candidates", scoreWdpBatch);
		auto rankWdpBatch = [&]()
		{
#pragma omp parallel for schedule(dynamic, 1)
			for (size_t s = 0; s < scans.size(); ++s)
			{
				MS2Scan *scan = scans[s];
				for (size_t i = 0; i < scanCandidates[s].size(); ++i)
				{
					const double wdp = candidateWdp[s][i];
					if (!std::isfinite(wdp))
						continue;
					const CandidateMatch &candidate = scanCandidates[s][i];
					const CascadeRecord &rec = batchRecords[candidate.recordIdx];
					std::vector<ScoringPsm> &topList = topPsms[s];
					if (!wouldKeepTopUniquePeptide(
							topList, rec.nakedPeptide, wdp,
							args.topPsmsPerScan))
						continue;
					ScoringPsm psm;
					psm.scanIdx = s;
					psm.libraryIndex = libraryIndex;
					psm.sfiRecordId = rec.sfiRecordId;
					psm.scanNumber = scan->iScanId;
					psm.label = spectraIndex.metadata().label;
					psm.wdp = wdp;
					psm.mvh = candidateMvh[s][i];
					psm.xcorrGeometryKey = scoringGeometryKey(scan);
					psm.ms2Pct = rec.ms2Pct;
					psm.nakedPeptide = rec.nakedPeptide;
					psm.proteins = normalizeProteinList(std::string(
						spectraIndex.proteins(rec.sfiRecordId)));
					psm.sipAtom = canonicalSipIsotope;
					pushTopUniquePeptide(
						topList, std::move(psm), args.topPsmsPerScan);
				}
			}
		};
		timing.run("Rank WDP candidates", "retain unique-peptide WDP winners",
			cascadeCandidates, "candidates", rankWdpBatch);
		totalWdpEnvelopeFailures += batchWdpEnvelopeFailures;
	}
	std::vector<std::vector<std::vector<double>>> featureYMass(nThreads);
	std::vector<std::vector<std::vector<double>>> featureYProb(nThreads);
	std::vector<std::vector<std::vector<double>>> featureBMass(nThreads);
	std::vector<std::vector<std::vector<double>>> featureBProb(nThreads);
	std::vector<std::vector<double>> wdpFeatureMz(nThreads);
	std::vector<std::vector<double>> wdpFeatureIntensity(nThreads);
	std::vector<std::vector<std::pair<double, double>>> featurePeaks(nThreads);
	std::vector<std::vector<double>> entropyAlignedP(nThreads);
	std::vector<std::vector<double>> entropyAlignedQ(nThreads);
	std::vector<std::vector<double>> cosineAlignedP(nThreads);
	std::vector<std::vector<double>> cosineAlignedQ(nThreads);
	size_t rankedPsms = 0;
	for (const auto &v : topPsms)
		rankedPsms += v.size();
	size_t postScorePrecursorMatches = 0;
	auto finalizeRetainedPsms = [&]()
	{
		std::vector<ScoringPsm *> retained;
		retained.reserve(rankedPsms);
		for (auto &scanPsms : topPsms)
			for (ScoringPsm &psm : scanPsms)
				retained.push_back(&psm);
		if (retained.empty())
			return size_t{0};
		std::sort(retained.begin(), retained.end(),
			[](const ScoringPsm *left, const ScoringPsm *right)
			{
				if (left->libraryIndex != right->libraryIndex)
					return left->libraryIndex < right->libraryIndex;
				if (left->sfiRecordId != right->sfiRecordId)
					return left->sfiRecordId < right->sfiRecordId;
				return left->xcorrGeometryKey < right->xcorrGeometryKey;
			});
		std::vector<size_t> groupOffsets;
		groupOffsets.reserve(retained.size() + 1U);
		groupOffsets.push_back(0);
		for (size_t i = 1; i < retained.size(); ++i)
		{
			if (retained[i]->libraryIndex != retained[i - 1U]->libraryIndex ||
				retained[i]->sfiRecordId != retained[i - 1U]->sfiRecordId ||
				retained[i]->xcorrGeometryKey !=
					retained[i - 1U]->xcorrGeometryKey)
			{
				groupOffsets.push_back(i);
			}
		}
		groupOffsets.push_back(retained.size());
		const size_t groupCount = groupOffsets.size() - 1U;

		long long failures = 0;
		long long matchedPrecursors = 0;
#pragma omp parallel for schedule(guided, 16) reduction(+ : failures, matchedPrecursors)
		for (size_t groupIndex = 0; groupIndex < groupCount; ++groupIndex)
		{
			const int tid = omp_get_thread_num();
			const size_t groupBegin = groupOffsets[groupIndex];
			const size_t groupEnd = groupOffsets[groupIndex + 1U];
			ScoringPsm &firstPsm = *retained[groupBegin];
			SipRecord rec = materializeSfiRecord(
				spectraIndices[firstPsm.libraryIndex], firstPsm.sfiRecordId);
			const auto abundance = std::lower_bound(
				wdpAbundances.begin(), wdpAbundances.end(), firstPsm.ms2Pct);
			if (abundance == wdpAbundances.end() ||
				*abundance != firstPsm.ms2Pct)
			{
				failures += static_cast<long long>(groupEnd - groupBegin);
				continue;
			}
			const size_t abundanceIndex = static_cast<size_t>(
				abundance - wdpAbundances.begin());
			if (!wdpIsotopologues[abundanceIndex].computeProductIon(
					rec.peptide,
					featureYMass[tid], featureYProb[tid],
					featureBMass[tid], featureBProb[tid]))
			{
				failures += static_cast<long long>(groupEnd - groupBegin);
				continue;
			}
			buildHighResolutionFeatureSpectrum(
				rec,
				featureYMass[tid], featureYProb[tid],
				featureBMass[tid], featureBProb[tid],
				featurePeaks[tid],
				wdpFeatureMz[tid], wdpFeatureIntensity[tid]);
			IsotopeDistribution precursorEnvelope;
			if (!buildPrecursorEnvelopeFromProductIons(
					wdpIsotopologues[abundanceIndex],
					featureYMass[tid], featureYProb[tid],
					featureBMass[tid], featureBProb[tid],
					precursorEnvelope))
			{
				failures += static_cast<long long>(groupEnd - groupBegin);
				continue;
			}

			MS2Scan *xcorrGeometryScan = scans[firstPsm.scanIdx];
			if (xcorrGeometryScan->pQuery != nullptr)
			{
				CometSearchMod::BinPeptideIonsSIPNoCancelOut(
					rec.vvdYionMass, rec.vvdYionProb,
					rec.vvdBionMass, rec.vvdBionProb, xcorrGeometryScan,
					tDuplSip[tid], tBinIonSip[tid], tBinSip[tid]);
			}

			for (size_t occurrence = groupBegin;
				 occurrence < groupEnd; ++occurrence)
			{
				ScoringPsm &psm = *retained[occurrence];
				MS2Scan *scan = scans[psm.scanIdx];
				if (scan->pQuery != nullptr)
				{
					CometSearchMod::ScoreBinnedPeptideSIPNoCancelOut(
						scan, tBinIonSip[tid], tBinSip[tid], psm.xcorr);
				}
				alignSpectraFeatures(
					wdpFeatureMz[tid], wdpFeatureIntensity[tid],
					scan->vdMZ, scan->vdIntensity, fragTol,
					entropyAlignedP[tid], entropyAlignedQ[tid],
					cosineAlignedP[tid], cosineAlignedQ[tid]);
				const double entropy = computeEntropy(
					entropyAlignedP[tid], entropyAlignedQ[tid]);
				const double cosine = computeCosine(
					cosineAlignedP[tid], cosineAlignedQ[tid]);
				const EnvCounts envelopeCounts = countMatchedEnvelopes(
					scan, featureYMass[tid], featureBMass[tid], fragTol);
				const sipros::SpectraIndex &psmIndex =
					spectraIndices[psm.libraryIndex];
				const sipros::SpectraIndexRecord &psmRecord =
					psmIndex.record(psm.sfiRecordId);
				PrecursorMatch featureMatch =
					matchSfiIsolationWindow(scan, psmRecord, rtTol);
				featureMatch = matchNearbyEnvelopePrecursor(
					ms1Data, scan->iParentScanID, scanRtMinutes(scan),
					rec.charge,
					precursorEnvelope.vMass, precursorEnvelope.vProb,
					mzTol, kPrecursorMatchScanRadius, featureMatch);
				if (featureMatch.matchedMs1ScanNumber > 0)
					++matchedPrecursors;
				rec.proteins = psm.proteins;
				psm.row = makeScoringRow(
					psm.scanIdx, scan, rec, psm.sipAtom, psm.ms2Pct,
					featureMatch,
					envelopeCounts, ms1DataPtr, mzTol,
					precursorEnvelope.vMass, precursorEnvelope.vProb,
					psm.wdp, psm.xcorr, psm.mvh, entropy, cosine);
				psm.featuresReady = true;
			}
		}
		auto removeFailures = [](std::vector<std::vector<ScoringPsm>> &perScan)
		{
#pragma omp parallel for schedule(static)
			for (size_t s = 0; s < perScan.size(); ++s)
			{
				auto &psms = perScan[s];
				psms.erase(std::remove_if(psms.begin(), psms.end(),
					[](const ScoringPsm &psm) { return !psm.featuresReady; }),
					psms.end());
			}
		};
		removeFailures(topPsms);
		postScorePrecursorMatches += static_cast<size_t>(matchedPrecursors);
		return static_cast<size_t>(failures);
	};
	timing.run("Finalize retained PSM features",
		"materialize and score only WDP winners",
		rankedPsms, "PSMs", [&]()
		{
			totalWdpEnvelopeFailures += finalizeRetainedPsms();
		});
	size_t retainedPsms = 0;
	for (const auto &v : topPsms)
		retainedPsms += v.size();

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

	timing.run("Build PIN columns", "materialize PIN rows",
			   retainedPsms, "PSMs", [&]()
			   {
				   std::vector<std::vector<MergedRow>> keptByScan(scans.size());
#pragma omp parallel for schedule(guided)
				   for (size_t s = 0; s < scans.size(); ++s)
				   {
					   std::vector<const ScoringPsm *> scanPsms;
					   scanPsms.reserve(topPsms[s].size());
					   appendPassing(scanPsms, topPsms[s]);
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
						   keptByScan[s].push_back(std::move(m));
					   }
				   }
				   size_t keptCount = 0;
				   for (const auto &scanRows : keptByScan)
					   keptCount += scanRows.size();
				   kept.reserve(keptCount);
				   for (auto &scanRows : keptByScan)
					   for (MergedRow &row : scanRows)
						   kept.push_back(std::move(row));
			   });

	const std::string pinSuffix = args.sfiLabel == 1
		? "_target.pin" : "_decoy.pin";
	std::filesystem::path pinPath = std::filesystem::path(args.outputDir) /
		(sampleBasename + pinSuffix);
	size_t writtenRows = 0;
	timing.run("Format and write PIN", "write PIN file",
			   0, pinPath.string(), [&]()
			   {
				   const std::vector<PinWriter::SearchSpectraPinRow> pinRows = makeSearchSpectraPinRows(kept);
				   writtenRows = PinWriter::writeSearchSpectraPin(pinPath.string(), sampleBasename, pinRows);
			   });
	// Cleanup
	timing.run("Cleanup scans", "cleanup scans",
			   scans.size(), "scans", [&]()
			   {
				   const long long scanCount =
					   static_cast<long long>(scans.size());
#pragma omp parallel for schedule(static)
				   for (long long i = 0; i < scanCount; ++i)
				   {
					   delete scans[static_cast<size_t>(i)];
				   }
			   });

	const sipros::PerformanceTiming fileTiming = processTimer.elapsed();
	double accountedWallSeconds = 0.0;
	double accountedCpuSeconds = 0.0;
	for (const TimingEntry &entry : timing.entries)
	{
		accountedWallSeconds += entry.timing.wallSeconds;
		accountedCpuSeconds += entry.timing.cpuSeconds;
	}
	// The file timer also covers initialization and coordination outside the
	// individually timed stages (SFI/profile setup, buffer allocation,
	// cascade pruning/ID compaction, and transient container cleanup).  Print
	// that remainder explicitly so the non-indented stage rows reconcile with
	// Search file (total).  The indented query rows are diagnostic sub-stages
	// and must not be added a second time.
	sipros::PerformanceTiming overheadTiming;
	overheadTiming.wallSeconds = std::max(
		0.0, fileTiming.wallSeconds - accountedWallSeconds);
	overheadTiming.cpuSeconds = std::max(
		0.0, fileTiming.cpuSeconds - accountedCpuSeconds);
	sipros::printPerformanceStage(
		std::cout, "Workflow setup/overhead", overheadTiming);
	sipros::printPerformanceStage(
		std::cout, "Search file (total)", fileTiming);
	timing.printFooter();
	std::cout << "\nSearch statistics\n"
			  << "  SFI library records      : "
			  << sipros::formatPerformanceCount(totalLibraryRecords) << "\n"
			  << "  Window+RT candidates     : "
			  << sipros::formatPerformanceCount(gateCounters.precursorCandidates)
			  << " (" << sipros::formatPerformanceCount(gateCounters.rtRejected)
			  << " RT rejected)\n"
			  << "  Product posting matches  : "
			  << sipros::formatPerformanceCount(gateCounters.productPostingMatches)
			  << "\n"
			  << "  Product RT segments      : "
			  << sipros::formatPerformanceCount(gateCounters.rtProductBins)
			  << "\n"
			  << "  Product mass intervals   : "
			  << sipros::formatPerformanceCount(gateCounters.productMassIntervals)
			  << "\n"
			  << "  Product range lookups    : "
			  << sipros::formatPerformanceCount(gateCounters.productRangeLookups)
			  << "\n"
			  << "  Product-gate survivors   : "
			  << sipros::formatPerformanceCount(totalAssignedCandidates) << "\n"
			  << "  MVH accepted             : "
			  << sipros::formatPerformanceCount(totalMvhAccepted) << "\n"
			  << "  WDP cascade candidates   : "
			  << sipros::formatPerformanceCount(totalCascadeCandidates) << "\n"
			  << "  WDP envelope failures    : "
			  << sipros::formatPerformanceCount(totalWdpEnvelopeFailures) << "\n"
			  << "  Finalized WDP winners    : "
			  << sipros::formatPerformanceCount(retainedPsms) << "\n"
			  << "  Post-score MS1 matches   : "
			  << sipros::formatPerformanceCount(postScorePrecursorMatches)
			  << "/" << sipros::formatPerformanceCount(retainedPsms) << "\n"
			  << "  PIN rows                 : "
			  << sipros::formatPerformanceCount(writtenRows) << "\n"
			  << "  PIN output               : " << pinPath.string() << "\n";
	return 0;
}

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

	std::vector<std::string> sfiFiles = listFiles(args.sfiDir, {".sfi", ".SFI"});
	sfiFiles.erase(std::remove_if(sfiFiles.begin(), sfiFiles.end(),
		[&](const std::string &path)
		{
			return (args.sfiLabel == -1) != isDecoySfi(path);
		}), sfiFiles.end());
	if (sfiFiles.empty())
	{
		std::cerr << "No matching SFI spectra library files in " << args.sfiDir
				  << "; HDF5 spectra libraries are not supported\n";
		return 1;
	}

	if (!std::filesystem::is_directory(args.outputDir))
		std::filesystem::create_directories(args.outputDir);

	omp_set_num_threads(std::max(1, args.threads));
	std::cout << "Sipros spectra search\n"
			  << "  Input files: " << scanFiles.size() << "\n"
			  << "  Mode       : search H5 scans against SFI\n"
			  << "  Label      : "
				  << (args.sfiLabel == -1 ? "decoy" : "target")
			  << "\n"
			  << "  Threads    : " << omp_get_max_threads() << "\n";
	const sipros::PerformanceTimer runTimer;
	int rc = 0;
	for (const auto &scanPath : scanFiles)
	{
		int r = processOneHdf5(args, scanPath, sfiFiles);
		if (r != 0)
			rc = r;
	}
	sipros::printPerformanceHeader(
		std::cout, "Run summary", omp_get_max_threads());
	sipros::printPerformanceStage(
		std::cout, "Complete spectra search", runTimer.elapsed(),
		sipros::formatPerformanceCount(scanFiles.size()) + " input file(s)");
	sipros::printPerformanceFooter(std::cout);
	std::cout << "\nSipros spectra search complete.\n";
	return rc;
}

} // namespace sipros
