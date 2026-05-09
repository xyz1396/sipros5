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
#include <map>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <H5Cpp.h>
#include "averagine.h"
#include "proNovoConfig.h"

#if !defined(_WIN32)
#include <cerrno>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace
{
struct Args
{
	std::string configPath;
	std::string inputPath;
	std::string ft2Path;
	std::string outputPath;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double fixedSipAbundancePct = 1.0;
	bool sipAbundanceRange = false;
	double sipAbundanceStartPct = 1.0;
	double sipAbundanceEndPct = 1.0;
	double sipAbundanceStepPct = 1.0;
	double probCutoff = 0.01;
	double ppmTolerance = 10.0;
	size_t minMatchedEnvelopes = 3;
	int threads = 0;
	bool writeDecoy = false;
	unsigned int decoySeed = 1;
};

struct PsmRow
{
	std::string psmId;
	std::string sample;
	int scanNumber = 0;
	std::string peptide;
	int precursorCharge = 1;
	std::string retentionText = "0";
	double probability = 0.0;
	double expectation = std::numeric_limits<double>::infinity();
	double hyperscore = -std::numeric_limits<double>::infinity();
	size_t order = 0;
};

struct ReadStats
{
	size_t totalRows = 0;
	size_t shortRows = 0;
	size_t invalidRows = 0;
	size_t unsupportedMods = 0;
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
	double targetSipAbundancePct = 0.0;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double probCutoff = 0.0;
	double ppmTolerance = 0.0;
	uint64_t minMatchedEnvelopes = 0;
};

struct Ft2SampleTask
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
	size_t written = 0;
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
	std::atomic<size_t> missingFt2{0};
	std::atomic<size_t> missingScan{0};
	std::atomic<size_t> computeFailed{0};
	std::atomic<size_t> precursorFailed{0};
	std::atomic<size_t> unmatched{0};
	std::atomic<size_t> targetFailed{0};
	std::atomic<size_t> decoyCollision{0};
	std::atomic<size_t> decoyComputeFailed{0};
};

int atomIndex(char sipAtom);

std::string trim(const std::string &s)
{
	size_t start = 0;
	while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])) != 0)
	{
		++start;
	}
	size_t end = s.size();
	while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0)
	{
		--end;
	}
	return s.substr(start, end - start);
}

std::string toLower(std::string s)
{
	for (char &ch : s)
	{
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	return s;
}

std::vector<std::string> splitTab(const std::string &line)
{
	std::vector<std::string> fields;
	std::stringstream ss(line);
	std::string field;
	while (std::getline(ss, field, '\t'))
	{
		fields.push_back(field);
	}
	if (!line.empty() && line.back() == '\t')
	{
		fields.emplace_back();
	}
	return fields;
}

void printUsage(const char *prog)
{
	std::cerr << "Usage: " << prog
			  << " -c <config.cfg> -i <psm.tsv|frag_dir> -f <ft2_file|sipros_dir> -o <output.h5|output_dir/>"
			  << " [-a <SIP atom/isotope, e.g. C13,H2,O18,N15,S34>]"
			  << " [-b <fixed SIP pct|lower-upper, default 1.0>] [-s|--step <pct, default 1.0>]"
			  << " [-p <prob cutoff, default 0.01>]"
			  << " [--ppm <match tolerance, default 10>] [--min-matched-envelopes <N, default 3>]"
			  << " [--decoy] [--decoy-seed <N, default 1>] [-t <threads>]\n";
	std::cerr << "FT2 matching is always performed at the baseline 1% C13 abundance; -b controls shifted output abundance(s).\n";
	std::cerr << "When one file is produced, -o is used as the output file and .h5 is appended if needed.\n";
	std::cerr << "When multiple files are produced, -o is used as an output directory unless it already names one.\n";
}

bool parseSipAtomSpec(const std::string &spec, char &sipAtom, int &sipIsotopeMassNumber)
{
	const std::string t = trim(spec);
	if (t.empty())
	{
		return false;
	}
	const char atom = static_cast<char>(std::toupper(static_cast<unsigned char>(t[0])));
	if (atomIndex(atom) < 0)
	{
		return false;
	}
	sipAtom = atom;

	if (t.size() == 1)
	{
		sipIsotopeMassNumber = -1;
		if (atom == 'O')
		{
			sipIsotopeMassNumber = 18;
		}
		else if (atom == 'S')
		{
			sipIsotopeMassNumber = 34;
		}
		return true;
	}

	for (size_t i = 1; i < t.size(); ++i)
	{
		if (!std::isdigit(static_cast<unsigned char>(t[i])))
		{
			return false;
		}
	}

	try
	{
		sipIsotopeMassNumber = std::stoi(t.substr(1));
	}
	catch (const std::exception &)
	{
		return false;
	}

	if (atom == 'C')
	{
		return sipIsotopeMassNumber == 13;
	}
	if (atom == 'H')
	{
		return sipIsotopeMassNumber == 2;
	}
	if (atom == 'N')
	{
		return sipIsotopeMassNumber == 15;
	}
	if (atom == 'O')
	{
		return sipIsotopeMassNumber == 18;
	}
	if (atom == 'S')
	{
		return sipIsotopeMassNumber == 34;
	}
	return false;
}

bool parseDoubleStrict(const std::string &text, double &value)
{
	try
	{
		const std::string t = trim(text);
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
	const std::string t = trim(spec);
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
		if (opt == "-c")
		{
			if (!requireValue(opt))
			{
				return false;
			}
			args.configPath = argv[++i];
		}
		else if (opt == "-i")
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
			args.ft2Path = argv[++i];
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
			if (!parseSipAtomSpec(sipSpec, args.sipAtom, args.sipIsotopeMassNumber))
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

	if (args.configPath.empty() || args.inputPath.empty() || args.ft2Path.empty() || args.outputPath.empty())
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
		const std::string t = trim(s);
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
		const std::string t = trim(s);
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

bool isApprox(double value, double target, double tolerance)
{
	return std::abs(value - target) <= tolerance;
}

bool parseModificationMass(const std::string &text, double &value)
{
	std::string t = trim(text);
	if (!t.empty() && t.front() == '+')
	{
		t.erase(t.begin());
	}
	return parseDoubleField(t, value);
}

bool applyResidueModification(char residue, const std::string &modText, std::string &body, std::string &reason)
{
	double mass = 0.0;
	if (!parseModificationMass(modText, mass))
	{
		reason = "invalid modification mass " + modText;
		return false;
	}

	if (residue == 'M' && (isApprox(mass, 147.0, 0.5) || isApprox(mass, 15.9949, 0.05)))
	{
		body.push_back('~');
		return true;
	}
	if (residue == 'C' && (isApprox(mass, 160.0, 0.5) || isApprox(mass, 57.0215, 0.05)))
	{
		return true;
	}
	if ((residue == 'N' || residue == 'Q') &&
		(isApprox(mass, residue == 'N' ? 115.0 : 129.0, 0.5) || isApprox(mass, 0.984, 0.05)))
	{
		body.push_back('!');
		return true;
	}

	reason = std::string("unsupported modification ") + residue + "[" + modText + "]";
	return false;
}

bool applyNTermModification(const std::string &modText, std::string &body, std::string &reason)
{
	double mass = 0.0;
	if (!parseModificationMass(modText, mass))
	{
		reason = "invalid N-terminal modification mass " + modText;
		return false;
	}
	if (isApprox(mass, 43.0, 0.5) || isApprox(mass, 42.0106, 0.05))
	{
		body.push_back('%');
		return true;
	}

	reason = "unsupported N-terminal modification n[" + modText + "]";
	return false;
}

bool convertModifiedPeptide(const std::string &plainPeptide,
							const std::string &modifiedPeptide,
							std::string &normalizedPeptide,
							std::string &reason)
{
	const std::string source = trim(modifiedPeptide).empty() ? trim(plainPeptide) : trim(modifiedPeptide);
	if (source.empty())
	{
		reason = "empty peptide";
		return false;
	}

	std::string body;
	for (size_t i = 0; i < source.size();)
	{
		const char raw = source[i];
		if (i == 0 && raw == 'n' && i + 1 < source.size() && source[i + 1] == '[')
		{
			const size_t rb = source.find(']', i + 2);
			if (rb == std::string::npos)
			{
				reason = "unterminated N-terminal modification";
				return false;
			}
			const std::string modText = source.substr(i + 2, rb - i - 2);
			if (!applyNTermModification(modText, body, reason))
			{
				return false;
			}
			i = rb + 1;
			continue;
		}
		if (std::isspace(static_cast<unsigned char>(raw)) != 0)
		{
			++i;
			continue;
		}
		if (!std::isalpha(static_cast<unsigned char>(raw)))
		{
			reason = "unsupported peptide character";
			return false;
		}

		const char residue = static_cast<char>(std::toupper(static_cast<unsigned char>(raw)));
		body.push_back(residue);
		++i;
		if (i < source.size() && source[i] == '[')
		{
			const size_t rb = source.find(']', i + 1);
			if (rb == std::string::npos)
			{
				reason = "unterminated modification";
				return false;
			}
			const std::string modText = source.substr(i + 1, rb - i - 1);
			if (!applyResidueModification(residue, modText, body, reason))
			{
				return false;
			}
			i = rb + 1;
		}
	}

	if (body.empty())
	{
		reason = "empty peptide";
		return false;
	}
	normalizedPeptide = "[" + body + "]";
	return true;
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
	const std::vector<std::string> parts = splitDot(trim(spectrum));
	if (parts.size() < 3)
	{
		return false;
	}
	sample = parts[0];
	if (sample.empty() || !parseIntField(parts[1], scanNumber))
	{
		return false;
	}
	if (parts.size() >= 4)
	{
		int parsedCharge = 0;
		if (parseIntField(parts.back(), parsedCharge) && parsedCharge > 0)
		{
			charge = parsedCharge;
		}
	}
	return scanNumber > 0;
}

std::vector<fs::path> collectPsmFiles(const std::string &inputPath)
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
			if (!entry.is_regular_file())
			{
				continue;
			}
			if (entry.path().filename() == "psm.tsv")
			{
				files.push_back(entry.path());
			}
		}
	}
	std::sort(files.begin(), files.end());
	if (files.empty())
	{
		throw std::runtime_error("No PSM TSV files found under: " + inputPath);
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

	const std::vector<std::string> headers = splitTab(headerLine);
	std::unordered_map<std::string, size_t> columns;
	for (size_t i = 0; i < headers.size(); ++i)
	{
		columns[trim(headers[i])] = i;
	}

	const size_t idxSpectrum = getRequiredColumn(columns, "Spectrum");
	const size_t idxPeptide = getRequiredColumn(columns, "Peptide");
	const size_t idxCharge = getRequiredColumn(columns, "Charge");
	const size_t idxRetention = getRequiredColumn(columns, "Retention");
	const size_t idxProbability = getRequiredColumn(columns, "Probability");
	const size_t idxModifiedPeptide = getOptionalColumn(columns, "Modified Peptide");
	const size_t idxExpectation = getOptionalColumn(columns, "Expectation");
	const size_t idxHyperscore = getOptionalColumn(columns, "Hyperscore");

	size_t requiredMax = std::max({idxSpectrum, idxPeptide, idxCharge, idxRetention, idxProbability});
	if (idxModifiedPeptide != std::string::npos)
	{
		requiredMax = std::max(requiredMax, idxModifiedPeptide);
	}

	std::string line;
	while (std::getline(in, line))
	{
		if (line.empty())
		{
			continue;
		}
		++stats.totalRows;
		const std::vector<std::string> fields = splitTab(line);
		if (fields.size() <= requiredMax)
		{
			++stats.shortRows;
			continue;
		}

		PsmRow row;
		row.psmId = trim(fields[idxSpectrum]);
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
		row.retentionText = trim(fields[idxRetention]);
		if (row.retentionText.empty())
		{
			row.retentionText = "0";
		}
		if (!parseDoubleField(fields[idxProbability], row.probability))
		{
			++stats.invalidRows;
			continue;
		}
		if (idxExpectation != std::string::npos && fields.size() > idxExpectation)
		{
			double value = 0.0;
			if (parseDoubleField(fields[idxExpectation], value))
			{
				row.expectation = value;
			}
		}
		if (idxHyperscore != std::string::npos && fields.size() > idxHyperscore)
		{
			double value = 0.0;
			if (parseDoubleField(fields[idxHyperscore], value))
			{
				row.hyperscore = value;
			}
		}

		const std::string modified = (idxModifiedPeptide != std::string::npos && fields.size() > idxModifiedPeptide)
										 ? fields[idxModifiedPeptide]
										 : std::string();
		std::string reason;
		if (!convertModifiedPeptide(fields[idxPeptide], modified, row.peptide, reason))
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

bool isBetterPsm(const PsmRow &candidate, const PsmRow &current)
{
	constexpr double eps = 1e-15;
	if (candidate.probability > current.probability + eps)
	{
		return true;
	}
	if (std::abs(candidate.probability - current.probability) <= eps)
	{
		if (candidate.expectation < current.expectation - eps)
		{
			return true;
		}
		if (std::abs(candidate.expectation - current.expectation) <= eps)
		{
			if (candidate.hyperscore > current.hyperscore + eps)
			{
				return true;
			}
			if (std::abs(candidate.hyperscore - current.hyperscore) <= eps)
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
		const std::string key = row.peptide + "\t" + std::to_string(row.precursorCharge);
		const auto it = bestIndex.find(key);
		if (it == bestIndex.end())
		{
			bestIndex[key] = selected.size();
			selected.push_back(row);
			continue;
		}
		PsmRow &current = selected[it->second];
		if (isBetterPsm(row, current))
		{
			current = row;
		}
	}

	std::sort(selected.begin(), selected.end(),
			  [](const PsmRow &a, const PsmRow &b)
			  { return a.order < b.order; });
	return selected;
}

std::unordered_map<std::string, fs::path> collectFt2Files(const std::string &ft2Path)
{
	const fs::path path(ft2Path);
	if (!fs::exists(path))
	{
		throw std::runtime_error("FT2 path does not exist: " + ft2Path);
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
			const std::string ext = toLower(entry.path().extension().string());
			if (ext == ".ft2")
			{
				files.push_back(entry.path());
			}
		}
	}
	std::sort(files.begin(), files.end());
	if (files.empty())
	{
		throw std::runtime_error("No FT2 files found under: " + ft2Path);
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

bool parseFt2PeakLine(const std::string &line, ObservedPeak &peak)
{
	std::istringstream iss(line);
	if (!(iss >> peak.mz >> peak.intensity))
	{
		return false;
	}

	std::string resolution;
	std::string baseline;
	std::string noise;
	std::string chargeText;
	int rawCharge = 1;
	if (iss >> resolution >> baseline >> noise >> chargeText)
	{
		parseIntField(chargeText, rawCharge);
	}
	peak.charge = 1;
	if (rawCharge == 1)
	{
		peak.charge = 1;
	}
	else if (rawCharge <= 0 || rawCharge > 1)
	{
		peak.charge = 1;
	}
	return true;
}

bool parseFt2ScanLine(const std::string &line, int &scanNumber)
{
	std::istringstream iss(line);
	char tag = '\0';
	return static_cast<bool>(iss >> tag >> scanNumber) && tag == 'S' && scanNumber > 0;
}

std::unordered_map<int, ObservedScan> readRequestedScans(const fs::path &ft2File,
														 const std::unordered_set<int> &requestedScans)
{
	std::unordered_map<int, ObservedScan> scans;
	if (requestedScans.empty())
	{
		return scans;
	}

	std::ifstream in;
	std::vector<char> buffer(2 * 1024 * 1024);
	in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
	in.open(ft2File);
	if (!in)
	{
		throw std::runtime_error("Cannot open FT2 file: " + ft2File.string());
	}

	bool capture = false;
	int currentScanNumber = 0;
	ObservedScan currentScan;
	std::string line;

	const auto flushScan = [&]()
	{
		if (capture && currentScanNumber > 0)
		{
			std::sort(currentScan.peaks.begin(), currentScan.peaks.end(),
					  [](const ObservedPeak &a, const ObservedPeak &b)
					  { return a.mz < b.mz; });
			scans[currentScanNumber] = std::move(currentScan);
			currentScan = ObservedScan();
		}
	};

	while (std::getline(in, line))
	{
		if (line.empty())
		{
			continue;
		}
		const char first = line[0];
		if (first == 'S')
		{
			flushScan();
			int scanNumber = 0;
			if (parseFt2ScanLine(line, scanNumber))
			{
				currentScanNumber = scanNumber;
				capture = requestedScans.find(scanNumber) != requestedScans.end();
				currentScan.peaks.clear();
			}
			else
			{
				currentScanNumber = 0;
				capture = false;
				currentScan.peaks.clear();
			}
		}
		else if (capture && first >= '0' && first <= '9')
		{
			ObservedPeak peak;
			if (parseFt2PeakLine(line, peak))
			{
				currentScan.peaks.push_back(peak);
			}
		}
	}
	flushScan();
	return scans;
}

std::vector<Ft2SampleTask> buildFt2SampleTasks(const std::vector<PsmRow> &rows,
											   const std::unordered_map<std::string, fs::path> &ft2Files,
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

	std::vector<Ft2SampleTask> tasks;
	tasks.reserve(groupedBySample.size());
	for (auto &kv : groupedBySample)
	{
		const auto ft2It = ft2Files.find(kv.first);
		if (ft2It == ft2Files.end())
		{
			std::cerr << "No FT2 file found for sample: " << kv.first << "\n";
			processingStats.missingFt2.fetch_add(kv.second.rowIndices.size(), std::memory_order_relaxed);
			continue;
		}

		Ft2SampleTask task;
		task.sample = kv.first;
		task.path = ft2It->second;
		task.requestedScans = std::move(kv.second.requestedScans);
		task.rowIndices = std::move(kv.second.rowIndices);
		tasks.push_back(std::move(task));
	}
	return tasks;
}

int atomIndex(char sipAtom)
{
	switch (sipAtom)
	{
	case 'C':
		return 0;
	case 'H':
		return 1;
	case 'O':
		return 2;
	case 'N':
		return 3;
	case 'P':
		return 4;
	case 'S':
		return 5;
	default:
		return -1;
	}
}

void refreshResidueDistributions(Isotopologue &iso)
{
	for (const auto &kv : iso.mResidueAtomicComposition)
	{
		IsotopeDistribution dist;
		iso.computeIsotopicDistribution(kv.second, dist);
		iso.vResidueIsotopicDistribution[kv.first] = dist;
	}
}

void ensureDefaultNTermAcetylation(Isotopologue &iso)
{
	if (iso.mResidueAtomicComposition.find("%") != iso.mResidueAtomicComposition.end())
	{
		return;
	}
	// Element_List is CHONPS in Sipros configs; % is acetylation C2H2O.
	iso.mResidueAtomicComposition["%"] = {2, 2, 1, 0, 0, 0};
}

int resolveSipIsotopeIndex(const Isotopologue &iso, char sipAtom, int isotopeMassNumber)
{
	const int ix = atomIndex(sipAtom);
	if (ix < 0)
	{
		throw std::runtime_error("Unsupported SIP atom. Use one of C,H,O,N,P,S.");
	}
	if (ix >= static_cast<int>(iso.vAtomIsotopicDistribution.size()) ||
		iso.vAtomIsotopicDistribution[ix].vProb.size() < 2 ||
		iso.vAtomIsotopicDistribution[ix].vMass.size() < 2)
	{
		throw std::runtime_error("Configured isotopic distribution for SIP atom is not usable.");
	}

	const int requestedMass = isotopeMassNumber > 0 ? isotopeMassNumber : 0;
	if (requestedMass <= 0)
	{
		return 1;
	}

	const auto &masses = iso.vAtomIsotopicDistribution[ix].vMass;
	for (size_t i = 1; i < masses.size(); ++i)
	{
		const int massInt = static_cast<int>(std::lround(masses[i]));
		if (massInt == requestedMass)
		{
			return static_cast<int>(i);
		}
	}
	throw std::runtime_error("Requested SIP isotope mass is not available for atom " + std::string(1, sipAtom));
}

void setSipAbundance(Isotopologue &iso, char sipAtom, int isotopeIndex, double sipPct)
{
	const int ix = atomIndex(sipAtom);
	if (ix < 0)
	{
		throw std::runtime_error("Unsupported SIP atom. Use one of C,H,O,N,P,S.");
	}
	if (ix >= static_cast<int>(iso.vAtomIsotopicDistribution.size()) ||
		isotopeIndex < 1 ||
		isotopeIndex >= static_cast<int>(iso.vAtomIsotopicDistribution[ix].vProb.size()))
	{
		throw std::runtime_error("Configured isotopic distribution for SIP atom/isotope is not usable.");
	}

	const int expectedIsotopeIndex = (sipAtom == 'O' || sipAtom == 'S') ? 2 : 1;
	if (isotopeIndex != expectedIsotopeIndex)
	{
		throw std::runtime_error("Requested SIP isotope is inconsistent with atom-specific SIP rule.");
	}

	auto &probs = iso.vAtomIsotopicDistribution[ix].vProb;
	if (!averagine::changeAtomProbability(probs, sipAtom, sipPct / 100.0))
	{
		throw std::runtime_error("Configured isotopic distribution for SIP atom/isotope is not usable.");
	}
	refreshResidueDistributions(iso);
}

void buildPrecursorChargePeaks(const IsotopeDistribution &dist,
							   int charge,
							   double probCutoff,
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

bool buildPrecursorDistributionFromProductIons(Isotopologue &iso,
											   const std::vector<std::vector<double>> &yMass,
											   const std::vector<std::vector<double>> &yProb,
											   const std::vector<std::vector<double>> &bMass,
											   const std::vector<std::vector<double>> &bProb,
											   IsotopeDistribution &precursorDist)
{
	if (bMass.empty() || bProb.empty() || yMass.empty() || yProb.empty())
	{
		return false;
	}
	const size_t bLast = bMass.size() - 1;
	const size_t yFirst = 0;
	IsotopeDistribution bLastDist(bMass[bLast], bProb[bLast]);
	IsotopeDistribution y1Dist(yMass[yFirst], yProb[yFirst]);
	precursorDist = iso.sum(bLastDist, y1Dist);
	return true;
}

bool findMatchedPeak(const std::vector<ObservedPeak> &peaks,
					 double targetMz,
					 double ppmTolerance,
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

			const bool matchedEnvelope = findMatchedPeak(scan.peaks, envelope[apexIndex].mz, ppmTolerance, matched);
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

void matchBaselineInFt2Batches(const std::vector<PsmRow> &rows,
							   const std::unordered_map<std::string, fs::path> &ft2Files,
							   const Isotopologue &baselineIso,
							   const Args &args,
							   int effectiveThreads,
							   std::vector<MatchedEnvelopeSet> &matchedEnvelopeSets,
							   std::vector<char> &baselineOk,
							   ProcessingStats &processingStats)
{
	const std::vector<Ft2SampleTask> tasks = buildFt2SampleTasks(rows, ft2Files, processingStats);
	const size_t batchSize = static_cast<size_t>(std::max(1, effectiveThreads));

	for (size_t batchStart = 0; batchStart < tasks.size(); batchStart += batchSize)
	{
		const size_t batchEnd = std::min(batchStart + batchSize, tasks.size());
		const size_t batchCount = batchEnd - batchStart;
		std::vector<std::unordered_map<int, ObservedScan>> loadedScans(batchCount);
		std::vector<std::string> errors(batchCount);

		{
			std::cout << "Timing parallel region: load FT2 batch "
					  << (batchStart / batchSize + 1) << " of "
					  << ((tasks.size() + batchSize - 1) / batchSize)
					  << " (" << batchCount << " FT2 files)" << std::endl;
			CLOCKSTART;
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < static_cast<int>(batchCount); ++i)
			{
				const size_t localTaskIndex = static_cast<size_t>(i);
				const Ft2SampleTask &task = tasks[batchStart + localTaskIndex];
#pragma omp critical(cerr_output)
				{
					std::cerr << "Reading " << task.requestedScans.size()
							  << " requested scans from " << task.path << "\n";
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
			CLOCKSTOP;
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
			const Ft2SampleTask &task = tasks[batchStart + localTaskIndex];
			for (size_t rowIndex : task.rowIndices)
			{
				batchRows.push_back({localTaskIndex, rowIndex});
			}
		}

		{
			std::cout << "Timing parallel region: baseline match FT2 batch "
					  << (batchStart / batchSize + 1) << " of "
					  << ((tasks.size() + batchSize - 1) / batchSize)
					  << " (" << batchRows.size() << " PSM rows)" << std::endl;
			CLOCKSTART;
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
			CLOCKSTOP;
		}
	}
}

bool buildShiftedChargeOneFragmentEntries(const std::vector<std::vector<double>> &bMass,
										  const std::vector<std::vector<double>> &bProb,
										  const std::vector<std::vector<double>> &yMass,
										  const std::vector<std::vector<double>> &yProb,
										  const MatchedEnvelopeSet &matchedSet,
										  double probCutoff,
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
						  const std::unordered_set<std::string> &experimentalPeptides,
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
		if (candidate != targetPeptide && experimentalPeptides.find(candidate) == experimentalPeptides.end())
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

fs::path appendHdf5ExtensionIfNeeded(const fs::path &path)
{
	if (toLower(path.extension().string()) == ".h5")
	{
		return path;
	}
	fs::path hdf5Path = path;
	hdf5Path += ".h5";
	return hdf5Path;
}

fs::path resolveOutputBasePath(const std::string &outputPath, bool multipleOutputFiles)
{
	const fs::path path(outputPath);
	if (hasTrailingPathSeparator(outputPath) ||
		(fs::exists(path) && fs::is_directory(path)))
	{
		return path / "spectra.h5";
	}

	const fs::path hdf5Path = appendHdf5ExtensionIfNeeded(path);
	if (!multipleOutputFiles)
	{
		return hdf5Path;
	}

	if (toLower(path.extension().string()) == ".h5")
	{
		throw std::runtime_error("Multiple output files requested, but -o names a single .h5 file. Use an output directory path, for example -o out or -o out/.");
	}

	return path / "spectra.h5";
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
	writeIntAttribute(file, "format_version", 1);
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

bool writeGeneratedSpectraFile(const fs::path &path,
							   const std::vector<SpectrumOutputRecord> &records,
							   const std::vector<char> &ok,
							   const Hdf5OutputMetadata &metadata,
							   size_t &written)
{
	Hdf5OutputData data = buildOutputDataFromRecords(records, ok);
	written = data.recordCount();
	return writeSpectraHdf5File(path, data, metadata);
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
								   OutputJobStats &jobStats)
{
	try
	{
		Isotopologue localIso = pristineIso;
		setSipAbundance(localIso, sipAtom, targetSipIsotopeIndex, job.metadata.targetSipAbundancePct);

		std::vector<SpectrumOutputRecord> records(rows.size());
		std::vector<char> ok(rows.size(), 0);

		for (size_t rowIndex = 0; rowIndex < rows.size(); ++rowIndex)
		{
			if (!baselineOk[rowIndex])
			{
				continue;
			}

			const PsmRow &row = rows[rowIndex];
			if (!job.decoy)
			{
				std::vector<std::vector<double>> yMass, yProb, bMass, bProb;
				if (!localIso.computeProductIon(row.peptide, yMass, yProb, bMass, bProb))
				{
					++jobStats.targetFailed;
					continue;
				}

				IsotopeDistribution precursorDist;
				if (!buildPrecursorDistributionFromProductIons(localIso, yMass, yProb, bMass, bProb, precursorDist))
				{
					++jobStats.targetFailed;
					continue;
				}

				SpectrumOutputRecord record;
				record.psmId = spectraHeaderId(row.psmId, job.metadata.targetSipAbundancePct, false);
				record.retention = row.retentionText;
				record.charge = row.precursorCharge;
				record.peptide = row.peptide;
				buildPrecursorChargePeaks(precursorDist, row.precursorCharge, args.probCutoff,
										  record.precursorMz, record.precursorIntensity);

				std::vector<FragmentEntry> entries;
				if (!buildShiftedChargeOneFragmentEntries(bMass, bProb, yMass, yProb,
														  matchedEnvelopeSets[rowIndex], args.probCutoff,
														  entries))
				{
					++jobStats.targetFailed;
					continue;
				}

				copyFragmentEntriesToRecord(entries, record);
				records[rowIndex] = std::move(record);
				ok[rowIndex] = 1;
				continue;
			}

			if (decoyPeptides[rowIndex].empty())
			{
				continue;
			}

			const std::vector<double> apexPool =
				shuffledApexIntensityPool(matchedEnvelopeSets[rowIndex], args.decoySeed, rowIndex);
			const double meanApex = meanIntensity(apexPool);
			if (apexPool.empty() || meanApex <= 0.0)
			{
				++jobStats.decoyComputeFailed;
				continue;
			}

			std::vector<std::vector<double>> decoyYMass, decoyYProb, decoyBMass, decoyBProb;
			if (!localIso.computeProductIon(decoyPeptides[rowIndex], decoyYMass, decoyYProb, decoyBMass, decoyBProb))
			{
				++jobStats.decoyComputeFailed;
				continue;
			}

			IsotopeDistribution decoyPrecursorDist;
			if (!buildPrecursorDistributionFromProductIons(localIso, decoyYMass, decoyYProb,
															decoyBMass, decoyBProb, decoyPrecursorDist))
			{
				++jobStats.decoyComputeFailed;
				continue;
			}

			SpectrumOutputRecord decoyRecord;
			decoyRecord.psmId = spectraHeaderId(row.psmId, job.metadata.targetSipAbundancePct, true);
			decoyRecord.retention = row.retentionText;
			decoyRecord.charge = row.precursorCharge;
			decoyRecord.peptide = decoyPeptides[rowIndex];
			buildPrecursorChargePeaks(decoyPrecursorDist, row.precursorCharge, args.probCutoff,
									  decoyRecord.precursorMz, decoyRecord.precursorIntensity);

			std::vector<FragmentEntry> decoyEntries;
			if (!buildDecoyChargeOneFragmentEntries(decoyBMass, decoyBProb,
													decoyYMass, decoyYProb,
													apexPool, meanApex,
													decoyAddedResidues[rowIndex] != 0,
													args.probCutoff,
													decoyEntries))
			{
				++jobStats.decoyComputeFailed;
				continue;
			}

			copyFragmentEntriesToRecord(decoyEntries, decoyRecord);
			records[rowIndex] = std::move(decoyRecord);
			ok[rowIndex] = 1;
		}

		job.success = writeGeneratedSpectraFile(job.path, records, ok, job.metadata, job.written);
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
			generateAndWriteOutputFileJob(outputJobs[jobIndex], rows, baselineOk, matchedEnvelopeSets,
										  decoyPeptides, decoyAddedResidues, pristineIso,
										  args, sipAtom, targetSipIsotopeIndex, childStats);
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
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < static_cast<int>(outputJobs.size()); ++i)
	{
		generateAndWriteOutputFileJob(outputJobs[static_cast<size_t>(i)], rows, baselineOk,
									  matchedEnvelopeSets, decoyPeptides, decoyAddedResidues,
									  pristineIso, args, sipAtom, targetSipIsotopeIndex,
									  jobStats[static_cast<size_t>(i)]);
	}
	for (const OutputJobStats &stats : jobStats)
	{
		addOutputJobStats(processingStats, stats);
	}
	return true;
#endif
}

} // namespace

int main(int argc, char **argv)
{
	std::setlocale(LC_ALL, "C");
	std::ios_base::sync_with_stdio(false);

	Args args;
	if (!parseArgs(argc, argv, args))
	{
		return 1;
	}

	if (!ProNovoConfig::setFilename(args.configPath))
	{
		std::cerr << "Could not load config file: " << args.configPath << "\n";
		return 1;
	}

	char sipAtom = args.sipAtom;
	int sipIsotopeMassNumber = args.sipIsotopeMassNumber;
	if (sipAtom == '\0')
	{
		const std::string cfgAtom = ProNovoConfig::getSetSIPelement();
		if (cfgAtom.size() != 1)
		{
			std::cerr << "Invalid SIP_Element in config. Pass -a explicitly.\n";
			return 1;
		}
		sipAtom = static_cast<char>(std::toupper(static_cast<unsigned char>(cfgAtom[0])));
		if (sipAtom == 'O')
		{
			sipIsotopeMassNumber = 18;
		}
		else if (sipAtom == 'S')
		{
			sipIsotopeMassNumber = 34;
		}
	}

	int targetSipIsotopeIndex = 1;
	try
	{
		targetSipIsotopeIndex = resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, sipAtom, sipIsotopeMassNumber);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	int baselineSipIsotopeIndex = 1;
	try
	{
		baselineSipIsotopeIndex = resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, 'C', 13);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Failed to resolve baseline C13 isotope: " << ex.what() << "\n";
		return 1;
	}

	Isotopologue pristineIso = ProNovoConfig::configIsotopologue;
	try
	{
		ensureDefaultNTermAcetylation(pristineIso);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Failed to initialize isotopologue configuration: " << ex.what() << "\n";
		return 1;
	}

	Isotopologue baselineIso = pristineIso;
	try
	{
		setSipAbundance(baselineIso, 'C', baselineSipIsotopeIndex, 1.0);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Failed to apply baseline 1% C13 abundance: " << ex.what() << "\n";
		return 1;
	}

	const std::vector<double> targetAbundances = makeTargetAbundances(args);

	if (args.threads > 0)
	{
		omp_set_num_threads(args.threads);
	}
	const int effectiveThreads = std::max(1, omp_get_max_threads());

	ReadStats readStats;
	std::vector<PsmRow> allRows;
	try
	{
		allRows = readInputRows(args.inputPath, readStats);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	std::unordered_set<std::string> experimentalPeptides;
	experimentalPeptides.reserve(allRows.size());
	for (const PsmRow &row : allRows)
	{
		experimentalPeptides.insert(row.peptide);
	}

	std::vector<PsmRow> rows = selectBestRowsByPeptideCharge(allRows);
	std::cerr << "Read " << readStats.totalRows << " PSM rows; retained " << rows.size()
			  << " highest-probability peptide/charge rows";
	if (readStats.shortRows > 0 || readStats.invalidRows > 0 || readStats.unsupportedMods > 0)
	{
		std::cerr << " (skipped short=" << readStats.shortRows
				  << ", invalid=" << readStats.invalidRows
				  << ", unsupported_mods=" << readStats.unsupportedMods << ")";
	}
	std::cerr << "\n";

	std::unordered_map<std::string, fs::path> ft2Files;
	try
	{
		ft2Files = collectFt2Files(args.ft2Path);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	std::vector<MatchedEnvelopeSet> matchedEnvelopeSets(rows.size());
	std::vector<char> baselineOk(rows.size(), 0);
	ProcessingStats processingStats;

	try
	{
		matchBaselineInFt2Batches(rows, ft2Files, baselineIso, args, effectiveThreads,
								  matchedEnvelopeSets, baselineOk, processingStats);
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
			std::cout << "Timing parallel region: generate decoy peptides ("
					  << rows.size() << " PSM rows)" << std::endl;
			CLOCKSTART;
#pragma omp parallel for schedule(dynamic) reduction(+ : generatedCount)
			for (int i = 0; i < static_cast<int>(rows.size()); ++i)
			{
				const size_t rowIndex = static_cast<size_t>(i);
				if (!baselineOk[rowIndex])
				{
					continue;
				}
				bool addedResidue = false;
				if (generateDecoyPeptide(rows[rowIndex].peptide, experimentalPeptides,
										 args.decoySeed, rowIndex, decoyPeptides[rowIndex],
										 addedResidue))
				{
					decoyAddedResidues[rowIndex] = addedResidue ? 1 : 0;
					++generatedCount;
				}
				else
				{
					++processingStats.decoyCollision;
				}
			}
			CLOCKSTOP;
		}
		decoyPeptidesGenerated = generatedCount;
	}

	const size_t plannedOutputFileCount = targetAbundances.size() * (args.writeDecoy ? 2 : 1);
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
	const std::string outputSipLabel = sipLabelForPath(sipAtom, sipIsotopeMassNumber);
	size_t outputFilesWritten = 0;
	size_t decoyFilesWritten = 0;
	std::vector<OutputFileJob> outputJobs;
	outputJobs.reserve(targetAbundances.size() * (args.writeDecoy ? 2 : 1));
	for (double targetAbundancePct : targetAbundances)
	{
		Hdf5OutputMetadata outputMetadata;
		outputMetadata.recordKind = "target";
		outputMetadata.targetSipAbundancePct = targetAbundancePct;
		outputMetadata.sipAtom = sipAtom;
		outputMetadata.sipIsotopeMassNumber = sipIsotopeMassNumber;
		outputMetadata.probCutoff = args.probCutoff;
		outputMetadata.ppmTolerance = args.ppmTolerance;
		outputMetadata.minMatchedEnvelopes = static_cast<uint64_t>(args.minMatchedEnvelopes);

		OutputFileJob targetJob;
		targetJob.path = multipleOutputFiles
							 ? outputPathForAbundance(outputBasePath, outputSipLabel, targetAbundancePct)
							 : outputBasePath;
		targetJob.metadata = outputMetadata;
		targetJob.decoy = false;
		outputJobs.push_back(std::move(targetJob));

		if (args.writeDecoy)
		{
			Hdf5OutputMetadata decoyMetadata = outputMetadata;
			decoyMetadata.recordKind = "decoy";
			OutputFileJob decoyJob;
			decoyJob.path = decoyOutputPathForAbundance(outputBasePath, outputSipLabel, targetAbundancePct);
			decoyJob.metadata = decoyMetadata;
			decoyJob.decoy = true;
			outputJobs.push_back(std::move(decoyJob));
		}
	}

	{
		std::cout << "Timing parallel region: process output HDF5 files with fork pool ("
				  << outputJobs.size() << " files, up to " << effectiveThreads
				  << " processes)" << std::endl;
		CLOCKSTART;
		if (!runOutputJobs(outputJobs, rows, baselineOk, matchedEnvelopeSets,
						   decoyPeptides, decoyAddedResidues, pristineIso, args,
						   sipAtom, targetSipIsotopeIndex, effectiveThreads, processingStats))
		{
			return 1;
		}
		CLOCKSTOP;
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

		if (job.decoy)
		{
			++decoyFilesWritten;
			std::cerr << "Wrote decoy shifted spectra for " << job.written
					  << " PSMs at " << job.metadata.targetSipAbundancePct
					  << "% to " << job.path << "\n";
		}
		else
		{
			++outputFilesWritten;
			std::cerr << "Wrote shifted matched experimental spectra for " << job.written
					  << " PSMs at " << job.metadata.targetSipAbundancePct
					  << "% to " << job.path << "\n";
		}
	}

	if (!allOutputSucceeded)
	{
		return 1;
	}

	std::cerr << "Input peptide spectra: " << rows.size()
			  << "; baseline matched peptide spectra extracted: " << baselineMatchedRows
			  << "; output files written: " << outputFilesWritten
			  << "; minimum matched envelopes: " << args.minMatchedEnvelopes;
	if (args.writeDecoy)
	{
		std::cerr << "; decoy peptides generated: " << decoyPeptidesGenerated
				  << "; decoy files written: " << decoyFilesWritten;
	}
	std::cerr << "\n";
	std::cerr << "Skipped during processing: missing_ft2=" << processingStats.missingFt2.load()
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
	return 0;
}
