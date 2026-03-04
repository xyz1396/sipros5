#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "averagine.h"
#include "proNovoConfig.h"

namespace
{
struct Args
{
	std::string configPath;
	std::string inputPath;
	std::string outputPath;
	char sipAtom = '\0';
	int sipIsotopeMassNumber = -1;
	double probCutoff = 0.01;
	int threads = 0;
};

struct PsmRow
{
	std::string psmId;
	double sipPct = 0.0;
	std::string peptide;
	int precursorCharge = 1;
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
			  << " -c <config.cfg> -i <input.tsv|input.pin> -o <output.txt> [-a <SIP atom/isotope, e.g. C13,H2,O18,N15,S34>] [-p <prob cutoff>] [-t <threads>]\n";
	std::cerr << "Required columns in input: (PSMId or SpecId), MS2IsotopicAbundances, Peptide\n";
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
		// Use canonical SIP isotopes for O and S when isotope mass is omitted.
		if (atom == 'O')
		{
			sipIsotopeMassNumber = 18;
		}
		else if (atom == 'S')
		{
			sipIsotopeMassNumber = 34;
		}
		else
		{
			sipIsotopeMassNumber = -1;
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
		if (i + 1 >= argc)
		{
			std::cerr << "Missing value for option: " << opt << "\n";
			return false;
		}
		if (opt == "-c")
		{
			args.configPath = argv[++i];
		}
		else if (opt == "-i")
		{
			args.inputPath = argv[++i];
		}
		else if (opt == "-o")
		{
			args.outputPath = argv[++i];
		}
		else if (opt == "-a" || opt == "--sip-atom")
		{
			const std::string sipSpec = trim(argv[++i]);
			if (!parseSipAtomSpec(sipSpec, args.sipAtom, args.sipIsotopeMassNumber))
			{
				std::cerr << "Invalid SIP atom/isotope: " << sipSpec
						  << ". Use C13,H2,O18,N15,S34 (or C/H/O/N/P/S).\n";
				return false;
			}
		}
		else if (opt == "-p" || opt == "--prob-cutoff")
		{
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
		else if (opt == "-t" || opt == "--threads")
		{
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

	if (args.configPath.empty() || args.inputPath.empty() || args.outputPath.empty())
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
	const std::string t = trim(s);
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

std::string normalizePeptide(const std::string &raw)
{
	const std::string p = trim(raw);
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

std::vector<PsmRow> readInputRows(const std::string &path)
{
	std::ifstream in(path);
	if (!in)
	{
		throw std::runtime_error("Cannot open input file: " + path);
	}

	std::string headerLine;
	if (!std::getline(in, headerLine))
	{
		throw std::runtime_error("Input file is empty: " + path);
	}
	const std::vector<std::string> headers = splitTab(headerLine);
	std::unordered_map<std::string, size_t> col;
	for (size_t i = 0; i < headers.size(); ++i)
	{
		col[trim(headers[i])] = i;
	}

	const auto getRequired = [&](const std::string &name) -> size_t
	{
		const auto it = col.find(name);
		if (it == col.end())
		{
			throw std::runtime_error("Missing required column: " + name);
		}
		return it->second;
	};

	size_t idxPSMId = std::string::npos;
	auto itPSMId = col.find("PSMId");
	if (itPSMId != col.end())
	{
		idxPSMId = itPSMId->second;
	}
	else
	{
		auto itSpecId = col.find("SpecId");
		if (itSpecId != col.end())
		{
			idxPSMId = itSpecId->second;
		}
	}
	if (idxPSMId == std::string::npos)
	{
		throw std::runtime_error("Missing required ID column: PSMId or SpecId");
	}
	const size_t idxMS2Pct = getRequired("MS2IsotopicAbundances");
	const size_t idxPeptide = getRequired("Peptide");
	size_t idxCharge = std::string::npos;
	const auto itCharge = col.find("parentCharges");
	if (itCharge != col.end())
	{
		idxCharge = itCharge->second;
	}

	std::vector<PsmRow> rows;
	std::string line;
	size_t lineNo = 1;
	while (std::getline(in, line))
	{
		++lineNo;
		if (line.empty())
		{
			continue;
		}
		const std::vector<std::string> f = splitTab(line);
		const size_t need = std::max({idxPSMId, idxMS2Pct, idxPeptide});
		if (f.size() <= need)
		{
			continue;
		}
		PsmRow row;
		row.psmId = trim(f[idxPSMId]);
		row.peptide = normalizePeptide(f[idxPeptide]);
		if (row.psmId.empty() || row.peptide.empty())
		{
			continue;
		}
			try
			{
				row.sipPct = parseFirstDouble(f[idxMS2Pct]);
				if (idxCharge != std::string::npos)
				{
					row.precursorCharge = std::max(1, std::stoi(trim(f[idxCharge])));
				}
			}
			catch (const std::exception &)
			{
				std::cerr << "Skipping line " << lineNo << ": invalid MS2IsotopicAbundances or parentCharges.\n";
				continue;
			}
		rows.push_back(std::move(row));
	}
	return rows;
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

	std::ostringstream supported;
	bool first = true;
	for (size_t i = 1; i < masses.size(); ++i)
	{
		if (!first)
		{
			supported << ',';
		}
		supported << static_cast<int>(std::lround(masses[i]));
		first = false;
	}
	throw std::runtime_error("Requested SIP isotope mass " + std::to_string(requestedMass) +
							 " is not available for atom " + std::string(1, sipAtom) +
							 ". Supported isotope masses: " + supported.str());
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
	const bool updated = averagine::changeAtomProbability(probs, sipAtom, sipPct / 100.0);
	if (!updated)
	{
		throw std::runtime_error("Configured isotopic distribution for SIP atom/isotope is not usable.");
	}

	refreshResidueDistributions(iso);
}

void appendChargeSeriesLine(std::ostream &out,
							const std::vector<std::vector<double>> &masses,
							const std::vector<std::vector<double>> &probs,
							int charge,
							double probCutoff)
{
	const double proton = ProNovoConfig::getProtonMass();
	std::vector<std::pair<double, double>> peaks;
	peaks.reserve(256);
	for (size_t i = 0; i < masses.size(); ++i)
	{
		for (size_t j = 0; j < masses[i].size() && j < probs[i].size(); ++j)
		{
			if (probs[i][j] < probCutoff)
			{
				continue;
			}
			const double mz = masses[i][j] / static_cast<double>(charge) + proton;
			peaks.emplace_back(mz, probs[i][j]);
		}
	}

	std::sort(peaks.begin(), peaks.end(),
			  [](const std::pair<double, double> &a, const std::pair<double, double> &b)
			  { return a.first < b.first; });

	bool first = true;
	for (const auto &peak : peaks)
	{
		if (!first)
		{
			out << ' ';
		}
		out << std::setprecision(10) << peak.first << ' ' << std::setprecision(10) << peak.second;
		first = false;
	}
	out << '\n';
}

void appendPrecursorHeaderLine(std::ostream &out,
							   const std::string &psmId,
							   const IsotopeDistribution &dist,
							   int charge,
							   double probCutoff)
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

	out << psmId;
	for (const auto &peak : peaks)
	{
		out << ' ' << std::setprecision(10) << peak.first << ' ' << std::setprecision(10) << peak.second;
	}
	out << '\n';
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
} // namespace

int main(int argc, char **argv)
{
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

	if (atomIndex(sipAtom) < 0)
	{
		std::cerr << "Invalid SIP atom '" << sipAtom << "'. Valid options: C,H,O,N,P,S\n";
		return 1;
	}
	int sipIsotopeIndex = 1;
	try
	{
		sipIsotopeIndex = resolveSipIsotopeIndex(ProNovoConfig::configIsotopologue, sipAtom, sipIsotopeMassNumber);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}
	std::cout << "SIP atom: " << sipAtom;
	if (sipIsotopeMassNumber > 0)
	{
		std::cout << sipIsotopeMassNumber;
	}
	std::cout << "\n";

	if (args.threads > 0)
	{
		omp_set_num_threads(args.threads);
	}

	std::vector<PsmRow> rows;
	try
	{
		rows = readInputRows(args.inputPath);
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		return 1;
	}

	std::ofstream out(args.outputPath);
	if (!out)
	{
		std::cerr << "Cannot open output file: " << args.outputPath << "\n";
		return 1;
	}

	const Isotopologue baseIso = ProNovoConfig::configIsotopologue;
	std::vector<std::string> blocks(rows.size());
	std::vector<char> ok(rows.size(), 0);
	size_t written = 0;

#pragma omp parallel
	{
		Isotopologue localIso = baseIso;
#pragma omp for schedule(dynamic)
		for (int i = 0; i < static_cast<int>(rows.size()); ++i)
		{
			const auto &row = rows[static_cast<size_t>(i)];
			try
			{
				setSipAbundance(localIso, sipAtom, sipIsotopeIndex, row.sipPct);
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
				if (!buildPrecursorDistributionFromProductIons(localIso, yMass, yProb, bMass, bProb, precursorDist))
				{
#pragma omp critical
					{
						std::cerr << "Skipping PSM " << row.psmId << ": precursor reconstruction from product ions failed.\n";
					}
					continue;
				}

				std::ostringstream ss;
				appendPrecursorHeaderLine(ss, row.psmId, precursorDist, row.precursorCharge, args.probCutoff);
				appendChargeSeriesLine(ss, bMass, bProb, 1, args.probCutoff);
				appendChargeSeriesLine(ss, bMass, bProb, 2, args.probCutoff);
				appendChargeSeriesLine(ss, bMass, bProb, 3, args.probCutoff);
				appendChargeSeriesLine(ss, yMass, yProb, 1, args.probCutoff);
				appendChargeSeriesLine(ss, yMass, yProb, 2, args.probCutoff);
				appendChargeSeriesLine(ss, yMass, yProb, 3, args.probCutoff);

				blocks[static_cast<size_t>(i)] = ss.str();
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

	for (size_t i = 0; i < blocks.size(); ++i)
	{
		if (!ok[i])
		{
			continue;
		}
		out << blocks[i];
		++written;
	}

	std::cerr << "Wrote theoretical spectra for " << written << " PSMs to " << args.outputPath << "\n";
	return 0;
}
