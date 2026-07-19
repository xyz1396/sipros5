#include "SiprosSearchRunner.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void require(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

bool parse(const std::vector<std::string> &values,
		   sipros::DatabaseSearchArguments &arguments,
		   std::string &error)
{
	std::vector<std::string> storage = values;
	std::vector<char *> argv;
	argv.reserve(storage.size());
	for (std::string &value : storage)
	{
		argv.push_back(value.data());
	}
	std::ostringstream output;
	std::ostringstream errors;
	sipros::SiprosSearchRunner runner;
	const bool result = runner.initializeArguments(
		static_cast<int>(argv.size()), argv.data(), arguments, output, errors);
	error = errors.str();
	return result;
}

} // namespace

int main()
{
	try
	{
		std::string error;
		sipros::DatabaseSearchArguments regular;
		require(parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "-f", "two.hdf5", "-o", "out"},
			regular, error), error);
		require(regular.scanFiles ==
			std::vector<std::string>({"one.h5", "two.hdf5"}),
			"repeatable -f did not preserve every input in order");
		require(regular.raxportReadOptions.precursorSource ==
			sipros::PrecursorSource::Ms1Neighborhood,
			"Regular search did not resolve to MS1-neighborhood precursors");

		sipros::DatabaseSearchArguments removedLegacyFlag;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "--indexed-search"}, removedLegacyFlag, error) &&
			error.find("Unknown option --indexed-search") != std::string::npos,
			"removed indexed-search compatibility flag was still accepted");

		sipros::DatabaseSearchArguments prepareOnly;
		require(parse({"sipros search-fasta", "-fasta", "target.fasta",
			"--fragment-index-cache", "target.sfi", "--prepare-only"},
			prepareOnly, error), error);
		require(prepareOnly.prepareOnly,
			"--prepare-only did not select cache preparation mode");
		require(prepareOnly.scanFiles.empty() &&
			prepareOnly.workingDirectory.empty(),
			"cache preparation unexpectedly discovered scan files");

		sipros::DatabaseSearchArguments prepareWithoutCache;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"--prepare-only"}, prepareWithoutCache, error) &&
			error.find("requires --fragment-index-cache") != std::string::npos,
			"cache preparation without a persistent cache was not rejected");

		sipros::DatabaseSearchArguments prepareWithScan;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "--fragment-index-cache", "target.sfi",
			"--prepare-only"}, prepareWithScan, error) &&
			error.find("does not accept -f or -w") != std::string::npos,
			"cache preparation with a scan input was not rejected");

		sipros::DatabaseSearchArguments sip;
		require(parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "-a", "C13", "-b", "1-5", "-s", "1"},
			sip, error), error);
		require(sip.raxportReadOptions.precursorSource ==
			sipros::PrecursorSource::RaxportCandidates,
			"SIP search did not retain Raxport precursor candidates");

		sipros::DatabaseSearchArguments badSipIndex;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "-a", "C13", "-b", "1-5", "-s", "1",
			"--fragment-index-cache", "sip.sfi"}, badSipIndex, error) &&
			error.find("legacy H5 precursor candidates") != std::string::npos,
			"SIP fragment-index conflict was not rejected clearly");

		sipros::DatabaseSearchArguments badRegularSource;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "--precursor-source", "raxport-candidates"},
			badRegularSource, error) &&
			error.find("requires --precursor-source ms1-neighborhood") !=
				std::string::npos,
			"Regular legacy precursor-source conflict was not rejected");

		sipros::DatabaseSearchArguments badPinOutput;
		require(!parse({"sipros search-fasta", "-fasta", "target.fasta",
			"-f", "one.h5", "-f", "two.h5", "--pin-output", "one.pin"},
			badPinOutput, error) &&
			error.find("exactly one -f") != std::string::npos,
			"multi-file --pin-output conflict was not rejected");

		return 0;
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
