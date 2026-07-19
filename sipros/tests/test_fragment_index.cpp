#include "MVH.h"
#include "fragmentindex.h"
#include "peptide.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace
{

void check(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

} // namespace

int main()
{
	try
	{
		const unsigned long processId =
#if defined(__unix__) || defined(__APPLE__)
			static_cast<unsigned long>(getpid());
#else
			0;
#endif
		const fs::path fasta = fs::temp_directory_path() /
			("sipros_fragment_index_test_" + std::to_string(processId) + ".fasta");
		const fs::path cache = fs::temp_directory_path() /
			("sipros_fragment_index_test_" + std::to_string(processId) + ".sfi");
		{
			std::ofstream out(fasta);
			out << ">protein_one\n"
				<< "MPEPTIDERKACDEFGHIKLMNPQRSTVWYRK\n"
				<< ">protein_two\n"
				<< "MKRHQNPEPTIDEKAAAAAAAK\n";
			// Cross several bounded precursor blocks so local posting ids and
			// block directories are exercised by the mmap round trip.
			for (int i = 0; i < 270; ++i)
			{
				out << ">bulk_protein_" << i << "\n"
					<< "MPEPTIDEKAAAAAAAKPEPTIDER\n";
			}
		}

		check(ProNovoConfig::load(ProNovoConfig::Profile::Regular),
			"could not load Regular profile");
		std::string configError;
		check(ProNovoConfig::configureVariablePtms({"none"}, 1, configError),
			"could not disable variable PTMs: " + configError);
		ProNovoConfig::setFASTAfilename(fasta.string());
		struct SerialPeptide
		{
			std::string identified;
			std::string original;
			std::string protein;
			int begin = 0;
			char identifyPrefix = '-';
			char originalPrefix = '-';
		};
		std::vector<SerialPeptide> serialPeptides;
		ProteinDatabase serialDatabase;
		serialDatabase.loadDatabase();
		if (serialDatabase.getFirstProtein())
		{
			Peptide peptide;
			while (serialDatabase.getNextPeptide(&peptide))
			{
				serialPeptides.push_back({
					peptide.sPeptide,
					peptide.sOriginalPeptide,
					peptide.sProteinName,
					peptide.ibeginPos,
					peptide.cIdentifyPrefix,
					peptide.cOriginalPrefix});
			}
		}

		sipros::FragmentIndex built;
		std::string error;
		check(built.loadOrBuild(cache.string(), true, error),
			"could not build fragment index: " + error);
		check(!built.stats().loadedFromCache, "forced build unexpectedly used cache");
		check(built.peptideCount() > 0, "fragment index has no peptides");
		check(built.peptideCount() == serialPeptides.size(),
			"parallel digest peptide count differs from serial ProteinDatabase");
		check(built.fragmentCount() > 0, "fragment index has no postings");
		check(built.fragmentBinCount() > 0 &&
			built.fragmentBinCount() <= built.fragmentCount(),
			"sparse fragment-bin directory has invalid cardinality");
		check(built.precursorBlockCount() > 1,
			"fragment-index test did not cross a bounded precursor block");
		check(fs::file_size(cache) == built.stats().cacheBytes,
			"reported cache size does not match file size");

		for (uint32_t peptideId = 0; peptideId < built.peptideCount(); ++peptideId)
		{
			const auto &record = built.peptide(peptideId);
			check(record.generationOrdinal < serialPeptides.size(),
				"parallel digest generation ordinal is invalid");
			const SerialPeptide &serial =
				serialPeptides[record.generationOrdinal];
			check(record.scoringOffset == record.peptideOffset &&
				record.scoringSize == record.peptideSize,
				"identical scoring sequence was not aliased in compact storage");
			check(built.peptideSequence(peptideId) == serial.identified &&
				built.originalSequence(peptideId) == serial.original &&
				built.proteinNames(peptideId) == serial.protein &&
				record.beginPosition == serial.begin &&
				record.identifyPrefix == serial.identifyPrefix &&
				record.originalPrefix == serial.originalPrefix,
				"parallel digest metadata differs from serial ProteinDatabase");
			std::string sequence(built.scoringSequence(peptideId));
			std::vector<double> charged;
			std::vector<double> forward;
			std::vector<double> reverse;
			std::vector<char> residues;
			check(MVH::CalculateSequenceIons(sequence, 2,
				MVH::bUseSmartPlusThreeModel, &charged, &forward, &reverse,
				&residues), "MVH rejected indexed scoring sequence");
			check(forward.size() == record.peptideLength &&
				reverse.size() == record.peptideLength,
				"universal index omitted a full-length neutral ion");
			const uint32_t block = built.precursorBlockForPeptide(peptideId);
			auto postingExists = [&](double mass)
			{
				const auto range = built.fragmentRange(block, mass, mass);
				for (const sipros::FragmentPosting *posting = range.first;
					 posting != range.second; ++posting)
				{
					if (posting != nullptr &&
						built.postingPeptideId(block, *posting) == peptideId)
					{
						return true;
					}
				}
				return false;
			};
			for (double mass : forward)
			{
				check(postingExists(mass), "indexed b-ion posting is missing");
			}
			for (double mass : reverse)
			{
				check(postingExists(mass), "indexed y-ion posting is missing");
			}

			Peptide materialized;
			built.materializePeptide(peptideId, materialized);
			check(materialized.sPeptide == built.peptideSequence(peptideId),
				"materialized peptide sequence differs from cache record");
			check(materialized.sNeutralLossPeptide == built.scoringSequence(peptideId),
				"materialized scoring sequence differs from cache record");
		}

		sipros::FragmentIndex loaded;
		check(loaded.loadOrBuild(cache.string(), false, error),
			"could not mmap fragment index: " + error);
		check(loaded.stats().loadedFromCache, "warm load rebuilt the cache");
		check(loaded.peptideCount() == built.peptideCount() &&
			loaded.fragmentCount() == built.fragmentCount(),
			"warm cache counts differ from built index");
		check(loaded.peptideSequence(0) == built.peptideSequence(0),
			"warm cache peptide differs from built index");

		{
			std::fstream corrupt(cache,
				std::ios::in | std::ios::out | std::ios::binary);
			check(static_cast<bool>(corrupt), "could not open cache for corruption test");
			corrupt.seekg(-1, std::ios::end);
			char byte = 0;
			corrupt.read(&byte, 1);
			check(static_cast<bool>(corrupt), "could not read cache corruption byte");
			byte ^= static_cast<char>(0x5a);
			corrupt.seekp(-1, std::ios::end);
			corrupt.write(&byte, 1);
			corrupt.close();
			check(static_cast<bool>(corrupt), "could not write cache corruption byte");
		}
		sipros::FragmentIndex repaired;
		check(repaired.loadOrBuild(cache.string(), false, error),
			"could not rebuild corrupted fragment index: " + error);
		check(!repaired.stats().loadedFromCache,
			"corrupted cache passed payload validation");

		std::error_code ignored;
		fs::remove(cache, ignored);
		fs::remove(fasta, ignored);
		std::cout << "ok: fragment index preserves universal neutral ions and mmap cache\n";
		return 0;
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
