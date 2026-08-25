#include "MVH.h"
#include "fragmentindex.h"
#include "peptide.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#elif defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif



void check(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

std::string canonicalNakedPeptide(std::string_view decorated)
{
	const size_t open = decorated.find('[');
	const size_t close = decorated.rfind(']');
	if (open != std::string_view::npos && close > open)
	{
		decorated = decorated.substr(open + 1, close - open - 1);
	}
	std::string result;
	for (char symbol : decorated)
	{
		if (std::isalpha(static_cast<unsigned char>(symbol)) == 0)
		{
			continue;
		}
		char residue = static_cast<char>(std::toupper(
			static_cast<unsigned char>(symbol)));
		if (residue == 'I' || residue == 'J')
		{
			residue = 'L';
		}
		result.push_back(residue);
	}
	return result;
}


int main()
{
	try
	{
		const unsigned long processId =
#if defined(_WIN32)
			static_cast<unsigned long>(_getpid());
#elif defined(__unix__) || defined(__APPLE__)
			static_cast<unsigned long>(getpid());
#else
			0;
#endif
		const std::filesystem::path fasta = std::filesystem::temp_directory_path() /
			("sipros_fragment_index_test_" + std::to_string(processId) + ".fasta");
		const std::filesystem::path cache = std::filesystem::temp_directory_path() /
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

		std::string error;
		{
			sipros::FragmentIndex built;
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
		check(std::filesystem::file_size(cache) == built.stats().cacheBytes,
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
		}

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

		const std::filesystem::path guardDirectory = std::filesystem::temp_directory_path() /
			("sipros_fragment_index_guard_test_" +
			 std::to_string(processId));
		std::filesystem::create_directories(guardDirectory);
		const std::filesystem::path targetFasta = guardDirectory / "target.fasta";
		const std::filesystem::path decoyFasta = guardDirectory / "decoy.fasta";
		const std::filesystem::path targetCache = guardDirectory / "target.sfi";
		const std::filesystem::path legacyDecoyCache = guardDirectory / "legacy.sfi";
		const std::filesystem::path decoyCache = guardDirectory / "decoy.sfi";
		{
			std::ofstream out(targetFasta);
			out << ">target_collision\n"
				<< "MPEPTIDERKACDEFGHIKLMNPQRSTVWYRK\n";
		}
		{
			std::ofstream out(decoyFasta);
			out << ">decoy_collision\n"
				<< "MPEPTIDERKACDEFGHIKLMNPQRSTVWYRK\n"
				<< ">decoy_unique\n"
				<< "MYYYYYYYYKGGGGGGGGR\n";
		}

		ProNovoConfig::setFASTAfilename(targetFasta.string());
		sipros::FragmentIndex targetGuardIndex;
		check(targetGuardIndex.loadOrBuild(
			targetCache.string(), true, error),
			"could not build collision-guard target.sfi: " + error);
		std::unordered_set<std::string> targetIdentities;
		for (uint32_t peptideId = 0;
			 peptideId < targetGuardIndex.peptideCount(); ++peptideId)
		{
			targetIdentities.insert(canonicalNakedPeptide(
				targetGuardIndex.peptideSequence(peptideId)));
		}

		ProNovoConfig::setFASTAfilename(decoyFasta.string());
		uint64_t legacyDecoyPeptideCount = 0;
		std::vector<std::string> expectedGuardedPeptides;
		{
			sipros::FragmentIndex legacyDecoyIndex;
			check(legacyDecoyIndex.loadOrBuild(
				legacyDecoyCache.string(), true, error),
				"could not build legacy unguarded decoy cache: " + error);
			legacyDecoyPeptideCount = legacyDecoyIndex.peptideCount();
			for (uint32_t peptideId = 0;
				 peptideId < legacyDecoyIndex.peptideCount(); ++peptideId)
			{
				const std::string peptide(
					legacyDecoyIndex.peptideSequence(peptideId));
				if (targetIdentities.find(canonicalNakedPeptide(peptide)) ==
					targetIdentities.end())
				{
					expectedGuardedPeptides.push_back(peptide);
				}
			}
		}
		std::sort(expectedGuardedPeptides.begin(),
			expectedGuardedPeptides.end());
		std::filesystem::rename(legacyDecoyCache, decoyCache);

		sipros::FragmentIndex guardedDecoyIndex;
		check(guardedDecoyIndex.loadOrBuild(
			decoyCache.string(), false, error),
			"could not build default guarded decoy.sfi: " + error);
		check(!guardedDecoyIndex.stats().loadedFromCache,
			"legacy unguarded decoy.sfi fingerprint was accepted");
		check(guardedDecoyIndex.stats().collisionGuardTargetPeptides ==
			targetIdentities.size(),
			"guarded decoy.sfi reports the wrong target identity count");
		check(guardedDecoyIndex.stats().collisionExcludedPeptides > 0,
			"default decoy.sfi guard did not exclude a target collision");
		check(legacyDecoyPeptideCount - guardedDecoyIndex.peptideCount() ==
			guardedDecoyIndex.stats().collisionExcludedPeptides,
			"guarded decoy.sfi peptide-count reduction does not match its audit count");
		std::vector<std::string> actualGuardedPeptides;
		actualGuardedPeptides.reserve(
			static_cast<size_t>(guardedDecoyIndex.peptideCount()));
		for (uint32_t peptideId = 0;
			 peptideId < guardedDecoyIndex.peptideCount(); ++peptideId)
		{
			const std::string peptide(
				guardedDecoyIndex.peptideSequence(peptideId));
			check(targetIdentities.find(canonicalNakedPeptide(peptide)) ==
				targetIdentities.end(),
				"guarded decoy.sfi still contains a target-equivalent peptide");
			actualGuardedPeptides.push_back(peptide);
		}
		std::sort(actualGuardedPeptides.begin(), actualGuardedPeptides.end());
		check(actualGuardedPeptides == expectedGuardedPeptides,
			"guarded decoy.sfi did not preserve the exact non-colliding peptide records");

		sipros::FragmentIndex warmGuardedDecoyIndex;
		check(warmGuardedDecoyIndex.loadOrBuild(
			decoyCache.string(), false, error),
			"could not reload guarded decoy.sfi: " + error);
		check(warmGuardedDecoyIndex.stats().loadedFromCache,
			"warm guarded decoy.sfi unexpectedly rebuilt");
		check(warmGuardedDecoyIndex.stats().collisionGuardSeconds == 0.0,
			"warm guarded decoy.sfi rebuilt the canonical target set");
		check(warmGuardedDecoyIndex.stats().collisionExcludedPeptides ==
			guardedDecoyIndex.stats().collisionExcludedPeptides,
			"warm guarded decoy.sfi lost its collision audit count");

		std::error_code ignored;
		std::filesystem::remove(cache, ignored);
		std::filesystem::remove(fasta, ignored);
		std::filesystem::remove_all(guardDirectory, ignored);
		std::cout << "ok: fragment index preserves universal neutral ions, mmap cache, and default decoy collision guard\n";
		return 0;
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
