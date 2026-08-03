#include "isotopologue.h"
#include "proNovoConfig.h"

#include <cctype>
#include <unordered_map>
#include <utility>

namespace
{

bool hasNegativeAtomCount(const sipros::SourcedComposition &composition)
{
	for (size_t source = 0; source < sipros::IsotopeSourceCount; ++source)
		for (size_t element = 0; element < sipros::ElementCount; ++element)
			if (composition.atoms[source][element] < 0)
				return true;
	return false;
}

} // namespace

IsotopeDistribution::IsotopeDistribution()
{
}

IsotopeDistribution::IsotopeDistribution(vector<double> vItsMass, vector<double> vItsProb)
{
	vMass = vItsMass;
	vProb = vItsProb;
}

double IsotopeDistribution::getMostAbundantMass()
{
	double dMaxProb = 0;
	double dMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		if (dMaxProb < vProb[i])
		{
			dMaxProb = vProb[i];
			dMass = vMass[i];
		}
	}
	return dMass;
}

Isotopologue::Isotopologue()
{
	AtomNumber = 0;
}

bool Isotopologue::setupIsotopologue(
	const map<string, sipros::SourcedComposition> &residueAtomicComposition,
	const vector<IsotopeDistribution> &atomIsotopicDistribution,
	const string &atomNames)
{
	AtomNumber = atomNames.size();
	if (AtomNumber == 0 ||
		atomIsotopicDistribution.size() != AtomNumber)
	{
		cerr << "ERROR: built-in atom names and isotope distributions do not match." << endl;
		return false;
	}
	if (AtomNumber != sipros::ElementCount)
	{
		cerr << "ERROR: built-in chemistry must define the six real CHONPS elements."
			 << endl;
		return false;
	}

	mResidueSourcedComposition = residueAtomicComposition;
	vAtomIsotopicDistribution = atomIsotopicDistribution;
	for (IsotopeDistribution &distribution : vAtomIsotopicDistribution)
	{
		if (distribution.vMass.size() != distribution.vProb.size() ||
			!CheckMass(distribution.vMass, distribution.vProb))
		{
			cerr << "ERROR: invalid built-in isotope distribution." << endl;
			return false;
		}
	}
	vNaturalAtomIsotopicDistribution = vAtomIsotopicDistribution;
	return refreshResidueStateCache({});
}

bool Isotopologue::refreshResidueStateCache(
	const map<string, string> &ptmSites,
	const vector<string> &nTermPtms,
	const vector<string> &cTermPtms)
{
	vResidueIsotopicDistribution.clear();
	residueStateCache.clear();
	bIonTerminalStateCache.clear();
	yIonTerminalStateCache.clear();
	for (const auto &entry : mResidueSourcedComposition)
	{
		IsotopeDistribution distribution;
		if (!computeIsotopicDistribution(entry.second, distribution))
		{
			cerr << "ERROR: cannot calculate the isotopic distribution for "
				 << entry.first << endl;
			return false;
		}
		vResidueIsotopicDistribution[entry.first] = distribution;
		if (entry.first.size() == 1 &&
			std::isalpha(static_cast<unsigned char>(entry.first.front())))
		{
			residueStateCache.emplace(
				entry.first,
				CachedResidueState{entry.second, std::move(distribution)});
		}
	}

	for (const auto &ptm : ptmSites)
	{
		if (ptm.first.size() != 1)
			continue;
		const auto ptmComposition =
			mResidueSourcedComposition.find(ptm.first);
		if (ptmComposition == mResidueSourcedComposition.end())
		{
			cerr << "ERROR: cannot find PTM composition for "
				 << ptm.first << endl;
			return false;
		}
		for (char site : ptm.second)
		{
			const string residueName(1, site);
			const auto residueComposition =
				mResidueSourcedComposition.find(residueName);
			if (residueComposition == mResidueSourcedComposition.end())
			{
				cerr << "ERROR: cannot find PTM site residue "
					 << residueName << endl;
				return false;
			}
			sipros::SourcedComposition combined =
				residueComposition->second + ptmComposition->second;
			if (hasNegativeAtomCount(combined))
			{
				cerr << "ERROR: negative sourced atom count for residue state "
					 << residueName + ptm.first << endl;
				return false;
			}
			IsotopeDistribution distribution;
			if (!computeIsotopicDistribution(combined, distribution))
				return false;
			residueStateCache[residueName + ptm.first] =
				CachedResidueState{combined, std::move(distribution)};
		}
	}

	const auto nTerm = mResidueSourcedComposition.find("Nterm");
	const auto cTerm = mResidueSourcedComposition.find("Cterm");
	if (nTerm == mResidueSourcedComposition.end() ||
		cTerm == mResidueSourcedComposition.end())
	{
		cerr << "ERROR: cannot find peptide terminal compositions" << endl;
		return false;
	}
	const sipros::SourcedComposition yIonBaseComposition =
		nTerm->second + cTerm->second;
	auto cacheTerminalState = [&](unordered_map<string, IsotopeDistribution> &cache,
								  const string &stateKey,
								  const sipros::SourcedComposition &composition,
								  const char *terminus)
		-> bool
	{
		if (hasNegativeAtomCount(composition))
		{
			cerr << "ERROR: negative sourced atom count for " << terminus
				 << " terminal state " << stateKey << endl;
			return false;
		}
		IsotopeDistribution distribution;
		if (!computeIsotopicDistribution(composition, distribution))
			return false;
		cache[stateKey] = std::move(distribution);
		return true;
	};
	if (!cacheTerminalState(bIonTerminalStateCache, "", {}, "b-ion") ||
		!cacheTerminalState(
			yIonTerminalStateCache, "", yIonBaseComposition, "y-ion"))
		return false;
	auto cachePtmStates = [&](const vector<string> &ptms,
							  unordered_map<string, IsotopeDistribution> &cache,
							  const sipros::SourcedComposition &baseComposition,
							  const char *terminus) -> bool
	{
		for (const string &ptm : ptms)
		{
			if (ptm.size() != 1)
			{
				cerr << "ERROR: terminal PTM tokens must be one character" << endl;
				return false;
			}
			const auto ptmComposition = mResidueSourcedComposition.find(ptm);
			if (ptmComposition == mResidueSourcedComposition.end())
			{
				cerr << "ERROR: cannot find terminal PTM composition for "
					 << ptm << endl;
				return false;
			}
			if (!cacheTerminalState(
					cache, ptm, baseComposition + ptmComposition->second,
					terminus))
				return false;
		}
		return true;
	};
	return cachePtmStates(nTermPtms, bIonTerminalStateCache, {}, "b-ion") &&
		cachePtmStates(
			cTermPtms, yIonTerminalStateCache, yIonBaseComposition, "y-ion");
}

bool Isotopologue::getCachedResidueState(
	const string &state,
	sipros::SourcedComposition &composition,
	IsotopeDistribution &distribution) const
{
	const auto cached = residueStateCache.find(state);
	if (cached == residueStateCache.end())
		return false;
	composition = cached->second.composition;
	distribution = cached->second.distribution;
	return true;
}

bool Isotopologue::getSingleResidueMostAbundantMasses(vector<string> &vsResidues, vector<double> &vdMostAbundantMasses, double &dTerminusMassN,
													  double &dTerminusMassC)
{
	vsResidues.clear();
	vdMostAbundantMasses.clear();

	map<string, IsotopeDistribution>::iterator ResidueIter;
	string sCurrentResidue;
	IsotopeDistribution currentDistribution;
	double dCurrentMostAbundantMasses;

	// for single amino acid

	for (ResidueIter = vResidueIsotopicDistribution.begin(); ResidueIter != vResidueIsotopicDistribution.end(); ResidueIter++)
	{
		sCurrentResidue = ResidueIter->first;
		currentDistribution = ResidueIter->second;
		dCurrentMostAbundantMasses = currentDistribution.getMostAbundantMass();
		if (sCurrentResidue.size() == 1)
		{
			vsResidues.push_back(sCurrentResidue);
			vdMostAbundantMasses.push_back(dCurrentMostAbundantMasses);
			//			cout << sCurrentResidue << "  " << dCurrentMostAbundantMasses << endl;
		}
		else if (sCurrentResidue == "NTerm" || sCurrentResidue == "Nterm")
		{
			dTerminusMassN = dCurrentMostAbundantMasses;
		}
		else if (sCurrentResidue == "CTerm" || sCurrentResidue == "Cterm")
		{
			dTerminusMassC = dCurrentMostAbundantMasses;
		}
		else
		{
			cerr << "ERROR: Unknown residue or PTM " << sCurrentResidue << endl;
		}
	}

	unsigned int i;

	// bubble sort the list by mass
	unsigned int n = vsResidues.size();
	unsigned int pass;
	double dCurrentMass;
	for (pass = 1; pass < n; pass++)
	{ // count how many times
		// This next loop becomes shorter and shorter
		for (i = 0; i < n - pass; i++)
		{
			if (vdMostAbundantMasses.at(i) > vdMostAbundantMasses.at(i + 1))
			{
				// exchange
				dCurrentMass = vdMostAbundantMasses.at(i);
				sCurrentResidue = vsResidues.at(i);

				vdMostAbundantMasses.at(i) = vdMostAbundantMasses.at(i + 1);
				vsResidues.at(i) = vsResidues.at(i + 1);

				vdMostAbundantMasses.at(i + 1) = dCurrentMass;
				vsResidues.at(i + 1) = sCurrentResidue;
			}
		}
	}

	return true;
}

bool Isotopologue::computeIsotopicDistribution(string sSequence, IsotopeDistribution &myIsotopeDistribution)
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;
	map<string, IsotopeDistribution>::iterator ResidueIter;

	ResidueIter = vResidueIsotopicDistribution.find("Nterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cerr << "ERROR: can't find the N-terminus" << endl;
		return false;
	}

	ResidueIter = vResidueIsotopicDistribution.find("Cterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = sum(currentDistribution, sumDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// add up all residues's isotopic distribution
	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		ResidueIter = vResidueIsotopicDistribution.find(currentResidue);
		if (ResidueIter != vResidueIsotopicDistribution.end())
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum(currentDistribution, sumDistribution);
		}
		else
		{
			cerr << "ERROR: can't find the residue " << currentResidue << endl;
			return false;
		}
	}

	myIsotopeDistribution = sumDistribution;

	return true;
}

bool Isotopologue::computePeptideIsotopicDistribution(
	string decoratedSequence,
	IsotopeDistribution &myIsotopeDistribution)
{
	if (decoratedSequence.empty() || decoratedSequence.front() != '[')
	{
		cerr << "ERROR: First character in a peptide sequence must be [."
			 << endl;
		return false;
	}
	const size_t closingBracket = decoratedSequence.find(']', 1);
	if (closingBracket == string::npos ||
		decoratedSequence.find('[', 1) != string::npos ||
		decoratedSequence.find(']', closingBracket + 1) != string::npos)
	{
		cerr << "ERROR: Invalid decorated peptide brackets." << endl;
		return false;
	}

	string symbols = decoratedSequence.substr(1, closingBracket - 1);
	symbols.append(decoratedSequence.substr(closingBracket + 1));
	if (symbols.empty())
	{
		cerr << "ERROR: Peptide sequence is empty." << endl;
		return false;
	}
	sipros::SourcedComposition composition;
	return computeSourcedComposition(symbols, composition) &&
		computeIsotopicDistribution(composition, myIsotopeDistribution);
}

bool Isotopologue::parseProductIonSequence(
	const string &sequence,
	ProductIonState &state) const
{
	state = ProductIonState{};
	if (sequence.empty() || sequence.front() != '[')
	{
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}

	const size_t closingBracket = sequence.find(']', 1);
	if (closingBracket == string::npos)
	{
		cerr << "ERROR: Peptide sequence must contain a closing ]." << endl;
		return false;
	}
	if (sequence.find(']', closingBracket + 1) != string::npos)
	{
		cerr << "ERROR: Peptide sequence contains more than one closing ]." << endl;
		return false;
	}

	auto addNamedComposition = [&](char symbol,
								   sipros::SourcedComposition &composition,
								   string &stateKey,
								   const char *kind) -> bool
	{
		const string name(1, symbol);
		const auto found = mResidueSourcedComposition.find(name);
		if (found == mResidueSourcedComposition.end())
		{
			cerr << "ERROR: cannot find " << kind
				 << " in the built-in chemistry: " << name << endl;
			return false;
		}
		composition += found->second;
		stateKey.push_back(symbol);
		return true;
	};
	const auto isResidue = [](char symbol)
	{
		return std::isalpha(static_cast<unsigned char>(symbol)) != 0;
	};

	size_t cursor = 1;
	sipros::SourcedComposition nTermPtmComposition;
	while (cursor < closingBracket && !isResidue(sequence[cursor]))
	{
		if (!state.nTermStateKey.empty())
		{
			cerr << "ERROR: multiple N-terminal PTMs are not supported." << endl;
			return false;
		}
		if (sequence[cursor] == '[' ||
			!addNamedComposition(sequence[cursor], nTermPtmComposition,
								 state.nTermStateKey,
								 "N-terminal PTM"))
			return false;
		++cursor;
	}

	while (cursor < closingBracket)
	{
		if (!isResidue(sequence[cursor]))
		{
			cerr << "ERROR: expected an amino-acid residue in peptide sequence."
				 << endl;
			return false;
		}

		sipros::SourcedComposition residueComposition;
		string residueStateKey;
		if (!addNamedComposition(sequence[cursor], residueComposition,
								 residueStateKey,
								 "residue"))
			return false;
		++cursor;
		while (cursor < closingBracket && !isResidue(sequence[cursor]))
		{
			if (residueStateKey.size() > 1)
			{
				cerr << "ERROR: multiple PTMs on one residue are not supported."
					 << endl;
				return false;
			}
			if (sequence[cursor] == '[' ||
				!addNamedComposition(sequence[cursor], residueComposition,
									 residueStateKey,
									 "PTM"))
				return false;
			++cursor;
		}
		if (hasNegativeAtomCount(residueComposition))
		{
			cerr << "ERROR: negative sourced atom count for residue state "
				 << residueStateKey << endl;
			return false;
		}
		state.residueStateKeys.push_back(std::move(residueStateKey));
	}

	sipros::SourcedComposition cTermPtmComposition;
	for (cursor = closingBracket + 1; cursor < sequence.size(); ++cursor)
	{
		if (!state.cTermStateKey.empty())
		{
			cerr << "ERROR: multiple C-terminal PTMs are not supported." << endl;
			return false;
		}
		if (isResidue(sequence[cursor]) || sequence[cursor] == '[' ||
			sequence[cursor] == ']' ||
			!addNamedComposition(sequence[cursor], cTermPtmComposition,
								 state.cTermStateKey,
								 "C-terminal PTM"))
		{
			cerr << "ERROR: invalid character after the peptide closing ]."
				 << endl;
			return false;
		}
	}

	const int peptideLength = static_cast<int>(state.residueStateKeys.size());
	if (peptideLength < ProNovoConfig::getMinPeptideLength())
	{
		cerr << "ERROR: Peptide sequence is too short " << sequence << endl;
		return false;
	}
	if (bIonTerminalStateCache.find(state.nTermStateKey) ==
			bIonTerminalStateCache.end())
	{
		cerr << "ERROR: unsupported N-terminal PTM state "
			 << state.nTermStateKey << endl;
		return false;
	}
	if (yIonTerminalStateCache.find(state.cTermStateKey) ==
			yIonTerminalStateCache.end())
	{
		cerr << "ERROR: unsupported C-terminal PTM state "
			 << state.cTermStateKey << endl;
		return false;
	}
	for (const string &stateKey : state.residueStateKeys)
	{
		if (residueStateCache.find(stateKey) == residueStateCache.end())
		{
			cerr << "ERROR: unsupported residue/PTM state " << stateKey << endl;
			return false;
		}
	}
	return true;
}

bool Isotopologue::computeProductIon(
	string sequence,
	vector<vector<double>> &yIonMass,
	vector<vector<double>> &yIonProb,
	vector<vector<double>> &bIonMass,
	vector<vector<double>> &bIonProb) const
{
	yIonMass.clear();
	yIonProb.clear();
	bIonMass.clear();
	bIonProb.clear();
	ProductIonState state;
	if (!parseProductIonSequence(sequence, state))
		return false;

	const size_t productIonCount = state.residueStateKeys.size() - 1;
	yIonMass.reserve(productIonCount);
	yIonProb.reserve(productIonCount);
	bIonMass.reserve(productIonCount);
	bIonProb.reserve(productIonCount);

	IsotopeDistribution cumulative =
		bIonTerminalStateCache.at(state.nTermStateKey);
	for (size_t residue = 0; residue < productIonCount; ++residue)
	{
		cumulative = sum(
			residueStateCache.at(state.residueStateKeys[residue]).distribution,
			cumulative);
		bIonMass.push_back(cumulative.vMass);
		bIonProb.push_back(cumulative.vProb);
	}

	cumulative = yIonTerminalStateCache.at(state.cTermStateKey);
	for (size_t offset = 0; offset < productIonCount; ++offset)
	{
		const size_t residue = state.residueStateKeys.size() - 1 - offset;
		cumulative = sum(
			residueStateCache.at(state.residueStateKeys[residue]).distribution,
			cumulative);
		yIonMass.push_back(cumulative.vMass);
		yIonProb.push_back(cumulative.vProb);
	}

	// Keep product envelopes as neutral-composition masses.  Scoring converts
	// each isotope peak to charge z with (neutralMass + z * protonMass) / z.

	return true;
}

bool Isotopologue::computeIsotopicDistribution(
	const sipros::SourcedComposition &composition,
	IsotopeDistribution &myIsotopeDistribution)
{
	if (vAtomIsotopicDistribution.size() != AtomNumber ||
		vNaturalAtomIsotopicDistribution.size() != AtomNumber)
		return false;

	IsotopeDistribution result({0.0}, {1.0});
	for (size_t element = 0; element < AtomNumber; ++element)
	{
		const int biosyntheticCount =
			composition[sipros::IsotopeSource::Biosynthetic][element];
		const int naturalCount =
			composition[sipros::IsotopeSource::ReagentNatural][element] +
			composition[sipros::IsotopeSource::DigestionSolvent][element];
		const IsotopeDistribution &active =
			vAtomIsotopicDistribution[element];
		const IsotopeDistribution &natural =
			vNaturalAtomIsotopicDistribution[element];

		// At natural abundance, grouping identical channels produces the same
		// convolution (and pruning order) as a conventional flat formula.
		if (active.vMass == natural.vMass && active.vProb == natural.vProb)
		{
			const int totalCount = biosyntheticCount + naturalCount;
			if (totalCount != 0)
				result = sum(multiply(natural, totalCount), result);
			continue;
		}
		if (biosyntheticCount != 0)
			result = sum(multiply(active, biosyntheticCount), result);
		if (naturalCount != 0)
			result = sum(multiply(natural, naturalCount), result);
	}
	myIsotopeDistribution = result;
	return true;
}

bool Isotopologue::computeSourcedComposition(
	string sSequence,
	sipros::SourcedComposition &composition)
{
	composition = {};
	auto addNamedComposition = [&](const string &name) -> bool
	{
		const auto residue = mResidueSourcedComposition.find(name);
		if (residue == mResidueSourcedComposition.end())
			return false;
		composition += residue->second;
		return true;
	};

	if (!addNamedComposition("Nterm"))
	{
		cerr << "ERROR: can't find the atomic composition for the N-terminus"
			 << endl;
		return false;
	}
	if (!addNamedComposition("Cterm"))
	{
		cerr << "ERROR: can't find the atomic composition for the C-terminus"
			 << endl;
		return false;
	}

	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		if (!addNamedComposition(currentResidue))
		{
			cerr << "ERROR: can't find the atomic composition for residue/PTM: " << currentResidue << endl;
			return false;
		}
	}
	return true;
}

IsotopeDistribution Isotopologue::sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1) const
{
	double ProbabilityCutoff_local = 0.000001;

	IsotopeDistribution sumDistribution;
	int iSizeDistribution0 = distribution0.vMass.size();
	int iSizeDistribution1 = distribution1.vMass.size();
	size_t newSize = iSizeDistribution0 + iSizeDistribution1 - 1;
	sumDistribution.vMass.reserve(newSize);
	sumDistribution.vProb.reserve(newSize);
	for (int k = 0; k < newSize; k++)
	{
		double sumweight = 0, summass = 0;
		int start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		int end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1;	   // min(g_n - 1, k)
		int iCount = 0;
		double dSum = 0;
		for (int i = start; i <= end; i++)
		{
			double weight = distribution0.vProb[i] * distribution1.vProb[k - i];
			double mass = distribution0.vMass[i] + distribution1.vMass[k - i];
			sumweight += weight;
			summass += weight * mass;
			iCount += 1;
			dSum += mass;
		}
		if (sumweight == 0)
		{
			sumDistribution.vMass.push_back(dSum / ((double)iCount));
		}
		else
		{
			sumDistribution.vMass.push_back(summass / sumweight);
		}
		sumDistribution.vProb.push_back(sumweight);
	}

	// prune small probabilities
	vector<double>::iterator iteProb = sumDistribution.vProb.begin();
	vector<double>::iterator iteMass = sumDistribution.vMass.begin();
	while (iteProb != sumDistribution.vProb.end())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb++;
		iteMass++;
	}
	if (iteProb != sumDistribution.vProb.begin())
	{
		sumDistribution.vProb.erase(sumDistribution.vProb.begin(), iteProb);
		sumDistribution.vMass.erase(sumDistribution.vMass.begin(), iteMass);
	}

	iteProb = sumDistribution.vProb.end() - 1;
	iteMass = sumDistribution.vMass.end() - 1;
	while (iteProb != sumDistribution.vProb.begin())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb--;
		iteMass--;
	}
	if (iteProb != sumDistribution.vProb.end() - 1)
	{
		sumDistribution.vProb.erase(iteProb + 1, sumDistribution.vProb.end());
		sumDistribution.vMass.erase(iteMass + 1, sumDistribution.vMass.end());
	}

	// normalize the probability space to 1
	int iSizeSumDistribution;
	int i;
	double sumProb = 0;
	iSizeSumDistribution = sumDistribution.vMass.size();
	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumProb += sumDistribution.vProb[i];
	}

	if (sumProb <= 0)
	{
		cerr << "Error: Sum of distribution is zero" << endl;
		exit(1);
		return sumDistribution;
	}

	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumDistribution.vProb[i] = sumDistribution.vProb[i] / sumProb;
	}

	return sumDistribution;
}

IsotopeDistribution Isotopologue::multiply(const IsotopeDistribution &distribution0, int count)
{
	if (count == 1)
		return distribution0;
	IsotopeDistribution productDistribution;
	productDistribution.vMass.push_back(0.0);
	productDistribution.vProb.push_back(1.0);

	if (count < 0)
	{
		IsotopeDistribution negativeDistribution = distribution0;
		for (unsigned int n = 0; n < negativeDistribution.vMass.size(); n++)
		{
			negativeDistribution.vMass[n] = -negativeDistribution.vMass[n];
		}

		for (int i = 0; i < abs(count); ++i)
		{
			productDistribution = sum(productDistribution, negativeDistribution);
		}
	}
	else
	{
		struct DistributionPowerCache
		{
			vector<double> masses;
			vector<double> probabilities;
			vector<IsotopeDistribution> powers;
		};
		// Build each atom-count convolution once per thread. Extending the
		// previous count with one atom retains the existing left-to-right
		// pruning and normalization order exactly.
		thread_local std::unordered_map<
			const IsotopeDistribution *, DistributionPowerCache> caches;
		DistributionPowerCache &cache = caches[&distribution0];
		if (cache.masses != distribution0.vMass ||
			cache.probabilities != distribution0.vProb)
		{
			cache.masses = distribution0.vMass;
			cache.probabilities = distribution0.vProb;
			cache.powers.clear();
			cache.powers.push_back(productDistribution);
			cache.powers.push_back(distribution0);
		}
		while (cache.powers.size() <= static_cast<size_t>(count))
			cache.powers.push_back(sum(cache.powers.back(), distribution0));
		productDistribution = cache.powers[static_cast<size_t>(count)];
	}
	return productDistribution;
}

bool Isotopologue::CheckMass(vector<double> &vdMass, vector<double> &vdNaturalCompositionTemp)
{
	double dDiff = 0;

	int iMassCount = vdMass.size();
	if (iMassCount < 2)
	{
		return true;
	}
	int pass, i;
	double dCurrentMass;
	double dCurrentComposition;
	// count how many times
	// This next loop becomes shorter and shorter
	for (pass = 1; pass < iMassCount; pass++)
	{
		for (i = 0; i < iMassCount - pass; i++)
		{
			if (vdMass.at(i) > vdMass.at(i + 1))
			{
				// exchange
				dCurrentMass = vdMass.at(i);
				dCurrentComposition = vdNaturalCompositionTemp.at(i);

				vdMass.at(i) = vdMass.at(i + 1);
				vdNaturalCompositionTemp.at(i) = vdNaturalCompositionTemp.at(i + 1);

				vdMass.at(i + 1) = dCurrentMass;
				vdNaturalCompositionTemp.at(i + 1) = dCurrentComposition;
			}
		}
	}

	for (i = 1; i < (int)vdMass.size(); ++i)
	{
		dDiff = round(vdMass.at(i) - vdMass.at(i - 1));
		if (dDiff != 1.0)
		{
			return false;
		}
	}
	return true;
}
