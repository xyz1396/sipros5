#include "isotopologue.h"
#include "proNovoConfig.h"

#include <cctype>
#include <utility>

IsotopeDistribution::IsotopeDistribution()
{
}

IsotopeDistribution::IsotopeDistribution(vector<double> vItsMass, vector<double> vItsProb)
{
	vMass = vItsMass;
	vProb = vItsProb;
}

IsotopeDistribution::~IsotopeDistribution()
{
	// destructor
}

void IsotopeDistribution::print()
{
	cout << "Mass " << '\t' << "Inten" << endl;
	for (unsigned int i = 0; i < vMass.size(); i++)
	{
		cout << setprecision(8) << vMass[i] << '\t' << vProb[i] << endl;
	}
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

double IsotopeDistribution::getAverageMass()
{
	double dSumProb = 0;
	double dSumMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		dSumProb = dSumProb + vProb[i];
		dSumMass = dSumMass + vProb[i] * vMass[i];
	}

	if (dSumProb <= 0)
		return 1.0;

	return (dSumMass / dSumProb);
}
void IsotopeDistribution::filterProbCutoff(double dProbCutoff)
{
	vector<double> vMassCopy = vMass;
	vector<double> vProbCopy = vProb;
	vMass.clear();
	vProb.clear();
	for (unsigned int i = 0; i < vProbCopy.size(); ++i)
	{
		if (vProbCopy[i] >= dProbCutoff)
		{
			vMass.push_back(vMassCopy[i]);
			vProb.push_back(vProbCopy[i]);
		}
	}
}
double IsotopeDistribution::getLowestMass()
{
	return *min_element(vMass.begin(), vMass.end());
}

Isotopologue::Isotopologue() : MassPrecision(0.01), ProbabilityCutoff(0.000000001)
{
	AtomNumber = 0;
}

Isotopologue::~Isotopologue()
{
	// destructor
}

bool Isotopologue::setupIsotopologue(
	const map<string, sipros::SourcedComposition> &residueAtomicComposition,
	const vector<IsotopeDistribution> &atomIsotopicDistribution,
	const string &atomNames)
{
	AtomName = atomNames;
	AtomNumber = AtomName.size();
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
	vResidueIsotopicDistribution.clear();
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

	// calculate Isotopic distribution for all residues
	map<string, sipros::SourcedComposition>::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for (ResidueIter = mResidueSourcedComposition.begin(); ResidueIter != mResidueSourcedComposition.end(); ResidueIter++)
	{
		if (!computeIsotopicDistribution(ResidueIter->second, tempIsotopeDistribution))
		{
			cerr << "ERROR: cannot calculate the isotopic distribution for residue " << ResidueIter->first << endl;
			return false;
		}

		vResidueIsotopicDistribution[ResidueIter->first] = tempIsotopeDistribution;
	}

	return true;
}

double Isotopologue::computeMostAbundantMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getMostAbundantMass();
}

double Isotopologue::computeAverageMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getAverageMass();
}

double Isotopologue::computeMonoisotopicMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getLowestMass();
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

bool Isotopologue::computeProductIon(string sSequence, vector<vector<double>> &vvdYionMass, vector<vector<double>> &vvdYionProb,
									 vector<vector<double>> &vvdBionMass, vector<vector<double>> &vvdBionProb)
{
	vvdYionMass.clear();
	vvdYionProb.clear();
	vvdBionMass.clear();
	vvdBionProb.clear();

	if (sSequence.empty() || sSequence.front() != '[')
	{
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}

	const size_t closingBracket = sSequence.find(']', 1);
	if (closingBracket == string::npos)
	{
		cerr << "ERROR: Peptide sequence must contain a closing ]." << endl;
		return false;
	}
	if (sSequence.find(']', closingBracket + 1) != string::npos)
	{
		cerr << "ERROR: Peptide sequence contains more than one closing ]." << endl;
		return false;
	}

	auto addNamedComposition = [&](char symbol,
								   sipros::SourcedComposition &composition,
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
		return true;
	};
	auto loadTerminalComposition = [&](const char *name,
									   sipros::SourcedComposition &composition) -> bool
	{
		const auto found = mResidueSourcedComposition.find(name);
		if (found == mResidueSourcedComposition.end())
		{
			cerr << "ERROR: can't find the " << name << " composition" << endl;
			return false;
		}
		composition = found->second;
		return true;
	};
	const auto isResidue = [](char symbol)
	{
		return std::isalpha(static_cast<unsigned char>(symbol)) != 0;
	};

	sipros::SourcedComposition nTermComposition;
	sipros::SourcedComposition cTermComposition;
	if (!loadTerminalComposition("Nterm", nTermComposition) ||
		!loadTerminalComposition("Cterm", cTermComposition))
		return false;
	const sipros::SourcedComposition baseNTermComposition =
		nTermComposition;

	// Decorations between '[' and the first amino acid belong to the
	// N-terminus. Decorations after ']' belong to the C-terminus. Each
	// decoration is a registered one-character PTM token.
	size_t cursor = 1;
	while (cursor < closingBracket && !isResidue(sSequence[cursor]))
	{
		if (sSequence[cursor] == '[' ||
			!addNamedComposition(sSequence[cursor], nTermComposition,
								 "N-terminal PTM"))
			return false;
		++cursor;
	}

	vector<sipros::SourcedComposition> residueCompositions;
	while (cursor < closingBracket)
	{
		if (!isResidue(sSequence[cursor]))
		{
			cerr << "ERROR: expected an amino-acid residue in peptide sequence."
				 << endl;
			return false;
		}

		sipros::SourcedComposition residueComposition;
		if (!addNamedComposition(sSequence[cursor], residueComposition,
								 "residue"))
			return false;
		++cursor;
		while (cursor < closingBracket && !isResidue(sSequence[cursor]))
		{
			if (sSequence[cursor] == '[' ||
				!addNamedComposition(sSequence[cursor], residueComposition,
									 "PTM"))
				return false;
			++cursor;
		}
		residueCompositions.push_back(residueComposition);
	}

	for (cursor = closingBracket + 1; cursor < sSequence.size(); ++cursor)
	{
		if (isResidue(sSequence[cursor]) || sSequence[cursor] == '[' ||
			sSequence[cursor] == ']' ||
			!addNamedComposition(sSequence[cursor], cTermComposition,
								 "C-terminal PTM"))
		{
			cerr << "ERROR: invalid character after the peptide closing ]."
				 << endl;
			return false;
		}
	}

	const int peptideLength = static_cast<int>(residueCompositions.size());
	if (peptideLength < ProNovoConfig::getMinPeptideLength())
	{
		cerr << "ERROR: Peptide sequence is too short " << sSequence << endl;
		return false;
	}

	const size_t productIonCount = residueCompositions.size() - 1;
	vvdYionMass.reserve(productIonCount);
	vvdYionProb.reserve(productIonCount);
	vvdBionMass.reserve(productIonCount);
	vvdBionProb.reserve(productIonCount);

	// PTMs may contain negative atom counts. Aggregate each fragment's net
	// sourced composition before convolution so those atoms cancel against
	// the residue (or fixed modification) that they replace.
	// Product-ion masses are stored as neutral masses; callers add one proton
	// when converting them to charge-one m/z. The base N-terminal hydrogen is
	// therefore transferred to the y ion together with the C-terminal OH.
	// N-terminal PTM deltas remain on b ions, and C-terminal PTM deltas remain
	// on y ions.
	sipros::SourcedComposition bComposition =
		nTermComposition - baseNTermComposition;
	for (size_t residue = 0; residue < productIonCount; ++residue)
	{
		bComposition += residueCompositions[residue];
		IsotopeDistribution distribution;
		if (!computeIsotopicDistribution(bComposition, distribution))
			return false;
		vvdBionMass.push_back(std::move(distribution.vMass));
		vvdBionProb.push_back(std::move(distribution.vProb));
	}

	sipros::SourcedComposition yComposition =
		baseNTermComposition + cTermComposition;
	for (size_t offset = 0; offset < productIonCount; ++offset)
	{
		const size_t residue = residueCompositions.size() - 1 - offset;
		yComposition += residueCompositions[residue];
		IsotopeDistribution distribution;
		if (!computeIsotopicDistribution(yComposition, distribution))
			return false;
		vvdYionMass.push_back(std::move(distribution.vMass));
		vvdYionProb.push_back(std::move(distribution.vProb));
	}

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

IsotopeDistribution Isotopologue::sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1)
{
	double ProbabilityCutoff_local = 0.000001;

	IsotopeDistribution sumDistribution;
	double currentMass;
	double currentProb;
	int iSizeDistribution0 = distribution0.vMass.size();
	int iSizeDistribution1 = distribution1.vMass.size();
	int iCount = 0;
	double dSum = 0;
	size_t newSize = iSizeDistribution0 + iSizeDistribution1 - 1;
	sumDistribution.vMass.reserve(newSize);
	sumDistribution.vProb.reserve(newSize);
#pragma omp simd
	for (int k = 0; k < newSize; k++)
	{
		double sumweight = 0, summass = 0;
		int start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		int end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1;	   // min(g_n - 1, k)
		iCount = 0;
		dSum = 0;
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
			currentMass = dSum / ((double)iCount);
		}
		else
		{
			currentMass = summass / sumweight;
		}
		currentProb = sumweight;
		sumDistribution.vMass.push_back(currentMass);
		sumDistribution.vProb.push_back(currentProb);
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
#pragma omp simd
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
		for (int i = 0; i < count; ++i)
		{
			productDistribution = sum(productDistribution, distribution0);
		}
	}
	return productDistribution;
}

void Isotopologue::shiftMass(IsotopeDistribution &distribution0, double dMass)
{
	for (unsigned int i = 0; i < distribution0.vMass.size(); ++i)
	{
		distribution0.vMass[i] = distribution0.vMass[i] + dMass;
	}
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
