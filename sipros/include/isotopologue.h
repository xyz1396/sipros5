#ifndef ISOTOPOLOGUE_H
#define ISOTOPOLOGUE_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include "proNovoConfig.h"
using namespace std;

/**/
class IsotopeDistribution
{
public:
	IsotopeDistribution();
	IsotopeDistribution(vector<double> vItsMass, vector<double> vItsComposition);

	vector<double> vMass;
	vector<double> vProb;

	double getMostAbundantMass();
};

class Isotopologue
{
public:
	Isotopologue();

	// Initialize chemistry from the built-in, strongly typed profile.
	bool setupIsotopologue(
		const map<string, sipros::SourcedComposition> &residueAtomicComposition,
		const vector<IsotopeDistribution> &atomIsotopicDistribution,
		const string &atomNames);
	bool refreshResidueStateCache(
		const map<string, string> &ptmSites,
		const vector<string> &nTermPtms = {},
		const vector<string> &cTermPtms = {});
	bool getCachedResidueState(
		const string &state,
		sipros::SourcedComposition &composition,
		IsotopeDistribution &distribution) const;

	// get the MostAbundant masses of  residues
	bool getSingleResidueMostAbundantMasses(vector<string> &vsResidues, vector<double> &vdMostAbundantMasses, double &dTerminusMassN,
											double &dTerminusMassC);

	// variables for this isotopologue
	map<string, sipros::SourcedComposition> mResidueSourcedComposition;
	// Active distributions apply only to biosynthetic atoms. Natural-source
	// atoms always use this immutable natural-distribution snapshot.
	vector<IsotopeDistribution> vAtomIsotopicDistribution;
	vector<IsotopeDistribution> vNaturalAtomIsotopicDistribution;
	map<string, IsotopeDistribution> vResidueIsotopicDistribution;

	// emass functions for IsotopeDistribution's arithmetic
	IsotopeDistribution sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1);

	// compute isotoptic distributions for all product ions of a sequence
	// this isotopologue class is modified by only adding this function
	// The first dimension of vvdYionMass and vvdYionProb is from y1, y2, ...
	// The first dimension of vvdBionMass and vvdBionProb is from b1, b2, ...
	// the mass is calculated assuming cleavage of the peptide bond
	bool computeProductIon(string sSequence, vector<vector<double>> &vvdYionMass, vector<vector<double>> &vvdYionProb,
						   vector<vector<double>> &vvdBionMass, vector<vector<double>> &vvdBionProb);

	// compute isotoptic distribution for an amino acid sequence
	bool computeIsotopicDistribution(string sSequence, IsotopeDistribution &myIsotopeDistribution);

	// Compute the full precursor distribution for Sipros' decorated peptide
	// syntax, including terminal PTMs outside the square brackets.
	bool computePeptideIsotopicDistribution(
		string decoratedSequence,
		IsotopeDistribution &myIsotopeDistribution);

	// Compute an isotope distribution while preserving atom provenance.
	bool computeIsotopicDistribution(const sipros::SourcedComposition &composition,
									 IsotopeDistribution &myIsotopeDistribution);

	bool computeSourcedComposition(string sSequence,
								  sipros::SourcedComposition &composition);

private:
	struct CachedResidueState
	{
		sipros::SourcedComposition composition;
		IsotopeDistribution distribution;
	};
	struct ProductIonState
	{
		string nTermStateKey;
		string cTermStateKey;
		vector<string> residueStateKeys;
	};
	unordered_map<string, CachedResidueState> residueStateCache;
	unordered_map<string, IsotopeDistribution> bIonTerminalStateCache;
	unordered_map<string, IsotopeDistribution> yIonTerminalStateCache;

	bool parseProductIonSequence(
		const string &sequence,
		ProductIonState &state);

	IsotopeDistribution multiply(const IsotopeDistribution &distribution0, int count);

	// Number of CHONPS elements in the active distributions.
	unsigned int AtomNumber;

	// Sipros Ensemble
	// emass needs the mass to be one nucleus difference
	bool CheckMass(vector<double> &vdMass, vector<double> &vdNaturalCompositionTemp);
};

#endif // ISOTOPOLOGUE_H
