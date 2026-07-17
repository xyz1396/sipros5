#pragma once
#include "isotopologue.h"

class PeptideIsotopeCalculator
{
private:
	static double lightMass(size_t element);
	static double baseMass(const sipros::SourcedComposition &composition);
public:
	struct PrecursorEstimate
	{
		double mass = 0.0;
		double neutronMass = 0.0;
		int nominalShift = 0;
	};

	// Biosynthetic counts are SIP-labelable. Reagent and digestion-solvent
	// counts retain their natural isotope distributions.
	sipros::SourcedComposition pepComposition;

    static bool changeAtomProbability(std::vector<double> &probs, char atom, const double pct);
    void calPepAtomCounts(const string &pepSeq);
    double calPrecursorBaseMass(const string &pepSeq);
	// Center precursor windows on the modal bin of the full source-aware
	// isotope convolution; neutronMass retains its composition-weighted value.
	PrecursorEstimate calPrecursorEstimate(const string &pepSeq);
    double calPrecursorMass(const string &pepSeq);
};
