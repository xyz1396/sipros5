#pragma once
#include "isotopologue.h"
#include <array>

class averagine
{
private:
	static double lightMass(size_t element);
	static double baseMass(const sipros::SourcedComposition &composition);
public:
	 const string averagineResidue = "a";
	 const string SIPatoms = "CHONS";
	 sipros::SourcedComposition averagineComposition;
	 vector<sipros::SourcedComposition> averaginePepCompositions;
	 // Biosynthetic counts are the SIP-labelable counts. The complete sourced
	 // composition is retained so natural reagent/solvent background is not lost.
	 sipros::SourcedComposition pepComposition;
	 // C,H,O,N,P,S Atom count of BYions
	 std::vector<sipros::SourcedComposition> BionsCompositions;
	 std::vector<sipros::SourcedComposition> YionsCompositions;
    std::vector<double> BionsBaseMasses;
    std::vector<double> YionsBaseMasses;
    // Atom count difference bettween averagine and peptide
	 IsotopeDistribution averagineSIPdistribution;
    vector<IsotopeDistribution> averaginePepSIPdistributions;
	 int minPepLen = 0, maxPepLen = 0, pepLenRange = 0;
    averagine(const int minPepLen, const int maxPepLen);
    averagine();
    ~averagine();
    static bool changeAtomProbability(std::vector<double> &probs, char atom, const double pct);
    void changeAtomSIPabundance(const char SIPatom, const double pct);
    // Effective mass per nominal-neutron shift for the estimated precursor
    // isotopologue. O18/S34 contribute their full +2 isotope delta before
    // this value is normalized; use calPrecursorMass when an exact peak mass
    // rather than a nominal mass-window spacing is required.
    double calNetronMass(const string &pepSeq);
    void calAveraginePepAtomCounts();
	 sipros::SourcedComposition *getAveraginePepComposition(const int pepLen);
    void calAveraginePepSIPdistributions();
    IsotopeDistribution *getAveraginePepSIPdistribution(const int pepLen);
    void calPepAtomCounts(const string &pepSeq);
    void calBYionsAtomCounts(const string &pepSeq);
    // for peptide base mass without isotope
    double calPrecursorBaseMass(const string &pepSeq);
    void calBYionBaseMasses(const string &pepSeq);
    double calPrecursorMass(const string &pepSeq);
    void calDiffAtomCounts(const string &pepSeq);
    void calPrecursorIsotopeDistribution(const string &pepSeq, IsotopeDistribution &tempSIPdistribution);
};
