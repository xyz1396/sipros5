#include "PeptideIsotopeCalculator.h"
#include "proNovoConfig.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

double PeptideIsotopeCalculator::lightMass(size_t element)
{
    const auto &atoms = ProNovoConfig::configIsotopologue
                            .vNaturalAtomIsotopicDistribution;
    if (element >= atoms.size() || atoms[element].vMass.empty())
        return 0.0;
    return atoms[element].vMass.front();
}

double PeptideIsotopeCalculator::baseMass(
    const sipros::SourcedComposition &composition)
{
    const sipros::AtomCounts total = composition.total();
    double mass = 0.0;
    for (size_t element = 0; element < total.size(); ++element)
        mass += static_cast<double>(total[element]) * lightMass(element);
    return mass;
}

bool PeptideIsotopeCalculator::changeAtomProbability(
    std::vector<double> &probabilities,
    char atom,
    const double fraction)
{
    const double targetFraction =
        std::max(0.0, std::min(1.0, fraction));
    if (targetFraction != fraction)
    {
        cerr << "Warning: SIP abundance " << fraction
             << " is outside [0,1], clamped to "
             << targetFraction << "." << endl;
    }

    const size_t targetIsotopeIndex =
        (atom == 'O' || atom == 'S') ? 2u : 1u;
    if (probabilities.empty() ||
        targetIsotopeIndex >= probabilities.size())
    {
        cout << atom
             << " target isotope is unavailable in atom distribution!"
             << endl;
        return false;
    }

    double nonTargetTotal = 0.0;
    for (size_t isotope = 0; isotope < probabilities.size(); ++isotope)
        if (isotope != targetIsotopeIndex)
            nonTargetTotal += probabilities[isotope];

    if (!(nonTargetTotal > 0.0))
    {
        std::fill(probabilities.begin(), probabilities.end(), 0.0);
        probabilities[0] = 1.0 - targetFraction;
    }
    else
    {
        const double scale =
            (1.0 - targetFraction) / nonTargetTotal;
        for (size_t isotope = 0; isotope < probabilities.size(); ++isotope)
            if (isotope != targetIsotopeIndex)
                probabilities[isotope] *= scale;
    }
    probabilities[targetIsotopeIndex] = targetFraction;
    return true;
}

void PeptideIsotopeCalculator::calPepAtomCounts(
    const string &peptideSequence)
{
    pepComposition = {};
    if (!ProNovoConfig::configIsotopologue.computeSourcedComposition(
            peptideSequence, pepComposition))
        return;
}

double PeptideIsotopeCalculator::calPrecursorBaseMass(
    const string &peptideSequence)
{
    calPepAtomCounts(peptideSequence);
    return baseMass(pepComposition);
}

PeptideIsotopeCalculator::PrecursorEstimate
PeptideIsotopeCalculator::calPrecursorEstimate(
    const string &peptideSequence)
{
    calPepAtomCounts(peptideSequence);
    const auto &active = ProNovoConfig::configIsotopologue
                             .vAtomIsotopicDistribution;
    const auto &natural = ProNovoConfig::configIsotopologue
                              .vNaturalAtomIsotopicDistribution;
    if (active.size() < sipros::ElementCount ||
        natural.size() < sipros::ElementCount)
        throw std::runtime_error(
            "Cannot estimate precursor mass: incomplete isotope distributions.");

    double expectedMassShift = 0.0;
    double expectedNominalShift = 0.0;
    for (size_t source = 0;
         source < sipros::IsotopeSourceCount;
         ++source)
    {
        const bool biosynthetic =
            source == static_cast<size_t>(
                          sipros::IsotopeSource::Biosynthetic);
        for (size_t element = 0;
             element < sipros::ElementCount;
             ++element)
        {
            const int atomCount = pepComposition.atoms[source][element];
            if (atomCount <= 0)
                continue;
            const IsotopeDistribution &distribution =
                biosynthetic ? active[element] : natural[element];
            const size_t isotopeCount = std::min(
                distribution.vMass.size(), distribution.vProb.size());
            if (isotopeCount == 0)
                throw std::runtime_error(
                    "Cannot estimate precursor mass: empty isotope distribution.");
            for (size_t isotope = 1;
                 isotope < isotopeCount;
                 ++isotope)
            {
                const double expectedCount =
                    static_cast<double>(atomCount) *
                    distribution.vProb[isotope];
                const double exactDelta =
                    distribution.vMass[isotope] -
                    distribution.vMass.front();
                expectedMassShift += expectedCount * exactDelta;
                // The built-in isotope-table index is its nominal neutron
                // shift: O18/S34 are +2 and S36 is +4.
                expectedNominalShift +=
                    expectedCount * static_cast<double>(isotope);
            }
        }
    }
    if (!(expectedNominalShift > 0.0) ||
        !std::isfinite(expectedMassShift) ||
        !std::isfinite(expectedNominalShift))
        throw std::runtime_error(
            "Cannot estimate precursor mass: invalid expected isotope shift.");

    PrecursorEstimate estimate;
    estimate.neutronMass =
        expectedMassShift / expectedNominalShift;
    if (!(estimate.neutronMass > 0.0) ||
        !std::isfinite(estimate.neutronMass))
        throw std::runtime_error(
            "Cannot estimate precursor mass: invalid neutron spacing.");

    IsotopeDistribution envelope;
    if (!ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
            pepComposition, envelope) ||
        envelope.vMass.empty() ||
        envelope.vMass.size() != envelope.vProb.size())
        throw std::runtime_error(
            "Cannot estimate precursor mass: invalid isotope envelope.");

    const auto apex = std::max_element(
        envelope.vProb.begin(), envelope.vProb.end());
    if (apex == envelope.vProb.end() || !std::isfinite(*apex))
        throw std::runtime_error(
            "Cannot estimate precursor mass: invalid isotope-envelope apex.");

    const size_t apexIndex = static_cast<size_t>(
        std::distance(envelope.vProb.begin(), apex));
    const double lightestMass = baseMass(pepComposition);
    estimate.nominalShift = static_cast<int>(std::lround(
        envelope.vMass[apexIndex] - lightestMass));
    estimate.mass = lightestMass +
        static_cast<double>(estimate.nominalShift) *
            estimate.neutronMass;
    if (!std::isfinite(estimate.mass))
        throw std::runtime_error(
            "Cannot estimate precursor mass: invalid modal estimate.");
    return estimate;
}

double PeptideIsotopeCalculator::calPrecursorMass(
    const string &peptideSequence)
{
	return calPrecursorEstimate(peptideSequence).mass;
}
