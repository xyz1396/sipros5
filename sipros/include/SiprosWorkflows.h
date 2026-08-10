#pragma once

namespace sipros
{

// Raxport uses charge 0 when it could not assign a peak charge. Such peaks
// remain eligible, but a peak assigned to another charge must not substitute
// for the requested theoretical ion.
constexpr bool observedPeakChargeMatches(int observedCharge,
									 int theoreticalCharge) noexcept
{
	return observedCharge == 0 || observedCharge == theoreticalCharge;
}

class SearchSpectraWorkflow
{
public:
	int run(int argc, char **argv);
};

} // namespace sipros

class DatabaseSearchWorkflow
{
public:
	int run(int argc, char **argv);
};

class TheoreticalSpectraWorkflow
{
public:
	int run(int argc, char **argv);
};

class ExperimentalSpectraWorkflow
{
public:
	int run(int argc, char **argv);
};
