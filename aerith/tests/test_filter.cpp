#include "filter.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    assert(aerith::stripped_peptide("K[PEPTM~IDE]R") == "PEPTMIDE");
    assert(aerith::stripped_peptide("-.ACD[+57.0]EF.-") == "ACDEF");
    assert(aerith::stripped_peptide("[%M~PEPTIDE]R") == "MPEPTIDE");
    assert(aerith::stripped_peptide("K[PEPN!IDE]") == "PEPNIDE");
    assert(aerith::stripped_peptide("[PEPTIDE]R") == "PEPTIDE");

    const std::vector<double> scores{10.0, 9.0, 8.0, 7.0, 6.0};
    const std::vector<int> labels{1, 1, -1, 1, -1};
    const auto q = aerith::target_decoy_qvalues(scores, labels);
    assert(q.size() == scores.size());
    assert(std::abs(q[0] - 0.5) < 1e-12);
    assert(std::abs(q[1] - 0.5) < 1e-12);
    assert(q[2] <= q[4]);
    assert(q[3] <= q[4]);

    double pi0 = 0.0;
    const auto mixmax = aerith::mixmax_qvalues(scores, labels, &pi0);
    assert(mixmax.size() == scores.size());
    assert(pi0 >= 0.0 && pi0 <= 1.0);
    assert(mixmax[0] <= mixmax[1]);
    assert(mixmax[1] <= mixmax[3]);

    std::cout << "aerith unit tests passed\n";
    return 0;
}
