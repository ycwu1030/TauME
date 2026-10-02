#include <cassert>

#include "tauamp/delphes/rho_family_event_loop_v2.h"

int main() {
    using tauamp::delphes::g7_v2_all_charged;
    using tauamp::delphes::g7_v2_hardest_charged;
    using tauamp::tautau::ReconstructedPion;
    const std::vector<ReconstructedPion> charged{
        {0U, 211, 0.4, 0.0, 0.0, {0.1, 0.0, 0.0, 0.2}},
        {1U, 211, 0.9, 0.0, 0.0, {0.2, 0.0, 0.0, 0.3}},
        {2U, -211, 0.7, 0.0, 0.0, {-0.1, 0.0, 0.0, 0.2}}};
    const auto plus = g7_v2_all_charged(charged, 211);
    const auto hardest = g7_v2_hardest_charged(charged, 211);
    assert(plus.size() == 2U);
    assert(hardest.size() == 1U && hardest.front().index == 1U);
    const auto minus = g7_v2_hardest_charged(charged, -211);
    assert(minus.size() == 1U && minus.front().index == 2U);
}
