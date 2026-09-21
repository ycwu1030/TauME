#include <cassert>
#include <cmath>
#include <complex>

#include "tauamp/tautau/rho_a1_currents.h"

int main() {
    using namespace tauamp::tautau;
    const FourMomentum charged(0.32, 0.10, 0.61, 0.72);
    const FourMomentum neutral(-0.15, 0.08, -0.48, 0.55);
    const auto rho = rho_kw_current(charged, neutral);
    assert(std::abs(minkowski_dot(charged + neutral, rho)) < 1.0e-12);
    assert(rho_kw_lineshape(0.64) != rho_taudecay_lineshape(0.64));

    const FourMomentum q1(0.22, 0.10, 0.47, 0.55);
    const FourMomentum q2(-0.18, 0.16, -0.39, 0.48);
    const FourMomentum qop(-0.04, -0.26, 0.51, 0.62);
    const auto a1 = a1_kw_current(q1, q2, qop);
    const auto swapped = a1_kw_current(q2, q1, qop);
    const auto conjugate = a1_kw_current(q1, q2, qop, true);
    assert(std::abs(minkowski_dot(q1 + q2 + qop, a1)) < 1.0e-12);
    for (std::size_t index = 0; index < a1.size(); ++index) {
        assert(std::abs(a1[index] - swapped[index]) < 1.0e-14);
        assert(conjugate[index] == std::conj(a1[index]));
    }
}
