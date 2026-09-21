#include <cassert>
#include <cmath>

#include "tauamp/tautau/hadronic_analyser.h"

int main() {
    using namespace tauamp::tautau;
    const double mass = 1.77686;
    const FourMomentum neutrino{0.2, -0.1, 0.3, 0.4};
    const FourMomentum charged{0.3, 0.1, -0.2, 0.5};
    const FourMomentum neutral{-0.1, 0.2, 0.1, 0.4};
    const auto current = rho_kw_current(charged, neutral);
    const auto coefficients = decay_rest_coefficients(mass, neutrino, current);
    assert(coefficients.rate > 0.0);
    const auto analyser = normalized_analyser(coefficients);
    const double norm = std::sqrt(analyser[0] * analyser[0] + analyser[1] * analyser[1] + analyser[2] * analyser[2]);
    assert(norm <= 1.0 + 1.0e-9);

    auto scaled_current = current;
    for (auto& component : scaled_current) component *= std::complex<double>(2.0, 1.0);
    const auto scaled = decay_rest_coefficients(mass, neutrino, scaled_current);
    assert(std::abs(scaled.rate / coefficients.rate - 5.0) < 1.0e-8);
    const auto scaled_analyser = normalized_analyser(scaled);
    for (int index = 0; index < 3; ++index) assert(std::abs(scaled_analyser[index] - analyser[index]) < 1.0e-8);

    const auto anti = decay_rest_coefficients(mass, neutrino, current, true);
    assert(std::abs(anti.rate - coefficients.rate) < 1.0e-8);
    for (int index = 0; index < 3; ++index) assert(std::abs(anti.spin[index] + coefficients.spin[index]) < 1.0e-8);

    const FourMomentum same_charge_1{0.2, 0.1, 0.0, 0.4};
    const FourMomentum same_charge_2{-0.1, 0.0, 0.2, 0.35};
    const FourMomentum opposite_charge{0.15, -0.1, -0.05, 0.3};
    const auto a1_current = a1_kw_current(same_charge_1, same_charge_2, opposite_charge);
    const auto a1_coefficients = decay_rest_coefficients(mass, neutrino, a1_current);
    assert(a1_coefficients.rate > 0.0);
    const auto a1_analyser = normalized_analyser(a1_coefficients);
    const double a1_norm = std::sqrt(a1_analyser[0] * a1_analyser[0] + a1_analyser[1] * a1_analyser[1] + a1_analyser[2] * a1_analyser[2]);
    assert(a1_norm <= 1.0 + 1.0e-9);
}
