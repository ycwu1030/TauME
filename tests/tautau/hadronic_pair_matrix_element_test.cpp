#include <array>
#include <cassert>
#include <cmath>

#include "tauamp/tautau/hadronic_pair_matrix_element.h"
#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_analyser.h"

namespace {
constexpr double tolerance = 2.0e-10;
bool close(double left, double right) { return std::abs(left - right) < tolerance; }

tauamp::tautau::FourMomentum pion_from_tnl_direction(const std::array<double, 3>& direction,
                                                      const tauamp::tautau::FourMomentum& tau_cm, double theta) {
    const std::array<double, 3> transverse{{-std::cos(theta), 0.0, std::sin(theta)}};
    const std::array<double, 3> normal{{0.0, -1.0, 0.0}};
    const std::array<double, 3> longitudinal{{std::sin(theta), 0.0, std::cos(theta)}};
    const tauamp::tautau::FourMomentum pion_rest(
        direction[0] * transverse[0] + direction[1] * normal[0] + direction[2] * longitudinal[0],
        direction[0] * transverse[1] + direction[1] * normal[1] + direction[2] * longitudinal[1],
        direction[0] * transverse[2] + direction[1] * normal[2] + direction[2] * longitudinal[2], 1.0);
    return pion_rest.boosted({{tau_cm.px() / tau_cm.energy(), tau_cm.py() / tau_cm.energy(),
                               tau_cm.pz() / tau_cm.energy()}});
}
}

int main() {
    using namespace tauamp::tautau;
    constexpr double mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double momentum = std::sqrt(energy * energy - mass * mass);
    const FourMomentum electron(0.0, 0.0, energy, energy);
    const FourMomentum positron(0.0, 0.0, -energy, energy);
    const FourMomentum tau_minus(momentum * std::sin(theta), 0.0, momentum * std::cos(theta), energy);
    const FourMomentum tau_plus(-tau_minus.px(), 0.0, -tau_minus.pz(), energy);
    const BeamState beams(electron, positron, 0.37, -0.41);
    const TauPairKinematicPoint kinematics(beams, tau_minus, tau_plus, mass);
    const std::array<double, 3> minus_direction{{0.2, -0.3, std::sqrt(0.87)}};
    const std::array<double, 3> plus_direction{{-0.4, 0.1, std::sqrt(0.83)}};
    const FourMomentum pion_minus = pion_from_tnl_direction(minus_direction, tau_minus, theta);
    const FourMomentum pion_plus = pion_from_tnl_direction(plus_direction, tau_plus, theta);
    const FourMomentum neutrino_minus = tau_minus - pion_minus;
    const FourMomentum neutrino_plus = tau_plus - pion_plus;
    const auto generic_minus = hadronic_analyser_tnl(kinematics, neutrino_minus, complex_scaled(pion_minus, 1.0), false);
    const auto generic_plus = hadronic_analyser_tnl(kinematics, neutrino_plus, complex_scaled(pion_plus, 1.0), true);
    const auto pion_minus_analyser = tau_minus_pion_analyser(kinematics, pion_minus);
    const auto pion_plus_analyser = tau_plus_pion_analyser(kinematics, pion_plus);
    for (int index = 0; index < 3; ++index) {
        assert(close(generic_minus[index], pion_minus_analyser[index]));
        assert(close(generic_plus[index], pion_plus_analyser[index]));
    }
    const ElectroweakParameters parameters{0.313, 0.48, 91.1876, 2.4952};
    const PionPairMatrixElement pion_matrix(parameters);
    const HadronicPairMatrixElement generic_matrix(parameters);
    const auto expected = pion_matrix.polynomial_components(kinematics, pion_minus, pion_plus);
    const auto actual = generic_matrix.polynomial_components(kinematics, generic_minus, generic_plus);
    const std::array<double, 15> expected_values{{expected.sm, expected.f2_real, expected.f2_imaginary, expected.f3_real, expected.f3_imaginary, expected.f2_real_f2_real, expected.f2_real_f2_imaginary, expected.f2_real_f3_real, expected.f2_real_f3_imaginary, expected.f2_imaginary_f2_imaginary, expected.f2_imaginary_f3_real, expected.f2_imaginary_f3_imaginary, expected.f3_real_f3_real, expected.f3_real_f3_imaginary, expected.f3_imaginary_f3_imaginary}};
    const std::array<double, 15> actual_values{{actual.sm, actual.f2_real, actual.f2_imaginary, actual.f3_real, actual.f3_imaginary, actual.f2_real_f2_real, actual.f2_real_f2_imaginary, actual.f2_real_f3_real, actual.f2_real_f3_imaginary, actual.f2_imaginary_f2_imaginary, actual.f2_imaginary_f3_real, actual.f2_imaginary_f3_imaginary, actual.f3_real_f3_real, actual.f3_real_f3_imaginary, actual.f3_imaginary_f3_imaginary}};
    for (int index = 0; index < 15; ++index) assert(close(expected_values[index], actual_values[index]));
}
