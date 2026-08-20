#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>

#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_analyser.h"

namespace {
constexpr double kTolerance = 3e-11;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

template <class Callback>
void expect_invalid_argument(Callback callback) {
    bool thrown = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        thrown = true;
    }
    assert(thrown);
}

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
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::ElectroweakParameters;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::PhotonDipoleFormFactors;
    using tauamp::tautau::PionPairMatrixElement;
    using tauamp::tautau::TauPairKinematicPoint;

    constexpr double mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double momentum = std::sqrt(energy * energy - mass * mass);
    const FourMomentum electron_cm(0.0, 0.0, energy, energy);
    const FourMomentum positron_cm(0.0, 0.0, -energy, energy);
    const FourMomentum tau_minus_cm(momentum * std::sin(theta), 0.0, momentum * std::cos(theta), energy);
    const FourMomentum tau_plus_cm(-tau_minus_cm.px(), 0.0, -tau_minus_cm.pz(), energy);
    const BeamState beams(electron_cm, positron_cm, 0.37, -0.41);
    const TauPairKinematicPoint kinematics(beams, tau_minus_cm, tau_plus_cm, mass);

    const std::array<double, 3> pion_minus_direction{{0.2, -0.3, std::sqrt(0.87)}};
    const std::array<double, 3> pion_plus_direction{{-0.4, 0.1, std::sqrt(0.83)}};
    const FourMomentum pion_minus_cm = pion_from_tnl_direction(pion_minus_direction, tau_minus_cm, theta);
    const FourMomentum pion_plus_cm = pion_from_tnl_direction(pion_plus_direction, tau_plus_cm, theta);
    const auto minus_analyser = tauamp::tautau::tau_minus_pion_analyser(kinematics, pion_minus_cm);
    const auto plus_analyser = tauamp::tautau::tau_plus_pion_analyser(kinematics, pion_plus_cm);
    for (std::size_t index = 0; index < 3; ++index) {
        assert(close(minus_analyser[index], pion_minus_direction[index]));
        assert(close(plus_analyser[index], -pion_plus_direction[index]));
    }
    expect_invalid_argument([&] { tauamp::tautau::tau_minus_pion_analyser(kinematics, FourMomentum()); });

    const PionPairMatrixElement matrix_element(ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952});
    const auto components = matrix_element.components(kinematics, pion_minus_cm, pion_plus_cm);
    assert(close(components.sm, 0.020226557273898237));
    assert(close(components.f2_real, 0.05119425075109686));
    assert(close(components.f2_imaginary, -0.0007211197775460229));
    assert(close(components.f3_real, -0.0002487522803563775));
    assert(close(components.f3_imaginary, -0.019025560829261927));
    assert(close(components.recompose(PhotonDipoleFormFactors{}), components.sm));

    const std::array<double, 3> lab_boost{{0.0, 0.0, std::tanh(0.35)}};
    const BeamState asymmetric_beams(electron_cm.boosted(lab_boost), positron_cm.boosted(lab_boost), 0.37, -0.41);
    const TauPairKinematicPoint asymmetric_kinematics(
        asymmetric_beams, tau_minus_cm.boosted(lab_boost), tau_plus_cm.boosted(lab_boost), mass);
    const auto asymmetric_components = matrix_element.components(
        asymmetric_kinematics, pion_minus_cm.boosted(lab_boost), pion_plus_cm.boosted(lab_boost));
    assert(close(asymmetric_components.sm, components.sm));
    assert(close(asymmetric_components.f2_real, components.f2_real));
    assert(close(asymmetric_components.f2_imaginary, components.f2_imaginary));
    assert(close(asymmetric_components.f3_real, components.f3_real));
    assert(close(asymmetric_components.f3_imaginary, components.f3_imaginary));
}
