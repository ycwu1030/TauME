#include <cassert>
#include <cmath>
#include <stdexcept>

#include "tauamp/tautau/kinematics.h"

namespace {
constexpr double kTolerance = 1e-10;

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
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::TauPairKinematicPoint;

    constexpr double mass = 1.77;
    constexpr double beam_energy = 5.0;
    constexpr double theta = 0.7;
    const double tau_momentum = std::sqrt(beam_energy * beam_energy - mass * mass);
    const FourMomentum electron_cm(0.0, 0.0, beam_energy, beam_energy);
    const FourMomentum positron_cm(0.0, 0.0, -beam_energy, beam_energy);
    const FourMomentum tau_minus_cm(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), beam_energy);
    const FourMomentum tau_plus_cm(-tau_minus_cm.px(), 0.0, -tau_minus_cm.pz(), beam_energy);
    const BeamState symmetric_beams(electron_cm, positron_cm, 0.37, -0.41);
    const TauPairKinematicPoint symmetric(symmetric_beams, tau_minus_cm, tau_plus_cm, mass);

    const std::array<double, 3> lab_boost{{0.0, 0.0, std::tanh(0.35)}};
    const BeamState asymmetric_beams(electron_cm.boosted(lab_boost), positron_cm.boosted(lab_boost), 0.37, -0.41);
    const TauPairKinematicPoint asymmetric(asymmetric_beams, tau_minus_cm.boosted(lab_boost), tau_plus_cm.boosted(lab_boost), mass);

    assert(close(symmetric.sqrt_s(), asymmetric.sqrt_s()));
    assert(close(symmetric.pair_cm().tau_minus.px(), asymmetric.pair_cm().tau_minus.px()));
    assert(close(symmetric.pair_cm().tau_minus.py(), asymmetric.pair_cm().tau_minus.py()));
    assert(close(symmetric.pair_cm().tau_minus.pz(), asymmetric.pair_cm().tau_minus.pz()));
    assert(close(symmetric.pair_cm().tau_plus.px(), asymmetric.pair_cm().tau_plus.px()));
    assert(close(symmetric.electron_polarization(), 0.37));
    assert(close(symmetric.positron_polarization(), -0.41));

    expect_invalid_argument([&] { BeamState(electron_cm, positron_cm, 1.01, 0.0); });
    expect_invalid_argument([&] { BeamState(electron_cm, positron_cm, 0.0, -1.01); });
    expect_invalid_argument([&] { TauPairKinematicPoint(symmetric_beams, tau_minus_cm, FourMomentum(0.0, 0.0, 0.0, mass), mass); });
    expect_invalid_argument([&] { BeamState(FourMomentum(0.2, 0.0, beam_energy, beam_energy), positron_cm, 0.0, 0.0); });
}
