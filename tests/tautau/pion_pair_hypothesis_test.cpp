#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_pair_hypothesis.h"

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

void expect_equal(const tauamp::tautau::LinearComponents& left, const tauamp::tautau::LinearComponents& right) {
    assert(close(left.sm, right.sm));
    assert(close(left.f2_real, right.f2_real));
    assert(close(left.f2_imaginary, right.f2_imaginary));
    assert(close(left.f3_real, right.f3_real));
    assert(close(left.f3_imaginary, right.f3_imaginary));
}
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::ElectroweakParameters;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::KinematicPointOrigin;
    using tauamp::tautau::PionPairKinematicEnsemble;
    using tauamp::tautau::PionPairKinematicHypothesis;
    using tauamp::tautau::PionPairMatrixElement;
    using tauamp::tautau::TauPairKinematicPoint;
    using tauamp::tautau::WeightedPionPairKinematicHypothesis;

    constexpr double mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double momentum = std::sqrt(energy * energy - mass * mass);
    const FourMomentum electron(0.0, 0.0, energy, energy);
    const FourMomentum positron(0.0, 0.0, -energy, energy);
    const FourMomentum tau_minus(momentum * std::sin(theta), 0.0, momentum * std::cos(theta), energy);
    const FourMomentum tau_plus(-tau_minus.px(), 0.0, -tau_minus.pz(), energy);
    const TauPairKinematicPoint point(BeamState(electron, positron, 0.37, -0.41), tau_minus, tau_plus, mass);
    const FourMomentum pion_minus(0.2, -0.3, 0.4, 1.0);
    const FourMomentum pion_plus(-0.4, 0.1, -0.2, 1.0);
    const PionPairKinematicHypothesis hypothesis{point, pion_minus, pion_plus, KinematicPointOrigin::truth};
    const PionPairMatrixElement matrix_element(ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952});

    expect_equal(matrix_element.components(hypothesis), matrix_element.components(point, pion_minus, pion_plus));

    expect_invalid_argument([&] { PionPairKinematicEnsemble({}); });
    expect_invalid_argument([&] {
        PionPairKinematicEnsemble({WeightedPionPairKinematicHypothesis{hypothesis, -1.0}});
    });
    expect_invalid_argument([&] {
        PionPairKinematicEnsemble({WeightedPionPairKinematicHypothesis{hypothesis, std::numeric_limits<double>::quiet_NaN()}});
    });
    expect_invalid_argument([&] {
        PionPairKinematicEnsemble({WeightedPionPairKinematicHypothesis{hypothesis, std::numeric_limits<double>::infinity()}});
    });
    expect_invalid_argument([&] {
        PionPairKinematicEnsemble({WeightedPionPairKinematicHypothesis{hypothesis, 0.0}});
    });

    const PionPairKinematicEnsemble ensemble({WeightedPionPairKinematicHypothesis{hypothesis, 0.25},
                                               WeightedPionPairKinematicHypothesis{hypothesis, 0.75}});
    assert(ensemble.entries().size() == 2U);
    assert(close(ensemble.total_weight(), 1.0));
}
