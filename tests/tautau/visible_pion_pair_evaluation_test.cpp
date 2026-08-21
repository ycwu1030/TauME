#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_pair_marginalization.h"
#include "tauamp/tautau/visible_pion_pair_evaluation.h"

namespace {
constexpr double kTolerance = 3e-11;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

void expect_equal(const tauamp::tautau::LinearComponents& left, const tauamp::tautau::LinearComponents& right) {
    assert(close(left.sm, right.sm));
    assert(close(left.f2_real, right.f2_real));
    assert(close(left.f2_imaginary, right.f2_imaginary));
    assert(close(left.f3_real, right.f3_real));
    assert(close(left.f3_imaginary, right.f3_imaginary));
}

void expect_equal(const tauamp::tautau::LinearObservableComponents& left,
                  const tauamp::tautau::LinearObservableComponents& right) {
    assert(close(left.f2_real, right.f2_real));
    assert(close(left.f2_imaginary, right.f2_imaginary));
    assert(close(left.f3_real, right.f3_real));
    assert(close(left.f3_imaginary, right.f3_imaginary));
}

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

tauamp::tautau::FourMomentum rotate_y(const tauamp::tautau::FourMomentum& momentum, double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return {cosine * momentum.px() + sine * momentum.pz(), momentum.py(),
            -sine * momentum.px() + cosine * momentum.pz(), momentum.energy()};
}

tauamp::tautau::FourMomentum visible_from_tau_rest(const std::array<double, 3>& direction,
                                                    const tauamp::tautau::FourMomentum& tau, double tau_mass,
                                                    double visible_mass) {
    const double rest_energy = (tau_mass * tau_mass + visible_mass * visible_mass) / (2.0 * tau_mass);
    const double rest_momentum = (tau_mass * tau_mass - visible_mass * visible_mass) / (2.0 * tau_mass);
    return tauamp::tautau::FourMomentum(rest_momentum * direction[0], rest_momentum * direction[1],
                                         rest_momentum * direction[2], rest_energy)
        .boosted({{tau.px() / tau.energy(), tau.py() / tau.energy(), tau.pz() / tau.energy()}});
}

std::vector<tauamp::tautau::WeightedLinearComponents> weighted_components(
    const tauamp::tautau::PionPairVisibleEventEvaluation& evaluation) {
    std::vector<tauamp::tautau::WeightedLinearComponents> result;
    for (const auto& entry : evaluation.entries()) result.push_back({entry.components, entry.weight});
    return result;
}
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::ElectroweakParameters;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::HadronicTauPairSolutionStatus;
    using tauamp::tautau::KinematicBranchCombination;
    using tauamp::tautau::KinematicPointOrigin;
    using tauamp::tautau::PionPairMatrixElement;
    using tauamp::tautau::PionPairVisibleEventEvaluation;
    using tauamp::tautau::PionPairVisibleEventEvaluator;

    constexpr double tau_mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double pion_mass = 0.13957;
    constexpr double theta = 0.91;
    const double tau_momentum = std::sqrt(energy * energy - tau_mass * tau_mass);
    const FourMomentum electron_cm(0.0, 0.0, energy, energy);
    const FourMomentum positron_cm(0.0, 0.0, -energy, energy);
    const BeamState beams_cm(electron_cm, positron_cm, 0.37, -0.41);
    const FourMomentum tau_minus_cm(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), energy);
    const FourMomentum tau_plus_cm(-tau_minus_cm.px(), -tau_minus_cm.py(), -tau_minus_cm.pz(), energy);
    const FourMomentum pion_minus =
        visible_from_tau_rest({{0.2, -0.3, std::sqrt(0.87)}}, tau_minus_cm, tau_mass, pion_mass);
    const FourMomentum pion_plus =
        visible_from_tau_rest({{-0.4, 0.1, std::sqrt(0.83)}}, tau_plus_cm, tau_mass, pion_mass);
    const ElectroweakParameters parameters{0.313, 0.48, 91.1876, 2.4952};
    const PionPairVisibleEventEvaluator evaluator(parameters);
    const PionPairMatrixElement direct_matrix_element(parameters);

    const auto twofold = evaluator.evaluate(beams_cm, pion_minus, pion_plus, tau_mass);
    assert(twofold.status() == HadronicTauPairSolutionStatus::two_solutions);
    assert(twofold.entries().size() == 2U);
    for (const auto& entry : twofold.entries()) {
        assert(entry.hypothesis.origin == KinematicPointOrigin::analytic_branch);
        assert(close(entry.weight, 0.5));
        expect_equal(entry.components, direct_matrix_element.components(entry.hypothesis));
    }
    const auto expected_average = tauamp::tautau::average_conditional_observables(weighted_components(twofold));
    const auto expected_ratio = tauamp::tautau::ratio_of_average_components(weighted_components(twofold));
    assert(!close(expected_average.f3_real, expected_ratio.f3_real));
    expect_equal(twofold.observable_components(KinematicBranchCombination::average_conditional_observables), expected_average);
    expect_equal(twofold.observable_components(KinematicBranchCombination::ratio_of_average_components), expected_ratio);

    auto reversed_entries = twofold.entries();
    std::reverse(reversed_entries.begin(), reversed_entries.end());
    const PionPairVisibleEventEvaluation reversed(twofold.status(), std::move(reversed_entries));
    expect_equal(reversed.observable_components(KinematicBranchCombination::average_conditional_observables), expected_average);
    expect_equal(reversed.observable_components(KinematicBranchCombination::ratio_of_average_components), expected_ratio);

    const std::array<double, 3> lab_boost{{0.0, 0.0, std::tanh(0.35)}};
    const auto boosted = evaluator.evaluate(BeamState(electron_cm.boosted(lab_boost), positron_cm.boosted(lab_boost), 0.37, -0.41),
                                            pion_minus.boosted(lab_boost), pion_plus.boosted(lab_boost), tau_mass);
    assert(boosted.status() == HadronicTauPairSolutionStatus::two_solutions);
    assert(boosted.entries().size() == twofold.entries().size());
    for (std::size_t index = 0; index < twofold.entries().size(); ++index)
        expect_equal(boosted.entries()[index].components, twofold.entries()[index].components);

    const FourMomentum tau_minus_forward(0.0, 0.0, tau_momentum, energy);
    const FourMomentum tau_plus_forward(0.0, 0.0, -tau_momentum, energy);
    constexpr double tangent_production_angle = 0.91;
    const auto tangent = evaluator.evaluate(
        BeamState(electron_cm, positron_cm, 0.0, 0.0),
        rotate_y(visible_from_tau_rest({{std::sin(0.7), 0.0, std::cos(0.7)}}, tau_minus_forward, tau_mass, pion_mass),
                 tangent_production_angle),
        rotate_y(visible_from_tau_rest({{std::sin(1.0), 0.0, std::cos(1.0)}}, tau_plus_forward, tau_mass, pion_mass),
                 tangent_production_angle), tau_mass);
    assert(tangent.status() == HadronicTauPairSolutionStatus::one_solution);
    assert(tangent.entries().size() == 1U);
    assert(close(tangent.entries().front().weight, 1.0));
    expect_equal(tangent.observable_components(KinematicBranchCombination::average_conditional_observables),
                 tangent.observable_components(KinematicBranchCombination::ratio_of_average_components));

    const auto no_solution = evaluator.evaluate(beams_cm, FourMomentum(0.0, 0.0, 0.0, 1.0), pion_plus, tau_mass);
    assert(no_solution.status() == HadronicTauPairSolutionStatus::no_solution);
    assert(no_solution.entries().empty());
    expect_invalid_argument([&] {
        no_solution.observable_components(KinematicBranchCombination::ratio_of_average_components);
    });

    const double cosine = 0.5;
    const double visible_minus_energy = tau_mass * tau_mass / (2.0 * (energy - tau_momentum * cosine));
    const double visible_plus_energy = tau_mass * tau_mass / (2.0 * (energy + tau_momentum * cosine));
    const auto non_unique = evaluator.evaluate(
        BeamState(electron_cm, positron_cm, 0.0, 0.0), FourMomentum(0.0, 0.0, visible_minus_energy, visible_minus_energy),
        FourMomentum(0.0, 0.0, visible_plus_energy, visible_plus_energy), tau_mass);
    assert(non_unique.status() == HadronicTauPairSolutionStatus::non_unique);
    assert(non_unique.entries().empty());
    expect_invalid_argument([&] {
        non_unique.observable_components(KinematicBranchCombination::average_conditional_observables);
    });
}
