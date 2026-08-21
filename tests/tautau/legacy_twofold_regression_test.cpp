#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "legacy_twofold_control.h"
#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"

namespace {
constexpr double kTolerance = 5e-9;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

bool same_momentum(const tauamp::tautau::FourMomentum& left, const tauamp::tautau::FourMomentum& right) {
    return close(left.px(), right.px()) && close(left.py(), right.py()) && close(left.pz(), right.pz()) &&
           close(left.energy(), right.energy());
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

bool matches(const legacy_twofold_control::NeutrinoPair& legacy,
             const tauamp::tautau::HadronicTauPairKinematicSolution& current,
             const tauamp::tautau::FourMomentum& visible_minus, const tauamp::tautau::FourMomentum& visible_plus) {
    return same_momentum(legacy.neutrino_minus, current.neutrino_minus_lab) &&
           same_momentum(legacy.neutrino_plus, current.neutrino_plus_lab) &&
           same_momentum(visible_minus + legacy.neutrino_minus, current.point.tau_minus_lab()) &&
           same_momentum(visible_plus + legacy.neutrino_plus, current.point.tau_plus_lab());
}
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::HadronicTauPairKinematicSolver;
    using tauamp::tautau::HadronicTauPairObservation;
    using tauamp::tautau::HadronicTauPairSolutionStatus;

    constexpr double tau_mass = 1.77;  // Legacy TauME MTAU in include/tauamp/constants.h.
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double tau_momentum = std::sqrt(energy * energy - tau_mass * tau_mass);
    const FourMomentum electron(0.0, 0.0, energy, energy);
    const FourMomentum positron(0.0, 0.0, -energy, energy);
    const FourMomentum tau_minus(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), energy);
    const FourMomentum tau_plus(-tau_minus.px(), -tau_minus.py(), -tau_minus.pz(), energy);
    const FourMomentum visible_minus =
        visible_from_tau_rest({{0.2, -0.3, std::sqrt(0.87)}}, tau_minus, tau_mass, 0.6);
    const FourMomentum visible_plus =
        visible_from_tau_rest({{-0.4, 0.1, std::sqrt(0.83)}}, tau_plus, tau_mass, 0.8);

    const auto current = HadronicTauPairKinematicSolver::solve(
        {BeamState(electron, positron, 0.0, 0.0), visible_minus, visible_plus, tau_mass});
    assert(current.status == HadronicTauPairSolutionStatus::two_solutions);
    assert(current.solutions.size() == 2U);

    std::vector<legacy_twofold_control::NeutrinoPair> legacy;
    assert(legacy_twofold_control::reconstruct_neutrinos(2.0 * energy, visible_minus, visible_plus, legacy));
    assert(legacy.size() == 2U);
    assert((matches(legacy[0], current.solutions[0], visible_minus, visible_plus) &&
            matches(legacy[1], current.solutions[1], visible_minus, visible_plus)) ||
           (matches(legacy[0], current.solutions[1], visible_minus, visible_plus) &&
            matches(legacy[1], current.solutions[0], visible_minus, visible_plus)));

    constexpr double second_theta = 1.2;
    const FourMomentum second_tau_minus(tau_momentum * std::sin(second_theta), 0.0,
                                        tau_momentum * std::cos(second_theta), energy);
    const FourMomentum second_tau_plus(-second_tau_minus.px(), -second_tau_minus.py(), -second_tau_minus.pz(), energy);
    const FourMomentum second_visible_minus =
        visible_from_tau_rest({{0.5, 0.4, std::sqrt(0.59)}}, second_tau_minus, tau_mass, 0.7);
    const FourMomentum second_visible_plus =
        visible_from_tau_rest({{-0.25, 0.35, std::sqrt(0.815)}}, second_tau_plus, tau_mass, 0.5);
    const auto second_current = HadronicTauPairKinematicSolver::solve(
        {BeamState(electron, positron, 0.0, 0.0), second_visible_minus, second_visible_plus, tau_mass});
    assert(second_current.status == HadronicTauPairSolutionStatus::two_solutions);
    assert(second_current.solutions.size() == 2U);

    std::vector<legacy_twofold_control::NeutrinoPair> second_legacy;
    assert(legacy_twofold_control::reconstruct_neutrinos(2.0 * energy, second_visible_minus, second_visible_plus, second_legacy));
    assert(second_legacy.size() == 2U);
    assert((matches(second_legacy[0], second_current.solutions[0], second_visible_minus, second_visible_plus) &&
            matches(second_legacy[1], second_current.solutions[1], second_visible_minus, second_visible_plus)) ||
           (matches(second_legacy[0], second_current.solutions[1], second_visible_minus, second_visible_plus) &&
            matches(second_legacy[1], second_current.solutions[0], second_visible_minus, second_visible_plus)));
}
