#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>

#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"

namespace {
constexpr double kTolerance = 5e-10;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

bool same_momentum(const tauamp::tautau::FourMomentum& left, const tauamp::tautau::FourMomentum& right) {
    return close(left.px(), right.px()) && close(left.py(), right.py()) && close(left.pz(), right.pz()) &&
           close(left.energy(), right.energy());
}

tauamp::tautau::FourMomentum visible_from_tau_rest(const std::array<double, 3>& direction,
                                                    const tauamp::tautau::FourMomentum& tau_cm,
                                                    double tau_mass, double visible_mass) {
    const double visible_energy = (tau_mass * tau_mass + visible_mass * visible_mass) / (2.0 * tau_mass);
    const double visible_momentum = (tau_mass * tau_mass - visible_mass * visible_mass) / (2.0 * tau_mass);
    const tauamp::tautau::FourMomentum visible_rest(direction[0] * visible_momentum, direction[1] * visible_momentum,
                                                     direction[2] * visible_momentum, visible_energy);
    return visible_rest.boosted({{tau_cm.px() / tau_cm.energy(), tau_cm.py() / tau_cm.energy(),
                                  tau_cm.pz() / tau_cm.energy()}});
}

void expect_solution_constraints(const tauamp::tautau::HadronicTauPairKinematicSolution& solution,
                                 const tauamp::tautau::HadronicTauPairObservation& observation) {
    const auto total = observation.beams.electron_lab() + observation.beams.positron_lab();
    const auto reconstructed_total = solution.point.tau_minus_lab() + solution.point.tau_plus_lab();
    assert(same_momentum(total, reconstructed_total));
    assert(close(solution.point.tau_minus_lab().mass_squared(), observation.tau_mass * observation.tau_mass));
    assert(close(solution.point.tau_plus_lab().mass_squared(), observation.tau_mass * observation.tau_mass));
    assert(same_momentum(solution.neutrino_minus_lab, solution.point.tau_minus_lab() - observation.visible_minus_lab));
    assert(same_momentum(solution.neutrino_plus_lab, solution.point.tau_plus_lab() - observation.visible_plus_lab));
    assert(close(solution.neutrino_minus_lab.mass_squared(), 0.0));
    assert(close(solution.neutrino_plus_lab.mass_squared(), 0.0));
    assert(solution.neutrino_minus_lab.energy() >= -kTolerance);
    assert(solution.neutrino_plus_lab.energy() >= -kTolerance);
}
}  // namespace

int main() {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::HadronicTauPairKinematicSolver;
    using tauamp::tautau::HadronicTauPairObservation;
    using tauamp::tautau::HadronicTauPairSolutionStatus;

    constexpr double tau_mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double tau_momentum = std::sqrt(energy * energy - tau_mass * tau_mass);
    const FourMomentum electron_cm(0.0, 0.0, energy, energy);
    const FourMomentum positron_cm(0.0, 0.0, -energy, energy);
    const FourMomentum tau_minus_cm(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), energy);
    const FourMomentum tau_plus_cm(-tau_minus_cm.px(), 0.0, -tau_minus_cm.pz(), energy);
    const FourMomentum visible_minus = visible_from_tau_rest({{0.2, -0.3, std::sqrt(0.87)}}, tau_minus_cm, tau_mass, 0.6);
    const FourMomentum visible_plus = visible_from_tau_rest({{-0.4, 0.1, std::sqrt(0.83)}}, tau_plus_cm, tau_mass, 0.8);
    const HadronicTauPairObservation observation{
        BeamState(electron_cm, positron_cm, 0.37, -0.41), visible_minus, visible_plus, tau_mass};

    const auto twofold = HadronicTauPairKinematicSolver::solve(observation);
    assert(twofold.status == HadronicTauPairSolutionStatus::two_solutions);
    assert(twofold.solutions.size() == 2U);
    for (const auto& solution : twofold.solutions) expect_solution_constraints(solution, observation);
    assert((same_momentum(twofold.solutions[0].point.tau_minus_lab(), tau_minus_cm) ||
            same_momentum(twofold.solutions[1].point.tau_minus_lab(), tau_minus_cm)));

    const std::array<double, 3> lab_boost{{0.0, 0.0, std::tanh(0.35)}};
    const HadronicTauPairObservation asymmetric{
        BeamState(electron_cm.boosted(lab_boost), positron_cm.boosted(lab_boost), 0.37, -0.41),
        visible_minus.boosted(lab_boost), visible_plus.boosted(lab_boost), tau_mass};
    const auto asymmetric_solutions = HadronicTauPairKinematicSolver::solve(asymmetric);
    assert(asymmetric_solutions.status == HadronicTauPairSolutionStatus::two_solutions);
    assert(asymmetric_solutions.solutions.size() == 2U);
    for (const auto& solution : asymmetric_solutions.solutions) expect_solution_constraints(solution, asymmetric);
    assert(same_momentum(asymmetric_solutions.solutions[0].point.pair_cm().tau_minus,
                         twofold.solutions[0].point.pair_cm().tau_minus));
    assert(same_momentum(asymmetric_solutions.solutions[1].point.pair_cm().tau_minus,
                         twofold.solutions[1].point.pair_cm().tau_minus));

    const FourMomentum tau_minus_forward(0.0, 0.0, tau_momentum, energy);
    const FourMomentum tau_plus_forward(0.0, 0.0, -tau_momentum, energy);
    const HadronicTauPairObservation tangent{
        BeamState(electron_cm, positron_cm, 0.0, 0.0),
        visible_from_tau_rest({{std::sin(0.7), 0.0, std::cos(0.7)}}, tau_minus_forward, tau_mass, 0.6),
        visible_from_tau_rest({{std::sin(1.0), 0.0, std::cos(1.0)}}, tau_plus_forward, tau_mass, 0.8), tau_mass};
    const auto tangent_solutions = HadronicTauPairKinematicSolver::solve(tangent);
    assert(tangent_solutions.status == HadronicTauPairSolutionStatus::one_solution);
    assert(tangent_solutions.solutions.size() == 1U);
    expect_solution_constraints(tangent_solutions.solutions.front(), tangent);
    assert(same_momentum(tangent_solutions.solutions.front().point.tau_minus_lab(), tau_minus_forward));

    const double c = 0.5;
    const double visible_minus_energy = tau_mass * tau_mass / (2.0 * (energy - tau_momentum * c));
    const double visible_plus_energy = tau_mass * tau_mass / (2.0 * (energy + tau_momentum * c));
    const HadronicTauPairObservation non_unique{
        BeamState(electron_cm, positron_cm, 0.0, 0.0), FourMomentum(0.0, 0.0, visible_minus_energy, visible_minus_energy),
        FourMomentum(0.0, 0.0, visible_plus_energy, visible_plus_energy), tau_mass};
    const auto non_unique_solutions = HadronicTauPairKinematicSolver::solve(non_unique);
    assert(non_unique_solutions.status == HadronicTauPairSolutionStatus::non_unique);
    assert(non_unique_solutions.solutions.empty());

    const HadronicTauPairObservation no_solution{
        BeamState(electron_cm, positron_cm, 0.0, 0.0), FourMomentum(0.0, 0.0, 0.0, 1.0), visible_plus, tau_mass};
    const auto no_solutions = HadronicTauPairKinematicSolver::solve(no_solution);
    assert(no_solutions.status == HadronicTauPairSolutionStatus::no_solution);
    assert(no_solutions.solutions.empty());
}
