#ifndef TAUAMP_TAUTAU_HADRONIC_TAU_PAIR_KINEMATICS_H_
#define TAUAMP_TAUTAU_HADRONIC_TAU_PAIR_KINEMATICS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/kinematics.h"

namespace tauamp::tautau {

struct HadronicTauPairObservation {
    BeamState beams;
    FourMomentum visible_minus_lab;
    FourMomentum visible_plus_lab;
    double tau_mass;
};

enum class HadronicTauPairSolutionStatus { no_solution, one_solution, two_solutions, non_unique };

struct HadronicTauPairKinematicSolution {
    TauPairKinematicPoint point;
    FourMomentum neutrino_minus_lab;
    FourMomentum neutrino_plus_lab;
};

struct HadronicTauPairSolutionSet {
    HadronicTauPairSolutionStatus status;
    std::vector<HadronicTauPairKinematicSolution> solutions;
};

class HadronicTauPairKinematicSolver {
public:
    static HadronicTauPairSolutionSet solve(const HadronicTauPairObservation& observation) {
        constexpr double tolerance = 1e-9;
        if (observation.tau_mass <= 0.0 || observation.beams.sqrt_s() <= 2.0 * observation.tau_mass)
            return no_solution();

        const FourMomentum beam_total = observation.beams.electron_lab() + observation.beams.positron_lab();
        const std::array<double, 3> to_pair_cm{{-beam_total.px() / beam_total.energy(), -beam_total.py() / beam_total.energy(),
                                                 -beam_total.pz() / beam_total.energy()}};
        const std::array<double, 3> to_lab{{-to_pair_cm[0], -to_pair_cm[1], -to_pair_cm[2]}};
        const FourMomentum visible_minus_cm = observation.visible_minus_lab.boosted(to_pair_cm);
        const FourMomentum visible_plus_cm = observation.visible_plus_lab.boosted(to_pair_cm);
        const std::array<double, 3> minus_spatial = detail::spatial_components(visible_minus_cm);
        const std::array<double, 3> plus_spatial = detail::spatial_components(visible_plus_cm);
        const double minus_momentum = detail::spatial_norm(minus_spatial);
        const double plus_momentum = detail::spatial_norm(plus_spatial);
        if (minus_momentum < tolerance || plus_momentum < tolerance) return no_solution();

        const double tau_energy = observation.beams.sqrt_s() / 2.0;
        const double tau_momentum_squared = tau_energy * tau_energy - observation.tau_mass * observation.tau_mass;
        if (tau_momentum_squared <= tolerance * tolerance) return no_solution();
        const double tau_momentum = std::sqrt(tau_momentum_squared);
        const double minus_mass_squared = visible_minus_cm.mass_squared();
        const double plus_mass_squared = visible_plus_cm.mass_squared();
        if (minus_mass_squared < -tolerance || plus_mass_squared < -tolerance) return no_solution();

        const std::array<double, 3> minus_direction{{minus_spatial[0] / minus_momentum, minus_spatial[1] / minus_momentum,
                                                      minus_spatial[2] / minus_momentum}};
        const std::array<double, 3> plus_direction{{plus_spatial[0] / plus_momentum, plus_spatial[1] / plus_momentum,
                                                    plus_spatial[2] / plus_momentum}};
        const double c_minus = (2.0 * tau_energy * visible_minus_cm.energy() - observation.tau_mass * observation.tau_mass -
                                minus_mass_squared) /
                               (2.0 * tau_momentum * minus_momentum);
        const double c_plus = (observation.tau_mass * observation.tau_mass + plus_mass_squared -
                               2.0 * tau_energy * visible_plus_cm.energy()) /
                              (2.0 * tau_momentum * plus_momentum);
        if (c_minus < -1.0 - tolerance || c_minus > 1.0 + tolerance || c_plus < -1.0 - tolerance ||
            c_plus > 1.0 + tolerance)
            return no_solution();
        const double bounded_c_minus = bound_unit(c_minus);
        const double bounded_c_plus = bound_unit(c_plus);
        const double dot = detail::spatial_dot(minus_direction, plus_direction);
        const double denominator = std::max(0.0, 1.0 - dot * dot);

        if (denominator <= tolerance * tolerance) {
            if (std::abs(bounded_c_plus - dot * bounded_c_minus) > tolerance) return no_solution();
            if (std::abs(std::abs(bounded_c_minus) - 1.0) <= tolerance) {
                const std::array<double, 3> direction{{bounded_c_minus * minus_direction[0], bounded_c_minus * minus_direction[1],
                                                        bounded_c_minus * minus_direction[2]}};
                return one_solution(make_solution(observation, visible_minus_cm, visible_plus_cm, direction, tau_energy,
                                                  tau_momentum, to_lab, tolerance));
            }
            return {HadronicTauPairSolutionStatus::non_unique, {}};
        }

        const std::array<double, 3> in_plane{{(bounded_c_minus - dot * bounded_c_plus) * minus_direction[0] +
                                                   (bounded_c_plus - dot * bounded_c_minus) * plus_direction[0],
                                               (bounded_c_minus - dot * bounded_c_plus) * minus_direction[1] +
                                                   (bounded_c_plus - dot * bounded_c_minus) * plus_direction[1],
                                               (bounded_c_minus - dot * bounded_c_plus) * minus_direction[2] +
                                                   (bounded_c_plus - dot * bounded_c_minus) * plus_direction[2]}};
        const std::array<double, 3> base{{in_plane[0] / denominator, in_plane[1] / denominator, in_plane[2] / denominator}};
        const double discriminant = 1.0 - detail::spatial_dot(base, base);
        if (discriminant < -tolerance) return no_solution();
        if (discriminant <= tolerance) {
            return one_solution(make_solution(observation, visible_minus_cm, visible_plus_cm,
                                              detail::normalized(base, "tangent tau direction cannot vanish"), tau_energy,
                                              tau_momentum, to_lab, tolerance));
        }

        const std::array<double, 3> cross = detail::spatial_cross(minus_direction, plus_direction);
        const std::array<double, 3> normal{{cross[0] / std::sqrt(denominator), cross[1] / std::sqrt(denominator),
                                             cross[2] / std::sqrt(denominator)}};
        const double scale = std::sqrt(discriminant);
        const std::array<double, 3> first{{base[0] + scale * normal[0], base[1] + scale * normal[1], base[2] + scale * normal[2]}};
        const std::array<double, 3> second{{base[0] - scale * normal[0], base[1] - scale * normal[1], base[2] - scale * normal[2]}};
        HadronicTauPairSolutionSet result{HadronicTauPairSolutionStatus::two_solutions, {}};
        result.solutions.push_back(
            make_solution(observation, visible_minus_cm, visible_plus_cm, first, tau_energy, tau_momentum, to_lab, tolerance));
        result.solutions.push_back(
            make_solution(observation, visible_minus_cm, visible_plus_cm, second, tau_energy, tau_momentum, to_lab, tolerance));
        return result;
    }

private:
    static double bound_unit(double value) { return value < -1.0 ? -1.0 : (value > 1.0 ? 1.0 : value); }

    static HadronicTauPairSolutionSet no_solution() { return {HadronicTauPairSolutionStatus::no_solution, {}}; }

    static HadronicTauPairSolutionSet one_solution(HadronicTauPairKinematicSolution solution) {
        return {HadronicTauPairSolutionStatus::one_solution, {solution}};
    }

    static HadronicTauPairKinematicSolution make_solution(const HadronicTauPairObservation& observation,
                                                           const FourMomentum& visible_minus_cm,
                                                           const FourMomentum& visible_plus_cm,
                                                           const std::array<double, 3>& direction, double tau_energy,
                                                           double tau_momentum, const std::array<double, 3>& to_lab,
                                                           double tolerance) {
        const FourMomentum tau_minus_cm(tau_momentum * direction[0], tau_momentum * direction[1], tau_momentum * direction[2],
                                         tau_energy);
        const FourMomentum tau_plus_cm(-tau_minus_cm.px(), -tau_minus_cm.py(), -tau_minus_cm.pz(), tau_energy);
        const FourMomentum neutrino_minus_cm = tau_minus_cm - visible_minus_cm;
        const FourMomentum neutrino_plus_cm = tau_plus_cm - visible_plus_cm;
        if (neutrino_minus_cm.energy() < -tolerance || neutrino_plus_cm.energy() < -tolerance ||
            std::abs(neutrino_minus_cm.mass_squared()) > tolerance || std::abs(neutrino_plus_cm.mass_squared()) > tolerance)
            throw std::runtime_error("cone intersection generated an invalid neutrino candidate");

        const FourMomentum tau_minus_lab = tau_minus_cm.boosted(to_lab);
        const FourMomentum tau_plus_lab = tau_plus_cm.boosted(to_lab);
        return {TauPairKinematicPoint(observation.beams, tau_minus_lab, tau_plus_lab, observation.tau_mass),
                neutrino_minus_cm.boosted(to_lab), neutrino_plus_cm.boosted(to_lab)};
    }
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_HADRONIC_TAU_PAIR_KINEMATICS_H_
