#ifndef TAUAMP_TAUTAU_VISIBLE_PION_PAIR_EVALUATION_H_
#define TAUAMP_TAUTAU_VISIBLE_PION_PAIR_EVALUATION_H_

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"
#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_pair_marginalization.h"

namespace tauamp::tautau {

enum class KinematicBranchCombination { average_conditional_observables, ratio_of_average_components };

struct PionPairHypothesisComponents {
    PionPairKinematicHypothesis hypothesis;
    LinearComponents components;
    PolynomialComponents polynomial;
    double weight;
};

class PionPairVisibleEventEvaluation {
public:
    PionPairVisibleEventEvaluation(HadronicTauPairSolutionStatus status, std::vector<PionPairHypothesisComponents> entries)
        : status_(status), entries_(std::move(entries)) {
        const std::size_t expected_entries =
            status_ == HadronicTauPairSolutionStatus::one_solution ? 1U
            : status_ == HadronicTauPairSolutionStatus::two_solutions ? 2U
                                                                        : 0U;
        if (entries_.size() != expected_entries)
            throw std::invalid_argument("visible pion-pair evaluation has incompatible solution entries");
        double total_weight = 0.0;
        for (const auto& entry : entries_) {
            if (!std::isfinite(entry.weight) || entry.weight < 0.0)
                throw std::invalid_argument("visible pion-pair evaluation weights must be finite and non-negative");
            total_weight += entry.weight;
        }
        if (!entries_.empty() && (!std::isfinite(total_weight) || total_weight <= 0.0))
            throw std::invalid_argument("visible pion-pair evaluation total weight must be finite and positive");
    }

    HadronicTauPairSolutionStatus status() const { return status_; }
    const std::vector<PionPairHypothesisComponents>& entries() const { return entries_; }

    PolynomialComponents marginalized_polynomial_components() const {
        std::vector<WeightedPolynomialComponents> components;
        components.reserve(entries_.size());
        for (const auto& entry : entries_) components.push_back({entry.polynomial, entry.weight});
        return average_polynomial_components(components);
    }

    LinearObservableComponents observable_components(KinematicBranchCombination combination) const {
        std::vector<WeightedLinearComponents> components;
        components.reserve(entries_.size());
        for (const auto& entry : entries_) components.push_back({entry.components, entry.weight});
        if (combination == KinematicBranchCombination::average_conditional_observables)
            return average_conditional_observables(components);
        return ratio_of_average_components(components);
    }

private:
    HadronicTauPairSolutionStatus status_;
    std::vector<PionPairHypothesisComponents> entries_;
};

class PionPairVisibleEventEvaluator {
public:
    explicit PionPairVisibleEventEvaluator(ElectroweakParameters parameters,
                                           ProductionBosons boson_selection = ProductionBosons::photon_and_z)
        : matrix_element_(parameters, boson_selection) {}

    PionPairVisibleEventEvaluation evaluate(const BeamState& beams, const FourMomentum& pion_minus_lab,
                                            const FourMomentum& pion_plus_lab, double tau_mass) const {
        const HadronicTauPairSolutionSet solution_set =
            HadronicTauPairKinematicSolver::solve({beams, pion_minus_lab, pion_plus_lab, tau_mass});
        std::vector<PionPairHypothesisComponents> entries;
        const double weight = solution_set.solutions.empty() ? 0.0 : 1.0 / solution_set.solutions.size();
        entries.reserve(solution_set.solutions.size());
        for (const auto& solution : solution_set.solutions) {
            const PionPairKinematicHypothesis hypothesis{solution.point, pion_minus_lab, pion_plus_lab,
                                                          KinematicPointOrigin::analytic_branch};
            entries.push_back({hypothesis, matrix_element_.components(hypothesis),
                               matrix_element_.polynomial_components(hypothesis), weight});
        }
        return {solution_set.status, std::move(entries)};
    }

private:
    PionPairMatrixElement matrix_element_;
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_VISIBLE_PION_PAIR_EVALUATION_H_
