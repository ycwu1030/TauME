#ifndef TAUAMP_TAUTAU_ORDERED_CHANNEL_EVALUATION_H_
#define TAUAMP_TAUTAU_ORDERED_CHANNEL_EVALUATION_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tauamp/tautau/common_channel_export.h"
#include "tauamp/tautau/hadronic_pair_matrix_element.h"
#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"

namespace tauamp::tautau {

struct OrderedDecayEvaluation {
    FourMomentum visible_lab;
    ComplexFourVector current_lab;
};

struct OrderedCandidateEvaluation {
    HadronicTauPairSolutionStatus status{HadronicTauPairSolutionStatus::no_solution};
    std::vector<CommonHypothesisRecord> branches;
};

struct OrderedEventEvaluation {
    std::size_t candidate_count{};
    std::size_t valid_candidate_count{};
    std::vector<HadronicTauPairSolutionStatus> candidate_statuses;
    std::vector<CommonHypothesisRecord> branches;
};

namespace detail {
inline const ReconstructedPion& object_by_index(const std::vector<ReconstructedPion>& objects,
                                                std::size_t index, const char* collection) {
    const auto found = std::find_if(objects.begin(), objects.end(),
                                    [index](const ReconstructedPion& object) { return object.index == index; });
    if (found == objects.end()) throw std::invalid_argument(std::string("missing object index in ") + collection);
    return *found;
}

inline std::array<double, g6_component_count> component_array(const PolynomialComponents& value) {
    return {{value.sm, value.f2_real, value.f2_imaginary, value.f3_real, value.f3_imaginary,
             value.f2_real_f2_real, value.f2_real_f2_imaginary, value.f2_real_f3_real,
             value.f2_real_f3_imaginary, value.f2_imaginary_f2_imaginary,
             value.f2_imaginary_f3_real, value.f2_imaginary_f3_imaginary,
             value.f3_real_f3_real, value.f3_real_f3_imaginary,
             value.f3_imaginary_f3_imaginary}};
}
}  // namespace detail

inline OrderedDecayEvaluation evaluate_ordered_decay_side(const OrderedDecaySide& side,
                                                           const std::vector<ReconstructedPion>& charged,
                                                           const std::vector<ReconstructedPion>& neutral,
                                                           bool tau_plus) {
    if (side.mode == "pi") {
        if (side.charged_indices.size() != 1U || !side.neutral_indices.empty())
            throw std::invalid_argument("pion side has invalid daughter content");
        const FourMomentum pion = detail::object_by_index(charged, side.charged_indices[0], "charged").momentum;
        return {pion, complex_scaled(pion, 1.0)};
    }
    if (side.mode == "rho") {
        if (side.charged_indices.size() != 1U || side.neutral_indices.size() != 1U)
            throw std::invalid_argument("rho side has invalid daughter content");
        const FourMomentum charged_pion =
            detail::object_by_index(charged, side.charged_indices[0], "charged").momentum;
        const FourMomentum neutral_pion =
            detail::object_by_index(neutral, side.neutral_indices[0], "neutral").momentum;
        return {charged_pion + neutral_pion, rho_kw_current(charged_pion, neutral_pion, tau_plus)};
    }
    if (side.mode == "a1") {
        if (side.same_charge_slots.size() != 2U || !side.has_opposite_charge_slot || !side.neutral_indices.empty())
            throw std::invalid_argument("a1 side has invalid daughter content");
        const FourMomentum same_charge_1 =
            detail::object_by_index(charged, side.same_charge_slots[0], "charged").momentum;
        const FourMomentum same_charge_2 =
            detail::object_by_index(charged, side.same_charge_slots[1], "charged").momentum;
        const FourMomentum opposite =
            detail::object_by_index(charged, side.opposite_charge_slot, "charged").momentum;
        return {same_charge_1 + same_charge_2 + opposite,
                a1_kw_current(same_charge_1, same_charge_2, opposite, tau_plus)};
    }
    throw std::invalid_argument("unknown ordered decay-side mode");
}

inline OrderedCandidateEvaluation evaluate_ordered_candidate(
    const OrderedChannelHypothesis& candidate, const std::vector<ReconstructedPion>& charged,
    const std::vector<ReconstructedPion>& neutral, const BeamState& beams, double tau_mass,
    const HadronicPairMatrixElement& matrix_element) {
    const OrderedDecayEvaluation tau_plus =
        evaluate_ordered_decay_side(candidate.tau_plus, charged, neutral, true);
    const OrderedDecayEvaluation tau_minus =
        evaluate_ordered_decay_side(candidate.tau_minus, charged, neutral, false);
    const HadronicTauPairSolutionSet solutions = HadronicTauPairKinematicSolver::solve(
        {beams, tau_minus.visible_lab, tau_plus.visible_lab, tau_mass});
    OrderedCandidateEvaluation result;
    result.status = solutions.status;
    if (solutions.solutions.empty()) return result;
    const double weight = 1.0 / static_cast<double>(solutions.solutions.size());
    for (std::size_t index = 0; index < solutions.solutions.size(); ++index) {
        const auto& solution = solutions.solutions[index];
        const auto minus_analyser = hadronic_analyser_tnl(
            solution.point, solution.neutrino_minus_lab, tau_minus.current_lab, false);
        const auto plus_analyser = hadronic_analyser_tnl(
            solution.point, solution.neutrino_plus_lab, tau_plus.current_lab, true);
        const PolynomialComponents polynomial =
            matrix_element.polynomial_components(solution.point, minus_analyser, plus_analyser);
        CommonHypothesisRecord branch;
        branch.hypothesis_index = index;
        branch.weight = weight;
        branch.valid = true;
        branch.components = detail::component_array(polynomial);
        result.branches.push_back(branch);
    }
    return result;
}

inline OrderedEventEvaluation evaluate_ordered_event(
    OrderedChannel channel, const std::vector<ReconstructedPion>& charged,
    const std::vector<ReconstructedPion>& neutral, const BeamState& beams, double tau_mass,
    const HadronicPairMatrixElement& matrix_element) {
    const auto candidates = enumerate_ordered_channel_hypotheses(channel, charged, neutral);
    OrderedEventEvaluation result;
    result.candidate_count = candidates.size();
    std::vector<OrderedCandidateEvaluation> valid;
    for (const auto& candidate : candidates) {
        OrderedCandidateEvaluation evaluation =
            evaluate_ordered_candidate(candidate, charged, neutral, beams, tau_mass, matrix_element);
        result.candidate_statuses.push_back(evaluation.status);
        if (!evaluation.branches.empty()) valid.push_back(std::move(evaluation));
    }
    result.valid_candidate_count = valid.size();
    if (valid.empty()) return result;
    const double candidate_weight = 1.0 / static_cast<double>(valid.size());
    std::size_t output_index = 0U;
    for (auto& candidate : valid) {
        for (auto& branch : candidate.branches) {
            branch.hypothesis_index = output_index++;
            branch.weight *= candidate_weight;
            result.branches.push_back(branch);
        }
    }
    return result;
}

inline CommonEventRecord make_ordered_common_event_record(long long source_event_id, OrderedChannel channel,
                                                          const OrderedEventEvaluation& evaluation) {
    if (evaluation.branches.empty())
        throw std::invalid_argument("ordered event has no valid physical TauOO branch");
    return {source_event_id, channel, std::string("selected::") + ordered_channel_name(channel),
            ordered_channel_name(channel), evaluation.candidate_count, evaluation.branches};
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_ORDERED_CHANNEL_EVALUATION_H_
