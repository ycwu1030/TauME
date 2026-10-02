#ifndef TAUAMP_DELPHES_RHO_FAMILY_EVENT_LOOP_V2_H_
#define TAUAMP_DELPHES_RHO_FAMILY_EVENT_LOOP_V2_H_

#include <array>
#include <cmath>
#include <optional>

#include "tauamp/delphes/rho_family_selection_v2.h"

namespace tauamp::delphes {

struct RhoFamilyEventEvaluationV2 {
    std::array<RhoV2TargetResult, 3> targets;
    std::optional<std::size_t> winner;
    bool pre_resolution_multi_target{false};
    bool tie_resolved_by_delta_r{false};
};

inline RhoFamilyEventEvaluationV2 evaluate_rho_family_event_v2(
    const CommonChannelPionCollections& collections, const tautau::BeamState& beams,
    double tau_mass, const tautau::HadronicPairMatrixElement& matrix_element,
    double width_multiplier) {
    RhoFamilyEventEvaluationV2 result;
    for (std::size_t index = 0; index < g7_rho_channels.size(); ++index)
        result.targets[index] = evaluate_g7_rho_target_v2(
            g7_rho_channels[index], collections, beams, tau_mass, matrix_element, width_multiplier);

    std::array<std::size_t, 3> valid{};
    std::size_t valid_count = 0U;
    for (std::size_t index = 0; index < result.targets.size(); ++index)
        if (result.targets[index].diagnostic.valid_candidate_count != 0U) valid[valid_count++] = index;
    result.pre_resolution_multi_target = valid_count > 1U;
    if (valid_count == 0U) return result;

    std::size_t winner = valid[0];
    for (std::size_t position = 1; position < valid_count; ++position) {
        const std::size_t candidate = valid[position];
        const double mass_difference = result.targets[candidate].diagnostic.mass_score - result.targets[winner].diagnostic.mass_score;
        if (mass_difference < -1.0e-12) {
            winner = candidate;
        } else if (std::abs(mass_difference) <= 1.0e-12) {
            result.tie_resolved_by_delta_r = true;
            const double delta_r_difference = result.targets[candidate].diagnostic.delta_r_score - result.targets[winner].diagnostic.delta_r_score;
            if (delta_r_difference < -1.0e-12) winner = candidate;
        }
    }
    result.winner = winner;
    return result;
}

}  // namespace tauamp::delphes

#endif
