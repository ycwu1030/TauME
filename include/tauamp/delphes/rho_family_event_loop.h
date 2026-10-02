#ifndef TAUAMP_DELPHES_RHO_FAMILY_EVENT_LOOP_H_
#define TAUAMP_DELPHES_RHO_FAMILY_EVENT_LOOP_H_

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "tauamp/delphes/common_channel_delphes_reader.h"
#include "tauamp/delphes/rho_family_cross_feed.h"
#include "tauamp/tautau/ordered_channel_evaluation.h"

namespace tauamp::delphes {

struct RhoFamilyEventEvaluation {
    std::array<std::optional<tautau::OrderedEventEvaluation>, g7_rho_channels.size()> targets;
    std::vector<RhoTargetDiagnostic> diagnostics;
};

inline RhoFamilyEventEvaluation evaluate_rho_family_event(
    const CommonChannelPionCollections& collections, const tautau::BeamState& beams, double tau_mass,
    const tautau::HadronicPairMatrixElement& matrix_element) {
    RhoFamilyEventEvaluation result;
    result.diagnostics.reserve(g7_rho_channels.size());
    for (std::size_t index = 0; index < g7_rho_channels.size(); ++index) {
        try {
            auto evaluation = tautau::evaluate_ordered_event(
                g7_rho_channels[index], collections.charged, collections.neutral, beams, tau_mass, matrix_element);
            result.diagnostics.push_back(
                {g7_rho_channels[index], evaluation.candidate_count, evaluation.valid_candidate_count, false});
            result.targets[index] = std::move(evaluation);
        } catch (...) {
            result.diagnostics.push_back({g7_rho_channels[index], 0U, 0U, true});
        }
    }
    return result;
}

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_RHO_FAMILY_EVENT_LOOP_H_
