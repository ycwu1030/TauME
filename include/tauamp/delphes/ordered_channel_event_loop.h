#ifndef TAUAMP_DELPHES_ORDERED_CHANNEL_EVENT_LOOP_H_
#define TAUAMP_DELPHES_ORDERED_CHANNEL_EVENT_LOOP_H_

#include <array>
#include <cstddef>

#include "tauamp/delphes/common_channel_delphes_reader.h"
#include "tauamp/tautau/ordered_channel_scoring.h"

namespace tauamp::delphes {

struct OrderedChannelCandidateAssignment {
    tautau::OrderedEventDecision decision;
    std::array<std::size_t, tautau::g9_ordered_channels.size()> candidate_counts{};
};

inline OrderedChannelCandidateAssignment evaluate_ordered_channel_candidate_assignment(
    const CommonChannelPionCollections& collections,
    const tautau::OrderedSelectionScoreConfig& configuration) {
    OrderedChannelCandidateAssignment result;
    for (std::size_t index = 0; index < tautau::g9_ordered_channels.size(); ++index)
        result.candidate_counts[index] = tautau::enumerate_ordered_channel_hypotheses(
            tautau::g9_ordered_channels[index], collections.charged, collections.neutral).size();
    result.decision = tautau::select_ordered_channel_candidate(
        collections.charged, collections.neutral, configuration);
    return result;
}

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_ORDERED_CHANNEL_EVENT_LOOP_H_
