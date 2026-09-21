#ifndef TAUAMP_DELPHES_COMMON_CHANNEL_EVENT_ADAPTER_H_
#define TAUAMP_DELPHES_COMMON_CHANNEL_EVENT_ADAPTER_H_

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tauamp/delphes/common_channel_delphes_reader.h"
#include "tauamp/tautau/common_channel_export.h"

namespace tauamp::delphes {

using CommonChannelContraction =
    std::function<tautau::CommonHypothesisRecord(const tautau::OrderedChannelHypothesis&)>;

struct CommonChannelEventBuildResult {
    bool selected{false};
    tautau::CommonEventRecord record{};
    std::vector<std::string> diagnostics;
};

inline std::vector<std::string> classify_common_channel_rejection(
    tautau::OrderedChannel channel, const CommonChannelPionCollections& collections) {
    const auto candidates = tautau::enumerate_ordered_channel_hypotheses(
        channel, collections.charged, collections.neutral);
    if (!candidates.empty()) return {};
    std::vector<std::string> diagnostics{"rejected::no_ownership_candidate"};
    const bool rho_channel = channel == tautau::OrderedChannel::pi_rho ||
                             channel == tautau::OrderedChannel::rho_pi ||
                             channel == tautau::OrderedChannel::rho_rho ||
                             channel == tautau::OrderedChannel::rho_a1 ||
                             channel == tautau::OrderedChannel::a1_rho;
    if (rho_channel) diagnostics.push_back("rejected::rho_mass_or_object_requirement");
    if (rho_channel && collections.neutral.size() > 2U)
        diagnostics.push_back("diagnostic::extra_neutral_pions_hardest_required_count");
    return diagnostics;
}

inline CommonChannelEventBuildResult try_build_common_channel_event_record(
    long long source_event_id, tautau::OrderedChannel channel,
    const CommonChannelPionCollections& collections, const std::string& category_code,
    CommonChannelContraction contract) {
    if (category_code.empty()) throw std::invalid_argument("common-channel category code must be non-empty");
    if (!contract) throw std::invalid_argument("common-channel contraction callback must be provided");

    const auto candidates = tautau::enumerate_ordered_channel_hypotheses(
        channel, collections.charged, collections.neutral);
    if (candidates.empty()) {
        return {false, {}, classify_common_channel_rejection(channel, collections)};
    }

    std::vector<tautau::CommonHypothesisRecord> hypotheses;
    hypotheses.reserve(candidates.size());
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        tautau::CommonHypothesisRecord record = contract(candidates[index]);
        record.hypothesis_index = index;
        record.valid = true;
        hypotheses.push_back(record);
    }
    std::vector<std::string> diagnostics;
    if (candidates.size() > 1U)
        diagnostics.push_back("diagnostic::cross_feed_or_ownership_ambiguity");
    return {true,
            {source_event_id, channel, category_code, tautau::ordered_channel_name(channel), candidates.size(),
             std::move(hypotheses)},
            std::move(diagnostics)};
}

inline tautau::CommonEventRecord build_common_channel_event_record(
    long long source_event_id, tautau::OrderedChannel channel,
    const CommonChannelPionCollections& collections, const std::string& category_code,
    CommonChannelContraction contract) {
    if (category_code.empty()) throw std::invalid_argument("common-channel category code must be non-empty");
    if (!contract) throw std::invalid_argument("common-channel contraction callback must be provided");
    const auto result = try_build_common_channel_event_record(
        source_event_id, channel, collections, category_code, std::move(contract));
    if (!result.selected) throw std::invalid_argument("common-channel event has no selected ownership hypothesis");
    return result.record;
}

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_COMMON_CHANNEL_EVENT_ADAPTER_H_
