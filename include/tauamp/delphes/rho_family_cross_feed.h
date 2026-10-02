#ifndef TAUAMP_DELPHES_RHO_FAMILY_CROSS_FEED_H_
#define TAUAMP_DELPHES_RHO_FAMILY_CROSS_FEED_H_

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "tauamp/tautau/ordered_channel_adapter.h"

namespace tauamp::delphes {

inline constexpr std::array<tautau::OrderedChannel, 3> g7_rho_channels{
    tautau::OrderedChannel::pi_rho, tautau::OrderedChannel::rho_pi, tautau::OrderedChannel::rho_rho};

inline std::size_t g7_rho_channel_index(tautau::OrderedChannel channel) {
    for (std::size_t index = 0; index < g7_rho_channels.size(); ++index)
        if (g7_rho_channels[index] == channel) return index;
    throw std::invalid_argument("channel is not in the G7 rho family");
}

struct RhoTargetDiagnostic {
    tautau::OrderedChannel channel{};
    std::size_t candidate_count{0};
    std::size_t valid_candidate_count{0};
    bool technical_failure{false};
};

struct RhoCrossFeedDecision {
    std::string category;
    tautau::OrderedChannel target{};
    bool has_target{false};
};

class RhoFamilyCrossFeedCounter {
public:
    explicit RhoFamilyCrossFeedCounter(tautau::OrderedChannel source) : source_(source) {
        (void)g7_rho_channel_index(source_);
    }

    RhoCrossFeedDecision record(long long event_id, const std::vector<RhoTargetDiagnostic>& diagnostics) {
        if (event_id < 0) throw std::invalid_argument("G7 rho source event IDs must be non-negative");
        if (event_ids_.count(event_id) != 0U) throw std::invalid_argument("duplicate G7 rho source event ID");
        event_ids_[event_id] = {};

        std::vector<tautau::OrderedChannel> valid_targets;
        bool has_candidate = false;
        bool technical = false;
        for (const auto& diagnostic : diagnostics) {
            (void)g7_rho_channel_index(diagnostic.channel);
            has_candidate = has_candidate || diagnostic.candidate_count != 0U;
            technical = technical || diagnostic.technical_failure;
            if (diagnostic.valid_candidate_count != 0U) valid_targets.push_back(diagnostic.channel);
        }

        RhoCrossFeedDecision decision{};
        if (technical) {
            decision.category = "technical_failure";
        } else if (valid_targets.size() > 1U) {
            decision.category = "ambiguous_multi_hypothesis";
        } else if (valid_targets.size() == 1U) {
            decision.target = valid_targets.front();
            decision.has_target = true;
            const auto source_name = std::string(tautau::ordered_channel_name(source_));
            const auto target_name = std::string(tautau::ordered_channel_name(decision.target));
            if (decision.target == source_)
                decision.category = "selected::" + target_name;
            else
                decision.category = "cross_feed::" + source_name + "::" + target_name;
            ++matrix_[g7_rho_channel_index(source_)][g7_rho_channel_index(decision.target)];
        } else if (has_candidate) {
            decision.category = "unsupported_or_no_solution";
        } else {
            decision.category = "rejected::no_ownership_candidate";
        }
        event_ids_[event_id] = decision.category;
        ++category_counts_[decision.category];
        return decision;
    }

    const std::array<std::array<long long, 3>, 3>& matrix() const { return matrix_; }
    const std::map<std::string, long long>& category_counts() const { return category_counts_; }
    const std::map<long long, std::string>& event_categories() const { return event_ids_; }

private:
    tautau::OrderedChannel source_{};
    std::array<std::array<long long, 3>, 3> matrix_{};
    std::map<std::string, long long> category_counts_;
    std::map<long long, std::string> event_ids_;
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_RHO_FAMILY_CROSS_FEED_H_
