#ifndef TAUAMP_TAUTAU_ORDERED_CHANNEL_SCORING_H_
#define TAUAMP_TAUTAU_ORDERED_CHANNEL_SCORING_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tauamp/tautau/ordered_channel_adapter.h"

namespace tauamp::tautau {

struct OrderedSelectionScoreConfig {
    double lambda_unused_charged{0.0};
    double lambda_unused_neutral{0.0};
    double lambda_resonant_extension{0.0};
    double lambda_topology_complexity{0.0};
    double charged_unused_weight{1.0};
    double neutral_unused_weight{1.0};
};

struct OrderedCandidateScore {
    double base_cost{0.0};
    double generic_unused_cost{0.0};
    double resonant_extension_evidence{0.0};
    double resonant_extension_cost{0.0};
    double topology_complexity_cost{0.0};
    double total_cost{0.0};
};

struct OrderedDecisionCandidate {
    OrderedChannel channel{};
    std::size_t candidate_index{0};
    OrderedCandidateScore score;
};

struct OrderedEventDecision {
    bool discarded{true};
    OrderedChannel channel{OrderedChannel::pi_pi};
    std::size_t candidate_index{0};
    double winner_cost{std::numeric_limits<double>::infinity()};
    bool has_runner_up{false};
    OrderedChannel runner_up_channel{OrderedChannel::pi_pi};
    std::size_t runner_up_candidate_index{0};
    double runner_up_cost{std::numeric_limits<double>::infinity()};

    double margin() const {
        if (!has_runner_up || discarded) return std::numeric_limits<double>::infinity();
        return runner_up_cost - winner_cost;
    }
};

inline constexpr std::array<OrderedChannel, 9> g9_ordered_channels{
    OrderedChannel::pi_pi,  OrderedChannel::pi_rho, OrderedChannel::rho_pi,
    OrderedChannel::rho_rho, OrderedChannel::pi_a1, OrderedChannel::a1_pi,
    OrderedChannel::rho_a1, OrderedChannel::a1_rho, OrderedChannel::a1_a1};

inline double ordered_candidate_resonant_extension_evidence(
    const OrderedChannelHypothesis& hypothesis, const std::vector<ReconstructedPion>& charged,
    const std::vector<ReconstructedPion>& neutral) {
    std::vector<std::size_t> used_charged;
    used_charged.insert(used_charged.end(), hypothesis.tau_plus.charged_indices.begin(),
                        hypothesis.tau_plus.charged_indices.end());
    used_charged.insert(used_charged.end(), hypothesis.tau_minus.charged_indices.begin(),
                        hypothesis.tau_minus.charged_indices.end());
    std::vector<std::size_t> used_neutral;
    used_neutral.insert(used_neutral.end(), hypothesis.tau_plus.neutral_indices.begin(),
                        hypothesis.tau_plus.neutral_indices.end());
    used_neutral.insert(used_neutral.end(), hypothesis.tau_minus.neutral_indices.begin(),
                        hypothesis.tau_minus.neutral_indices.end());
    const auto is_used = [](const std::vector<std::size_t>& used, std::size_t index) {
        return std::find(used.begin(), used.end(), index) != used.end();
    };
    const auto find_charged = [&](std::size_t index) -> const ReconstructedPion& {
        const auto found = std::find_if(charged.begin(), charged.end(),
                                        [index](const ReconstructedPion& object) { return object.index == index; });
        if (found == charged.end()) throw std::invalid_argument("extension references missing charged object");
        return *found;
    };
    double best = 0.0;
    for (const auto& pi0 : neutral) {
        if (is_used(used_neutral, pi0.index)) continue;
        for (const auto index : used_charged) {
            const auto& pion = find_charged(index);
            const double mass = detail::pair_mass(pion, pi0);
            if (detail::rho_pass(mass)) best = std::max(best, std::exp(-0.5 * detail::rho_mass_cost(mass)));
        }
    }
    const auto a1_extension = [&](const OrderedDecaySide& side, bool tau_plus) {
        if (side.mode != "pi" || side.charged_indices.size() != 1U) return 0.0;
        const auto& direct = find_charged(side.charged_indices.front());
        const int required_same_charge = tau_plus ? 1 : -1;
        const int required_opposite_charge = tau_plus ? -1 : 1;
        double evidence = 0.0;
        for (const auto& first : charged) {
            if (is_used(used_charged, first.index) || first.charge() != required_same_charge) continue;
            for (const auto& second : charged) {
                if (is_used(used_charged, second.index) || second.index == first.index ||
                    second.charge() != required_opposite_charge)
                    continue;
                const double mass = std::sqrt(std::max(
                    (direct.momentum + first.momentum + second.momentum).mass_squared(), 0.0));
                if (detail::a1_pass(mass)) evidence = std::max(evidence, std::exp(-0.5 * detail::a1_mass_cost(mass)));
            }
        }
        return evidence;
    };
    best = std::max(best, a1_extension(hypothesis.tau_plus, true));
    best = std::max(best, a1_extension(hypothesis.tau_minus, false));
    return best;
}

namespace detail {

inline bool score_is_better(const OrderedDecisionCandidate& left,
                            const OrderedDecisionCandidate& right) {
    constexpr double tolerance = 1.0e-12;
    if (left.score.total_cost + tolerance < right.score.total_cost) return true;
    if (std::abs(left.score.total_cost - right.score.total_cost) > tolerance) return false;
    const auto left_name = ordered_channel_name(left.channel);
    const auto right_name = ordered_channel_name(right.channel);
    if (std::string(left_name) != std::string(right_name))
        return std::string(left_name) < std::string(right_name);
    return left.candidate_index < right.candidate_index;
}

}  // namespace detail

inline OrderedCandidateScore score_ordered_candidate(
    const OrderedChannelHypothesis& hypothesis, double resonant_extension_evidence,
    const OrderedSelectionScoreConfig& configuration) {
    if (!std::isfinite(hypothesis.normalized_mass_cost) || hypothesis.normalized_mass_cost < 0.0)
        throw std::invalid_argument("ordered candidate has invalid base cost");
    if (!std::isfinite(resonant_extension_evidence) || resonant_extension_evidence < 0.0)
        throw std::invalid_argument("ordered candidate has invalid extension evidence");
    const double generic_unused_cost =
        configuration.lambda_unused_charged * configuration.charged_unused_weight *
            static_cast<double>(hypothesis.unused_charged_count) +
        configuration.lambda_unused_neutral * configuration.neutral_unused_weight *
            static_cast<double>(hypothesis.unused_neutral_count);
    const double extension_cost = configuration.lambda_resonant_extension * resonant_extension_evidence;
    const double topology_complexity_cost = configuration.lambda_topology_complexity * static_cast<double>(
        hypothesis.rho_masses.size() + hypothesis.a1_masses.size());
    return {hypothesis.normalized_mass_cost, generic_unused_cost, resonant_extension_evidence,
            extension_cost, topology_complexity_cost,
            hypothesis.normalized_mass_cost + generic_unused_cost + extension_cost + topology_complexity_cost};
}

inline OrderedEventDecision select_ordered_event_winner(
    std::vector<OrderedDecisionCandidate> candidates) {
    OrderedEventDecision result;
    if (candidates.empty()) return result;
    std::sort(candidates.begin(), candidates.end(), detail::score_is_better);
    const auto& winner = candidates.front();
    result.discarded = false;
    result.channel = winner.channel;
    result.candidate_index = winner.candidate_index;
    result.winner_cost = winner.score.total_cost;
    if (candidates.size() > 1U) {
        const auto& runner_up = candidates[1];
        result.has_runner_up = true;
        result.runner_up_channel = runner_up.channel;
        result.runner_up_candidate_index = runner_up.candidate_index;
        result.runner_up_cost = runner_up.score.total_cost;
    }
    return result;
}

inline OrderedEventDecision select_ordered_channel_candidate(
    const std::vector<ReconstructedPion>& charged, const std::vector<ReconstructedPion>& neutral,
    const OrderedSelectionScoreConfig& configuration) {
    std::vector<OrderedDecisionCandidate> candidates;
    for (const auto channel : g9_ordered_channels) {
        const auto hypotheses = enumerate_ordered_channel_hypotheses(channel, charged, neutral);
        for (const auto& hypothesis : hypotheses) {
            const double evidence = ordered_candidate_resonant_extension_evidence(hypothesis, charged, neutral);
            candidates.push_back({channel, hypothesis.hypothesis_index,
                                  score_ordered_candidate(hypothesis, evidence, configuration)});
        }
    }
    return select_ordered_event_winner(std::move(candidates));
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_ORDERED_CHANNEL_SCORING_H_
