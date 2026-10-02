#ifndef TAUAMP_TAUTAU_ORDERED_CHANNEL_EXPLANATION_SCORING_H_
#define TAUAMP_TAUTAU_ORDERED_CHANNEL_EXPLANATION_SCORING_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/ordered_channel_scoring.h"

namespace tauamp::tautau {

struct OrderedExplanationScoreConfig {
    double lambda_unexplained_object{0.01};
    double lambda_unexplained_resonance{0.10};
    double lambda_used_resonance_reward{0.10};
    double lambda_topology_complexity{0.01};
};

struct OrderedExplanationScore {
    double base_mass_cost{0.0};
    double used_resonance_evidence{0.0};
    double unexplained_object_cost{0.0};
    double unexplained_resonance_evidence{0.0};
    double used_resonance_reward{0.0};
    double topology_complexity_cost{0.0};
    double total_cost{0.0};
};

namespace explanation_detail {

inline bool contains_index(const std::vector<std::size_t>& indices, std::size_t index) {
    return std::find(indices.begin(), indices.end(), index) != indices.end();
}

inline std::vector<std::size_t> used_charged_indices(const OrderedChannelHypothesis& hypothesis) {
    std::vector<std::size_t> result = hypothesis.tau_plus.charged_indices;
    result.insert(result.end(), hypothesis.tau_minus.charged_indices.begin(), hypothesis.tau_minus.charged_indices.end());
    return result;
}

inline std::vector<std::size_t> used_neutral_indices(const OrderedChannelHypothesis& hypothesis) {
    std::vector<std::size_t> result = hypothesis.tau_plus.neutral_indices;
    result.insert(result.end(), hypothesis.tau_minus.neutral_indices.begin(), hypothesis.tau_minus.neutral_indices.end());
    return result;
}

inline double best_neutral_rho_opportunity(const ReconstructedPion& neutral,
                                           const std::vector<std::size_t>& used_charged_indices,
                                           const std::vector<ReconstructedPion>& charged) {
    double best = 0.0;
    for (const auto index : used_charged_indices) {
        const auto found = std::find_if(charged.begin(), charged.end(),
                                        [index](const ReconstructedPion& object) { return object.index == index; });
        if (found == charged.end()) throw std::invalid_argument("missing used charged object");
        const double mass = detail::pair_mass(*found, neutral);
        if (detail::rho_pass(mass)) best = std::max(best, std::exp(-0.5 * detail::rho_mass_cost(mass)));
    }
    return best;
}

inline double best_charged_a1_opportunity(const ReconstructedPion& object,
                                          const std::vector<ReconstructedPion>& charged) {
    double best = 0.0;
    for (std::size_t first = 0; first < charged.size(); ++first) {
        for (std::size_t second = first + 1U; second < charged.size(); ++second) {
            if (charged[first].index == object.index || charged[second].index == object.index) continue;
            const int charge_sum = object.charge() + charged[first].charge() + charged[second].charge();
            if (std::abs(charge_sum) != 1) continue;
            const double mass = std::sqrt(std::max(
                (object.momentum + charged[first].momentum + charged[second].momentum).mass_squared(), 0.0));
            if (detail::a1_pass(mass)) best = std::max(best, std::exp(-0.5 * detail::a1_mass_cost(mass)));
        }
    }
    return best;
}

inline double used_resonance_evidence(const OrderedChannelHypothesis& hypothesis) {
    double evidence = 0.0;
    for (const double mass : hypothesis.rho_masses)
        evidence += std::exp(-0.5 * detail::rho_mass_cost(mass));
    for (const double mass : hypothesis.a1_masses)
        evidence += std::exp(-0.5 * detail::a1_mass_cost(mass));
    return evidence;
}

}  // namespace explanation_detail

inline OrderedExplanationScore score_ordered_candidate_explanation(
    const OrderedChannelHypothesis& hypothesis, const std::vector<ReconstructedPion>& charged,
    const std::vector<ReconstructedPion>& neutral, const OrderedExplanationScoreConfig& configuration) {
    if (!std::isfinite(hypothesis.normalized_mass_cost) || hypothesis.normalized_mass_cost < 0.0)
        throw std::invalid_argument("ordered explanation candidate has invalid mass cost");
    const auto used_charged = explanation_detail::used_charged_indices(hypothesis);
    const auto used_neutral = explanation_detail::used_neutral_indices(hypothesis);
    double unexplained_object_cost = 0.0;
    double unexplained_resonance_evidence = 0.0;
    for (const auto& object : charged) {
        if (explanation_detail::contains_index(used_charged, object.index)) continue;
        const double opportunity = explanation_detail::best_charged_a1_opportunity(object, charged);
        unexplained_resonance_evidence += opportunity;
        unexplained_object_cost += configuration.lambda_unexplained_object *
                                   (1.0 + configuration.lambda_unexplained_resonance * opportunity);
    }
    for (const auto& object : neutral) {
        if (explanation_detail::contains_index(used_neutral, object.index)) continue;
        const double opportunity = explanation_detail::best_neutral_rho_opportunity(object, used_charged, charged);
        unexplained_resonance_evidence += opportunity;
        unexplained_object_cost += configuration.lambda_unexplained_object *
                                   (1.0 + configuration.lambda_unexplained_resonance * opportunity);
    }
    const double used_evidence = explanation_detail::used_resonance_evidence(hypothesis);
    const double used_reward = -configuration.lambda_used_resonance_reward * used_evidence;
    const double complexity = configuration.lambda_topology_complexity *
                              std::log1p(static_cast<double>(hypothesis.rho_masses.size() + hypothesis.a1_masses.size()));
    return {hypothesis.normalized_mass_cost, used_evidence, unexplained_object_cost,
            unexplained_resonance_evidence, used_reward, complexity,
            hypothesis.normalized_mass_cost + unexplained_object_cost + used_reward + complexity};
}

inline OrderedEventDecision select_ordered_channel_explanation_candidate(
    const std::vector<ReconstructedPion>& charged, const std::vector<ReconstructedPion>& neutral,
    const OrderedExplanationScoreConfig& configuration) {
    std::vector<OrderedDecisionCandidate> candidates;
    for (const auto channel : g9_ordered_channels) {
        const auto hypotheses = enumerate_ordered_channel_hypotheses(channel, charged, neutral);
        for (const auto& hypothesis : hypotheses) {
            const auto explanation = score_ordered_candidate_explanation(hypothesis, charged, neutral, configuration);
            OrderedCandidateScore score;
            score.base_cost = explanation.base_mass_cost;
            score.generic_unused_cost = explanation.unexplained_object_cost;
            score.resonant_extension_evidence = explanation.unexplained_resonance_evidence;
            score.resonant_extension_cost = explanation.used_resonance_reward;
            score.topology_complexity_cost = explanation.topology_complexity_cost;
            score.total_cost = explanation.total_cost;
            candidates.push_back({channel, hypothesis.hypothesis_index, score});
        }
    }
    return select_ordered_event_winner(std::move(candidates));
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_ORDERED_CHANNEL_EXPLANATION_SCORING_H_
