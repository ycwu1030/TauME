#ifndef TAUAMP_DELPHES_RHO_FAMILY_SELECTION_V2_H_
#define TAUAMP_DELPHES_RHO_FAMILY_SELECTION_V2_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include "tauamp/delphes/common_channel_delphes_reader.h"
#include "tauamp/delphes/rho_family_cross_feed.h"
#include "tauamp/tautau/ordered_channel_evaluation.h"

namespace tauamp::delphes {

struct RhoV2TargetDiagnostic {
    tautau::OrderedChannel channel{};
    std::size_t candidate_count{0};
    std::size_t valid_candidate_count{0};
    bool technical_failure{false};
    double mass_score{0.0};
    double delta_r_score{0.0};
};

struct RhoV2TargetResult {
    RhoV2TargetDiagnostic diagnostic;
    std::optional<tautau::OrderedEventEvaluation> evaluation;
};

inline std::vector<tautau::ReconstructedPion> g7_v2_hardest_charged(
    const std::vector<tautau::ReconstructedPion>& charged, int pid) {
    std::vector<tautau::ReconstructedPion> result;
    for (const auto& object : charged)
        if (object.pid == pid) result.push_back(object);
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        if (left.pt != right.pt) return left.pt > right.pt;
        return left.index < right.index;
    });
    if (result.size() > 1U) result.resize(1U);
    return result;
}

inline std::vector<tautau::ReconstructedPion> g7_v2_all_charged(
    const std::vector<tautau::ReconstructedPion>& charged, int pid) {
    std::vector<tautau::ReconstructedPion> result;
    for (const auto& object : charged)
        if (object.pid == pid) result.push_back(object);
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        if (left.pt != right.pt) return left.pt > right.pt;
        return left.index < right.index;
    });
    return result;
}

inline RhoV2TargetResult evaluate_g7_rho_target_v2(
    tautau::OrderedChannel channel, const CommonChannelPionCollections& collections,
    const tautau::BeamState& beams, double tau_mass,
    const tautau::HadronicPairMatrixElement& matrix_element, double width_multiplier) {
    using namespace tautau;
    RhoV2TargetResult result;
    result.diagnostic.channel = channel;
    const auto plus = g7_v2_all_charged(collections.charged, 211);
    const auto minus = g7_v2_all_charged(collections.charged, -211);
    const auto neutral = detail::hardest_neutrals(collections.neutral, channel == OrderedChannel::rho_rho ? 2U : 1U);
    if (plus.empty() || minus.empty() || neutral.empty()) return result;

    OrderedChannelHypothesis hypothesis;
    hypothesis.channel = channel;
    hypothesis.hypothesis_index = 0U;
    const double window = width_multiplier * g6_rho_width;
    if (channel == OrderedChannel::pi_rho || channel == OrderedChannel::rho_pi) {
        const auto& rho_candidates = channel == OrderedChannel::pi_rho ? minus : plus;
        const auto& pion_side_charged = channel == OrderedChannel::pi_rho ? plus.front() : minus.front();
        std::size_t best_index = 0U;
        double best_mass = 0.0;
        double best_score = 1.0e300;
        double best_dr = 1.0e300;
        for (std::size_t index = 0; index < rho_candidates.size(); ++index) {
            const double mass = detail::pair_mass(rho_candidates[index], neutral.front());
            const double score = std::pow(mass - g6_rho_mass, 2);
            const double dr = detail::delta_r(rho_candidates[index], neutral.front());
            if (score < best_score || (std::abs(score - best_score) <= 1.0e-12 && dr < best_dr)) {
                best_index = index;
                best_mass = mass;
                best_score = score;
                best_dr = dr;
            }
        }
        if (std::abs(best_mass - g6_rho_mass) >= window) return result;
        result.diagnostic.candidate_count = 1U;
        result.diagnostic.mass_score = best_score;
        result.diagnostic.delta_r_score = best_dr;
        hypothesis.rho_masses = {best_mass};
        hypothesis.tau_plus = channel == OrderedChannel::pi_rho ? detail::pion_side(pion_side_charged) : detail::rho_side(rho_candidates[best_index], neutral.front());
        hypothesis.tau_minus = channel == OrderedChannel::pi_rho ? detail::rho_side(rho_candidates[best_index], neutral.front()) : detail::pion_side(pion_side_charged);
    } else if (channel == OrderedChannel::rho_rho) {
        struct Assignment { double score; double delta_r; std::size_t plus_index; std::size_t minus_index; double plus_mass; double minus_mass; };
        Assignment best{1.0e300, 1.0e300, 0U, 0U, 0.0, 0.0};
        std::size_t best_plus_neutral = 0U;
        std::size_t best_minus_neutral = 1U;
        for (std::size_t pi = 0; pi < plus.size(); ++pi) for (std::size_t mi = 0; mi < minus.size(); ++mi) {
            for (const auto ordering : {std::pair<std::size_t, std::size_t>{0U, 1U}, {1U, 0U}}) {
                const double plus_mass = detail::pair_mass(plus[pi], neutral[ordering.first]);
                const double minus_mass = detail::pair_mass(minus[mi], neutral[ordering.second]);
                const Assignment candidate{
                    std::pow(plus_mass - g6_rho_mass, 2) + std::pow(minus_mass - g6_rho_mass, 2),
                    detail::delta_r(plus[pi], neutral[ordering.first]) + detail::delta_r(minus[mi], neutral[ordering.second]),
                    pi, mi, plus_mass, minus_mass};
                if (candidate.score < best.score || (std::abs(candidate.score - best.score) <= 1.0e-4 && candidate.delta_r < best.delta_r)) {
                    best = candidate;
                    best_plus_neutral = ordering.first;
                    best_minus_neutral = ordering.second;
                }
            }
        }
        if (std::abs(best.plus_mass - g6_rho_mass) >= window || std::abs(best.minus_mass - g6_rho_mass) >= window) return result;
        result.diagnostic.candidate_count = 1U;
        result.diagnostic.mass_score = 0.5 * best.score;
        result.diagnostic.delta_r_score = 0.5 * best.delta_r;
        hypothesis.rho_masses = {best.plus_mass, best.minus_mass};
        hypothesis.tau_plus = detail::rho_side(plus[best.plus_index], neutral[best_plus_neutral]);
        hypothesis.tau_minus = detail::rho_side(minus[best.minus_index], neutral[best_minus_neutral]);
    } else {
        return result;
    }

    try {
        const auto evaluated = evaluate_ordered_candidate(
            hypothesis, collections.charged, collections.neutral, beams, tau_mass, matrix_element);
        if (!evaluated.branches.empty()) {
            result.diagnostic.valid_candidate_count = 1U;
            result.evaluation = OrderedEventEvaluation{1U, 1U, {evaluated.status}, evaluated.branches};
        }
    } catch (...) {
        result.diagnostic.technical_failure = true;
    }
    return result;
}

}  // namespace tauamp::delphes

#endif
