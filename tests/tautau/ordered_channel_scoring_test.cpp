#include <cassert>
#include <cmath>
#include <vector>

#include "tauamp/tautau/ordered_channel_scoring.h"

namespace {
using tauamp::tautau::FourMomentum;
using tauamp::tautau::ReconstructedPion;

ReconstructedPion object(std::size_t index, int pid, double pt, double phi) {
    const double mass = pid == 111 ? 0.1349766 : 0.13957018;
    return {index, pid, pt, 0.0, phi,
            FourMomentum(pt * std::cos(phi), pt * std::sin(phi), 0.0, std::sqrt(pt * pt + mass * mass))};
}
}

int main() {
    using namespace tauamp::tautau;
    constexpr double pi = 3.14159265358979323846;
    const std::vector<ReconstructedPion> charged{object(1, 211, 0.4, 0.0), object(2, -211, 0.4, 0.0)};
    const std::vector<ReconstructedPion> neutral{object(10, 111, 0.4, pi)};
    const auto pi_pi = enumerate_ordered_channel_hypotheses(OrderedChannel::pi_pi, charged, neutral);
    const auto pi_rho = enumerate_ordered_channel_hypotheses(OrderedChannel::pi_rho, charged, neutral);
    assert(pi_pi.size() == 1U);
    assert(pi_rho.size() == 1U);
    const auto evidence = ordered_candidate_resonant_extension_evidence(pi_pi.front(), charged, neutral);
    assert(evidence > 0.0 && evidence <= 1.0);

    OrderedSelectionScoreConfig configuration;
    configuration.lambda_unused_neutral = 0.05;
    configuration.lambda_resonant_extension = 0.10;
    const auto pi_pi_score = score_ordered_candidate(pi_pi.front(), evidence, configuration);
    const auto pi_rho_score = score_ordered_candidate(pi_rho.front(), 0.0, configuration);
    assert(pi_pi_score.generic_unused_cost > 0.0);
    assert(pi_pi_score.resonant_extension_cost > 0.0);
    assert(pi_rho_score.generic_unused_cost == 0.0);

    const auto decision = select_ordered_event_winner({
        {OrderedChannel::pi_pi, 0U, pi_pi_score},
        {OrderedChannel::pi_rho, 0U, pi_rho_score}});
    assert(!decision.discarded);
    assert(decision.has_runner_up);
    assert(decision.channel == OrderedChannel::pi_rho);
    assert(decision.margin() > 0.0);

    const auto global_decision = select_ordered_channel_candidate(charged, neutral, configuration);
    assert(!global_decision.discarded);
    assert(global_decision.channel == OrderedChannel::pi_rho);
}
