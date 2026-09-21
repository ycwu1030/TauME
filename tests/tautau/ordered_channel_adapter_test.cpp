#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "tauamp/tautau/ordered_channel_adapter.h"

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
    const std::vector<ReconstructedPion> pi_pair{object(1, 211, 0.4, 0.0), object(2, -211, 0.4, 0.0)};
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::pi_pi, pi_pair, {}).size() == 1U);
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::pi_rho, pi_pair,
                                                 {object(10, 111, 0.4, pi)}).size() == 1U);
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::rho_pi, pi_pair,
                                                 {object(10, 111, 0.4, pi)}).size() == 1U);
    const std::vector<ReconstructedPion> rho_pair{object(1, 211, 0.4, 0.0), object(2, -211, 0.4, pi)};
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::rho_rho, rho_pair,
                                                 {object(10, 111, 0.4, pi), object(11, 111, 0.4, 0.0)}).size() == 1U);
    const std::vector<ReconstructedPion> mixed{
        object(1, 211, 0.4, 0.0), object(2, 211, 0.4, 0.0),
        object(3, -211, 0.4, 0.0), object(4, -211, 0.4, 0.0)};
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::pi_a1, mixed, {}).size() == 2U);
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::a1_pi, mixed, {}).size() == 2U);
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::rho_a1, mixed,
                                                 {object(10, 111, 0.4, pi)}).size() == 2U);
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::a1_rho, mixed,
                                                 {object(10, 111, 0.4, pi)}).size() == 2U);
    std::vector<ReconstructedPion> a1_pair;
    for (std::size_t index = 1; index <= 3; ++index) a1_pair.push_back(object(index, 211, 0.3 + 0.05 * index, 0.0));
    for (std::size_t index = 4; index <= 6; ++index) a1_pair.push_back(object(index, -211, 0.3 + 0.03 * index, 0.0));
    assert(enumerate_ordered_channel_hypotheses(OrderedChannel::a1_a1, a1_pair, {}).size() == 9U);
}
