#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "tauamp/tautau/ordered_channel_evaluation.h"

namespace {
using namespace tauamp::tautau;
constexpr double tau_mass = 1.77686;
constexpr double charged_mass = 0.13957018;
constexpr double neutral_mass = 0.1349766;
constexpr double a1_mass = 1.2;

std::array<double, 3> unit(std::array<double, 3> value) {
    const double norm = std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
    for (double& component : value) component /= norm;
    return value;
}

FourMomentum tau_rest_visible(double mass, const std::array<double, 3>& direction) {
    const double momentum = (tau_mass * tau_mass - mass * mass) / (2.0 * tau_mass);
    const double energy = (tau_mass * tau_mass + mass * mass) / (2.0 * tau_mass);
    const auto direction_unit = unit(direction);
    return {momentum * direction_unit[0], momentum * direction_unit[1], momentum * direction_unit[2], energy};
}

FourMomentum boost_from_parent_rest(const FourMomentum& daughter, const FourMomentum& parent) {
    return daughter.boosted({{parent.px() / parent.energy(), parent.py() / parent.energy(), parent.pz() / parent.energy()}});
}

ReconstructedPion reco(std::size_t index, int pid, const FourMomentum& momentum) {
    const double pt = std::hypot(momentum.px(), momentum.py());
    const double total = std::sqrt(momentum.px() * momentum.px() + momentum.py() * momentum.py() + momentum.pz() * momentum.pz());
    const double eta = total == std::abs(momentum.pz()) ? (momentum.pz() >= 0.0 ? 1.0e9 : -1.0e9)
                                                        : 0.5 * std::log((total + momentum.pz()) / (total - momentum.pz()));
    return {index, pid, pt, eta, std::atan2(momentum.py(), momentum.px()), momentum};
}

std::pair<FourMomentum, FourMomentum> rho_daughters(const FourMomentum& rho_tau_rest,
                                                    const FourMomentum& tau_lab, double shift) {
    const double first = g6_rho_mass * g6_rho_mass - std::pow(charged_mass + neutral_mass, 2);
    const double second = g6_rho_mass * g6_rho_mass - std::pow(charged_mass - neutral_mass, 2);
    const double momentum = std::sqrt(first * second) / (2.0 * g6_rho_mass);
    const auto direction = unit({{0.3 + shift, -0.4, 0.8}});
    const FourMomentum charged_rest(momentum * direction[0], momentum * direction[1], momentum * direction[2],
                                    std::sqrt(momentum * momentum + charged_mass * charged_mass));
    const FourMomentum neutral_rest(-momentum * direction[0], -momentum * direction[1], -momentum * direction[2],
                                    std::sqrt(momentum * momentum + neutral_mass * neutral_mass));
    return {boost_from_parent_rest(boost_from_parent_rest(charged_rest, rho_tau_rest), tau_lab),
            boost_from_parent_rest(boost_from_parent_rest(neutral_rest, rho_tau_rest), tau_lab)};
}

std::array<FourMomentum, 3> a1_daughters(const FourMomentum& a1_tau_rest, const FourMomentum& tau_lab,
                                         double shift) {
    const double energy = a1_mass / 3.0;
    const double momentum = std::sqrt(energy * energy - charged_mass * charged_mass);
    const double c = 0.45 + shift;
    const double s = std::sqrt(1.0 - c * c);
    const FourMomentum q1(momentum * c, momentum * s, 0.0, energy);
    const FourMomentum q2(-0.5 * momentum * c, -0.5 * momentum * s, std::sqrt(0.75) * momentum, energy);
    const FourMomentum q3(-0.5 * momentum * c, -0.5 * momentum * s, -std::sqrt(0.75) * momentum, energy);
    return {{boost_from_parent_rest(boost_from_parent_rest(q1, a1_tau_rest), tau_lab),
             boost_from_parent_rest(boost_from_parent_rest(q2, a1_tau_rest), tau_lab),
             boost_from_parent_rest(boost_from_parent_rest(q3, a1_tau_rest), tau_lab)}};
}

struct Fixture { std::vector<ReconstructedPion> charged; std::vector<ReconstructedPion> neutral; };

Fixture fixture(OrderedChannel channel, const FourMomentum& tau_plus, const FourMomentum& tau_minus) {
    Fixture result;
    std::size_t charged_index = 0;
    const auto plus = [&](const FourMomentum& p) { result.charged.push_back(reco(charged_index++, 211, p)); };
    const auto minus = [&](const FourMomentum& p) { result.charged.push_back(reco(charged_index++, -211, p)); };
    const auto pion_plus = [&] { plus(boost_from_parent_rest(tau_rest_visible(charged_mass, {{0.2, 0.4, 0.8}}), tau_plus)); };
    const auto pion_minus = [&] { minus(boost_from_parent_rest(tau_rest_visible(charged_mass, {{-0.2, 0.4, 0.8}}), tau_minus)); };
    const auto rho_plus = [&](double shift) {
        const auto daughters = rho_daughters(tau_rest_visible(g6_rho_mass, {{0.2 + shift, 0.5, 0.7}}), tau_plus, shift);
        plus(daughters.first); result.neutral.push_back(reco(result.neutral.size(), 111, daughters.second));
    };
    const auto rho_minus = [&](double shift) {
        const auto daughters = rho_daughters(tau_rest_visible(g6_rho_mass, {{-0.2 + shift, 0.5, 0.7}}), tau_minus, shift);
        minus(daughters.first); result.neutral.push_back(reco(result.neutral.size(), 111, daughters.second));
    };
    const auto a1_plus = [&] {
        const auto daughters = a1_daughters(tau_rest_visible(a1_mass, {{-0.3, 0.4, 0.7}}), tau_plus, 0.0);
        plus(daughters[0]); plus(daughters[1]); minus(daughters[2]);
    };
    const auto a1_minus = [&] {
        const auto daughters = a1_daughters(tau_rest_visible(a1_mass, {{0.3, 0.4, 0.7}}), tau_minus, 0.0);
        minus(daughters[0]); minus(daughters[1]); plus(daughters[2]);
    };
    switch (channel) {
        case OrderedChannel::pi_pi: pion_plus(); pion_minus(); break;
        case OrderedChannel::pi_rho: pion_plus(); rho_minus(0.0); break;
        case OrderedChannel::rho_pi: rho_plus(0.0); pion_minus(); break;
        case OrderedChannel::rho_rho: rho_plus(0.0); rho_minus(0.12); break;
        case OrderedChannel::pi_a1: pion_plus(); a1_minus(); break;
        case OrderedChannel::a1_pi: a1_plus(); pion_minus(); break;
        case OrderedChannel::rho_a1: rho_minus(0.0); a1_plus(); break;
        case OrderedChannel::a1_rho: a1_minus(); rho_plus(0.0); break;
        case OrderedChannel::a1_a1: a1_plus(); a1_minus(); break;
    }
    return result;
}

void check(const OrderedEventEvaluation& evaluation) {
    assert(evaluation.candidate_count > 0U);
    assert(evaluation.valid_candidate_count > 0U);
    double total = 0.0;
    for (const auto& branch : evaluation.branches) {
        assert(branch.valid && std::isfinite(branch.weight) && branch.weight > 0.0);
        assert(std::isfinite(branch.components[0]) && branch.components[0] > 0.0);
        total += branch.weight;
    }
    assert(std::abs(total - 1.0) < 1.0e-12);
}
}

int main() {
    constexpr double beam_energy = 2.13;
    constexpr double theta = 0.83;
    const double tau_momentum = std::sqrt(beam_energy * beam_energy - tau_mass * tau_mass);
    const FourMomentum electron(0.0, 0.0, -beam_energy, beam_energy);
    const FourMomentum positron(0.0, 0.0, beam_energy, beam_energy);
    const FourMomentum tau_minus(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), beam_energy);
    const FourMomentum tau_plus(-tau_minus.px(), 0.0, -tau_minus.pz(), beam_energy);
    const BeamState beams(electron, positron, 0.0, 0.0);
    const HadronicPairMatrixElement matrix(ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952}, ProductionBosons::photon_only);
    const std::array<OrderedChannel, 9> channels{{OrderedChannel::pi_pi, OrderedChannel::pi_rho, OrderedChannel::rho_pi,
        OrderedChannel::rho_rho, OrderedChannel::pi_a1, OrderedChannel::a1_pi, OrderedChannel::rho_a1,
        OrderedChannel::a1_rho, OrderedChannel::a1_a1}};
    for (const auto channel : channels) {
        const auto objects = fixture(channel, tau_plus, tau_minus);
        check(evaluate_ordered_event(channel, objects.charged, objects.neutral, beams, tau_mass, matrix));
    }
}
