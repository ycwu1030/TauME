#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "tauamp/delphes/common_channel_root_writer.h"
#include "tauamp/tautau/ordered_channel_evaluation.h"

namespace {
using tauamp::tautau::FourMomentum;
using tauamp::tautau::ReconstructedPion;
constexpr double tau_mass = 1.77686;
constexpr double charged_mass = 0.13957018;
constexpr double neutral_mass = 0.1349766;

std::array<double, 3> unit(std::array<double, 3> value) {
    const double norm = std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
    for (double& component : value) component /= norm;
    return value;
}

FourMomentum visible_in_tau_rest(double mass, const std::array<double, 3>& direction) {
    const double momentum = (tau_mass * tau_mass - mass * mass) / (2.0 * tau_mass);
    const double energy = (tau_mass * tau_mass + mass * mass) / (2.0 * tau_mass);
    return {momentum * direction[0], momentum * direction[1], momentum * direction[2], energy};
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
                                                    const FourMomentum& tau_cm) {
    const double first = tauamp::tautau::g6_rho_mass * tauamp::tautau::g6_rho_mass -
                         std::pow(charged_mass + neutral_mass, 2);
    const double second = tauamp::tautau::g6_rho_mass * tauamp::tautau::g6_rho_mass -
                          std::pow(charged_mass - neutral_mass, 2);
    const double momentum = std::sqrt(first * second) / (2.0 * tauamp::tautau::g6_rho_mass);
    const auto direction = unit({{0.3, -0.4, 0.8}});
    const FourMomentum charged_rest(momentum * direction[0], momentum * direction[1], momentum * direction[2],
                                    std::sqrt(momentum * momentum + charged_mass * charged_mass));
    const FourMomentum neutral_rest(-momentum * direction[0], -momentum * direction[1], -momentum * direction[2],
                                    std::sqrt(momentum * momentum + neutral_mass * neutral_mass));
    const FourMomentum charged_tau_rest = boost_from_parent_rest(charged_rest, rho_tau_rest);
    const FourMomentum neutral_tau_rest = boost_from_parent_rest(neutral_rest, rho_tau_rest);
    return {boost_from_parent_rest(charged_tau_rest, tau_cm), boost_from_parent_rest(neutral_tau_rest, tau_cm)};
}

void check_weights(const tauamp::tautau::OrderedEventEvaluation& evaluation) {
    double sum = 0.0;
    for (const auto& branch : evaluation.branches) {
        assert(branch.valid);
        assert(branch.weight > 0.0 && std::isfinite(branch.weight));
        assert(branch.components[0] > 0.0 && std::isfinite(branch.components[0]));
        for (double value : branch.components) assert(std::isfinite(value));
        sum += branch.weight;
    }
    assert(std::abs(sum - 1.0) < 1.0e-12);
}

tauamp::tautau::CommonChannelMetadata rho_pi_metadata() {
    tauamp::tautau::CommonChannelMetadata metadata;
    metadata.ordered_channel = tauamp::tautau::OrderedChannel::rho_pi;
    metadata.momentum_frame = "lab";
    metadata.object_source = "reconstructed";
    metadata.branch_weight_rule = "candidate_then_tauoo_branch_marginalization_v1";
    metadata.tau_plus_decay = "rho+ nu";
    metadata.tau_minus_decay = "pi- nu~";
    metadata.decay_model_revision = "g6-rho-a1-primary-v1";
    metadata.current_revision = "g6-rho-transverse-rt-v1";
    metadata.lineshape_revision = "g6-rho-kw-constant-width-v1";
    metadata.phase_space_measure_revision = "g6-lorentz-invariant-two-body-v1";
    metadata.current_normalization_revision = "g6-rho-reduced-normalization-v1";
    metadata.identical_slot_rule = "not_applicable";
    metadata.charge_convention_revision = "g6-charge-labelled-contraction-v1";
    metadata.analyser_convention_revision = "g6-p43-p44-primary-v1";
    metadata.pi0_object_source = "NeutralPion";
    metadata.selection_revision = "g6-rho-a1-candidate-v1";
    metadata.tauoo_revision = "modern-g6-physical-evaluator-v1";
    metadata.exporter_revision = "common-root-writer-v1";
    metadata.source_hashes = {{"ordered_channel_evaluation_test", std::string(64, 'c')}};
    metadata.cross_feed_policy = "retain-overlap-record-v1";
    return metadata;
}
}

int main() {
    using namespace tauamp::tautau;
    constexpr double beam_energy = 2.13;
    constexpr double theta = 0.83;
    const double tau_momentum = std::sqrt(beam_energy * beam_energy - tau_mass * tau_mass);
    const FourMomentum electron(0.0, 0.0, -beam_energy, beam_energy);
    const FourMomentum positron(0.0, 0.0, beam_energy, beam_energy);
    const FourMomentum tau_minus(tau_momentum * std::sin(theta), 0.0, tau_momentum * std::cos(theta), beam_energy);
    const FourMomentum tau_plus(-tau_minus.px(), 0.0, -tau_minus.pz(), beam_energy);
    const BeamState beams(electron, positron, 0.0, 0.0);
    const HadronicPairMatrixElement matrix(ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952},
                                           ProductionBosons::photon_only);

    const FourMomentum rho_tau_rest = visible_in_tau_rest(g6_rho_mass, unit({{0.2, 0.5, 0.7}}));
    const auto [rho_charged, rho_neutral] = rho_daughters(rho_tau_rest, tau_plus);
    const FourMomentum pion_minus = boost_from_parent_rest(
        visible_in_tau_rest(charged_mass, unit({{-0.4, 0.3, 0.8}})), tau_minus);
    std::vector<ReconstructedPion> rho_pi_charged{reco(0, 211, rho_charged), reco(1, -211, pion_minus)};
    std::vector<ReconstructedPion> rho_pi_neutral{reco(0, 111, rho_neutral)};
    const auto rho_pi = evaluate_ordered_event(OrderedChannel::rho_pi, rho_pi_charged, rho_pi_neutral,
                                               beams, tau_mass, matrix);
    assert(rho_pi.candidate_count == 1);
    assert(rho_pi.valid_candidate_count == 1);
    assert(!rho_pi.branches.empty());
    check_weights(rho_pi);
    CommonEventRecord rho_pi_record = make_ordered_common_event_record(1001, OrderedChannel::rho_pi, rho_pi);
    assert(rho_pi_record.hypotheses.size() == rho_pi.branches.size());

    const std::string output = "/private/tmp/tauamp_ordered_physical_evaluation_roundtrip.root";
    std::remove(output.c_str());
    tauamp::delphes::CommonChannelRootWriter writer(output, rho_pi_metadata());
    writer.write(rho_pi_record);
    writer.finalize();
    std::unique_ptr<TFile> output_file(TFile::Open(output.c_str(), "READ"));
    assert(output_file && !output_file->IsZombie());
    auto* events = dynamic_cast<TTree*>(output_file->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(output_file->Get("Hypotheses"));
    assert(events && hypotheses);
    assert(events->GetEntries() == 1);
    assert(hypotheses->GetEntries() == static_cast<Long64_t>(rho_pi_record.hypotheses.size()));
    std::remove(output.c_str());

    constexpr double a1_mass = 1.2;
    const FourMomentum a1_tau_rest = visible_in_tau_rest(a1_mass, unit({{-0.3, 0.4, 0.7}}));
    const double daughter_energy = a1_mass / 3.0;
    const double daughter_momentum = std::sqrt(daughter_energy * daughter_energy - charged_mass * charged_mass);
    const FourMomentum q1_rest(daughter_momentum, 0.0, 0.0, daughter_energy);
    const FourMomentum q2_rest(-0.5 * daughter_momentum, std::sqrt(0.75) * daughter_momentum, 0.0, daughter_energy);
    const FourMomentum qop_rest(-0.5 * daughter_momentum, -std::sqrt(0.75) * daughter_momentum, 0.0, daughter_energy);
    const FourMomentum q1 = boost_from_parent_rest(boost_from_parent_rest(q1_rest, a1_tau_rest), tau_plus);
    const FourMomentum q2 = boost_from_parent_rest(boost_from_parent_rest(q2_rest, a1_tau_rest), tau_plus);
    const FourMomentum qop = boost_from_parent_rest(boost_from_parent_rest(qop_rest, a1_tau_rest), tau_plus);
    const FourMomentum direct_minus = boost_from_parent_rest(
        visible_in_tau_rest(charged_mass, unit({{0.5, -0.2, 0.7}})), tau_minus);
    std::vector<ReconstructedPion> a1_pi_charged{
        reco(0, 211, q1), reco(1, 211, q2), reco(2, -211, qop), reco(3, -211, direct_minus)};
    const auto a1_pi = evaluate_ordered_event(OrderedChannel::a1_pi, a1_pi_charged, {}, beams, tau_mass, matrix);
    assert(a1_pi.candidate_count == 2);
    assert(a1_pi.valid_candidate_count >= 1);
    assert(!a1_pi.branches.empty());
    check_weights(a1_pi);
    CommonEventRecord a1_pi_record = make_ordered_common_event_record(1002, OrderedChannel::a1_pi, a1_pi);
    assert(a1_pi_record.hypotheses.size() == a1_pi.branches.size());
}
