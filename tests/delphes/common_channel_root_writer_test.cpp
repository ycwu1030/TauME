#include <cassert>
#include <cstdio>
#include <memory>
#include <string>

#include <TFile.h>
#include <TKey.h>
#include <TObjString.h>
#include <TTree.h>

#include "tauamp/delphes/common_channel_root_writer.h"

namespace {

tauamp::tautau::CommonChannelMetadata metadata() {
    tauamp::tautau::CommonChannelMetadata value;
    value.sqrt_s_GeV = 4.26;
    value.momentum_frame = "lab";
    value.object_source = "reconstructed";
    value.branch_weight_rule = "candidate_then_tauoo_branch_marginalization_v1";
    value.ordered_channel = tauamp::tautau::OrderedChannel::rho_a1;
    value.tau_plus_decay = "a1+ nu~";
    value.tau_minus_decay = "rho- nu";
    value.decay_model_revision = "g6-rho-a1-primary-v1";
    value.current_revision = "g6-rho-a1-current-v1";
    value.lineshape_revision = "g6-rho-kw-constant-width-v1";
    value.phase_space_measure_revision = "g6-lorentz-invariant-four-body-massless-neutrino-v1";
    value.current_normalization_revision = "g6-a1-reduced-normalization-GKW1-Na1-1-v1";
    value.identical_slot_rule = "descending-reconstructed-pT;stable-object-index-tiebreak";
    value.charge_convention_revision = "g6-charge-labelled-contraction-v1";
    value.analyser_convention_revision = "g6-p43-p44-primary-v1";
    value.pi0_object_source = "NeutralPion";
    value.selection_revision = "g6-rho-a1-candidate-v1";
    value.tauoo_revision = "fixture";
    value.exporter_revision = "fixture";
    value.source_hashes = {{"fixture", std::string(64, 'a')}};
    value.cross_feed_policy = "retain-overlap-record-v1";
    value.failure_record = "none";
    return value;
}

}  // namespace

int main() {
    const std::string path = "/private/tmp/tauamp_common_channel_root_writer_test.root";
    std::remove(path.c_str());
    tauamp::delphes::CommonChannelRootWriter writer(path, metadata());

    tauamp::tautau::CommonHypothesisRecord first;
    first.hypothesis_index = 0;
    first.weight = 0.5;
    first.valid = true;
    first.components[0] = 1.0;
    first.components[1] = 2.0;
    tauamp::tautau::CommonHypothesisRecord second = first;
    second.hypothesis_index = 1;
    second.weight = 0.5;
    second.components[0] = 3.0;
    tauamp::tautau::CommonEventRecord event{7, tauamp::tautau::OrderedChannel::rho_a1,
                                            "selected::rho_a1", "rho_a1", 2, {first, second}};
    writer.write(event);
    writer.finalize();

    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    assert(file && !file->IsZombie());
    const auto* manifest = dynamic_cast<TObjString*>(file->Get("TauAmpCommonChannelMetadata"));
    assert(manifest != nullptr);
    const std::string text = manifest->GetString().Data();
    assert(text.find("schema_version=2\n") != std::string::npos);
    assert(text.find("ordered_channel=rho_a1\n") != std::string::npos);
    assert(text.find("sqrt_s_GeV=4.26\n") != std::string::npos);
    assert(text.find("event_count=1\n") != std::string::npos);
    assert(text.find("hypothesis_count=2\n") != std::string::npos);
    assert(text.find("source_hash.fixture=" + std::string(64, 'a') + "\n") != std::string::npos);

    auto* events = dynamic_cast<TTree*>(file->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(file->Get("Hypotheses"));
    assert(events && hypotheses);
    assert(events->GetEntries() == 1);
    assert(hypotheses->GetEntries() == 2);

    Long64_t source_event_id = -1;
    Int_t channel = -1;
    Int_t hypothesis_index = -1;
    Double_t weight = 0.0;
    Double_t components[15]{};
    events->SetBranchAddress("source_event_id", &source_event_id);
    events->SetBranchAddress("ordered_channel", &channel);
    events->GetEntry(0);
    assert(source_event_id == 7);
    assert(channel == static_cast<Int_t>(tauamp::tautau::OrderedChannel::rho_a1));
    hypotheses->SetBranchAddress("hypothesis_index", &hypothesis_index);
    hypotheses->SetBranchAddress("kinematic_weight", &weight);
    hypotheses->SetBranchAddress("components", components);
    hypotheses->GetEntry(1);
    assert(hypothesis_index == 1);
    assert(weight == 0.5);
    assert(components[0] == 3.0);
}
