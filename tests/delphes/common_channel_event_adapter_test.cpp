#include <cassert>
#include <cstdio>
#include <memory>
#include <string>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "tauamp/delphes/common_channel_event_adapter.h"
#include "tauamp/delphes/common_channel_root_writer.h"

namespace {
void append_particle(TClonesArray& collection, int index, int pid, int charge, double px, double py, double pz,
                     double energy) {
    auto* particle = new (collection[index]) GenParticle();
    particle->PID = pid;
    particle->Charge = charge;
    particle->Px = static_cast<Float_t>(px);
    particle->Py = static_cast<Float_t>(py);
    particle->Pz = static_cast<Float_t>(pz);
    particle->E = static_cast<Float_t>(energy);
}

void write_rho_pi_fixture(const std::string& path) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
    assert(file && !file->IsZombie());
    TTree tree("Delphes", "rho-pi adapter fixture");
    TClonesArray* charged = new TClonesArray("GenParticle");
    TClonesArray* neutral = new TClonesArray("GenParticle");
    tree.Branch("ChargedPion", &charged, 32000, 0);
    tree.Branch("NeutralPion", &neutral, 32000, 0);
    const double p = 0.362;
    const double ec = std::sqrt(p * p + 0.13957018 * 0.13957018);
    const double en = std::sqrt(p * p + 0.1349766 * 0.1349766);
    append_particle(*charged, 0, 211, 1, p, 0.0, 0.0, ec);
    append_particle(*charged, 1, -211, -1, 0.0, 0.0, 0.0, 0.13957018);
    append_particle(*neutral, 0, 111, 0, -p, 0.0, 0.0, en);
    tree.Fill();
    tree.Write();
    file->Write();
    file->Close();
    delete charged;
    delete neutral;
}
}

int main() {
    const std::string path = "/private/tmp/tauamp_common_channel_event_adapter_fixture.root";
    std::remove(path.c_str());
    write_rho_pi_fixture(path);
    tauamp::delphes::CommonChannelDelphesReader reader(path, {{}, {}});
    const auto collections = reader.read(0);
    const auto record = tauamp::delphes::build_common_channel_event_record(
        17, tauamp::tautau::OrderedChannel::rho_pi, collections, "selected::rho_pi",
        [](const tauamp::tautau::OrderedChannelHypothesis&) {
            tauamp::tautau::CommonHypothesisRecord result;
            result.weight = 1.0;
            result.components[0] = 1.0;
            return result;
        });
    assert(record.candidate_hypothesis_count == 1);
    assert(record.hypotheses.size() == 1);
    assert(record.hypotheses.front().hypothesis_index == 0);
    assert(record.truth_ordered_channel == "rho_pi");

    auto rejected = collections;
    rejected.neutral.clear();
    const auto rejection = tauamp::delphes::try_build_common_channel_event_record(
        18, tauamp::tautau::OrderedChannel::rho_pi, rejected, "diagnostic::rho_pi",
        [](const tauamp::tautau::OrderedChannelHypothesis&) {
            tauamp::tautau::CommonHypothesisRecord result;
            result.weight = 1.0;
            result.components[0] = 1.0;
            return result;
        });
    assert(!rejection.selected);
    assert(rejection.record.hypotheses.empty());
    assert(!rejection.diagnostics.empty());
    assert(rejection.diagnostics.front() == "rejected::no_ownership_candidate");
    assert(rejection.diagnostics.back() == "rejected::rho_mass_or_object_requirement");

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
    metadata.tauoo_revision = "fixture";
    metadata.exporter_revision = "fixture";
    metadata.source_hashes = {{"fixture", std::string(64, 'b')}};
    metadata.cross_feed_policy = "retain-overlap-record-v1";
    const std::string output = "/private/tmp/tauamp_common_channel_event_adapter_output.root";
    std::remove(output.c_str());
    tauamp::delphes::CommonChannelRootWriter writer(output, metadata);
    writer.write(record);
    writer.finalize();
    std::unique_ptr<TFile> output_file(TFile::Open(output.c_str(), "READ"));
    assert(output_file && !output_file->IsZombie());
    auto* events = dynamic_cast<TTree*>(output_file->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(output_file->Get("Hypotheses"));
    assert(events && hypotheses);
    assert(events->GetEntries() == 1);
    assert(hypotheses->GetEntries() == 1);
    std::remove(output.c_str());
    std::remove(path.c_str());
}
