#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <string>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "tauamp/delphes/pion_pair_root_writer.h"

namespace {
constexpr double kTolerance = 5e-6;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }
}  // namespace

int main(int argc, char* argv[]) {
    assert(argc == 2);
    const std::string output_path = "/private/tmp/tauamp_pion_pair_root_writer_test.root";
    std::remove(output_path.c_str());

    const tauamp::delphes::PionPairReaderConfig reader_config{0.0, 0.0, 1.77686,
                                                               tauamp::delphes::PionObjectSource::generated};
    const tauamp::tautau::ElectroweakParameters electroweak_parameters{0.313, 0.48, 91.1876, 2.4952};
    const tauamp::delphes::PionPairDelphesReader reader(argv[1], reader_config);
    const tauamp::delphes::PionPairDelphesEventEvaluator evaluator(electroweak_parameters);
    const tauamp::delphes::PionPairDelphesEventLoop loop(reader, evaluator);
    const tauamp::delphes::PionPairRootOutputMetadata metadata{
        "a34d645cb9df4dc9a9527263277ebf3b14eb86aa9e4d74feab60744d28e84daa", "Delphes",
        tauamp::delphes::PionObjectSource::generated, reader_config.electron_polarization,
        reader_config.positron_polarization, reader_config.tau_mass, electroweak_parameters, "writer-test-revision"};
    tauamp::delphes::PionPairRootWriter writer(output_path, metadata);

    std::optional<tauamp::delphes::PionPairEventLoopRecord> raw_record;
    loop.for_each(0, 1, [&](const tauamp::delphes::PionPairEventLoopRecord& record) {
        raw_record = record;
        writer.write(record);
    });
    assert(raw_record.has_value());
    writer.write({1, tauamp::delphes::TauPairDecayClassification::recognized_unsupported,
                  tauamp::delphes::TauMomentumProvenance::unavailable,
                  tauamp::delphes::PionPairSelectionStatus::unsupported_input, std::nullopt, std::nullopt});

    bool rejected_duplicate = false;
    try {
        writer.write(*raw_record);
    } catch (const std::invalid_argument&) {
        rejected_duplicate = true;
    }
    assert(rejected_duplicate);
    writer.finalize();

    std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "READ"));
    assert(output && !output->IsZombie());
    const auto* manifest = dynamic_cast<TObjString*>(output->Get("TauAmpPionPairMetadata"));
    assert(manifest != nullptr);
    const std::string manifest_text = manifest->GetString().Data();
    assert(manifest_text.find("schema_version=1\n") != std::string::npos);
    assert(manifest_text.find("input_file_sha256=a34d645cb9df4dc9a9527263277ebf3b14eb86aa9e4d74feab60744d28e84daa\n") !=
           std::string::npos);
    assert(manifest_text.find("object_source=generated\n") != std::string::npos);
    assert(manifest_text.find("truth_tau_input_used=false\n") != std::string::npos);
    assert(manifest_text.find("implementation_revision=writer-test-revision\n") != std::string::npos);

    auto* events = dynamic_cast<TTree*>(output->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(output->Get("Hypotheses"));
    assert(events != nullptr);
    assert(hypotheses != nullptr);
    assert(events->GetEntries() == 2);
    const Long64_t expected_hypothesis_count = static_cast<Long64_t>(raw_record->evaluation->entries().size());
    assert(hypotheses->GetEntries() == expected_hypothesis_count);

    Long64_t source_entry = -1;
    Int_t classification = -1;
    Int_t provenance = -1;
    Int_t eligibility = -1;
    Int_t hypothesis_count = -1;
    Bool_t has_lab_observation = false;
    Double_t electron_pz = 0.0;
    Double_t pion_minus_py = 0.0;
    events->SetBranchAddress("source_entry", &source_entry);
    events->SetBranchAddress("classification", &classification);
    events->SetBranchAddress("tau_momentum_provenance", &provenance);
    events->SetBranchAddress("eligibility_status", &eligibility);
    events->SetBranchAddress("hypothesis_count", &hypothesis_count);
    events->SetBranchAddress("has_lab_observation", &has_lab_observation);
    events->SetBranchAddress("electron_pz", &electron_pz);
    events->SetBranchAddress("pion_minus_py", &pion_minus_py);

    events->GetEntry(0);
    assert(source_entry == 0);
    assert(classification == static_cast<Int_t>(raw_record->classification));
    assert(provenance == static_cast<Int_t>(raw_record->tau_momentum_provenance));
    assert(eligibility == static_cast<Int_t>(raw_record->selection_status));
    assert(hypothesis_count == expected_hypothesis_count);
    assert(has_lab_observation);
    assert(close(electron_pz, raw_record->observation->beams.electron_lab().pz()));
    assert(close(pion_minus_py, raw_record->observation->pion_minus_lab.py()));

    events->GetEntry(1);
    assert(source_entry == 1);
    assert(classification == static_cast<Int_t>(tauamp::delphes::TauPairDecayClassification::recognized_unsupported));
    assert(!has_lab_observation);
    assert(hypothesis_count == 0);
    assert(std::isnan(electron_pz));
    assert(std::isnan(pion_minus_py));

    if (expected_hypothesis_count > 0) {
        Int_t hypothesis_index = -1;
        Double_t kinematic_weight = 0.0;
        Double_t component_sm = 0.0;
        Double_t component_f2_real = 0.0;
        hypotheses->SetBranchAddress("source_entry", &source_entry);
        hypotheses->SetBranchAddress("hypothesis_index", &hypothesis_index);
        hypotheses->SetBranchAddress("kinematic_weight", &kinematic_weight);
        hypotheses->SetBranchAddress("component_sm", &component_sm);
        hypotheses->SetBranchAddress("component_f2_real", &component_f2_real);
        hypotheses->GetEntry(0);
        const auto& expected = raw_record->evaluation->entries().front();
        assert(source_entry == 0);
        assert(hypothesis_index == 0);
        assert(close(kinematic_weight, expected.weight));
        assert(close(component_sm, expected.components.sm));
        assert(close(component_f2_real, expected.components.f2_real));
    }
}
