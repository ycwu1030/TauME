#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "pion_pair_delphes_fixture.h"
#include "tauamp/delphes/pion_pair_event_loop.h"
#include "tauamp/delphes/pion_pair_reader.h"

namespace {

std::string shell_quote(const std::string& value) {
    std::string quoted{"'"};
    for (const char character : value) {
        if (character == '\'') {
            quoted += "'\\''";
        } else {
            quoted += character;
        }
    }
    return quoted + "'";
}

std::string read_file(const std::string& path) {
    std::ifstream input(path);
    std::ostringstream content;
    content << input.rdbuf();
    return content.str();
}

std::string caller_command(const std::string& caller_path, const std::string& input_path, const std::string& output_path,
                           const std::string& standard_output_path, const std::string& standard_error_path) {
    return shell_quote(caller_path) + " --input " + shell_quote(input_path) + " --output " + shell_quote(output_path) +
           " --begin 0 --end 1 --pion-source generated --electron-polarization 0 --positron-polarization 0"
           " --tau-mass 1.77686 --electric-charge 0.313 --weak-mixing-sine 0.48 --z-mass 91.1876 --z-width 2.4952"
           " --input-sha256 ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
           " --implementation-revision portable-fixture-test > " + shell_quote(standard_output_path) +
           " 2> " + shell_quote(standard_error_path);
}

}  // namespace

int main(int argc, char* argv[]) {
    assert(argc == 2);
    const std::string input_path = "/private/tmp/tauamp_pion_pair_delphes_fixture_input.root";
    const std::string output_path = "/private/tmp/tauamp_pion_pair_delphes_fixture_output.root";
    const std::string standard_output_path = "/private/tmp/tauamp_pion_pair_delphes_fixture.stdout";
    const std::string standard_error_path = "/private/tmp/tauamp_pion_pair_delphes_fixture.stderr";
    std::remove(input_path.c_str());
    std::remove(output_path.c_str());
    std::remove(standard_output_path.c_str());
    std::remove(standard_error_path.c_str());

    tauamp::delphes::test::write_pion_pair_delphes_fixture(input_path);

    const tauamp::delphes::PionPairReaderConfig generated_config{0.0, 0.0, 1.77686,
                                                                  tauamp::delphes::PionObjectSource::generated};
    const tauamp::delphes::PionPairDelphesReader generated_reader(input_path, generated_config);
    assert(generated_reader.entry_count() == 1);
    const auto generated = generated_reader.read(0);
    assert(generated.classification == tauamp::delphes::TauPairDecayClassification::supported_pion_pair);
    assert(generated.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::generator_truth_available);
    assert(generated.observation.has_value());
    assert(generated.observation->beams.electron_lab().pz() < 0.0);
    assert(generated.observation->beams.positron_lab().pz() > 0.0);
    assert(generated.observation->pion_minus_lab.py() > 0.0);
    assert(generated.observation->pion_plus_lab.pz() > 0.0);

    const tauamp::delphes::PionPairReaderConfig reconstructed_config{0.0, 0.0, 1.77686,
                                                                      tauamp::delphes::PionObjectSource::reconstructed};
    const tauamp::delphes::PionPairDelphesReader reconstructed_reader(input_path, reconstructed_config);
    const auto reconstructed = reconstructed_reader.read(0);
    assert(reconstructed.classification == tauamp::delphes::TauPairDecayClassification::supported_pion_pair);
    assert(reconstructed.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::unavailable);
    assert(reconstructed.observation.has_value());
    assert(reconstructed.observation->pion_minus_lab.py() > 0.0);
    assert(reconstructed.observation->pion_plus_lab.pz() > 0.0);

    const tauamp::delphes::PionPairDelphesEventEvaluator evaluator(
        tauamp::tautau::ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952});
    const tauamp::delphes::PionPairDelphesEventLoop loop(generated_reader, evaluator);
    std::vector<tauamp::delphes::PionPairEventLoopRecord> records;
    const auto summary = loop.for_each(0, 1, [&](const tauamp::delphes::PionPairEventLoopRecord& record) {
        records.push_back(record);
    });
    assert(summary.entries_visited == 1);
    assert(summary.evaluable_two_solutions == 1);
    assert(records.size() == 1U);
    assert(records.front().evaluation.has_value());
    assert(records.front().evaluation->entries().size() == 2U);

    const int caller_result = std::system(
        caller_command(argv[1], input_path, output_path, standard_output_path, standard_error_path).c_str());
    assert(caller_result == 0);
    assert(read_file(standard_output_path) ==
           "entries_visited=1\n"
           "supported_input=1\n"
           "unsupported_input=0\n"
           "kinematic_no_solution=0\n"
           "kinematic_non_unique=0\n"
           "evaluable_one_solution=0\n"
           "evaluable_two_solutions=1\n");

    std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "READ"));
    assert(output && !output->IsZombie());
    auto* events = dynamic_cast<TTree*>(output->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(output->Get("Hypotheses"));
    auto* manifest = dynamic_cast<TObjString*>(output->Get("TauAmpPionPairMetadata"));
    assert(events != nullptr && events->GetEntries() == 1);
    assert(hypotheses != nullptr && hypotheses->GetEntries() == 2);
    assert(manifest != nullptr);
    assert(std::string(manifest->GetString().Data()).find("input_file_sha256=ffffffff") != std::string::npos);
    assert(std::string(manifest->GetString().Data()).find("implementation_revision=portable-fixture-test") !=
           std::string::npos);

    output.reset();
    std::remove(input_path.c_str());
    std::remove(output_path.c_str());
    std::remove(standard_output_path.c_str());
    std::remove(standard_error_path.c_str());
}
