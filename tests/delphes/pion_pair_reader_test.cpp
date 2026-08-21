#include <cassert>
#include <cmath>
#include <string>

#include "tauamp/delphes/pion_pair_reader.h"

namespace {
constexpr double kTolerance = 5e-6;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }
}  // namespace

int main(int argc, char* argv[]) {
    assert(argc == 2);
    const std::string path = argv[1];
    const tauamp::delphes::PionPairReaderConfig generated_config{0.0, 0.0, 1.77686,
                                                                  tauamp::delphes::PionObjectSource::generated};
    const tauamp::delphes::PionPairDelphesReader generated_reader(path, generated_config);
    const auto generated = generated_reader.read(0);
    assert(generated.classification == tauamp::delphes::TauPairDecayClassification::supported_pion_pair);
    assert(generated.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::generator_truth_available);
    assert(generated.observation.has_value());
    assert(close(generated.observation->beams.electron_lab().pz(), -2.13));
    assert(close(generated.observation->beams.positron_lab().pz(), 2.13));
    assert(close(generated.observation->pion_minus_lab.py(), -1.11791));
    assert(close(generated.observation->pion_plus_lab.py(), 0.79137));

    const tauamp::delphes::PionPairReaderConfig reconstructed_config{0.0, 0.0, 1.77686,
                                                                      tauamp::delphes::PionObjectSource::reconstructed};
    const tauamp::delphes::PionPairDelphesReader reconstructed_reader(path, reconstructed_config);
    const auto reconstructed = reconstructed_reader.read(0);
    assert(reconstructed.classification == tauamp::delphes::TauPairDecayClassification::supported_pion_pair);
    assert(reconstructed.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::unavailable);
    assert(reconstructed.observation.has_value());
    assert(close(reconstructed.observation->pion_minus_lab.py(), -1.11143));
    assert(close(reconstructed.observation->pion_plus_lab.py(), 0.785371));

    using tauamp::delphes::PionCandidate;
    using tauamp::delphes::TauPairDecayClassification;
    assert(tauamp::delphes::classify_pion_pair({PionCandidate{-211, -1}, PionCandidate{211, 1}}, 0) ==
           TauPairDecayClassification::supported_pion_pair);
    assert(tauamp::delphes::classify_pion_pair({PionCandidate{-211, -1}, PionCandidate{211, 1}}, 2) ==
           TauPairDecayClassification::recognized_unsupported);
    assert(tauamp::delphes::classify_pion_pair({PionCandidate{-211, -1}, PionCandidate{-211, -1}}, 0) ==
           TauPairDecayClassification::unclassified);
    assert(tauamp::delphes::classify_pion_pair({PionCandidate{-211, -1}, PionCandidate{211, 1}, PionCandidate{211, 1}}, 0) ==
           TauPairDecayClassification::ambiguous);
}
