#include <cassert>
#include <stdexcept>
#include <string>

#include "tauamp/delphes/pion_pair_evaluation.h"
#include "tauamp/delphes/pion_pair_reader.h"

int main(int argc, char* argv[]) {
    assert(argc == 2);
    const tauamp::delphes::PionPairReaderConfig config{0.0, 0.0, 1.77686,
                                                        tauamp::delphes::PionObjectSource::generated};
    const tauamp::delphes::PionPairDelphesReader reader(argv[1], config);
    const tauamp::delphes::PionPairDelphesEventEvaluator evaluator(
        tauamp::tautau::ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952});

    const auto supported = evaluator.evaluate(reader.read(0));
    assert(supported.classification == tauamp::delphes::TauPairDecayClassification::supported_pion_pair);
    assert(supported.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::generator_truth_available);
    assert(supported.evaluation.has_value());

    const auto unsupported = evaluator.evaluate({tauamp::delphes::TauPairDecayClassification::recognized_unsupported,
                                                 tauamp::delphes::TauMomentumProvenance::unavailable, std::nullopt});
    assert(unsupported.classification == tauamp::delphes::TauPairDecayClassification::recognized_unsupported);
    assert(unsupported.tau_momentum_provenance == tauamp::delphes::TauMomentumProvenance::unavailable);
    assert(!unsupported.evaluation.has_value());

    bool rejected_inconsistent_supported_event = false;
    try {
        evaluator.evaluate({tauamp::delphes::TauPairDecayClassification::supported_pion_pair,
                            tauamp::delphes::TauMomentumProvenance::unavailable, std::nullopt});
    } catch (const std::invalid_argument&) {
        rejected_inconsistent_supported_event = true;
    }
    assert(rejected_inconsistent_supported_event);
}
