#ifndef TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_
#define TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_

#include <optional>
#include <stdexcept>
#include <utility>

#include "tauamp/delphes/pion_pair_reader.h"
#include "tauamp/tautau/visible_pion_pair_evaluation.h"

namespace tauamp::delphes {

struct PionPairDelphesEvaluation {
    TauPairDecayClassification classification;
    TauMomentumProvenance tau_momentum_provenance;
    std::optional<tautau::PionPairVisibleEventEvaluation> evaluation;
};

class PionPairDelphesEventEvaluator {
public:
    explicit PionPairDelphesEventEvaluator(tautau::ElectroweakParameters parameters) : evaluator_(parameters) {}

    PionPairDelphesEvaluation evaluate(const PionPairDelphesEvent& event) const {
        if (event.classification != TauPairDecayClassification::supported_pion_pair)
            return {event.classification, event.tau_momentum_provenance, std::nullopt};
        if (!event.observation.has_value())
            throw std::invalid_argument("supported pion-pair event lacks a lab-frame observation");
        const PionPairLabObservation& observation = *event.observation;
        return {event.classification, event.tau_momentum_provenance,
                evaluator_.evaluate(observation.beams, observation.pion_minus_lab, observation.pion_plus_lab,
                                    observation.tau_mass)};
    }

private:
    tautau::PionPairVisibleEventEvaluator evaluator_;
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_
