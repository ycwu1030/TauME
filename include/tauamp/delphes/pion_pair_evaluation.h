#ifndef TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_
#define TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_

#include <optional>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "tauamp/delphes/pion_pair_reader.h"
#include "tauamp/tautau/visible_pion_pair_evaluation.h"

namespace tauamp::delphes {

struct PionPairDelphesEvaluation {
    TauPairDecayClassification classification;
    TauMomentumProvenance tau_momentum_provenance;
    std::optional<tautau::PionPairVisibleEventEvaluation> evaluation;
    std::optional<PionPairLabObservation> observation{};
};

class PionPairDelphesEventEvaluator {
public:
    explicit PionPairDelphesEventEvaluator(
        tautau::ElectroweakParameters parameters,
        tautau::ProductionBosons boson_selection = tautau::ProductionBosons::photon_and_z)
        : evaluator_(parameters, boson_selection) {}

    PionPairDelphesEvaluation evaluate(const PionPairDelphesEvent& event) const {
        if (event.classification != TauPairDecayClassification::supported_pion_pair)
            return {event.classification, event.tau_momentum_provenance, std::nullopt, std::nullopt};
        if (!event.observation.has_value())
            throw std::invalid_argument("supported pion-pair event lacks a lab-frame observation");
        const PionPairLabObservation& observation = *event.observation;
        if (observation.truth_tau_minus_lab.has_value() && observation.truth_tau_plus_lab.has_value()) {
            const double truth_minus_mass = std::sqrt(observation.truth_tau_minus_lab->mass_squared());
            const double truth_plus_mass = std::sqrt(observation.truth_tau_plus_lab->mass_squared());
            const tautau::TauPairKinematicPoint point(
                observation.beams, *observation.truth_tau_minus_lab, *observation.truth_tau_plus_lab,
                truth_minus_mass, truth_plus_mass, false);
            const tautau::PionPairKinematicHypothesis hypothesis{
                point, observation.pion_minus_lab, observation.pion_plus_lab, tautau::KinematicPointOrigin::truth};
            return {event.classification, event.tau_momentum_provenance,
                    evaluator_.evaluate_truth(point, observation.pion_minus_lab, observation.pion_plus_lab),
                    event.observation};
        }
        return {event.classification, event.tau_momentum_provenance,
                evaluator_.evaluate(observation.beams, observation.pion_minus_lab, observation.pion_plus_lab,
                                    observation.tau_mass),
                event.observation};
    }


private:
    tautau::PionPairVisibleEventEvaluator evaluator_;
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_EVALUATION_H_
