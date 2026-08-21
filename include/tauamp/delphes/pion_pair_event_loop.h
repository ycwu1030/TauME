#ifndef TAUAMP_DELPHES_PION_PAIR_EVENT_LOOP_H_
#define TAUAMP_DELPHES_PION_PAIR_EVENT_LOOP_H_

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "tauamp/delphes/pion_pair_evaluation.h"

namespace tauamp::delphes {

enum class PionPairSelectionStatus {
    unsupported_input,
    kinematic_no_solution,
    kinematic_non_unique,
    evaluable_one_solution,
    evaluable_two_solutions
};

struct PionPairEventLoopRecord {
    Long64_t entry;
    TauPairDecayClassification classification;
    TauMomentumProvenance tau_momentum_provenance;
    PionPairSelectionStatus selection_status;
    std::optional<tautau::PionPairVisibleEventEvaluation> evaluation;
    std::optional<PionPairLabObservation> observation;
};

struct PionPairEventLoopSummary {
    Long64_t entries_visited{};
    Long64_t supported_input{};
    Long64_t unsupported_input{};
    Long64_t kinematic_no_solution{};
    Long64_t kinematic_non_unique{};
    Long64_t evaluable_one_solution{};
    Long64_t evaluable_two_solutions{};
};

inline PionPairEventLoopRecord make_pion_pair_event_loop_record(Long64_t entry,
                                                                 PionPairDelphesEvaluation evaluation) {
    if (evaluation.classification != TauPairDecayClassification::supported_pion_pair) {
        if (evaluation.evaluation.has_value())
            throw std::invalid_argument("unsupported pion-pair event must not carry an evaluation");
        return {entry, evaluation.classification, evaluation.tau_momentum_provenance,
                PionPairSelectionStatus::unsupported_input, std::nullopt, std::nullopt};
    }
    if (!evaluation.evaluation.has_value())
        throw std::invalid_argument("supported pion-pair event must carry a visible-event evaluation");

    const tautau::HadronicTauPairSolutionStatus kinematic_status = evaluation.evaluation->status();
    PionPairSelectionStatus selection_status;
    switch (kinematic_status) {
        case tautau::HadronicTauPairSolutionStatus::no_solution:
            selection_status = PionPairSelectionStatus::kinematic_no_solution;
            break;
        case tautau::HadronicTauPairSolutionStatus::non_unique:
            selection_status = PionPairSelectionStatus::kinematic_non_unique;
            break;
        case tautau::HadronicTauPairSolutionStatus::one_solution:
            selection_status = PionPairSelectionStatus::evaluable_one_solution;
            break;
        case tautau::HadronicTauPairSolutionStatus::two_solutions:
            selection_status = PionPairSelectionStatus::evaluable_two_solutions;
            break;
    }
    return {entry, evaluation.classification, evaluation.tau_momentum_provenance, selection_status,
            std::move(evaluation.evaluation), std::move(evaluation.observation)};
}

inline void account_pion_pair_event(PionPairEventLoopSummary& summary, const PionPairEventLoopRecord& record) {
    ++summary.entries_visited;
    if (record.classification == TauPairDecayClassification::supported_pion_pair) {
        ++summary.supported_input;
    } else {
        ++summary.unsupported_input;
    }
    switch (record.selection_status) {
        case PionPairSelectionStatus::unsupported_input:
            break;
        case PionPairSelectionStatus::kinematic_no_solution:
            ++summary.kinematic_no_solution;
            break;
        case PionPairSelectionStatus::kinematic_non_unique:
            ++summary.kinematic_non_unique;
            break;
        case PionPairSelectionStatus::evaluable_one_solution:
            ++summary.evaluable_one_solution;
            break;
        case PionPairSelectionStatus::evaluable_two_solutions:
            ++summary.evaluable_two_solutions;
            break;
    }
}

class PionPairDelphesEventLoop {
public:
    PionPairDelphesEventLoop(const PionPairDelphesReader& reader, const PionPairDelphesEventEvaluator& evaluator)
        : reader_(reader), evaluator_(evaluator) {}

    template <class Consumer>
    PionPairEventLoopSummary for_each(Long64_t begin, Long64_t end, Consumer consumer) const {
        if (begin < 0 || end < begin || end > reader_.entry_count())
            throw std::out_of_range("Delphes event-loop range is outside reader entries");
        PionPairEventLoopSummary summary;
        for (Long64_t entry = begin; entry < end; ++entry) {
            PionPairEventLoopRecord record = process_entry(entry);
            account_pion_pair_event(summary, record);
            consumer(record);
        }
        return summary;
    }

private:
    PionPairEventLoopRecord process_entry(Long64_t entry) const {
        try {
            return make_pion_pair_event_loop_record(entry, evaluator_.evaluate(reader_.read(entry)));
        } catch (const std::exception& error) {
            throw std::runtime_error("cannot process Delphes entry " + std::to_string(entry) + ": " + error.what());
        }
    }

    const PionPairDelphesReader& reader_;
    const PionPairDelphesEventEvaluator& evaluator_;
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_EVENT_LOOP_H_
