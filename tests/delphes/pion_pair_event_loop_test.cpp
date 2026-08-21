#include <cassert>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tauamp/delphes/pion_pair_event_loop.h"
#include "tauamp/delphes/pion_pair_reader.h"

int main(int argc, char* argv[]) {
    assert(argc == 2);

    using tauamp::delphes::PionPairDelphesEvaluation;
    using tauamp::delphes::PionPairEventLoopRecord;
    using tauamp::delphes::PionPairSelectionStatus;
    using tauamp::delphes::TauMomentumProvenance;
    using tauamp::delphes::TauPairDecayClassification;
    using tauamp::tautau::HadronicTauPairSolutionStatus;
    using tauamp::tautau::PionPairVisibleEventEvaluation;

    const PionPairDelphesEvaluation unsupported{TauPairDecayClassification::recognized_unsupported,
                                                TauMomentumProvenance::unavailable, std::nullopt};
    const PionPairEventLoopRecord rejected = tauamp::delphes::make_pion_pair_event_loop_record(7, unsupported);
    assert(rejected.entry == 7);
    assert(rejected.selection_status == PionPairSelectionStatus::unsupported_input);
    assert(!rejected.evaluation.has_value());

    const PionPairDelphesEvaluation no_solution{TauPairDecayClassification::supported_pion_pair,
                                                TauMomentumProvenance::unavailable,
                                                PionPairVisibleEventEvaluation(HadronicTauPairSolutionStatus::no_solution, {})};
    assert(tauamp::delphes::make_pion_pair_event_loop_record(8, no_solution).selection_status ==
           PionPairSelectionStatus::kinematic_no_solution);

    const PionPairDelphesEvaluation non_unique{
        TauPairDecayClassification::supported_pion_pair, TauMomentumProvenance::unavailable,
        PionPairVisibleEventEvaluation(HadronicTauPairSolutionStatus::non_unique, {})};
    assert(tauamp::delphes::make_pion_pair_event_loop_record(9, non_unique).selection_status ==
           PionPairSelectionStatus::kinematic_non_unique);

    const tauamp::delphes::PionPairReaderConfig config{0.0, 0.0, 1.77686,
                                                        tauamp::delphes::PionObjectSource::generated};
    const tauamp::delphes::PionPairDelphesReader reader(argv[1], config);
    assert(reader.entry_count() > 0);
    const tauamp::delphes::PionPairDelphesEventEvaluator evaluator(
        tauamp::tautau::ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952});
    const tauamp::delphes::PionPairDelphesEventLoop loop(reader, evaluator);

    std::vector<PionPairEventLoopRecord> records;
    const auto summary = loop.for_each(0, 1, [&](const PionPairEventLoopRecord& record) { records.push_back(record); });
    assert(records.size() == 1U);
    assert(records.front().entry == 0);
    assert(records.front().classification == TauPairDecayClassification::supported_pion_pair);
    assert(records.front().tau_momentum_provenance == TauMomentumProvenance::generator_truth_available);
    assert(records.front().observation.has_value());
    assert(records.front().observation->beams.electron_lab().pz() < 0.0);
    assert(records.front().evaluation.has_value());
    assert(summary.entries_visited == 1);
    assert(summary.supported_input == 1);
    assert(summary.unsupported_input == 0);
    assert(summary.kinematic_no_solution + summary.kinematic_non_unique + summary.evaluable_one_solution +
               summary.evaluable_two_solutions ==
           1);

    bool rejected_invalid_range = false;
    try {
        loop.for_each(1, 0, [](const PionPairEventLoopRecord&) {});
    } catch (const std::out_of_range&) {
        rejected_invalid_range = true;
    }
    assert(rejected_invalid_range);
}
