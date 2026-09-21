#include <cassert>
#include <stdexcept>

#include "tauamp/tautau/common_channel_export.h"

int main() {
    using namespace tauamp::tautau;
    CommonChannelMetadata metadata;
    metadata.momentum_frame = "lab";
    metadata.object_source = "reconstructed";
    metadata.branch_weight_rule = "candidate_then_tauoo_branch_marginalization_v1";
    metadata.ordered_channel = OrderedChannel::rho_a1;
    metadata.tau_plus_decay = "a1+ nu~";
    metadata.tau_minus_decay = "rho- nu";
    metadata.decay_model_revision = "g6-rho-a1-primary-v1";
    metadata.current_revision = "g6-rho-a1-current-v1";
    metadata.lineshape_revision = "g6-rho-kw-constant-width-v1";
    metadata.phase_space_measure_revision = "g6-lorentz-invariant-four-body-massless-neutrino-v1";
    metadata.current_normalization_revision = "g6-a1-reduced-normalization-GKW1-Na1-1-v1";
    metadata.identical_slot_rule = "descending-reconstructed-pT;stable-object-index-tiebreak";
    metadata.charge_convention_revision = "g6-charge-labelled-contraction-v1";
    metadata.analyser_convention_revision = "g6-p43-p44-primary-v1";
    metadata.pi0_object_source = "NeutralPion";
    metadata.selection_revision = "g6-rho-a1-candidate-v1";
    metadata.tauoo_revision = "fixture";
    metadata.exporter_revision = "fixture";
    metadata.source_hashes = {{"fixture", std::string(64, '0')}};
    metadata.cross_feed_policy = "retain-overlap-record-v1";
    metadata.failure_record = "none";

    CommonHypothesisRecord hypothesis;
    hypothesis.valid = true;
    hypothesis.weight = 1.0;
    CommonEventRecord event{101, OrderedChannel::rho_a1, "selected::rho_a1", "rho_a1", 2U, {hypothesis}};
    validate_common_event_record(metadata, event);

    metadata.truth_tau_input = true;
    bool rejected = false;
    try {
        validate_common_event_record(metadata, event);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
}
