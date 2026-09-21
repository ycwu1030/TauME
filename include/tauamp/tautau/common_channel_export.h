#ifndef TAUAMP_TAUTAU_COMMON_CHANNEL_EXPORT_H_
#define TAUAMP_TAUTAU_COMMON_CHANNEL_EXPORT_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "tauamp/tautau/ordered_channel_adapter.h"

namespace tauamp::tautau {

inline constexpr std::size_t g6_component_count = 15U;

struct CommonChannelMetadata {
    int schema_version{2};
    std::string model_interpretation{"linear_plus_quadratic"};
    std::string quadratic_cross_term_convention{"stored_coefficient_includes_both_amplitude_orderings"};
    std::string production_bosons{"photon_only"};
    double sqrt_s_GeV{4.26};
    std::string momentum_frame;
    std::string object_source;
    std::string branch_weight_rule;
    bool truth_tau_input{false};
    bool truth_neutrino_input{false};
    OrderedChannel ordered_channel{};
    std::string tau_plus_decay;
    std::string tau_minus_decay;
    std::string decay_model_revision;
    std::string current_revision;
    std::string lineshape_revision;
    std::string phase_space_measure_revision;
    std::string current_normalization_revision;
    std::string identical_slot_rule;
    std::string charge_convention_revision;
    std::string analyser_convention_revision;
    std::string pi0_object_source;
    std::string selection_revision;
    std::string tauoo_revision;
    std::string exporter_revision;
    std::map<std::string, std::string> source_hashes;
    long long event_count{0};
    long long hypothesis_count{0};
    std::string cross_feed_policy;
    std::string failure_record;
};

struct CommonHypothesisRecord {
    std::size_t hypothesis_index{};
    double weight{};
    bool valid{};
    std::array<double, g6_component_count> components{};
};

struct CommonEventRecord {
    long long source_event_id{};
    OrderedChannel ordered_channel{};
    std::string category_code;
    std::string truth_ordered_channel;
    std::size_t candidate_hypothesis_count{};
    std::vector<CommonHypothesisRecord> hypotheses;
};

inline bool g6_sha256(const std::string& value) {
    if (value.size() != 64U) return false;
    for (const char character : value) {
        const bool digit = character >= '0' && character <= '9';
        const bool lower = character >= 'a' && character <= 'f';
        const bool upper = character >= 'A' && character <= 'F';
        if (!digit && !lower && !upper) return false;
    }
    return true;
}

inline void validate_common_channel_metadata(const CommonChannelMetadata& metadata) {
    if (metadata.schema_version != 2) throw std::invalid_argument("G6 common export requires schema version 2");
    if (metadata.model_interpretation != "linear_plus_quadratic")
        throw std::invalid_argument("G6 common export requires linear_plus_quadratic model interpretation");
    if (!std::isfinite(metadata.sqrt_s_GeV) || metadata.sqrt_s_GeV <= 0.0)
        throw std::invalid_argument("G6 common export requires a finite positive collision energy");
    if (metadata.production_bosons != "photon_only")
        throw std::invalid_argument("G6 common export requires photon_only production");
    if (metadata.object_source == "reconstructed" && (metadata.truth_tau_input || metadata.truth_neutrino_input))
        throw std::invalid_argument("reconstructed common export cannot use truth tau or neutrino input");
    const std::array<const std::string*, 18> required{{
        &metadata.momentum_frame, &metadata.object_source, &metadata.branch_weight_rule,
        &metadata.tau_plus_decay, &metadata.tau_minus_decay, &metadata.decay_model_revision,
        &metadata.current_revision, &metadata.lineshape_revision, &metadata.phase_space_measure_revision,
        &metadata.current_normalization_revision, &metadata.identical_slot_rule,
        &metadata.charge_convention_revision, &metadata.analyser_convention_revision,
        &metadata.pi0_object_source, &metadata.selection_revision, &metadata.tauoo_revision,
        &metadata.exporter_revision, &metadata.cross_feed_policy}};
    for (const auto* value : required)
        if (value->empty()) throw std::invalid_argument("G6 common metadata contains an empty required field");
    if (metadata.source_hashes.empty()) throw std::invalid_argument("G6 common metadata requires source hashes");
    for (const auto& [name, digest] : metadata.source_hashes)
        if (name.empty() || !g6_sha256(digest)) throw std::invalid_argument("G6 common metadata has an invalid source hash");
}

inline void validate_common_event_record(const CommonChannelMetadata& metadata,
                                         const CommonEventRecord& event) {
    validate_common_channel_metadata(metadata);
    if (event.ordered_channel != metadata.ordered_channel)
        throw std::invalid_argument("event ordered channel does not match metadata");
    if (event.category_code.empty()) throw std::invalid_argument("event category code must be non-empty");
    if (event.candidate_hypothesis_count == 0U)
        throw std::invalid_argument("selected common event requires at least one candidate hypothesis");
    std::size_t expected_index = 0U;
    double total_weight = 0.0;
    for (const auto& hypothesis : event.hypotheses) {
        if (!hypothesis.valid) continue;
        if (hypothesis.hypothesis_index != expected_index++)
            throw std::invalid_argument("valid common hypothesis indices must be contiguous from zero");
        if (!std::isfinite(hypothesis.weight) || hypothesis.weight < 0.0)
            throw std::invalid_argument("valid common hypothesis weight must be finite and non-negative");
        total_weight += hypothesis.weight;
        for (const double component : hypothesis.components)
            if (!std::isfinite(component)) throw std::invalid_argument("common hypothesis components must be finite");
    }
    if (expected_index == 0U || !std::isfinite(total_weight) || total_weight <= 0.0)
        throw std::invalid_argument("common event requires positive finite total valid weight");
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_COMMON_CHANNEL_EXPORT_H_
