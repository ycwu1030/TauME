#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tauamp/delphes/common_channel_root_writer.h"
#include "tauamp/delphes/rho_family_event_loop_v2.h"

namespace {

using tauamp::delphes::g7_rho_channels;
using tauamp::tautau::OrderedChannel;

struct Configuration {
    std::string input;
    std::string output_dir;
    OrderedChannel truth_source;
    std::string input_sha256;
    std::string implementation_sha256;
    double width_multiplier{3.0};
    long long entry_begin{0};
    long long entry_end{-1};
    long long event_id_offset{0};
};

const std::string& value(const std::map<std::string, std::string>& values, const char* name) {
    const auto found = values.find(name);
    if (found == values.end()) throw std::invalid_argument(std::string("missing option ") + name);
    return found->second;
}

OrderedChannel parse_channel(const std::string& name) {
    for (const auto channel : g7_rho_channels)
        if (name == tauamp::tautau::ordered_channel_name(channel)) return channel;
    throw std::invalid_argument("--truth-source must be pi_rho, rho_pi, or rho_rho");
}

void validate_hash(const std::string& hash, const char* name) {
    if (hash.size() != 64U) throw std::invalid_argument(std::string(name) + " must be a 64-character SHA-256");
    for (const unsigned char character : hash)
        if (std::isxdigit(character) == 0) throw std::invalid_argument(std::string(name) + " must be hexadecimal");
}

Configuration parse_configuration(int argc, char* argv[]) {
    if (argc < 2 || argc % 2 == 0) throw std::invalid_argument("options must be supplied as name/value pairs");
    std::map<std::string, std::string> values;
    for (int index = 1; index < argc; index += 2) {
        if (!values.emplace(argv[index], argv[index + 1]).second)
            throw std::invalid_argument("option supplied more than once");
    }
    constexpr const char* required[] = {
        "--input", "--output-dir", "--truth-source", "--input-sha256", "--implementation-sha256"};
    for (const char* option : required) (void)value(values, option);
    constexpr const char* optional[] = {
        "--width-multiplier", "--entry-begin", "--entry-end", "--event-id-offset"};
    for (const auto& [name, ignored] : values) {
        bool known = false;
        for (const char* option : required) known = known || name == option;
        for (const char* option : optional) known = known || name == option;
        if (!known) throw std::invalid_argument("unknown option " + name);
    }
    Configuration result{value(values, "--input"), value(values, "--output-dir"),
                          parse_channel(value(values, "--truth-source")),
                          value(values, "--input-sha256"), value(values, "--implementation-sha256"),
                          values.count("--width-multiplier") ? std::stod(values.at("--width-multiplier")) : 3.0,
                          values.count("--entry-begin") ? std::stoll(values.at("--entry-begin")) : 0,
                          values.count("--entry-end") ? std::stoll(values.at("--entry-end")) : -1,
                          values.count("--event-id-offset") ? std::stoll(values.at("--event-id-offset")) : 0};
    validate_hash(result.input_sha256, "--input-sha256");
    validate_hash(result.implementation_sha256, "--implementation-sha256");
    if (!(result.width_multiplier > 0.0)) throw std::invalid_argument("--width-multiplier must be positive");
    if (result.entry_begin < 0) throw std::invalid_argument("--entry-begin must be non-negative");
    if (result.entry_end >= 0 && result.entry_end < result.entry_begin)
        throw std::invalid_argument("--entry-end must not precede --entry-begin");
    if (result.event_id_offset < 0) throw std::invalid_argument("--event-id-offset must be non-negative");
    return result;
}

std::string decay_name(OrderedChannel channel, bool tau_plus) {
    switch (channel) {
        case OrderedChannel::pi_rho: return tau_plus ? "pi+ nu~" : "rho- -> pi- pi0 nu";
        case OrderedChannel::rho_pi: return tau_plus ? "rho+ -> pi+ pi0 nu~" : "pi- nu";
        case OrderedChannel::rho_rho: return tau_plus ? "rho+ -> pi+ pi0 nu~" : "rho- -> pi- pi0 nu";
        default: throw std::invalid_argument("invalid G7 rho channel");
    }
}

tauamp::tautau::CommonChannelMetadata metadata(
    OrderedChannel channel, const Configuration& configuration) {
    tauamp::tautau::CommonChannelMetadata result;
    result.ordered_channel = channel;
    result.momentum_frame = "lab";
    result.object_source = "reconstructed";
    result.branch_weight_rule = "candidate_then_tauoo_branch_marginalization_v2";
    result.tau_plus_decay = decay_name(channel, true);
    result.tau_minus_decay = decay_name(channel, false);
    result.decay_model_revision = "g6-rho-a1-primary-v1";
    result.current_revision = "g6-rho-transverse-rt-v1";
    result.lineshape_revision = "g6-rho-kw-constant-width-v1";
    result.phase_space_measure_revision = "g6-lorentz-invariant-two-body-v1";
    result.current_normalization_revision = "g6-rho-reduced-normalization-v1";
    result.identical_slot_rule = "not_applicable";
    result.charge_convention_revision = "g6-charge-labelled-contraction-v1";
    result.analyser_convention_revision = "g6-p43-p44-primary-v1";
    result.pi0_object_source = "NeutralPion";
    result.selection_revision = "g7-rho-mass-score-unique-target-v2-3gamma";
    result.tauoo_revision = "48289fea12b4198ddc1c230048ad72919215d29c";
    result.exporter_revision = "g7-rho-family-event-loop-v2-chunked-v1";
    result.source_hashes = {{"input", configuration.input_sha256}, {"implementation", configuration.implementation_sha256}};
    result.cross_feed_policy = "g7-rho-unique-target-mass-score-deltar-v2";
    result.failure_record = "none";
    return result;
}


}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Configuration configuration = parse_configuration(argc, argv);
        const std::filesystem::path output_dir(configuration.output_dir);
        if (std::filesystem::exists(output_dir)) throw std::invalid_argument("output directory already exists");
        std::filesystem::create_directories(output_dir);

        const tauamp::tautau::FourMomentum electron(0.0, 0.0, -2.13, 2.13);
        const tauamp::tautau::FourMomentum positron(0.0, 0.0, 2.13, 2.13);
        const tauamp::tautau::BeamState beams(electron, positron, 0.0, 0.0);
        const tauamp::tautau::HadronicPairMatrixElement matrix(
            tauamp::tautau::ElectroweakParameters{0.313, 0.48, 91.1876, 2.4952},
            tauamp::tautau::ProductionBosons::photon_only);
        const tauamp::delphes::CommonChannelDelphesReader reader(
            configuration.input, {electron, positron});
        const Long64_t input_entries = reader.entry_count();
        const Long64_t entry_begin = static_cast<Long64_t>(configuration.entry_begin);
        const Long64_t entry_end = configuration.entry_end < 0
            ? input_entries
            : static_cast<Long64_t>(configuration.entry_end);
        if (entry_begin > input_entries || entry_end > input_entries)
            throw std::invalid_argument("requested entry range exceeds input entry count");

        std::array<std::unique_ptr<tauamp::delphes::CommonChannelRootWriter>, 3> writers;
        for (std::size_t index = 0; index < g7_rho_channels.size(); ++index) {
            const auto channel_name = std::string(tauamp::tautau::ordered_channel_name(g7_rho_channels[index]));
            writers[index] = std::make_unique<tauamp::delphes::CommonChannelRootWriter>(
                (output_dir / (std::string("rho_family_") + channel_name + ".root")).string(),
                metadata(g7_rho_channels[index], configuration));
        }
        std::array<long long, 3> candidate_events{};
        std::array<long long, 3> evaluable_events{};
        std::array<std::array<long long, 3>, 3> matrix_counts{};
        std::map<std::string, long long> category_counts;
        std::map<long long, std::string> event_categories;
        long long technical_events = 0;
        long long pre_resolution_multi_events = 0;

        for (Long64_t entry = entry_begin; entry < entry_end; ++entry) {
            const long long source_event_id = configuration.event_id_offset + entry;
            const auto collections = reader.read(entry);
            const auto evaluated = tauamp::delphes::evaluate_rho_family_event_v2(
                collections, beams, 1.77686, matrix, configuration.width_multiplier);
            for (std::size_t index = 0; index < g7_rho_channels.size(); ++index) {
                candidate_events[index] += evaluated.targets[index].diagnostic.candidate_count != 0U;
                evaluable_events[index] += evaluated.targets[index].diagnostic.valid_candidate_count != 0U;
            }
            if (evaluated.pre_resolution_multi_target) ++pre_resolution_multi_events;
            std::string decision_category;
            std::optional<std::size_t> target_index = evaluated.winner;
            bool event_technical = false;
            bool event_candidate = false;
            for (const auto& target : evaluated.targets) event_technical = event_technical || target.diagnostic.technical_failure;
            for (const auto& target : evaluated.targets) event_candidate = event_candidate || target.diagnostic.candidate_count != 0U;
            if (event_technical) ++technical_events;
            if (target_index.has_value()) {
                const auto target = g7_rho_channels[target_index.value()];
                decision_category = target == configuration.truth_source
                    ? std::string("selected::") + tauamp::tautau::ordered_channel_name(target)
                    : std::string("cross_feed::") + tauamp::tautau::ordered_channel_name(configuration.truth_source) + "::" + tauamp::tautau::ordered_channel_name(target);
                ++matrix_counts[tauamp::delphes::g7_rho_channel_index(configuration.truth_source)][target_index.value()];
                if (!evaluated.targets[target_index.value()].evaluation.has_value()) throw std::runtime_error("selected v2 target evaluation was not retained");
                auto record = tauamp::tautau::make_ordered_common_event_record(
                    source_event_id, target, evaluated.targets[target_index.value()].evaluation.value());
                record.category_code = decision_category;
                record.truth_ordered_channel = tauamp::tautau::ordered_channel_name(configuration.truth_source);
                writers[target_index.value()]->write(record);
            } else if (event_technical) {
                decision_category = "technical_failure";
            } else if (!event_candidate) {
                decision_category = "rejected::no_ownership_candidate";
            } else {
                decision_category = "unsupported_or_no_solution";
            }
            ++category_counts[decision_category];
            event_categories[source_event_id] = decision_category;
        }
        for (auto& writer : writers) writer->finalize();

        std::ofstream summary(output_dir / "rho_family_cross_feed_summary.json");
        if (!summary) throw std::runtime_error("cannot create G7 rho summary");
        summary << "{\n  \"result_identity\": \"detector_event_accounting_selection_scan\",\n";
        summary << "  \"truth_source_channel\": \"" << tauamp::tautau::ordered_channel_name(configuration.truth_source) << "\",\n";
        summary << "  \"input_entries_total\": " << input_entries << ",\n";
        summary << "  \"entry_begin\": " << entry_begin << ",\n";
        summary << "  \"entry_end\": " << entry_end << ",\n";
        summary << "  \"event_id_offset\": " << configuration.event_id_offset << ",\n";
        summary << "  \"generated_or_delphes_entries\": " << (entry_end - entry_begin) << ",\n";
        summary << "  \"technical_failure_events\": " << technical_events << ",\n";
        summary << "  \"pre_resolution_multi_target_events\": " << pre_resolution_multi_events << ",\n";
        summary << "  \"width_multiplier\": " << configuration.width_multiplier << ",\n";
        summary << "  \"candidate_events\": [" << candidate_events[0] << ", " << candidate_events[1] << ", " << candidate_events[2] << "],\n";
        summary << "  \"evaluable_events\": [" << evaluable_events[0] << ", " << evaluable_events[1] << ", " << evaluable_events[2] << "],\n";
        summary << "  \"target_order\": [\"pi_rho\", \"rho_pi\", \"rho_rho\"],\n  \"matrix\": [\n";
        for (std::size_t row = 0; row < 3; ++row) {
            summary << "    [" << matrix_counts[row][0] << ", " << matrix_counts[row][1] << ", " << matrix_counts[row][2] << "]" << (row + 1 == 3 ? "\n" : ",\n");
        }
        summary << "  ],\n  \"category_counts\": {\n";
        std::size_t category_index = 0;
        for (const auto& [category, count] : category_counts) {
            summary << "    \"" << category << "\": " << count << (++category_index == category_counts.size() ? "\n" : ",\n");
        }
        summary << "  },\n  \"event_categories\": {\n";
        std::size_t event_index = 0;
        for (const auto& [event_id, category] : event_categories) {
            summary << "    \"" << event_id << "\": \"" << category << "\"" << (++event_index == event_categories.size() ? "\n" : ",\n");
        }
        summary << "  }\n}\n";
        std::cout << "entries_visited=" << (entry_end - entry_begin) << "\n";
        std::cout << "summary=" << (output_dir / "rho_family_cross_feed_summary.json") << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rho_family_delphes_export: " << error.what() << '\n';
        return 1;
    }
}
