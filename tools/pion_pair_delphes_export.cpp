#include <cmath>
#include <cctype>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

#include "tauamp/delphes/pion_pair_event_loop.h"
#include "tauamp/delphes/pion_pair_root_writer.h"

namespace {

struct Configuration {
    std::string input_path;
    std::string output_path;
    Long64_t begin;
    Long64_t end;
    tauamp::delphes::PionObjectSource pion_source;
    double electron_polarization;
    double positron_polarization;
    double tau_mass;
    tauamp::tautau::ElectroweakParameters electroweak_parameters;
    std::string input_sha256;
    std::string implementation_revision;
};

const std::string& required_value(const std::map<std::string, std::string>& values, const char* name) {
    const auto found = values.find(name);
    if (found == values.end()) throw std::invalid_argument(std::string("missing required option ") + name);
    return found->second;
}

Long64_t parse_entry(const std::string& value, const char* name) {
    std::size_t consumed = 0;
    long long parsed{};
    try {
        parsed = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string(name) + " must be an integer");
    }
    if (consumed != value.size() || parsed < std::numeric_limits<Long64_t>::min() ||
        parsed > std::numeric_limits<Long64_t>::max())
        throw std::invalid_argument(std::string(name) + " must be an integer");
    return static_cast<Long64_t>(parsed);
}

double parse_finite_number(const std::string& value, const char* name) {
    std::size_t consumed = 0;
    double parsed{};
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string(name) + " must be a finite number");
    }
    if (consumed != value.size() || !std::isfinite(parsed))
        throw std::invalid_argument(std::string(name) + " must be a finite number");
    return parsed;
}

tauamp::delphes::PionObjectSource parse_pion_source(const std::string& value) {
    if (value == "generated") return tauamp::delphes::PionObjectSource::generated;
    if (value == "reconstructed") return tauamp::delphes::PionObjectSource::reconstructed;
    throw std::invalid_argument("--pion-source must be generated or reconstructed");
}

void validate_sha256(const std::string& value) {
    if (value.size() != 64U)
        throw std::invalid_argument("--input-sha256 must contain exactly 64 hexadecimal characters");
    for (const unsigned char character : value) {
        if (std::isxdigit(character) == 0)
            throw std::invalid_argument("--input-sha256 must contain exactly 64 hexadecimal characters");
    }
}

Configuration parse_configuration(int argc, char* argv[]) {
    if (argc < 2 || (argc % 2) == 0)
        throw std::invalid_argument("options must be supplied as exactly one value for each required flag");

    std::map<std::string, std::string> values;
    for (int index = 1; index < argc; index += 2) {
        const std::string option = argv[index];
        const auto inserted = values.emplace(option, argv[index + 1]);
        if (!inserted.second) throw std::invalid_argument("option supplied more than once: " + option);
    }

    constexpr const char* required_options[] = {
        "--input",          "--output",              "--begin",          "--end",     "--pion-source",
        "--electron-polarization", "--positron-polarization", "--tau-mass",       "--electric-charge",
        "--weak-mixing-sine",      "--z-mass",                 "--z-width",       "--input-sha256",
        "--implementation-revision"};
    for (const char* option : required_options) (void)required_value(values, option);
    if (values.size() != std::size(required_options)) throw std::invalid_argument("unknown command-line option");

    Configuration configuration{
        required_value(values, "--input"),
        required_value(values, "--output"),
        parse_entry(required_value(values, "--begin"), "--begin"),
        parse_entry(required_value(values, "--end"), "--end"),
        parse_pion_source(required_value(values, "--pion-source")),
        parse_finite_number(required_value(values, "--electron-polarization"), "--electron-polarization"),
        parse_finite_number(required_value(values, "--positron-polarization"), "--positron-polarization"),
        parse_finite_number(required_value(values, "--tau-mass"), "--tau-mass"),
        {parse_finite_number(required_value(values, "--electric-charge"), "--electric-charge"),
         parse_finite_number(required_value(values, "--weak-mixing-sine"), "--weak-mixing-sine"),
         parse_finite_number(required_value(values, "--z-mass"), "--z-mass"),
         parse_finite_number(required_value(values, "--z-width"), "--z-width")},
        required_value(values, "--input-sha256"),
        required_value(values, "--implementation-revision")};
    validate_sha256(configuration.input_sha256);
    return configuration;
}

void validate_range(const Configuration& configuration, Long64_t entry_count) {
    if (configuration.begin < 0 || configuration.end < configuration.begin || configuration.end > entry_count)
        throw std::out_of_range("Delphes event-loop range is outside reader entries");
}

void print_summary(const tauamp::delphes::PionPairEventLoopSummary& summary) {
    std::cout << "entries_visited=" << summary.entries_visited << '\n';
    std::cout << "supported_input=" << summary.supported_input << '\n';
    std::cout << "unsupported_input=" << summary.unsupported_input << '\n';
    std::cout << "kinematic_no_solution=" << summary.kinematic_no_solution << '\n';
    std::cout << "kinematic_non_unique=" << summary.kinematic_non_unique << '\n';
    std::cout << "evaluable_one_solution=" << summary.evaluable_one_solution << '\n';
    std::cout << "evaluable_two_solutions=" << summary.evaluable_two_solutions << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Configuration configuration = parse_configuration(argc, argv);
        const tauamp::delphes::PionPairReaderConfig reader_configuration{
            configuration.electron_polarization, configuration.positron_polarization, configuration.tau_mass,
            configuration.pion_source};
        const tauamp::delphes::PionPairDelphesReader reader(configuration.input_path, reader_configuration);
        validate_range(configuration, reader.entry_count());

        const tauamp::delphes::PionPairDelphesEventEvaluator evaluator(configuration.electroweak_parameters);
        const tauamp::delphes::PionPairDelphesEventLoop loop(reader, evaluator);
        const tauamp::delphes::PionPairRootOutputMetadata metadata{
            configuration.input_sha256,
            "Delphes",
            configuration.pion_source,
            configuration.electron_polarization,
            configuration.positron_polarization,
            configuration.tau_mass,
            configuration.electroweak_parameters,
            configuration.implementation_revision};
        tauamp::delphes::PionPairRootWriter writer(configuration.output_path, metadata);
        const tauamp::delphes::PionPairEventLoopSummary summary = loop.for_each(
            configuration.begin, configuration.end,
            [&writer](const tauamp::delphes::PionPairEventLoopRecord& record) { writer.write(record); });
        writer.finalize();
        print_summary(summary);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "pion_pair_delphes_export: " << error.what() << '\n';
        return 1;
    }
}
