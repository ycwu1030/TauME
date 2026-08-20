#ifndef TAUAMP_TAUTAU_TEST_REFERENCE_FIXTURE_H_
#define TAUAMP_TAUTAU_TEST_REFERENCE_FIXTURE_H_

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_pair_marginalization.h"

namespace tautau_test {

struct PionPairReferenceCase {
    std::string name;
    tauamp::tautau::ElectroweakParameters parameters{};
    tauamp::tautau::FourMomentum electron;
    tauamp::tautau::FourMomentum positron;
    double electron_polarization{};
    double positron_polarization{};
    double tau_mass{};
    tauamp::tautau::FourMomentum tau_minus;
    tauamp::tautau::FourMomentum tau_plus;
    tauamp::tautau::FourMomentum pion_minus;
    tauamp::tautau::FourMomentum pion_plus;
    tauamp::tautau::LinearComponents expected{};
};

struct PionPairEnsembleCase {
    std::string name;
    std::vector<tauamp::tautau::WeightedLinearComponents> components;
    tauamp::tautau::LinearObservableComponents expected_average_of_ratios{};
    tauamp::tautau::LinearObservableComponents expected_ratio_after_average{};
};

inline void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

inline tauamp::tautau::FourMomentum read_four_momentum(std::istream& input, const std::string& field) {
    double px{};
    double py{};
    double pz{};
    double energy{};
    require(static_cast<bool>(input >> px >> py >> pz >> energy), "invalid four-momentum in " + field);
    return {px, py, pz, energy};
}

inline std::vector<PionPairReferenceCase> read_pion_pair_reference_cases(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pion-pair reference fixture: " + path);

    std::vector<PionPairReferenceCase> result;
    PionPairReferenceCase current;
    bool in_case = false;
    std::string field;
    while (input >> field) {
        if (field == "case") {
            require(!in_case, "nested reference-fixture case");
            current = PionPairReferenceCase{};
            require(static_cast<bool>(input >> current.name), "missing reference-fixture case name");
            in_case = true;
        } else if (field == "parameters") {
            require(in_case, "parameters outside reference-fixture case");
            require(static_cast<bool>(input >> current.parameters.electron_charge >> current.parameters.sin_theta_w >>
                                     current.parameters.z_mass >> current.parameters.z_width),
                    "invalid reference-fixture parameters");
        } else if (field == "electron") {
            require(in_case, "electron outside reference-fixture case");
            current.electron = read_four_momentum(input, field);
        } else if (field == "positron") {
            require(in_case, "positron outside reference-fixture case");
            current.positron = read_four_momentum(input, field);
        } else if (field == "polarizations") {
            require(in_case, "polarizations outside reference-fixture case");
            require(static_cast<bool>(input >> current.electron_polarization >> current.positron_polarization),
                    "invalid reference-fixture polarizations");
        } else if (field == "tau_mass") {
            require(in_case, "tau_mass outside reference-fixture case");
            require(static_cast<bool>(input >> current.tau_mass), "invalid reference-fixture tau mass");
        } else if (field == "tau_minus") {
            require(in_case, "tau_minus outside reference-fixture case");
            current.tau_minus = read_four_momentum(input, field);
        } else if (field == "tau_plus") {
            require(in_case, "tau_plus outside reference-fixture case");
            current.tau_plus = read_four_momentum(input, field);
        } else if (field == "pion_minus") {
            require(in_case, "pion_minus outside reference-fixture case");
            current.pion_minus = read_four_momentum(input, field);
        } else if (field == "pion_plus") {
            require(in_case, "pion_plus outside reference-fixture case");
            current.pion_plus = read_four_momentum(input, field);
        } else if (field == "expected") {
            require(in_case, "expected outside reference-fixture case");
            require(static_cast<bool>(input >> current.expected.sm >> current.expected.f2_real >>
                                     current.expected.f2_imaginary >> current.expected.f3_real >>
                                     current.expected.f3_imaginary),
                    "invalid reference-fixture expected components");
        } else if (field == "end") {
            require(in_case, "end outside reference-fixture case");
            result.push_back(current);
            in_case = false;
        } else {
            throw std::runtime_error("unknown reference-fixture field: " + field);
        }
    }
    require(!in_case, "unterminated reference-fixture case");
    require(!result.empty(), "reference fixture has no cases");
    return result;
}

inline std::vector<PionPairEnsembleCase> read_pion_pair_ensemble_cases(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pion-pair ensemble fixture: " + path);

    std::vector<PionPairEnsembleCase> result;
    PionPairEnsembleCase current;
    bool in_case = false;
    std::string field;
    while (input >> field) {
        if (field == "case") {
            require(!in_case, "nested ensemble-fixture case");
            current = PionPairEnsembleCase{};
            require(static_cast<bool>(input >> current.name), "missing ensemble-fixture case name");
            in_case = true;
        } else if (field == "component") {
            require(in_case, "component outside ensemble-fixture case");
            tauamp::tautau::WeightedLinearComponents value;
            require(static_cast<bool>(input >> value.components.sm >> value.components.f2_real >>
                                     value.components.f2_imaginary >> value.components.f3_real >>
                                     value.components.f3_imaginary >> value.weight),
                    "invalid ensemble-fixture component");
            current.components.push_back(value);
        } else if (field == "expected_average_of_ratios") {
            require(in_case, "expected_average_of_ratios outside ensemble-fixture case");
            require(static_cast<bool>(input >> current.expected_average_of_ratios.f2_real >>
                                     current.expected_average_of_ratios.f2_imaginary >>
                                     current.expected_average_of_ratios.f3_real >>
                                     current.expected_average_of_ratios.f3_imaginary),
                    "invalid ensemble-fixture average-of-ratios result");
        } else if (field == "expected_ratio_after_average") {
            require(in_case, "expected_ratio_after_average outside ensemble-fixture case");
            require(static_cast<bool>(input >> current.expected_ratio_after_average.f2_real >>
                                     current.expected_ratio_after_average.f2_imaginary >>
                                     current.expected_ratio_after_average.f3_real >>
                                     current.expected_ratio_after_average.f3_imaginary),
                    "invalid ensemble-fixture ratio-after-average result");
        } else if (field == "end") {
            require(in_case, "end outside ensemble-fixture case");
            require(!current.components.empty(), "ensemble-fixture case has no components");
            result.push_back(std::move(current));
            in_case = false;
        } else {
            throw std::runtime_error("unknown ensemble-fixture field: " + field);
        }
    }
    require(!in_case, "unterminated ensemble-fixture case");
    require(!result.empty(), "ensemble fixture has no cases");
    return result;
}

}  // namespace tautau_test

#endif  // TAUAMP_TAUTAU_TEST_REFERENCE_FIXTURE_H_
