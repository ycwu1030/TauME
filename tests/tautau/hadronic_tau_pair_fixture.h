#ifndef TAUAMP_TAUTAU_TEST_HADRONIC_TAU_PAIR_FIXTURE_H_
#define TAUAMP_TAUTAU_TEST_HADRONIC_TAU_PAIR_FIXTURE_H_

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"

namespace tautau_test {

struct HadronicTauPairFixtureCase {
    std::string name;
    tauamp::tautau::HadronicTauPairObservation observation{tauamp::tautau::BeamState(
        tauamp::tautau::FourMomentum(0.0, 0.0, 1.0, 1.0), tauamp::tautau::FourMomentum(0.0, 0.0, -1.0, 1.0), 0.0, 0.0),
                                                               {}, {}, 1.0};
    tauamp::tautau::HadronicTauPairSolutionStatus status{};
    std::vector<tauamp::tautau::HadronicTauPairKinematicSolution> solutions;
};

inline void require_hadronic_fixture(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

inline tauamp::tautau::FourMomentum read_hadronic_four_momentum(std::istream& input, const std::string& field) {
    double px{};
    double py{};
    double pz{};
    double energy{};
    require_hadronic_fixture(static_cast<bool>(input >> px >> py >> pz >> energy), "invalid four-momentum in " + field);
    return {px, py, pz, energy};
}

inline tauamp::tautau::HadronicTauPairSolutionStatus read_hadronic_status(const std::string& value) {
    using Status = tauamp::tautau::HadronicTauPairSolutionStatus;
    if (value == "no_solution") return Status::no_solution;
    if (value == "one_solution") return Status::one_solution;
    if (value == "two_solutions") return Status::two_solutions;
    if (value == "non_unique") return Status::non_unique;
    throw std::runtime_error("unknown hadronic fixture status: " + value);
}

inline std::vector<HadronicTauPairFixtureCase> read_hadronic_tau_pair_cases(const std::string& path) {
    using namespace tauamp::tautau;
    std::ifstream input(path);
    require_hadronic_fixture(static_cast<bool>(input), "cannot open hadronic twofold fixture: " + path);

    std::vector<HadronicTauPairFixtureCase> result;
    std::string field;
    std::string name;
    bool in_case = false;
    FourMomentum electron;
    FourMomentum positron;
    FourMomentum visible_minus;
    FourMomentum visible_plus;
    double electron_polarization{};
    double positron_polarization{};
    double tau_mass{};
    HadronicTauPairSolutionStatus status{};
    std::vector<HadronicTauPairKinematicSolution> solutions;

    while (input >> field) {
        if (field == "case") {
            require_hadronic_fixture(!in_case, "nested hadronic fixture case");
            require_hadronic_fixture(static_cast<bool>(input >> name), "missing hadronic fixture case name");
            solutions.clear();
            in_case = true;
        } else if (field == "electron") {
            require_hadronic_fixture(in_case, "electron outside hadronic fixture case");
            electron = read_hadronic_four_momentum(input, field);
        } else if (field == "positron") {
            require_hadronic_fixture(in_case, "positron outside hadronic fixture case");
            positron = read_hadronic_four_momentum(input, field);
        } else if (field == "polarizations") {
            require_hadronic_fixture(in_case, "polarizations outside hadronic fixture case");
            require_hadronic_fixture(static_cast<bool>(input >> electron_polarization >> positron_polarization),
                                     "invalid hadronic fixture polarizations");
        } else if (field == "tau_mass") {
            require_hadronic_fixture(in_case, "tau_mass outside hadronic fixture case");
            require_hadronic_fixture(static_cast<bool>(input >> tau_mass), "invalid hadronic fixture tau mass");
        } else if (field == "visible_minus") {
            require_hadronic_fixture(in_case, "visible_minus outside hadronic fixture case");
            visible_minus = read_hadronic_four_momentum(input, field);
        } else if (field == "visible_plus") {
            require_hadronic_fixture(in_case, "visible_plus outside hadronic fixture case");
            visible_plus = read_hadronic_four_momentum(input, field);
        } else if (field == "status") {
            require_hadronic_fixture(in_case, "status outside hadronic fixture case");
            std::string value;
            require_hadronic_fixture(static_cast<bool>(input >> value), "missing hadronic fixture status");
            status = read_hadronic_status(value);
        } else if (field == "solution") {
            require_hadronic_fixture(in_case, "solution outside hadronic fixture case");
            const FourMomentum tau_minus = read_hadronic_four_momentum(input, field);
            const FourMomentum tau_plus = read_hadronic_four_momentum(input, field);
            const FourMomentum neutrino_minus = read_hadronic_four_momentum(input, field);
            const FourMomentum neutrino_plus = read_hadronic_four_momentum(input, field);
            const BeamState beams(electron, positron, electron_polarization, positron_polarization);
            solutions.push_back({TauPairKinematicPoint(beams, tau_minus, tau_plus, tau_mass), neutrino_minus, neutrino_plus});
        } else if (field == "end") {
            require_hadronic_fixture(in_case, "end outside hadronic fixture case");
            const BeamState beams(electron, positron, electron_polarization, positron_polarization);
            result.push_back({name, {beams, visible_minus, visible_plus, tau_mass}, status, std::move(solutions)});
            in_case = false;
        } else {
            throw std::runtime_error("unknown hadronic fixture field: " + field);
        }
    }
    require_hadronic_fixture(!in_case, "unterminated hadronic fixture case");
    require_hadronic_fixture(!result.empty(), "hadronic twofold fixture has no cases");
    return result;
}

}  // namespace tautau_test

#endif  // TAUAMP_TAUTAU_TEST_HADRONIC_TAU_PAIR_FIXTURE_H_
