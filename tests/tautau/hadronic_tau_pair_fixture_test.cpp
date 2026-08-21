#include <cassert>
#include <cmath>
#include <string>

#include "hadronic_tau_pair_fixture.h"
#include "tauamp/tautau/hadronic_tau_pair_kinematics.h"

namespace {
constexpr double kTolerance = 5e-10;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

void expect_equal(const tauamp::tautau::FourMomentum& left, const tauamp::tautau::FourMomentum& right) {
    assert(close(left.px(), right.px()));
    assert(close(left.py(), right.py()));
    assert(close(left.pz(), right.pz()));
    assert(close(left.energy(), right.energy()));
}
}  // namespace

int main(int argc, char* argv[]) {
    const std::string fixture_path = argc == 2 ? argv[1] : "tests/fixtures/tautau/hadronic_tau_pair_twofold_cases.dat";
    for (const auto& fixture : tautau_test::read_hadronic_tau_pair_cases(fixture_path)) {
        const auto result = tauamp::tautau::HadronicTauPairKinematicSolver::solve(fixture.observation);
        assert(result.status == fixture.status);
        assert(result.solutions.size() == fixture.solutions.size());
        for (std::size_t index = 0; index < result.solutions.size(); ++index) {
            expect_equal(result.solutions[index].point.tau_minus_lab(), fixture.solutions[index].point.tau_minus_lab());
            expect_equal(result.solutions[index].point.tau_plus_lab(), fixture.solutions[index].point.tau_plus_lab());
            expect_equal(result.solutions[index].neutrino_minus_lab, fixture.solutions[index].neutrino_minus_lab);
            expect_equal(result.solutions[index].neutrino_plus_lab, fixture.solutions[index].neutrino_plus_lab);
        }
    }
}
