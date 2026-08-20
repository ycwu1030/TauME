#include <cassert>
#include <cmath>
#include <string>

#include "reference_fixture.h"
#include "tauamp/tautau/matrix_element.h"
#include "tauamp/tautau/pion_pair_marginalization.h"

namespace {
constexpr double kTolerance = 3e-11;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

void expect_equal(const tauamp::tautau::LinearComponents& left, const tauamp::tautau::LinearComponents& right) {
    assert(close(left.sm, right.sm));
    assert(close(left.f2_real, right.f2_real));
    assert(close(left.f2_imaginary, right.f2_imaginary));
    assert(close(left.f3_real, right.f3_real));
    assert(close(left.f3_imaginary, right.f3_imaginary));
}

void expect_equal(const tauamp::tautau::LinearObservableComponents& left,
                  const tauamp::tautau::LinearObservableComponents& right) {
    assert(close(left.f2_real, right.f2_real));
    assert(close(left.f2_imaginary, right.f2_imaginary));
    assert(close(left.f3_real, right.f3_real));
    assert(close(left.f3_imaginary, right.f3_imaginary));
}
}  // namespace

int main(int argc, char* argv[]) {
    using tauamp::tautau::BeamState;
    using tauamp::tautau::ElectroweakParameters;
    using tauamp::tautau::PionPairMatrixElement;
    using tauamp::tautau::TauPairKinematicPoint;

    const std::string fixture_directory = argc == 2 ? argv[1] : "tests/fixtures/tautau";
    for (const auto& fixture : tautau_test::read_pion_pair_reference_cases(fixture_directory + "/pion_pair_reference_cases.dat")) {
        const BeamState beams(fixture.electron, fixture.positron, fixture.electron_polarization,
                              fixture.positron_polarization);
        const TauPairKinematicPoint point(beams, fixture.tau_minus, fixture.tau_plus, fixture.tau_mass);
        const PionPairMatrixElement matrix_element(fixture.parameters);
        expect_equal(matrix_element.components(point, fixture.pion_minus, fixture.pion_plus), fixture.expected);
    }

    for (const auto& fixture : tautau_test::read_pion_pair_ensemble_cases(fixture_directory + "/pion_pair_ensemble_cases.dat")) {
        expect_equal(tauamp::tautau::average_conditional_observables(fixture.components), fixture.expected_average_of_ratios);
        expect_equal(tauamp::tautau::ratio_of_average_components(fixture.components), fixture.expected_ratio_after_average);
    }
}
