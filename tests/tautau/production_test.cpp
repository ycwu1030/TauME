#include <array>
#include <cassert>
#include <cmath>

#include "tauamp/tautau/production.h"

namespace {
constexpr double kTolerance = 3e-11;

bool close(double left, double right) {
    return std::abs(left - right) < kTolerance;
}

void assert_matrix_close(const tauamp::tautau::PauliBasisMatrix& matrix,
                         const std::array<double, 16>& expected) {
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            assert(close(matrix(row, column), expected.at(4 * row + column)));
        }
    }
}

constexpr std::array<double, 16> kSm = {
    0.020141724142036844, 0.010001549243338076, 1.5303951932885371e-07, 0.0091837045472659175,
    0.010001549243338076, 0.011823717055622626, 9.0311600780165996e-08, 0.0090591275539739398,
    1.5303951932885371e-07, 9.0311600780166843e-08, -0.0019574364278415546, 5.9403462832338488e-08,
    0.0091837045472659209, 0.0090591275539739433, 5.9403462832338488e-08, 0.010275443514255765,
};
constexpr std::array<double, 16> kF2Real = {
    0.044207019995651839, 0.023971846086302824, 2.1029740955117085e-07, 0.018373173642305243,
    0.023971846086302831, 0.027561341928923639, 9.0311600780166843e-08, 0.021717674559940128,
    2.1029740955117085e-07, 9.0311600780166843e-08, -3.4694469519536142e-18, 7.4041968013015758e-08,
    0.018373173642305271, 0.021717674559940121, 7.4041968013014911e-08, 0.016645678066728196,
};
constexpr std::array<double, 16> kF2Imaginary = {
    -2.3858701516180081e-07, -2.9019061643237309e-08, 0.0035888944614727645, -1.5810114267419961e-07,
    -2.9019061643237309e-08, 2.6568893864500898e-08, -3.2934544673904805e-06, -1.9282142323800988e-07,
    0.0035888944614727645, -3.2934544673913508e-06, -4.3368086899420177e-19, 0.0039660565386280032,
    -1.5810114267419961e-07, -1.9282142323800988e-07, 0.0039660565386280041, -2.6515590902803643e-07,
};
constexpr std::array<double, 16> kF3Real = {
    0.0, -1.0742842704476729e-07, 0.0074452610242628217, -7.0804295013854324e-09,
    1.0742842704476729e-07, 0.0, 0.007344910160030446, 2.74650634127932e-08,
    -0.0074452610242628217, -0.007344910160030446, 0.0, -0.0067440285164716693,
    7.0804295013854324e-09, -2.74650634127932e-08, 0.0067440285164716693, 0.0,
};
constexpr std::array<double, 16> kF3Imaginary = {
    0.0, -0.0067440285164716719, -2.74650634127878e-08, 0.0073449101600304469,
    0.0067440285164716719, 0.0, 7.0804295011648862e-09, 0.0074452610242628191,
    2.74650634127878e-08, -7.0804295011652038e-09, 0.0, 1.0742842704420073e-07,
    -0.0073449101600304434, -0.0074452610242628208, -1.0742842704420073e-07, -1.7347234759768071e-18,
};
constexpr std::array<double, 16> kUnpolarizedSm = {
    0.017486799937152842, 1.7516105429769053e-06, 1.4033267616726074e-07, 2.2091693784754307e-06,
    1.7516105429769053e-06, 0.010265500416838554, -1.0265855465588651e-08, 0.0078649770765990626,
    1.4033267616726074e-07, -1.0265855465588651e-08, -0.0016994709338559132, -6.7524809473445982e-09,
    2.2091693784823696e-06, 0.0078649770765990556, -6.7524809473445982e-09, 0.0089207704541702099,
};
}  // namespace

int main() {
    using tauamp::tautau::FourMomentum;
    using tauamp::tautau::BeamState;
    using tauamp::tautau::ElectronPositronToTauPair;
    using tauamp::tautau::ElectroweakParameters;
    using tauamp::tautau::TauPairKinematicPoint;

    constexpr double mass = 1.77686;
    constexpr double energy = 2.1;
    constexpr double theta = 0.91;
    const double momentum = std::sqrt(energy * energy - mass * mass);
    const BeamState beams(FourMomentum(0.0, 0.0, energy, energy), FourMomentum(0.0, 0.0, -energy, energy), 0.37, -0.41);
    const TauPairKinematicPoint symmetric(
        beams, FourMomentum(momentum * std::sin(theta), 0.0, momentum * std::cos(theta), energy),
        FourMomentum(-momentum * std::sin(theta), 0.0, -momentum * std::cos(theta), energy), mass);
    const std::array<double, 3> lab_boost{{0.0, 0.0, std::tanh(0.35)}};
    const BeamState asymmetric_beams(
        FourMomentum(0.0, 0.0, energy, energy).boosted(lab_boost),
        FourMomentum(0.0, 0.0, -energy, energy).boosted(lab_boost), 0.37, -0.41);
    const TauPairKinematicPoint asymmetric(
        asymmetric_beams, symmetric.pair_cm().tau_minus.boosted(lab_boost),
        symmetric.pair_cm().tau_plus.boosted(lab_boost), mass);

    const ElectroweakParameters parameters{0.313, 0.48, 91.1876, 2.4952};
    const ElectronPositronToTauPair production(parameters);
    const auto components = production.components(symmetric);
    const auto polynomial = production.polynomial_components(symmetric);
    const auto asymmetric_components = production.components(asymmetric);

    assert_matrix_close(components.sm, kSm);
    assert_matrix_close(components.f2_real, kF2Real);
    assert_matrix_close(components.f2_imaginary, kF2Imaginary);
    assert_matrix_close(components.f3_real, kF3Real);
    assert_matrix_close(components.f3_imaginary, kF3Imaginary);
    assert_matrix_close(components.recompose({}), kSm);
    assert_matrix_close(polynomial.recompose_linear({}), kSm);
    assert_matrix_close(polynomial.sm, kSm);
    assert(!close(polynomial.f2_real_f2_real(0, 0), 0.0));
    bool has_f2_f3_cross_term = false;
    for (std::size_t row = 0; row < 4; ++row)
        for (std::size_t column = 0; column < 4; ++column)
            has_f2_f3_cross_term = has_f2_f3_cross_term || !close(polynomial.f2_real_f3_real(row, column), 0.0);
    assert(has_f2_f3_cross_term);
    assert(close(polynomial.f2_real_f2_imaginary(0, 0), 0.0));
    assert(close(polynomial.f3_real_f3_imaginary(0, 0), 0.0));

    const tauamp::tautau::PhotonDipoleFormFactors finite_form_factors{{0.002, -0.004}, {0.003, -0.001}};
    const auto full = polynomial.recompose_full(finite_form_factors);
    const auto reversed = polynomial.recompose_full(
        tauamp::tautau::PhotonDipoleFormFactors{{-0.002, 0.004}, {-0.003, 0.001}});
    const auto linear = polynomial.recompose_linear(finite_form_factors);
    const auto reversed_linear = polynomial.recompose_linear(
        tauamp::tautau::PhotonDipoleFormFactors{{-0.002, 0.004}, {-0.003, 0.001}});
    bool observed_quadratic_effect = false;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            assert(close(full(row, column) + reversed(row, column) - linear(row, column) -
                             reversed_linear(row, column),
                         2.0 * (full(row, column) - linear(row, column))));
            if (!close(full(row, column), linear(row, column))) observed_quadratic_effect = true;
        }
    }
    assert(observed_quadratic_effect);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            assert(close(components.sm(row, column), asymmetric_components.sm(row, column)));
            assert(close(components.f2_real(row, column), asymmetric_components.f2_real(row, column)));
            assert(close(components.f2_imaginary(row, column), asymmetric_components.f2_imaginary(row, column)));
            assert(close(components.f3_real(row, column), asymmetric_components.f3_real(row, column)));
            assert(close(components.f3_imaginary(row, column), asymmetric_components.f3_imaginary(row, column)));
        }
    }

    const BeamState positron_unpolarized_beams(
        FourMomentum(0.0, 0.0, energy, energy), FourMomentum(0.0, 0.0, -energy, energy), 0.37, 0.0);
    const TauPairKinematicPoint positron_unpolarized(
        positron_unpolarized_beams, symmetric.pair_cm().tau_minus, symmetric.pair_cm().tau_plus, mass);
    const auto positron_unpolarized_components = production.components(positron_unpolarized);
    assert(!close(components.sm(0, 1), positron_unpolarized_components.sm(0, 1)));

    const BeamState unpolarized_beams(
        FourMomentum(0.0, 0.0, energy, energy), FourMomentum(0.0, 0.0, -energy, energy), 0.0, 0.0);
    const TauPairKinematicPoint unpolarized(
        unpolarized_beams, symmetric.pair_cm().tau_minus, symmetric.pair_cm().tau_plus, mass);
    assert_matrix_close(production.components(unpolarized).sm, kUnpolarizedSm);
}
