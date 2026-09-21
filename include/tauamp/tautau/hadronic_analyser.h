#ifndef TAUAMP_TAUTAU_HADRONIC_ANALYSER_H_
#define TAUAMP_TAUTAU_HADRONIC_ANALYSER_H_

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "tauamp/tautau/rho_a1_currents.h"

namespace tauamp::tautau {

struct DecayRestCoefficients {
    double rate{};
    std::array<double, 3> spin{};
};

namespace detail {
using Complex = std::complex<double>;
using Matrix4 = std::array<Complex, 16>;

inline Matrix4 matrix_zero() { return {}; }
inline Matrix4 matrix_identity() {
    Matrix4 result{};
    for (int index = 0; index < 4; ++index) result[4 * index + index] = 1.0;
    return result;
}
inline Matrix4 matrix_add(const Matrix4& left, const Matrix4& right) {
    Matrix4 result{};
    for (int index = 0; index < 16; ++index) result[index] = left[index] + right[index];
    return result;
}
inline Matrix4 matrix_subtract(const Matrix4& left, const Matrix4& right) {
    Matrix4 result{};
    for (int index = 0; index < 16; ++index) result[index] = left[index] - right[index];
    return result;
}
inline Matrix4 matrix_scale(const Matrix4& value, Complex factor) {
    Matrix4 result{};
    for (int index = 0; index < 16; ++index) result[index] = factor * value[index];
    return result;
}
inline Matrix4 matrix_product(const Matrix4& left, const Matrix4& right) {
    Matrix4 result{};
    for (int row = 0; row < 4; ++row)
        for (int column = 0; column < 4; ++column)
            for (int inner = 0; inner < 4; ++inner) result[4 * row + column] += left[4 * row + inner] * right[4 * inner + column];
    return result;
}
inline Complex matrix_trace(const Matrix4& value) {
    Complex result{};
    for (int index = 0; index < 4; ++index) result += value[4 * index + index];
    return result;
}
inline Matrix4 gamma0() {
    Matrix4 result{};
    result[0] = 1.0; result[5] = 1.0; result[10] = -1.0; result[15] = -1.0;
    return result;
}
inline Matrix4 gamma_spatial(int axis) {
    const std::array<Matrix4, 3> values{{
        Matrix4{0,0,0,1, 0,0,1,0, 0,-1,0,0, -1,0,0,0},
        Matrix4{0,0,0,Complex(0,-1), 0,0,Complex(0,1),0, 0,Complex(0,1),0,0, Complex(0,-1),0,0,0},
        Matrix4{0,0,1,0, 0,0,0,-1, -1,0,0,0, 0,1,0,0}}};
    if (axis < 0 || axis >= 3) throw std::invalid_argument("spatial gamma index must lie in [0,2]");
    return values[axis];
}
inline Matrix4 gamma5() {
    return matrix_scale(matrix_product(matrix_product(matrix_product(gamma0(), gamma_spatial(0)), gamma_spatial(1)), gamma_spatial(2)), Complex(0, 1));
}
inline Matrix4 slash(const std::array<Complex, 4>& vector) {
    Matrix4 result = matrix_scale(gamma0(), vector[0]);
    for (int axis = 0; axis < 3; ++axis) result = matrix_subtract(result, matrix_scale(gamma_spatial(axis), vector[axis + 1]));
    return result;
}
inline std::array<Complex, 4> four_vector(const FourMomentum& value) {
    return {value.energy(), value.px(), value.py(), value.pz()};
}
inline std::array<Complex, 4> four_vector(const ComplexFourVector& value) { return value; }
inline std::array<Complex, 4> conjugate_vector(const ComplexFourVector& value) {
    return {std::conj(value[0]), std::conj(value[1]), std::conj(value[2]), std::conj(value[3])};
}
inline Complex trace_rate(double mass, const FourMomentum& neutrino, const ComplexFourVector& current,
                          const std::array<double, 3>& spin, bool anti) {
    const Matrix4 identity = matrix_identity();
    const Matrix4 left = matrix_scale(matrix_subtract(identity, gamma5()), 0.5);
    const Matrix4 spin_slash = slash({0.0, spin[0], spin[1], spin[2]});
    const Matrix4 p_slash = slash({mass, 0.0, 0.0, 0.0});
    const Matrix4 projector = matrix_scale(
        matrix_product(matrix_add(p_slash, matrix_scale(identity, anti ? -mass : mass)),
                       matrix_add(identity, matrix_product(gamma5(), spin_slash))), 0.5);
    const Matrix4 current_slash = slash(four_vector(current));
    const Matrix4 conjugate_slash = slash(conjugate_vector(current));
    const Matrix4 neutrino_slash = slash(four_vector(neutrino));
    Matrix4 trace = matrix_product(neutrino_slash, current_slash);
    if (anti) {
        trace = matrix_product(projector, current_slash);
        trace = matrix_product(trace, left);
        trace = matrix_product(trace, neutrino_slash);
        trace = matrix_product(trace, conjugate_slash);
        trace = matrix_product(trace, left);
    } else {
        trace = matrix_product(trace, left);
        trace = matrix_product(trace, projector);
        trace = matrix_product(trace, conjugate_slash);
        trace = matrix_product(trace, left);
    }
    return matrix_trace(trace);
}
}  // namespace detail

inline DecayRestCoefficients decay_rest_coefficients(double tau_mass, const FourMomentum& neutrino,
                                                     const ComplexFourVector& current, bool anti = false) {
    if (!std::isfinite(tau_mass) || tau_mass <= 0.0) throw std::invalid_argument("tau mass must be finite and positive");
    const std::array<double, 3> zero{{0.0, 0.0, 0.0}};
    DecayRestCoefficients result{std::real(detail::trace_rate(tau_mass, neutrino, current, zero, anti)), {}};
    for (int index = 0; index < 3; ++index) {
        auto spin = zero;
        spin[index] = 1.0;
        result.spin[index] = std::real(detail::trace_rate(tau_mass, neutrino, current, spin, anti)) - result.rate;
    }
    if (!std::isfinite(result.rate)) throw std::runtime_error("hadronic decay trace produced a non-finite rate");
    return result;
}

inline std::array<double, 3> normalized_analyser(const DecayRestCoefficients& coefficients) {
    if (!(coefficients.rate > 0.0)) throw std::invalid_argument("hadronic analyser requires a positive decay rate");
    return {{coefficients.spin[0] / coefficients.rate, coefficients.spin[1] / coefficients.rate,
             coefficients.spin[2] / coefficients.rate}};
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_HADRONIC_ANALYSER_H_
