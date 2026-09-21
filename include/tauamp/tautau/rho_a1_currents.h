#ifndef TAUAMP_TAUTAU_RHO_A1_CURRENTS_H_
#define TAUAMP_TAUTAU_RHO_A1_CURRENTS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "tauamp/tautau/four_momentum.h"

namespace tauamp::tautau {

inline constexpr double g6_rho_mass = 0.77549;
inline constexpr double g6_rho_width = 0.1491;
inline constexpr double g6_rho_prime_mass = 1.465;
inline constexpr double g6_rho_prime_width = 0.40;
inline constexpr double g6_rho_prime_beta = -0.145;
inline constexpr double g6_charged_pion_mass = 0.13957018;
inline constexpr double g6_neutral_pion_mass = 0.1349766;

using ComplexFourVector = std::array<std::complex<double>, 4>;

inline double minkowski_dot(const FourMomentum& left, const FourMomentum& right) {
    return left.energy() * right.energy() - left.px() * right.px() - left.py() * right.py() -
           left.pz() * right.pz();
}

inline std::complex<double> minkowski_dot(const FourMomentum& left, const ComplexFourVector& right) {
    return left.energy() * right[0] - left.px() * right[1] - left.py() * right[2] - left.pz() * right[3];
}

inline FourMomentum scaled(const FourMomentum& value, double factor) {
    return {value.px() * factor, value.py() * factor, value.pz() * factor, value.energy() * factor};
}

inline FourMomentum transverse(const FourMomentum& vector, const FourMomentum& total) {
    const double total_squared = total.mass_squared();
    if (!std::isfinite(total_squared) || total_squared <= 0.0)
        throw std::invalid_argument("transverse projection requires a finite time-like total momentum");
    return vector - scaled(total, minkowski_dot(total, vector) / total_squared);
}

inline ComplexFourVector complex_scaled(const FourMomentum& value, std::complex<double> factor) {
    return {factor * value.energy(), factor * value.px(), factor * value.py(), factor * value.pz()};
}

inline ComplexFourVector operator+(const ComplexFourVector& left, const ComplexFourVector& right) {
    ComplexFourVector result{};
    for (std::size_t index = 0; index < result.size(); ++index) result[index] = left[index] + right[index];
    return result;
}

inline std::complex<double> rho_kw_lineshape(double s, double mass = g6_rho_mass,
                                             double width = g6_rho_width) {
    if (!std::isfinite(s)) throw std::invalid_argument("rho invariant mass squared must be finite");
    return 1.0 / std::complex<double>(s - mass * mass, mass * width);
}

inline double g6_klambda(double a, double b) { return 1.0 + a * a + b * b - 2.0 * (a * b + b + a); }

inline std::complex<double> rho_taudecay_breit_wigner(double s, double mass, double width,
                                                       bool charged_mode = true) {
    if (!std::isfinite(s) || s < 0.0)
        throw std::invalid_argument("running-width rho invariant mass squared must be finite and non-negative");
    const double pion1 = g6_charged_pion_mass;
    const double pion2 = charged_mode ? g6_neutral_pion_mass : g6_charged_pion_mass;
    double running_width = 0.0;
    if (s > (pion1 + pion2) * (pion1 + pion2)) {
        const double root_s = std::sqrt(s);
        const double qs = std::sqrt(std::max(g6_klambda(pion1 * pion1 / s, pion2 * pion2 / s), 0.0));
        const double qm = std::sqrt(std::max(
            g6_klambda(pion1 * pion1 / (mass * mass), pion2 * pion2 / (mass * mass)), 0.0));
        running_width = width * (root_s / mass) * std::pow(qs / qm, 3);
    }
    return -mass * mass / std::complex<double>(s - mass * mass, std::sqrt(s) * running_width);
}

inline std::complex<double> rho_taudecay_lineshape(double s, bool charged_mode = true) {
    const auto rho = rho_taudecay_breit_wigner(s, g6_rho_mass, g6_rho_width, charged_mode);
    const auto rho_prime =
        rho_taudecay_breit_wigner(s, g6_rho_prime_mass, g6_rho_prime_width, charged_mode);
    return (rho + g6_rho_prime_beta * rho_prime) / (1.0 + g6_rho_prime_beta);
}

inline ComplexFourVector rho_kw_current(const FourMomentum& charged_pion,
                                        const FourMomentum& neutral_pion,
                                        bool charge_conjugate = false) {
    const FourMomentum total = charged_pion + neutral_pion;
    const FourMomentum basis = transverse(charged_pion - neutral_pion, total);
    ComplexFourVector current = complex_scaled(basis, rho_kw_lineshape(total.mass_squared()));
    if (charge_conjugate) {
        for (auto& component : current) component = std::conj(component);
    }
    return current;
}

inline ComplexFourVector a1_kw_current(const FourMomentum& same_charge_1,
                                       const FourMomentum& same_charge_2,
                                       const FourMomentum& opposite_charge,
                                       bool charge_conjugate = false) {
    const FourMomentum total = same_charge_1 + same_charge_2 + opposite_charge;
    const FourMomentum t1 = transverse(same_charge_1 - opposite_charge, total);
    const FourMomentum t2 = transverse(same_charge_2 - opposite_charge, total);
    const double s1 = (same_charge_2 + opposite_charge).mass_squared();
    const double s2 = (same_charge_1 + opposite_charge).mass_squared();
    ComplexFourVector current = complex_scaled(t1, rho_kw_lineshape(s2)) +
                                complex_scaled(t2, rho_kw_lineshape(s1));
    if (charge_conjugate) {
        for (auto& component : current) component = std::conj(component);
    }
    return current;
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_RHO_A1_CURRENTS_H_
