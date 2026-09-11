#ifndef TAUAMP_TAUTAU_SPIN_DENSITY_H_
#define TAUAMP_TAUTAU_SPIN_DENSITY_H_

#include <array>
#include <complex>
#include <cstddef>

namespace tauamp::tautau {

struct PhotonDipoleFormFactors {
    std::complex<double> f2_gamma{};
    std::complex<double> f3_gamma{};
};

class PauliBasisMatrix {
public:
    double& operator()(std::size_t row, std::size_t column) { return entries_.at(4 * row + column); }
    double operator()(std::size_t row, std::size_t column) const { return entries_.at(4 * row + column); }
    const std::array<double, 16>& entries() const { return entries_; }

private:
    std::array<double, 16> entries_{};
};

struct LinearComponents {
    double sm{};
    double f2_real{};
    double f2_imaginary{};
    double f3_real{};
    double f3_imaginary{};

    double recompose(const PhotonDipoleFormFactors& form_factors) const {
        return sm + form_factors.f2_gamma.real() * f2_real + form_factors.f2_gamma.imag() * f2_imaginary +
               form_factors.f3_gamma.real() * f3_real + form_factors.f3_gamma.imag() * f3_imaginary;
    }
};

struct SpinDensityComponents {
    PauliBasisMatrix sm;
    PauliBasisMatrix f2_real;
    PauliBasisMatrix f2_imaginary;
    PauliBasisMatrix f3_real;
    PauliBasisMatrix f3_imaginary;

    PauliBasisMatrix recompose(const PhotonDipoleFormFactors& form_factors) const {
        PauliBasisMatrix result;
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                result(row, column) = sm(row, column) + form_factors.f2_gamma.real() * f2_real(row, column) +
                                      form_factors.f2_gamma.imag() * f2_imaginary(row, column) +
                                      form_factors.f3_gamma.real() * f3_real(row, column) +
                                      form_factors.f3_gamma.imag() * f3_imaginary(row, column);
            }
        }
        return result;
    }
};

struct PolynomialComponents {
    double sm{};
    double f2_real{};
    double f2_imaginary{};
    double f3_real{};
    double f3_imaginary{};
    double f2_real_f2_real{};
    double f2_real_f2_imaginary{};
    double f2_real_f3_real{};
    double f2_real_f3_imaginary{};
    double f2_imaginary_f2_imaginary{};
    double f2_imaginary_f3_real{};
    double f2_imaginary_f3_imaginary{};
    double f3_real_f3_real{};
    double f3_real_f3_imaginary{};
    double f3_imaginary_f3_imaginary{};

    double recompose_linear(const PhotonDipoleFormFactors& form_factors) const {
        return sm + form_factors.f2_gamma.real() * f2_real + form_factors.f2_gamma.imag() * f2_imaginary +
               form_factors.f3_gamma.real() * f3_real + form_factors.f3_gamma.imag() * f3_imaginary;
    }

    double recompose_full(const PhotonDipoleFormFactors& form_factors) const {
        const double f2_real_value = form_factors.f2_gamma.real();
        const double f2_imaginary_value = form_factors.f2_gamma.imag();
        const double f3_real_value = form_factors.f3_gamma.real();
        const double f3_imaginary_value = form_factors.f3_gamma.imag();
        return recompose_linear(form_factors) + f2_real_value * f2_real_value * f2_real_f2_real +
               f2_real_value * f2_imaginary_value * f2_real_f2_imaginary +
               f2_real_value * f3_real_value * f2_real_f3_real +
               f2_real_value * f3_imaginary_value * f2_real_f3_imaginary +
               f2_imaginary_value * f2_imaginary_value * f2_imaginary_f2_imaginary +
               f2_imaginary_value * f3_real_value * f2_imaginary_f3_real +
               f2_imaginary_value * f3_imaginary_value * f2_imaginary_f3_imaginary +
               f3_real_value * f3_real_value * f3_real_f3_real +
               f3_real_value * f3_imaginary_value * f3_real_f3_imaginary +
               f3_imaginary_value * f3_imaginary_value * f3_imaginary_f3_imaginary;
    }
};

struct SpinDensityPolynomialComponents {
    PauliBasisMatrix sm;
    PauliBasisMatrix f2_real;
    PauliBasisMatrix f2_imaginary;
    PauliBasisMatrix f3_real;
    PauliBasisMatrix f3_imaginary;
    PauliBasisMatrix f2_real_f2_real;
    PauliBasisMatrix f2_real_f2_imaginary;
    PauliBasisMatrix f2_real_f3_real;
    PauliBasisMatrix f2_real_f3_imaginary;
    PauliBasisMatrix f2_imaginary_f2_imaginary;
    PauliBasisMatrix f2_imaginary_f3_real;
    PauliBasisMatrix f2_imaginary_f3_imaginary;
    PauliBasisMatrix f3_real_f3_real;
    PauliBasisMatrix f3_real_f3_imaginary;
    PauliBasisMatrix f3_imaginary_f3_imaginary;

    SpinDensityPolynomialComponents() = default;

    PauliBasisMatrix recompose_linear(const PhotonDipoleFormFactors& form_factors) const {
        return recompose(form_factors);
    }

    PauliBasisMatrix recompose(const PhotonDipoleFormFactors& form_factors) const {
        PauliBasisMatrix result;
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                result(row, column) = sm(row, column) + form_factors.f2_gamma.real() * f2_real(row, column) +
                                      form_factors.f2_gamma.imag() * f2_imaginary(row, column) +
                                      form_factors.f3_gamma.real() * f3_real(row, column) +
                                      form_factors.f3_gamma.imag() * f3_imaginary(row, column);
            }
        }
        return result;
    }

    PauliBasisMatrix recompose_full(const PhotonDipoleFormFactors& form_factors) const {
        PauliBasisMatrix result = recompose_linear(form_factors);
        const std::array<double, 4> values{{form_factors.f2_gamma.real(), form_factors.f2_gamma.imag(),
                                             form_factors.f3_gamma.real(), form_factors.f3_gamma.imag()}};
        const std::array<const PauliBasisMatrix*, 10> terms{{
            &f2_real_f2_real, &f2_real_f2_imaginary, &f2_real_f3_real, &f2_real_f3_imaginary,
            &f2_imaginary_f2_imaginary, &f2_imaginary_f3_real, &f2_imaginary_f3_imaginary,
            &f3_real_f3_real, &f3_real_f3_imaginary, &f3_imaginary_f3_imaginary,
        }};
        const std::array<std::array<std::size_t, 2>, 10> indices{{
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 1}},
            {{1, 2}}, {{1, 3}}, {{2, 2}}, {{2, 3}}, {{3, 3}},
        }};
        for (std::size_t term = 0; term < terms.size(); ++term) {
            for (std::size_t row = 0; row < 4; ++row) {
                for (std::size_t column = 0; column < 4; ++column)
                    result(row, column) += values[indices[term][0]] * values[indices[term][1]] *
                                           (*terms[term])(row, column);
            }
        }
        return result;
    }
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_SPIN_DENSITY_H_
