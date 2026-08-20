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

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_SPIN_DENSITY_H_
