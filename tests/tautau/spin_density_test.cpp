#include <cassert>
#include <cmath>
#include <complex>

#include "tauamp/tautau/spin_density.h"

int main() {
    using tauamp::tautau::PauliBasisMatrix;
    using tauamp::tautau::PhotonDipoleFormFactors;
    using tauamp::tautau::SpinDensityComponents;

    SpinDensityComponents components;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double base = static_cast<double>(4 * row + column + 1);
            components.sm(row, column) = base;
            components.f2_real(row, column) = 2.0 * base;
            components.f2_imaginary(row, column) = -3.0 * base;
            components.f3_real(row, column) = 5.0 * base;
            components.f3_imaginary(row, column) = -7.0 * base;
        }
    }
    const PhotonDipoleFormFactors form_factors{{0.2, -0.4}, {-0.3, 0.6}};
    const PauliBasisMatrix result = components.recompose(form_factors);
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double base = static_cast<double>(4 * row + column + 1);
            assert(std::abs(result(row, column) - base * (1.0 + 0.4 + 1.2 - 1.5 - 4.2)) < 1e-12);
        }
    }
}
