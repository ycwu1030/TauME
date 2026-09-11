#ifndef TAUAMP_TAUTAU_PRODUCTION_H_
#define TAUAMP_TAUTAU_PRODUCTION_H_

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>

#include "tauamp/tautau/kinematics.h"
#include "tauamp/tautau/spin_density.h"

namespace tauamp::tautau {

struct ElectroweakParameters {
    double electron_charge;
    double sin_theta_w;
    double z_mass;
    double z_width;
};

enum class ProductionBosons { photon_and_z, photon_only };

namespace detail {

using Complex = std::complex<double>;

class DiracMatrix {
public:
    Complex& operator()(std::size_t row, std::size_t column) { return entries_.at(4 * row + column); }
    Complex operator()(std::size_t row, std::size_t column) const { return entries_.at(4 * row + column); }

    static DiracMatrix identity() {
        DiracMatrix result;
        for (std::size_t index = 0; index < 4; ++index) result(index, index) = 1.0;
        return result;
    }

private:
    std::array<Complex, 16> entries_{};
};

inline DiracMatrix operator+(const DiracMatrix& left, const DiracMatrix& right) {
    DiracMatrix result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) result(row, column) = left(row, column) + right(row, column);
    }
    return result;
}

inline DiracMatrix operator-(const DiracMatrix& left, const DiracMatrix& right) {
    DiracMatrix result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) result(row, column) = left(row, column) - right(row, column);
    }
    return result;
}

inline DiracMatrix operator*(Complex factor, const DiracMatrix& matrix) {
    DiracMatrix result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) result(row, column) = factor * matrix(row, column);
    }
    return result;
}

inline DiracMatrix operator*(const DiracMatrix& matrix, Complex factor) { return factor * matrix; }

inline DiracMatrix operator*(const DiracMatrix& left, const DiracMatrix& right) {
    DiracMatrix result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            for (std::size_t index = 0; index < 4; ++index) result(row, column) += left(row, index) * right(index, column);
        }
    }
    return result;
}

inline Complex trace(const DiracMatrix& matrix) {
    Complex result{};
    for (std::size_t index = 0; index < 4; ++index) result += matrix(index, index);
    return result;
}

inline DiracMatrix bar(const DiracMatrix& matrix, const DiracMatrix& gamma_zero) {
    DiracMatrix conjugate_transpose;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column)
            conjugate_transpose(row, column) = std::conj(matrix(column, row));
    }
    return gamma_zero * conjugate_transpose * gamma_zero;
}

inline const std::array<DiracMatrix, 4>& gamma_matrices() {
    static const std::array<DiracMatrix, 4> gamma = [] {
        const Complex imaginary{0.0, 1.0};
        DiracMatrix gamma_zero;
        gamma_zero(0, 0) = 1.0;
        gamma_zero(1, 1) = 1.0;
        gamma_zero(2, 2) = -1.0;
        gamma_zero(3, 3) = -1.0;

        const std::array<std::array<Complex, 4>, 3> pauli = {{
            {{0.0, 1.0, 1.0, 0.0}},
            {{0.0, -imaginary, imaginary, 0.0}},
            {{1.0, 0.0, 0.0, -1.0}},
        }};

        std::array<DiracMatrix, 4> result{};
        result[0] = gamma_zero;
        for (std::size_t spatial = 0; spatial < 3; ++spatial) {
            for (std::size_t row = 0; row < 2; ++row) {
                for (std::size_t column = 0; column < 2; ++column) {
                    const Complex entry = pauli[spatial][2 * row + column];
                    result[spatial + 1](row, column + 2) = entry;
                    result[spatial + 1](row + 2, column) = -entry;
                }
            }
        }
        return result;
    }();
    return gamma;
}

inline const DiracMatrix& gamma_five() {
    static const DiracMatrix gamma5 = Complex{0.0, 1.0} * gamma_matrices()[0] * gamma_matrices()[1] *
                                      gamma_matrices()[2] * gamma_matrices()[3];
    return gamma5;
}

inline DiracMatrix slash(const std::array<double, 4>& vector) {
    const std::array<double, 4> metric{{1.0, -1.0, -1.0, -1.0}};
    DiracMatrix result;
    for (std::size_t index = 0; index < 4; ++index) result = result + metric[index] * vector[index] * gamma_matrices()[index];
    return result;
}

inline std::array<double, 4> components(const FourMomentum& vector) {
    return {{vector.energy(), vector.px(), vector.py(), vector.pz()}};
}

inline std::array<double, 4> spin_vector(const FourMomentum& momentum, const std::array<double, 3>& axis, double mass) {
    const double spatial_projection = momentum.px() * axis[0] + momentum.py() * axis[1] + momentum.pz() * axis[2];
    const double factor = spatial_projection / (mass * (momentum.energy() + mass));
    return {{spatial_projection / mass, axis[0] + momentum.px() * factor, axis[1] + momentum.py() * factor,
             axis[2] + momentum.pz() * factor}};
}


struct BosonCouplings {
    double mass;
    double width;
    double electron_vector;
    double electron_axial;
    double tau_vector;
    double tau_axial;
};

inline std::array<BosonCouplings, 2> bosons(const ElectroweakParameters& parameters) {
    if (!(parameters.electron_charge > 0.0) || !(parameters.sin_theta_w > 0.0) ||
        !(parameters.sin_theta_w < 1.0) || !(parameters.z_mass > 0.0) || parameters.z_width < 0.0)
        throw std::invalid_argument("invalid electroweak parameters");
    const double cosine = std::sqrt(1.0 - parameters.sin_theta_w * parameters.sin_theta_w);
    const double scale = parameters.electron_charge / (2.0 * parameters.sin_theta_w * cosine);
    const double vector = scale * (-0.5 + 2.0 * parameters.sin_theta_w * parameters.sin_theta_w);
    const double axial = -0.5 * scale;
    return {{{0.0, 0.0, -parameters.electron_charge, 0.0, -parameters.electron_charge, 0.0},
             {parameters.z_mass, parameters.z_width, vector, axial, vector, axial}}};
}

inline std::array<DiracMatrix, 4> standard_vertices(double vector, double axial) {
    std::array<DiracMatrix, 4> result{};
    for (std::size_t mu = 0; mu < 4; ++mu)
        result[mu] = gamma_matrices()[mu] * (vector * DiracMatrix::identity() - axial * gamma_five());
    return result;
}

inline std::array<DiracMatrix, 4> barred_vertices(const std::array<DiracMatrix, 4>& vertices) {
    std::array<DiracMatrix, 4> result{};
    for (std::size_t mu = 0; mu < 4; ++mu) result[mu] = bar(vertices[mu], gamma_matrices()[0]);
    return result;
}

inline DiracMatrix sigma_q(std::size_t mu, const std::array<double, 4>& q) {
    const std::array<double, 4> metric{{1.0, -1.0, -1.0, -1.0}};
    DiracMatrix result;
    for (std::size_t nu = 0; nu < 4; ++nu) {
        const DiracMatrix sigma = Complex{0.0, 0.5} * (gamma_matrices()[mu] * gamma_matrices()[nu] -
                                                         gamma_matrices()[nu] * gamma_matrices()[mu]);
        result = result + metric[nu] * q[nu] * sigma;
    }
    return result;
}

enum class TauVertexSector { standard_model, f2, f3 };

inline std::array<DiracMatrix, 4> tau_vertices(const BosonCouplings& boson, bool photon, double tau_mass,
                                                 const std::array<double, 4>& q, TauVertexSector sector,
                                                 Complex form_factor) {
    if (!photon || sector == TauVertexSector::standard_model)
        return standard_vertices(boson.tau_vector, boson.tau_axial);
    std::array<DiracMatrix, 4> result{};
    for (std::size_t mu = 0; mu < 4; ++mu) {
        const DiracMatrix dipole = sigma_q(mu, q);
        if (sector == TauVertexSector::f2)
            result[mu] = boson.tau_vector * dipole * (Complex{0.0, 1.0} * form_factor / (2.0 * tau_mass));
        else
            result[mu] = boson.tau_vector * dipole * (form_factor / (2.0 * tau_mass)) * gamma_five();
    }
    return result;
}

inline Complex propagator_denominator(double s, const BosonCouplings& boson) {
    return {s - boson.mass * boson.mass, boson.mass * boson.width};
}

inline Complex propagator_weight(double s, const BosonCouplings& left, const BosonCouplings& right) {
    return 1.0 / (propagator_denominator(s, left) * std::conj(propagator_denominator(s, right)));
}

struct PairContract {
    std::array<std::array<Complex, 4>, 4> lepton_tensor;
    std::array<DiracMatrix, 4> left_tau;
    std::array<DiracMatrix, 4> right_tau_bar;
    Complex weight;

    Complex operator()(const DiracMatrix& minus_projector, const DiracMatrix& plus_projector) const {
        Complex result{};
        for (std::size_t mu = 0; mu < 4; ++mu) {
            for (std::size_t nu = 0; nu < 4; ++nu)
                result += lepton_tensor[mu][nu] * trace(minus_projector * left_tau[mu] * plus_projector * right_tau_bar[nu]);
        }
        return weight * result;
    }
};

inline PairContract pair_contract(const TauPairKinematicPoint& kinematics, const BosonCouplings& left_boson,
                                  const BosonCouplings& right_boson, bool left_is_photon, bool right_is_photon,
                                  TauVertexSector left_sector, TauVertexSector right_sector, Complex left_form_factor,
                                  Complex right_form_factor) {
    const PairCenterOfMass& pair_cm = kinematics.pair_cm();
    const std::array<double, 4> electron = components(pair_cm.electron);
    const std::array<double, 4> positron = components(pair_cm.positron);
    const std::array<double, 4> tau_minus = components(pair_cm.tau_minus);
    const std::array<double, 4> tau_plus = components(pair_cm.tau_plus);
    const std::array<double, 4> q{{tau_minus[0] + tau_plus[0], tau_minus[1] + tau_plus[1], tau_minus[2] + tau_plus[2],
                                    tau_minus[3] + tau_plus[3]}};
    const DiracMatrix electron_density =
        (DiracMatrix::identity() + kinematics.electron_polarization() * gamma_five()) * slash(electron) * 0.5;
    const DiracMatrix positron_density =
        (DiracMatrix::identity() - kinematics.positron_polarization() * gamma_five()) * slash(positron) * 0.5;
    const auto left_electron = standard_vertices(left_boson.electron_vector, left_boson.electron_axial);
    const auto right_electron_bar = barred_vertices(standard_vertices(right_boson.electron_vector, right_boson.electron_axial));
    const auto left_tau_vertices = tau_vertices(left_boson, left_is_photon, kinematics.tau_mass(), q, left_sector, left_form_factor);
    const auto right_tau_bar_vertices =
        barred_vertices(tau_vertices(right_boson, right_is_photon, kinematics.tau_mass(), q, right_sector, right_form_factor));
    const std::array<double, 4> metric{{1.0, -1.0, -1.0, -1.0}};

    PairContract result{};
    for (std::size_t mu = 0; mu < 4; ++mu) {
        result.left_tau[mu] = metric[mu] * left_tau_vertices[mu];
        result.right_tau_bar[mu] = metric[mu] * right_tau_bar_vertices[mu];
        for (std::size_t nu = 0; nu < 4; ++nu)
            result.lepton_tensor[mu][nu] = trace(positron_density * left_electron[mu] * electron_density * right_electron_bar[nu]);
    }
    const double s = kinematics.sqrt_s() * kinematics.sqrt_s();
    result.weight = propagator_weight(s, left_boson, right_boson);
    return result;
}

class ComplexPauliBasisMatrix {
public:
    Complex& operator()(std::size_t row, std::size_t column) { return entries_.at(4 * row + column); }
    Complex operator()(std::size_t row, std::size_t column) const { return entries_.at(4 * row + column); }

private:
    std::array<Complex, 16> entries_{};
};

inline ComplexPauliBasisMatrix coefficients(const TauPairKinematicPoint& kinematics, const PairContract& contract) {
    const PairCenterOfMass& pair_cm = kinematics.pair_cm();
    const DiracMatrix base_minus = slash(components(pair_cm.tau_minus)) + kinematics.tau_mass() * DiracMatrix::identity();
    const DiracMatrix base_plus = slash(components(pair_cm.tau_plus)) - kinematics.tau_mass() * DiracMatrix::identity();
    const auto axes = spin_axes(pair_cm.electron, pair_cm.tau_minus);
    std::array<std::array<Complex, 4>, 4> values{};

    for (std::size_t minus_index = 0; minus_index < 4; ++minus_index) {
        for (std::size_t plus_index = 0; plus_index < 4; ++plus_index) {
            DiracMatrix minus_projector = base_minus;
            DiracMatrix plus_projector = base_plus;
            if (minus_index)
                minus_projector = base_minus * (DiracMatrix::identity() + gamma_five() * slash(spin_vector(pair_cm.tau_minus, ((minus_index == 1) ? axes.transverse : (minus_index == 2) ? axes.normal : axes.longitudinal), kinematics.tau_mass()))) * 0.5;
            if (plus_index)
                plus_projector = base_plus * (DiracMatrix::identity() + gamma_five() * slash(spin_vector(pair_cm.tau_plus, ((plus_index == 1) ? axes.transverse : (plus_index == 2) ? axes.normal : axes.longitudinal), kinematics.tau_mass()))) * 0.5;
            values[minus_index][plus_index] = contract(minus_projector, plus_projector);
        }
    }

    ComplexPauliBasisMatrix result;
    result(0, 0) = values[0][0];
    for (std::size_t minus_index = 1; minus_index < 4; ++minus_index) result(minus_index, 0) = 2.0 * values[minus_index][0] - result(0, 0);
    for (std::size_t plus_index = 1; plus_index < 4; ++plus_index) result(0, plus_index) = 2.0 * values[0][plus_index] - result(0, 0);
    for (std::size_t minus_index = 1; minus_index < 4; ++minus_index) {
        for (std::size_t plus_index = 1; plus_index < 4; ++plus_index)
            result(minus_index, plus_index) = 4.0 * values[minus_index][plus_index] - result(0, 0) -
                                               result(minus_index, 0) - result(0, plus_index);
    }
    return result;
}

inline PauliBasisMatrix real_coefficients(const ComplexPauliBasisMatrix& complex_matrix) {
    PauliBasisMatrix result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            if (std::abs(complex_matrix(row, column).imag()) > 2e-9)
                throw std::runtime_error("summed production coefficient has an unexpected imaginary residue");
            result(row, column) = complex_matrix(row, column).real();
        }
    }
    return result;
}

inline PauliBasisMatrix sum_coefficients(const TauPairKinematicPoint& kinematics, const std::array<BosonCouplings, 2>& model,
                                         TauVertexSector sector, Complex form_factor, ProductionBosons boson_selection) {
    ComplexPauliBasisMatrix total;
    const auto add = [&](std::size_t left, std::size_t right, TauVertexSector left_sector, TauVertexSector right_sector,
                         Complex left_factor, Complex right_factor) {
        const ComplexPauliBasisMatrix contribution = coefficients(
            kinematics, pair_contract(kinematics, model[left], model[right], left == 0, right == 0, left_sector,
                                       right_sector, left_factor, right_factor));
        for (std::size_t row = 0; row < 4; ++row)
            for (std::size_t column = 0; column < 4; ++column) total(row, column) += contribution(row, column);
    };

    if (sector == TauVertexSector::standard_model) {
        const std::size_t boson_count = boson_selection == ProductionBosons::photon_only ? 1U : 2U;
        for (std::size_t left = 0; left < boson_count; ++left) {
            for (std::size_t right = 0; right < boson_count; ++right)
                add(left, right, TauVertexSector::standard_model, TauVertexSector::standard_model, {}, {});
        }
    } else {
        add(0, 0, sector, TauVertexSector::standard_model, form_factor, {});
        add(0, 0, TauVertexSector::standard_model, sector, {}, form_factor);
        if (boson_selection == ProductionBosons::photon_and_z) {
            add(0, 1, sector, TauVertexSector::standard_model, form_factor, {});
            add(1, 0, TauVertexSector::standard_model, sector, {}, form_factor);
        }
    }
    return real_coefficients(total);
}

}  // namespace detail

class ElectronPositronToTauPair {
public:
    explicit ElectronPositronToTauPair(ElectroweakParameters parameters,
                                       ProductionBosons boson_selection = ProductionBosons::photon_and_z)
        : parameters_(parameters), boson_selection_(boson_selection) {}

    SpinDensityComponents components(const TauPairKinematicPoint& kinematics) const {
        const auto model = detail::bosons(parameters_);
        SpinDensityComponents result;
        result.sm = detail::sum_coefficients(kinematics, model, detail::TauVertexSector::standard_model, {}, boson_selection_);
        result.f2_real = detail::sum_coefficients(kinematics, model, detail::TauVertexSector::f2, {1.0, 0.0}, boson_selection_);
        result.f2_imaginary = detail::sum_coefficients(kinematics, model, detail::TauVertexSector::f2, {0.0, 1.0}, boson_selection_);
        result.f3_real = detail::sum_coefficients(kinematics, model, detail::TauVertexSector::f3, {1.0, 0.0}, boson_selection_);
        result.f3_imaginary = detail::sum_coefficients(kinematics, model, detail::TauVertexSector::f3, {0.0, 1.0}, boson_selection_);
        return result;
    }

    SpinDensityPolynomialComponents polynomial_components(const TauPairKinematicPoint& kinematics) const {
        const auto model = detail::bosons(parameters_);
        const SpinDensityComponents linear = components(kinematics);
        SpinDensityPolynomialComponents result;
        result.sm = linear.sm;
        result.f2_real = linear.f2_real;
        result.f2_imaginary = linear.f2_imaginary;
        result.f3_real = linear.f3_real;
        result.f3_imaginary = linear.f3_imaginary;

        const std::array<detail::TauVertexSector, 4> sectors{{detail::TauVertexSector::f2,
                                                                detail::TauVertexSector::f2,
                                                                detail::TauVertexSector::f3,
                                                                detail::TauVertexSector::f3}};
        const std::array<detail::Complex, 4> factors{{detail::Complex{1.0, 0.0}, detail::Complex{0.0, 1.0},
                                                       detail::Complex{1.0, 0.0}, detail::Complex{0.0, 1.0}}};
        const std::array<PauliBasisMatrix*, 10> terms{{
            &result.f2_real_f2_real, &result.f2_real_f2_imaginary, &result.f2_real_f3_real,
            &result.f2_real_f3_imaginary, &result.f2_imaginary_f2_imaginary, &result.f2_imaginary_f3_real,
            &result.f2_imaginary_f3_imaginary, &result.f3_real_f3_real, &result.f3_real_f3_imaginary,
            &result.f3_imaginary_f3_imaginary,
        }};
        const std::array<std::array<std::size_t, 2>, 10> indices{{
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 1}},
            {{1, 2}}, {{1, 3}}, {{2, 2}}, {{2, 3}}, {{3, 3}},
        }};
        for (std::size_t term = 0; term < terms.size(); ++term) {
            const std::size_t left_direction = indices[term][0];
            const std::size_t right_direction = indices[term][1];
            detail::ComplexPauliBasisMatrix total;
            const auto add = [&](std::size_t left, std::size_t right) {
                const detail::ComplexPauliBasisMatrix contribution = detail::coefficients(
                    kinematics, detail::pair_contract(kinematics, model[0], model[0], true, true,
                                                       sectors[left], sectors[right], factors[left], factors[right]));
                for (std::size_t row = 0; row < 4; ++row)
                    for (std::size_t column = 0; column < 4; ++column)
                        total(row, column) += contribution(row, column);
            };
            add(left_direction, right_direction);
            if (left_direction != right_direction) add(right_direction, left_direction);
            *terms[term] = detail::real_coefficients(total);
        }
        return result;
    }

private:
    ElectroweakParameters parameters_;
    ProductionBosons boson_selection_;
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_PRODUCTION_H_
