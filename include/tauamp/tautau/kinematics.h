#ifndef TAUAMP_TAUTAU_KINEMATICS_H_
#define TAUAMP_TAUTAU_KINEMATICS_H_

#include <array>
#include <cmath>
#include <stdexcept>

#include "tauamp/tautau/four_momentum.h"

namespace tauamp::tautau {

namespace detail {
constexpr double kKinematicTolerance = 1e-9;

inline std::array<double, 3> spatial_components(const FourMomentum& momentum) {
    return {{momentum.px(), momentum.py(), momentum.pz()}};
}

inline double spatial_dot(const std::array<double, 3>& left, const std::array<double, 3>& right) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

inline std::array<double, 3> spatial_cross(const std::array<double, 3>& left, const std::array<double, 3>& right) {
    return {{left[1] * right[2] - left[2] * right[1], left[2] * right[0] - left[0] * right[2],
             left[0] * right[1] - left[1] * right[0]}};
}

inline double spatial_norm(const std::array<double, 3>& vector) {
    return std::sqrt(spatial_dot(vector, vector));
}

inline std::array<double, 3> normalized(const std::array<double, 3>& vector, const char* message) {
    const double norm = spatial_norm(vector);
    if (norm < kKinematicTolerance) throw std::invalid_argument(message);
    return {{vector[0] / norm, vector[1] / norm, vector[2] / norm}};
}

inline FourMomentum rotate_to_electron_axis(const FourMomentum& value, const FourMomentum& electron) {
    const std::array<double, 3> z_axis = normalized(spatial_components(electron), "electron momentum cannot vanish in pair CM");
    std::array<double, 3> reference{{0.0, 1.0, 0.0}};
    if (std::abs(spatial_dot(reference, z_axis)) > 0.9) reference = {{1.0, 0.0, 0.0}};
    const std::array<double, 3> x_axis = normalized(spatial_cross(reference, z_axis), "electron-axis rotation is undefined");
    const std::array<double, 3> y_axis = spatial_cross(z_axis, x_axis);
    const std::array<double, 3> spatial = spatial_components(value);
    return {spatial_dot(spatial, x_axis), spatial_dot(spatial, y_axis), spatial_dot(spatial, z_axis), value.energy()};
}

struct SpinAxes {
    std::array<double, 3> transverse;
    std::array<double, 3> normal;
    std::array<double, 3> longitudinal;
};

inline SpinAxes spin_axes(const FourMomentum& electron, const FourMomentum& tau_minus) {
    const std::array<double, 3> longitudinal = normalized(spatial_components(tau_minus), "spin basis is undefined at tau-pair threshold");
    const std::array<double, 3> beam_cross_tau = spatial_cross(spatial_components(electron), spatial_components(tau_minus));
    const std::array<double, 3> normal = normalized(
        {{-beam_cross_tau[0], -beam_cross_tau[1], -beam_cross_tau[2]}},
        "spin basis is undefined for exactly forward or backward tau production");
    return {spatial_cross(normal, longitudinal), normal, longitudinal};
}
}  // namespace detail

class BeamState {
public:
    BeamState(const FourMomentum& electron_lab, const FourMomentum& positron_lab, double electron_polarization,
              double positron_polarization)
        : electron_lab_(electron_lab), positron_lab_(positron_lab), electron_polarization_(electron_polarization),
          positron_polarization_(positron_polarization) {
        if (std::abs(electron_polarization_) > 1.0 || std::abs(positron_polarization_) > 1.0)
            throw std::invalid_argument("longitudinal polarization must lie in [-1, 1]");
        const FourMomentum total = electron_lab_ + positron_lab_;
        if (total.mass_squared() <= 0.0 || total.energy() <= 0.0)
            throw std::invalid_argument("beam pair must have positive time-like invariant mass");
        const std::array<double, 3> electron_spatial = detail::spatial_components(electron_lab_);
        const std::array<double, 3> positron_spatial = detail::spatial_components(positron_lab_);
        const double electron_norm = detail::spatial_norm(electron_spatial);
        const double positron_norm = detail::spatial_norm(positron_spatial);
        if (electron_norm < detail::kKinematicTolerance || positron_norm < detail::kKinematicTolerance)
            throw std::invalid_argument("beam momentum cannot vanish");
        if (detail::spatial_norm(detail::spatial_cross(electron_spatial, positron_spatial)) >
            detail::kKinematicTolerance * electron_norm * positron_norm)
            throw std::invalid_argument("non-collinear beam crossing angle is not implemented");
    }

    const FourMomentum& electron_lab() const { return electron_lab_; }
    const FourMomentum& positron_lab() const { return positron_lab_; }
    double electron_polarization() const { return electron_polarization_; }
    double positron_polarization() const { return positron_polarization_; }
    double sqrt_s() const { return std::sqrt((electron_lab_ + positron_lab_).mass_squared()); }

private:
    FourMomentum electron_lab_;
    FourMomentum positron_lab_;
    double electron_polarization_;
    double positron_polarization_;
};

struct PairCenterOfMass {
    FourMomentum electron;
    FourMomentum positron;
    FourMomentum tau_minus;
    FourMomentum tau_plus;
};

// A complete, conditional lab-frame tau-pair point for the matrix-element kernel.
// This value does not represent a detector observation and never infers parent-tau
// momenta from visible decay products.
class TauPairKinematicPoint {
public:
    TauPairKinematicPoint(const BeamState& beams, const FourMomentum& tau_minus_lab, const FourMomentum& tau_plus_lab,
                      double tau_mass)
        : beams_(beams), tau_minus_lab_(tau_minus_lab), tau_plus_lab_(tau_plus_lab), tau_mass_(tau_mass) {
        if (tau_mass_ <= 0.0) throw std::invalid_argument("tau mass must be positive");
        const FourMomentum total = beams_.electron_lab() + beams_.positron_lab();
        const FourMomentum final_total = tau_minus_lab_ + tau_plus_lab_;
        if (std::abs(total.px() - final_total.px()) > detail::kKinematicTolerance ||
            std::abs(total.py() - final_total.py()) > detail::kKinematicTolerance ||
            std::abs(total.pz() - final_total.pz()) > detail::kKinematicTolerance ||
            std::abs(total.energy() - final_total.energy()) > detail::kKinematicTolerance)
            throw std::invalid_argument("tau-pair momenta must conserve the beam four-momentum");
        if (sqrt_s() <= 2.0 * tau_mass_) throw std::invalid_argument("pair invariant mass must exceed tau-pair threshold");

        to_pair_cm_beta_ = {{-total.px() / total.energy(), -total.py() / total.energy(), -total.pz() / total.energy()}};
        electron_cm_before_rotation_ = beams_.electron_lab().boosted(to_pair_cm_beta_);
        pair_cm_ = {detail::rotate_to_electron_axis(electron_cm_before_rotation_, electron_cm_before_rotation_),
                    detail::rotate_to_electron_axis(beams_.positron_lab().boosted(to_pair_cm_beta_), electron_cm_before_rotation_),
                    detail::rotate_to_electron_axis(tau_minus_lab_.boosted(to_pair_cm_beta_), electron_cm_before_rotation_),
                    detail::rotate_to_electron_axis(tau_plus_lab_.boosted(to_pair_cm_beta_), electron_cm_before_rotation_)};
    }

    const PairCenterOfMass& pair_cm() const { return pair_cm_; }
    const FourMomentum& tau_minus_lab() const { return tau_minus_lab_; }
    const FourMomentum& tau_plus_lab() const { return tau_plus_lab_; }
    FourMomentum to_pair_cm(const FourMomentum& lab_momentum) const {
        return detail::rotate_to_electron_axis(lab_momentum.boosted(to_pair_cm_beta_), electron_cm_before_rotation_);
    }
    double sqrt_s() const { return beams_.sqrt_s(); }
    double tau_mass() const { return tau_mass_; }
    double electron_polarization() const { return beams_.electron_polarization(); }
    double positron_polarization() const { return beams_.positron_polarization(); }

private:
    BeamState beams_;
    FourMomentum tau_minus_lab_;
    FourMomentum tau_plus_lab_;
    double tau_mass_;
    std::array<double, 3> to_pair_cm_beta_{};
    FourMomentum electron_cm_before_rotation_;
    PairCenterOfMass pair_cm_;
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_KINEMATICS_H_
