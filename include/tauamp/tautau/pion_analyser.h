#ifndef TAUAMP_TAUTAU_PION_ANALYSER_H_
#define TAUAMP_TAUTAU_PION_ANALYSER_H_

#include <array>
#include <stdexcept>

#include "tauamp/tautau/kinematics.h"

namespace tauamp::tautau {

namespace detail {
inline std::array<double, 3> pion_direction_in_tau_rest_tnl(const TauPairKinematicPoint& kinematics,
                                                              const FourMomentum& pion_lab, bool tau_minus) {
    const PairCenterOfMass& pair_cm = kinematics.pair_cm();
    const FourMomentum& tau_cm = tau_minus ? pair_cm.tau_minus : pair_cm.tau_plus;
    const FourMomentum pion_cm = kinematics.to_pair_cm(pion_lab);
    const std::array<double, 3> to_tau_rest{{-tau_cm.px() / tau_cm.energy(), -tau_cm.py() / tau_cm.energy(),
                                               -tau_cm.pz() / tau_cm.energy()}};
    const std::array<double, 3> pion_direction = normalized(
        spatial_components(pion_cm.boosted(to_tau_rest)), "pion momentum cannot vanish in the tau rest frame");
    const SpinAxes axes = spin_axes(pair_cm.electron, pair_cm.tau_minus);
    return {{spatial_dot(pion_direction, axes.transverse), spatial_dot(pion_direction, axes.normal),
             spatial_dot(pion_direction, axes.longitudinal)}};
}
}  // namespace detail

inline std::array<double, 3> tau_minus_pion_analyser(const TauPairKinematicPoint& kinematics,
                                                       const FourMomentum& pion_minus_lab) {
    return detail::pion_direction_in_tau_rest_tnl(kinematics, pion_minus_lab, true);
}

inline std::array<double, 3> tau_plus_pion_analyser(const TauPairKinematicPoint& kinematics,
                                                      const FourMomentum& pion_plus_lab) {
    const std::array<double, 3> pion_direction = detail::pion_direction_in_tau_rest_tnl(kinematics, pion_plus_lab, false);
    return {{-pion_direction[0], -pion_direction[1], -pion_direction[2]}};
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_PION_ANALYSER_H_
