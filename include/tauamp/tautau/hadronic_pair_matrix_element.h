#ifndef TAUAMP_TAUTAU_HADRONIC_PAIR_MATRIX_ELEMENT_H_
#define TAUAMP_TAUTAU_HADRONIC_PAIR_MATRIX_ELEMENT_H_

#include <array>
#include <complex>

#include "tauamp/tautau/hadronic_analyser.h"
#include "tauamp/tautau/production.h"

namespace tauamp::tautau {

namespace detail {
inline ComplexFourVector boost_current(const ComplexFourVector& current, const std::array<double, 3>& beta) {
    const FourMomentum real_part(std::real(current[1]), std::real(current[2]), std::real(current[3]), std::real(current[0]));
    const FourMomentum imaginary_part(std::imag(current[1]), std::imag(current[2]), std::imag(current[3]), std::imag(current[0]));
    const FourMomentum real_boosted = real_part.boosted(beta);
    const FourMomentum imaginary_boosted = imaginary_part.boosted(beta);
    return {{std::complex<double>(real_boosted.energy(), imaginary_boosted.energy()),
             std::complex<double>(real_boosted.px(), imaginary_boosted.px()),
             std::complex<double>(real_boosted.py(), imaginary_boosted.py()),
             std::complex<double>(real_boosted.pz(), imaginary_boosted.pz())}};
}

inline ComplexFourVector to_pair_cm(const TauPairKinematicPoint& kinematics, const ComplexFourVector& current) {
    const FourMomentum real_part(std::real(current[1]), std::real(current[2]), std::real(current[3]), std::real(current[0]));
    const FourMomentum imaginary_part(std::imag(current[1]), std::imag(current[2]), std::imag(current[3]), std::imag(current[0]));
    const FourMomentum real_cm = kinematics.to_pair_cm(real_part);
    const FourMomentum imaginary_cm = kinematics.to_pair_cm(imaginary_part);
    return {{std::complex<double>(real_cm.energy(), imaginary_cm.energy()),
             std::complex<double>(real_cm.px(), imaginary_cm.px()),
             std::complex<double>(real_cm.py(), imaginary_cm.py()),
             std::complex<double>(real_cm.pz(), imaginary_cm.pz())}};
}

inline double contract_decay_analysers(const PauliBasisMatrix& matrix,
                                       const std::array<double, 3>& tau_minus_analyser,
                                       const std::array<double, 3>& tau_plus_analyser) {
    double result = matrix(0, 0);
    for (std::size_t index = 0; index < 3; ++index) {
        result += matrix(index + 1, 0) * tau_minus_analyser[index];
        result += matrix(0, index + 1) * tau_plus_analyser[index];
    }
    for (std::size_t minus_index = 0; minus_index < 3; ++minus_index)
        for (std::size_t plus_index = 0; plus_index < 3; ++plus_index)
            result += matrix(minus_index + 1, plus_index + 1) * tau_minus_analyser[minus_index] *
                      tau_plus_analyser[plus_index];
    return result;
}
}  // namespace detail

inline std::array<double, 3> hadronic_analyser_tnl(const TauPairKinematicPoint& kinematics,
                                                   const FourMomentum& neutrino_lab,
                                                   const ComplexFourVector& current_lab, bool tau_plus) {
    const PairCenterOfMass& pair_cm = kinematics.pair_cm();
    const FourMomentum& tau_cm = tau_plus ? pair_cm.tau_plus : pair_cm.tau_minus;
    const std::array<double, 3> to_tau_rest{{-tau_cm.px() / tau_cm.energy(), -tau_cm.py() / tau_cm.energy(),
                                             -tau_cm.pz() / tau_cm.energy()}};
    const FourMomentum neutrino_rest = kinematics.to_pair_cm(neutrino_lab).boosted(to_tau_rest);
    const ComplexFourVector current_rest = detail::boost_current(detail::to_pair_cm(kinematics, current_lab), to_tau_rest);
    const DecayRestCoefficients coefficients =
        decay_rest_coefficients(tau_plus ? kinematics.tau_plus_mass() : kinematics.tau_minus_mass(), neutrino_rest,
                                current_rest, tau_plus);
    const std::array<double, 3> analyser = normalized_analyser(coefficients);
    const detail::SpinAxes axes = detail::spin_axes(pair_cm.electron, pair_cm.tau_minus);
    return {{detail::spatial_dot(analyser, axes.transverse), detail::spatial_dot(analyser, axes.normal),
             detail::spatial_dot(analyser, axes.longitudinal)}};
}

class HadronicPairMatrixElement {
public:
    explicit HadronicPairMatrixElement(ElectroweakParameters parameters,
                                       ProductionBosons boson_selection = ProductionBosons::photon_and_z)
        : production_(parameters, boson_selection) {}

    PolynomialComponents polynomial_components(const TauPairKinematicPoint& kinematics,
                                               const std::array<double, 3>& tau_minus_analyser,
                                               const std::array<double, 3>& tau_plus_analyser) const {
        const SpinDensityPolynomialComponents production = production_.polynomial_components(kinematics);
        return {
            detail::contract_decay_analysers(production.sm, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f3_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f3_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_real_f2_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_real_f2_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_real_f3_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_real_f3_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_imaginary_f2_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_imaginary_f3_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f2_imaginary_f3_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f3_real_f3_real, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f3_real_f3_imaginary, tau_minus_analyser, tau_plus_analyser),
            detail::contract_decay_analysers(production.f3_imaginary_f3_imaginary, tau_minus_analyser, tau_plus_analyser)};
    }

private:
    ElectronPositronToTauPair production_;
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_HADRONIC_PAIR_MATRIX_ELEMENT_H_
