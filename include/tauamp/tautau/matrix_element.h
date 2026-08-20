#ifndef TAUAMP_TAUTAU_MATRIX_ELEMENT_H_
#define TAUAMP_TAUTAU_MATRIX_ELEMENT_H_

#include <array>

#include "tauamp/tautau/pion_analyser.h"
#include "tauamp/tautau/pion_pair_hypothesis.h"
#include "tauamp/tautau/production.h"

namespace tauamp::tautau {

inline double contract_pion_analysers(const PauliBasisMatrix& matrix, const std::array<double, 3>& tau_minus_analyser,
                                      const std::array<double, 3>& tau_plus_analyser) {
    double result = matrix(0, 0);
    for (std::size_t index = 0; index < 3; ++index) {
        result += matrix(index + 1, 0) * tau_minus_analyser[index];
        result += matrix(0, index + 1) * tau_plus_analyser[index];
    }
    for (std::size_t minus_index = 0; minus_index < 3; ++minus_index) {
        for (std::size_t plus_index = 0; plus_index < 3; ++plus_index)
            result += matrix(minus_index + 1, plus_index + 1) * tau_minus_analyser[minus_index] *
                      tau_plus_analyser[plus_index];
    }
    return result;
}

class PionPairMatrixElement {
public:
    explicit PionPairMatrixElement(ElectroweakParameters parameters) : production_(parameters) {}

    LinearComponents components(const PionPairKinematicHypothesis& hypothesis) const {
        return components(hypothesis.point, hypothesis.pion_minus_lab, hypothesis.pion_plus_lab);
    }

    LinearComponents components(const TauPairKinematicPoint& kinematics, const FourMomentum& pion_minus_lab,
                                const FourMomentum& pion_plus_lab) const {
        const auto tau_minus_analyser = tau_minus_pion_analyser(kinematics, pion_minus_lab);
        const auto tau_plus_analyser = tau_plus_pion_analyser(kinematics, pion_plus_lab);
        const SpinDensityComponents production_components = production_.components(kinematics);
        return {contract_pion_analysers(production_components.sm, tau_minus_analyser, tau_plus_analyser),
                contract_pion_analysers(production_components.f2_real, tau_minus_analyser, tau_plus_analyser),
                contract_pion_analysers(production_components.f2_imaginary, tau_minus_analyser, tau_plus_analyser),
                contract_pion_analysers(production_components.f3_real, tau_minus_analyser, tau_plus_analyser),
                contract_pion_analysers(production_components.f3_imaginary, tau_minus_analyser, tau_plus_analyser)};
    }

private:
    ElectronPositronToTauPair production_;
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_MATRIX_ELEMENT_H_
