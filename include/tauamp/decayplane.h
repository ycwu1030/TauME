#ifndef TAU_OO_DECAYPLANE_H_
#define TAU_OO_DECAYPLANE_H_

#include "constants.h"
#include "utilities.h"

namespace tauamp {
class DecayPlaneRhoRho_t {
public:
    void set_momenta(std::vector<dlv_t> p_taum_decays_list, std::vector<dlv_t> p_taup_decays_list) {
        p_pim_lab = p_taum_decays_list[1];
        p_pip_lab = p_taup_decays_list[1];
        p_pi0m_lab = p_taum_decays_list[2];
        p_pi0p_lab = p_taup_decays_list[2];
        boost_to_zmf();
        calc_Phi_CP();
    }
    double get_Phi_CP() const { return Phi_CP; }

private:
    void boost_to_zmf() {
        dlv_t p_pic_tot = p_pim_lab + p_pip_lab;
        auto bst = p_pic_tot.BoostToCM();

        p_pim_zmf = ROOT::Math::VectorUtil::boost(p_pim_lab, bst);
        p_pip_zmf = ROOT::Math::VectorUtil::boost(p_pip_lab, bst);
        p_pi0m_zmf = ROOT::Math::VectorUtil::boost(p_pi0m_lab, bst);
        p_pi0p_zmf = ROOT::Math::VectorUtil::boost(p_pi0p_lab, bst);
        // std::cout << "p_pim_zmf" << p_pim_zmf << std::endl;
        // std::cout << "p_pip_zmf" << p_pip_zmf << std::endl;
    }

    void calc_Phi_CP() {
        double Epip = p_pip_lab.E();
        double Epim = p_pim_lab.E();
        double Epi0p = p_pi0p_lab.E();
        double Epi0m = p_pi0m_lab.E();

        double Y = (Epip - Epi0p) * (Epim - Epi0m);
        auto pip3 = p_pip_zmf.Vect().Unit();
        auto pim3 = p_pim_zmf.Vect().Unit();
        auto pi0p3 = p_pi0p_zmf.Vect().Unit();
        auto pi0m3 = p_pi0m_zmf.Vect().Unit();

        auto pi0p3_perp = ROOT::Math::VectorUtil::PerpVector(pi0p3, pip3);
        auto pi0m3_perp = ROOT::Math::VectorUtil::PerpVector(pi0m3, pim3);

        // std::cout << "-------------" << std::endl;
        // std::cout << "pip3: " << pip3.x() << "," << pip3.y() << "," << pip3.z() << std::endl;
        // std::cout << "pim3: " << pim3.x() << "," << pim3.y() << "," << pim3.z() << std::endl;
        // std::cout << "pi0p3: " << pi0p3.x() << "," << pi0p3.y() << "," << pi0p3.z() << std::endl;
        // std::cout << "pi0m3: " << pi0m3.x() << "," << pi0m3.y() << "," << pi0m3.z() << std::endl;
        // std::cout << "pi0p3_perp: " << pi0p3_perp.x() << "," << pi0p3_perp.y() << "," << pi0p3_perp.z() << std::endl;
        // std::cout << "pi0m3_perp: " << pi0m3_perp.x() << "," << pi0m3_perp.y() << "," << pi0m3_perp.z() << std::endl;
        // std::cout << "-------------" << std::endl;

        double sign = (pi0m3_perp.Cross(pi0p3_perp)).Dot(pim3) >= 0 ? 1 : -1;
        double X = acos(pi0p3_perp.Unit().Dot(pi0m3_perp.Unit())) * sign;

        Phi_CP = X;
        if (Y < 0) {
            if (Phi_CP > 0) {
                Phi_CP -= M_PI;
            } else {
                Phi_CP += M_PI;
            }
        }
    }

    dlv_t p_pim_lab;
    dlv_t p_pip_lab;

    dlv_t p_pi0m_lab;
    dlv_t p_pi0p_lab;

    dlv_t p_pim_zmf;
    dlv_t p_pip_zmf;
    dlv_t p_pi0m_zmf;
    dlv_t p_pi0p_zmf;

    double Phi_CP;
};
}  // namespace tauamp
#endif  // TAU_OO_DECAYPLANE_H_
