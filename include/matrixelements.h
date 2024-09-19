#ifndef TAU_OO_MATRIX_ELEMENTS_H_
#define TAU_OO_MATRIX_ELEMENTS_H_

#include <iostream>

#include "constants.h"
#include "utilities.h"

template <class TAUM_DECAY_t, class TAUP_DECAY_t>
class ME_Base_t {
public:
    ME_Base_t() {}
    virtual ~ME_Base_t() {}

    void set_momenta(std::vector<dlv_t> p_production_list, std::vector<dlv_t> p_taum_decays_list,
                     std::vector<dlv_t> p_taup_decays_list) {
        set_taum_decay_momenta(p_taum_decays_list);
        set_taup_decay_momenta(p_taup_decays_list);
        set_production_momenta(p_production_list);
    }
    // * Suppose the BSM contribution is proportional to some parameter c1;
    // * ME2 = ME2_SM + c1 * ME2_Interference + c1^2 ME2_BSM;
    double get_ME2_SM() const { return m_ME2[0]; }
    double get_ME2_Interference() const { return m_ME2[1]; }
    double get_ME2_BSM() const { return m_ME2[2]; }

protected:
    virtual void set_production_momenta(
        std::vector<dlv_t> p_list) = 0;  // * User should make sure setting momenta for taupm decay first
    void set_taum_decay_momenta(std::vector<dlv_t> p_list) {
        m_taum_decay.set_momenta(p_list);
        m_H_m = m_taum_decay.get_H();
        m_omega_m = m_taum_decay.get_omega();
    }
    void set_taup_decay_momenta(std::vector<dlv_t> p_list) {
        m_taup_decay.set_momenta(p_list);
        m_H_p = m_taup_decay.get_H();
        m_omega_p = m_taup_decay.get_omega();
    }

    TAUM_DECAY_t m_taum_decay;
    TAUP_DECAY_t m_taup_decay;

    clv_t m_H_m;
    clv_t m_H_p;

    cd_t m_omega_m;
    cd_t m_omega_p;

    double m_ME2[3];
};

template <class TAUM_DECAY_t, class TAUP_DECAY_t>
class ME_EDM_Re_t : public ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t> {
    // * For e+ e- > gamma* > tau+ tau-
    // * Including EDM operator between gamma and tau lepton
public:
    ME_EDM_Re_t() {}
    ~ME_EDM_Re_t() {}

protected:
    virtual void set_production_momenta(std::vector<dlv_t> p_list) override {
        clv_t _p_ep = ToComplex(p_list[0]);
        clv_t _p_em = ToComplex(p_list[1]);
        clv_t _p_taup = ToComplex(p_list[2]);
        clv_t _p_taum = ToComplex(p_list[3]);
        clv_t _p_s = _p_ep + _p_em;
        clv_t _H_m = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_H_m;
        clv_t _H_p = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_H_p;

        cd_t SP56 = _H_m.Dot(_H_p);
        cd_t SP50 = _H_m.Dot(_p_em);
        cd_t SP51 = _H_m.Dot(_p_ep);
        cd_t SP54 = _H_m.Dot(_p_taup);
        cd_t SP52 = _H_m.Dot(_p_s);
        cd_t SP60 = _H_p.Dot(_p_em);
        cd_t SP61 = _H_p.Dot(_p_ep);
        cd_t SP63 = _H_p.Dot(_p_taum);
        cd_t SP62 = _H_p.Dot(_p_s);
        cd_t SP01 = _p_em.Dot(_p_ep);
        cd_t SP03 = _p_em.Dot(_p_taum);
        cd_t SP04 = _p_em.Dot(_p_taup);
        cd_t SP02 = _p_em.Dot(_p_s);
        cd_t SP13 = _p_ep.Dot(_p_taum);
        cd_t SP14 = _p_ep.Dot(_p_taup);
        cd_t SP12 = _p_ep.Dot(_p_s);
        cd_t SP34 = _p_taum.Dot(_p_taup);
        cd_t SP32 = _p_taum.Dot(_p_s);
        cd_t SP42 = _p_taup.Dot(_p_s);
        cd_t SP22 = _p_s.Dot(_p_s);
        // std::cout << "q.q = " << SP22 << std::endl;

        cd_t EP5603 = epsilon(_H_m, _H_p, _p_em, _p_taum);
        cd_t EP5602 = epsilon(_H_m, _H_p, _p_em, _p_s);
        cd_t EP5613 = epsilon(_H_m, _H_p, _p_ep, _p_taum);
        cd_t EP5612 = epsilon(_H_m, _H_p, _p_ep, _p_s);
        cd_t EP5642 = epsilon(_H_m, _H_p, _p_taup, _p_s);
        cd_t EP5032 = epsilon(_H_m, _p_em, _p_taum, _p_s);
        cd_t EP5042 = epsilon(_H_m, _p_em, _p_taup, _p_s);
        cd_t EP5132 = epsilon(_H_m, _p_ep, _p_taum, _p_s);
        cd_t EP5142 = epsilon(_H_m, _p_ep, _p_taup, _p_s);
        cd_t EP6032 = epsilon(_H_p, _p_em, _p_taum, _p_s);
        cd_t EP6042 = epsilon(_H_p, _p_em, _p_taup, _p_s);
        cd_t EP6132 = epsilon(_H_p, _p_ep, _p_taum, _p_s);
        cd_t EP6142 = epsilon(_H_p, _p_ep, _p_taup, _p_s);

        // * The calculation procedure is as following
        // * For e+ e- > gamma* > tau+ tau- we have
        // * ME2 ~ M0 + M1.sp + M2.sm + M3.sp.sm
        // * For tau+ > nu + xxxxx
        // * ME2 ~ omega_p + H_p.sp
        // * For tau- > nu + xxxxx
        // * ME2 ~ omega_m + H_m.sp
        // * Then the total ME2 reads
        // * ME2 ~ M0 omega_p omega_m - M1.H_p - M2.H_m + M3.H_p.H_m
        // * And the dependence on the BSM coefficient is given by
        // * M0 ~ M01 + M02 c1 + M03 c1^2
        // * M1 ~ M11 + M12 c1 + M13 c1^2
        // * M2 ~ M21 + M22 c1 + M23 c1^2
        // * M3 ~ M31 + M32 c1 + M33 c1^2
        // * ME2 ~     (M01 omega_p omega_m - M11.Hp - M21.Hm + M31.Hp.Hm)  // -> _ME2[0]
        // *      + c1 (M02 omega_p omega_m - M12.Hp - M22.Hm + M32.Hp.Hm)  // -> _ME2[1]
        // *    + c1^2 (M03 omega_p omega_m - M13.Hp - M23.Hm + M33.Hp.Hm); // -> _ME2[2]
        // * But in the following, we will use the following convention
        // * ME2 ~     (M01 omega_p omega_m - M11 - M21 + M31) // -> _ME2[0]
        // *      + c1 (M02 omega_p omega_m - M12 - M22 + M32) // -> _ME2[1]
        // *    + c1^2 (M03 omega_p omega_m - M13 - M23 + M33) // -> _ME2[2]
        // * Which mean in the calculation, we already dotted with the Hp/Hm;

        cd_t M01 = 8 * pow(MTAU, 2) * SP01 + 8 * (SP04 * SP13 + SP03 * SP14);
        cd_t M02 = 0;
        cd_t M03 = -8 * SP04 * SP13 * SP22 - 8 * SP03 * SP14 * SP22 -
                   4 * pow(MTAU, 2) * (2 * SP02 * SP12 + SP01 * SP22) + 8 * SP04 * SP12 * SP32 +
                   8 * SP02 * SP14 * SP32 - 8 * SP02 * SP12 * SP34 + 4 * SP01 * SP22 * SP34 + 8 * SP03 * SP12 * SP42 +
                   8 * SP02 * SP13 * SP42;
        cd_t M11 = 0.0;
        cd_t M12 = 0.0;
        cd_t M13 = 0.0;
        cd_t M21 = 0.0;
        cd_t M22 = 0.0;
        cd_t M23 = 0.0;
        cd_t M31 = -8 * (SP04 * SP13 * SP56 + SP03 * SP14 * SP56 - SP01 * SP34 * SP56 + pow(MTAU, 2) * SP51 * SP60 +
                         SP34 * SP51 * SP60 - SP13 * SP54 * SP60 + pow(MTAU, 2) * SP50 * SP61 + SP34 * SP50 * SP61 -
                         SP03 * SP54 * SP61 - SP14 * SP50 * SP63 - SP04 * SP51 * SP63 + SP01 * SP54 * SP63);
        cd_t M32 = -4 * MTAU *
                   (-2 * EP5642 * SP01 - EP5613 * SP02 - EP5612 * SP03 + 2 * EP5612 * SP04 - EP5603 * SP12 -
                    EP5602 * SP13 + 2 * EP5602 * SP14 + EP6132 * SP50 + 2 * EP6142 * SP50 + EP6032 * SP51 +
                    2 * EP6042 * SP51 + 3 * EP5132 * SP60 + 2 * EP5142 * SP60 + 3 * EP5032 * SP61 + 2 * EP5042 * SP61);
        cd_t M33 = 4 * (-2 * SP03 * SP14 * SP22 * SP56 + 2 * SP02 * SP14 * SP32 * SP56 - 2 * SP02 * SP12 * SP34 * SP56 +
                        3 * SP01 * SP22 * SP34 * SP56 + 2 * SP03 * SP12 * SP42 * SP56 + 2 * SP02 * SP13 * SP42 * SP56 -
                        4 * SP01 * SP32 * SP42 * SP56 - 2 * SP22 * SP34 * SP51 * SP60 + 4 * SP32 * SP42 * SP51 * SP60 +
                        2 * SP12 * SP34 * SP52 * SP60 - 4 * SP13 * SP42 * SP52 * SP60 + 2 * SP13 * SP22 * SP54 * SP60 -
                        2 * SP12 * SP32 * SP54 * SP60 - 2 * SP22 * SP34 * SP50 * SP61 + 4 * SP32 * SP42 * SP50 * SP61 +
                        2 * SP02 * SP34 * SP52 * SP61 - 4 * SP03 * SP42 * SP52 * SP61 + 2 * SP03 * SP22 * SP54 * SP61 -
                        2 * SP02 * SP32 * SP54 * SP61 - 4 * SP14 * SP32 * SP50 * SP62 + 2 * SP12 * SP34 * SP50 * SP62 +
                        2 * SP02 * SP34 * SP51 * SP62 + 4 * SP03 * SP14 * SP52 * SP62 - 4 * SP01 * SP34 * SP52 * SP62 -
                        2 * SP03 * SP12 * SP54 * SP62 - 2 * SP02 * SP13 * SP54 * SP62 + 4 * SP01 * SP32 * SP54 * SP62 +
                        pow(MTAU, 2) * (SP01 * SP22 * SP56 - 2 * SP22 * SP51 * SP60 + 2 * SP12 * SP52 * SP60 -
                                        2 * SP22 * SP50 * SP61 + 2 * SP12 * SP50 * SP62 +
                                        SP02 * (-2 * SP12 * SP56 + 2 * SP52 * SP61 + 2 * SP51 * SP62)) +
                        2 * SP14 * SP22 * SP50 * SP63 - 2 * SP12 * SP42 * SP50 * SP63 - 2 * SP02 * SP42 * SP51 * SP63 -
                        2 * SP02 * SP14 * SP52 * SP63 + 4 * SP01 * SP42 * SP52 * SP63 + 2 * SP02 * SP12 * SP54 * SP63 -
                        3 * SP01 * SP22 * SP54 * SP63 -
                        2 * SP04 *
                            (SP13 * SP22 * SP56 - SP12 * SP32 * SP56 + 2 * SP32 * SP51 * SP62 - 2 * SP13 * SP52 * SP62 -
                             SP22 * SP51 * SP63 + SP12 * SP52 * SP63));

        cd_t _omega_m = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_omega_m;
        cd_t _omega_p = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_omega_p;
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_ME2[0] = real(M01 * _omega_m * _omega_p - M11 - M21 + M31);
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_ME2[1] = real(M02 * _omega_m * _omega_p - M12 - M22 + M32);
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::m_ME2[2] = real(M03 * _omega_m * _omega_p - M13 - M23 + M33);
    }
};

#endif  // TAU_OO_MATRIX_ELEMENTS_H_
