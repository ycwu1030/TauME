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

    virtual void set_production_momenta(
        std::vector<vd_t> p_list) = 0;  // * User should make sure setting momenta for taupm decay first
    void set_taum_decay_momenta(std::vector<vd_t> p_list) {
        _taum_decay.set_momenta(p_list);
        _H_m = _taum_decay.get_H();
        _omega_m = _taum_decay.get_omega();
    }
    void set_taup_decay_momenta(std::vector<vd_t> p_list) {
        _taup_decay.set_momenta(p_list);
        _H_p = _taup_decay.get_H();
        _omega_p = _taup_decay.get_omega();
    }
    // * Suppose the BSM contribution is proportional to some parameter c1;
    // * ME2 = ME2_SM + c1 * ME2_Interference + c1^2 ME2_BSM;
    double get_ME2_SM() const { return _ME2[0]; }
    double get_ME2_Interference() const { return _ME2[1]; }
    double get_ME2_BSM() const { return _ME2[2]; }

protected:
    TAUM_DECAY_t _taum_decay;
    TAUP_DECAY_t _taup_decay;

    vd_t _H_m;
    vd_t _H_p;

    double _omega_m;
    double _omega_p;

    double _ME2[3];
};

template <class TAUM_DECAY_t, class TAUP_DECAY_t>
class ME_EDM_Re_t : public ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t> {
    // * For e+ e- > gamma* > tau+ tau-
    // * Including EDM operator between gamma and tau lepton
public:
    ME_EDM_Re_t() {}
    ~ME_EDM_Re_t() {}

    virtual void set_production_momenta(std::vector<vd_t> p_list) override {
        vd_t _p_ep = p_list[0];
        vd_t _p_em = p_list[1];
        vd_t _p_taup = p_list[2];
        vd_t _p_taum = p_list[3];
        vd_t _p_s = _p_ep + _p_em;
        vd_t _H_m = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_H_m;
        vd_t _H_p = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_H_p;

        double SP56 = _H_m * _H_p;
        double SP50 = _H_m * _p_em;
        double SP51 = _H_m * _p_ep;
        double SP54 = _H_m * _p_taup;
        double SP52 = _H_m * _p_s;
        double SP60 = _H_p * _p_em;
        double SP61 = _H_p * _p_ep;
        double SP63 = _H_p * _p_taum;
        double SP62 = _H_p * _p_s;
        double SP01 = _p_em * _p_ep;
        double SP03 = _p_em * _p_taum;
        double SP04 = _p_em * _p_taup;
        double SP02 = _p_em * _p_s;
        double SP13 = _p_ep * _p_taum;
        double SP14 = _p_ep * _p_taup;
        double SP12 = _p_ep * _p_s;
        double SP34 = _p_taum * _p_taup;
        double SP32 = _p_taum * _p_s;
        double SP42 = _p_taup * _p_s;
        double SP22 = _p_s * _p_s;
        // std::cout << "q.q = " << SP22 << std::endl;

        double EP5603 = epsilon(_H_m, _H_p, _p_em, _p_taum);
        double EP5602 = epsilon(_H_m, _H_p, _p_em, _p_s);
        double EP5613 = epsilon(_H_m, _H_p, _p_ep, _p_taum);
        double EP5612 = epsilon(_H_m, _H_p, _p_ep, _p_s);
        double EP5642 = epsilon(_H_m, _H_p, _p_taup, _p_s);
        double EP5032 = epsilon(_H_m, _p_em, _p_taum, _p_s);
        double EP5042 = epsilon(_H_m, _p_em, _p_taup, _p_s);
        double EP5132 = epsilon(_H_m, _p_ep, _p_taum, _p_s);
        double EP5142 = epsilon(_H_m, _p_ep, _p_taup, _p_s);
        double EP6032 = epsilon(_H_p, _p_em, _p_taum, _p_s);
        double EP6042 = epsilon(_H_p, _p_em, _p_taup, _p_s);
        double EP6132 = epsilon(_H_p, _p_ep, _p_taum, _p_s);
        double EP6142 = epsilon(_H_p, _p_ep, _p_taup, _p_s);

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

        double M01 = 8 * pow(MTAU, 2) * SP01 + 8 * (SP04 * SP13 + SP03 * SP14);
        double M02 = 0;
        double M03 = -8 * SP04 * SP13 * SP22 - 8 * SP03 * SP14 * SP22 -
                     4 * pow(MTAU, 2) * (2 * SP02 * SP12 + SP01 * SP22) + 8 * SP04 * SP12 * SP32 +
                     8 * SP02 * SP14 * SP32 - 8 * SP02 * SP12 * SP34 + 4 * SP01 * SP22 * SP34 + 8 * SP03 * SP12 * SP42 +
                     8 * SP02 * SP13 * SP42;
        double M11 = 0;
        double M12 = 0;
        double M13 = 0;
        double M21 = 0;
        double M22 = 0;
        double M23 = 0;
        double M31 = -8 * (SP04 * SP13 * SP56 + SP03 * SP14 * SP56 - SP01 * SP34 * SP56 + pow(MTAU, 2) * SP51 * SP60 +
                           SP34 * SP51 * SP60 - SP13 * SP54 * SP60 + pow(MTAU, 2) * SP50 * SP61 + SP34 * SP50 * SP61 -
                           SP03 * SP54 * SP61 - SP14 * SP50 * SP63 - SP04 * SP51 * SP63 + SP01 * SP54 * SP63);
        double M32 =
            -4 * MTAU *
            (-2 * EP5642 * SP01 - EP5613 * SP02 - EP5612 * SP03 + 2 * EP5612 * SP04 - EP5603 * SP12 - EP5602 * SP13 +
             2 * EP5602 * SP14 + EP6132 * SP50 + 2 * EP6142 * SP50 + EP6032 * SP51 + 2 * EP6042 * SP51 +
             3 * EP5132 * SP60 + 2 * EP5142 * SP60 + 3 * EP5032 * SP61 + 2 * EP5042 * SP61);
        double M33 =
            4 * (-2 * SP03 * SP14 * SP22 * SP56 + 2 * SP02 * SP14 * SP32 * SP56 - 2 * SP02 * SP12 * SP34 * SP56 +
                 3 * SP01 * SP22 * SP34 * SP56 + 2 * SP03 * SP12 * SP42 * SP56 + 2 * SP02 * SP13 * SP42 * SP56 -
                 4 * SP01 * SP32 * SP42 * SP56 - 2 * SP22 * SP34 * SP51 * SP60 + 4 * SP32 * SP42 * SP51 * SP60 +
                 2 * SP12 * SP34 * SP52 * SP60 - 4 * SP13 * SP42 * SP52 * SP60 + 2 * SP13 * SP22 * SP54 * SP60 -
                 2 * SP12 * SP32 * SP54 * SP60 - 2 * SP22 * SP34 * SP50 * SP61 + 4 * SP32 * SP42 * SP50 * SP61 +
                 2 * SP02 * SP34 * SP52 * SP61 - 4 * SP03 * SP42 * SP52 * SP61 + 2 * SP03 * SP22 * SP54 * SP61 -
                 2 * SP02 * SP32 * SP54 * SP61 - 4 * SP14 * SP32 * SP50 * SP62 + 2 * SP12 * SP34 * SP50 * SP62 +
                 2 * SP02 * SP34 * SP51 * SP62 + 4 * SP03 * SP14 * SP52 * SP62 - 4 * SP01 * SP34 * SP52 * SP62 -
                 2 * SP03 * SP12 * SP54 * SP62 - 2 * SP02 * SP13 * SP54 * SP62 + 4 * SP01 * SP32 * SP54 * SP62 +
                 pow(MTAU, 2) *
                     (SP01 * SP22 * SP56 - 2 * SP22 * SP51 * SP60 + 2 * SP12 * SP52 * SP60 - 2 * SP22 * SP50 * SP61 +
                      2 * SP12 * SP50 * SP62 + SP02 * (-2 * SP12 * SP56 + 2 * SP52 * SP61 + 2 * SP51 * SP62)) +
                 2 * SP14 * SP22 * SP50 * SP63 - 2 * SP12 * SP42 * SP50 * SP63 - 2 * SP02 * SP42 * SP51 * SP63 -
                 2 * SP02 * SP14 * SP52 * SP63 + 4 * SP01 * SP42 * SP52 * SP63 + 2 * SP02 * SP12 * SP54 * SP63 -
                 3 * SP01 * SP22 * SP54 * SP63 -
                 2 * SP04 *
                     (SP13 * SP22 * SP56 - SP12 * SP32 * SP56 + 2 * SP32 * SP51 * SP62 - 2 * SP13 * SP52 * SP62 -
                      SP22 * SP51 * SP63 + SP12 * SP52 * SP63));

        double _omega_m = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_omega_m;
        double _omega_p = ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_omega_p;
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_ME2[0] = M01 * _omega_m * _omega_p - M11 - M21 + M31;
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_ME2[1] = M02 * _omega_m * _omega_p - M12 - M22 + M32;
        ME_Base_t<TAUM_DECAY_t, TAUP_DECAY_t>::_ME2[2] = M03 * _omega_m * _omega_p - M13 - M23 + M33;
    }
};

#endif  // TAU_OO_MATRIX_ELEMENTS_H_
