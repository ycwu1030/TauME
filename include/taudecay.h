#ifndef TAU_OO_TAUDECAY_H_
#define TAU_OO_TAUDECAY_H_

#include <vector>

#include "constants.h"
#include "utilities.h"

template <bool is_anti = false>
class TauDecay_t {
public:
    TauDecay_t() : _gamma_va(1.0) { _sign = is_anti ? -1.0 : 1.0; }
    // * is_anti: false for taum;
    // * is_anti: true for taup;
    virtual ~TauDecay_t() {}

    virtual void set_momenta(std::vector<dlv_t> p_list) {
        _set_momenta(p_list);
        _calc_Pi();
        _calc_omgega_H();
    }
    cd_t get_omega() const { return _omega; }
    clv_t get_H() const { return _H; }

protected:
    virtual void _set_momenta(
        std::vector<dlv_t> p_list) = 0;  // * Shall be implemented in each mode, where the J will be calculated

    virtual void _calc_Pi() {
        clv_t _Jc = conjugate(_J);
        clv_t tmp_Pi = 2.0 * ((_Jc.Dot(_p_nu_tau)) * _J + (_J.Dot(_p_nu_tau)) * _Jc - (_Jc.Dot(_J)) * _p_nu_tau);
        _Pi = TrimImaginary(tmp_Pi);
        _Pi5 = ToComplex(2.0 * imag(epsilon(_Jc, _J, _p_nu_tau)));
    }
    virtual void _calc_omgega_H() {
        _omega = (_Pi + _gamma_va * _Pi5).Dot(_p_tau);
        clv_t tmp = -_Pi5 - _gamma_va * _Pi;
        _H = _sign * (MTAU * tmp - (_p_tau.Dot(tmp)) * _p_tau * (1.0 / MTAU));
        // _omega = _Pi * _p_tau;
        // _H = _sign * 1.0 * (-1.0 * MTAU * _gamma_va * _Pi + _gamma_va * (_Pi * _p_tau) * _p_tau / MTAU);
    }

    cd_t _gamma_va;

    clv_t _p_tau;
    clv_t _p_nu_tau;

    cd_t _sign;

    clv_t _J;

    cd_t _omega;
    clv_t _H;

    clv_t _Pi;
    clv_t _Pi5;
};

template <bool is_anti = false>
class TauDecay_pi : public TauDecay_t<is_anti> {
public:
    TauDecay_pi() {}
    ~TauDecay_pi() {}

protected:
    virtual void _set_momenta(std::vector<dlv_t> p_list) override {
        TauDecay_t<is_anti>::_p_nu_tau = ToComplex(p_list[0]);
        _p_pi = ToComplex(p_list[1]);
        TauDecay_t<is_anti>::_p_tau = TauDecay_t<is_anti>::_p_nu_tau + _p_pi;
        TauDecay_t<is_anti>::_J = _p_pi;
    }

private:
    clv_t _p_pi;
};

template <bool is_anti = false>
class TauDecay_rho : public TauDecay_t<is_anti> {
public:
    TauDecay_rho() {}
    ~TauDecay_rho() {}

protected:
    virtual void _set_momenta(std::vector<dlv_t> p_list) override {
        TauDecay_t<is_anti>::_p_nu_tau = ToComplex(p_list[0]);
        _p_pi0 = ToComplex(p_list[1]);
        _p_pic = ToComplex(p_list[2]);
        TauDecay_t<is_anti>::_p_tau = TauDecay_t<is_anti>::_p_nu_tau + _p_pi0 + _p_pic;
        TauDecay_t<is_anti>::_J = _p_pic - _p_pi0;
    }

private:
    clv_t _p_pi0;
    clv_t _p_pic;
};

namespace {
inline double _lam_Kallen(double x, double y, double z) {
    return x * x + y * y + z * z - 2.0 * x * y - 2.0 * y * z - 2.0 * z * x;
}

inline double _Pi_2body_decay(double m0, double mi, double mj) {
    double m02 = m0 * m0;
    double mi2 = mi * mi;
    double mj2 = mj * mj;
    double lam = _lam_Kallen(1.0, mi2 / m02, mj2 / m02);
    return m0 / 2.0 * sqrt(lam);
}

inline double _Gamma_Rho(double Q2, double mass, double width) {
    double Q = sqrt(Q2);
    double P_pi_on_shell =
        _Pi_2body_decay(mass, MPI0, MPI0);  // * ignore the difference of the mass between pi+- and pi0
    double P_pi_off_shell = _Pi_2body_decay(Q, MPI0, MPI0);
    return width * mass / Q * pow(P_pi_off_shell / P_pi_on_shell, 3);
}

inline cd_t _BW_RHO(double Q2, double mass, double width) {
    cd_t I{0.0, 1.0};
    return mass * mass / (mass * mass - Q2 - I * mass * _Gamma_Rho(Q2, mass, width));
}

inline cd_t _F_PI(double Q2) {
    cd_t BW_RHO = _BW_RHO(Q2, MRHO, GAMMARHO);
    cd_t BW_RHOP = _BW_RHO(Q2, MRHOP, GAMMARHOP);
    return (BW_RHO + BETA * BW_RHOP) / (1.0 + BETA);
}
}  // namespace

template <bool is_anti = false>
class TauDecay_a1 : public TauDecay_t<is_anti> {
public:
    TauDecay_a1() {}
    ~TauDecay_a1() {}

protected:
    virtual void _set_momenta(std::vector<dlv_t> p_list) override {
        TauDecay_t<is_anti>::_p_nu_tau = ToComplex(p_list[0]);
        _p_pi1 = ToComplex(p_list[1]);
        _p_pi2 = ToComplex(p_list[2]);
        _p_pi3 = ToComplex(p_list[3]);
        _p_Q = _p_pi1 + _p_pi2 + _p_pi3;
        _p_q13 = _p_pi1 + _p_pi3;
        _p_q23 = _p_pi2 + _p_pi3;

        double Q2 = real(_p_Q.Dot(_p_Q));

        clv_t _q13 = _p_pi1 - _p_pi3;
        double Q13sq = real(_p_q13.Dot(_p_q13));
        cd_t F13 = _F_PI(Q13sq);

        clv_t _q23 = _p_pi2 - _p_pi3;
        double Q23sq = real(_p_q23.Dot(_p_q23));
        cd_t F23 = _F_PI(Q23sq);

        TauDecay_t<is_anti>::_p_tau = TauDecay_t<is_anti>::_p_nu_tau + _p_pi1 + _p_pi2 + _p_pi3;
        TauDecay_t<is_anti>::_J =
            (_q13 - (_p_Q.Dot(_q13) / Q2) * _p_Q) * F13 + (_q23 - (_p_Q.Dot(_q23) / Q2) * _p_Q) * F23;
        // * When calculating J, I ignore some factor that won't affect the final results
    }

private:
    // * tau -> nu pi1 pi2 pi3
    // * with either
    // * tau -> nu pi- pi- pi+
    // * or
    // * tau -> nu pi0 pi0 pi-
    clv_t _p_pi1;
    clv_t _p_pi2;
    clv_t _p_pi3;
    clv_t _p_q13;
    clv_t _p_q23;
    clv_t _p_Q;
};

#endif  // TAU_OO_TAUDECAY_FINALSTATES_H_
