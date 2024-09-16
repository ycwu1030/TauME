#ifndef TAU_OO_TAUDECAY_H_
#define TAU_OO_TAUDECAY_H_

#include <vector>

#include "constants.h"
#include "utilities.h"

template <bool is_anti = false>
class TauDecay_t {
public:
    TauDecay_t() : _gamma_va(1.0) { _sign = is_anti ? -1 : 1; }
    // * is_anti: false for taum;
    // * is_anti: true for taup;
    virtual ~TauDecay_t() {}

    virtual void set_momenta(std::vector<vd_t> p_list) {
        _set_momenta(p_list);
        _calc_Pi();
        _calc_omgega_H();
    }
    double get_omega() const { return _omega; }
    vd_t get_H() const { return _H; }

protected:
    virtual void _set_momenta(
        std::vector<vd_t> p_list) = 0;  // * Shall be implemented in each mode, where the J will be calculated

    virtual void _calc_Pi() {
        vcd_t _Jc = conjugate(_J);
        vcd_t tmp_Pi = 2.0 * ((_Jc * _p_nu_tau) * _J + (_J * _p_nu_tau) * _Jc - (_Jc * _J) * _p_nu_tau);
        _Pi = real(tmp_Pi);

        _Pi5 = 2.0 * imag(epsilon(_Jc, _J, _p_nu_tau));
    }
    virtual void _calc_omgega_H() {
        // _omega = (_Pi + _gamma_va * _Pi5) * _p_tau;
        // vd_t tmp = -1.0 * _Pi5 - _gamma_va * _Pi;
        // _H = _sign * 1.0 * (MTAU * tmp - (_p_tau * tmp) * _p_tau / MTAU);
        _omega = _Pi * _p_tau;
        _H = _sign * 1.0 * (-1.0 * MTAU * _gamma_va * _Pi + _gamma_va * (_Pi * _p_tau) * _p_tau / MTAU);
    }

    double _gamma_va;

    vd_t _p_tau;
    vd_t _p_nu_tau;

    int _sign;

    vcd_t _J;

    double _omega;
    vd_t _H;

    vd_t _Pi;
    vd_t _Pi5;
};

template <bool is_anti = false>
class TauDecay_pi : public TauDecay_t<is_anti> {
public:
    TauDecay_pi() {}
    ~TauDecay_pi() {}

protected:
    virtual void _set_momenta(std::vector<vd_t> p_list) override {
        TauDecay_t<is_anti>::_p_nu_tau = p_list[0];
        _p_pi = p_list[1];
        TauDecay_t<is_anti>::_p_tau = TauDecay_t<is_anti>::_p_nu_tau + _p_pi;
        TauDecay_t<is_anti>::_J = tocomplex(_p_pi);
    }

private:
    vd_t _p_pi;
};

template <bool is_anti = false>
class TauDecay_rho : public TauDecay_t<is_anti> {
public:
    TauDecay_rho() {}
    ~TauDecay_rho() {}

protected:
    virtual void _set_momenta(std::vector<vd_t> p_list) override {
        TauDecay_t<is_anti>::_p_nu_tau = p_list[0];
        _p_pi0 = p_list[1];
        _p_pic = p_list[2];
        TauDecay_t<is_anti>::_p_tau = TauDecay_t<is_anti>::_p_nu_tau + _p_pi0 + _p_pic;
        TauDecay_t<is_anti>::_J = tocomplex(_p_pic - _p_pi0);
    }

private:
    vd_t _p_pi0;
    vd_t _p_pic;
};

#endif  // TAU_OO_TAUDECAY_FINALSTATES_H_
