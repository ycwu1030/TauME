#ifndef TAUAMP_TAUTAU_TEST_LEGACY_TWOFOLD_CONTROL_H_
#define TAUAMP_TAUTAU_TEST_LEGACY_TWOFOLD_CONTROL_H_

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "tauamp/tautau/four_momentum.h"

// Test-only extraction of tauamp::EventReader_t::reconstruct_neutrinos from
// TauOO/include/tauamp/eventreader.h (legacy main, lines 218-333). It retains
// the legacy p_miss-axis rotations and quadratic construction, while replacing
// EventReader member-state writes with an explicit return vector. The original
// macro MTAU is 1.77; this control deliberately preserves that fixed value.
namespace legacy_twofold_control {

struct NeutrinoPair {
    tauamp::tautau::FourMomentum neutrino_minus;
    tauamp::tautau::FourMomentum neutrino_plus;
};

inline TLorentzVector root_momentum(const tauamp::tautau::FourMomentum& momentum) {
    return {momentum.px(), momentum.py(), momentum.pz(), momentum.energy()};
}

inline tauamp::tautau::FourMomentum core_momentum(const TLorentzVector& momentum) {
    return {momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()};
}

inline bool reconstruct_neutrinos(double Ecm, const tauamp::tautau::FourMomentum& visible_minus,
                                  const tauamp::tautau::FourMomentum& visible_plus,
                                  std::vector<NeutrinoPair>& solutions) {
    constexpr double kLegacyTauMass = 1.77;
    solutions.clear();

    // Literal algorithmic control adapted only at its input/output boundary.
    double Etau = Ecm / 2.0;
    double ptau = Ecm / 2.0 * std::sqrt(1.0 - 4.0 * kLegacyTauMass * kLegacyTauMass / Ecm / Ecm);
    TLorentzVector pcm(0, 0, 0, Ecm);
    TLorentzVector p_h_m = root_momentum(visible_minus);
    TLorentzVector p_h_p = root_momentum(visible_plus);
    TLorentzVector p_miss = pcm - (p_h_m + p_h_p);
    double pm = p_miss.P();
    double thetam = p_miss.Theta();
    double phim = p_miss.Phi();

    double Enu = Ecm / 2.0 - p_h_m.E();
    double Enubar = Ecm / 2.0 - p_h_p.E();
    if (Enu < 0 || Enubar < 0) return false;

    p_h_m.RotateZ(-phim);
    p_h_m.RotateY(-thetam);
    p_h_p.RotateZ(-phim);
    p_h_p.RotateY(-thetam);
    p_miss.RotateZ(-phim);
    p_miss.RotateY(-thetam);

    double cth = (pm * pm + Enu * Enu - Enubar * Enubar) / (2.0 * pm * Enu);
    if (std::abs(cth) > 1.0) return false;
    double theta = std::acos(cth);
    double sth = std::sin(theta);

    double z = Enu * cth;
    double r0 = Enu * sth;
    double x0 = p_h_m.Px();
    double y0 = p_h_m.Py();
    double z02 = ptau * ptau - std::pow(p_h_m.Pz() + z, 2);

    double ax = 4.0 * (x0 * x0 + y0 * y0);
    double bx = 4.0 * x0 * (r0 * r0 + x0 * x0 + y0 * y0 - z02);
    double cx = r0 * r0 * r0 * r0 + 2.0 * r0 * r0 * (x0 * x0 - y0 * y0 - z02) +
                std::pow(x0 * x0 + y0 * y0 - z02, 2);
    double Delta = bx * bx - 4.0 * ax * cx;
    if (Delta < 0) return false;
    double x_sol_1 = (-bx + std::sqrt(Delta)) / 2.0 / ax;
    double x_sol_2 = (-bx - std::sqrt(Delta)) / 2.0 / ax;

    double y_sol_1 = (z02 - r0 * r0 - y0 * y0 - x0 * x0 - 2.0 * x0 * x_sol_1) / 2.0 / y0;
    double y_sol_2 = (z02 - r0 * r0 - y0 * y0 - x0 * x0 - 2.0 * x0 * x_sol_2) / 2.0 / y0;

    TLorentzVector p_nu_sol_1;
    p_nu_sol_1.SetXYZM(x_sol_1, y_sol_1, z, 0.0);
    TLorentzVector p_nubar_sol_1 = p_miss - p_nu_sol_1;

    TLorentzVector p_nu_sol_2;
    p_nu_sol_2.SetXYZM(x_sol_2, y_sol_2, z, 0.0);
    TLorentzVector p_nubar_sol_2 = p_miss - p_nu_sol_2;

    p_h_m.RotateY(thetam);
    p_h_m.RotateZ(phim);
    p_h_p.RotateY(thetam);
    p_h_p.RotateZ(phim);

    p_nu_sol_1.RotateY(thetam);
    p_nu_sol_1.RotateZ(phim);
    p_nubar_sol_1.RotateY(thetam);
    p_nubar_sol_1.RotateZ(phim);
    p_nu_sol_2.RotateY(thetam);
    p_nu_sol_2.RotateZ(phim);
    p_nubar_sol_2.RotateY(thetam);
    p_nubar_sol_2.RotateZ(phim);

    solutions.push_back({core_momentum(p_nu_sol_1), core_momentum(p_nubar_sol_1)});
    solutions.push_back({core_momentum(p_nu_sol_2), core_momentum(p_nubar_sol_2)});
    return true;
}

}  // namespace legacy_twofold_control

#endif  // TAUAMP_TAUTAU_TEST_LEGACY_TWOFOLD_CONTROL_H_
