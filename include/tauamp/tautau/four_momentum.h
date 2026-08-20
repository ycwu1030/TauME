#ifndef TAUAMP_TAUTAU_FOUR_MOMENTUM_H_
#define TAUAMP_TAUTAU_FOUR_MOMENTUM_H_

#include <array>
#include <cmath>
#include <stdexcept>

namespace tauamp::tautau {

class FourMomentum {
public:
    FourMomentum() = default;
    FourMomentum(double px, double py, double pz, double energy) : px_(px), py_(py), pz_(pz), energy_(energy) {}

    double px() const { return px_; }
    double py() const { return py_; }
    double pz() const { return pz_; }
    double energy() const { return energy_; }
    double mass_squared() const { return energy_ * energy_ - px_ * px_ - py_ * py_ - pz_ * pz_; }

    FourMomentum boosted(const std::array<double, 3>& beta) const {
        const double beta_squared = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
        if (beta_squared >= 1.0) throw std::invalid_argument("boost velocity must be subluminal");
        if (beta_squared < 1e-18) return *this;
        const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
        const double beta_dot_p = beta[0] * px_ + beta[1] * py_ + beta[2] * pz_;
        const double factor = (gamma - 1.0) * beta_dot_p / beta_squared + gamma * energy_;
        return {px_ + factor * beta[0], py_ + factor * beta[1], pz_ + factor * beta[2],
                gamma * (energy_ + beta_dot_p)};
    }

private:
    double px_{};
    double py_{};
    double pz_{};
    double energy_{};
};

inline FourMomentum operator+(const FourMomentum& left, const FourMomentum& right) {
    return {left.px() + right.px(), left.py() + right.py(), left.pz() + right.pz(),
            left.energy() + right.energy()};
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_FOUR_MOMENTUM_H_
