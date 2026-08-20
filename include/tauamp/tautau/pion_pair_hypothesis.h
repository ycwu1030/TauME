#ifndef TAUAMP_TAUTAU_PION_PAIR_HYPOTHESIS_H_
#define TAUAMP_TAUTAU_PION_PAIR_HYPOTHESIS_H_

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tauamp/tautau/kinematics.h"

namespace tauamp::tautau {

enum class KinematicPointOrigin { truth, analytic_branch, constrained_fit, sampled_trial, external };

// A latent complete pion-pair point supplied by a declared construction method.
// It is not a detector observation and does not reconstruct invisible momenta.
struct PionPairKinematicHypothesis {
    TauPairKinematicPoint point;
    FourMomentum pion_minus_lab;
    FourMomentum pion_plus_lab;
    KinematicPointOrigin origin;
};

struct WeightedPionPairKinematicHypothesis {
    PionPairKinematicHypothesis hypothesis;
    double weight;
};

class PionPairKinematicEnsemble {
public:
    explicit PionPairKinematicEnsemble(std::vector<WeightedPionPairKinematicHypothesis> entries)
        : entries_(std::move(entries)) {
        if (entries_.empty()) throw std::invalid_argument("pion-pair kinematic ensemble cannot be empty");
        for (const auto& entry : entries_) {
            if (!std::isfinite(entry.weight) || entry.weight < 0.0)
                throw std::invalid_argument("pion-pair kinematic ensemble weights must be finite and non-negative");
            total_weight_ += entry.weight;
        }
        if (!std::isfinite(total_weight_) || total_weight_ <= 0.0)
            throw std::invalid_argument("pion-pair kinematic ensemble total weight must be finite and positive");
    }

    const std::vector<WeightedPionPairKinematicHypothesis>& entries() const { return entries_; }
    double total_weight() const { return total_weight_; }

private:
    std::vector<WeightedPionPairKinematicHypothesis> entries_;
    double total_weight_{};
};

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_PION_PAIR_HYPOTHESIS_H_
