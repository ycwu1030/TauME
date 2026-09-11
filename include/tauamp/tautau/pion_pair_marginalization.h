#ifndef TAUAMP_TAUTAU_PION_PAIR_MARGINALIZATION_H_
#define TAUAMP_TAUTAU_PION_PAIR_MARGINALIZATION_H_

#include <cmath>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/spin_density.h"

namespace tauamp::tautau {

struct LinearObservableComponents {
    double f2_real{};
    double f2_imaginary{};
    double f3_real{};
    double f3_imaginary{};
};

struct WeightedLinearComponents {
    LinearComponents components;
    double weight;
};

struct WeightedPolynomialComponents {
    PolynomialComponents components;
    double weight;
};

namespace detail {
inline double checked_total_weight(const std::vector<WeightedLinearComponents>& entries) {
    if (entries.empty()) throw std::invalid_argument("linear-component ensemble cannot be empty");
    double total_weight = 0.0;
    for (const auto& entry : entries) {
        if (!std::isfinite(entry.weight) || entry.weight < 0.0)
            throw std::invalid_argument("linear-component ensemble weights must be finite and non-negative");
        total_weight += entry.weight;
    }
    if (!std::isfinite(total_weight) || total_weight <= 0.0)
        throw std::invalid_argument("linear-component ensemble total weight must be finite and positive");
    return total_weight;
}
}  // namespace detail

inline LinearObservableComponents conditional_observables(const LinearComponents& components) {
    if (components.sm == 0.0) throw std::invalid_argument("conditional observable requires a nonzero SM component");
    return {components.f2_real / components.sm, components.f2_imaginary / components.sm, components.f3_real / components.sm,
            components.f3_imaginary / components.sm};
}

// This is the average of the conditional component-to-SM ratios at supplied points.
inline LinearObservableComponents average_conditional_observables(const std::vector<WeightedLinearComponents>& entries) {
    const double total_weight = detail::checked_total_weight(entries);
    LinearObservableComponents result;
    for (const auto& entry : entries) {
        const LinearObservableComponents conditional = conditional_observables(entry.components);
        result.f2_real += entry.weight * conditional.f2_real;
        result.f2_imaginary += entry.weight * conditional.f2_imaginary;
        result.f3_real += entry.weight * conditional.f3_real;
        result.f3_imaginary += entry.weight * conditional.f3_imaginary;
    }
    result.f2_real /= total_weight;
    result.f2_imaginary /= total_weight;
    result.f3_real /= total_weight;
    result.f3_imaginary /= total_weight;
    return result;
}

// This is the ratio formed after separately averaging the stripped components.
inline LinearObservableComponents ratio_of_average_components(const std::vector<WeightedLinearComponents>& entries) {
    detail::checked_total_weight(entries);
    LinearComponents sum;
    for (const auto& entry : entries) {
        sum.sm += entry.weight * entry.components.sm;
        sum.f2_real += entry.weight * entry.components.f2_real;
        sum.f2_imaginary += entry.weight * entry.components.f2_imaginary;
        sum.f3_real += entry.weight * entry.components.f3_real;
        sum.f3_imaginary += entry.weight * entry.components.f3_imaginary;
    }
    return conditional_observables(sum);
}

inline PolynomialComponents average_polynomial_components(const std::vector<WeightedPolynomialComponents>& entries) {
    double total_weight = 0.0;
    if (entries.empty()) throw std::invalid_argument("polynomial-component ensemble cannot be empty");
    for (const auto& entry : entries) {
        if (!std::isfinite(entry.weight) || entry.weight < 0.0)
            throw std::invalid_argument("polynomial-component ensemble weights must be finite and non-negative");
        total_weight += entry.weight;
    }
    if (!std::isfinite(total_weight) || total_weight <= 0.0)
        throw std::invalid_argument("polynomial-component ensemble total weight must be finite and positive");

    PolynomialComponents result;
    for (const auto& entry : entries) {
        const double weight = entry.weight;
        const PolynomialComponents& components = entry.components;
        result.sm += weight * components.sm;
        result.f2_real += weight * components.f2_real;
        result.f2_imaginary += weight * components.f2_imaginary;
        result.f3_real += weight * components.f3_real;
        result.f3_imaginary += weight * components.f3_imaginary;
        result.f2_real_f2_real += weight * components.f2_real_f2_real;
        result.f2_real_f2_imaginary += weight * components.f2_real_f2_imaginary;
        result.f2_real_f3_real += weight * components.f2_real_f3_real;
        result.f2_real_f3_imaginary += weight * components.f2_real_f3_imaginary;
        result.f2_imaginary_f2_imaginary += weight * components.f2_imaginary_f2_imaginary;
        result.f2_imaginary_f3_real += weight * components.f2_imaginary_f3_real;
        result.f2_imaginary_f3_imaginary += weight * components.f2_imaginary_f3_imaginary;
        result.f3_real_f3_real += weight * components.f3_real_f3_real;
        result.f3_real_f3_imaginary += weight * components.f3_real_f3_imaginary;
        result.f3_imaginary_f3_imaginary += weight * components.f3_imaginary_f3_imaginary;
    }
    result.sm /= total_weight;
    result.f2_real /= total_weight;
    result.f2_imaginary /= total_weight;
    result.f3_real /= total_weight;
    result.f3_imaginary /= total_weight;
    result.f2_real_f2_real /= total_weight;
    result.f2_real_f2_imaginary /= total_weight;
    result.f2_real_f3_real /= total_weight;
    result.f2_real_f3_imaginary /= total_weight;
    result.f2_imaginary_f2_imaginary /= total_weight;
    result.f2_imaginary_f3_real /= total_weight;
    result.f2_imaginary_f3_imaginary /= total_weight;
    result.f3_real_f3_real /= total_weight;
    result.f3_real_f3_imaginary /= total_weight;
    result.f3_imaginary_f3_imaginary /= total_weight;
    return result;
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_PION_PAIR_MARGINALIZATION_H_
