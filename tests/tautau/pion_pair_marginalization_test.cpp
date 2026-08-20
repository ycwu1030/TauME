#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "tauamp/tautau/pion_pair_marginalization.h"

namespace {
constexpr double kTolerance = 1e-14;

bool close(double left, double right) { return std::abs(left - right) < kTolerance; }

template <class Callback>
void expect_invalid_argument(Callback callback) {
    bool thrown = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        thrown = true;
    }
    assert(thrown);
}
}  // namespace

int main() {
    using tauamp::tautau::LinearComponents;
    using tauamp::tautau::WeightedLinearComponents;

    const LinearComponents first{2.0, 4.0, -2.0, 2.0, -6.0};
    const LinearComponents second{8.0, 8.0, 16.0, 0.0, 24.0};
    const std::vector<WeightedLinearComponents> components{{first, 1.0}, {second, 3.0}};

    const auto average_of_ratios = tauamp::tautau::average_conditional_observables(components);
    const auto ratio_after_average = tauamp::tautau::ratio_of_average_components(components);
    assert(close(average_of_ratios.f2_real, 1.25));
    assert(close(average_of_ratios.f2_imaginary, 1.25));
    assert(close(average_of_ratios.f3_real, 0.25));
    assert(close(average_of_ratios.f3_imaginary, 1.5));
    assert(close(ratio_after_average.f2_real, 28.0 / 26.0));
    assert(close(ratio_after_average.f2_imaginary, 46.0 / 26.0));
    assert(close(ratio_after_average.f3_real, 2.0 / 26.0));
    assert(close(ratio_after_average.f3_imaginary, 66.0 / 26.0));
    assert(!close(average_of_ratios.f3_real, ratio_after_average.f3_real));

    const std::vector<WeightedLinearComponents> one_candidate{{first, 2.0}};
    const auto conditional = tauamp::tautau::conditional_observables(first);
    const auto one_average_of_ratios = tauamp::tautau::average_conditional_observables(one_candidate);
    const auto one_ratio_after_average = tauamp::tautau::ratio_of_average_components(one_candidate);
    assert(close(conditional.f3_real, one_average_of_ratios.f3_real));
    assert(close(conditional.f3_real, one_ratio_after_average.f3_real));

    expect_invalid_argument([&] { tauamp::tautau::conditional_observables(LinearComponents{}); });
    expect_invalid_argument([&] { tauamp::tautau::average_conditional_observables({}); });
    expect_invalid_argument([&] {
        tauamp::tautau::average_conditional_observables({WeightedLinearComponents{first, -1.0}});
    });
    expect_invalid_argument([&] {
        tauamp::tautau::ratio_of_average_components(
            {WeightedLinearComponents{first, std::numeric_limits<double>::infinity()}});
    });
    expect_invalid_argument([&] {
        tauamp::tautau::average_conditional_observables({WeightedLinearComponents{first, 0.0}});
    });
    expect_invalid_argument([&] {
        tauamp::tautau::ratio_of_average_components({WeightedLinearComponents{LinearComponents{1.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
                                                      WeightedLinearComponents{LinearComponents{-1.0, 0.0, 0.0, 0.0, 0.0}, 1.0}});
    });
}
