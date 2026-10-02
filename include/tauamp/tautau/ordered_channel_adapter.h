#ifndef TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_
#define TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "tauamp/tautau/four_momentum.h"
#include "tauamp/tautau/rho_a1_currents.h"

namespace tauamp::tautau {

enum class OrderedChannel { pi_pi, pi_rho, rho_pi, rho_rho, pi_a1, a1_pi, rho_a1, a1_rho, a1_a1 };

inline const char* ordered_channel_name(OrderedChannel channel) {
    switch (channel) {
        case OrderedChannel::pi_pi: return "pi_pi";
        case OrderedChannel::pi_rho: return "pi_rho";
        case OrderedChannel::rho_pi: return "rho_pi";
        case OrderedChannel::rho_rho: return "rho_rho";
        case OrderedChannel::pi_a1: return "pi_a1";
        case OrderedChannel::a1_pi: return "a1_pi";
        case OrderedChannel::rho_a1: return "rho_a1";
        case OrderedChannel::a1_rho: return "a1_rho";
        case OrderedChannel::a1_a1: return "a1_a1";
    }
    throw std::invalid_argument("invalid ordered channel");
}

struct ReconstructedPion {
    std::size_t index{};
    int pid{};
    double pt{};
    double eta{};
    double phi{};
    FourMomentum momentum{};

    int charge() const {
        if (pid == 211) return 1;
        if (pid == -211) return -1;
        if (pid == 111) return 0;
        throw std::invalid_argument("G6 channel adapter received an unsupported pion PID");
    }
};

struct OrderedDecaySide {
    std::string mode;
    std::vector<std::size_t> charged_indices;
    std::vector<std::size_t> neutral_indices;
    std::vector<std::size_t> same_charge_slots;
    std::size_t opposite_charge_slot{};
    bool has_opposite_charge_slot{false};
};

struct OrderedChannelHypothesis {
    OrderedChannel channel{};
    std::size_t hypothesis_index{};
    OrderedDecaySide tau_plus;
    OrderedDecaySide tau_minus;
    std::vector<double> rho_masses;
    std::vector<double> a1_masses;
    double normalized_mass_cost{0.0};
    std::size_t unused_charged_count{0};
    std::size_t unused_neutral_count{0};
};

namespace detail {

inline std::pair<std::vector<ReconstructedPion>, std::vector<ReconstructedPion>> split_charges(
    const std::vector<ReconstructedPion>& charged) {
    std::vector<ReconstructedPion> plus;
    std::vector<ReconstructedPion> minus;
    for (const auto& object : charged) {
        if (object.charge() == 1) plus.push_back(object);
        else if (object.charge() == -1) minus.push_back(object);
        else throw std::invalid_argument("charged collection contains a neutral pion");
    }
    const auto by_index = [](const ReconstructedPion& left, const ReconstructedPion& right) {
        return left.index < right.index;
    };
    std::sort(plus.begin(), plus.end(), by_index);
    std::sort(minus.begin(), minus.end(), by_index);
    return {plus, minus};
}

inline std::vector<ReconstructedPion> all_neutrals(const std::vector<ReconstructedPion>& neutral) {
    for (const auto& object : neutral) {
        if (object.charge() != 0) throw std::invalid_argument("neutral collection contains a charged pion");
    }
    std::vector<ReconstructedPion> result = neutral;
    std::sort(result.begin(), result.end(), [](const ReconstructedPion& left, const ReconstructedPion& right) {
        return left.index < right.index;
    });
    return result;
}

template <typename T, typename Callback>
inline void choose_subsets(const std::vector<T>& objects, std::size_t required, Callback callback) {
    if (required > objects.size()) return;
    std::vector<T> selected;
    std::function<void(std::size_t)> visit = [&](std::size_t begin) {
        if (selected.size() == required) {
            callback(selected);
            return;
        }
        const std::size_t remaining = required - selected.size();
        for (std::size_t index = begin; index + remaining <= objects.size(); ++index) {
            selected.push_back(objects[index]);
            visit(index + 1U);
            selected.pop_back();
        }
    };
    visit(0U);
}

inline double pair_mass(const ReconstructedPion& charged, const ReconstructedPion& neutral) {
    return std::sqrt(std::max((charged.momentum + neutral.momentum).mass_squared(), 0.0));
}

inline constexpr double g9_rho_window_multiplier = 3.0;
inline constexpr double g9_a1_mass = 1.23;
inline constexpr double g9_a1_width = 0.42;
inline constexpr double g9_a1_window_multiplier = 3.0;

inline bool rho_pass(double mass) {
    return mass + 1.0e-12 >= g6_charged_pion_mass + g6_neutral_pion_mass &&
           std::abs(mass - g6_rho_mass) < g9_rho_window_multiplier * g6_rho_width;
}

inline bool a1_pass(double mass) {
    return mass + 1.0e-12 >= 3.0 * g6_charged_pion_mass &&
           std::abs(mass - g9_a1_mass) < g9_a1_window_multiplier * g9_a1_width;
}

inline double rho_mass_cost(double mass) {
    const double z = (mass - g6_rho_mass) / g6_rho_width;
    return z * z;
}

inline double a1_mass_cost(double mass) {
    const double z = (mass - g9_a1_mass) / g9_a1_width;
    return z * z;
}

inline double delta_phi(double left, double right) {
    constexpr double pi = 3.14159265358979323846;
    double difference = std::fmod(left - right + pi, 2.0 * pi);
    if (difference < 0.0) difference += 2.0 * pi;
    return difference - pi;
}

inline double delta_r(const ReconstructedPion& left, const ReconstructedPion& right) {
    return std::hypot(left.eta - right.eta, delta_phi(left.phi, right.phi));
}

inline OrderedDecaySide pion_side(const ReconstructedPion& pion) {
    return {"pi", {pion.index}, {}, {}, 0U, false};
}

inline OrderedDecaySide rho_side(const ReconstructedPion& charged, const ReconstructedPion& neutral) {
    return {"rho", {charged.index}, {neutral.index}, {}, 0U, false};
}

inline OrderedDecaySide a1_side(const std::vector<ReconstructedPion>& same_charge,
                                const ReconstructedPion& opposite) {
    std::vector<ReconstructedPion> ordered = same_charge;
    std::sort(ordered.begin(), ordered.end(), [](const ReconstructedPion& left, const ReconstructedPion& right) {
        if (left.pt != right.pt) return left.pt > right.pt;
        return left.index < right.index;
    });
    std::vector<std::size_t> slots;
    for (const auto& object : ordered) slots.push_back(object.index);
    std::vector<std::size_t> charged = slots;
    charged.push_back(opposite.index);
    return {"a1", charged, {}, slots, opposite.index, true};
}

}  // namespace detail

inline std::vector<OrderedChannelHypothesis> enumerate_ordered_channel_hypotheses(
    OrderedChannel channel, const std::vector<ReconstructedPion>& charged,
    const std::vector<ReconstructedPion>& neutral) {
    const auto charge_split = detail::split_charges(charged);
    const auto& plus = charge_split.first;
    const auto& minus = charge_split.second;
    const auto neutrals = detail::all_neutrals(neutral);
    std::vector<OrderedChannelHypothesis> result;
    const auto add = [&](OrderedDecaySide tau_plus, OrderedDecaySide tau_minus,
                         std::vector<double> rho_masses = {}, std::vector<double> a1_masses = {}) {
        const double mass_cost =
            std::accumulate(rho_masses.begin(), rho_masses.end(), 0.0, [](double sum, double mass) {
                return sum + detail::rho_mass_cost(mass);
            }) +
            std::accumulate(a1_masses.begin(), a1_masses.end(), 0.0, [](double sum, double mass) {
                return sum + detail::a1_mass_cost(mass);
            });
        const std::size_t slot_count = 2U + rho_masses.size() + a1_masses.size();
        std::vector<std::size_t> used_charged;
        used_charged.insert(used_charged.end(), tau_plus.charged_indices.begin(), tau_plus.charged_indices.end());
        used_charged.insert(used_charged.end(), tau_minus.charged_indices.begin(), tau_minus.charged_indices.end());
        std::vector<std::size_t> used_neutral;
        used_neutral.insert(used_neutral.end(), tau_plus.neutral_indices.begin(), tau_plus.neutral_indices.end());
        used_neutral.insert(used_neutral.end(), tau_minus.neutral_indices.begin(), tau_minus.neutral_indices.end());
        const auto is_used = [](const std::vector<std::size_t>& used, std::size_t index) {
            return std::find(used.begin(), used.end(), index) != used.end();
        };
        std::size_t unused_charged = 0U;
        for (const auto& object : charged) if (!is_used(used_charged, object.index)) ++unused_charged;
        std::size_t unused_neutral = 0U;
        for (const auto& object : neutral) if (!is_used(used_neutral, object.index)) ++unused_neutral;
        result.push_back({channel, result.size(), std::move(tau_plus), std::move(tau_minus),
                          std::move(rho_masses), std::move(a1_masses),
                          mass_cost / static_cast<double>(slot_count), unused_charged, unused_neutral});
    };

    if (channel == OrderedChannel::pi_pi) {
        for (const auto& plus_pion : plus)
            for (const auto& minus_pion : minus)
                add(detail::pion_side(plus_pion), detail::pion_side(minus_pion));
    } else if (channel == OrderedChannel::pi_rho || channel == OrderedChannel::rho_pi) {
        for (const auto& plus_pion : plus) for (const auto& minus_pion : minus)
            for (const auto& pi0 : neutrals) {
                const auto& rho_charged = channel == OrderedChannel::pi_rho ? minus_pion : plus_pion;
                const double mass = detail::pair_mass(rho_charged, pi0);
                if (!detail::rho_pass(mass)) continue;
                if (channel == OrderedChannel::pi_rho)
                    add(detail::pion_side(plus_pion), detail::rho_side(minus_pion, pi0), {mass});
                else
                    add(detail::rho_side(plus_pion, pi0), detail::pion_side(minus_pion), {mass});
            }
    } else if (channel == OrderedChannel::rho_rho) {
        for (const auto& plus_pion : plus) for (const auto& minus_pion : minus)
            for (std::size_t first = 0; first < neutrals.size(); ++first)
                for (std::size_t second = 0; second < neutrals.size(); ++second) {
                    if (first == second) continue;
                    const double plus_mass = detail::pair_mass(plus_pion, neutrals[first]);
                    const double minus_mass = detail::pair_mass(minus_pion, neutrals[second]);
                    if (detail::rho_pass(plus_mass) && detail::rho_pass(minus_mass))
                        add(detail::rho_side(plus_pion, neutrals[first]),
                            detail::rho_side(minus_pion, neutrals[second]), {plus_mass, minus_mass});
                }
    } else if (channel == OrderedChannel::pi_a1) {
        for (const auto& direct : plus) {
            std::vector<ReconstructedPion> remaining_plus;
            for (const auto& object : plus) if (object.index != direct.index) remaining_plus.push_back(object);
            detail::choose_subsets(minus, 2U, [&](const auto& same_minus) {
                for (const auto& opposite : remaining_plus) {
                    const double mass = std::sqrt(std::max(
                        (same_minus[0].momentum + same_minus[1].momentum + opposite.momentum).mass_squared(), 0.0));
                    if (detail::a1_pass(mass))
                        add(detail::pion_side(direct), detail::a1_side(same_minus, opposite), {}, {mass});
                }
            });
        }
    } else if (channel == OrderedChannel::a1_pi) {
        for (const auto& direct : minus) {
            std::vector<ReconstructedPion> remaining_minus;
            for (const auto& object : minus) if (object.index != direct.index) remaining_minus.push_back(object);
            detail::choose_subsets(plus, 2U, [&](const auto& same_plus) {
                for (const auto& opposite : remaining_minus) {
                    const double mass = std::sqrt(std::max(
                        (same_plus[0].momentum + same_plus[1].momentum + opposite.momentum).mass_squared(), 0.0));
                    if (detail::a1_pass(mass))
                        add(detail::a1_side(same_plus, opposite), detail::pion_side(direct), {}, {mass});
                }
            });
        }
    } else if (channel == OrderedChannel::rho_a1 || channel == OrderedChannel::a1_rho) {
        if (channel == OrderedChannel::rho_a1) {
            for (const auto& rho_charged : plus) for (const auto& pi0 : neutrals) {
                const double rho_mass = detail::pair_mass(rho_charged, pi0);
                if (!detail::rho_pass(rho_mass)) continue;
                std::vector<ReconstructedPion> remaining_plus;
                for (const auto& object : plus) if (object.index != rho_charged.index) remaining_plus.push_back(object);
                detail::choose_subsets(minus, 2U, [&](const auto& same_minus) {
                    for (const auto& opposite : remaining_plus) {
                        const double a1_mass = std::sqrt(std::max(
                            (same_minus[0].momentum + same_minus[1].momentum + opposite.momentum).mass_squared(), 0.0));
                        if (detail::a1_pass(a1_mass))
                            add(detail::rho_side(rho_charged, pi0), detail::a1_side(same_minus, opposite),
                                {rho_mass}, {a1_mass});
                    }
                });
            }
        } else {
            for (const auto& rho_charged : minus) for (const auto& pi0 : neutrals) {
                const double rho_mass = detail::pair_mass(rho_charged, pi0);
                if (!detail::rho_pass(rho_mass)) continue;
                std::vector<ReconstructedPion> remaining_minus;
                for (const auto& object : minus) if (object.index != rho_charged.index) remaining_minus.push_back(object);
                detail::choose_subsets(plus, 2U, [&](const auto& same_plus) {
                    for (const auto& opposite : remaining_minus) {
                        const double a1_mass = std::sqrt(std::max(
                            (same_plus[0].momentum + same_plus[1].momentum + opposite.momentum).mass_squared(), 0.0));
                        if (detail::a1_pass(a1_mass))
                            add(detail::a1_side(same_plus, opposite), detail::rho_side(rho_charged, pi0),
                                {rho_mass}, {a1_mass});
                    }
                });
            }
        }
    } else if (channel == OrderedChannel::a1_a1) {
        detail::choose_subsets(plus, 3U, [&](const auto& plus_triplet) {
            detail::choose_subsets(minus, 3U, [&](const auto& minus_triplet) {
                for (const auto& plus_opposite : plus_triplet) for (const auto& minus_opposite : minus_triplet) {
                    std::vector<ReconstructedPion> plus_same, minus_same;
                    for (const auto& object : plus_triplet) if (object.index != plus_opposite.index) plus_same.push_back(object);
                    for (const auto& object : minus_triplet) if (object.index != minus_opposite.index) minus_same.push_back(object);
                    const double plus_mass = std::sqrt(std::max(
                        (plus_same[0].momentum + plus_same[1].momentum + minus_opposite.momentum).mass_squared(), 0.0));
                    const double minus_mass = std::sqrt(std::max(
                        (minus_same[0].momentum + minus_same[1].momentum + plus_opposite.momentum).mass_squared(), 0.0));
                    if (detail::a1_pass(plus_mass) && detail::a1_pass(minus_mass))
                        add(detail::a1_side(plus_same, minus_opposite), detail::a1_side(minus_same, plus_opposite),
                            {}, {plus_mass, minus_mass});
                }
            });
        });
    }
    return result;
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_
