#ifndef TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_
#define TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
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

inline std::vector<ReconstructedPion> hardest_neutrals(const std::vector<ReconstructedPion>& neutral,
                                                        std::size_t required) {
    for (const auto& object : neutral) {
        if (object.charge() != 0) throw std::invalid_argument("neutral collection contains a charged pion");
    }
    if (neutral.size() < required) return {};
    std::vector<ReconstructedPion> result = neutral;
    std::sort(result.begin(), result.end(), [](const ReconstructedPion& left, const ReconstructedPion& right) {
        if (left.pt != right.pt) return left.pt > right.pt;
        return left.index < right.index;
    });
    result.resize(required);
    return result;
}

inline double pair_mass(const ReconstructedPion& charged, const ReconstructedPion& neutral) {
    return std::sqrt(std::max((charged.momentum + neutral.momentum).mass_squared(), 0.0));
}

inline bool rho_pass(double mass) { return std::abs(mass - g6_rho_mass) < g6_rho_width; }

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
    const auto [plus, minus] = detail::split_charges(charged);
    std::vector<OrderedChannelHypothesis> result;
    const auto add = [&](OrderedDecaySide tau_plus, OrderedDecaySide tau_minus,
                         std::vector<double> masses = {}) {
        result.push_back({channel, result.size(), std::move(tau_plus), std::move(tau_minus), std::move(masses)});
    };

    if (channel == OrderedChannel::pi_pi) {
        if (plus.size() == 1U && minus.size() == 1U && neutral.empty())
            add(detail::pion_side(plus[0]), detail::pion_side(minus[0]));
    } else if (channel == OrderedChannel::pi_rho || channel == OrderedChannel::rho_pi) {
        const auto pi0 = detail::hardest_neutrals(neutral, 1U);
        if (plus.size() == 1U && minus.size() == 1U && pi0.size() == 1U) {
            const auto& rho_charged = channel == OrderedChannel::pi_rho ? minus[0] : plus[0];
            const double mass = detail::pair_mass(rho_charged, pi0[0]);
            if (detail::rho_pass(mass)) {
                if (channel == OrderedChannel::pi_rho)
                    add(detail::pion_side(plus[0]), detail::rho_side(minus[0], pi0[0]), {mass});
                else
                    add(detail::rho_side(plus[0], pi0[0]), detail::pion_side(minus[0]), {mass});
            }
        }
    } else if (channel == OrderedChannel::rho_rho) {
        const auto pi0 = detail::hardest_neutrals(neutral, 2U);
        if (plus.size() == 1U && minus.size() == 1U && pi0.size() == 2U) {
            struct Assignment { double cost; double sum_delta_r; std::size_t plus_index; std::size_t minus_index; double plus_mass; double minus_mass; };
            std::vector<Assignment> assignments;
            for (const auto ordering : {std::pair<std::size_t, std::size_t>{0U, 1U}, {1U, 0U}}) {
                const double plus_mass = detail::pair_mass(plus[0], pi0[ordering.first]);
                const double minus_mass = detail::pair_mass(minus[0], pi0[ordering.second]);
                assignments.push_back({
                    std::pow(plus_mass - g6_rho_mass, 2) + std::pow(minus_mass - g6_rho_mass, 2),
                    detail::delta_r(plus[0], pi0[ordering.first]) + detail::delta_r(minus[0], pi0[ordering.second]),
                    ordering.first, ordering.second, plus_mass, minus_mass});
            }
            auto best = assignments[0];
            const auto& other = assignments[1];
            if (other.cost < best.cost) best = other;
            if (std::abs(assignments[0].cost - assignments[1].cost) <= 1.0e-4)
                best = assignments[0].sum_delta_r <= assignments[1].sum_delta_r ? assignments[0] : assignments[1];
            if (detail::rho_pass(best.plus_mass) && detail::rho_pass(best.minus_mass))
                add(detail::rho_side(plus[0], pi0[best.plus_index]),
                    detail::rho_side(minus[0], pi0[best.minus_index]), {best.plus_mass, best.minus_mass});
        }
    } else if (channel == OrderedChannel::pi_a1 && plus.size() == 2U && minus.size() == 2U && neutral.empty()) {
        for (std::size_t direct = 0; direct < plus.size(); ++direct)
            add(detail::pion_side(plus[direct]), detail::a1_side(minus, plus[1U - direct]));
    } else if (channel == OrderedChannel::a1_pi && plus.size() == 2U && minus.size() == 2U && neutral.empty()) {
        for (std::size_t direct = 0; direct < minus.size(); ++direct)
            add(detail::a1_side(plus, minus[1U - direct]), detail::pion_side(minus[direct]));
    } else if (channel == OrderedChannel::rho_a1) {
        const auto pi0 = detail::hardest_neutrals(neutral, 1U);
        if (plus.size() == 2U && minus.size() == 2U && pi0.size() == 1U) {
            for (std::size_t rho = 0; rho < minus.size(); ++rho) {
                const double mass = detail::pair_mass(minus[rho], pi0[0]);
                if (detail::rho_pass(mass))
                    add(detail::a1_side(plus, minus[1U - rho]), detail::rho_side(minus[rho], pi0[0]), {mass});
            }
        }
    } else if (channel == OrderedChannel::a1_rho) {
        const auto pi0 = detail::hardest_neutrals(neutral, 1U);
        if (plus.size() == 2U && minus.size() == 2U && pi0.size() == 1U) {
            for (std::size_t rho = 0; rho < plus.size(); ++rho) {
                const double mass = detail::pair_mass(plus[rho], pi0[0]);
                if (detail::rho_pass(mass))
                    add(detail::rho_side(plus[rho], pi0[0]), detail::a1_side(minus, plus[1U - rho]), {mass});
            }
        }
    } else if (channel == OrderedChannel::a1_a1 && plus.size() == 3U && minus.size() == 3U && neutral.empty()) {
        for (std::size_t plus_opposite = 0; plus_opposite < plus.size(); ++plus_opposite) {
            for (std::size_t minus_opposite = 0; minus_opposite < minus.size(); ++minus_opposite) {
                std::vector<ReconstructedPion> plus_same;
                std::vector<ReconstructedPion> minus_same;
                for (std::size_t index = 0; index < plus.size(); ++index)
                    if (index != plus_opposite) plus_same.push_back(plus[index]);
                for (std::size_t index = 0; index < minus.size(); ++index)
                    if (index != minus_opposite) minus_same.push_back(minus[index]);
                add(detail::a1_side(plus_same, minus[minus_opposite]),
                    detail::a1_side(minus_same, plus[plus_opposite]));
            }
        }
    }
    return result;
}

}  // namespace tauamp::tautau

#endif  // TAUAMP_TAUTAU_ORDERED_CHANNEL_ADAPTER_H_
