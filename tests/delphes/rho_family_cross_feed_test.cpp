#include <cassert>

#include "tauamp/delphes/rho_family_cross_feed.h"

int main() {
    using namespace tauamp::delphes;
    using tauamp::tautau::OrderedChannel;
    RhoFamilyCrossFeedCounter counter(OrderedChannel::pi_rho);
    assert(counter.record(7, {{OrderedChannel::pi_rho, 1, 1, false}, {OrderedChannel::rho_pi, 0, 0, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "selected::pi_rho");
    assert(counter.record(8, {{OrderedChannel::pi_rho, 0, 0, false}, {OrderedChannel::rho_pi, 1, 1, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "cross_feed::pi_rho::rho_pi");
    assert(counter.record(9, {{OrderedChannel::pi_rho, 1, 1, false}, {OrderedChannel::rho_pi, 1, 1, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "ambiguous_multi_hypothesis");
    assert(counter.record(10, {{OrderedChannel::pi_rho, 1, 0, false}, {OrderedChannel::rho_pi, 0, 0, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "unsupported_or_no_solution");
    assert(counter.record(11, {{OrderedChannel::pi_rho, 0, 0, false}, {OrderedChannel::rho_pi, 0, 0, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "rejected::no_ownership_candidate");
    assert(counter.record(12, {{OrderedChannel::pi_rho, 0, 0, true}, {OrderedChannel::rho_pi, 0, 0, false}, {OrderedChannel::rho_rho, 0, 0, false}}).category == "technical_failure");
    assert(counter.matrix()[0][0] == 1 && counter.matrix()[0][1] == 1);
    assert(counter.event_categories().size() == 6);
}
