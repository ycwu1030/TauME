#include <cassert>
#include <cstdio>
#include <string>

#include "pion_pair_delphes_fixture.h"
#include "tauamp/delphes/common_channel_delphes_reader.h"

int main() {
    const std::string path = "/private/tmp/tauamp_common_channel_reader_fixture.root";
    std::remove(path.c_str());
    tauamp::delphes::test::write_pion_pair_delphes_fixture(path);

    const tauamp::delphes::CommonChannelDelphesReaderConfig config{
        {0.0, 0.0, -2.13, 2.13}, {0.0, 0.0, 2.13, 2.13}};
    const tauamp::delphes::CommonChannelDelphesReader reader(path, config);
    assert(reader.entry_count() == 1);
    const auto collections = reader.read(0);
    assert(collections.charged.size() == 2);
    assert(collections.neutral.empty());
    assert(collections.charged[0].pid == 211);
    assert(collections.charged[1].pid == -211);
    assert(collections.charged[0].charge() == 1);
    assert(collections.charged[1].charge() == -1);
    assert(collections.charged[0].momentum.energy() > 0.0);
    std::remove(path.c_str());
}
