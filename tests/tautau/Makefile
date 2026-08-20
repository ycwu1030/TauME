CXX ?= c++
CXXFLAGS := -std=c++17 -Wall -Wextra -Werror -I../../include
ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --libs)

TESTS := kinematics_test.x spin_density_test.x production_test.x pion_pair_test.x pion_pair_hypothesis_test.x pion_pair_marginalization_test.x reference_fixture_test.x lab_frame_pion_test.x legacy_compile_test.x

all: $(TESTS)

HEADERS := $(shell find ../../include/tauamp/tautau -type f -name "*.h")

# Existing legacy declarations trigger unused-variable warnings under -Werror.
legacy_compile_test.x: CXXFLAGS += $(ROOT_CXXFLAGS) -Wno-unused-variable
legacy_compile_test.x: LDFLAGS += $(ROOT_LDFLAGS)

%.x: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TESTS)
