CXX = g++


CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -I../../include -I./

LDFLAGS += $(shell root-config --libs) -lEG

all: test_matrixelement.x

LHEF.o: LHEF.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.x: %.cpp LHEF.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
