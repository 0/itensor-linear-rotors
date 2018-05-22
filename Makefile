CXX=$(ITENSOR_CXX)
CXXFLAGS=$(ITENSOR_CXXFLAGS) -Wall -pedantic
LDLIBS=$(ITENSOR_LDLIBS)

OBJECTS=linrot/operators.o
BINARIES=dipoles_dmrg


.PHONY: all clean

all: $(BINARIES)

clean:
	rm -rf $(OBJECTS) $(BINARIES)


dipoles_dmrg: dipoles_dmrg.cc $(OBJECTS)
