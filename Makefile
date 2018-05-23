CXX=$(ITENSOR_CXX)
CXXFLAGS=$(ITENSOR_CXXFLAGS) -Wall -pedantic
LDLIBS=$(ITENSOR_LDLIBS)

OBJECTS=linrot/operators.o dmrg.o
BINARIES=dipoles_dmrg dipoles_dmrg_field dipoles_dmrg_nonlinear


.PHONY: all clean

all: $(BINARIES)

clean:
	rm -rf $(OBJECTS) $(BINARIES)


dipoles_dmrg: dipoles_dmrg.cc $(OBJECTS)

dipoles_dmrg_field: dipoles_dmrg_field.cc $(OBJECTS)

dipoles_dmrg_nonlinear: dipoles_dmrg_nonlinear.cc $(OBJECTS)
