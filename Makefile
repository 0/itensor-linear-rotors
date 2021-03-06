CXX=$(ITENSOR_CXX)
CXXFLAGS=$(ITENSOR_CXXFLAGS) -I. -Wall -pedantic
LDLIBS=$(ITENSOR_LDLIBS)

OBJECTS=linrot/linrigrot.o linrot/operators.o dmrg.o

SOURCES=$(wildcard bin/*.cc)
BINARIES=$(SOURCES:.cc=)


.PHONY: all clean

all: $(BINARIES)

clean:
	rm -rf $(OBJECTS) $(BINARIES)


linrot/linrigrot.o: linrot/linrigrot.cc linrot/linrigrot.h

linrot/operators.o: linrot/operators.cc linrot/operators.h

dmrg.o: dmrg.cc dmrg.h linrot/linrigrot.h linrot/operators.h observer.h


bin/%: bin/%.cc $(OBJECTS) argparse.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECTS) $(LDLIBS)
