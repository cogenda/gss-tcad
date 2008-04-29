ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    mix2.o semiequ2mix.o insulequ2mix.o 
	cp *.o ..

mix2.o         : mix2.cc
semiequ2mix.o  : semiequ2mix.cc
insulequ2mix.o : insulequ2mix.cc

.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<

clean:
	rm -f *.o

