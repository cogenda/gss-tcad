ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    mix1.o semiequ1mix.o insulequ1mix.o mixsocket.o
	cp *.o ..

mix1.o         : mix1.cc
semiequ1mix.o  : semiequ1mix.cc
insulequ1mix.o : insulequ1mix.cc
mixsocket.o    : mixsocket.cc
.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<

clean:
	rm -f *.o

