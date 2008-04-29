ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    ddm_nt2e.o semiequ2e.o insulequ2.o elecequ2.o
	cp *.o ..
	
ddm_nt2e.o   : ddm_nt2e.cc 
semiequ2e.o  : semiequ2e.cc jflux2e.h
insulequ2.o  : insulequ2.cc
elecequ2.o   : elecequ2.cc
.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<
		
clean:
	rm -f *.o

