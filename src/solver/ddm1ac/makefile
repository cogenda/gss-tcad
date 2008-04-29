ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    ddm_nt1ac.o semiequ1ac.o insulequ1ac.o elecequ1ac.o 
	cp *.o ..

ddm_nt1ac.o   : ddm_nt1ac.cc
semiequ1ac.o  : semiequ1ac.cc
insulequ1ac.o : insulequ1ac.cc
elecequ1ac.o  : elecequ1ac.cc 

.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<

clean:
	rm -f *.o

