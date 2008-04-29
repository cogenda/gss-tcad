ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    ddm_nt1e.o semiequ1e.o insulequ1.o elecequ1.o
	cp *.o ..
	
ddm_nt1e.o  : ddm_nt1e.cc 
semiequ1e.o : semiequ1e.cc jflux1e.h
insulequ1.o : insulequ1.cc
elecequ1.o  : elecequ1.cc
.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<
		
clean:
	rm -f *.o

