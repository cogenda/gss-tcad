ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    qddm_nt1e.o semiequ1q.o insulequ1q.o elecequ1q.o
	cp *.o ..
	
qddm_nt1e.o   : qddm_nt1e.cc 
semiequ1q.o   : semiequ1q.cc jflux1q.h
insulequ1q.o  : insulequ1q.cc
elecequ1q.o   : elecequ1q.cc

.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<
		
clean:
	rm -f *.o

