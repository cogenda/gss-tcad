ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    ebm_nt3e.o semiequ3e.o insulequ3.o elecequ3.o
	cp *.o ..
	
ebm_nt3e.o   : ebm_nt3e.cc 
semiequ3e.o  : semiequ3e.cc flux3e.h
insulequ3.o  : insulequ3.cc
elecequ3.o   : elecequ3.cc

.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<
		
clean:
	rm -f *.o

