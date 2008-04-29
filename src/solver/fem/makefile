ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:    em_fem.o
	cp *.o ..
	
em_fem.o: em_fem.cc

.cc.o:
	$(CXX)   $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<
		
clean:
	rm -f *.o

