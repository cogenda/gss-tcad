ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libCu.so
	cp libCu.so ..

OBJS = 	Cu_basic.o Cu_thermal.o Cu_Optical.o
libCu.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libCu.so $(OBJS) -lm

Cu_basic.o    : Cu_basic.cc
Cu_thermal.o  : Cu_thermal.cc
Cu_Optical.o  : Cu_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
