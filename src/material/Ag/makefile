ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libAg.so
	cp libAg.so ..

OBJS = 	Ag_basic.o Ag_thermal.o Ag_Optical.o
libAg.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libAg.so $(OBJS) -lm

Ag_basic.o    : Ag_basic.cc
Ag_thermal.o  : Ag_thermal.cc
Ag_Optical.o  : Ag_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
