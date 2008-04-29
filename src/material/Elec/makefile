ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libElec.so
	cp libElec.so ..

OBJS = 	Elec_basic.o Elec_thermal.o Elec_Optical.o
libElec.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libElec.so $(OBJS) -lm

Elec_basic.o    : Elec_basic.cc
Elec_thermal.o  : Elec_thermal.cc
Elec_Optical.o  : Elec_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
