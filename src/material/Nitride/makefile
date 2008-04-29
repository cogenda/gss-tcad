ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libNitride.so
	cp libNitride.so ..

OBJS = 	Nit_basic.o Nit_thermal.o Nit_Optical.o
libNitride.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libNitride.so $(OBJS) -lm

Nit_basic.o    : Nit_basic.cc
Nit_thermal.o  : Nit_thermal.cc
Nit_Optical.o  : Nit_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
