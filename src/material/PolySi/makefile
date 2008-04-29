ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libPolySi.so
	cp libPolySi.so ..

OBJS = 	PolySi_basic.o PolySi_thermal.o PolySi_Optical.o
libPolySi.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libPolySi.so $(OBJS) -lm

PolySi_basic.o    : PolySi_basic.cc
PolySi_thermal.o  : PolySi_thermal.cc
PolySi_Optical.o  : PolySi_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
