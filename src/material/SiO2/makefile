ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libSiO2.so
	cp libSiO2.so ..

OBJS = 	SiO2_basic.o SiO2_thermal.o SiO2_Optical.o
libSiO2.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libSiO2.so $(OBJS) -lm

SiO2_basic.o    : SiO2_basic.cc
SiO2_thermal.o  : SiO2_thermal.cc
SiO2_Optical.o  : SiO2_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
