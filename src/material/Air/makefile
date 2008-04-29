ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libAir.so
	cp libAir.so ..

OBJS = 	Air_basic.o Air_thermal.o Air_Optical.o
libAir.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libAir.so $(OBJS) -lm

Air_basic.o    : Air_basic.cc
Air_thermal.o  : Air_thermal.cc
Air_Optical.o  : Air_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
