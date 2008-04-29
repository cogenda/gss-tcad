ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libAl.so
	cp libAl.so ..

OBJS = 	Al_basic.o Al_thermal.o Al_Optical.o
libAl.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libAl.so $(OBJS) -lm

Al_basic.o    : Al_basic.cc
Al_thermal.o  : Al_thermal.cc
Al_Optical.o  : Al_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
