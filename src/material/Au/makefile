ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libAu.so
	cp libAu.so ..

OBJS = 	Au_basic.o Au_thermal.o Au_Optical.o
libAu.so:  $(OBJS)  
	$(CXX)  $(LDFSHARED) -o libAu.so $(OBJS) -lm

Au_basic.o    : Au_basic.cc
Au_thermal.o  : Au_thermal.cc
Au_Optical.o  : Au_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
