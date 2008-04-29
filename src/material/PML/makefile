ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libPML.so
	cp libPML.so ..

OBJS = 	PML_basic.o PML_thermal.o
libPML.so:  $(OBJS)  
	$(CXX)  -shared -o libPML.so $(OBJS) -lm

PML_basic.o    : PML_basic.cc
PML_thermal.o  : PML_thermal.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
