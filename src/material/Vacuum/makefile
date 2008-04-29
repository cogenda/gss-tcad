ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libVacuum.so
	cp libVacuum.so ..

OBJS = 	Vacuum_basic.o Vacuum_thermal.o
libVacuum.so:  $(OBJS)  
	$(CXX)  -shared -o libVacuum.so $(OBJS) -lm

Vacuum_basic.o    : Vacuum_basic.cc
Vacuum_thermal.o  : Vacuum_thermal.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
