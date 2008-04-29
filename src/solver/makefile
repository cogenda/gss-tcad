ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../make.defs

all:    initial.o doping.o source.o data_io.o probe.o data_init.o refine.o bc.o
	cd ddm1e;       make;
	cd ddm1ac;      make;
	cd mix1;        make;
	cd ddm2e;       make;
	cd mix2;        make;
	cd ebm3e;       make;
	cd qddm1e;      make;
	cd fem;         make;
	cp *.o ..

initial.o   : initial.cc
source.o    : source.cc
doping.o    : doping.cc
data_init.o : data_init.cc
data_io.o   : data_io.cc
probe.o     : probe.cc
refine.o    : refine.cc ../include/bsolver.h ../include/mesh.h
bc.o        : bc.cc

.cc.o:
	$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<

clean:
	cd ddm1e;       make clean;
	cd ddm1ac;      make clean;
	cd mix1;        make clean;
	cd ddm2e;	make clean;
	cd mix2;        make clean;
	cd ebm3e;       make clean;
	cd qddm1e;      make clean;
	cd fem;         make clean;
	rm -f *.o

