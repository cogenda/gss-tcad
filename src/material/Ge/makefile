ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libGe.so
	cp libGe.so ..

OBJS = Ge_basic.o Ge_band.o Ge_gen.o Ge_thermal.o \
       Ge_mob_Analytic.o \
       Ge_Optical.o \
       ../adolc_init.o

libGe.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libGe.so $(OBJS) -lm

Ge_band.o    : Ge_band.cc
Ge_basic.o   : Ge_basic.cc
Ge_gen.o     : Ge_gen.cc
Ge_thermal.o : Ge_thermal.cc

Ge_mob_Analytic.o : Ge_mob_Analytic.cc

Ge_Optical.o      : Ge_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
