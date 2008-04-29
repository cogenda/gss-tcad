ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libSiC3C.so
	cp libSiC3C.so ..

OBJS = 3C-SiC_basic.o 3C-SiC_band.o 3C-SiC_gen.o  3C-SiC_thermal.o \
       3C-SiC_mob_Analytic.o   3C-SiC_Optical.o \
       ../adolc_init.o

libSiC3C.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libSiC3C.so $(OBJS) -lm

3C-SiC_band.o          : 3C-SiC_band.cc
3C-SiC_basic.o         : 3C-SiC_basic.cc
3C-SiC_gen.o           : 3C-SiC_gen.cc
3C-SiC_thermal.o       : 3C-SiC_thermal.cc
3C-SiC_mob_Analytic.o  : 3C-SiC_mob_Analytic.cc
3C-SiC_Optical.o       : 3C-SiC_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
