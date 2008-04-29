ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libInN.so
	cp libInN.so ..
OBJS = InN_basic.o InN_band.o InN_Optical.o InN_thermal.o InN_gen.o\
       InN_mob_Analytic.o \
       ../adolc_init.o

libInN.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libInN.so $(OBJS) -lm

InN_band.o    : InN_band.cc
InN_basic.o   : InN_basic.cc
InN_thermal.o : InN_thermal.cc 
InN_Optical.o : InN_Optical.cc
InN_gen.o     : InN_gen.cc

InN_mob_Analytic.o :InN_mob_Analytic.cc
 
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
