ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libGaAs.so
	cp libGaAs.so ..
OBJS = GaAs_basic.o GaAs_band.o GaAs_gen.o GaAs_thermal.o \
       GaAs_mob_Analytic.o  GaAs_mob_Hypertang.o \
       GaAs_Optical.o \
       ../adolc_init.o

libGaAs.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libGaAs.so $(OBJS) -lm

GaAs_band.o    : GaAs_band.cc
GaAs_basic.o   : GaAs_basic.cc
GaAs_gen.o     : GaAs_gen.cc
GaAs_thermal.o : GaAs_thermal.cc

GaAs_mob_Analytic.o  : GaAs_mob_Analytic.cc
GaAs_mob_Hypertang.o : GaAs_mob_Hypertang.cc

GaAs_Optical.o       : GaAs_Optical.cc 
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
