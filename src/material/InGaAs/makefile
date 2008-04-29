ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libInGaAs.so
	cp libInGaAs.so ..
OBJS = InGaAs_basic.o InGaAs_band.o InGaAs_gen.o InGaAs_thermal.o \
       InGaAs_mob_Analytic.o \
       InGaAs_Optical.o \
       ../adolc_init.o

libInGaAs.so:$(OBJS)
	$(CXX) -shared -o libInGaAs.so $(OBJS) -lm

InGaAs_band.o    : InGaAs_band.cc
InGaAs_basic.o   : InGaAs_basic.cc
InGaAs_gen.o     : InGaAs_gen.cc
InGaAs_thermal.o : InGaAs_thermal.cc

InGaAs_mob_Analytic.o : InGaAs_mob_Analytic.cc

InGaAs_Optical.o      : InGaAs_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
