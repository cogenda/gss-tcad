ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libAlGaAs.so
	cp libAlGaAs.so ..
OBJS = AlGaAs_basic.o AlGaAs_band.o AlGaAs_gen.o AlGaAs_thermal.o \
       AlGaAs_mob_Analytic.o  \
       AlGaAs_Optical.o \
       ../adolc_init.o

libAlGaAs.so:$(OBJS)
	$(CXX) -shared -o libAlGaAs.so $(OBJS) -lm

AlGaAs_band.o    : AlGaAs_band.cc
AlGaAs_basic.o   : AlGaAs_basic.cc
AlGaAs_gen.o     : AlGaAs_gen.cc
AlGaAs_thermal.o : AlGaAs_thermal.cc

AlGaAs_mob_Analytic.o : AlGaAs_mob_Analytic.cc

AlGaAs_Optical.o : AlGaAs_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
