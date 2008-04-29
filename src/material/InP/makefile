ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libInP.so
	cp libInP.so ..
OBJS = InP_basic.o InP_band.o InP_thermal.o InP_gen.o InP_Optical.o \
       InP_mob_Analytic.o \
       ../adolc_init.o

libInP.so:$(OBJS)
	$(CXX) -shared -o libInP.so $(OBJS) -lm

InP_band.o    : InP_band.cc
InP_basic.o   : InP_basic.cc
InP_thermal.o : InP_thermal.cc
InP_gen.o     : InP_gen.cc
InP_Optical.o : InP_Optical.cc

InP_mob_Analytic.o : InP_mob_Analytic.cc
 
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
