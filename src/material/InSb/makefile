ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libInSb.so
	cp libInSb.so ..
OBJS = InSb_basic.o InSb_band.o InSb_Optical.o InSb_thermal.o InSb_gen.o\
       InSb_mob_Analytic.o \
       ../adolc_init.o

libInSb.so:$(OBJS)
	$(CXX) -shared -o libInSb.so $(OBJS) -lm

InSb_band.o    : InSb_band.cc
InSb_basic.o   : InSb_basic.cc
InSb_thermal.o : InSb_thermal.cc 
InSb_Optical.o : InSb_Optical.cc
InSb_gen.o     : InSb_gen.cc

InSb_mob_Analytic.o :InSb_mob_Analytic.cc
 
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
