ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libSiGe.so
	cp libSiGe.so ..
OBJS = SiGe_basic.o SiGe_band.o SiGe_gen.o SiGe_thermal.o \
       SiGe_mob_Constant.o SiGe_mob_Analytic.o SiGe_mob_Philips.o SiGe_mob_HP.o SiGe_mob_Lombardi.o SiGe_mob_Lucent.o\
       SiGe_Optical.o  \
       ../adolc_init.o

libSiGe.so:$(OBJS)
	$(CXX) -shared -o libSiGe.so $(OBJS) -lm

SiGe_band.o    : SiGe_band.cc
SiGe_basic.o   : SiGe_basic.cc
SiGe_gen.o     : SiGe_gen.cc
SiGe_thermal.o : SiGe_thermal.cc
SiGe_mob_Constant.o : SiGe_mob_Constant.cc
SiGe_mob_Analytic.o : SiGe_mob_Analytic.cc
SiGe_mob_Philips.o  : SiGe_mob_Philips.cc
SiGe_mob_Lombardi.o : SiGe_mob_Lombardi.cc
SiGe_mob_Lucent.o   : SiGe_mob_Lucent.cc
SiGe_mob_HP.o       : SiGe_mob_HP.cc

SiGe_Optical.o      : SiGe_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
