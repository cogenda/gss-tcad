ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libSi.so
	cp libSi.so ..

       
OBJS = Si_basic.o Si_band.o Si_gen.o Si_gen_Valdinoci.o Si_thermal.o \
       Si_mob_Constant.o Si_mob_Analytic.o Si_mob_Philips.o Si_mob_HP.o Si_mob_Lombardi.o Si_mob_Lucent.o\
       Si_Optical.o \
       ../adolc_init.o

libSi.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libSi.so $(OBJS) -lm

Si_band.o          : Si_band.cc
Si_basic.o         : Si_basic.cc
Si_gen.o           : Si_gen.cc
Si_gen_Valdinoci.o : Si_gen_Valdinoci.cc
Si_thermal.o       : Si_thermal.cc
Si_mob_Constant.o  : Si_mob_Constant.cc
Si_mob_Analytic.o  : Si_mob_Analytic.cc
Si_mob_Philips.o   : Si_mob_Philips.cc
Si_mob_Lucent.o    : Si_mob_Lucent.cc
Si_mob_HP.o        : Si_mob_HP.cc
Si_mob_Lombardi.o  : Si_mob_Lombardi.cc
Si_Optical.o       : Si_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
