ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libSiO2S.so
	cp libSiO2S.so ..

       
OBJS = SiO2S_basic.o SiO2S_band.o SiO2S_gen.o  SiO2S_thermal.o \
       SiO2S_mob_Constant.o\
       SiO2S_Optical.o \
       ../adolc_init.o

libSiO2S.so:$(OBJS)
	$(CXX) $(LDFSHARED) -o libSiO2S.so $(OBJS) -lm

SiO2S_band.o          : SiO2S_band.cc
SiO2S_basic.o         : SiO2S_basic.cc
SiO2S_gen.o           : SiO2S_gen.cc
SiO2S_thermal.o       : SiO2S_thermal.cc
SiO2S_mob_Constant.o  : SiO2S_mob_Constant.cc
SiO2S_Optical.o       : SiO2S_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
