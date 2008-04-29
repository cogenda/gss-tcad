ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libInAs.so
	cp libInAs.so ..
OBJS = InAs_basic.o InAs_band.o InAs_Optical.o InAs_thermal.o InAs_gen.o\
       InAs_mob_Analytic.o \
       ../adolc_init.o

libInAs.so:$(OBJS)
	$(CXX) -shared -o libInAs.so $(OBJS) -lm

InAs_band.o    : InAs_band.cc
InAs_basic.o   : InAs_basic.cc
InAs_thermal.o : InAs_thermal.cc 
InAs_Optical.o : InAs_Optical.cc
InAs_gen.o     : InAs_gen.cc

InAs_mob_Analytic.o :InAs_mob_Analytic.cc
 
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
