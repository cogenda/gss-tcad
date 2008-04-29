ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libHgCdTe.so
	cp libHgCdTe.so ..
OBJS = HgCdTe_basic.o HgCdTe_band.o HgCdTe_gen.o HgCdTe_thermal.o \
       HgCdTe_mob_Analytic.o  \
       HgCdTe_Optical.o \
       ../adolc_init.o

libHgCdTe.so:$(OBJS)
	$(CXX) -shared -o libHgCdTe.so $(OBJS) -lm

HgCdTe_band.o    : HgCdTe_band.cc
HgCdTe_basic.o   : HgCdTe_basic.cc
HgCdTe_gen.o     : HgCdTe_gen.cc
HgCdTe_thermal.o : HgCdTe_thermal.cc

HgCdTe_mob_Analytic.o : HgCdTe_mob_Analytic.cc


HgCdTe_Optical.o : HgCdTe_Optical.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
