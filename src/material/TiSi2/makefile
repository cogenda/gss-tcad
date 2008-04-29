ALL: all


include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all:  libTiSi2.so
	cp libTiSi2.so ..

OBJS = 	TiSi2_basic.o TiSi2_thermal.o TiSi2_Optical.o
libTiSi2.so:  $(OBJS)  
	$(CXX)  -shared -o libTiSi2.so $(OBJS) -lm

TiSi2_basic.o    : TiSi2_basic.cc
TiSi2_thermal.o  : TiSi2_thermal.cc
TiSi2_Optical.o  : TiSi2_Optical.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c  -fPIC $< $(INCLUDE) $(PETSC_INCLUDE)

clean:
	rm -f *.o *.so
