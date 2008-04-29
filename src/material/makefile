ALL: all

include ${PETSC_DIR}/bmake/common/base
include ../../make.defs

all:  phy_scale.o material.o adolc_init.o
	cd Al;      make;
	cd Ag;      make;
	cd Au;      make;
	cd Cu;      make;
	cd TiSi2;   make;
	cd PolySi;  make;
	cd Elec;    make;
	cd Si;      make;
	cd GaAs;    make;
	cd Ge;      make;
	cd SiGe;    make;
	cd AlGaAs;  make;
	cd InGaAs;  make;
	cd InP;     make;
	cd InAs;    make;
	cd InSb;    make;
	cd InN;     make;
	cd HgCdTe;  make;
	cd 3C-SiC;  make;
	cd S-SiO2;  make;
	cd Air;     make;
	cd SiO2;    make;
	cd Nitride; make;
	cd Vacuum;  make;
	cd PML;     make;	
	cp phy_scale.o material.o ..


phy_scale.o: phy_scale.cc
material.o : material.cc matdefine.h
adolc_init.o : adolc_init.cc
	$(CXX) $(CXXFLAGS) -fPIC -c $< $(INCLUDE) $(PETSC_INCLUDE)

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDE) $(PETSC_INCLUDE)


clean:
	cd Al;      make clean;
	cd Ag;      make clean;
	cd Au;      make clean;
	cd Cu;      make clean;
	cd TiSi2;   make clean;
	cd PolySi;  make clean;
	cd Elec;    make clean;
	cd Si;      make clean;
	cd GaAs;    make clean;
	cd Ge;      make clean;
	cd SiGe;    make clean;
	cd AlGaAs;  make clean;
	cd InGaAs;  make clean;
	cd InP;     make clean;
	cd InAs;    make clean;
	cd InSb;    make clean;
	cd InN;     make clean;
	cd HgCdTe;  make clean;
	cd 3C-SiC;  make clean;
	cd S-SiO2;  make clean;
	cd Air;     make clean;
	cd SiO2;    make clean;
	cd Nitride; make clean;
	cd Vacuum;  make clean;
	cd PML;     make clean;	
	rm -f *.o  *.so
