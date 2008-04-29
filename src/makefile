All:all

include ${PETSC_DIR}/bmake/common/base
include ../make.defs

OBJ =  main.o  initial.o data_init.o data_io.o probe.o source.o doping.o  bc.o   \
       material.o phy_scale.o \
       cmdbuf.o lex.o parser.o \
       grid.o trimesh.o mesh.o zone.o element.o interface.o triangle.o geom.o refine.o \
       vtk_plot.o plot.o plot3d.o xgraph.o wgraph.o grafix3d.o showmesh.o\
       log.o  ctrl.o  \
       ddm_nt1e.o semiequ1e.o insulequ1.o elecequ1.o \
       ddm_nt1ac.o semiequ1ac.o insulequ1ac.o elecequ1ac.o \
       mix1.o semiequ1mix.o insulequ1mix.o mixsocket.o\
       ddm_nt2e.o  semiequ2e.o insulequ2.o elecequ2.o\
       mix2.o semiequ2mix.o insulequ2mix.o\
       ebm_nt3e.o semiequ3e.o insulequ3.o elecequ3.o\
       qddm_nt1e.o semiequ1q.o insulequ1q.o elecequ1q.o\
       em_fem.o


all:    objs main.o
	$(LINK) -O -o  gss$(EXE) $(OBJ)  ${PETSC_LIB} $(CLIBS)

objs:
	cd utils;    make;
	cd cmd;      make;
	cd grid;     make;
	cd mesh;     make;
	cd solver;   make;
	cd graphic;  make;
	cd runctrl;  make;
	cd material; make;

main.o: main.cc
	$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c main.cc

clean:
	cd cmd;      make clean;
	cd grid;     make clean;
	cd mesh;     make clean;
	cd solver;   make clean;
	cd graphic;  make clean;
	cd runctrl;  make clean;
	cd material; make clean;
	cd utils;    make clean;
	rm -f *.o
	rm -f gss$(EXE)


