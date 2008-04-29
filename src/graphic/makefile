ALL: GSS
include ${PETSC_DIR}/bmake/common/base
include ../../make.defs



OBJS = xgraph.o wgraph.o grafix3d.o plot.o plot3d.o vtk_plot.o showmesh.o

GSS:$(OBJS)
	cp *.o ..

xgraph.o   : xgraph.cpp ../include/xgraph.h
wgraph.o   : wgraph.cpp ../include/wgraph.h
grafix3d.o : grafix3d.cpp ../include/grafix3d.h
plot.o     : plot.cc
plot3d.o   : plot3d.cc
vtk_plot.o : vtk_plot.cc
showmesh.o : showmesh.cc

.cc.o:
	$(CXX)  $(CXXFLAGS) -c $< $(INCLUDE) $(PETSC_INCLUDE)
.cpp.o:
	$(CXX)  $(CXXFLAGS) -c $< $(INCLUDE)

clean:
	rm -f *.o
