ALL: GSS

include ${PETSC_DIR}/bmake/common/base
include ../../make.defs

GSS:   grid.o trimesh.o
	cp *.o ..
	
grid.o: grid.cc skeleton.h
trimesh.o: trimesh.cc skeleton.h
.cc.o:
	$(CXX) $(CXXFLAGS) -c $< $(PETSC_INCLUDE) $(INCLUDE)
	
	
clean:
	rm -f *.o


	
	
	
