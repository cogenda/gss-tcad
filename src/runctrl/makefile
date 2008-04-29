ALL: GSS

include ${PETSC_DIR}/bmake/common/base
include ../../make.defs

GSS:   ctrl.o 
	cp *.o ..

ctrl.o: ctrl.cc \
        ../include/ctrl.h \
        ../include/ddm_nt1e.h \
        ../include/qddm_nt1e.h \
        ../include/ddm_nt2e.h \
        ../include/ebm_nt3e.h \
        ../include/mix1.h\
        ../include/mix2.h\
        ../include/em_fem.h
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDE) $(PETSC_INCLUDE)


clean:
	rm -f *.o

