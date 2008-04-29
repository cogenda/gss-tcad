All:all

include ${PETSC_DIR}/bmake/common/base
include ../../../make.defs

all: dumptif$(EXE)

dumptif$(EXE): main.o tocgns.o  tif_parser.o tif_lex.o
	$(LINK) -O -o dumptif$(EXE)   main.o tocgns.o tif_lex.o tif_parser.o ${PETSC_LIB} $(CLIBS)

tif_parser.o : tif_parser.cc
tif_lex.o    : tif_lex.cc
tocgns.o     : tocgns.cc
main.o       : main.cc

.cc.o:
	$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(INCLUDE) -c $<

tif_lex.cc:
	$(LEX) -l tif_lex.l
	mv lex.yy.c tif_lex.cc

tif_parser.cc:
	$(YACC) -d tif_parser.y 
	mv y.tab.h tif_parser.tab.h
	mv y.tab.c tif_parser.cc

clean:
	rm  -f dumptif$(EXE)
	rm  -f *.o tif_lex.cc tif_parser.cc tif_parser.tab.h

install: 
	cp dumptif$(EXE) ../../../bin
