ALL: brkpnts.h

include ../../../make.defs


brkpnts.h: brkpnts$(EXE)
	./brkpnts$(EXE) 
	cp brkpnts.h ../../include

brkpnts$(EXE): brkpnts.c	
	$(CC) -o brkpnts$(EXE) brkpnts.c -lm
	

clean: 
	rm -f *.o brkpnts$(EXE) brkpnts.txt brkpnts.h ../../include/brkpnts.h

