include make.defs

All:
	cd src;    make;

install: src/gss$(EXE)
	cp src/gss$(EXE) bin/
	cp src/material/lib*.so lib/
	cp src/utils/tif2cgns/dumptif$(EXE) bin/


clean:
	cd src;      make clean;
	rm -f bin/gss$(EXE)
	rm -f bin/dumptif$(EXE)
	rm -f lib/*.so


distclean:
	rm -f config.h
	rm -f make.defs
	rm -f config.log
	rm -f config.status

