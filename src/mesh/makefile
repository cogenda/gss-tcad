include ../../make.defs

OBJ =  mesh.o \
       triangle.o \
       geom.o \
       interface.o \
       element.o   \
       zone.o

gss:    $(OBJ)
	cp *.o ..


mesh.o:   mesh.cc
geom.o:   geom.cc
interface.o: interface.cc 
element.o: element.cc
zone.o: zone.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDE)

triangle.o: triangle.c ../include/triangle.h
	$(CC) -O1 -c -DTRILIBRARY  triangle.c $(INCLUDE)

clean:
	rm -f *.o

