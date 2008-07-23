ALL: all
	
include ../../../make.defs

LIBSRC =  csa.c\
          minell.c\
          svd.c\
          nan.c

LIBOBJS = csa.o\
          minell.o\
          svd.o\
          nan.o

all: libcsa.a


libcsa.a: $(LIBOBJS)
	$(AR) cru libcsa.a $(LIBOBJS)
	chmod go+r libcsa.a

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean:	
	rm -f *.o libcsa.a 

