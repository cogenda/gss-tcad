include ../../make.defs

GSS:   	log.o lex.o parser.o cmdbuf.o
	cp *.o ..

lex.o   : lex.cc ../include/parser.h
parser.o: parser.cc ../include/parser.h
cmdbuf.o: cmdbuf.cc ../include/cmdbuf.h
log.o   : log.cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDE)


lex.cc: lex.l
	$(LEX) -l -i lex.l
	mv lex.yy.c lex.cc

parser.cc ../include/parser.h: parser.y
	$(YACC) -d parser.y
	mv y.tab.h ../include/parser.h
	mv y.tab.c parser.cc

clean:
	rm -f *.o
	rm -f lex.cc parser.cc ../include/parser.h




