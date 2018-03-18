CXX	= g++
CXXOBJS = main.o simplex.o 

INCL	= Makefile
CXXINCL = opt.h simplex.h position.h	params.h utils.h function.h
CXXFLAGS = -O3 -Wall 

DEST	= /usr/local/bin
LIBS	= -L/usr/lib 
LIBS	+= -lm 

EXE	= Optimize 


EXE: ${CXXOBJS}
	${CXX} ${CXXFLAGS} ${CXXOBJS} ${LIBS} -o ${EXE} 

CXXOBJS: ${INCL} ${CXXINCL}

clean:	
	rm -f *.o *.gnu ${CXXOBJS} ${EXE}
