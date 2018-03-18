CXX	= g++
CXXOBJS = main.o simplex.o 

CXXINCL = simplex.h Makefile
CXXFLAGS = -O3 -pg 

DEST	= /usr/local/bin
LIBS	= -L/usr/lib 
LIBS	+= -lm 

EXE	= simplex 


EXE: ${CXXOBJS}
	${CXX} ${CXXFLAGS} ${CXXOBJS} ${LIBS} -o ${EXE} 

CXXOBJS: ${CXXINCL}

clean:	
	rm -f *.o *~ ${CXXOBJS} ${EXE}
