CXX	= g++
CXXOBJS = main.o io.o km2ap.o initialize.o movingmesh.o delaunay.o voronoi.o tree.o pm.o time_evolve.o solve_hydro.o geometric.o exArith.o utilities/helperFunctions.o utilities/timer.o treepm.o hilbert.o

INCL	= Makefile
CXXINCL = params.h utils.h particle.h mesh.h def.h tree.h exArith.h utilities/helperFunctions.h utilities/timer.h vec.h
CXXFLAGS = -O3 -Wall -mtune=native -march=native -I/usr/local/include 

DEST	= /usr/local/bin
LIBS	= -L/usr/lib -I/usr/local/include -L/usr/local/lib 
LIBS	+= -lm -lfftw -lrfftw  

EXE	= EXE 


EXE: ${CXXOBJS}
	${CXX} ${CXXFLAGS} ${LIBS} ${CXXOBJS} -o ${EXE} 

CXXOBJS: ${INCL} ${CXXINCL} 

clean:	
	rm -f *.o utilities/*.o *.gnu ${CXXOBJS} ${EXE}
