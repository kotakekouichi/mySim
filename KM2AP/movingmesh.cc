#include <iostream>
#include <vector>
#include <map>
#include "params.h"
#include "mesh.h"
#include "particle.h"
#include "pm.h"
#include "tree.h"
#include "def.h"

using namespace std;

/* --- function --- */
void exactinit();
int initializeDelaunay(particleClass *Particle, meshClass *mesh, delaunayClass *Delaunay, paramsClass &paramsClass);
int delaunaygenerator2d(particleClass *Particle, meshClass *mesh, delaunayClass *delaunay, paramsClass &params);
int delaunaygenerator3d(particleClass *Particle, meshClass *mesh, delaunayClass *delaunay, paramsClass &params);
int voronoigenerator2d(particleClass *Particle, meshClass *Mesh, paramsClass &params);
int voronoigenerator3d(particleClass *Particle, meshClass *Mesh, paramsClass &params);

int clear(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, paramsClass &params);
void outputVoronoi(meshClass *Mesh, paramsClass &params);
void outputParticle(string dmfilename, string gasfilename, particleClass *Particle, paramsClass &params);

int time_evolve(particleClass *Particle, pmMeshClass* PMMesh, meshClass *Mesh, paramsClass &params);
int solve_hydro(meshClass *Mesh, paramsClass &params);
/* ---------------- */

int movingmesh(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, pmMeshClass *PMMesh, paramsClass &params)
{
  int rval = 0;

  if(params.dim == 2)
  {
    // delaynay generator
    for(int loop = 0; loop < 100; ++loop)
    {
      rval = initializeDelaunay(Particle, Mesh, delaunay, params);
      clock_t time= clock();
      rval = delaunaygenerator2d(Particle, Mesh, delaunay, params);
      // voronoi generator
      rval = voronoigenerator2d(Particle, Mesh, params);
      time = clock() - time;
      cout << (double) time / CLOCKS_PER_SEC << endl;
      getchar();

      // solve hydro
      solve_hydro(Mesh, params);

      outputVoronoi(Mesh, params);

      rval = clear(Particle, Mesh, delaunay, params);
      outputParticle("output/dm.dat", "output/gas.dat", Particle, params);
    }
  } else if(params.dim == 3)
  {
#if 1
    for(int loop = 0; loop < 10; ++loop) 
    {
      rval = initializeDelaunay(Particle, Mesh, delaunay, params);
      rval = delaunaygenerator3d(Particle, Mesh, delaunay, params);
      rval = voronoigenerator3d(Particle, Mesh, params);
//      outputVoronoi(Mesh, params);
//      outputParticle("output/dm.dat", "output/gas.dat", Particle, params);
      cout << "end voronoi" << endl;
      solve_hydro(Mesh, params);
      rval = clear(Particle, Mesh, delaunay, params);

    } 

    return rval;
#endif 
 
    const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE; 
    tree<particleClass> Root(params.Ndm + params.Ngas, (LEFTSIDE + RIGHTSIDE) / 2.0, (LEFTSIDE + RIGHTSIDE) / 2.0, (LEFTSIDE + RIGHTSIDE) / 2.0, boxsize, 1);
    Root.setdim(3);
    Root.initialize(Particle);

    while(a < 1.0)
    {
      //delaunay generator 

      Root.construct();

      rval = initializeDelaunay(Particle, Mesh, delaunay, params);
      rval = delaunaygenerator3d(Particle, Mesh, delaunay, params);
      rval = voronoigenerator3d(Particle, Mesh, params);

      outputParticle("output/dm.dat", "output/gas.dat", Particle, params);
      time_evolve(Particle, PMMesh, Mesh, params);

      rval = clear(Particle, Mesh, delaunay, params);
      Root.free();

      //voronoi generator
      //rval = delaunaygenerator3d(Particle, Mesh, delaunay, params);
      //time_evolve(Particle, PMMesh, Mesh, params);
      //break;
    }
    outputParticle("output/dm.dat", "output/gas.dat", Particle, params);
  }


  return rval;
}

 int clear(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, paramsClass &params)
{
  cout << params.Ndelaunay << endl;
  for(int i = 0; i < params.Ndelaunay; ++i)
  {
    delaunay[i].nNeighbor = 0;
    
   for(int j = 0; j < 4; ++j)
    if(delaunay[i].faces[j] != NULL)
    {
      //cout << delaunay[i].indd << " " <<  delaunay[i].faces[j]->interDel[0]->indd << " " <<
	      //delaunay[i].faces[j]->interDel[1]->indd << " " << delaunay[i].flag << endl;

      int idel = 0;
      if(delaunay[i].faces[j]->interDel[idel]->indd == delaunay[i].indd)
      {
        idel = 1; 
      }

      delaunayClass *tmp = delaunay[i].faces[j]->interDel[idel];
      int idx = delaunay[i].faces[j]->interidx[idel];
      //cout << " "<<  j << " " << idx << endl;

      //delete tmp->faces[idx];
      //if(tmp->faces[idx] != NULL)
	//delete delaunay[i].faces[j];
      delaunay[i].faces[j] = NULL;
    }
  }
  params.Ndelaunay = 1;

  const int buf = params.dim == 2 ? 3 : 4;
  for(int imesh = 0; imesh < params.Ngas + buf + params.Nghost; ++imesh)
  {
    Mesh[imesh].delaunay.clear();
    for(int in = 0; in < Mesh[imesh].neighbors.size(); ++in)
      Mesh[imesh].neighbors[in].delaunay.clear();
    Mesh[imesh].neighbors.clear();
    Mesh[imesh].next = NULL;
  }

  return 0;
}

