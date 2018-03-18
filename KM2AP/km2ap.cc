#include <iostream>
#include <vector>
#include "params.h"
#include "mesh.h"
#include "particle.h"
#include "pm.h"

using namespace std;

double a, t;
/* --- function --- */
int initialize2d(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params);
int initialize3d(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params);
int initialize_pmmesh(pmMeshClass *PMMesh, paramsClass &params);

int movingmesh(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, pmMeshClass *PmMesh, paramsClass &params);
int free(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay);
int freepm(pmMeshClass *PMMesh);
/* --------------- */

int km2ap(paramsClass &params)
{
  int rval = 0;
  const int Ngas = params.Ngas;
  const int Ndm = params.Ndm;
  const int Npart = Ngas + Ndm + 4;

  //
  a = 1.0 / (1.0 + 50.0);
  //
  particleClass *Particle = new particleClass[Npart*3]; 
  meshClass *Mesh = new meshClass[Ngas*2 + 4];
  pmMeshClass *PmMesh = NULL;
  delaunayClass *Delaunay = new delaunayClass[2000000];

  cout << "start simulation " << endl;
  
  if(params.dim == 2)
  {
    rval = initialize2d(Particle, Mesh, Delaunay, params); 
    rval = movingmesh(Particle, Mesh, Delaunay, PmMesh, params);
  } else if(params.dim == 3) 
  {
    rval = initialize3d(Particle, Mesh, Delaunay, params);

    if(params.periodic && params.selfgravity)
    {
      const int NPMGRID = params.NPMGRID;
      PmMesh = new pmMeshClass[NPMGRID*NPMGRID*NPMGRID];	
      rval = initialize_pmmesh(PmMesh, params);
    }

    rval = movingmesh(Particle, Mesh, Delaunay, PmMesh, params);
  }

  rval = free(Particle, Mesh, Delaunay);
  rval = freepm(PmMesh);

  return rval;
}

int free(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay)
{
  int rval = 0;
  
  cout << "free" << endl;
  delete [] Mesh;
  delete [] Particle;
  delete [] Delaunay;

  return rval;
}


