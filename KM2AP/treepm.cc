#include <iostream>
#include "tree.h"
#include "pm.h"
#include "particle.h"

int gravityforce_pm(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);

int gravityforce_treepm(particleClass *Particle, pmMeshClass *PMMesh, tree<particleClass> root, paramsClass &params)
{
  int rval = 0;
  const int Npart = params.Ndm + params.Ngas;

  //short range
  root.gravity_shortrange_force(Particle);

  //lang range
  gravityforce_pm(Particle, PMMesh, params);


  return rval;
}
