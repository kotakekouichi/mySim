#include <iostream>
#include <math.h>
#include "particle.h"
#include "mesh.h"
#include "params.h"
#include "pm.h"
#include "def.h"

/*--- function --- */
void periodic_dx(double &dx);
void solve_friedmanneq(double *da, double Coeff, double dt, paramsClass &params);
int initializeParticleForce(particleClass *particleClass, paramsClass &params);
int gravityforce_pm(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);
/*----------------*/

int time_evolve(particleClass *Particle, pmMeshClass *PMMesh, meshClass *Mesh, paramsClass &params)
{
  cout << "time evolve" << endl;
  int rval = 0;
  const int Ndm = params.Ndm, Ngas = params.Ngas;
  const int Npart = Ndm + Ngas;
  const int NPMGRID = params.NPMGRID;
  const int mass = 0, px = 1, py = 2, pz = 3, energy = params.dim == 3 ? 4 : 3;
  double dt = 1.0e+32;
  double *vx_half = new double[Npart+4], *vy_half = new double[Npart+4], *vz_half = new double[Npart+4];
  double da, da_half, Ha = 0.0, Ha_n = 0.0, Ha_n1 = 0.0; 

  if(params.selfgravity)
  {
    //初期k必要 Fは重力加速度
    initializeParticleForce(Particle, params);
    gravityforce_pm(Particle, PMMesh, params); 
    double dx = 1.0 / static_cast<double>(NPMGRID);
    double F = 0.0;
    
    for(int ipart = 4; ipart < Npart + 4; ++ipart)
    {
      double fx = Particle[ipart].getFx();
      double fy = Particle[ipart].getFy();
      double fz = Particle[ipart].getFz();  

      F = F < (fx * fx + fy * fy + fz * fz) ?
	(fx * fx + fy * fy + fz * fz) : F;
    }

    F = sqrt(F);
    dt = min(0.33 / sqrt(F / dx), dt);
  }
  //if(hydro) solve_hydro;

  for(int ipart = 4; ipart < Npart + 4; ++ipart)
  {
    particleClass *part = &Particle[ipart];
    double xp = part->getPosx(), yp = part->getPosy(), zp = part->getPosz();
    if(part->getParticleType() == dark_matter)
    {
      vx_half[ipart] = part->getVelx() + part->getFx() * dt / 2.0; 
      vy_half[ipart] = part->getVely() + part->getFy() * dt / 2.0; 
      vz_half[ipart] = part->getVelz() + part->getFz() * dt / 2.0; 

      xp += vx_half[ipart] * dt;
      yp += vy_half[ipart] * dt;
      zp += vz_half[ipart] * dt;
    } else if(part->getParticleType() == gas)
    {
      vx_half[ipart] = part->getVelx() + part->getFx() * dt / 2.0;
      vy_half[ipart] = part->getVely() + part->getFy() * dt / 2.0;
      vz_half[ipart] = part->getVelz() + part->getFz() * dt / 2.0;

      part->setPx( part->getPx() + 0.5 * part->getMass() * part->getFx() * dt);
      part->setPy( part->getPy() + 0.5 * part->getMass() * part->getFy() * dt);
      part->setPz( part->getPz() + 0.5 * part->getMass() * part->getFz() * dt); 

      part->setE( part->getE() + 0.5 * part->getMass() * 
	  (part->getVelx() * part->getFx() 
	 + part->getVely() * part->getFy() 
	 + part->getVelz() * part->getFz()) * dt );

      xp += vx_half[ipart] * dt;
      yp += vy_half[ipart] * dt;
      zp += vz_half[ipart] * dt;

    }

    periodic_dx(xp);
    periodic_dx(yp);
    periodic_dx(zp);

    part->setPosx(xp);
    part->setPosy(yp);
    part->setPosz(zp);
  }

  if(params.cosmological) 
  {
    solve_friedmanneq(&da, 1.0, dt, params);
    da_half = da * 0.5; 
    Ha = da / dt / a;
    Ha_n = Ha;
    a += da_half;
  }

  cout << a << " " << 1.0 / a - 1.0 << endl;

  if(params.selfgravity)
  {
    initializeParticleForce(Particle, params);
    gravityforce_pm(Particle, PMMesh, params);
  }
  
  if(params.cosmological)
  {
    a += da_half;
    solve_friedmanneq(&da, 1.0, dt, params);
    Ha_n1 = da / dt / a;
  }

  for(int ipart = 4; ipart < Npart + 4; ++ipart)
  {
    particleClass *part = &Particle[ipart]; 
    
    if(part->getParticleType() == dark_matter)
    {
      part->setVelx(vx_half[ipart] +  part->getFx() * dt * 0.5);
      part->setVely(vy_half[ipart] +  part->getFy() * dt * 0.5);
      part->setVelz(vz_half[ipart] +  part->getFz() * dt * 0.5);
    }
    else if(part->getParticleType() == gas)
    {
       meshClass *ptrMgp = part->getMeshGenerationPoint();
       double oldPx = part->getPx(), oldPy = part->getPy(), oldPz = part->getPz();
       
       part->setMass( part->getMass() - ptrMgp->dQi[mass] * dt);
       part->setPx( oldPx + 0.5 * part->getMass() * part->getFx() * dt - ptrMgp->dQi[px] * dt );
       part->setPy( oldPy + 0.5 * part->getMass() * part->getFy() * dt - ptrMgp->dQi[py] * dt );
       part->setPz( oldPz + 0.5 * part->getMass() *part->getFz() * dt - ptrMgp->dQi[pz] * dt ); 

       part->setE( part->getE() + 0.5 * part->getMass() * 
	   (part->getVelx() * part->getFx() 
	    + part->getVely() * part->getFy() 
	    + part->getVelz() * part->getFz()) * dt - ptrMgp->dQi[energy] * dt );

#if 0
       if(params.cosmological)
       {
	 part->setPx((oldPx - 0.5 * Ha_n * oldPx *dt) / (1.0 + Ha_n1 * dt * 0.5));
	 part->setPy((oldPy - 0.5 * Ha_n * oldPy *dt) / (1.0 + Ha_n1 * dt * 0.5));
	 part->setPz((oldPz - 0.5 * Ha_n * oldPz *dt) / (1.0 + Ha_n1 * dt * 0.5));
	 part->setE((part->getE() - 1.0 * Ha * part->getE() * dt) / (1.0 + 1.0 * Ha_n1 * dt));
       }

#endif
     
       // 正しい？
       part->setVelx(part->getPx() / part->getMass());
       part->setVely(part->getPy() / part->getMass());
       part->setVelz(part->getPz() / part->getMass());
    }

  }

  t += dt;

  delete [] vx_half;
  delete [] vy_half;
  delete [] vz_half;

  cout << "end time evolve " << endl;

  return rval;
}

void periodic_dx(double &dx)
{
  if(dx > RIGHTSIDE)
    dx -= BOXSIZE;
  else if(dx < LEFTSIDE)
    dx += BOXSIZE;
}

void solve_friedmanneq(double *da, double Coeff, double dt, paramsClass &params){

  double a_pre, da_half;
  const double OmegaMatter = params.OmegaMatter, OmegaLambda = params.OmegaLambda;

  da_half = sqrt(OmegaMatter / a + OmegaLambda * a * a) * dt / 2.0 * a* a * Coeff;
  a_pre = a + da_half;
  da_half = sqrt(OmegaMatter /a_pre + OmegaLambda * a_pre * a_pre) * dt*a_pre*a_pre * Coeff;
  *da = da_half;

}


int initializeParticleForce(particleClass *Particle, paramsClass &params)
{
  int rval = 0;
  const int Npart = 4 + params.Ndm + params.Ngas;

  for(int ipart = 4; ipart < Npart; ++ipart)
  {
    Particle[ipart].setFx(0.0);
    Particle[ipart].setFy(0.0);
    Particle[ipart].setFz(0.0);
  }

  return 0;
}

