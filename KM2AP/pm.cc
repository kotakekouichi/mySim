#include <iostream>
#include <math.h>
#include <vector>
#include "complex.h"
#include "fftw.h"
#include "rfftw.h"
#include "particle.h"
#include "pm.h"
#include "params.h"
#include "def.h"
#include "utils.h"
#define IND(a, b, c) ( (a) + NPMGRID * ( (b) + (c) * NPMGRID) )
#define INDF(a,b,c) ((c) + (NPMGRID+2) * ((b) + (a)*NPMGRID))
#define INDK(a,b,c) ((c) + (NPMGRID/2+1) * ((b) + (a)*NPMGRID))
#define sink(a,b,c) ( sq(sin(M_PI * (double)a /(double)NPMGRID)) + sq(sin(M_PI * (double)b /(double)NPMGRID)) + sq(sin(M_PI * (double)c /(double)NPMGRID)))  

/* --- function --- */ 
int gravityforce_pm(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);
int initialize_pmmesh(pmMeshClass *PMMesh, paramsClass &params);
int compute_densityfield(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);
int compute_phix(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);
int getunit(double &scale_d, double &scale_t, double &H0, paramsClass &params);
int compute_force_from_2PFDA(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params);
/* ---------------*/

int initialize_pmmesh(pmMeshClass *PMMesh, paramsClass &params)
{
  int rval = 0, index = 0;
  const int NPMGRID = params.NPMGRID;
  const int nx = NPMGRID, ny = NPMGRID, nz = NPMGRID ;
  const double dx = BOXSIZE / static_cast<double>(NPMGRID); 

  for(int i = 0; i < nx; ++i)
  {
    for(int j = 0; j < ny; ++j)
    {
      for(int k = 0; k < nz; ++k)
      {
	index = IND(i, j , k); 
	PMMesh[index].setdx(dx);
	PMMesh[index].setxc((i + 0.5) * dx);
	PMMesh[index].setyc((j + 0.5) * dx);
	PMMesh[index].setzc((k + 0.5) * dx);
      }
    }
  } 

  return rval;
}

int compute_densityfield(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params)
{
  cout << "deinsty" << endl;
  int rval = 0, i, j,k;
  const int NPMGRID = params.NPMGRID;
  const int Npart = params.Ndm + params.Ngas;
  const int nx = NPMGRID, ny = NPMGRID, nz = NPMGRID;
  const double dx = BOXSIZE / static_cast<double>(NPMGRID);

  for(k = 0; k < nz; ++k)
  {
    for(j = 0; j < ny; ++j)
    {
      for(i = 0; i < nx; ++i)
      {
	PMMesh[IND(i, j, k)].setrho(0.0);
      }
    }
  }

  for(int ipart = 4; ipart < Npart + 4; ++ipart)
  {
    particleClass *part = &Particle[ipart];
    i = static_cast<int>(part->getPosx() / dx);
    j = static_cast<int>(part->getPosy() / dx);
    k = static_cast<int>(part->getPosz() / dx);

    int is = part->getPosx() < PMMesh[IND(i, j, k)].getxc() ? -1 : 0; 
    int ie = PMMesh[IND(i,j, k)].getxc() < part->getPosx() ? 1 : 0;

    int js = part->getPosy() < PMMesh[IND(i, j, k)].getyc() ? -1 : 0; 
    int je = PMMesh[IND(i,j, k)].getyc() < part->getPosy() ? 1 : 0;

    int ks = part->getPosz() < PMMesh[IND(i, j, k)].getzc() ? -1 : 0; 
    int ke = PMMesh[IND(i,j, k)].getzc() < part->getPosz() ? 1 : 0;

    for(int kk = k + ks; kk <= k + ke; ++kk)
    {
      for(int jj = j + js; jj <= j + je; ++jj)
      {
	for(int ii = i + is; ii <= i + ie; ++ii)
	{
	  int ix = ii, iy = jj, iz = kk;	  
	  double Lx = 0.0, Ly = 0.0, Lz = 0.0;
	  if(ii == NPMGRID || ii == -1)
	  {
	    ix = ii - (double) ii / fabs(ii) * NPMGRID;
	    Lx = (double) ii / fabs(ii);
	  }
	  if(jj == NPMGRID || jj == -1)
	  {
	    iy = jj - (double) jj / fabs(jj) * NPMGRID;
	    Ly = (double) jj / fabs(jj);
	  }
	  if(kk == NPMGRID || kk == -1)
	  {
	    iz = kk - (double) kk / fabs(kk) * NPMGRID;
	    Lz = (double) kk / fabs(kk);
	  }

	  int ind = IND(ix, iy, iz);
	  double Wx = 0.0, Wy = 0.0, Wz = 0.0;
	  pmMeshClass *pmmesh = &PMMesh[ind]; 
	  if(fabs(pmmesh->getxc() + Lx - part->getPosx()) < pmmesh->getdx())
	    Wx = 1.0 - fabs(part->getPosx() - pmmesh->getxc() - Lx) / pmmesh->getdx();
	  if(fabs(pmmesh->getyc() + Ly - part->getPosy()) < pmmesh->getdx())
	    Wy = 1.0 - fabs(part->getPosy() - pmmesh->getyc() - Ly) / pmmesh->getdx();
	  if(fabs(pmmesh->getzc() + Lz - part->getPosz()) < pmmesh->getdx())
	    Wz = 1.0 - fabs(part->getPosz() - pmmesh->getzc() - Lz) / pmmesh->getdx();

	  double W = Wx * Wy * Wz;

	  pmmesh->setrho(pmmesh->getrho() + part->getMass() * W / cube(pmmesh->getdx()));
	}
      }
    }
  }

  cout << "end density" << endl;
  return rval;
}

int getunit(double &scale_d, double &scale_t, double &H0, paramsClass &params)
{
  int rval = 0;
  scale_d = params.OmegaMatter * 1.88880e-29 * (params.Hubble0/100.0) * (params.Hubble0/100.0) / a / a / a;
  scale_t = a * a / (params.Hubble0*1.0e5/3.08e24);
  H0 = params.Hubble0 / 3.08e24 * 1.0e5;

  return rval;
}

int compute_phix(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params)
{
  int rval = 0;
  const int NPMGRID = params.NPMGRID;
  const int N3 = cube(NPMGRID); 
  const double G = 6.67e-8;
  double rho_c = 0.0, scale_d, scale_t, H0;

  fftw_real *density, *phi_x;
  fftw_complex *phi_k;
  rfftwnd_plan p, pinv;

  getunit(scale_d, scale_t, H0, params);


  //G = 6.67e-8 * scale_d * scale_t * scale_t;
  rho_c = (3.0 * H0 * H0) / (8.0 * G * M_PI);
  //rho_c = 1.8888e-29 *(Hubble0/100.0) *(Hubble0/100.0)/scale_d;
  p = rfftw3d_create_plan(NPMGRID, NPMGRID, NPMGRID, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  pinv = rfftw3d_create_plan(NPMGRID, NPMGRID, NPMGRID, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  density = (fftw_real*) fftw_malloc(sizeof(fftw_real) * NPMGRID * NPMGRID *(NPMGRID + 2));

  for(int i = 0; i < NPMGRID; ++i)
    for(int j = 0; j < NPMGRID; ++j)
    {
      for(int k = 0; k < NPMGRID + 2; ++k)
	density[INDF(i,j,k)] = 0.0; 

      // ∆xphibar = 1.5aOmega(rhobar - 1)
      for(int k = 0; k < NPMGRID; ++k)
	density[INDF(i, j, k)] = 1.5 * a * params.OmegaMatter * (PMMesh[IND(i, j, k)].getrho() * 
	    scale_d / rho_c / params.OmegaMatter * cube(a) - 1.0);
    }

  rfftwnd_one_real_to_complex(p, density, NULL);
  phi_k = (fftw_complex *) density;

  for(int i = 0; i < NPMGRID / 2 ; ++i)
  {
    int ixp = NPMGRID / 2 + i;
    int ixm = NPMGRID / 2 - i;
    for(int j = 0; j < NPMGRID / 2; ++j)
    {
      int iyp = NPMGRID / 2 + j; 
      int iym = NPMGRID / 2 - j;
      for(int k = 0; k < NPMGRID / 2 + 1; ++k)
      {

	double coeff = -1.0 / 4.0 * sq(1.0 / (double) NPMGRID);
	phi_k[INDK(i,j,k)].re *= coeff / sink(i,j,k); 
	phi_k[INDK(ixp,j,k)].re *= coeff / sink(ixm,j,k); 
	phi_k[INDK(i,iyp,k)].re *= coeff / sink(i,iym,k); 
	phi_k[INDK(ixp,iyp,k)].re *= coeff / sink(ixm,iym,k);

	phi_k[INDK(i,j,k)].im *= coeff / sink(i,j,k); 
	phi_k[INDK(ixp,j,k)].im *= coeff / sink(ixm,j,k); 
	phi_k[INDK(i,iyp,k)].im *=  coeff / sink(i,iym,k);
	phi_k[INDK(ixp,iyp,k)].im *= coeff / sink(ixm,iym,k);

	if (sink(i,j,k) == 0.0){
	  phi_k[INDK(i,j,k)].re = 0.0;
	  phi_k[INDK(i,j,k)].im = 0.0;
	}
      }
    }
  }

  rfftwnd_one_complex_to_real(pinv, phi_k, NULL);

  phi_x = (fftw_real *)phi_k;

  for(int i = 0; i < NPMGRID; ++i)
    for(int j = 0; j < NPMGRID;++j)
      for(int k = 0; k < NPMGRID; ++k)
	PMMesh[IND(i,j,k)].setphix(phi_x[INDF(i,j,k)] / N3); 

  fftwnd_destroy_plan(p);
  fftwnd_destroy_plan(pinv);

  return rval;
}

int compute_force_from_2PFDA(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params)
{
  int rval = 0;
  const int Npart = params.Ndm + params.Ngas + 4;
  const int NPMGRID = params.NPMGRID;
  const double dx = BOXSIZE / static_cast<double>(NPMGRID);
/*
  for(int ipart = 4; ipart < Npart; ++ipart)
  {
    Particle[ipart].setFx(0.0);
    Particle[ipart].setFy(0.0);
    Particle[ipart].setFz(0.0);
  }
*/
  for(int k = 0; k < NPMGRID; ++k)
  {
    for(int j = 0; j < NPMGRID; ++j)
    {
      for(int i = 0; i < NPMGRID; ++i)
      {
	int i1 = i + 1 == NPMGRID ? i1 = 0 : i + 1; 
	int i0 = i - 1 == -1 ? i0 = NPMGRID - 1 : i - 1; 

	int j1 = j + 1 == NPMGRID ? j1 = 0 : j + 1; 
	int j0 = j - 1 == -1 ? j0 = NPMGRID - 1 : j - 1; 

	int k1 = k + 1 == NPMGRID ? k1 = 0 : k + 1; 
	int k0 = k - 1 == -1 ? k0 = NPMGRID - 1 : k - 1; 

	PMMesh[IND(i,j,k)].setax((PMMesh[IND(i1,j,k)].getphix() - PMMesh[IND(i0,j,k)].getphix()) / (-2.0 * dx));
	PMMesh[IND(i,j,k)].setay((PMMesh[IND(i,j1,k)].getphix() - PMMesh[IND(i,j0,k)].getphix()) / (-2.0 * dx));
	PMMesh[IND(i,j,k)].setaz((PMMesh[IND(i,j,k1)].getphix() - PMMesh[IND(i,j,k0)].getphix()) / (-2.0 * dx));
      }
    } 
  }

  for(int ipart = 4;ipart < Npart; ++ipart)
  {
    particleClass *part = &Particle[ipart];
    int i = static_cast<int>(part->getPosx() / dx);
    int j = static_cast<int>(part->getPosy() / dx);
    int k = static_cast<int>(part->getPosz() / dx);

    int is = part->getPosx() < PMMesh[IND(i, j, k)].getxc() ? -1 : 0; 
    int ie = PMMesh[IND(i,j, k)].getxc() < part->getPosx() ? 1 : 0;

    int js = part->getPosy() < PMMesh[IND(i, j, k)].getyc() ? -1 : 0; 
    int je = PMMesh[IND(i,j, k)].getyc() < part->getPosy() ? 1 : 0;

    int ks = part->getPosz() < PMMesh[IND(i, j, k)].getzc() ? -1 : 0; 
    int ke = PMMesh[IND(i,j, k)].getzc() < part->getPosz() ? 1 : 0;

    for(int kk = k + ks; kk <= k + ke; ++kk)
    {
      for(int jj = j + js; jj <= j + je; ++jj)
      {
	for(int ii = i + is; ii <= i + ie; ++ii)
	{

	  int ix = ii, iy = jj, iz = kk;
	  double Lx = 0.0,Ly =0.0,Lz = 0.0;

	  if(ii == NPMGRID || ii == -1){
	    ix = ii - (double)ii / fabs(ii) * NPMGRID; 
	    Lx = (double)ii / fabs(ii);
	  }
	  if(jj == NPMGRID || jj == -1){
	    iy = jj - (double)jj / fabs(jj) * NPMGRID;
	    Ly = (double)jj / fabs(jj);
	  }
	  if(kk == NPMGRID || kk == -1){
	    iz = kk - (double)kk / fabs(kk) * NPMGRID;
	    Lz = (double)kk / fabs(kk);
	  }
	  int ind = IND(ix,iy,iz);
	  double Wx = 0.0, Wy = 0.0, Wz = 0.0, W = 0.0;
	  pmMeshClass *pmmesh = &PMMesh[ind];
	  if(fabs(pmmesh->getxc() + Lx - part->getPosx()) < pmmesh->getdx())
	    Wx = 1.0 - fabs(part->getPosx() - pmmesh->getxc() - Lx) / pmmesh->getdx();
	  if(fabs(pmmesh->getyc() + Ly - part->getPosy()) < pmmesh->getdx())
	    Wy = 1.0 - fabs(part->getPosy() - pmmesh->getyc() - Ly) / pmmesh->getdx();
	  if(fabs(pmmesh->getzc() + Lz - part->getPosz()) < pmmesh->getdx())
	    Wz = 1.0 - fabs(part->getPosz() - pmmesh->getzc() - Lz) / pmmesh->getdx();

	  W = Wx * Wy * Wz;
	  part->setFx(part->getFx() + pmmesh->getax() * W);
	  part->setFy(part->getFy() + pmmesh->getay() * W);
	  part->setFz(part->getFz() + pmmesh->getaz() * W);
	}
      }
    }
  }


  return rval;
}


int gravityforce_pm(particleClass *Particle, pmMeshClass *PMMesh, paramsClass &params)
{
  int rval = 0;

  rval = compute_densityfield(Particle, PMMesh, params);
  rval = compute_phix(Particle, PMMesh, params);
  rval = compute_force_from_2PFDA(Particle, PMMesh, params);

  return rval;

}

int freepm(pmMeshClass *PMMesh)
{
  int rval = 0;

  delete [] PMMesh;

  return rval;
}
