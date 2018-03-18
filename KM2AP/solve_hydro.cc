#include <iostream>
#include <math.h>
#include "mesh.h"
#include "utils.h"
#include "params.h"

using namespace std;

struct variables{
  double rho;
  double u;
  double p;
};

struct Input_vec_Riemann{
  struct variables L;
  struct variables R;
};

struct Riemann_outputs{
  double P_M;
  double S_M;
  variables var; 
};


void compute_w(double *OoldL, double *QoldR, double *f, double *w, double *n, paramsClass &params, meshClass *MeshL, meshClass *MeshR);
void compute_wp(double *OoldL, double *QoldR, meshClass *MeshL, meshClass *MeshR, double *WpL, double *WpR, double *n, double *w, paramsClass &parms);

double compute_timestep(meshClass* Mesh, paramsClass &params);

int compute_gradient(int imesh, double dfdr[5][3], double *w, meshClass *Mesh, int nhydro, int ndim);
int compute_gradient_t(double *Wp, double dfdr[5][3], double *dfdt, int ndim, int nhydro);

void riemann_solver(double *WF, double *uleft, double *uright, paramsClass &params);
double guess_for_pressure(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, 
    double v_line_L, double v_line_R, double cs_L, double cs_R);
int iterative_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_standard(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R);

void sample_reimann_vaccum_internal(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_vaccum_left(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_vaccum_right(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R);

int solve_hydro(meshClass *Mesh, paramsClass &params)
{
  int rval = 0;
  double f[3], w[5], norm, gamma = 5.0 / 3.0, Aij, dt;
  double WpL[5], WppL[5], WpppL[5], WpR[5], WppR[5], WpppR[5];
  double dfdrL[5][3], dfdtL[5], dfdrR[5][3], dfdtR[5];
  double Flux[5], Wlab[5], WF[5];
  meshClass *MeshL, *MeshR; 

  dt = compute_timestep(Mesh, params);

  cout << "sovle hydro dt = " << dt << endl;

  const int buf = params.dim == 3 ? 4 : 3;
  const int Ngas = params.Ngas + buf;
  const int nhydro = params.dim == 3 ? 5 : 4;
  const int mass = 0, px = 1, py = 2, pz = 3, energy = params.dim == 3 ? 4 : 3;
  const int ndim = params.dim;

  double **Qold = new double*[Ngas], **Qnew = new double*[Ngas];

  //allocate
  for(int imesh = 0; imesh < Ngas; ++imesh)
  {
    Qold[imesh] = new double[nhydro];
    Qnew[imesh] = new double[nhydro];
  }

  //initialize 
  for(int imesh = buf; imesh < Ngas; ++imesh)
  {
    meshClass *mesh = &Mesh[imesh];
    Qold[imesh][mass] = Qnew[imesh][mass] = mesh->getMass();
    Qold[imesh][px] = Qnew[imesh][px] = mesh->getPx();
    Qold[imesh][py] = Qnew[imesh][py] = mesh->getPy();
    Qold[imesh][pz] = Qnew[imesh][pz] = mesh->getPz();
    Qold[imesh][energy] = Qnew[imesh][energy] = mesh->getE();
    //Qold[imesh][4] = Qnew[imesh][4] = 0.0;
  }

  for(int lmesh = buf; lmesh < Ngas; ++lmesh)
  {
    MeshL = &Mesh[lmesh]; 

    double QoldL[nhydro], dQi[nhydro];

    for(int i = 0; i < nhydro; ++i)
    {
      QoldL[i] = Qold[lmesh][i]; 
      dQi[i] = 0.0;
      MeshL->dQi[i] = 0.0;
    }

    if(!(0.1 < MeshL->getPosx() && MeshL->getPosx() < .9 && 0.1 < MeshL->getPosy() && MeshL->getPosy() < .9 && 0.1 < MeshL->getPosz() && MeshL->getPosz() < 0.9))
      continue;
    
    for(int inbr = 0; inbr < MeshL->npoly; ++inbr)
    {
      int rmesh = MeshL->neighbors[inbr].index; 
      if(rmesh < buf) continue;

      MeshR = &Mesh[rmesh];

      if(!(0.1 < MeshR->getPosx() && MeshR->getPosx() < .9 && 0.1 < MeshR->getPosy() && MeshR->getPosy() < .9 && 0.1 < MeshR->getPosz() && MeshR->getPosz() < 0.9))
      continue;
      
      double sL[3] = {MeshL->getxg(), MeshL->getyg(), MeshL->getzg()}, sR[3] = {MeshR->getxg(), MeshR->getyg(), MeshR->getzg()};
      double QoldR[nhydro];

      for(int i = 0; i < nhydro; ++i) QoldR[i] = Qold[rmesh][i];

      for(int i = 0; i < ndim; ++i) f[i] = MeshL->neighbors[inbr].f[i];

      if(!(0.1 < f[0] < .9 && 0.1 < f[1]< .9 && 0.1 < f[2] < 0.9))
      {
	//cout << MeshR->index << endl;
	//cout << f[0] << " " << f[1] << "  " << f[2] << endl;
	//getchar();
      } 
      double n[ndim];
      compute_w(QoldL, QoldR, f, w, n, params, MeshL, MeshR);
      //mesh_rounder();
      
      compute_wp(QoldL, QoldR, MeshL, MeshR, WpL, WpR, n, w, params);

      compute_gradient(lmesh, dfdrL, w, Mesh, nhydro, ndim);
      compute_gradient(rmesh, dfdrR, w, Mesh, nhydro, ndim);

      ///////////
     compute_gradient_t(WpL, dfdrL, dfdtL, ndim, nhydro);
     compute_gradient_t(WpR, dfdrR, dfdtR, ndim, nhydro);

      /*dfdtL[0] = -(WpL[1]*dfdrL[0][0] + WpL[2]*dfdrL[0][1] + WpL[0]*(dfdrL[1][0]+dfdrL[2][1]));
      dfdtL[1] = -(WpL[1]*dfdrL[1][0] + WpL[2]*dfdrL[1][1] + dfdrL[3][0] / WpL[0]);
      dfdtL[2] = -(WpL[1]*dfdrL[2][0] + WpL[2]*dfdrL[2][1] + dfdrL[3][1] / WpL[0]);
      dfdtL[3] = -(gamma*WpL[3]*(dfdrL[1][0] + dfdrL[2][1]) + WpL[1] * dfdrL[3][0] + WpL[2] * dfdrL[3][1]);

      dfdtR[0] = -(WpR[1]*dfdrR[0][0] + WpR[2]*dfdrR[0][1] + WpR[0]*(dfdrR[1][0]+dfdrR[2][1]));
      dfdtR[1] = -(WpR[1]*dfdrR[1][0] + WpR[2]*dfdrR[1][1] + dfdrR[3][0] / WpR[0]);
      dfdtR[2] = -(WpR[1]*dfdrR[2][0] + WpR[2]*dfdrR[2][1] + dfdrR[3][1] / WpR[0]);
      dfdtR[3] = -(gamma*WpR[3]*(dfdrR[1][0] + dfdrR[2][1]) + WpR[1] * dfdrR[3][0] + WpR[2] * dfdrR[3][1]);
       */
      //MeshL.getForce()...
      ///////////

      for(int i = 0; i < nhydro; ++i)
      {
	WppL[i] = WpL[i] + dfdtL[i] * dt / 2.0 * 0.0;
	WppR[i] = WpR[i] + dfdtR[i] * dt / 2.0 * 0.0;

	for(int idim = 0; idim < ndim; ++idim)
	{
           //WppL[i] += dfdrL[i][idim] * (f[idim] - sL[idim]);	
	   //WppR[i] += dfdrR[i][idim] * (f[idim] - sR[idim]);
	}
      }

      double cos = n[0], sin = n[1], cos1 = 0.0, cos2 = 0.0, sin1 = 0.0, sin2 = 0.0, n2[3] ={n[0],n[1],n[2]};
      double Wpp2L[5], Wpp2R[5];
      double nvec[3] = {n[0], n[1], n[2]};

      if(ndim == 2)
      {
	WpppL[0] = WppL[0];
	WpppR[0] = WppR[0];
	WpppL[1] = cos * WppL[1] + sin * WppL[2];
	WpppR[1] = cos * WppR[1] + sin * WppR[2]; 
	WpppL[2] = -sin * WppL[1] + cos * WppL[2];
	WpppR[2] = -sin * WppR[1] + cos * WppR[2]; 
	WpppL[3] = WppL[3];
	WpppR[3] = WppR[3];
      } else if(ndim == 3)
      {

	WpppL[0] = WppL[0];
	WpppR[0] = WppR[0];
	// rotation to xz
	norm = sqrt(n[0]*n[0] +n[1]*n[1]);
	if(norm != 0.0){
	  cos1 = n[0]/norm; sin1 = n[1]/norm;	
	  Wpp2L[1] = cos1 * WppL[1] + sin1 * WppL[2];
	  Wpp2L[2] = -sin1 * WppL[1] + cos1 * WppL[2];
	  Wpp2L[3] = WppL[3];
	  Wpp2R[1] = cos1 * WppR[1] + sin1 * WppR[2];
	  Wpp2R[2] = -sin1 * WppR[1] + cos1 * WppR[2];
	  Wpp2R[3] = WppR[3];
	  n2[0]  = cos1 * n[0] + sin1 * n[1];
	  n2[1] = -sin1 * n[0] + cos1 * n[1];
	  n2[2] = n[2];
	  //for(int i=0;i<3;i++)
	    //n[i] = n2[i];
	}
	else{
	  for(int i=1;i<4;i++){
	    Wpp2L[i] = WppL[i];
	    Wpp2R[i] = WppR[i];
	  }
	}
	// rotation to x
	norm = sqrt(n2[0]*n2[0] + n2[2]*n2[2]);
	if(norm != 0.0){
	  cos2 = n2[0]/norm; sin2 = n2[2]/norm;
	  WpppL[1] = cos2 * Wpp2L[1] + sin2 * Wpp2L[3];
	  WpppL[2] = Wpp2L[2];
	  WpppL[3] = -sin2 * Wpp2L[1] + cos2 * Wpp2L[3];
	  WpppR[1] = cos2 * Wpp2R[1] + sin2 * Wpp2R[3];
	  WpppR[2] = Wpp2R[2];
	  WpppR[3] = -sin2 * Wpp2R[1] + cos2 * Wpp2R[3];
	}
	else{
	  for(int i=1;i<4;i++){
	    WpppL[i] = Wpp2L[i];
	    WpppR[i] = Wpp2R[i];
	  }
	}

	WpppL[4] = WppL[4];
	WpppR[4] = WppR[4];

      }

      //riemann 
      riemann_solver(WF, WpppL, WpppR, params);

      // reverse //
      if(ndim == 2)
      {
	double ff[2] = {WF[1], WF[2]};

	WF[1] = cos * ff[0] - sin * ff[1];
	WF[2] = sin * ff[0] + cos * ff[1];
      } else if(ndim == 3)
      {
	double ff[5];
	for(int i=0;i<5;i++)
	  ff[i] = WF[i];

	if(cos2 <= 1.0 && sin2 <= 1.0 ){
	  WF[1] = cos2 * ff[1] - sin2 * ff[3];
	  WF[2] = ff[2];
	  WF[3] = sin2 * ff[1] + cos2 * ff[3];
	}
	else{
	  for(int i=0;i<5;i++)
	    WF[i] = ff[i];
	}
	for(int i=0;i<5;i++)
	  ff[i] = WF[i];
	if(cos1 <= 1.0 && sin1 <= 1.0){
	  WF[1] = cos1 * ff[1] - sin1 * ff[2];
	  WF[2] = sin1 * ff[1] + cos1 * ff[2];
	  WF[3] = ff[3];
	}
	else{
	  for(int i=0;i<5;i++)
	    WF[i] = ff[i];
	}
      }

      for(int i = 0; i < nhydro; ++i) Wlab[i] = WF[i] + w[i];

      Aij = MeshL->neighbors[inbr].Aij;

      double elab = 0.0, WF_n = 0.0;

      for(int i = px, in = 0; i < energy; ++i, ++in)
      {
	elab += Wlab[mass] * sq(Wlab[i]) / 2.0; 
	WF_n += WF[i] * n[in];
      } 
      elab +=  Wlab[energy] / (gamma - 1.0);

      Flux[mass] = Wlab[mass] * WF_n;

      for(int i = px, in = 0; i < energy; ++i, ++in)
	Flux[i] = Wlab[mass] * Wlab[i] * WF_n + Wlab[energy] * n[in];

      Flux[energy] = 0.0;
      for(int i = px, in = 0; i < energy; ++i, ++in)
	Flux[energy] += (elab * WF[i] + Wlab[energy] * Wlab[i]) * n[in];

      for(int i = 0; i < nhydro; ++i)
      {
	dQi[i] += - Flux[i] * Aij;
	MeshL->dQi[i] += -Flux[i] * Aij;
      }

    }

    for(int i = 0; i < nhydro; ++i) Qnew[lmesh][i] += dt * dQi[i];

    //cout << MeshL->getPosx() << " " << MeshL->getPosy() << " " <<  MeshL->getPosz() << " " << MeshL->V << endl;

  }

  for(int lmesh = buf; lmesh < Ngas; ++lmesh)
  {
    meshClass *mesh = &Mesh[lmesh];
    mesh->setMass(Qnew[lmesh][mass]);
    mesh->setPx(Qnew[lmesh][px]);
    mesh->setPy(Qnew[lmesh][py]);
    mesh->setPz(Qnew[lmesh][pz]);
    mesh->setE(Qnew[lmesh][energy]);

    mesh->setPosx( mesh->getPosx() + dt * (Qold[lmesh][px] / Qold[lmesh][mass]) / 2.0);
    mesh->setPosy( mesh->getPosy() + dt * (Qold[lmesh][py] / Qold[lmesh][mass]) / 2.0);
    mesh->setPosz( mesh->getPosz() + dt * (Qold[lmesh][pz] / Qold[lmesh][mass]) / 2.0);

    mesh->setPosx( mesh->getPosx() + dt * (mesh->getPx() / mesh->getMass()) / 2.0);
    mesh->setPosy( mesh->getPosy() + dt * (mesh->getPy() / mesh->getMass()) / 2.0);
    mesh->setPosz( mesh->getPosz() + dt * (mesh->getPz() / mesh->getMass()) / 2.0);
    if(!(0.0 < mesh->getPosx() && mesh->getPosx() < 1.0))
    { cout << "warning x " <<  mesh->getPosx() << " " << mesh->V << endl; getchar();}

    if(!(0.0 < mesh->getPosy() && mesh->getPosy() < 1.0))
    { cout << "warning y " << mesh->getPosy() << " " << mesh->V << endl; getchar();}

    if(!(0.0 < mesh->getPosz() && mesh->getPosz() < 1.0))
    { cout << "warning z" << mesh->getPosz() << " " << mesh->V << endl; getchar();}
  }

  t += dt;

  //free
  for(int imesh = 0; imesh < Ngas; ++imesh)
  {
    delete Qold[imesh];
    delete Qnew[imesh];
  }

  delete [] Qold;
  delete [] Qnew;

  return  rval;
}

double compute_timestep(meshClass *Mesh, paramsClass &params)
{
  double dbl = 1.0e+32, Ri, ci, vi, P, rho;
  const double gamma = 5.0 / 3.0, Ccfl = 0.3;
  const int buf = params.dim == 3 ? 4 : 3;
  const int Ngas  = params.Ngas + buf;

  switch(params.dim)
  {
    case(2):
      for(int i = buf; i < Ngas; ++i)
      {

	Ri = sqrt(Mesh[i].S / M_PI); 
	rho = Mesh[i].getMass() / Mesh[i].S;
	P = (gamma - 1.0) * (Mesh[i].getE() - 1.0 / 2.0 / Mesh[i].getMass() * (Mesh[i].getPx() * Mesh[i].getPx() + Mesh[i].getPy() * Mesh[i].getPy()))/ Mesh[i].S;
	ci = sqrt(gamma * P / rho);
	vi = sqrt(Mesh[i].getPx() * Mesh[i].getPx()  + Mesh[i].getPy() * Mesh[i].getPy()) / Mesh[i].getMass();
	dbl = min(dbl, Ccfl* Ri / (ci + vi));
      } 
      break;

    case(3):
      for(int i = buf; i < Ngas; ++i)
      {
	Ri = pow(Mesh[i].V * 3.0 / 4.0, 1.0 / 3.0);
	rho = Mesh[i].getMass() / Mesh[i].V;	 
	P = Mesh[i].pressure();
	//cout << rho <<   " " << Mesh[i].getMass() << " " << Mesh[i].V << " " << P << " / " << Mesh[i].getPosx() << " " << Mesh[i].getPosy() << " " << Mesh[i].getPosz() << endl;
	if(P < 0.0) continue;
	ci = sqrt(gamma * P / rho);
	vi = sqrt(sq(Mesh[i].getPx()) + sq(Mesh[i].getPy()) + sq(Mesh[i].getPz())) / Mesh[i].getMass();  
	dbl = min(dbl, Ccfl * Ri / (ci + vi));
      }
      break; 
  }
  return dbl;
}

int compute_gradient(int imesh, double dfdr[5][3], double *w, meshClass *Mesh, int nhydro, int ndim)
{
  const int mass = 0, px = 1, energy = ndim == 3 ? 4 : 3;
  const int buf = ndim == 3 ? 4 : 3;
  const double gamma = 5.0 / 3.0;
  int rval = 0;
  int jmesh = 0;
  double f[ndim], fi[nhydro], ri[ndim], gi[ndim], vi[ndim];
  meshClass *iMesh = &Mesh[imesh], *jMesh = NULL;
 
  ri[0] = iMesh->getPosx();
  gi[0] = iMesh->getxg();
  vi[0] = iMesh->getPx() / iMesh->getMass();
  if(ndim >= 2){ri[1] = iMesh->getPosy(); gi[1] = iMesh->getyg(); vi[1] = iMesh->getPy() / iMesh->getMass();}
  if(ndim >= 3){ri[2] = iMesh->getPosz(); gi[2] = iMesh->getzg(); vi[2] = iMesh->getPz() / iMesh->getMass();} 

  // initialize fi; 
  double vol = ndim == 3 ? iMesh->V : iMesh->S, vi2 = 0.0, kinetic = 0.0, m = iMesh->getMass();
  for(int idim = 0; idim < ndim; ++idim)
  { vi2 += sq(vi[idim]); }
  kinetic = 0.5 * m * vi2 / vol;
  
  fi[mass] = iMesh->getMass();
  for(int idim = 0, ihydro = px; idim < ndim; ++idim, ++ihydro) fi[ihydro] = vi[idim];
  fi[energy] = (gamma - 1.0) * (iMesh->getE() - kinetic);

  // initialize dfdr
  double fmax[nhydro], fmin[nhydro], alpha[nhydro], phi2[nhydro];

  for(int i = 0; i < nhydro; ++i)
  {
    fi[i] = fi[i] - w[i];
    fmax[i] = fi[i];
    fmin[i] = fi[i];
    alpha[i] = 1.0;
  }

  for(int i = 0; i < nhydro; ++i)
    for(int j = 0; j < ndim; ++j)
      dfdr[i][j] = 0.0;

  for(int j = 0; j < iMesh->npoly; ++j)
  {
    jmesh = iMesh->neighbors[j].index; 
    if(jmesh < buf) continue;
    jMesh = &Mesh[jmesh];

    double rj[ndim], gj[ndim], vj[ndim];
    double Aij = iMesh->neighbors[j].Aij, rij[ndim], r = 0.0, cij[ndim];
    double fj[nhydro], vj2 = 0.0;

    // initialize r, g, v
    rj[0] = jMesh->getPosx();
    gj[0] = jMesh->getxg();
    vj[0] = jMesh->getPx() / jMesh->getMass();
    if(ndim >= 2){rj[1] = jMesh->getPosy(); gj[1] = jMesh->getyg(); vj[1] = jMesh->getPy() / jMesh->getMass();}
    if(ndim >= 3){rj[2] = jMesh->getPosz(); gj[2] = jMesh->getzg(); vj[2] = jMesh->getPz() / jMesh->getMass();} 

    vol = ndim == 3 ? jMesh->V : jMesh->S;
    m = jMesh->getMass();

    // calc kinetic
    for(int idim = 0; idim < ndim; ++idim)
    { vj2 += sq(vj[idim]); }
    kinetic = 0.5 * m * vj2 / vol; 
 
    // set fj
    fj[mass] = jMesh->getMass();
    for(int idim = 0, ihydro = px; idim < ndim; ++idim, ++ihydro) fj[ihydro] = vj[idim];
    fj[energy] = (gamma - 1.0) * (jMesh->getE() - kinetic); 

    // calc r
    for(int idim = 0; idim < ndim; ++idim) rij[idim] = ri[idim] - rj[idim];
    
    for(int idim = 0; idim < ndim; ++idim) r += sq(rij[idim]);
    r = sqrt(r);

    for(int idim = 0; idim < ndim; ++idim) f[idim] = jMesh->f[idim];

    // set cij
    for(int idim = 0; idim < ndim; ++idim) cij[idim] = f[idim] - (ri[idim] + rj[idim]) / 2.0;

    // set max min
    for(int i = 0; i < nhydro; ++i)
    {
      fj[i] = fj[i] - w[i];
      fmax[i] = max(fmax[i], fj[i]);
      fmin[i] = min(fmin[i], fj[i]);
    }
    // calc dfdr
    vol = ndim == 3 ? iMesh->V : iMesh->S;
    for(int i = 0; i < nhydro; ++i)
      for(int k = 0; k < ndim; ++k)
	dfdr[i][k] += (Aij * ((fj[i] - fi[i]) * cij[k] / r - (fi[i] + fj[i]) / 2.0 * rij[k] / r)) / vol;

  }
  // calc weight 
  for(int i = 0; i < nhydro; ++i)
  {
    for(int j = 0; j < iMesh->npoly; ++j)
    {
      for(int idim = 0; idim < ndim; ++idim) f[idim] = iMesh->neighbors[j].f[idim];

      phi2[i] = 0.0;
      for(int idim = 0; idim < ndim; ++idim) 
	phi2[i] += dfdr[i][idim] * ( f[idim] - gi[idim] );

      if(phi2[i] > 0.0){
	alpha[i] = min(alpha[i], (fmax[i] - fi[i]) / phi2[i]);
      }
      else if(phi2[i] < 0.0){
	alpha[i] = min(alpha[i], (fmin[i] - fi[i]) / phi2[i]);
      }
      else if(phi2[i] == 0.0)
	alpha[i] = min(alpha[i],1.0);
    }
  }

  // weight alpha
  for(int i = 0; i < nhydro; ++i)
    for(int j = 0; j < ndim; ++j)
    {
      dfdr[i][j] *= alpha[i]*1.0;
    }

  return rval;
}

// riemann_solver //

#define GAMMA    (5.0 / 3.0)
#define GAMMA_G1 ((GAMMA-1.0)/(2.0*GAMMA))
#define GAMMA_G2 ((GAMMA+1.0)/(2.0*GAMMA))
#define GAMMA_G3 ((2.0*GAMMA/(GAMMA-1.0)))
#define GAMMA_G4 (2.0/(GAMMA-1.0))
#define GAMMA_G5 (2.0/(GAMMA+1.0))
#define GAMMA_G6 ((GAMMA-1.0)/(GAMMA+1.0))
#define GAMMA_G7 (0.5*(GAMMA-1.0))
#define GAMMA_G8 (1.0/GAMMA)
#define GAMMA_G9 (GAMMA-1.0)
#define Gamma  5.0 / 3.0
#define phi1(x) sqrt((Gamma - 1.0) * 0.5 + (Gamma + 1.0) * 0.5 * (x))
#define phi2(x) ((Gamma - 1.0) / (2.0 * sqrt(Gamma)) * (1.0 - (x)) / (1.0 - pow((x), (Gamma - 1.0) / (2.0 * Gamma))))
#define TOL_ITER 1.e-6
#define NMAX_ITER 1000
#if 0
void riemann_solver(double *WF, double *uleft, double *uright)
{
  double c, gamma = 5.0 / 3.0;
  double eps, w = 1.25, ML, MR;
  double Pold = 0.0, Pnew, unew, rhonew;
  double Q = 2.0;
  double aL = sqrt(gamma * uleft[3] / uleft[0]), aR = sqrt(gamma * uright[3] / uright[0]);
  double z = (gamma - 1.0) * 0.5 / gamma;
  double PLR = pow(uleft[3] / uright[3], z);
  double CL = uleft[0] * aL, CR = uright[0] * aR;
  double pmin = min(uleft[3], uright[3]), pmax = max(uleft[3], uright[3]);
  double S, ST, SH;
  if(uleft[0] < 0.0 || uleft[3] < 0.0 || uright[0] < 0.0 || uright[3] < 0.0)
    printf("war%f	%f	%f	%f	\n",uleft[0], uleft[3], uright[0], uright[3]);

  Pnew = (CR * uleft[3] + CL * uright[3] + CL * CR * (uleft[1] - uright[1])) / (CL + CR);
  eps = 1.0;

  //while(fabs(eps) > 1.0e-4){
  if(Q > pmax / pmin && pmin <= Pnew && Pnew <= pmax){
    Pold = Pnew;
    Pnew = (CR * uleft[3] + CL * uright[3] + CL * CR * (uleft[1] - uright[1])) / (CL + CR);
    unew = (CL * uleft[1] + CR * uright[1] + uleft[3] - uright[3]) / (CL + CR); 
    if(0.0 <= unew){
      S = uleft[1] - aL;
      if(S <= 0.0)
	rhonew = uleft[0] + (Pnew - uleft[3]) / aL / aL;
      else{
	rhonew = uleft[0];
	unew = uleft[1];
	Pnew = uleft[3];
      }
    }
    if(unew <= 0.0){
      S = uright[1] + aR;
      if(0.0 <= S)
	rhonew = uright[0] + (Pnew - uright[3]) / aR / aR; 
      else{
	rhonew = uright[0];
	unew = uright[1];
	Pnew = uright[3];
      }
    }
    //break;
  }
  else{
    double p0 = max(Pnew, 0.0);
    if(p0 > pmin){
      while(eps >= 1.0e-6){
	ML = sqrt(uleft[0]*uleft[3]) * phi1(p0/uleft[3]);
	MR = sqrt(uright[0]*uright[3]) * phi1(p0/uright[3]);
	Pold = Pnew;
	Pnew = (uleft[1] - uright[1] + uleft[3] / ML + uright[3] / MR) / (1.0 / ML + 1.0 / MR);
	p0 = Pnew;
	eps = fabs(Pnew - Pold);
      }

      unew = 0.5 * (uleft[1] + uright[1] + (Pnew - uright[3]) / MR - (Pnew - uleft[3]) / ML);

      if(0.0 <= unew){
	S = uleft[1] - aL * sqrt((gamma + 1.0) * Pnew / 2.0 / gamma / uleft[3] + (gamma - 1.0) / 2.0 / gamma) ;
	if(S <= 0.0)
	  rhonew = uleft[0] * (Pnew / uleft[3] + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * Pnew / uleft[3] + 1.0) ;
	else{
	  rhonew = uleft[0];
	  unew = uleft[1];
	  Pnew = uleft[3];
	}
      }
      else{
	S = uright[1] + aR * sqrt((gamma + 1.0) * Pnew / 2.0 / gamma / uright[3] + (gamma - 1.0) / 2.0 / gamma);
	if(0.0 <= S)
	  rhonew = uright[0] * (Pnew / uright[3] + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * Pnew / uright[3] + 1.0) ;
	else{
	  rhonew = uright[0];
	  unew = uright[1];
	  Pnew = uright[3];
	}
      }

    }

    else if (p0 <= pmin){
      Pnew = pow((aL + aR - (gamma - 1.0) * (uright[1] - uleft[1]) / 2.0) / (aL / pow(uleft[3], z) + aR / pow(uright[3], z)), 1.0 / z);
      unew = (PLR * uleft[1] / aL + uright[1] / aR + 2.0 * (PLR - 1.0) / (gamma - 1.0)) / (PLR / aL + 1.0 / aR);
      Pold = Pnew;

      if(0.0 <= unew){
	ST = uleft[1] - aL ; 
	SH = unew - aL*pow(Pnew/uleft[3], (gamma - 1.0) / 2.0 / gamma);
	if(0.0 <= SH){
	  rhonew = uleft[0];
	  unew = uleft[1];
	  Pnew = uleft[3];
	}
	else if(SH <= 0.0 && 0.0 <= ST){
	  rhonew = uleft[0] * pow(2.0 / (gamma + 1.0) + (gamma - 1.0) * uleft[1] / (gamma + 1.0) / aL, 2.0 / (gamma - 1.0));
	  unew = 2.0 / (gamma + 1.0) * (aL + (gamma - 1.0) / 2.0 * uleft[1]);
	  Pnew = uleft[3] * pow(2.0 / (gamma + 1.0) + (gamma - 1.0) * uleft[1] / (gamma + 1.0) / aL, 2.0 * gamma / (gamma - 1.0));
	}
	else
	  rhonew = uleft[0] * pow(Pnew / uleft[3], z);
      }
      else{
	ST =  uright[1] + aR;
	SH = unew + aR*pow(Pnew/uright[3], (gamma - 1.0) / 2.0 / gamma); 
	if(0.0 <= ST)
	  rhonew = uright[0] * pow(Pnew / uright[3], z);
	else if(ST <= 0.0 && 0.0 <= SH){
	  rhonew = uright[0] * pow(2.0 / (gamma + 1.0) - (gamma - 1.0) * uright[1] / (gamma + 1.0) / aR, 2.0 / (gamma - 1.0));
	  unew = 2.0 / (gamma + 1.0) * (-aR + (gamma - 1.0) / 2.0 * uright[1]);
	  Pnew = uright[3] * pow(2.0 / (gamma + 1.0) - (gamma - 1.0) * uright[1] / (gamma + 1.0) / aR, 2.0 * gamma / (gamma - 1.0));
	}
	else{
	  rhonew = uright[0];
	  unew = uright[1];
	  Pnew = uright[3];
	}
      }
      //  break;
    }
  } 
  //  }

  WF[0] = rhonew;
  WF[1] = unew; 
  WF[2] = 0.0;
  WF[3] = Pnew;

}
#endif
#if 1
void riemann_solver(double *WF, double *uleft, double *uright, paramsClass &params)
{
  const int density = 0, px = 1, py = 2, pz = 3, energy = params.dim == 3 ? 4 : 3;

  Input_vec_Riemann Riemann_vec;
  Riemann_outputs Riemann_out;
  Riemann_vec.L.rho = uleft[density];
  Riemann_vec.L.u   = uleft[px];
  Riemann_vec.L.p   = uleft[energy]; 
  Riemann_vec.R.rho = uright[density];
  Riemann_vec.R.u   = uright[px];
  Riemann_vec.R.p   = uright[energy];
  double cs_L = sqrt(GAMMA * uleft[energy] / uleft[density]);
  double cs_R = sqrt(GAMMA * uright[energy] / uright[density]);
  double v_line_L = uleft[px];
  double v_line_R = uright[px];

  if((Riemann_vec.L.rho > 0) && (Riemann_vec.R.rho > 0))
  {
    if(iterative_Riemann_solver(Riemann_vec, &Riemann_out, v_line_L, v_line_R, cs_L, cs_R))
    {
      /* this is the 'normal' Reimann solution */
      sample_reimann_standard(0.0,Riemann_vec,&Riemann_out,v_line_L,v_line_R,cs_L,cs_R);
    }
    else
    {
      /* ICs lead to vacuum, need to sample vacuum solution */
      sample_reimann_vaccum_internal(0.0,Riemann_vec,&Riemann_out,v_line_L,v_line_R,cs_L,cs_R);

      if(isnan(WF[energy]) || isnan(WF[px]) || isnan(WF[density]))
	exit(1); 
    }
  } 
  /*  else {
  // one of the densities is zero or negative //
  if((Riemann_vec.L.rho<0)||(Riemann_vec.R.rho<0))
  exit(1);
  if(Riemann_vec.L.rho>0)
  sample_reimann_vaccum_right(0.0,Riemann_vec,&Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
  if(Riemann_vec.R.rho>0)
  sample_reimann_vaccum_left(0.0,Riemann_vec,&Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
  }*/
  WF[density] = Riemann_out.var.rho;
  WF[px] = Riemann_out.var.u;
  WF[py] = WF[pz] = 0.0;
  WF[energy] = Riemann_out.var.p;
}
#endif

double guess_for_pressure(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  double pmin, pmax;
  /* start with the usual lowest-order guess for the contact wave pressure */
  double pv = 0.5*(Riemann_vec.L.p+Riemann_vec.R.p) - 0.125*(v_line_R-v_line_L)*(Riemann_vec.L.p+Riemann_vec.R.p)*(cs_L+cs_R);
  pmin = min(Riemann_vec.L.p,Riemann_vec.R.p);
  pmax = max(Riemann_vec.L.p,Riemann_vec.R.p);

  /* if one side is vacuum, guess half the mean */
  if(pmin<=0)
    return 0.5*(pmin+pmax);

  /* if the two are sufficiently close, and pv is between both values, return it */
  double qrat = pmax / pmin;
  if(qrat <= 2.0 && (pmin <= pv && pv <= pmax))
    return pv;

  if(pv < pmin)
  {
    /* use two-rarefaction solution */
    double pnu = (cs_L+cs_R) - GAMMA_G7 * (v_line_R - v_line_L);
    double pde = cs_L / pow(Riemann_vec.L.p, GAMMA_G1) + cs_R / pow(Riemann_vec.R.p, GAMMA_G1);
    return pow(pnu / pde, GAMMA_G3);
  }
  else
  {
    /* two-shock approximation  */
    double gel = sqrt((GAMMA_G5 / Riemann_vec.L.rho) / (GAMMA_G6 * Riemann_vec.L.p + pv));
    double ger = sqrt((GAMMA_G5 / Riemann_vec.R.rho) / (GAMMA_G6 * Riemann_vec.R.p + pv));
    double x = (gel * Riemann_vec.L.p + ger * Riemann_vec.R.p - (v_line_R - v_line_L)) / (gel + ger);
    if(x < pmin || x > pmax)
      x = pmin;
    return x;
  }
}

int iterative_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  /* before going on, let's compare this to an exact Riemann solution calculated iteratively */
  double Pg,Pg_prev,W_L,W_R,Z_L,Z_R,tol,pratio; int niter_Riemann=0;
  double a0,a1,a2,dvel,check_vel;
  dvel = v_line_R - v_line_L;
  check_vel = GAMMA_G4 * (cs_R + cs_L) - dvel;
  /* if check_vel<0, this will produce a vacuum: need to use vacuum-specific subroutine */
  if(check_vel < 0) return 0;
  if(isnan(cs_L) || isnan(cs_R)) return 0;
  tol=100.0;
  Pg = guess_for_pressure(Riemann_vec, Riemann_out, v_line_L, v_line_R, cs_L, cs_R);
  while((tol>TOL_ITER)&&(niter_Riemann<NMAX_ITER))
  {
    Pg_prev=Pg;
    if(Pg>Riemann_vec.L.p)
    {
      /* shock wave */
      a0 = GAMMA_G5 / Riemann_vec.L.rho;
      a1 = GAMMA_G6 * Riemann_vec.L.p;
      a2 = sqrt(a0 / (Pg+a1));
      W_L = (Pg-Riemann_vec.L.p) * a2;
      Z_L = a2 * (1.0 - 0.5*(Pg-Riemann_vec.L.p)/(a1+Pg));
    } else {
      /* rarefaction wave */
      pratio = Pg / Riemann_vec.L.p;
      W_L = GAMMA_G4 * cs_L * (pow(pratio, GAMMA_G1)-1);
      Z_L = 1 / (Riemann_vec.L.rho*cs_L) * pow(Pg/Riemann_vec.L.p, -GAMMA_G2);
    }
    if(Pg>Riemann_vec.R.p)
    {
      /* shock wave */
      a0 = GAMMA_G5 / Riemann_vec.R.rho;
      a1 = GAMMA_G6 * Riemann_vec.R.p;
      a2 = sqrt(a0 / (Pg+a1));
      W_R = (Pg-Riemann_vec.R.p) * a2;
      Z_R = a2 * (1.0 - 0.5*(Pg-Riemann_vec.R.p)/(a1+Pg));
    } else {
      /* rarefaction wave */
      pratio = Pg / Riemann_vec.R.p;
      W_R = GAMMA_G4 * cs_R * (pow(pratio, GAMMA_G1)-1);
      Z_R = 1 / (Riemann_vec.R.rho*cs_R) * pow(pratio, -GAMMA_G2);
    }
    if(niter_Riemann < NMAX_ITER / 2)
      Pg -= (W_L + W_R + dvel) / (Z_L + Z_R);
    else
      Pg -= 0.5 * (W_L + W_R + dvel) / (Z_L + Z_R);

    if(Pg < 0.1 * Pg_prev)
      Pg = 0.1 * Pg_prev;

    tol = 2.0 * fabs((Pg-Pg_prev)/(Pg+Pg_prev));
    niter_Riemann++;
  }
  if(niter_Riemann<NMAX_ITER)
  {
    Riemann_out->P_M = Pg;
    Riemann_out->S_M = 0.5*(v_line_L+v_line_R) + 0.5*(W_R-W_L);
    return 1;
  } else {
    return 0;
  }
}

void sample_reimann_standard(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  double C_eff,S_eff;
  if(S <= Riemann_out->S_M)  /* sample point is left of contact discontinuity */
  {
    if(Riemann_out->P_M <= Riemann_vec.L.p)	/* left fan (rarefaction) */
    {
      double S_check_L = v_line_L - cs_L;
      if(S <= S_check_L) /* left data state */
      {
	Riemann_out->var.p = Riemann_vec.L.p;
	Riemann_out->var.rho = Riemann_vec.L.rho;
	Riemann_out->var.u = Riemann_vec.L.u;
	return;
      }
      else
      {
	double C_eff_L = cs_L * pow(Riemann_out->P_M / Riemann_vec.L.p, GAMMA_G1);
	double S_tmp_L = Riemann_out->S_M - C_eff_L;

	if(S > S_tmp_L)	/* middle left state */
	{
	  Riemann_out->var.rho = Riemann_vec.L.rho * pow(Riemann_out->P_M / Riemann_vec.L.p, GAMMA_G8);
	  Riemann_out->var.p = Riemann_out->P_M;
	  Riemann_out->var.u = Riemann_vec.L.u + (Riemann_out->S_M-v_line_L);
	  return;
	}
	else		/* left state inside fan */
	{
	  S_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
	  C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
	  Riemann_out->var.rho = Riemann_vec.L.rho * pow(C_eff / cs_L, GAMMA_G4);
	  Riemann_out->var.p = Riemann_vec.L.p * pow(C_eff / cs_L, GAMMA_G3);
	  Riemann_out->var.u = Riemann_vec.L.u + (S_eff-v_line_L);
	  return;
	}
      }
    }
    else			/* left shock */
    {
      if(Riemann_vec.L.p > 0)
      {
	double pml = Riemann_out->P_M / Riemann_vec.L.p;
	double S_L = v_line_L - cs_L * sqrt(GAMMA_G2 * pml + GAMMA_G1);

	if(S <= S_L)	/* left data state */
	{
	  Riemann_out->var.p = Riemann_vec.L.p;
	  Riemann_out->var.rho = Riemann_vec.L.rho;
	  Riemann_out->var.u = Riemann_vec.L.u;
	  return;
	}
	else		/* middle left state behind shock */
	{
	  Riemann_out->var.rho = Riemann_vec.L.rho * (pml + GAMMA_G6) / (pml * GAMMA_G6 + 1.0);
	  Riemann_out->var.p = Riemann_out->P_M;
	  Riemann_out->var.u = Riemann_vec.L.u + (Riemann_out->S_M-v_line_L);
	  return;
	}
      }
      else
      {
	Riemann_out->var.rho = Riemann_vec.L.rho / GAMMA_G6;
	Riemann_out->var.p = Riemann_out->P_M;
	Riemann_out->var.u = Riemann_vec.L.u + (Riemann_out->S_M-v_line_L);
	return;
      }
    }
  }
  else    /* sample point is right of contact discontinuity */
  {
    if(Riemann_out->P_M > Riemann_vec.R.p)	/* right shock */
    {
      if(Riemann_vec.R.p > 0)
      {
	double pmr = Riemann_out->P_M / Riemann_vec.R.p;
	double S_R = v_line_R + cs_R * sqrt(GAMMA_G2 * pmr + GAMMA_G1);

	if(S >= S_R)	/* right data state */
	{
	  Riemann_out->var.p = Riemann_vec.R.p;
	  Riemann_out->var.rho = Riemann_vec.R.rho;
	  Riemann_out->var.u = Riemann_vec.R.u;
	  return;
	}
	else		/* middle right state behind shock */
	{
	  Riemann_out->var.rho = Riemann_vec.R.rho * (pmr + GAMMA_G6) / (pmr * GAMMA_G6 + 1.0);
	  Riemann_out->var.p = Riemann_out->P_M;
	  Riemann_out->var.u = Riemann_vec.R.u + (Riemann_out->S_M-v_line_R);
	  return;
	}
      }
      else
      {
	Riemann_out->var.rho = Riemann_vec.R.rho / GAMMA_G6;
	Riemann_out->var.p = Riemann_out->P_M;
	Riemann_out->var.u = Riemann_vec.R.u + (Riemann_out->S_M-v_line_R);
	return;
      }
    }
    else			/* right fan */
    {
      double S_check_R = v_line_R + cs_R;
      if(S >= S_check_R)		/* right data state */
      {
	Riemann_out->var.p = Riemann_vec.R.p;
	Riemann_out->var.rho = Riemann_vec.R.rho;
	Riemann_out->var.u = Riemann_vec.R.u;
	return;
      }
      else
      {
	double C_eff_R = cs_R * pow(Riemann_out->P_M / Riemann_vec.R.p, GAMMA_G1);
	double S_tmp_R = Riemann_out->S_M + C_eff_R;

	if(S <= S_tmp_R)	/* middle right state */
	{
	  Riemann_out->var.rho = Riemann_vec.R.rho * pow(Riemann_out->P_M / Riemann_vec.R.p, GAMMA_G8);
	  Riemann_out->var.p = Riemann_out->P_M;
	  Riemann_out->var.u = Riemann_vec.R.u + (Riemann_out->S_M-v_line_R);
	  return;
	}
	else		/* fan right state */
	{
	  S_eff = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
	  C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
	  Riemann_out->var.rho = Riemann_vec.R.rho * pow(C_eff / cs_R, GAMMA_G4);
	  Riemann_out->var.p = Riemann_vec.R.p * pow(C_eff / cs_R, GAMMA_G3);
	  Riemann_out->var.u = Riemann_vec.R.u + (S_eff-v_line_R);
	  return;
	}
      }
    }
  }
}

void sample_reimann_vaccum_internal(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  double S_L = v_line_L + GAMMA_G4 * cs_L;
  double S_R = v_line_R - GAMMA_G4 * cs_R;
  if(S <= S_L)
  {
    /* left fan */
    sample_reimann_vaccum_right(S,Riemann_vec,Riemann_out,v_line_L,v_line_R,cs_L,cs_R);
  }
  else if(S >= S_R)
  {
    /* right fan */
    sample_reimann_vaccum_left(S,Riemann_vec,Riemann_out,v_line_L,v_line_R,cs_L,cs_R);
  }
  else
  {
    /* vacuum in between */
    Riemann_out->P_M = 0;
    Riemann_out->S_M = S;
    Riemann_out->var.rho = 0;
    Riemann_out->var.p = Riemann_out->P_M;
    //Riemann_out->var.u = (Riemann_vec.L.u + (Riemann_vec.R.u-Riemann_vec.L.u) * (S-S_L)/(S_R-S_L)) *
    //(1-1) + S * 1;
    Riemann_out->var.u = 0.0;
  }
}


/* --------------------------------------------------------------------------------- */
/* part of exact Riemann solver: */
/* left state is a vacuum, but right state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_left(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  double S_R = v_line_R - GAMMA_G4 * cs_R;

  if(S_R > S)
  {
    /* vacuum */
    Riemann_out->P_M = 0;
    Riemann_out->S_M = S_R;
    Riemann_out->var.rho = 0;
  } else {
    /* right fan */
    double S_R_check = v_line_R + cs_R;
    if(S_R_check > S)
    {
      /* rarefaction fan right state */
      double C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
      Riemann_out->P_M = Riemann_vec.R.p * pow(C_eff / cs_R, GAMMA_G3);
      Riemann_out->S_M = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
      Riemann_out->var.rho = Riemann_vec.R.rho * pow(C_eff / cs_R, GAMMA_G4);
    } else {
      /* right data state */
      Riemann_out->P_M = Riemann_vec.R.p;
      Riemann_out->S_M = v_line_R;
      Riemann_out->var.rho = Riemann_vec.R.rho;
    }
  }
  Riemann_out->var.p = Riemann_out->P_M;
  Riemann_out->var.u = Riemann_vec.R.u + (Riemann_out->S_M - v_line_R);
  return;
}


/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/* right state is a vacuum, but left state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_right(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
    double v_line_L, double v_line_R, double cs_L, double cs_R)
{
  double S_L = v_line_L - GAMMA_G4 * cs_L;

  if(S_L < S)
  {
    /* vacuum */
    Riemann_out->P_M = 0;
    Riemann_out->S_M = S_L;
    Riemann_out->var.rho = 0;
  } else {
    /* left fan */
    double S_L_check = v_line_L - cs_L;
    if(S_L_check < S)
    {
      /* rarefaction fan left state */
      double C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
      Riemann_out->P_M = Riemann_vec.L.p * pow(C_eff / cs_L, GAMMA_G3);
      Riemann_out->S_M = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
      Riemann_out->var.rho = Riemann_vec.L.rho * pow(C_eff / cs_L, GAMMA_G4);
    } else {
      /* left data state */
      Riemann_out->P_M = Riemann_vec.L.p;
      Riemann_out->S_M = v_line_L;
      Riemann_out->var.rho = Riemann_vec.L.rho;
    }
  }
  Riemann_out->var.p = Riemann_out->P_M;
  Riemann_out->var.u = Riemann_vec.L.u + (Riemann_out->S_M - v_line_L);
  return;
}

void compute_w(double *QoldL, double *QoldR, double *f, double *w, double *n, paramsClass &params, meshClass *MeshL, meshClass *MeshR)
{
  const int nhydro = params.dim == 3 ? 5 : 4;
  const int mass = 0, px = 1, energy = params.dim == 3 ? 4 : 3;
  const int ndim = params.dim;
  double rL[3] = {MeshL->getPosx(), MeshL->getPosy(), MeshL->getPosz()}, rR[3] = {MeshR->getPosx(), MeshR->getPosy(), MeshR->getPosz()};
  double vL[3] = {MeshL->getPx() / MeshL->getMass(), MeshL->getPy() / MeshL->getMass(), MeshL->getPz() / MeshL->getMass()}; 
  double vR[3] = {MeshR->getPx() / MeshR->getMass(), MeshR->getPy() / MeshR->getMass(), MeshR->getPz() / MeshR->getMass()};
  //double sL[3] = {MeshL->getxg(), MeshL->getyg(), MeshL->getzg()}, sR[3] = {MeshR->getxg(), MeshR->getyg(), MeshR->getzg()};
  double wL[nhydro], wR[nhydro], wp[nhydro];
  double norm = 0.0;

  for(int i = 0; i < ndim; ++i)
  {
    n[i] = rR[i] - rL[i]; 
    norm += sq(n[i]);
  }

  for(int i = 0; i < ndim; ++i)
  {
    n[i] = n[i] / sqrt(norm);
  }

  double dot  = 0.0, length = 0.0;

  wL[mass] = 0.0;
  for(int i = 0; i < ndim; ++i) 
    dot += vL[i] * n[i];
  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
    wL[ihydro] = dot * n[i];
  wL[energy] = 0.0;

  dot = 0.0;

  wR[mass] = 0.0;
  for(int i = 0; i < ndim; ++i) 
    dot += vR[i] * n[i];
  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
    wR[ihydro] = dot * n[i];
  wR[energy] = 0.0;

  w[mass] = w[energy] = 0.0;

  dot = 0.0;
  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
  {
    dot += (wL[ihydro] - wR[ihydro]) * (f[i] - (rR[i] + rL[i]) / 2.0); 
    length += sq(rR[i] - rL[i]);
  }

  //length = sqrt(length);

  for (int i = px; i < energy; ++i)
  {
    //wp[i] = dot * (rR[i-1] - rL[i-1]) / ( sq(rR[0] - rL[0]) + sq(rR[1] - rL[1]) ) ;
    wp[i] = dot * (rR[i - 1] - rL[i - 1]) / length;
    w[i] = (wR[i] + wL[i]) / 2.0 + wp[i];
  }

  return;
}

void compute_wp(double *QoldL, double *QoldR, meshClass *MeshL, meshClass *MeshR, double *WpL, double *WpR, double *n, double *w, paramsClass &params)
{
  //const int nhydro = params.dim == 3 ? 4 : 3;
  const int mass = 0, px = 1,  energy = params.dim == 3 ? 4 : 3;
  const int ndim = params.dim;
  const double Vl = ndim == 2 ? MeshL->S : MeshL->V;
  const double Vr = ndim == 2 ? MeshR->S : MeshR->V;
  const double gamma = 5.0 / 3.0;
  double kinetic_l = 0.0, kinetic_r= 0.0;
  
  WpL[mass] = QoldL[mass] / Vl;
  WpR[mass] = QoldR[mass] / Vr;

  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
  {
    kinetic_l += sq(QoldL[ihydro]);
    kinetic_r += sq(QoldR[ihydro]);
  }

  kinetic_l = 0.5 * kinetic_l / QoldL[mass];
  kinetic_r = 0.5 * kinetic_r / QoldR[mass];

  WpL[energy] = (gamma - 1.0) * (QoldL[energy] - kinetic_l) / Vl;
  WpR[energy] = (gamma - 1.0) * (QoldR[energy] - kinetic_r) / Vr;

  for(int i = px; i < energy; ++i)
  {
    WpL[i] = QoldL[i] / QoldL[mass] - w[i];
    WpR[i] = QoldR[i] / QoldR[mass] - w[i];
  }

  double dotl = 0.0, dotr = 0.0;

  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
  {
    dotl += n[i] * WpL[ihydro];
    dotr += n[i] * WpR[ihydro];
  }

  for(int i = 0, ihydro = px; i < ndim; ++i, ++ihydro)
  {
    WpL[ihydro] = dotl * n[i];
    WpR[ihydro] = dotr * n[i];
  }

  return;
}

int compute_gradient_t(double *Wp, double dfdr[5][3], double *dfdt, int ndim, int nhydro)
{
  int rval = 0;
  const int num = ndim * nhydro;
  const int mass = 0, px = 1,  energy = ndim == 3 ? 4 : 3;
  const double gamma = 5.0 / 3.0;
  double gradient[num];

  for(int i = 0; i < nhydro; ++i)
    for(int idim = 0; idim < ndim; ++idim)
      gradient[ndim*i+idim] = dfdr[i][idim];

  int igrad = 0;
  double dot1 = 0.0, dot2 = 0.0;
  for(int idim = 0, ihydro = px; idim < ndim; ++idim, ++ihydro, ++igrad)
  {
    dot1 += Wp[ihydro] * gradient[igrad]; 
    dot2 += Wp[mass] * dfdr[ihydro][idim];
  } 
  dfdt[mass] = -dot1 - dot2;
  
  for(int idim = 0, ihydro = px; idim < ndim; ++idim, ++ihydro)
  {
    dot1 = 0.0;
    dot2 = 0.0;

    for(int jdim = 0, jhydro = px; jdim < ndim; ++jdim, ++jhydro, ++igrad)
    {
      dot1 += Wp[jhydro] * gradient[igrad];	
    } 

    dot2 = dfdr[energy][idim] / Wp[mass];

    dfdt[ihydro] = -dot1 - dot2;
  }
  
  dot1 = 0.0; dot2 = 0.0;
  for(int idim = 0, ihydro = px; idim < ndim; ++idim, ++ihydro, ++igrad)
  {
   dot1 += gamma * Wp[energy] * dfdr[ihydro][idim];
   dot2 += Wp[ihydro] * gradient[igrad];
  }
  dfdt[energy] = -dot1 - dot2;

  //dfdt[0] = -(inner_product(1, Wp, 0, gradient) + Wp[0]*(dfdr[1][0]+dfdr[2][1]+dfdr[3][2]));
  //dfdt[1] = -(inner_product(1, Wp, 3, gradient) + dfdr[4][0] / Wp[0]);
  //dfdt[2] = -(inner_product(1, Wp, 6, gradient) + dfdr[4][1] / Wp[0]);
  //dfdt[3] = -(inner_product(1, Wp, 9, gradient) + dfdr[4][2] / Wp[0]);
  //dfdt[4] = -(Gamma*Wp[4]*(dfdr[1][0] + dfdr[2][1] + dfdr[3][2]) + inner_product(1, Wp, 12, gradient));

  return rval;
}



