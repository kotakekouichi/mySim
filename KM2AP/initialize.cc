#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "mesh.h"
#include "particle.h"
#include "params.h"
#include "def.h"
#include "utils.h"
#include "Vec.h"
#include "exArith.h"


int initializeDelaunay(delaunayClass *Delaunay, paramsClass &params);
double orient3d(Vec& a, Vec& b, Vec& c, Vec& d);

int initialize2d(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{
  const int Ngas = params.Ngas;
  const int Ndm = params.Ndm;
  const int Npart = Ngas + Ndm + 3;
  int rval = 0;

  //rval = initializeDelaunay(Delaunay, params);
  // dark matter particle
  for(int i = 3; i < Ndm + 3; ++i) 
  {
    Particle[i].setParticleType(dark_matter);
  }

  // gas particle 
#if 0
  const int Ngrid = 140;
  for(int i = 0, imesh = 3, ipart = 3; i < Ngrid; ++i)
  {
    for(int j = 0; j < Ngrid; ++j, ++imesh, ++ipart)
    {
      double dx = LEFTSIDE + i * BOXSIZE / Ngrid + BOXSIZE / Ngrid / 2.0; 
      double dy = LEFTSIDE + j * BOXSIZE / Ngrid + BOXSIZE / Ngrid / 2.0;

      Particle[ipart].setParticleType(gas);
      Particle[ipart].setPosx(dx);
      Particle[ipart].setPosy(dy);
      Mesh[imesh].setParticle(&Particle[ipart]);
      Mesh[imesh].delaunay.reserve(100);
      Mesh[imesh].index = imesh;
      Mesh[imesh].next = NULL;

      Mesh[imesh].setMass( 1.0 / Ngrid / Ngrid);
      Mesh[imesh].setPx(0.0);
      Mesh[imesh].setPy(0.0);

      double xc = (LEFTSIDE + RIGHTSIDE) / 2.0, yc = (LEFTSIDE + RIGHTSIDE) / 2.0;
      double dblE = sq(Mesh[imesh].getPosx() - xc) + sq(Mesh[imesh].getPosy() - yc) < sq(16.0 / (double) Ngrid) ? 
	exp(-10.0 * (sq(Mesh[imesh].getPosx() - xc) + sq(Mesh[imesh].getPosy() - yc))) : 1.0 / Ngrid / Ngrid; 
      Mesh[imesh].setE(dblE); 
    }
  }
#endif
#if 1
  for(int i = 3 + Ndm,  imesh = 3; i < Npart; ++i, ++imesh)
  {
    Particle[i].setParticleType(gas);
    double rnd = (double) rand() / RAND_MAX;
    Particle[i].setPosx(rnd);
    rnd = (double) rand() / RAND_MAX;
    Particle[i].setPosy(rnd);
    Mesh[imesh].setParticle(&Particle[i]); 
    Mesh[imesh].delaunay.reserve(100);
    Mesh[imesh].index = imesh;
    Mesh[imesh].next = NULL;
  }  
#endif

  if(params.periodic)
  {
    int nghost = 0;
    const int Nmesh = 3 + params.Ngas;

    for(int ipart = 3 + Ndm, imesh = 3; ipart < Npart; ++ipart, ++imesh)
    {
      for(int i = -1; i <= 1; ++i)
      {
	for(int j = -1; j <= 1; ++j)
	{
	  if(i == 0 && j == 0) continue;

	  double xp = i * 1.0 + Mesh[imesh].getPosx();
	  double yp = j * 1.0 + Mesh[imesh].getPosy();

	  if(!(-DX < xp && xp < BOXSIZE + DX && -DX < yp && yp < BOXSIZE + DX ) ) continue;

	  Particle[Npart + nghost].setParticleType(dummy);
	  Particle[Npart + nghost].setPosx(xp);
	  Particle[Npart + nghost].setPosy(yp);

	  Mesh[Nmesh + nghost].setParticle(&Particle[Npart + nghost]);
	  Mesh[Nmesh + nghost].index = Nmesh + nghost;
	  Mesh[Nmesh + nghost].next = NULL;

	  ++nghost;
	}
      } 
    }

    params.Nghost = nghost;
  }

  return rval;

}

int initialize3d(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{
  int rval = 0;
  const int Ngas = params.Ngas;
  const int Ndm = params.Ndm;
  const int Npart = Ngas + Ndm + 4;
  //const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE;

  cout << "initialize3d " << endl;
  // dummy particle
  for(int i = 0, imesh = 0; i < 4; ++i, ++imesh)
  {
    Particle[i].setParticleType(dummy);
  }

  // dark matter particle
  if(params.dmfile != "")
  {
    int i = 4, id = 0;
    double mp = 0.0, xp = 0.0 , yp = 0.0, zp = 0.0, vx =0.0 , vy = 0.0, vz = 0.0;
    FILE *fp = fopen(params.dmfile.c_str(), "r");

    while(fscanf(fp, "%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf", &id,&mp, &xp, &yp, &zp, &vx, &vy, &vz) != EOF)
    {
      Particle[i].setParticleType(dark_matter);
      Particle[i].setMass(mp); 
      Particle[i].setPosx(xp);
      Particle[i].setPosy(yp);
      Particle[i].setPosz(zp);
      Particle[i].setVelx(vx);
      Particle[i].setVely(vy);
      Particle[i].setVelz(vz);
      i++;
    }

    fclose(fp);

  } else
  {
    for(int i = 4; i < Ndm + 4; ++i) 
    {
      //...
      Particle[i].setParticleType(dark_matter);
    }
  }

  if(params.gasfile != "")
  {
    int i = params.Ndm + 4, id = 0, imesh = 4;
    double mp = 0.0, xp = 0.0 , yp = 0.0, zp = 0.0, vx =0.0 , vy = 0.0, vz = 0.0, U, E;
    FILE *fp = fopen(params.gasfile.c_str(), "r");

    //while(fscanf(fp, "%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf", &id,&mp, &xp, &yp, &zp, &vx, &vy, &vz, &E, &U) != EOF)
    while(fscanf(fp, "%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf", &id,&mp, &xp, &yp, &zp, &vx, &vy, &vz, &E, &U) != EOF)

    {
      Particle[i].setParticleType(gas);
      Particle[i].setMass(mp); 
      Particle[i].setPosx(xp);
      Particle[i].setPosy(yp);
      Particle[i].setPosz(zp);
      Particle[i].setVelx(vx);
      Particle[i].setVely(vy);
      Particle[i].setVelz(vz);
      Particle[i].setPx(mp * vx);
      Particle[i].setPy(mp * vy);
      Particle[i].setPz(mp * vz);
      Particle[i].setMeshGenerationPoint(&Mesh[imesh]);
      Mesh[imesh].setParticle(&Particle[i]);
      Mesh[imesh].index = imesh;
      Mesh[imesh].delaunay.reserve(100);
      Mesh[imesh].next = NULL;
      i++; ++imesh;
    }

    fclose(fp);

  }else 
  {
#if 0  
    //RANDOM
    for(int i = 4 + Ndm, imesh = 4; i < Ngas + Ndm + 4; ++i, ++imesh) 
    {
      Particle[i].setParticleType(gas);
      
      double rnd = (double) rand() / RAND_MAX;
      Particle[i].setPosx(rnd);
      
      rnd = (double) rand() / RAND_MAX;
      Particle[i].setPosy(rnd);

      rnd = (double) rand() / RAND_MAX;
      Particle[i].setPosz(rnd);

      Particel[i].setMeshGenerationPoint(&Mesh[imesh]);
      Mesh[imesh].setParticle(&Particle[i]); 
      Mesh[imesh].delaunay.reserve(200);
      Mesh[imesh].index = imesh;
      Mesh[imesh].next = NULL;
      Mesh[imesh].neighbors.reserve(200);
    }
#endif
#if 1
    //GRID
    int Ngrid = params.Ngrid;
    for(int i = 0, imesh = 4, ipart = 4; i < Ngrid; ++i)
    {
      for(int j = 0; j < Ngrid; ++j)
      {
	for(int k = 0; k < Ngrid; ++k, ++imesh, ++ipart)
	{
	  double dx = LEFTSIDE + i * BOXSIZE / Ngrid + BOXSIZE / Ngrid / 2.0; 
	  double dy = LEFTSIDE + j * BOXSIZE / Ngrid + BOXSIZE / Ngrid / 2.0;
	  double dz = LEFTSIDE + k * BOXSIZE / Ngrid + BOXSIZE / Ngrid / 2.0;

	  Particle[ipart].setParticleType(gas);

	  Particle[ipart].setMass( 1.0 / (double) cube(Ngrid) );
	  Particle[ipart].setPosx(dx);
	  Particle[ipart].setPosy(dy);
	  Particle[ipart].setPosz(dz);
	  Particle[ipart].setMeshGenerationPoint(&Mesh[imesh]);
      
	  Mesh[imesh].setParticle(&Particle[ipart]);
	  Mesh[imesh].setMass(1.0 / (double) cube(Ngrid));
	  Mesh[imesh].delaunay.reserve(100);
	  Mesh[imesh].index = imesh;
	  Mesh[imesh].next = NULL;
	  
	  double xc = (LEFTSIDE + RIGHTSIDE) / 2.0, yc = (LEFTSIDE + RIGHTSIDE) / 2.0, zc = (LEFTSIDE + RIGHTSIDE) / 2.0;
	  double length2 = sq(Mesh[imesh].getPosx() - xc) + sq(Mesh[imesh].getPosy() - yc) + sq(Mesh[imesh].getPosz() - zc);
	  //double dblE = length2 < sq(16.0 / (double) Ngrid) ? exp(-10.0 * (length2)): 1.0 / Ngrid / Ngrid / Ngrid; 
	  double dblE = length2 < sq(10.0 / 32.0) ?  exp(-30.0 * length2) : 1.0 / cube(32.0);
	  Mesh[imesh].setE(dblE);
	}
      }
    }
#endif
  }

  if(params.periodic)
  {
    int nghost = 0;
    const int Nmesh = 4 + params.Ngas;

    for(int ipart = 4 + Ndm, imesh = 4; ipart < Npart; ++ipart, ++imesh)
    {
      for(int i = -1; i <= 1; ++i)
      {
	for(int j = -1; j <= 1; ++j)
	{
	  for(int k = -1; k <= 1; ++k)
	  {
	    if(i == 0 && j == 0 && k == 0) continue;

	    double xp = i * 1.0 + Mesh[imesh].getPosx();
	    double yp = j * 1.0 + Mesh[imesh].getPosy();
	    double zp = k * 1.0 + Mesh[imesh].getPosz();

	    if(!(-DX < xp && xp < BOXSIZE + DX && -DX < yp && yp < BOXSIZE + DX && -DX < zp && zp < BOXSIZE + DX ) ) continue;
	  
	    Particle[Npart + nghost].setParticleType(dummy);
	    Particle[Npart + nghost].setPosx(xp);
	    Particle[Npart + nghost].setPosy(yp);
	    Particle[Npart + nghost].setPosz(zp);

	    Mesh[Nmesh + nghost].setParticle(&Particle[Npart + nghost]);
	    Mesh[Nmesh + nghost].index = Nmesh + nghost;
	    Mesh[Nmesh + nghost].next = NULL;

	    ++nghost;

	  }
	}
      } 
    }

    params.Nghost = nghost;
  }

  cout << "end initialize " << endl;

  return rval;

}

int initializeDelaunay(particleClass *Particle, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{

  int rval = 0;
  const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE;

  delaunayClass TopDelaynay;
  
  if(params.dim == 2)
  {
    double r_c = sqrt(boxsize * boxsize / 2.0 / 2.0 * 2.0);
    double cx = (RIGHTSIDE + LEFTSIDE) / 2.0;
    double cy = (RIGHTSIDE + LEFTSIDE) / 2.0;

    TopDelaynay.setPosx(0, cx - sqrt(3.0) * r_c);
    TopDelaynay.setPosy(0, cy - r_c);
    TopDelaynay.setPosx(1, cx + sqrt(3.0) * r_c);
    TopDelaynay.setPosy(1, cy - r_c);
    TopDelaynay.setPosx(2, cx);
    TopDelaynay.setPosy(2, cy + 2.0 * r_c);
    TopDelaynay.indd = 0;
    TopDelaynay.indm[0] = 0;
    TopDelaynay.indm[1] = 1;
    TopDelaynay.indm[2] = 2;

    //Delaunay.reserve(10000000);
    Delaunay[0] = TopDelaynay;

  } else if(params.dim == 3)
  {
    const double r = sqrt(pow(boxsize / 2.0, 2.0) * 3.0);
    const double R = 12.0 * r / sqrt(6.0);
    const double xc = (LEFTSIDE + RIGHTSIDE) / 2.0, yc = (LEFTSIDE + RIGHTSIDE) / 2.0, zc = (LEFTSIDE + RIGHTSIDE) / 2.0;

    TopDelaynay.setPosx(0, xc); TopDelaynay.setPosy(0, yc); TopDelaynay.setPosz(0, sqrt(6.0) / 4.0 * R + zc);
    TopDelaynay.setPosx(1, 1.0 / 2.0 * R + xc); TopDelaynay.setPosy(1, -sqrt(3.0) / 6.0 * R + yc); TopDelaynay.setPosz(1, -sqrt(6.0) / 12.0 * R + zc);
    TopDelaynay.setPosx(2, -1.0 / 2.0 * R + xc); TopDelaynay.setPosy(2, -sqrt(3.0) / 6.0 * R + yc); TopDelaynay.setPosz(2, -sqrt(6.0) / 12.0 * R + zc);
    TopDelaynay.setPosx(3, xc); TopDelaynay.setPosy(3, 1.0 / sqrt(3.0) * R + yc); TopDelaynay.setPosz(3, -sqrt(6.0) / 12.0 * R + zc);
    
    TopDelaynay.indd = 0;
    TopDelaynay.indm[0] = 0;
    TopDelaynay.indm[1] = 1;
    TopDelaynay.indm[2] = 2;
    TopDelaynay.indm[3] = 3;
    TopDelaynay.flag = true;
    
    for(int i = 0; i < 4; ++i)
    {
      TopDelaynay.neighbors[i] = NULL;
      TopDelaynay.faces[i] = NULL;
    }

    double a[3] = {TopDelaynay.getPosx(0), TopDelaynay.getPosy(0), TopDelaynay.getPosz(0)};
    double b[3] = {TopDelaynay.getPosx(1), TopDelaynay.getPosy(1), TopDelaynay.getPosz(1)};
    double c[3] = {TopDelaynay.getPosx(2), TopDelaynay.getPosy(2), TopDelaynay.getPosz(2)};
    double d[3] = {TopDelaynay.getPosx(3), TopDelaynay.getPosy(3), TopDelaynay.getPosz(3)}; 
    Vec avec12; 
    Vec bvec12 ;
    Vec cvec12 ;
    Vec dvec12 ;
    Vec evec12 ;
    avec12.setVec12(a[0], a[1], a[2]);
    bvec12.setVec12(b[0], b[1], b[2]);
    cvec12.setVec12(c[0], c[1], c[2]);
    dvec12.setVec12(d[0], d[1], d[2]);
    TopDelaynay.orient = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

    Delaunay[0] = TopDelaynay;
  }

  // dummy particle
  const int buf = params.dim == 2 ? 3 : 4;
  for(int i = 0, imesh = 0; i < buf; ++i, ++imesh)
  {
    Particle[i].setParticleType(dummy);
    Particle[i].setPosx(TopDelaynay.getPosx(i));
    Particle[i].setPosy(TopDelaynay.getPosy(i));
    if(params.dim == 3)
    {
      Particle[i].setPosz(TopDelaynay.getPosz(i));
    }
    Mesh[imesh].setParticle(&Particle[i]);
    Mesh[imesh].delaunay.reserve(200);
    Mesh[imesh].delaunay.push_back(&Delaunay[0]);
    Mesh[imesh].delaunaymap[0] = &Delaunay[0];
  }

  return rval;

}
