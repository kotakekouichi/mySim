#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include "params.h"
#include "mesh.h"
#include "particle.h"
#include "def.h"
#include "utils.h"

using namespace std;

int voronoigenerator2d(particleClass *Particle, meshClass *Mesh, paramsClass &params)
{
  int rval = 0;
  const int Nmesh = params.Ngas + 3;

  cout << "voronoi " << endl;

  for(int imesh = 3; imesh < Nmesh; ++imesh)
  {
    map<int, vector<delaunayClass*> > neighbormap;

    for(int id = 0; id < Mesh[imesh].delaunay.size(); ++id)
    {
      for(int im = 0; im < 3; ++im)
      {
	meshClass voronoi;
	int jmesh = Mesh[imesh].delaunay[id]->indm[im];

	if(imesh == jmesh) continue;

	neighbormap[jmesh].push_back(Mesh[imesh].delaunay[id]);
      }	
    } 
	
    int npoly = neighbormap.size();
    Mesh[imesh].npoly = npoly;
    for(map<int, vector<delaunayClass*> >::iterator ite = neighbormap.begin(); 
	ite != neighbormap.end(); ++ite)
    {
      meshClass voronoi;
      voronoi.index = ite->first;
      vector<delaunayClass*> tmpvec = ite->second;
      for(int k = 0; k < tmpvec.size(); ++k)
      {
	voronoi.delaunay.push_back(tmpvec[k]); 
      }
      Mesh[imesh].neighbors.push_back(voronoi);
    }

    Mesh[imesh].calcS2d(); // calculate S, Aij, f;
    Mesh[imesh].calcGravityCenter();
  } 

  return rval;
}

void meshClass::calcGravityCenter()
{
  this->xg = 0.0;
  this->yg = 0.0;
  this->zg = 0.0;
  const int nsize = this->delaunay.size();
  double dblx = 0.0, dbly = 0.0, dblz = 0.0;

  for(int id = 0; id < nsize;++id)
  {
    dblx += this->delaunay[id]->getxc();
    dbly += this->delaunay[id]->getyc();
    dblz += this->delaunay[id]->getzc();
  }
  this->xg = dblx / nsize;
  this->yg = dbly / nsize;
  this->zg = dblz / nsize;
}

void meshClass::calcS2d()
{
  const int Nnbr = this->neighbors.size();
  this->S = 0.0;

  for(int inbr = 0; inbr < Nnbr; ++inbr)
  {
    meshClass *jmesh = &this->neighbors[inbr]; 
    delaunayClass *d1 = jmesh->delaunay[0];
    delaunayClass *d2 = jmesh->delaunay[1];

    double x1 = this->getPosx() - d1->getxc();
    double x2 = this->getPosx() - d2->getxc();
    double y1 = this->getPosy() - d1->getyc();
    double y2 = this->getPosy() - d2->getyc();
     
    jmesh->Aij = sqrt(sq(d1->getxc() - d2->getxc()) + sq(d1->getyc() - d2->getyc()));

    jmesh->f[0] = (d1->getxc() + d2->getxc()) / 2.0;
    jmesh->f[1] = (d1->getyc() + d2->getyc()) / 2.0;

    this->S += 0.5 * fabs(x1 * y2 - x2 * y1);
  }

}

void meshClass::calcV3d()
{

  this->V = this->xg = this->yg = this->zg = 0.0;

  for(int j = 0; j < this->neighbors.size(); ++j)
  {
    meshClass *jmesh = &(this->neighbors[j]);
    const int nsize = jmesh->delaunay.size();

    for(int id = 0; id < nsize; ++id)
    {
      delaunayClass *tmpDelaunay = jmesh->delaunay[id];
#if 0 
      cout << j << " " << id << " ->  " << tmpDelaunay->getxc() << " " << tmpDelaunay->getyc() << " " << tmpDelaunay->getzc() << " " << tmpDelaunay->getRadius()<< endl;

      cout << tmpDelaunay->indm[0] << " " << tmpDelaunay->getPosx(0) << " " << tmpDelaunay->getPosy(0) << " " << tmpDelaunay->getPosz(0) << endl;
      cout << tmpDelaunay->indm[1] << " " << tmpDelaunay->getPosx(1) << " " << tmpDelaunay->getPosy(1) << " " << tmpDelaunay->getPosz(1) << endl;
      cout << tmpDelaunay->indm[2] << " " << tmpDelaunay->getPosx(2) << " " << tmpDelaunay->getPosy(2) << " " << tmpDelaunay->getPosz(2) << endl;
      cout << tmpDelaunay->indm[3] << " " << tmpDelaunay->getPosx(3) << " " << tmpDelaunay->getPosy(3) << " " << tmpDelaunay->getPosz(3) << endl;
    
#endif 
      
    }
  }

  for(int j = 0; j < this->neighbors.size(); ++j)
  {
     meshClass *jmesh = &(this->neighbors[j]);
     jmesh->Aij = 0.0;
     jmesh->f[0] = 0.0;
     jmesh->f[1] = 0.0;
     jmesh->f[2] = 0.0;
        	 
     const int nsize = jmesh->delaunay.size();

     if(nsize == 0) continue;

     for(int id = 0; id < nsize; ++id)
     {
       delaunayClass *tmpDelaunay = jmesh->delaunay[id];
        
       jmesh->f[0] += tmpDelaunay->getxc();
       jmesh->f[1] += tmpDelaunay->getyc();
       jmesh->f[2] += tmpDelaunay->getzc();

     }

     jmesh->f[0] /= nsize; jmesh->f[1] /= nsize; jmesh->f[2] /= nsize;

     double r1[3] = {jmesh->delaunay[0]->getxc() - jmesh->f[0], jmesh->delaunay[0]->getyc() - jmesh->f[1], jmesh->delaunay[0]->getzc() - jmesh->f[2]};

     data_tClass data_t[nsize];
     double nvec[3] = {this->getPosx() - jmesh->getPosx(), this->getPosy() - jmesh->getPosy(), this->getPosz() - jmesh->getPosz()};

     for(int id = 0; id < nsize; ++id)
     {
       double r2[3] = {jmesh->delaunay[id]->getxc() - jmesh->f[0], jmesh->delaunay[id]->getyc() - jmesh->f[1], jmesh->delaunay[id]->getzc() - jmesh->f[2]};
#if 1
       //cout << "r1 " << r1[0] << " " << r1[1] << " " << r1[2] << endl; //cout << "r" <<  id << " " << r2[0] << " " << r2[1] << " " << r2[2] << endl;
       //cout << "center" << jmesh->delaunay[id]->indd << " " << jmesh->delaunay[id]->getxc() << " " << jmesh->delaunay[id]->getyc() << " " << jmesh->delaunay[id]->getzc() << endl;
       //cout << " f " << jmesh->f[0] << " " << jmesh->f[1] << " " << jmesh->f[2] << endl;
       //cout << nsize << endl;
       //cout << index << " " << jmesh->index << endl;
#endif 
       data_t[id].index = id;
       data_t[id].compute_theta(r1, r2, nvec);
       data_t[id].xp = r2[0];
       data_t[id].yp = r2[1];
       data_t[id].zp = r2[2];

     }

     sort(data_t, data_t + nsize);

     double S = 0.0;
     
     for(int ip = 0; ip < nsize; ++ip)
     {
        int jp = ip == nsize - 1 ? 0 :  ip + 1; 

	double lx = data_t[ip].yp * data_t[jp].zp - data_t[ip].zp * data_t[jp].yp; 
        double ly = data_t[ip].zp * data_t[jp].xp - data_t[ip].xp * data_t[jp].zp; 
        double lz = data_t[ip].xp * data_t[jp].yp - data_t[ip].yp * data_t[jp].xp; 
        S += 0.5 * sqrt(lx * lx + ly * ly + lz * lz);	
     }

     double ll = sqrt(sq(nvec[0]) + sq(nvec[1]) + sq(nvec[2]));
     double height = fabs(nvec[0] * (getPosx() - jmesh->delaunay[0]->getxc()) 
	                + nvec[1] * (getPosy() - jmesh->delaunay[0]->getyc())
		        + nvec[2] * (getPosz() - jmesh->delaunay[0]->getzc())) / ll;	
     double Vi = 1.0 / 3.0 * S * height;

     //cout << "heigh and S " << height << " " << S << endl;
     xg += (getPosx() + 3.0 / 4.0 * height / ll * nvec[0]) * Vi;
     yg += (getPosy() + 3.0 / 4.0 * height / ll * nvec[1]) * Vi;
     zg += (getPosz() + 3.0 / 4.0 * height / ll * nvec[2]) * Vi;

     for(int idim = 0; idim < 3; ++idim)
       jmesh->f[idim] = 0.0;
     
     int nsize2 = 0;
     for(int ip = 0; ip < nsize; ++ip)
     {
       int ip1 = ip, ip2 = ip == nsize -1 ? 0 : ip + 1; 
       int id1  = data_t[ip1].index, id2 = data_t[ip2].index;
       delaunayClass *Point1 = jmesh->delaunay[id1], *Point2 = jmesh->delaunay[id2];

       if((fabs(Point1->getxc() - Point2->getxc()) < 1.0e-10) && (fabs(Point1->getyc() - Point2->getyc()) < 1.0e-10) && (fabs(Point1->getzc() - Point2->getzc()) < 1.0e-10))
	 continue;

       jmesh->f[0] += Point1->getxc();
       jmesh->f[1] += Point1->getyc();
       jmesh->f[2] += Point1->getzc();

       ++nsize2;
     }

     if(nsize2 != 0)
     {
       jmesh->f[0] = jmesh->f[0] / nsize2;
       jmesh->f[1] = jmesh->f[1] / nsize2;
       jmesh->f[2] = jmesh->f[2] / nsize2;
     }
#if 0 
     if(0.015625 < this->getPosx() && this->getPosx() < 1.0 - 0.015625 && 0.015625 < this->getPosy() && this->getPosy() < 1. - 0.015625 &&
	0.015625 < this->getPosz() && this->getPosz() < 1.0 - 0.015625)
     {
       cout << jmesh->f[0] << " " << jmesh->f[1] << " " << jmesh->f[2] << " " << S << " " << nsize2 << " " << nsize << " " << jmesh->index << endl;
       cout << jmesh->getPosx() << " " << jmesh->getPosy() << " " << jmesh->getPosz() << endl;
       cout << this->getPosx() << " " << this->getPosy() << " " << this->getPosz() << endl;

       for(int id = 0; id < nsize; ++id)
       {
	 delaunayClass *tmpDelaunay = jmesh->delaunay[id];
	 cout << j << " " << id << " ->  " << tmpDelaunay->getxc() << " " << tmpDelaunay->getyc() << " " << tmpDelaunay->getzc() << " " << tmpDelaunay->getRadius()<< endl;
       }


       cout << endl;
       //getchar();
     }
#endif  
     jmesh->Aij = S;
     V += Vi;
  }

  //cout << "V " << V << " " << cube(1.0 / 16.0) << " " 
    // << getPosx() << " " << getPosy() << " " << getPosz() << endl;

  xg /= V;
  yg /= V;
  zg /= V;
}

int voronoigenerator3d(particleClass *Particle, meshClass *Mesh, paramsClass &params)
{
  int rval = 0;

  const int Nmesh = params.Ngas + 4 + params.Nghost;

  cout << "voronoi " << endl;

  for(int imesh = 4; imesh < Nmesh; ++imesh)
  {

    map<int, vector<delaunayClass*> > neighbormap;

    for(int id = 0; id < Mesh[imesh].delaunay.size(); ++id)
    {
      for(int im = 0; im < 4; ++im)
      {
	int jmesh = Mesh[imesh].delaunay[id]->indm[im];

	//if(jmesh < 4) continue;
	if(imesh == jmesh) continue;

	neighbormap[jmesh].reserve(100);
	neighbormap[jmesh].push_back(Mesh[imesh].delaunay[id]);

	if(!Mesh[imesh].delaunay[id]->flag) getchar();

	//cout << imesh << " " << jmesh << " " << Mesh[imesh].delaunay[id]->indd << endl;
      }	
    } 

    int npoly = neighbormap.size();
    Mesh[imesh].npoly = npoly;
    Mesh[imesh].neighbors.reserve(1000);
    for(map<int, vector<delaunayClass*> >::iterator ite = neighbormap.begin(); 
	ite != neighbormap.end(); ++ite)
    {
      meshClass voronoi;
      voronoi.delaunay.reserve(100);
      voronoi.index = ite->first;
      voronoi.setParticle(Mesh[ite->first].getParticle());

      vector<delaunayClass*> tmpvec = ite->second;
      for(int k = 0; k < tmpvec.size(); ++k)
      {
	voronoi.delaunay.push_back(tmpvec[k]); 
      }
      voronoi.npoly = voronoi.delaunay.size();
      Mesh[imesh].neighbors.push_back(voronoi);

    }

    Mesh[imesh].calcV3d(); // calculate V Aij f xg;

    if(Mesh[imesh].V > 1.0) 
    {
       //cout << "imesh , Vol = " << imesh << " " << Mesh[imesh].V << endl;
       //cout << Mesh[imesh].npoly << endl;
       for(int j = 0; j < Mesh[imesh].neighbors.size(); ++j)
       {
	 meshClass *jmesh = &(Mesh[imesh].neighbors[j]);
	 const int nsize = jmesh->delaunay.size();

	 for(int id = 0; id < nsize; ++id)
	 {
	   delaunayClass *tmpDelaunay = jmesh->delaunay[id];
#if 0 
	   cout << j << " " << id << " ->  " << tmpDelaunay->getxc() << " " << tmpDelaunay->getyc() << " " << tmpDelaunay->getzc() << " " << tmpDelaunay->getRadius()<< " " << tmpDelaunay->birthtype << endl;
	   cout << tmpDelaunay->indm[0] << " " << tmpDelaunay->getPosx(0) << " " << tmpDelaunay->getPosy(0) << " " << tmpDelaunay->getPosz(0) << endl;
	   cout << tmpDelaunay->indm[1] << " " << tmpDelaunay->getPosx(1) << " " << tmpDelaunay->getPosy(1) << " " << tmpDelaunay->getPosz(1) << endl;
	   cout << tmpDelaunay->indm[2] << " " << tmpDelaunay->getPosx(2) << " " << tmpDelaunay->getPosy(2) << " " << tmpDelaunay->getPosz(2) << endl;
	   cout << tmpDelaunay->indm[3] << " " << tmpDelaunay->getPosx(3) << " " << tmpDelaunay->getPosy(3) << " " << tmpDelaunay->getPosz(3) << endl;

	  cout << "s,t,u = " << tmpDelaunay->sss << " " << tmpDelaunay->ttt << " " << tmpDelaunay->uuu << " " << tmpDelaunay->vvv << endl; 
#endif 

	 }
       }
       cout << "XYZ = " << Mesh[imesh].getPosx() << " " << Mesh[imesh].getPosy() << " " << Mesh[imesh].getPosz() << " " << Mesh[imesh].V << endl; 

    }
    //Mesh[imesh].calcGravityCenter();
  } 

  //getchar();

  return rval;
}
