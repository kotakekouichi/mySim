#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include "params.h"
#include "mesh.h"
#include "particle.h"
#include "def.h"
#include "utils.h"
#include "tree.h"
#include "Vec.h"
#include "exArith.h"

#define _DEBUG
using namespace std;

clock_t time3, time4, time5, time6, time7, time8, time9, gtmp;
static Vec avec12; 
static Vec bvec12 ;
static Vec cvec12 ;
static Vec dvec12 ;
static Vec evec12 ;
static delaunayClass oldDelaunay[50];
static vector<delaunayClass*> tmpStackDelaunays;
static delaunayClass** stackDelaunays;
static int istack;
static int lastDel;
static int men[2][3];
  
static vector<faceClass*> facevec;
faceClass *staticfaces;
static int iiface;

/* ------ function ------- */
void hilbertC(Vec *hilbertCurves, int &n, int s, double x, double y, double z, double dx, double dy, double dz, double dx2, double dy2, double dz2, double dx3, double dy3, double dz3);
void settingNextMesh(meshClass *Mesh, Vec *hilbertCurves, int n, paramsClass &params);

/* --- 2dimension --- */
int delaunaysplit2d(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params);
int delaunayflip2d(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params);
/* ------------------  */

/* --- 3dimension --- */
int delaunaysplit3d1to4(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params);
int delaunaysplit3dnto2n(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params);
int delaunayflip3d(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params);
int delaunayflip3d2to3(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params);
int delaunayflip3d3to2(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params);
int delaunayflip3d4to4(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params, int no_indm[2]);
int ch4to4(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params);

double orient3d(Vec& a, Vec& b, Vec& c, Vec& d);
double insphere(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e);

bool coplanar(int A, int B, int imesh, int jmesh, meshClass *Mesh);
bool onEdge(int A, int B, int C, meshClass *Mesh, double &eps);
int getDelaunaySharingEdge(int *no_indm, int &nsize, int imesh, int it, delaunayClass *oldDelaunay, delaunayClass *Delaunay, meshClass *Mesh);

/* ------------------  */

//int check_flip(int it, int imesh, meshClass *Mesh, delaunayClass    *Delaunay, delaunayClass *delaunayNew[], int Num, paramsClass &params);
int check_flip(int it, int imesh, meshClass *Mesh, delaunayClass    *Delaunay, int *newIndex, int Num, paramsClass &params);

void outputMesh(meshClass *Mesh, paramsClass &params);
void outputDelaunay(delaunayClass *Delaunay, paramsClass &params);
/* ----------------------------------------  */

/* --- test global parameter -------- */
int iii;
double _sss;
double _ttt;
double _uuu;
double _vvv;
bool _2to6Flag;

/*--- 2 dimenstion ---*/
int delaunaygenerator2d(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, paramsClass &params)
{
  cout << "start delaunay generator " << endl;
  int rval = 0;
  const int Nmesh = params.Ngas + 3 + params.Nghost;
  const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE;

  clock_t time = clock(), time1=0, time2=0;

  tree<meshClass> Root(0, (LEFTSIDE + RIGHTSIDE) / 2.0, (LEFTSIDE + RIGHTSIDE) / 2.0, boxsize, 10);
  Root.setdim(2);

  params.Ndelaunay = 1;

  int conetr = 0;
  for(int imesh = 3; imesh < Nmesh; ++imesh, ++conetr)
  {
    int it = 0;
    //cout << "imesh = "  << imesh << endl;
   
    splitType dummy;
    it = Root.delaunayfinder(imesh, Mesh, dummy);
    //it = -1;
    
    if(it < 0) 
    {
      for(it = 0; it < params.Ndelaunay; ++it)
      {
	if(delaunay[it].crossdeterminant(Mesh[imesh].getPosx(), Mesh[imesh].getPosy()) == 0)
	  break;
      }
    }

    rval = delaunaysplit2d(it, imesh, Mesh, delaunay, params);

  } 

  cout << (double) time1 / CLOCKS_PER_SEC << " " << (double) time2 / CLOCKS_PER_SEC << endl;
  time = clock() - time;
  cout << "clock time[sec] = " << (double) time / CLOCKS_PER_SEC << endl;
  
  //outputDelaunay(delaunay, params);
  //outputMesh(Mesh, params);
  Root.free();

  return rval;
} 

int delaunayClass::crossdeterminant(double dx, double dy)
{
  int rval = 0;
  int a, b, c;
  double det, A, B, C;
  double s,t;
  double x0 = this->getPosx(0), x1 = this->getPosx(1), x2 = this->getPosx(2);
  double y0 = this->getPosy(0), y1 = this->getPosy(1), y2 = this->getPosy(2);

  det = DET2D((x0 - x2) ,(x1 - x2) ,(y0 - y2) ,(y1 - y2));
  s = ((y1 - y2) * (dx - x2) - (x1 - x2) * (dy - y2)) / det;
  t = (-(y0 - y2) * (dx - x2) + (x0 - x2) * (dy - y2)) / det;

  if(0.0 <= s && s < 1.0 && 0.0 <= t && t < 1.0 && s + t <= 1.0)
  {
    rval = 0; 
    return rval;
  }

  A = DET2D(dx-x0,dy-y0,x2-x0,y2-y0);
  a = int(A / fabs(A));
  B = DET2D(dx-x1,dy-y1,x0-x1,y0-y1);
  b = int(B / fabs(B));
  C = DET2D(dx-x2,dy-y2,x1-x2,y1-y2);
  c = int(C / fabs(C));
  
  if(a == b && b == c)
    //det = -1.0;
    rval = 0;
  else 
    //det = 1.0; 
    rval = 1;
  return rval;
}

int delaunaysplit2d(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{
  //cout << "start split" << endl;
  int rval = 0;
  //vector<int> NewIndex;
  int NewIndex[3];
  //vector<delaunayClass> delaunayNew;
  delaunayClass *delaunayNew[3];
  delaunayClass aOldDelaunay = Delaunay[it];
  delaunayNew[0] = &Delaunay[it];
  delaunayNew[1] = &Delaunay[params.Ndelaunay];
  delaunayNew[2] = &Delaunay[params.Ndelaunay+1];

  delaunayNew[0]->nNeighbor = 0;

  for(int i = 0, j = 0; j < 2; ++j)
  {
    for(int k = j + 1; k < 3; ++k, ++i)
    {
      //int index = i == 0 ? Delaunay[it].indd : Delaunay.size() + i - 1;
      int index = i == 0 ? Delaunay[it].indd :  params.Ndelaunay + i - 1;
      //delaunayClass tmpDelaunay(index, Delaunay[it].indm[j], Delaunay[it].indm[k], imesh, Mesh); 
      //NewIndex.push_back(index);
      NewIndex[i] = index;
      //delaunayNew.push_back(tmpDelaunay);
      //delaunayNew[i] = tmpDelaunay;
      delaunayNew[i]->set(index, aOldDelaunay.indm[j], aOldDelaunay.indm[k], imesh, Mesh);
    }
  }

  params.Ndelaunay += 2;
/*
  for(int i = 0; i < 2; ++i)
  {
    delaunayClass tmp;
    Delaunay.push_back(tmp);
  }
  */
/*
  delaunayNew[0]->neighbors.push_back(&Delaunay[NewIndex[1]]); delaunayNew[0]->neighbors.push_back(&Delaunay[NewIndex[2]]);
  delaunayNew[1]->neighbors.push_back(&Delaunay[NewIndex[0]]); delaunayNew[1]->neighbors.push_back(&Delaunay[NewIndex[2]]);
  delaunayNew[2]->neighbors.push_back(&Delaunay[NewIndex[0]]); delaunayNew[2]->neighbors.push_back(&Delaunay[NewIndex[1]]);
*/
  delaunayNew[0]->neighbors[0] = &Delaunay[NewIndex[1]]; delaunayNew[0]->neighbors[1] = &Delaunay[NewIndex[2]];
  delaunayNew[0]->nNeighbor = 2;

  delaunayNew[1]->neighbors[0] = &Delaunay[NewIndex[0]]; delaunayNew[1]->neighbors[1] = & Delaunay[NewIndex[2]];
  delaunayNew[1]->nNeighbor = 2;

  delaunayNew[2]->neighbors[0] = &Delaunay[NewIndex[0]]; delaunayNew[2]->neighbors[1] = & Delaunay[NewIndex[1]];
  delaunayNew[2]->nNeighbor = 2;

  for(int i = 0, j = 0; j < 2; ++j)
  {
    int jmesh = aOldDelaunay.indm[j];
    for(int k = j + 1; k < 3; ++k, ++i)
    {

      delaunayClass *newDel = delaunayNew[i];
      int kmesh = aOldDelaunay.indm[k];
      for(int inbr = 0, nsize = aOldDelaunay.nNeighbor ; inbr < nsize; ++inbr)
      {
	delaunayClass *tmpNeighbor = aOldDelaunay.neighbors[inbr];
        if(tmpNeighbor->has_mesh(jmesh) && tmpNeighbor->has_mesh(kmesh))
	{
	  int index = 0;

	  for(int l = 0; l < tmpNeighbor->nNeighbor; ++l)
	  {
	    if(tmpNeighbor->neighbors[l]->indd == aOldDelaunay.indd) 
	    {
	      index = l; 
	      break;
	    } 
	  }
	  //delaunayNew[i]->neighbors.push_back(tmpNeighbor);
	  newDel->neighbors[ newDel->nNeighbor ] = tmpNeighbor; 
	  newDel->nNeighbor += 1;
	  tmpNeighbor->neighbors[index] = &Delaunay[NewIndex[i]]; 
	  //delaunayNew[i]->neighbors.push_back(oldDelaunay.neighbors[inbr]);
	  break;
	}	  
      }
    }
  }

  rval = aOldDelaunay.eraseDelaunaryOfMesh(Mesh, params);
/*
  for(int i = 0; i < 3; ++i)
  {
    Delaunay[NewIndex[i]] = delaunayNew[i];
  }
*/
  for(int i = 0; i < 3; ++i)
  {
    //int index = delaunayNew[i].indd;
    //rval = Delaunay[index].addDelaunayToMesh(Mesh);
    rval = delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }

  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 3, params);

  rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 3, params);
  return rval;
}

int delaunayClass::changeToNewDelaunay()
{
  int rval = 0;
  return rval;
} 

void delaunayClass::set(int index, int imesh, int jmesh, int kmesh, meshClass *Mesh)
{
  this->indd = index;

  this->indm[0] = imesh;
  this->indm[1] = jmesh;
  this->indm[2] = kmesh;
  this->indm[3] = -1;

  this->xp[0] = Mesh[imesh].getPosx();
  this->yp[0] = Mesh[imesh].getPosy();

  this->xp[1] = Mesh[jmesh].getPosx();
  this->yp[1] = Mesh[jmesh].getPosy();

  this->xp[2] = Mesh[kmesh].getPosx();
  this->yp[2] = Mesh[kmesh].getPosy();
  //this->neighbors.reserve(16);

  this->nNeighbor = 0;

  this->flag = true;
  this->checkFlag = false;

}

void delaunayClass::set(int index, int imesh, int jmesh, int kmesh, int lmesh, meshClass *Mesh)
{
  this->indd = index;

  this->indm[0] = imesh;
  this->indm[1] = jmesh;
  this->indm[2] = kmesh;
  this->indm[3] = lmesh;

  this->xp[0] = Mesh[imesh].getPosx();
  this->yp[0] = Mesh[imesh].getPosy();
  this->zp[0] = Mesh[imesh].getPosz();

  this->xp[1] = Mesh[jmesh].getPosx();
  this->yp[1] = Mesh[jmesh].getPosy();
  this->zp[1] = Mesh[jmesh].getPosz();

  this->xp[2] = Mesh[kmesh].getPosx();
  this->yp[2] = Mesh[kmesh].getPosy();
  this->zp[2] = Mesh[kmesh].getPosz();

  this->xp[3] = Mesh[lmesh].getPosx();
  this->yp[3] = Mesh[lmesh].getPosy();
  this->zp[3] = Mesh[lmesh].getPosz();
  //this->neighbors.reserve(16);

  //this->neighbors.clear();
  for(int i = 0; i < 4; ++i) this->neighbors[i] = NULL;
  for(int i = 0; i < 4; ++i) this->faces[i] = NULL;

  this->flag = true;
  this->checkFlag = false;
  this->nNeighbor = 0;

  this->radius = 0.0;

  double a[3] = {this->getPosx(0), this->getPosy(0), this->getPosz(0)};
  double b[3] = {this->getPosx(1), this->getPosy(1), this->getPosz(1)};
  double c[3] = {this->getPosx(2), this->getPosy(2), this->getPosz(2)};
  double d[3] = {this->getPosx(3), this->getPosy(3), this->getPosz(3)}; 
  avec12.setVec12(a[0], a[1], a[2]);
  bvec12.setVec12(b[0], b[1], b[2]);
  cvec12.setVec12(c[0], c[1], c[2]);
  dvec12.setVec12(d[0], d[1], d[2]);
  this->orient = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
}

//int check_flip(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, delaunayClass *delaunayNew[], int Num, paramsClass &params) 
int check_flip(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, int *NewIndex, int Num, paramsClass &params) 
{
  /*
  int rval = 0, idl, jmesh;
  const int nth = params.dim == 2 ? 3 : 4;
  clock_t tmp1 = 0, tmp2 = 0, a = 0;
  meshClass *jm_ptr = NULL;
  tmp1 = clock();
  
  for(int i = 0; i < Num; ++i)
  {
    //delaunayClass *tmpDelaunay = &Delaunay[delaunayNew[i]->indd];
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];

    if(!tmpDelaunay->flag || tmpDelaunay->checkFlag) continue;
    if(tmpDelaunay->nNeighbor != nth) continue;

      faceClass *face3 = tmpDelaunay->faces[3];
      delaunayClass *tmpNbr = face3->interDel[0];
      
      idl = 0;
      if(tmpNbr->indd == tmpDelaunay->indd)
      {tmpNbr = face3->interDel[1]; idl = 1;}

      jmesh = face3->interidx[idl];
      jmesh = tmpNbr->indm[jmesh];
      jm_ptr = &Mesh[jmesh];

      switch(params.dim)
      {
	case(2):
	  rval = Delaunay[tmpDelaunay->indd].determinant2d(Mesh[jmesh].getPosx(), Mesh[jmesh].getPosy());
	  if(rval == 0) continue;
	  rval = delaunayflip2d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, Delaunay, params);
	  goto NEXT;

	  break;

	case(3):

	  //tmp2 = clock();
	  rval = Delaunay[tmpDelaunay->indd].determinant3d(jm_ptr->getPosx(), jm_ptr->getPosy(), jm_ptr->getPosz());
	  //tmp2 = clock() - tmp2;
	  //a += tmp2;

	  //Delaunay[tmpDelaunay->indd].checkFlag = true;

	  if(rval == 0) continue;

	  tmp2 = clock();
	  rval = delaunayflip3d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, Delaunay, params);

	  tmp2 = clock() - tmp2;
	  a += tmp2;
	  goto NEXT;

	  break;
      }
NEXT:;

  }

  tmp1 = clock() - tmp1;

  tmp1 = tmp1 - a;

  //time5 += tmp1;
    
  return rval;

  */

  int rval = 0;
  const int nth = params.dim == 2 ? 3 : 4;
#if 0
 map<double, int> sortedMap;

  data_tClass *data_t = new data_tClass[Num];

  for(int i = 0; i < Num; ++i)
  {
    data_t[i].index = Delaunay[NewIndex[i]].indd; 
    //if(calcFlag)
      data_t[i].theta = 1.0 / (Delaunay[NewIndex[i]].getRadius());
    data_t[i].theta = i;
  }
  ///calcFlag = false;
  sort(data_t, data_t + Num);

  for(int i = 0; i < Num; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
    stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
  }
  return 0;
#endif 

  for(int i = 0; i < Num; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
    //delaunayClass *tmpDelaunay = &Delaunay[data_t[i].index];

    if(!tmpDelaunay->flag) continue;
    
    if(tmpDelaunay->nNeighbor != nth) continue;

    for(int inbr = 0; inbr < tmpDelaunay->nNeighbor; ++inbr)
    {
      delaunayClass *tmpNbr = tmpDelaunay->neighbors[inbr];
      
      if(tmpNbr->has_mesh(imesh)) continue;
      for(int j = 0; j < nth; ++j)
      {
	int jmesh = tmpNbr->indm[j]; 
	if(tmpDelaunay->has_mesh(jmesh)) continue;

	switch(params.dim)
	{
	  case(2):
	    rval = Delaunay[tmpDelaunay->indd].determinant2d(Mesh[jmesh].getPosx(), Mesh[jmesh].getPosy());
	    if(rval == 0) continue;
	    rval = delaunayflip2d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, Delaunay, params);
	    goto NEXT;
	    
	    break;

	  case(3):
	    //cout << "check  "<< i << " " << tmpDelaunay->indd << " " <<  jmesh << " " << tmpNbr->birthtype << endl;
	    rval = Delaunay[tmpDelaunay->indd].determinant3d(Mesh[jmesh].getPosx(), Mesh[jmesh].getPosy(), Mesh[jmesh].getPosz());
	    //cout << "endcheck" << endl;

	    if(rval == 0) continue;
	    rval = delaunayflip3d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, Delaunay, params);

	    //stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
	    goto NEXT;

	    break;
	}
	//goto NEXT;

      }
    }	
NEXT:;

  }

    
  return rval;

}

int delaunayflip2d(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params)
{
  int rval = 0;
  delaunayClass *delaunayNew[2];
  //delaunayClass oldDelaunay[2];

  oldDelaunay[0] = Delaunay[it];
  oldDelaunay[1] = Delaunay[in];
  delaunayNew[0] = &Delaunay[it];
  delaunayNew[1] = &Delaunay[in];

  delaunayNew[0]->nNeighbor = 0;
  delaunayNew[1]->nNeighbor = 0; 

  for(int i = 0; i < 2; ++i)
  {
    int index = i == 0 ? Delaunay[it].indd : Delaunay[in].indd;
    delaunayNew[i]->set(index, oldDelaunay[0].indm[i], jmesh,imesh, Mesh);
  }
  
  delaunayNew[0]->neighbors[0] = &Delaunay[in];
  delaunayNew[0]->nNeighbor = 1;
  delaunayNew[1]->neighbors[0] = &Delaunay[it];
  delaunayNew[1]->nNeighbor = 1;

  for(int i = 0; i < 2; ++i)
  {
    int idxd = i;

    delaunayClass *tmpDelaunay = &oldDelaunay[idxd];
    for(int j = 0; j < tmpDelaunay->nNeighbor; ++j)
    {
      delaunayClass *tmpNbr = tmpDelaunay->neighbors[j];
      if(tmpNbr->indd == delaunayNew[0]->indd || 
         tmpNbr->indd == delaunayNew[1]->indd) continue;
      
      int index = 0;
      for(int l = 0; l < tmpNbr->nNeighbor; ++l)
      {
	if(tmpNbr->neighbors[l]->indd == tmpDelaunay->indd) 
	{
	  index = l; break;
	} 
      }

      if(tmpNbr->has_mesh(oldDelaunay[0].indm[0]))
      {
	tmpNbr->neighbors[index] = delaunayNew[0];
	delaunayNew[0]->neighbors[delaunayNew[0]->nNeighbor] = tmpNbr;
	delaunayNew[0]->nNeighbor += 1;
      }
      else 
      {
	tmpNbr->neighbors[index] = delaunayNew[1];
	delaunayNew[1]->neighbors[delaunayNew[1]->nNeighbor] = tmpNbr;
	delaunayNew[1]->nNeighbor += 1;
      }
    }
  }

  /*
  for(int i = 0; i < 2; ++i)
  {
    oldDelaunay[i].eraseDelaunaryOfMesh(Mesh, params);
  }
  for(int i = 0; i < 2; ++i)
  {
    delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }*/

  int NewIndex[2] = {delaunayNew[0]->indd, delaunayNew[1]->indd};
  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 2, params);
  rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 2, params);

  return rval;
}

int delaunayClass::determinant2d(double xp, double yp)
{
  int rval = 0;
  double x0 = this->getPosx(0), x1 = this->getPosx(1), x2 = this->getPosx(2);
  double y0 = this->getPosy(0), y1 = this->getPosy(1), y2 = this->getPosy(2);

  rval = circumcenter(x0, y0, x1, y1, x2, y2, &xc, &yc, &radius);

  if(sq(xc - xp) + sq(yc - yp) < sq(radius))
    rval = 1;
  else 
    return rval = 0;
  return rval;
}

int delaunayClass::determinant3d(double _xp, double _yp, double _zp)
{
  int rval = 0;
  double a[3] = {this->getPosx(0), this->getPosy(0), this->getPosz(0)};
  double b[3] = {this->getPosx(1), this->getPosy(1), this->getPosz(1)};
  double c[3] = {this->getPosx(2), this->getPosy(2), this->getPosz(2)};
  double d[3] = {this->getPosx(3), this->getPosy(3), this->getPosz(3)}; 
  double e[3] = {_xp, _yp, _zp}; 
  double Orient, Insphere; 
  //clock_t tmptime = 0;

  if(this->radius == 1.0e+10){return 1; }

  int old = -1;
  
    //double eps2 = radius*radius - (sq(xc - _xp)+ sq(yc - _yp) + sq(zc - _zp));
  if(this->V > 1.0e-13)
  {

    double eps = radius - sqrt(sq(xc - _xp) + sq(yc - _yp) + sq(zc - _zp));
    if(eps == 0.0 || eps > 1.0e-6)
    {
      rval = 1;
      //old = 1;
      goto CLEANUP;
    }
    else if(eps < -1.0e-6)
    {
      rval = 0;
      //old = 0;
      goto CLEANUP; 
    }
  }

  //tmptime = clock();
  avec12.setVec12(a[0], a[1], a[2]);
  bvec12.setVec12(b[0], b[1], b[2]);
  cvec12.setVec12(c[0], c[1], c[2]);
  dvec12.setVec12(d[0], d[1], d[2]);
  evec12.setVec12(e[0], e[1], e[2]);
  //Orient = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
  Orient = this->orient;

  if(0 < Orient)
  {
    Insphere = ExactArithmetic::insphere(avec12, bvec12, cvec12, dvec12, evec12);
    if(0 <= Insphere)
      rval = 1;
    else 
      rval = 0;

    /*if(old != rval && old != -1 && !(fabs(eps) < 1.0e-4))
    {
      cout << Insphere << " " << Orient << " " << eps << " " << V << endl;
    }*/
  } else if(0 > Orient)
  {
    Insphere = ExactArithmetic::insphere(avec12, bvec12, cvec12, dvec12, evec12);
    if(0 >= Insphere)
      rval = 1;
    else 
      rval = 0;
/*
    if(old != rval && old != -1 && !(fabs(eps) < 1.0e-4))
    {
      cout << Insphere << " " << Orient << " " << eps << " " << V << endl;
    }
*/
  } else 
  {
    rval = 1;
/*
    if(old != rval && old != -1 && !(fabs(eps) < 1.0e-4))
    {
      cout << Insphere << " " << Orient << " " << eps << " " << V << endl;
    }
*/
  }

  //tmptime = clock() - tmptime;
  //time8 += tmptime;
CLEANUP:

  return rval;
}


int delaunayClass::circumcenter3d()
{

#define L(a, b,c) ((a) * (a) + (b) * (b) + (c) *(c))
  double det;
  double r2[3], a[3], b[3], c[3];
  //a_inv[3], b_inv[3], c_inv[3];
  for(int i=0;i<3;i++){
    a[i] = xp[i+1] - xp[0];
    b[i] = yp[i+1] - yp[0];
    c[i] = zp[i+1] - zp[0];
    r2[i] = L(xp[i+1], yp[i+1], zp[i+1]) - L(xp[0], yp[0], zp[0]);
  }

  det = DET3D(a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2]);

/* circumcenter */
#if 0
  xc = 0.0;yc =0.0;zc=0.0;
  a_inv[0] = (b[1]*c[2] - b[2]*c[1]); a_inv[1] = (b[0]*c[2]-b[2]*c[0]) * (-1.0) ; a_inv[2] = (b[0]*c[1] - b[1]*c[0]);
  b_inv[0] = -(a[1]*c[2] - a[2]*c[1]); b_inv[1] = (a[0]*c[2]-a[2]*c[0])         ; b_inv[2] = -(a[0]*c[1] - a[1]*c[0]);
  c_inv[0] = (a[1]*b[2] - a[2]*b[1]); c_inv[1] = (a[0]*b[2]-a[2]*b[0]) * (-1.0) ; c_inv[2] = (a[0]*b[1] - a[1]*b[0]);

  for(int i=0;i<3;i++){
    xc += r2[i] * a_inv[i] / (det * 2.0);  
    yc += r2[i] * b_inv[i] / (det * 2.0);  
    zc += r2[i] * c_inv[i] / (det * 2.0);  
  }

  radius = sqrt(L((xc - xp[0]), (yc - yp[0]) , (zc - zp[0])));

  if(fabs((det)) == 0.0){
    this->radius = 1.0e+10;
  }

  if(fabs(det) < 1.0e-10 && _2to6Flag) 
  {
    this->radius = 1.0e+10;
    //_2to6Flag = false;
  }
#endif

  double xc2 = 0.0, yc2 = 0.0, zc2 = 0.0;
  double r1[3], r3[3];

  r1[0] = a[0];
  r1[1] = b[0];
  r1[2] = c[0];

  r2[0] = a[1];
  r2[1] = b[1];
  r2[2] = c[1];

  r3[0] = a[2];
  r3[1] = b[2];
  r3[2] = c[2];

  double fac1 = r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2];
  double fac2 = r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2];
  double fac3 = r3[0] * r3[0] + r3[1] * r3[1] + r3[2] * r3[2];
  double R[3];
  R[0] = fac3 * (r1[1] * r2[2] - r1[2] * r2[1]) +
    fac2 * (r3[1] * r1[2] - r3[2] * r1[1]) +
    fac1 * (r2[1] * r3[2] - r2[2] * r3[1]);
  R[1] = fac3 * (r1[2] * r2[0] - r1[0] * r2[2]) +
    fac2 * (r3[2] * r1[0] - r3[0] * r1[2]) +
    fac1 * (r2[2] * r3[0] - r2[0] * r3[2]);
  R[2] = fac3 * (r1[0] * r2[1] - r1[1] * r2[0]) +
    fac2 * (r3[0] * r1[1] - r3[1] * r1[0]) +
    fac1 * (r2[0] * r3[1] - r2[1] * r3[0]);
  double _V = 2. * (r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
      r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
      r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
  if(fabs(_V) < 1.0e-13)
  {
    xc2 = 0.25 * (xp[0] + xp[1] + xp[2] + xp[3]);
    yc2 = 0.25 * (yp[0] + yp[1] + yp[2] + yp[3]);
    zc2 = 0.25 * (zp[0] + zp[1] + zp[2] + zp[3]); 
  } else {
    double Vinv = 1.0 / _V;
    xc2 = xp[0] + R[0] * Vinv;
    yc2 = yp[0] + R[1] * Vinv;
    zc2 = zp[0] + R[2] * Vinv;
  }
    
  xc = xc2;
  yc = yc2;
  zc = zc2;
  this->V = _V;
  radius = sqrt(L((xc - xp[0]), (yc - yp[0]) , (zc - zp[0])));

  double test1 = this->orient;
  if(_2to6Flag && test1 == 0.0) 
    this->radius = 1.0e+10;
  
  return 0;

}

int delaunayClass::circumcenter(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double *_xc, double *_yc, double *_radius)
{
  double det,c1,c2;
  double MATRIX[2][2];

  det = (Bx - Ax) * (Cy - Ay) - (Cx - Ax) * (By - Ay);
   if(det == 0.0){
     *_xc = 0.0;
     *_yc = 0.0;
     *_radius = 1.0e+2;
     return 0;
   }
   c1 = Bx * Bx + By * By - Ax * Ax - Ay * Ay;
   c2 = Cx * Cx + Cy * Cy - Ax * Ax - Ay * Ay;
 
   MATRIX[0][0] = 1.0 / 2.0 / det * (Cy - Ay);
   MATRIX[0][1] = 1.0 / 2.0 / det * (Ay - By);
   MATRIX[1][0] = 1.0 / 2.0 / det * (Ax - Cx); 
   MATRIX[1][1] = 1.0 / 2.0 / det * (Bx - Ax);
 
   *_xc = MATRIX[0][0] * c1 + MATRIX[0][1] * c2;
   *_yc = MATRIX[1][0] * c1 + MATRIX[1][1] * c2;
   *_radius = sqrt((*_xc - Ax)*(*_xc - Ax)+(*_yc -Ay)*(-Ay+*_yc));
   return 0;
 }

int delaunayClass::eraseDelaunaryOfMesh(meshClass *Mesh, paramsClass &params)
{
  int rval = 0;
  const int buf = params.dim == 2 ? 3 : 4;
#if 0
  for(int im = 0; im < buf; ++im)
  {
    int idxmesh = this->indm[im], id = 0;
    //if(idxmesh < buf) continue;
    bool findflag = false;
    meshClass *mesh = &Mesh[idxmesh];
    for(id = mesh->delaunay.size() - 1; id >= 0; --id)
    {
      if(mesh->delaunay[id]->indd == this->indd)
      {findflag = true; break;}
    } 
    if(findflag) mesh->delaunay.erase(mesh->delaunay.begin() + id);   
  }
#endif
  for(int im = 0; im < buf; ++im)
  {
    int idxmesh = this->indm[im];
    if(idxmesh < buf) continue;
    meshClass *mesh = &Mesh[idxmesh];
    unordered_map<int, delaunayClass*>::iterator itr = mesh->delaunaymap.find(this->indd);
    if(itr != mesh->delaunaymap.end())
      mesh->delaunaymap.erase(itr); 
  }
  return rval;
}

int delaunayClass::addDelaunayToMesh(meshClass *Mesh, paramsClass &params)
{
  int rval = 0;

  const int buf = params.dim == 2 ? 3 : 4;
  for(int im = 0; im < buf; ++im)
  {
    int idxmesh = this->indm[im]; 
    //Mesh[idxmesh].delaunay[this->indd] = this;
    Mesh[idxmesh].delaunay.push_back(this);
  } 

  return rval;

}
/*-------------------*/

/*--- 3 dimenstion ---*/
int delaunaygenerator3d(particleClass *Particle, meshClass *Mesh, delaunayClass *delaunay, paramsClass &params)
{
  int rval = 0;

  cout << "start delaunay generator " << endl;
  const int Nmesh = params.Ngas + 4 + params.Nghost;
  const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE;

  clock_t time1 = 0;
  clock_t time2 = 0;
  time3 = 0.0; time4 = 0.0, time5 = 0.0, time6 = 0, time7 = 0, time8 = 0;
  //clock_t time1=0, time2=0;

  params.Ndelaunay = 1;

  // == test hilbert == //
#define HILBERT
#ifdef HILBERT
  const int nGrid = 16;
  const int nGrid3 = cube(nGrid);
/*
  for(int i = 0;i < facevec.size(); ++i) 
  {
    delete facevec[i];
  }
  facevec.clear();
  facevec.reserve(1000000);
 */
  tmpStackDelaunays.reserve(100);
  staticfaces = new faceClass[10000000];
  iiface = 0; 
  lastDel = 0;
  int n = 0;
  Vec *hilbertCurves = new Vec[nGrid3];
  double x = 0, y = 0, z = 0;
  double dx = 1, dy = 0, dz = 0;
  double dx2 = 0, dy2 = 1, dz2 = 0;
  double dx3 = 0, dy3 = 0 ,dz3 = 1;
  
  hilbertC(hilbertCurves, n, nGrid, x, y, z, dx, dy, dz, dx2, dy2, dz2, dx3, dy3, dz3);
  settingNextMesh(Mesh, hilbertCurves, nGrid, params);

  clock_t time = clock();
  int meshCounter = 0;
  //stackDelaunays.reserve(10000);
  stackDelaunays = new delaunayClass*[10000];
  
  for(int i = 0; i < nGrid3; ++i)
  {
     for(meshClass *mptr = hilbertCurves[i].first; mptr != NULL; mptr = mptr->next, ++meshCounter)
     {
       int it = -1;
       //clock_t tmptime = clock(), tmptime2 = clock();
       istack = 0;

       //cout << "start" << endl;
       _2to6Flag = false;
       splitType splittype = _nosplit;// 
       //splittype = delaunay[params.Ndelaunay - 1].getSplitType_WalkNeighbor(it, mptr, params);
       if(delaunay[lastDel].flag)
	 splittype = delaunay[lastDel].getSplitType_WalkNeighbor(it, mptr, params);

       if(splittype == _nosplit)
       {
	 for(it = params.Ndelaunay -1; it >= 0; --it)
	 {
	   if(!delaunay[it].flag) continue;
	   splittype = delaunay[it].getSplitType(mptr->getPosx(), mptr->getPosy(), mptr->getPosz());
	   if(splittype != _nosplit) 
	   {
	     break;
	   }
	 }
       }

       //tmptime2 = clock() - tmptime2;
       //time1 += tmptime2;

       int imesh = mptr->index;
       if(meshCounter % 5000 == 0)
       cout  << "mesh => " << meshCounter << " " << Nmesh << " " << 
	 params.Ngas << " " << splittype << " " << params.Ndelaunay << " " << iiface << endl;
       //tmptime = clock();
       if(!delaunay[it].flag ) 
       {cout << "warning flag not found[" << it << "] " << endl;getchar();}

       if(splittype == _1to4)
	 rval = delaunaysplit3d1to4(it, imesh, Mesh, delaunay, params);
       else if(splittype == _nto2n)
	 rval = delaunaysplit3dnto2n(it, imesh, Mesh, delaunay, params);

       for(int is = 0; is < tmpStackDelaunays.size(); ++is, ++istack)
	 stackDelaunays[istack] = tmpStackDelaunays[is];
       tmpStackDelaunays.clear();

       //while(stackDelaunays.size() > 0)
       //tmptime = clock();
       while(istack > 0)
       {
	 //int s = stackDelaunays.size() - 1; 
	 int s = istack - 1;
	 delaunayClass *tmpDelaunay = stackDelaunays[s];

	 //clock_t tmp2 = clock();
	 if(tmpDelaunay->flag && tmpDelaunay->nNeighbor == 4) 
	 {
	   faceClass *face3 = tmpDelaunay->faces[3];
	   delaunayClass *tmpNbr = face3->interDel[0];

	   int idl = 0;
	   if(tmpNbr->indd == tmpDelaunay->indd)
	   {tmpNbr = face3->interDel[1]; idl = 1;}

	   int jmesh = face3->interidx[idl];
	   jmesh = tmpNbr->indm[jmesh];
	   meshClass* jm_ptr = &Mesh[jmesh];

	   rval = delaunay[tmpDelaunay->indd].determinant3d(jm_ptr->getPosx(), jm_ptr->getPosy(), jm_ptr->getPosz());

	 //tmpDelaunay->checkFlag = true;
	   if(rval != 0)
	   {

	     rval = delaunayflip3d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, delaunay, params);
	   }

	 }
	 //a += clock() - tmp2;
	 //stackDelaunays.pop_back();

	 istack--;

	 for(int is = 0; is < tmpStackDelaunays.size(); ++is, ++istack)
	 {
	   stackDelaunays[istack] = tmpStackDelaunays[is];
	 }
	 tmpStackDelaunays.clear();

       }

       //tmptime = clock() - tmptime;
       //time8 += tmptime - a;

     }
  } 
/*
	for(int inbr = 0; inbr < tmpDelaunay->nNeighbor; ++inbr)
	{
	delaunayClass *tmpNbr = tmpDelaunay->neighbors[inbr];

	if(tmpNbr->has_mesh(imesh)) continue;
	for(int j = 0; j < 4; ++j)
	{
	int jmesh = tmpNbr->indm[j]; 
	if(tmpDelaunay->has_mesh(jmesh)) continue;

	rval = delaunay[tmpDelaunay->indd].determinant3d(Mesh[jmesh].getPosx(), Mesh[jmesh].getPosy(), Mesh[jmesh].getPosz());

	if(rval != 0) 
	rval = delaunayflip3d(tmpDelaunay->indd, imesh, jmesh, Mesh, tmpNbr->indd, delaunay, params);

	goto NEXT;

	}
	}*/

  delete [] hilbertCurves;

  time = clock() - time;
  cout << "clock time[sec] = " << (double) time / CLOCKS_PER_SEC << endl;
  cout << "clock time1[sec] = " << (double) time1 / CLOCKS_PER_SEC << endl;
  cout << "clock time2[sec] = " << (double) time2 / CLOCKS_PER_SEC << endl;
  cout << "clock time3[sec] = " << (double) time3 / CLOCKS_PER_SEC << endl;
  cout << "clock time4[sec] = " << (double) time4 / CLOCKS_PER_SEC << endl;
  cout << "clock time5[sec] = " << (double) time5 / CLOCKS_PER_SEC << endl;
  cout << "clock time6[sec] = " << (double) time6 / CLOCKS_PER_SEC << endl;
  cout << "clock time7[sec] = " << (double) time7 / CLOCKS_PER_SEC << endl;
  cout << "clock time8[sec] = " << (double) time8 / CLOCKS_PER_SEC << endl;

  for(int i = 0; i < params.Ndelaunay; ++i)
  { 
    if(!delaunay[i].flag) continue;
    rval = delaunay[i].addDelaunayToMesh(Mesh, params);
  }

  getchar();
  
  delete [] staticfaces;
  delete [] stackDelaunays;
#endif
  // ================ //

/*

  tree<meshClass> Root(0, (LEFTSIDE + RIGHTSIDE) / 2.0, (LEFTSIDE + RIGHTSIDE) / 2.0, (LEFTSIDE + RIGHTSIDE) / 2.0, boxsize, 10);
  Root.setdim(3);
  for(int imesh = 4; imesh < Nmesh; ++imesh)
  {
    int it = -1;
 
    clock_t tmptime = clock();
    
    _2to6Flag = false;
    splitType splittype;
    iii = imesh;
    it = Root.delaunayfinder(imesh, Mesh, splittype);
      
    //it = -1;

    if(it < 0 || delaunay[it].flag == false) 
    {
      for(it = 0; it < params.Ndelaunay; ++it)
      {
	if(!delaunay[it].flag) continue;
	
	splittype = delaunay[it].getSplitType(Mesh[imesh].getPosx(), Mesh[imesh].getPosy(), Mesh[imesh].getPosz());
	
	if(splittype != _nosplit) 
	{
	  break;
	}
      }
    }

    tmptime = clock() - tmptime;
    time1 += tmptime;
    

    cout  << " " << imesh << " " << Nmesh << " " << params.Ngas << " " << splittype << endl;
    tmptime = clock();
    if(!delaunay[it].flag ) {cout << "warning flag not found[" << it << "] " << endl;getchar();}

    if(it == params.Ndelaunay)
    { 
      cout << "not found " << endl;
      getchar();
      continue;
    }

    if(splittype == _1to4)
      rval = delaunaysplit3d1to4(it, imesh, Mesh, delaunay, params);
    else if(splittype == _nto2n)
      rval = delaunaysplit3dnto2n(it, imesh, Mesh, delaunay, params);
    tmptime = clock() - tmptime;
    time2 += tmptime;
  } 

  cout << params.Ndelaunay << endl;
  //getchar();

  outputMesh(Mesh, params);
  
  for(int id = 0; id < params.Ndelaunay; ++id)
  {
    delaunay[id].circumcenter3d();
  }

  time = clock() - time;
  cout << "clock time[sec] = " << (double) time / CLOCKS_PER_SEC << endl;
  cout << "clock time1[sec] = " << (double) time1 / CLOCKS_PER_SEC << endl;
  cout << "clock time2[sec] = " << (double) time2 / CLOCKS_PER_SEC << endl;
  
  cout << "clock time3[sec] = " << (double) time3 / CLOCKS_PER_SEC << endl;
  cout << "clock time4[sec] = " << (double) time4 / CLOCKS_PER_SEC << endl;
  cout << "clock time5[sec] = " << (double) time5 / CLOCKS_PER_SEC << endl;
  cout << "clock time6[sec] = " << (double) time6 / CLOCKS_PER_SEC << endl;
  cout << "clock time7[sec] = " << (double) time7 / CLOCKS_PER_SEC << endl;

  getchar();
  Root.free();
  */
  //outputDelaunay(delaunay, params);
  cout << "output file " << endl;
  //outputMesh(Mesh, params);
  return rval;
}

splitType delaunayClass::getSplitType(double dx, double dy, double dz)
{
  splitType stype = _nosplit;
  double s = 0.0, t = 0.0, u = 0.0, p[3] = {dx - getPosx(0), dy - getPosy(0), dz - getPosz(0)};
  int _sign = 1;
#if 1 
  
  double det;
  double a[3], b[3], c[3];
  double a_inv[3], b_inv[3], c_inv[3];

  p[0] = dx - xp[0];
  p[1] = dy - yp[0];
  p[2] = dz - zp[0];
  for(int i=0;i<3;i++){
    a[i] = xp[i+1] - xp[0];
    b[i] = yp[i+1] - yp[0];
    c[i] = zp[i+1] - zp[0];
  }
  det = DET3D(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
  a_inv[0] = (b[1]*c[2] - b[2]*c[1]); a_inv[1] = (b[0]*c[2]-b[2]*c[0]) * (-1.0)      ; a_inv[2] = (b[0]*c[1] - b[1]*c[0]);
  b_inv[0] = -(a[1]*c[2] - a[2]*c[1]); b_inv[1] = (a[0]*c[2]-a[2]*c[0])              ; b_inv[2] = -(a[0]*c[1] - a[1]*c[0]);
  c_inv[0] = (a[1]*b[2] - a[2]*b[1]); c_inv[1] = (a[0]*b[2]-a[2]*b[0]) * (-1.0)      ; c_inv[2] = (a[0]*b[1] - a[1]*b[0]);
  s = (a_inv[0] * p[0] + b_inv[0] * p[1] + c_inv[0] * p[2]) / det;
  t = (a_inv[1] * p[0] + b_inv[1] * p[1] + c_inv[1] * p[2]) / det;
  u = (a_inv[2] * p[0] + b_inv[2] * p[1] + c_inv[2] * p[2]) / det;
  
  if(0.0 <= s && s <= 1.0 && 0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0 && s + t + u <= 1.0 + 1.0e-14)
  {
    _sss = s;
    _ttt = t;
    _uuu = u;
    _vvv = 1.0 - s - t- u;
    //stype = _1to4; 
  } 

  int count1 = 0;
  int count2 = 0;
  //(0,1,2,3) (0,1,3,2), (0,2,3,1) (1,2,3,0)
  for(int i = 0, idx = 3; i < 4; ++i)
      {
        for(int j = i + 1;j < 4; ++j)
	{
           for(int k = j + 1; k < 4; ++k, --idx)
	   {
	   
	     avec12.setVec12(xp[i], yp[i], zp[i]);
	     bvec12.setVec12(xp[j], yp[j], zp[j]);
	     cvec12.setVec12(xp[k], yp[k], zp[k]);
	     dvec12.setVec12(dx, dy, dz);

	     double test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

	     dvec12.setVec12(xp[idx], yp[idx], zp[idx]);

	     //double test2 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
	     double test2 = this->orient * _sign; _sign *= -1;

	     //cout << idx << " " << test2 << " " << this->orient << endl;
	     //getchar();

	     if(sign(test1) == sign(test2) && test1 != 0)
	       ++count1; 

	     if(test1 == 0)
	     {
	       men[count2][0] = indm[i]; men[count2][1] = indm[j]; men[count2][2] = indm[k];
	       ++count2;
	     }

	   }
	}
      }

  if(count1 == 4){
    stype = _1to4; 
    //cout << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s -t -u) << endl;
    //cout << "1to4" << endl;
    //getchar();
  }

  if(count2 == 1 && count1 == 3)
  {
    stype = _1to4;
    _2to6Flag = true;
    //cout << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s -t -u) << endl;
    //cout << "2to6" << endl;
    //getchar();
  } else if(count2 == 2 && count1 == 2)
  {
    stype = _nto2n;
    //cout << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s -t -u) << endl;

    //cout << "nto2n " << endl;
    //getchar();
  }
  
/*
  //const double eps2 = 1.0e-8;
  const double eps2 = 1.0e-5;
  if(-eps2 <= s && -eps2 <= t && -eps2 <= u && -eps2 <= 1.0-(s+t+u))
      if(stype != _nto2n && ((-eps2 <= s && s < eps2) || (-eps2 <= t && t < eps2)
	    || (-eps2 <= u && u < eps2) || (1.0-eps2 < s+t+u && s+t+u <= 1.0 + eps2)))
      {

	//cout <<  "2to6? " << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s -t -u) << endl;

	int Num0 = 0;

	for(int i = 0; i < 4; ++i)
	{
	  for(int j = i + 1;j < 4; ++j)
	  {
	    for(int k = j + 1; k < 4; ++k)
	    {

	      avec12.setVec12(xp[i], yp[i], zp[i]);
	      bvec12.setVec12(xp[j], yp[j], zp[j]);
	      cvec12.setVec12(xp[k], yp[k], zp[k]);
	      dvec12.setVec12(dx, dy, dz);

	      double test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
	      if(test1 == 0.0)
		++Num0;
	    }
	  }
	}
	//cout <<  "2to6? " << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s - t -u) << endl;
	if(Num0 == 1) 
	{
	  //stype = _1to4;
	  _2to6Flag = true;
	} 
      }
	
  const double epsnto2n = 1.0e-4;
  if((fabs(s) < epsnto2n && fabs(t) < epsnto2n && 0.0 < u && u < 1.0) ||
	(fabs(t) < epsnto2n && fabs(u) < epsnto2n && 0.0 < s && s < 1.0) ||
	(fabs(u) < epsnto2n && fabs(s) < epsnto2n && 0.0 < t && t < 1.0) ||
	(fabs(s) < epsnto2n && 0.0 < t && t < 1.0 && 0.0 < u && u < 1.0 && fabs(t + u - 1.0) < epsnto2n) ||
	(fabs(t) < epsnto2n && 0.0 < s && s < 1.0 && 0.0 < u && u < 1.0 && fabs(s + u - 1.0) < epsnto2n) || 
	(fabs(u) < epsnto2n && 0.0 < s && s < 1.0 && 0.0 < t && t < 1.0 && fabs(s + t - 1.0) < epsnto2n))
    {
      //stype = _nto2n;

      int Num0 = 0;

      for(int i = 0; i < 4; ++i)
      {
        for(int j = i + 1;j < 4; ++j)
	{
           for(int k = j + 1; k < 4; ++k)
	   {
	   
	     avec12.setVec12(xp[i], yp[i], zp[i]);
	     bvec12.setVec12(xp[j], yp[j], zp[j]);
	     cvec12.setVec12(xp[k], yp[k], zp[k]);
	     dvec12.setVec12(dx, dy, dz);

	     double test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
	     if(test1 == 0.0)
	       ++Num0;

	   }
	}
      }

      cout <<  "2to6? " << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s - t -u) << endl;
      cout << "Num0" = Num0 << endl;
      if(Num0 == 1) 
      {
	//stype = _1to4;
	_2to6Flag = true;
      } else if(Num0 == 2)
      {
        stype = _nto2n;
      }
      //getchar();

      //cout <<  "2to6? " << s << " " << t << " " << u << " " << s + t + u << " " << fabs(1.0 - s -t -u) << endl;
    }

*/
#endif
  return stype;
}

int delaunaysplit3d1to4(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{

  //clock_t tmptime = clock();

  int rval = 0;
  int NewIndex[4];
  delaunayClass *delaunayNew[4];
  //delaunayClass aOldDelaunay = Delaunay[it];
  oldDelaunay[0].clone(Delaunay[it]);
  delaunayNew[0] = &Delaunay[it];
  delaunayNew[1] = &Delaunay[params.Ndelaunay];
  delaunayNew[2] = &Delaunay[params.Ndelaunay+1];
  delaunayNew[3] = &Delaunay[params.Ndelaunay+2];

  //cout << "start 1to4 " << it << endl;
  delaunayNew[0]->nNeighbor = 0;

  for(int ia = 0, i = 0; i < 4; ++i)
  {
    for(int j = i + 1; j < 4; ++j)
    {
      for(int k = j + 1; k < 4; ++k, ++ia)
      {
	int index = ia == 0 ? Delaunay[it].indd : params.Ndelaunay + ia - 1;
	NewIndex[ia] = index;
	delaunayNew[ia]->set(index, oldDelaunay[0].indm[i], oldDelaunay[0].indm[j], oldDelaunay[0].indm[k], imesh, Mesh);
	delaunayNew[ia]->circumcenter3d();
	delaunayNew[ia]->birthtype = "1to4";
	delaunayNew[ia]->sss = _sss;
	delaunayNew[ia]->ttt = _ttt;
	delaunayNew[ia]->uuu = _uuu;
	delaunayNew[ia]->vvv = _vvv;
      }
    }
  }
  params.Ndelaunay += 3;

  int interidx[3] = {2, 1, 0};
  for(int ia = 0, iface = 2; ia < 4; ++ia, --iface)
  {
     for(int inbr = 0, i = ia + 1; i < 4; ++i)
     {
       if(ia == i) continue;

       if(delaunayNew[ia]->nNeighbor >= 4)
	 continue;
       
       //faceClass *aFace = new faceClass();
       faceClass *aFace = &staticfaces[iiface]; ++iiface;

       delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = &Delaunay[NewIndex[i]];
       
       delaunayNew[ia]->faces[interidx[delaunayNew[ia]->nNeighbor]] = aFace;
       aFace->interDel[0] = delaunayNew[ia];
       aFace->interidx[0] = interidx[delaunayNew[ia]->nNeighbor];
       aFace->inbrs[0] = delaunayNew[ia]->nNeighbor;
      
       ++delaunayNew[ia]->nNeighbor;

       if(Delaunay[NewIndex[i]].nNeighbor >= 4) 
	 continue;
       Delaunay[NewIndex[i]].neighbors[Delaunay[NewIndex[i]].nNeighbor] = delaunayNew[ia];
       
       Delaunay[NewIndex[i]].faces[interidx[Delaunay[NewIndex[i]].nNeighbor]] = aFace;
       aFace->interDel[1] = &Delaunay[NewIndex[i]];
       aFace->interidx[1] = interidx[Delaunay[NewIndex[i]].nNeighbor]; 
       aFace->inbrs[1] = Delaunay[NewIndex[i]].nNeighbor;
       
       ++Delaunay[NewIndex[i]].nNeighbor;
       ++inbr;
     } 
     //delaunayNew[ia]->nNeighbor = 3;
  }

  //tmptime = clock() - tmptime;
  //time3 += tmptime;  //1to4

  //tmptime = clock();
  //change neighbor;
#if 0

  vector<int> indexs;
  vector<int> indexofnbr;
  vector<int> indexofa;
  vector<int> indexs2, indexofnbr2, indexofa2;
  vector<int> nbrVec;

  indexs.reserve(16);
  indexofnbr.reserve(16);
  indexofa.reserve(16);
  indexs2.reserve(16);
  indexofnbr2.reserve(16);
  indexofa2.reserve(16);
  nbrVec.reserve(16);
  
  for(int inbr = 0; inbr < oldDelaunay[0].nNeighbor; ++inbr)
  {
    nbrVec.push_back(oldDelaunay[0].neighbors[inbr]->indd);
  }

  for(int i = 0, j = 0; j < 4; ++j)
  {
    int jmesh = oldDelaunay[0].indm[j];
    for(int k = j + 1; k < 4; ++k)
    {
      int kmesh = oldDelaunay[0].indm[k];
      for(int l = k + 1; l < 4; ++l, ++i)
      {
	delaunayClass *newDel = delaunayNew[i];
	int lmesh = oldDelaunay[0].indm[l];

	for(int inbr = 0, nsize = oldDelaunay[0].nNeighbor; inbr < nsize; ++inbr)
	{
	  int index = -1;

	  //delaunayClass *tmpNeighbor = &Delaunay[oldDelaunay[0].neighbors[inbr]->indd];
	  delaunayClass *tmpNeighbor = &Delaunay[nbrVec[inbr]];
	  if(tmpNeighbor->has_mesh(jmesh) && tmpNeighbor->has_mesh(kmesh) && tmpNeighbor->has_mesh(lmesh))
	  {

	    for(int m = 0; m < tmpNeighbor->nNeighbor; ++m)
	    {
	      if(tmpNeighbor->neighbors[m]->indd == oldDelaunay[0].indd) 
	      {
		index = m; 
		break;
	      } 
	    }

	    //newDel->neighbors[ newDel->nNeighbor ] = tmpNeighbor; 
	    indexs2.push_back(tmpNeighbor->indd);
	    indexofnbr2.push_back(newDel->nNeighbor);
	    indexofa2.push_back(newDel->indd);
	    newDel->nNeighbor += 1;
	    //tmpNeighbor->neighbors[index] = &Delaunay[NewIndex[i]]; 
	    if(newDel->nNeighbor > 4) {cout << "warning" << endl;getchar();}
	    indexs.push_back(tmpNeighbor->indd);
	    indexofnbr.push_back(index);
	    indexofa.push_back(newDel->indd);
 
	    break;
	  }	
	}
      }
    }
  }

#ifndef _DEBUG

  for(int i = 0; i < indexs2.size(); ++i)
  {
    Delaunay[indexofa2[i]].neighbors[indexofnbr2[i]] = &Delaunay[indexs2[i]];
  }

   for(int i = 0; i < indexs.size(); ++i) 
  {
     Delaunay[indexs[i]].neighbors[indexofnbr[i]] = &Delaunay[indexofa[i]];
  }
#endif

#endif

#ifdef _DEBUG
   //cout << "1to4" << endl;

  int nNeighbor[4] = {3,3,3,3};

  for(int ia = 0, ia2 = 3; ia < 4; ++ia, --ia2)
  {
    int idl = 0;
    faceClass *face = oldDelaunay[0].faces[ia2]; 
    
    if(face == NULL) continue;

    delaunayClass *interDel = face->interDel[0];

    //対面となる四面体を取得する
    if(interDel->indd == oldDelaunay[0].indd) 
    {interDel = face->interDel[1]; idl = 1;}
     
    int inbr = face->inbrs[idl];

    interDel->neighbors[inbr] = delaunayNew[ia]; 

    delaunayNew[ia]->neighbors[3] = interDel;

    idl =  idl^1;

    face->interDel[idl] = delaunayNew[ia];
    face->interidx[idl] = 3;
    face->inbrs[idl] = nNeighbor[ia];// delaunayNew[ia]->nNeighbor;
    delaunayNew[ia]->faces[3] = face;
    //++delaunayNew[ia]->nNeighbor;
    ++nNeighbor[ia];
  }

  for(int i = 0; i < 4; ++i)
    delaunayNew[i]->nNeighbor = nNeighbor[i];

#endif
 
  //tmptime = clock() - tmptime;
  //time4 += tmptime;
  
  //tmptime = clock();

  //rval = oldDelaunay[0].eraseDelaunaryOfMesh(Mesh, params);
  //for(int i = 0; i < 4; ++i){ rval = delaunayNew[i]->addDelaunayToMesh(Mesh, params);}

  //tmptime = clock() - tmptime;
  //time3 += tmptime;
  
  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 4, params);

  //rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 4, params);
  for(int i = 0; i < 4; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
//    stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
    //stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
    tmpStackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
  }

  lastDel = NewIndex[3];

  return rval;

}

int delaunayflip3d(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params)
{
  int rval = 0;
  int no_indm[2];

  //clock_t tmptime = clock();
  flipType fliptype = Delaunay[it].getFlipType(imesh, jmesh, Mesh, true, no_indm);
  //tmptime = clock() - tmptime;
  //time8 += tmptime;

  switch(fliptype)
  {
    case(_2to3):
      rval = delaunayflip3d2to3(it, imesh, jmesh, Mesh, in, Delaunay, params);
      break;
    case(_3to2):
      rval = delaunayflip3d3to2(it, imesh, jmesh, Mesh, in, Delaunay, params);
      break; 
    case (_4to4):
      rval = delaunayflip3d4to4(it, imesh, jmesh, Mesh, in, Delaunay, params, no_indm);
      break;
    case(_noflip):
      break;
  }

  return rval;
}

flipType delaunayClass::getFlipType(int imesh, int jmesh, meshClass *Mesh, bool flag, int *no_indm)
{
  flipType fliptype = _noflip;
  meshClass *ptrMesh = &Mesh[jmesh];
  const double xpp = ptrMesh->getPosx(), ypp = ptrMesh->getPosy(), zpp = ptrMesh->getPosz();

  /*
  double n, m ,s, u;
  double a[3], b[3], c[3], a_inv[3], b_inv[3], c_inv[3];
  double det2;
  double xx, yy, zz;

  const double eps6 = 1.0e-4;
  
  a[0] = xp[0] - xp[2]; b[0] = xp[1] - xp[2]; c[0] = -xp[3] + xpp;
  a[1] = yp[0] - yp[2]; b[1] = yp[1] - yp[2]; c[1] = -yp[3] + ypp;
  a[2] = zp[0] - zp[2]; b[2] = zp[1] - zp[2]; c[2] = -zp[3] + zpp;
  det2 = DET3D(a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2]);

  a_inv[0] = (b[1]*c[2] - b[2]*c[1]); a_inv[1] = (b[0]*c[2]-b[2]*c[0]) * (-1.0) ; a_inv[2] = (b[0]*c[1] - b[1]*c[0]);
  b_inv[0] = -(a[1]*c[2] - a[2]*c[1]); b_inv[1] = (a[0]*c[2]-a[2]*c[0])         ; b_inv[2] = -(a[0]*c[1] - a[1]*c[0]);
  c_inv[0] = (a[1]*b[2] - a[2]*b[1]); c_inv[1] = (a[0]*b[2]-a[2]*b[0]) * (-1.0) ; c_inv[2] = (a[0]*b[1] - a[1]*b[0]);

  n = ((xpp - xp[2]) * a_inv[0] + (ypp - yp[2]) * a_inv[1] + (zpp - zp[2]) * a_inv[2]) / det2;  
  m = ((xpp - xp[2]) * b_inv[0] + (ypp - yp[2]) * b_inv[1] + (zpp - zp[2]) * b_inv[2]) / det2;  
  s = ((xpp - xp[2]) * c_inv[0] + (ypp - yp[2]) * c_inv[1] + (zpp - zp[2]) * c_inv[2]) / det2;
  u = 1.0 - n - m;
  xx = xpp + (s) * (-xpp + xp[3]);
  yy = ypp + (s) * (-ypp + yp[3]);
  zz = zpp + (s) * (-zpp + zp[3]);

  _sss = n; _ttt = m; _uuu = 1.0 - n - m; _vvv = s;
  if(!(fabs(n) <= eps6 || fabs(m) <= eps6 || fabs(1.0 - n -m) <= eps6))
  {
    if(0.0 < n && n < 1.0 && 0.0 < m && m < 1.0 && 0.0 < 1.0 - n - m && 1.0 - n - m < 1.0 && 0.0 < s && s <= 1.0 + 1.0e-14)
    { 
      return fliptype =  _2to3;
    } else
    {
      //if(this->skipflipFlag(xx, yy, zz))
      {
//	fliptype = _noflip;
	//return fliptype;
      }
      //else
      {
	return fliptype = _3to2;
      }

    }

    //return fliptype;
  }
*/

  int i = 0, j = 0, k = 0;
  double test1, test2, ans = 0;

  //case 1 3 0 1 2
  i = 0; j = 1; k = 2;
  avec12.setVec12(xp[3], yp[3], zp[3]);
  bvec12.setVec12(xp[i], yp[i], zp[i]);
  cvec12.setVec12(xp[j], yp[j], zp[j]);
  dvec12.setVec12(xp[k], yp[k], zp[k]);

  //test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
  test1 = -this->orient;

  if(test1 == 0.0) goto CLEANUP;

  dvec12.setVec12(xpp, ypp, zpp);
  test2 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

  if(sign(test1) != sign(test2) && test2 != 0)
  {
    ans = 2;
    goto CLEANUP;
  }
  if(test2 == 0)
  {
    ans = 3;
    no_indm[0] = this->indm[i];
    no_indm[1] = this->indm[j];
    goto CLEANUP;
  }
  
  //case 2 3 0 1 2
  i = 0; j = 2; k = 1;
  avec12.setVec12(xp[3], yp[3], zp[3]);
  bvec12.setVec12(xp[i], yp[i], zp[i]);
  cvec12.setVec12(xp[j], yp[j], zp[j]);
  dvec12.setVec12(xp[k], yp[k], zp[k]);
  
  //test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
  test1 = this->orient;
  if(test1 == 0.0) goto CLEANUP;

  dvec12.setVec12(xpp, ypp, zpp);
  //dvec12 = Vec::getVec12(xpp, ypp, zpp);
  test2 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

  if(sign(test1) != sign(test2) && test2 != 0)
  {
    ans = 2;
    goto CLEANUP;

  }
  if(test2 == 0)
  {
    ans = 3;
    no_indm[0] = this->indm[i];
    no_indm[1] = this->indm[j];
    goto CLEANUP;
  }
  
  //case 3 3 1 2 0 
  i = 1; j = 2; k = 0;

  avec12.setVec12(xp[3], yp[3], zp[3]);
  bvec12.setVec12(xp[i], yp[i], zp[i]);
  cvec12.setVec12(xp[j], yp[j], zp[j]);
  dvec12.setVec12(xp[k], yp[k], zp[k]);
  //test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);
  test1 = -this->orient;
  if(test1 == 0.0) goto CLEANUP;

  dvec12.setVec12(xpp, ypp, zpp);
  test2 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

  if(sign(test1) != sign(test2) && test2 != 0)
  {
    ans = 2;
    goto CLEANUP;
  }
  if(test2 == 0)
  {
    ans = 3;
    no_indm[0] = this->indm[i];
    no_indm[1] = this->indm[j];
    goto CLEANUP;
  }
CLEANUP:

  if(test1 == 0.0)
  {
    fliptype = _2to3;
    _2to6Flag = false;
    return fliptype;
  }
  fliptype = _2to3;
  if(ans == 2)
    fliptype = _3to2;
  else if (ans == 3)
    fliptype = _4to4;
  
  return fliptype;
}

bool delaunayClass::skipflipFlag(double xp, double yp, double zp)
{
  double a[3] = {this->getPosx(1) - this->getPosx(0), this->getPosy(1) - this->getPosy(0), this->getPosz(1) - this->getPosz(0)};
  double b[3] = {this->getPosx(2) - this->getPosx(0), this->getPosy(2) - this->getPosy(0), this->getPosz(2) -this->getPosz(0)};
  double nvec[3] = {a[1]*b[2] - a[2]*b[1],
  		    a[2]*b[0] - a[0]*b[2],
  		    a[0]*b[1] - a[1]*b[0]};
  double leng = sqrt(sq(nvec[0]) +sq(nvec[1]) + sq(nvec[2]));

  nvec[0] = nvec[0] / leng;
  nvec[1] = nvec[1] / leng;
  nvec[2] = nvec[2] / leng;
  double A[3] = {this->getPosx(0), this->getPosy(0), this->getPosz(0)};
  double B[3] = {this->getPosx(1), this->getPosy(1), this->getPosz(1)};
  double C[3] = {this->getPosx(2), this->getPosy(2), this->getPosz(2)};
  double D[3] = {xp, yp, zp};
  double A2[3], B2[3], C2[3], D2[3];
  
  double norm = sqrt(nvec[0]*nvec[0] +nvec[1]*nvec[1]);
  double nvec2[3]={nvec[0], nvec[1], nvec[2]};
  double cos,sin;

  if(norm != 0.0)
  {
    cos = nvec[0] / norm; sin = nvec[1] / norm;     
    A2[0] = cos * A[0] + sin * A[1];
    A2[1] = -sin * A[0] + cos * A[1];
    A2[2] = A[2];

    B2[0] = cos * B[0] + sin * B[1];
    B2[1] = -sin * B[0] + cos * B[1];
    B2[2] = B[2];
 
    C2[0] = cos * C[0] + sin * C[1];
    C2[1] = -sin * C[0] + cos * C[1];
    C2[2] = C[2];
 
    D2[0] = cos * D[0] + sin * D[1];
    D2[1] = -sin * D[0] + cos * D[1];
    D2[2] = D[2];
 
    nvec2[0]  = cos * nvec[0] + sin * nvec[1];
    nvec2[1] = -sin * nvec[0]+ sin * nvec[1];
    nvec2[2] = nvec[2];
  } else
  {
    for(int i = 0; i < 3; ++i){
      A2[i] = A[i];
      B2[i] = B[i];
      C2[i] = C[i];
      D2[i] = D[i];
    }
  }

  norm = sqrt(nvec2[0]*nvec2[0] + nvec2[2]*nvec2[2]);
  if(norm != 0.0){
    cos = nvec2[0]/norm; sin = nvec2[2]/norm;
    A[0] = cos * A2[0] + sin * A2[2];
    A[1] = A2[1];
    A[2] = -sin * A2[0] + cos * A2[2];

    B[0] = cos * B2[0] + sin * B2[2];
    B[1] = B2[1];
    B[2] = -sin * B2[0] + cos * B2[2];
 
    C[0] = cos * C2[0] + sin * C2[2];
    C[1] = C2[1];
    C[2] = -sin * C2[0] + cos * C2[2];
 
    D[0] = cos * D2[0] + sin * D2[2];
    D[1] = D2[1];
    D[2] = -sin * D2[0] + cos * D2[2];
 
  }
  else{
    for(int i = 0; i < 3; ++i)
    {
      A[i] = A2[i];
      B[i] = B2[i];
      C[i] = C2[i];
      D[i] = D2[i];
    }
  }

  double s ,u;
  double p[2] = {D[1] - A[1], D[2] - A[2]};
  double X[2] = {B[1] - A[1], C[1] - A[1]}; 
  double Y[2] = {B[2] - A[2], C[2] - A[2]};
  double det = X[0]*Y[1] - X[1]*Y[0];
  const double eps2 = 0.0;

  if(det != 0.0){
    s = (Y[1]*p[0] - X[1]*p[1]) / det;
    u = (-Y[0]*p[0] + X[0]*p[1]) / det;
    if(s > eps2 && u > eps2)
      return false;
  }

  p[0] = D[1] - B[1]; p[1] = D[2] - B[2];
  X[0] = A[1] - B[1]; X[1] = C[1] - B[1]; 
  Y[0] = A[2] - B[2]; Y[1] = C[2] - B[2];
  det = X[0]*Y[1] - X[1]*Y[0];

  if(det != 0.0){
    s = (Y[1]*p[0] - X[1]*p[1]) / det;
    u = (-Y[0]*p[0] + X[0]*p[1]) / det;

    if(s > eps2 && u > eps2)
      return false;
  }

  p[0] = D[1] - C[1]; p[1] = D[2] - C[2];
  X[0] = A[1] - C[1]; X[1] = B[1] - C[1]; 
  Y[0] = A[2] - C[2]; Y[1] = B[2] - C[2];
  det = X[0]*Y[1] - X[1]*Y[0];
    
  if(det != 0.0){
    s = (Y[1]*p[0] - X[1]*p[1]) / det;
    u = (-Y[0]*p[0] + X[0]*p[1]) / det;

    if(s > eps2 && u > eps2)
      return false;
  }

  return true;

}

int delaunayflip3d2to3(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params)
{
  //clock_t tmptime = clock();
  int rval = 0, ntop[3];
  //cout << "2to3" << endl;

  int NewIndex[3] = {Delaunay[it].indd, Delaunay[in].indd, params.Ndelaunay};
  //delaunayClass oldDelaunay[2];
  delaunayClass *delaunayNew[3];
  //oldDelaunay[0] = Delaunay[it];
  //oldDelaunay[1] = Delaunay[in]; 
  oldDelaunay[0].clone(Delaunay[it]);
  oldDelaunay[1].clone(Delaunay[in]);
  delaunayNew[0] = &Delaunay[it];
  delaunayNew[1] = &Delaunay[in];
  delaunayNew[2] = &Delaunay[params.Ndelaunay];

  //cout << "2to3 " << Delaunay[it].getRadius() << endl;
  for(int i = 0, j = 0; i < 4; ++i)
  {
    if(oldDelaunay[0].indm[i] != imesh)
    {
      ntop[j] = oldDelaunay[0].indm[i]; ++j;
    } 
  }

  for(int ia = 0, i = 0; i < 3; ++i)
  {
    for(int j = i + 1; j < 3; ++j, ++ia)
    {
      delaunayNew[ia]->set(NewIndex[ia], ntop[i], ntop[j], jmesh, imesh,  Mesh);
      delaunayNew[ia]->circumcenter3d();
      delaunayNew[ia]->birthtype = "2to3";
    }
  }

  params.Ndelaunay += 1;

  //0. ntop0 ntop1 jmesh imesh (1, 0)
  //1. ntop0 ntop2 jmesh imesh (1, 0)
  //2. ntop1 ntop2 jmesh imesh (1, 0)
  int interidx[2] = {1, 0};
  
  for(int ia = 0; ia < 3; ++ia)
  {
     for(int inbr = 0, i = ia + 1; i < 3; ++i)
     {
       if(ia == i) continue;

       if(delaunayNew[ia]->nNeighbor >= 4)
	 continue;

       if(Delaunay[NewIndex[i]].nNeighbor >= 4)
	 continue;

       //faceClass *face = new faceClass();
       faceClass *face = &staticfaces[iiface]; ++iiface;

       delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = &Delaunay[NewIndex[i]];
       
       delaunayNew[ia]->faces[interidx[delaunayNew[ia]->nNeighbor]] = face;
       face->interDel[0] = delaunayNew[ia];
       face->interidx[0] = interidx[delaunayNew[ia]->nNeighbor];
       face->inbrs[0] = delaunayNew[ia]->nNeighbor;
       
       ++delaunayNew[ia]->nNeighbor;
       
       Delaunay[NewIndex[i]].neighbors[Delaunay[NewIndex[i]].nNeighbor] = delaunayNew[ia];
       
       Delaunay[NewIndex[i]].faces[interidx[Delaunay[NewIndex[i]].nNeighbor]] = face;
       face->interDel[1] = &Delaunay[NewIndex[i]];
       face->interidx[1] = interidx[Delaunay[NewIndex[i]].nNeighbor];
       face->inbrs[1] = Delaunay[NewIndex[i]].nNeighbor;
       
       ++Delaunay[NewIndex[i]].nNeighbor;

       ++inbr;
     } 
     delaunayNew[ia]->nNeighbor = 2;
  }
  
  //tmptime = clock() - tmptime;
  //time3 += tmptime; //2to3
  //tmptime = clock();

#if 0

  vector<int> indexs, indexs2;
  vector<int> indexsofnbr, indexofnbr2;
  vector<int> indexofa, indexofa2;
  vector<vector<int> > oldnbrVec;

  indexs.reserve(16); indexs2.reserve(16);
  indexsofnbr.reserve(16); indexofnbr2.reserve(16);
  indexofa.reserve(16); indexofa2.reserve(16);
  oldnbrVec.reserve(16);


  for(int ia = 0; ia < 2; ++ia)
  {
    vector<int> tmpVec; tmpVec.reserve(10);
    delaunayClass *tmpOld = &oldDelaunay[ia];
    for(int inbr = 0; inbr < tmpOld->nNeighbor; ++inbr)
    {
      tmpVec.push_back(tmpOld->neighbors[inbr]->indd);
    }
    oldnbrVec.push_back(tmpVec);
  }
  for(int ia = 0; ia < 2; ++ia)
  {
    delaunayClass *tmpOld = &oldDelaunay[ia];
    
    for(int inbr = 0; inbr < tmpOld->nNeighbor; ++inbr)
    {
      //delaunayClass *tmpNbr = tmpOld->neighbors[inbr]; 
      delaunayClass *tmpNbr = &Delaunay[oldnbrVec[ia][inbr]];
      if(tmpNbr->indd == oldDelaunay[0].indd || tmpNbr->indd == oldDelaunay[1].indd) continue;
    
      int index = 0;

      for(int l = 0; l < tmpNbr->nNeighbor; ++l)
	if(tmpNbr->neighbors[l]->indd == tmpOld->indd){ index = l; break; } 

      int ia2 = -1; 
      //if(tmpNbr->has_mesh(ntop[0]) && tmpNbr->has_mesh(ntop[1]) && delaunayNew[0]->neighbors[delaunayNew[0]->nNeighbor-1] != tmpNbr) ia2 = 0;
      //else if (tmpNbr->has_mesh(ntop[0]) && tmpNbr->has_mesh(ntop[2]) && delaunayNew[1]->neighbors[delaunayNew[1]->nNeighbor-1] != tmpNbr) ia2 = 1;
      //else if (tmpNbr->has_mesh(ntop[1]) && tmpNbr->has_mesh(ntop[2]) && delaunayNew[2]->neighbors[delaunayNew[2]->nNeighbor - 1] != tmpNbr) ia2 = 2; 
      if(tmpNbr->has_mesh(ntop[0]) && tmpNbr->has_mesh(ntop[1])) ia2 = 0;
      else if (tmpNbr->has_mesh(ntop[0]) && tmpNbr->has_mesh(ntop[2])) ia2 = 1;
      else if (tmpNbr->has_mesh(ntop[1]) && tmpNbr->has_mesh(ntop[2])) ia2 = 2; 
     
      if(ia2 != -1) 
      {
        //tmpNbr->neighbors[index] = delaunayNew[ia2];	
	indexs.push_back(tmpNbr->indd);
	indexsofnbr.push_back(index);
	indexofa.push_back(ia2);

	//delaunayNew[ia2]->neighbors[delaunayNew[ia2]->nNeighbor] = tmpNbr;
	indexofa2.push_back(delaunayNew[ia2]->indd);
	indexofnbr2.push_back(delaunayNew[ia2]->nNeighbor);
	indexs2.push_back(tmpNbr->indd);
	delaunayNew[ia2]->nNeighbor += 1;
	//break;

      }
    }
  }
#ifndef _DEBUG
   for(int i = 0; i < indexs2.size(); ++i)
  {
    Delaunay[indexofa2[i]].neighbors[indexofnbr2[i]] = &Delaunay[indexs2[i]];
  }

  for(int i = 0; i < indexs.size(); ++i) 
  {
     Delaunay[indexs[i]].neighbors[indexsofnbr[i]] = delaunayNew[indexofa[i]];
  }
#endif

#endif

  int ntopIdx[2][3];
  for(int ia = 0; ia < 2; ++ia)
  {
    int ii = 0;
    for(int j = 0; j < 4; ++j)
    {
      if(ntop[0] == oldDelaunay[ia].indm[j])
      {
        ntopIdx[ia][0] = j;
	++ii;
      } else if(ntop[1] == oldDelaunay[ia].indm[j])
      {
        ntopIdx[ia][1] = j;
	++ii;
      } else if(ntop[2] == oldDelaunay[ia].indm[j])
      {
        ntopIdx[ia][2] = j;
	++ii;
      }

      if(ii == 3) break;
    }
  }

 #ifdef _DEBUG
  //cout << "start 2to3 " << endl;

  int nNeighbor[3] = {2,2,2};

  for(int ia = 0, ii = 2; ia < 3; ++ia, --ii)
  {
    for(int j = 2; j < 4; ++j)
    {
      int idl = 0;
      faceClass *face = oldDelaunay[j-2].faces[ntopIdx[j-2][ii]]; 

      if(face == NULL) continue;

      delaunayClass *interDel = face->interDel[0];
      if(interDel->indd == oldDelaunay[j - 2].indd)
      {
        interDel = face->interDel[1]; idl = 1; 
      }

      int inbr = face->inbrs[idl];
      interDel->neighbors[inbr] = delaunayNew[ia];

      idl = idl^1;
      
      delaunayNew[ia]->neighbors[j] = interDel;
      delaunayNew[ia]->faces[j] = face;
      
      face->interidx[idl] = j;
      //face->inbrs[idl] = delaunayNew[ia]->nNeighbor;
      face->inbrs[idl] = nNeighbor[ia];
      face->interDel[idl] = delaunayNew[ia];
      nNeighbor[ia]++;
      //++delaunayNew[ia]->nNeighbor;
    
    }
  }
  //cout << "end" << endl;

  for(int i = 0; i < 3; ++i)
    delaunayNew[i]->nNeighbor = nNeighbor[i];
  
#endif
   
  //tmptime = clock() - tmptime;
  //time4 += tmptime;
  //tmptime = clock();
  /*
  for(int i = 0; i < 2; ++i)
  {
    oldDelaunay[i].eraseDelaunaryOfMesh(Mesh, params);
  }
  for(int i = 0; i < 3; ++i)
  {
    delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }
  */

  //tmptime = clock() - tmptime;
  //time5 += tmptime;

  //cout << "end 2to3 " << imesh << endl;

  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 3, params);

  //rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 3, params);
  for(int i = 0; i < 3; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
//    stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
    //stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
    tmpStackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);

  }

  lastDel = NewIndex[2];

  return rval;
}

int delaunayflip3d3to2(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params)
{
  //clock_t tmptime = clock();

  int rval = 0;
  delaunayClass *delaunayNew[2];
  int NewIndex[2] = {it, in}, no_indm[2];
  int sharedIndexes[3][2];
  delaunayNew[0] = &Delaunay[it];
  delaunayNew[1] = &Delaunay[in];

  //cout << "3to2" << endl;
  bool getnbrFlag = false;
  for(int inbr = 0; inbr < Delaunay[it].nNeighbor; ++inbr)
  {
    delaunayClass *pnbr = Delaunay[it].neighbors[inbr];
    if(pnbr->has_mesh(imesh) && pnbr->has_mesh(jmesh)) 
    {
       //oldDelaunay[2] = Delaunay[pnbr->indd]; 
      oldDelaunay[2].clone(Delaunay[pnbr->indd]);
      getnbrFlag = true;
    }
  }
  if(!getnbrFlag){rval = 2;/*tmptime = clock() - tmptime; time6 += tmptime*/; return rval;}

  //oldDelaunay[0] = Delaunay[it];
  //oldDelaunay[1] = Delaunay[in];
  oldDelaunay[0].clone(Delaunay[it]);
  oldDelaunay[1].clone(Delaunay[in]);

  for(int i = 1, k = 0; i < 3; ++i)
  {
    if(oldDelaunay[i].has_mesh(imesh))
    {
      for(int j = 0; j < 4; ++j)
      {
        if(oldDelaunay[i].indm[j] != imesh && oldDelaunay[i].indm[j] != jmesh)
	{
           no_indm[k] = oldDelaunay[i].indm[j];	
	   ++k;
	} 
      }	
    } 
  }

  delaunayClass *pdel = &oldDelaunay[0];
  for(int i = 0, ia = 0; i < 4; ++i)
  {
    for(int j = i + 1; j < 4; ++j)
    {
      if(pdel->indm[i] != imesh && pdel->indm[j] != imesh)
      {
	if((pdel->indm[i] != no_indm[0] || pdel->indm[j] != no_indm[1]) &&
	    (pdel->indm[i] != no_indm[1] || pdel->indm[j] != no_indm[0]))
	{
	  if(pdel->indm[i] == no_indm[0] || pdel->indm[j] == no_indm[0]) ia = 0;
	  else if(pdel->indm[i] == no_indm[1] || pdel->indm[j] == no_indm[1]) ia = 1;

	  if(pdel->indm[i] == no_indm[0] || pdel->indm[i] == no_indm[1])
	    delaunayNew[ia]->set(NewIndex[ia], pdel->indm[i], pdel->indm[j], jmesh, imesh, Mesh); 
	  else 
	    delaunayNew[ia]->set(NewIndex[ia], pdel->indm[j], pdel->indm[i], jmesh, imesh, Mesh);
	  delaunayNew[ia]->birthtype = "3to2";
	  delaunayNew[ia]->circumcenter3d();
	} 
      }
    } 
  }

  //get shard edge indexes
  for(int ia = 0; ia < 3; ++ia)
  {
    int ii = 0;
    for(int j = 0; j < 4; ++j)
    {
      if(no_indm[0] == oldDelaunay[ia].indm[j] || no_indm[1] == oldDelaunay[ia].indm[j])
      {
	if(no_indm[0] == oldDelaunay[ia].indm[j])
	  sharedIndexes[ia][0] = j; 
	else if(no_indm[1] == oldDelaunay[ia].indm[j])
	  sharedIndexes[ia][1] = j;
	++ii;
	if(ii == 2) 
	  break;
      }
    }
  }

  //0. no_indm[0] ? jmesh imesh
  //1. no_indm[1] ? jmesh imesh
  delaunayNew[0]->neighbors[0] = &Delaunay[NewIndex[1]];
  delaunayNew[1]->neighbors[0] = &Delaunay[NewIndex[0]];

  //faceClass *face = new faceClass();
  faceClass *face = &staticfaces[iiface]; ++iiface;
  delaunayNew[0]->faces[0] = face;
  face->interDel[0] = delaunayNew[0];
  face->interidx[0] = 0;
  face->inbrs[0] = 0;

  delaunayNew[1]->faces[0] = face;
  face->interDel[1] = delaunayNew[1];
  face->interidx[1] = 0;
  face->inbrs[1] = 0;

  delaunayNew[0]->nNeighbor = 1;
  delaunayNew[1]->nNeighbor = 1;

  Delaunay[oldDelaunay[2].indd].flag = false;

#if 0
  vector<int> indexs, indexs2;
  vector<int> indexofnbr, indexofnbr2;
  vector<int> indexofa, indexofa2;

  vector<vector<int> > oldnbrVec;

  indexs.reserve(16); indexs2.reserve(16);
  indexofnbr.reserve(16); indexofnbr2.reserve(16);
  indexofa.reserve(16); indexofa2.reserve(16);
  oldnbrVec.reserve(16);

  //tmptime = clock() - tmptime;
  //time3 += tmptime; //3to2
  //tmptime = clock();
  for(int ia = 0; ia < 3; ++ia)
  {
    vector<int> tmpVec;
    delaunayClass *tmpOld = &oldDelaunay[ia];
    for(int inbr = 0; inbr < tmpOld->nNeighbor; ++inbr)
    {
      tmpVec.push_back(tmpOld->neighbors[inbr]->indd);
    }
    oldnbrVec.push_back(tmpVec);
  }

  for(int ia = 0; ia < 3; ++ia)
  {
    pdel = &oldDelaunay[ia]; 

    for(int inbr = 0; inbr < pdel->nNeighbor; ++inbr)
    {
      //delaunayClass *tmpNbr = pdel->neighbors[inbr];
      delaunayClass *tmpNbr = &Delaunay[oldnbrVec[ia][inbr]];
      if(tmpNbr->indd == oldDelaunay[0].indd || tmpNbr->indd == oldDelaunay[1].indd|| tmpNbr->indd == oldDelaunay[2].indd) continue; 

      int index = 0;
      for(int l = 0; l < tmpNbr->nNeighbor; ++l)
      {
        if(tmpNbr->neighbors[l]->indd == pdel->indd) {index = l;break;} 
      }

      delaunayClass *newDel = NULL;
      
      //  if(tmpNbr->has_mesh(no_indm[0]) && delaunayNew[0]->neighbors[delaunayNew[0]->nNeighbor-1] != tmpNbr) newDel = delaunayNew[0];
      //else if(tmpNbr->has_mesh(no_indm[1]) && delaunayNew[1]->neighbors[delaunayNew[1]->nNeighbor-1] != tmpNbr) newDel = delaunayNew[1]; 

      if(tmpNbr->has_mesh(no_indm[0]) ) newDel = delaunayNew[0];
      else if(tmpNbr->has_mesh(no_indm[1]) ) newDel = delaunayNew[1]; 

      //if(tmpNbr->has_mesh(no_indm[0])) newDel = delaunayNew[0];
      //else if(tmpNbr->has_mesh(no_indm[1])) newDel = delaunayNew[1];
      //else {cout << "asdlfjkasd;f" << endl;getchar();}
      
      //newDel->neighbors[ newDel->nNeighbor ] = tmpNbr; 

      indexofa2.push_back(newDel->indd);
      indexofnbr2.push_back(newDel->nNeighbor);
      indexs2.push_back(tmpNbr->indd);

      newDel->nNeighbor += 1;
      if(newDel->nNeighbor >= 5) 
      {cout << "waring" << endl;getchar();}
      //tmpNbr->neighbors[index] = &Delaunay[newDel->indd]; 

      indexs.push_back(tmpNbr->indd);
      indexofnbr.push_back(index);
      indexofa.push_back(newDel->indd);
 
    } 
    
  }

#ifndef _DEBUG
  for(int i = 0; i < indexs2.size(); ++i)
  {
    Delaunay[indexofa2[i]].neighbors[indexofnbr2[i]] = &Delaunay[indexs2[i]]; 
  }

  for(int i = 0; i < indexs.size(); ++i) 
  {
    Delaunay[indexs[i]].neighbors[indexofnbr[i]] = &Delaunay[indexofa[i]];
  }
#endif

#endif
  
#ifdef _DEBUG 
  //cout << "3to2 " << endl;
  int nNeighbor[2] = {1, 1};

  for(int ianew = 0, idx = 1; ianew < 2; ++ianew, --idx)
  {
    for(int ia = 0; ia < 3; ++ia)
    {
      int idl = 0;
      faceClass *face = oldDelaunay[ia].faces[sharedIndexes[ia][idx]];

      if(face == NULL) continue;
      delaunayClass *interDel = face->interDel[0];
      
      if(interDel->indd == oldDelaunay[ia].indd)
      {interDel = face->interDel[1]; idl = 1;} 

      int inbr = face->inbrs[idl];
      interDel->neighbors[inbr] = delaunayNew[ianew];

      idl = idl^1; 

      face->interDel[idl] = delaunayNew[ianew];
      if(oldDelaunay[ia].has_mesh(imesh) && oldDelaunay[ia].has_mesh(jmesh))
      {
        //delaunayNew[ianew]->neighbors[delaunayNew[ianew]->nNeighbor] = interDel;
	delaunayNew[ianew]->neighbors[nNeighbor[ianew]] = interDel;
	delaunayNew[ianew]->faces[1] = face;
	face->interidx[idl] = 1;
	face->inbrs[idl] = nNeighbor[ianew];
	//face->inbrs[idl] = delaunayNew[ianew]->nNeighbor;
	//++delaunayNew[ianew]->nNeighbor;
	++nNeighbor[ianew];
      }
      else if (oldDelaunay[ia].has_mesh(imesh) && !oldDelaunay[ia].has_mesh(jmesh))
      {
        //delaunayNew[ianew]->neighbors[delaunayNew[ianew]->nNeighbor] = interDel;
	delaunayNew[ianew]->neighbors[nNeighbor[ianew]] = interDel;
	delaunayNew[ianew]->faces[2] = face;
	face->interidx[idl] = 2;
	face->inbrs[idl] = nNeighbor[ianew];
	//face->inbrs[idl] = delaunayNew[ianew]->nNeighbor;
	//++delaunayNew[ianew]->nNeighbor;
	++nNeighbor[ianew];
      }
      else if (!oldDelaunay[ia].has_mesh(imesh) && oldDelaunay[ia].has_mesh(jmesh))
      {
         //delaunayNew[ianew]->neighbors[delaunayNew[ianew]->nNeighbor] = interDel;
	 delaunayNew[ianew]->neighbors[nNeighbor[ianew]] = interDel;
	
	 delaunayNew[ianew]->faces[3] = face;
	 face->interidx[idl] = 3;
	 face->inbrs[idl] = nNeighbor[ianew];
	 //face->inbrs[idl] = delaunayNew[ianew]->nNeighbor;
	 //++delaunayNew[ianew]->nNeighbor;
	 ++nNeighbor[ianew];
      } 
    }
  }

  for(int i = 0; i < 2; ++i)
    delaunayNew[i]->nNeighbor = nNeighbor[i];

  //cout << "end 3to2" << endl;
   
#endif
    
  //tmptime = clock() - tmptime;
  //time4 += tmptime;

  //tmptime = clock();

  /*
  for(int i = 0; i < 3; ++i)
  {
    oldDelaunay[i].eraseDelaunaryOfMesh(Mesh, params);
  }
  for(int i = 0; i < 2; ++i)
  {
    delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }
  */

//  tmptime = clock() - tmptime;
//  time6 += tmptime;
//cout << "end 3to2 " << endl;
  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 2, params);
  //rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 2, params);
  for(int i = 0; i < 2; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
//    stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
    //stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
    tmpStackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
  }

    lastDel = NewIndex[1];

  return rval;
}

int delaunayflip3d4to4(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params, int no_indm[2])
{
  // need check coplane ? 
  //clock_t tmptime = clock();
  int rval = 0;
  int sharedEdgeIdx[4][2];
  delaunayClass *delaunayNew[4];
  delaunayClass *tmpOld0;

  oldDelaunay[0].clone(Delaunay[it]); oldDelaunay[2].clone(Delaunay[in]);
  //oldDelaunay[0] = Delaunay[it];
  //oldDelaunay[2] = Delaunay[in];
  tmpOld0 = &oldDelaunay[0];

  //cout << "start 4to4 " <<  it << endl;

  bool flag = false;
  delaunayClass *Nbr;

  for(int inbr = 0; inbr < tmpOld0->nNeighbor; ++inbr)
  {
    Nbr = tmpOld0->neighbors[inbr];
    if(!Nbr->has_mesh(jmesh) && Nbr->has_mesh(imesh))
    {
      if(Nbr->has_mesh(no_indm[0]) && Nbr->has_mesh(no_indm[1]))
      {
	oldDelaunay[1].clone(Delaunay[Nbr->indd]);
	//oldDelaunay[1] = Delaunay[Nbr->indd];
	flag = true;	   
	break;
      }
    } 
  }
  //tmptime = clock() - tmptime;

  if(!flag){ rval = 3;/* tmptime = clock() - tmptime; time7 += tmptime;*/ return rval;}
  // get oldDelaunay[3]
  flag = false;
  delaunayClass *tmpOld2 = &oldDelaunay[2];
  for(int inbr = 0; inbr < tmpOld2->nNeighbor; ++inbr)
  {
    Nbr = tmpOld2->neighbors[inbr]; 
    
    if(!Nbr->has_mesh(imesh) && Nbr->has_mesh(jmesh) && Nbr->has_mesh(no_indm[0]) && Nbr->has_mesh(no_indm[1]))
    {
      oldDelaunay[3].clone(Delaunay[Nbr->indd]);
      //oldDelaunay[3] = Delaunay[Nbr->indd];
      flag = true;
      break;
    }
  }
  
  if(!flag) { rval = 3;/*tmptime = clock() - tmptime; time7 += tmptime;*/ return rval;}

  flag = false;
  delaunayClass *tmpOld1 = &oldDelaunay[1];
  delaunayClass *tmpOld3 = &oldDelaunay[3];
  for(int inbr = 0; inbr < tmpOld1->nNeighbor; ++inbr)
  {
     if(tmpOld3->indd == tmpOld1->neighbors[inbr]->indd) 
     {flag = true; break;} 
  }

  if(!flag){ rval = 3;/* tmptime = clock() - tmptime; time7 += tmptime;*/  return rval;}

  int NewIndex[4] = {oldDelaunay[0].indd,oldDelaunay[1].indd, oldDelaunay[2].indd, oldDelaunay[3].indd};
  
  //time3 += tmptime; //4to4
  //tmptime = clock();
  //cout << "generate 4to4 " << endl;
  // create new delaunay;
  for(int i = 0, ic = 0; i < 2; ++i)
  {
    delaunayClass *tmpOld = &oldDelaunay[i];
    for(int j = 0; j < 3; ++j)
    {
      for(int k = j + 1; k < 3; ++k)
      {
      
	if( (tmpOld->indm[j] != no_indm[0] || tmpOld->indm[k] != no_indm[1]) && 
	    (tmpOld->indm[j] != no_indm[1] || tmpOld->indm[k] != no_indm[0]) )
	{
	  delaunayNew[ic] = &Delaunay[NewIndex[ic]];
	  if(tmpOld->indm[j] == no_indm[0] || tmpOld->indm[j] == no_indm[1])
	    delaunayNew[ic]->set(NewIndex[ic], tmpOld->indm[j], tmpOld->indm[k], jmesh, imesh, Mesh);
	  else 
	    delaunayNew[ic]->set(NewIndex[ic], tmpOld->indm[k], tmpOld->indm[j], jmesh, imesh, Mesh);
	  delaunayNew[ic]->circumcenter3d();
	  delaunayNew[ic]->birthtype = "4to4";
	  ++ic;
	}
      
      }
    }
  } 

  for(int ia = 0, ii = 0; ia < 4; ++ia)
  {
    ii = 0;
    for(int j = 0; j < 4; ++j)
    {
      if(no_indm[0] == oldDelaunay[ia].indm[j] || no_indm[1] == oldDelaunay[ia].indm[j])
      {

	if(no_indm[0] == oldDelaunay[ia].indm[j])
	  sharedEdgeIdx[ia][0] = j;
	else if(no_indm[1] == oldDelaunay[ia].indm[j])
	  sharedEdgeIdx[ia][1] = j;
	 ++ii;

	 if(ii == 2) 
	   break;
      }	
    }
  }

  // otonarisann
  delaunayNew[0]->neighbors[0] = delaunayNew[1]; 
  delaunayNew[1]->neighbors[0] = delaunayNew[0];
  delaunayNew[2]->neighbors[0] = delaunayNew[3];
  delaunayNew[3]->neighbors[0] = delaunayNew[2];

  for(int ia1 = 0, ia2 = 1; ia1 < 4; ia1 += 2, ia2 += 2)
  {
    //faceClass *face = new faceClass();
    faceClass *face = &staticfaces[iiface]; ++iiface;

    delaunayNew[ia1]->faces[0] = face;
    face->interDel[0] = delaunayNew[ia1];
    face->interidx[0] = 0;
    face->inbrs[0] = 0;
    
    delaunayNew[ia2]->faces[0] = face;
    face->interDel[1] = delaunayNew[ia2];
    face->interidx[1] = 0;
    face->inbrs[1] = 0;
  }


  if(delaunayNew[0]->indm[0] == delaunayNew[2]->indm[0] || delaunayNew[0]->indm[1] == delaunayNew[2]->indm[1] ||
     delaunayNew[0]->indm[0] == delaunayNew[2]->indm[1] || delaunayNew[0]->indm[1] == delaunayNew[2]->indm[0])
  {
  
    delaunayNew[0]->neighbors[1] = delaunayNew[2];
    delaunayNew[2]->neighbors[1] = delaunayNew[0];

//    faceClass *face = new faceClass();
    faceClass *face = &staticfaces[iiface]; ++iiface;

    delaunayNew[0]->faces[1] = face;
    delaunayNew[2]->faces[1] = face;
    face->interDel[0] = delaunayNew[0];
    face->interDel[1] = delaunayNew[2];
    face->interidx[0] = 1;
    face->interidx[1] = 1;
    face->inbrs[0] = 1;
    face->inbrs[1] = 1;

    delaunayNew[1]->neighbors[1] = delaunayNew[3];
    delaunayNew[3]->neighbors[1] = delaunayNew[1];

//    face = new faceClass();
    face = &staticfaces[iiface]; ++iiface;

    delaunayNew[1]->faces[1] = face;
    delaunayNew[3]->faces[1] = face;
    face->interDel[0] = delaunayNew[1];
    face->interDel[1] = delaunayNew[3];
    face->interidx[0] = 1;
    face->interidx[1] = 1;
    face->inbrs[0] = 1;
    face->inbrs[1] = 1;

  } else 
  {
    delaunayNew[0]->neighbors[1] = delaunayNew[3];
    delaunayNew[3]->neighbors[1] = delaunayNew[0];

    //faceClass *face = new faceClass();
    faceClass *face = &staticfaces[iiface]; ++iiface;

    delaunayNew[0]->faces[1] = face;
    delaunayNew[3]->faces[1] = face;
    face->interDel[0] = delaunayNew[0];
    face->interDel[1] = delaunayNew[3];
    face->interidx[0] = 1;
    face->interidx[1] = 1;
    face->inbrs[0] = 1;
    face->inbrs[1] = 1;
    
    delaunayNew[1]->neighbors[1] = delaunayNew[2];
    delaunayNew[2]->neighbors[1] = delaunayNew[1]; 

    //face = new faceClass();
    face = &staticfaces[iiface]; ++iiface;

    delaunayNew[1]->faces[1] = face;
    delaunayNew[2]->faces[1] = face;
    face->interDel[0] = delaunayNew[1];
    face->interDel[1] = delaunayNew[2];
    face->interidx[0] = 1;
    face->interidx[1] = 1;
    face->inbrs[0] = 1;
    face->inbrs[1] = 1;

  }
 
  delaunayNew[0]->nNeighbor = 2; delaunayNew[1]->nNeighbor = 2;
  delaunayNew[2]->nNeighbor = 2; delaunayNew[3]->nNeighbor = 2;

#if 0

  vector<int> indexs, indexs2;
  vector<int> indexofnbr, indexofnbr2;
  vector<int> indexofa, indexofa2;
  vector<vector<int> > oldnbrVec;

  for(int ia = 0; ia < 4; ++ia)
  {
    vector<int> tmpVec;
    delaunayClass *tmpOld = &oldDelaunay[ia];
    for(int inbr = 0; inbr < tmpOld->nNeighbor; ++inbr)
    {
      tmpVec.push_back(tmpOld->neighbors[inbr]->indd);
    }
    oldnbrVec.push_back(tmpVec);
  }

  for(int i = 0; i < 4; ++i)
  {
    for(int inbr = 0; inbr < oldDelaunay[i].nNeighbor; ++inbr)
    {
      //delaunayClass *tmpNbr = oldDelaunay[i].neighbors[inbr]; 
      delaunayClass *tmpNbr = &Delaunay[oldnbrVec[i][inbr]];
      //cout << oldnbrVec[i][inbr] << endl;
      if(tmpNbr->indd == oldDelaunay[0].indd || tmpNbr->indd == oldDelaunay[1].indd || tmpNbr->indd == oldDelaunay[2].indd || tmpNbr->indd == oldDelaunay[3].indd ) continue;
      
      if(!(tmpNbr->has_mesh(no_indm[0]) && tmpNbr->has_mesh(no_indm[1])))
      {
	int index = -1;
	if(tmpNbr->nNeighbor >= 5) getchar();
	for(int l = 0; l < tmpNbr->nNeighbor; ++l)
	{
	  if(tmpNbr->neighbors[l]->indd == oldDelaunay[i].indd)
	  {
	    index = l; 
	  }
	} 

	int ia = 0;

	if(delaunayNew[ia]->nNeighbor >= 5) 
	{
	  cout << "warning " <<  delaunayNew[ia]->nNeighbor << endl;
	  getchar();
	}
	
	for(ia = 0; ia < 4; ++ia)
	{
	  if(delaunayNew[ia]->nNeighbor == 4) continue;
           if(tmpNbr->has_mesh(delaunayNew[ia]->indm[0]) && tmpNbr->has_mesh(delaunayNew[ia]->indm[1]) && 
	       delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor - 1] != tmpNbr)  break;
	}

	delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = tmpNbr;

	indexofa2.push_back(delaunayNew[ia]->indd);
	indexofnbr2.push_back(delaunayNew[ia]->nNeighbor);
	indexs2.push_back(tmpNbr->indd);
	++delaunayNew[ia]->nNeighbor;

	indexs.push_back(tmpNbr->indd);
	indexofnbr.push_back(index);
	indexofa.push_back(delaunayNew[ia]->indd);

      }
    }
  }
   
  //cout << "end gene" << endl;

#ifndef _DEBUG
  for(int i = 0; i < indexs2.size(); ++i)
  {
    Delaunay[indexofa2[i]].neighbors[indexofnbr2[i]] = &Delaunay[indexs2[i]]; 
  }

  for(int i = 0; i < indexs.size(); ++i) 
  {
     Delaunay[indexs[i]].neighbors[indexofnbr[i]] = &Delaunay[indexofa[i]];
  }
  
#endif

  #endif 

#ifdef _DEBUG 
  int nNeighbor[4] = {2,2,2,2};

  //cout << "4to4" << endl;
  for(int ia = 0; ia < 4; ++ia)
  {
    for(int j = (ia < 2 ? 0 : 1); j < 4; j += 2)
    {
      int idl = 0;
      faceClass *face = NULL;    

      if(delaunayNew[ia]->has_mesh(no_indm[0])) 
      {
	face = oldDelaunay[j].faces[sharedEdgeIdx[j][1]];
      } else
      {
	face = oldDelaunay[j].faces[sharedEdgeIdx[j][0]];
      }
      
      if(face == NULL){continue;}
   
      delaunayClass *interDel = face->interDel[0];
      
      if(interDel->indd == oldDelaunay[j].indd)
      {interDel = face->interDel[1]; idl = 1;} 
     
      int inbr = face->inbrs[idl];

      interDel->neighbors[inbr] = delaunayNew[ia];

      idl = idl^1; 
      if(oldDelaunay[j].has_mesh(imesh))
      {
	//delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = interDel; 
	delaunayNew[ia]->neighbors[nNeighbor[ia]] = interDel;
	delaunayNew[ia]->faces[2] = face;

	face->interDel[idl] = delaunayNew[ia];
	face->interidx[idl] = 2;
	//face->inbrs[idl] = delaunayNew[ia]->nNeighbor;
	//++delaunayNew[ia]->nNeighbor;
	face->inbrs[idl] = nNeighbor[ia];
	++nNeighbor[ia];
      } else if(oldDelaunay[j].has_mesh(jmesh))
      {

   	//delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = interDel;
	delaunayNew[ia]->neighbors[nNeighbor[ia]] = interDel;
	delaunayNew[ia]->faces[3] = face;

	face->interDel[idl] = delaunayNew[ia];
	face->interidx[idl] = 3;
	//face->inbrs[idl] = delaunayNew[ia]->nNeighbor;
	//++delaunayNew[ia]->nNeighbor;
	face->inbrs[idl] = nNeighbor[ia];
	++nNeighbor[ia];
      }
      
    }
  }

  for(int i = 0; i < 4; ++i)
    delaunayNew[i]->nNeighbor = nNeighbor[i];
#endif

    //tmptime = clock() - tmptime;
  //time4 += tmptime;
  
  //tmptime = clock();
  /*
  for(int i = 0; i < 4; ++i)
  {
    oldDelaunay[i].eraseDelaunaryOfMesh(Mesh, params);
  }
  for(int i = 0; i < 4; ++i)
  {
    delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }
   */

  //tmptime = clock() - tmptime;
  //time7 += tmptime;
  
  //cout << "end 4to4" << endl;
  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 4, params);

  //rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 4, params);
  for(int i = 0; i < 4; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
    //stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
    //stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
    tmpStackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
  }
  lastDel = NewIndex[3];
  
  return rval;
}

bool coplanar(int A, int B, int imesh, int jmesh, meshClass *Mesh)
{

  bool flag = false;
  double a[3] = {Mesh[A].getPosx(), Mesh[A].getPosy(), Mesh[A].getPosz()};
  double b[3] = {Mesh[B].getPosx(), Mesh[B].getPosy(), Mesh[B].getPosz()};
  double c[3] = {Mesh[jmesh].getPosx(), Mesh[jmesh].getPosy(), Mesh[jmesh].getPosz()};
  double d[3] = {Mesh[imesh].getPosx(), Mesh[imesh].getPosy(), Mesh[imesh].getPosz()};

  //Orient = orient3d(a,b,c,d);  
#if 0
  Vec avec12 = Vec::getVec12(a[0], a[1], a[2]);
  Vec bvec12 = Vec::getVec12(b[0], b[1], b[2]);
  Vec cvec12 = Vec::getVec12(c[0], c[1], c[2]);
  Vec dvec12 = Vec::getVec12(d[0], d[1], d[2]);
#endif
  avec12.setVec12(a[0], a[1], a[2]);
  bvec12.setVec12(b[0], b[1], b[2]);
  cvec12.setVec12(c[0], c[1], c[2]);
  dvec12.setVec12(d[0], d[1], d[2]);

  double test1 = ExactArithmetic::orient3d(avec12, bvec12, cvec12, dvec12);

  //if(fabs(Orient) <= 1.0e-16) flag = true;
  if(test1 == 0) flag = true;
  return flag;
}

int delaunaysplit3dnto2n(int it, int imesh, meshClass *Mesh, delaunayClass *Delaunay, paramsClass &params)
{
  int rval = 0, nsize = 0, no_indm[2]={-1, -1};

  //clock_t tmptime = clock();
//  cout << "nto2n" << endl;
 
  //1. 
  rval = getDelaunaySharingEdge(no_indm, nsize, imesh, it, oldDelaunay, Delaunay, Mesh);
  
  //2
  delaunayClass *delaunayNew[2*nsize];
  int NewIndex[2*nsize];

  for(int i = 0, ia = 0; i < nsize; ++i)
  {
    delaunayNew[ia] = &Delaunay[oldDelaunay[i].indd];  
    delaunayNew[ia+1] = &Delaunay[params.Ndelaunay + i];
    ia+=2;
  }

  for(int ia = 0, i = 0; i < nsize; ++i)
  {
   for(int j = 0; j < 4; ++j)
   {
     for(int k = j + 1; k < 4; ++k)
     {
       if((oldDelaunay[i].indm[j] == no_indm[0] || oldDelaunay[i].indm[k] == no_indm[1]) || 
	  (oldDelaunay[i].indm[j] == no_indm[1] || oldDelaunay[i].indm[k] == no_indm[0]))  
	 continue;

       for(int l = 0; l < 2; ++l )
       {
	 NewIndex[ia] = l == 0 ? oldDelaunay[i].indd : params.Ndelaunay + i;
	 delaunayNew[ia]->set(NewIndex[ia], oldDelaunay[i].indm[j], oldDelaunay[i].indm[k], no_indm[l], imesh, Mesh);
	 delaunayNew[ia]->circumcenter3d();
#if 0
	 if(oldDelaunay[i].indm[j] == no_indm[0] || oldDelaunay[i].indm[j] == no_indm[0])
	   delaunayNew[ia]->set(NewIndex[ia], oldDelaunay[i].indm[j], oldDelaunay[i].indm[k], no_indm[l], imesh, Mesh);
	 else if(oldDelaunay[i].indm[k] == no_indm[1] || oldDelaunay[i].indm[k] == no_indm[1])
	   delaunayNew[ia]->set(NewIndex[ia], oldDelaunay[i].indm[k], oldDelaunay[i].indm[j], no_indm[l], imesh, Mesh);
#endif	

	 delaunayNew[ia]->birthtype = "nto2n";
	 ++ia;
       }
     }
   } 
  } 

  params.Ndelaunay += nsize;

  //face2
  for(int ia = 0; ia < nsize; ++ia)
  {
    //faceClass *face = new faceClass();
    faceClass *face = &staticfaces[iiface]; ++iiface;

    delaunayNew[2*ia]->neighbors[delaunayNew[2*ia]->nNeighbor] = delaunayNew[2*ia + 1];
    delaunayNew[2*ia + 1]->neighbors[delaunayNew[2*ia + 1]->nNeighbor] = delaunayNew[2*ia];

    delaunayNew[2*ia]->faces[2] = face;
    delaunayNew[2*ia + 1]->faces[2] = face;
    
    face->interDel[0] = delaunayNew[2*ia];
    face->interDel[1] = delaunayNew[2*ia+1];
    face->interidx[0] = 2;
    face->interidx[1] = 2;
    face->inbrs[0] = delaunayNew[2*ia]->nNeighbor;
    face->inbrs[1] = delaunayNew[2*ia + 1]->nNeighbor; 
    
    delaunayNew[2*ia]->nNeighbor = 1;
    delaunayNew[2*ia + 1]->nNeighbor = 1;
  }

  //face 0 or 1
  for(int i = 0; i < nsize; ++i)
  {
    for(int j = i + 1; j < nsize; ++j)
    {
      
      if(delaunayNew[2*j]->has_mesh(delaunayNew[2*i]->indm[0]) || delaunayNew[2*j]->has_mesh(delaunayNew[2*i]->indm[1]))
      {
        delaunayNew[2*i]->neighbors[delaunayNew[2*i]->nNeighbor] = delaunayNew[2*j]; 
	delaunayNew[2*j]->neighbors[delaunayNew[2*j]->nNeighbor] = delaunayNew[2*i];

	//faceClass *face = new faceClass();

	faceClass *face = &staticfaces[iiface]; ++iiface;

	if(delaunayNew[2*i]->indm[0] == delaunayNew[2*j]->indm[0] || 
	   delaunayNew[2*i]->indm[0] == delaunayNew[2*j]->indm[1])
	  {
	    if(delaunayNew[2*i]->indm[0] == delaunayNew[2*j]->indm[0])
	    {
	      
	      delaunayNew[2*i]->faces[1] = face;
	      face->interDel[0] = delaunayNew[2*i];
	      face->interidx[0] = 1;
	      face->inbrs[0] = delaunayNew[2*i]->nNeighbor;

	      delaunayNew[2*j]->faces[1] = face;
	      face->interDel[1] = delaunayNew[2*j];
	      face->interidx[1] = 1;
	      face->inbrs[1] = delaunayNew[2*j]->nNeighbor;
	      
	    } else if(delaunayNew[2*i]->indm[0] == delaunayNew[2*j]->indm[1])
	    {
	      
	      delaunayNew[2*i]->faces[1] = face;  
	      face->interDel[0] = delaunayNew[2*i];
	      face->interidx[0] = 1;
	      face->inbrs[0] = delaunayNew[2*i]->nNeighbor;

	      delaunayNew[2*j]->faces[0] = face;
	      face->interDel[1] = delaunayNew[2*j];
	      face->interidx[1] = 0;
	      face->inbrs[1] = delaunayNew[2*j]->nNeighbor; 
	    }
	  } else if(delaunayNew[2*i]->indm[1] == delaunayNew[2*j]->indm[0] ||
	            delaunayNew[2*i]->indm[1] == delaunayNew[2*j]->indm[1])
	  {
	     if(delaunayNew[2*i]->indm[1] == delaunayNew[2*j]->indm[0])
	    {
	      delaunayNew[2*i]->faces[0] = face;
	      face->interDel[0] = delaunayNew[2*i];
	      face->interidx[0] = 0;
	      face->inbrs[0] = delaunayNew[2*i]->nNeighbor;

	      delaunayNew[2*j]->faces[1] = face;
	      face->interDel[1] = delaunayNew[2*j];
	      face->interidx[1] = 1;
	      face->inbrs[1] = delaunayNew[2*j]->nNeighbor;
	      
	    } else if(delaunayNew[2*i]->indm[1] == delaunayNew[2*j]->indm[1])
	    {
	      delaunayNew[2*i]->faces[0] = face;  
	      face->interDel[0] = delaunayNew[2*i];
	      face->interidx[0] = 0;
	      face->inbrs[0] = delaunayNew[2*i]->nNeighbor;

	      delaunayNew[2*j]->faces[0] = face;
	      face->interDel[1] = delaunayNew[2*j];
	      face->interidx[1] = 0;
	      face->inbrs[1] = delaunayNew[2*j]->nNeighbor; 
	    }
	  }

	delaunayNew[2*i]->nNeighbor += 1;
	delaunayNew[2*j]->nNeighbor += 1;
      }	

      if(delaunayNew[2*j+1]->has_mesh(delaunayNew[2*i+1]->indm[0]) || delaunayNew[2*j+1]->has_mesh(delaunayNew[2*i+1]->indm[1]))
      {
        delaunayNew[2*i+1]->neighbors[delaunayNew[2*i+1]->nNeighbor] = delaunayNew[2*j+1]; 
	delaunayNew[2*j+1]->neighbors[delaunayNew[2*j+1]->nNeighbor] = delaunayNew[2*i+1];

	//faceClass *face = new faceClass();
	faceClass *face = &staticfaces[iiface]; ++iiface;

	if(delaunayNew[2*i + 1]->indm[0] == delaunayNew[2*j + 1]->indm[0] || 
	   delaunayNew[2*i + 1]->indm[0] == delaunayNew[2*j + 1]->indm[1])
	  {
	    if(delaunayNew[2*i + 1]->indm[0] == delaunayNew[2*j + 1]->indm[0])
	    {
	      delaunayNew[2*i + 1]->faces[1] = face;
	      face->interDel[0] = delaunayNew[2*i + 1];
	      face->interidx[0] = 1;
	      face->inbrs[0] = delaunayNew[2*i + 1]->nNeighbor;

	      delaunayNew[2*j + 1]->faces[1] = face;
	      face->interDel[1] = delaunayNew[2*j + 1];
	      face->interidx[1] = 1;
	      face->inbrs[1] = delaunayNew[2*j + 1]->nNeighbor;
	      
	    } else if(delaunayNew[2*i + 1]->indm[0] == delaunayNew[2*j + 1]->indm[1])
	    {
	      delaunayNew[2*i + 1]->faces[1] = face;  
	      face->interDel[0] = delaunayNew[2*i + 1];
	      face->interidx[0] = 1;
	      face->inbrs[0] = delaunayNew[2*i + 1]->nNeighbor;

	      delaunayNew[2*j + 1]->faces[0] = face;
	      face->interDel[1] = delaunayNew[2*j + 1];
	      face->interidx[1] = 0;
	      face->inbrs[1] = delaunayNew[2*j + 1]->nNeighbor; 
	    }
	  } else if(delaunayNew[2*i + 1]->indm[1] == delaunayNew[2*j + 1]->indm[0] ||
	            delaunayNew[2*i + 1]->indm[1] == delaunayNew[2*j + 1]->indm[1])
	  {
	     if(delaunayNew[2*i + 1]->indm[1] == delaunayNew[2*j + 1]->indm[0])
	    {
	      delaunayNew[2*i + 1]->faces[0] = face;
	      face->interDel[0] = delaunayNew[2*i + 1];
	      face->interidx[0] = 0;
	      face->inbrs[0] = delaunayNew[2*i + 1]->nNeighbor;

	      delaunayNew[2*j + 1]->faces[1] = face;
	      face->interDel[1] = delaunayNew[2*j + 1];
	      face->interidx[1] = 1;
	      face->inbrs[1] = delaunayNew[2*j + 1]->nNeighbor;
	      
	    } else if(delaunayNew[2*i + 1]->indm[1] == delaunayNew[2*j  +1]->indm[1])
	    {
	      delaunayNew[2*i + 1]->faces[0] = face;  
	      face->interDel[0] = delaunayNew[2*i + 1];
	      face->interidx[0] = 0;
	      face->inbrs[0] = delaunayNew[2*i + 1]->nNeighbor;

	      delaunayNew[2*j + 1]->faces[0] = face;
	      face->interDel[1] = delaunayNew[2*j + 1];
	      face->interidx[1] = 0;
	      face->inbrs[1] = delaunayNew[2*j + 1]->nNeighbor; 
	    }
	  }

	delaunayNew[2*i+1]->nNeighbor += 1;
	delaunayNew[2*j+1]->nNeighbor += 1;

      }	
    } 
  }

  //tmptime = clock() - tmptime;
  //time3 += tmptime;
  //tmptime = clock();
#if 0

  vector<int> indexs, index2;
  vector<int> indexofnbr, indexofnbr2;
  vector<int> indexofa, indexofa2;

  map<int, int> tmpIndexMap;
  for(int i = 0; i < nsize; ++i)
  {
    tmpIndexMap[oldDelaunay[i].indd] = oldDelaunay[i].indd; 
  }

  for(int i = 0; i < nsize; ++i)
  {
    for(int inbr = 0; inbr < oldDelaunay[i].nNeighbor; ++inbr)
    {
      delaunayClass *tmpNbr = oldDelaunay[i].neighbors[inbr];
      if(tmpIndexMap.find(tmpNbr->indd) != tmpIndexMap.end()) continue;

      int index = 0;
      for(int l = 0; l < tmpNbr->nNeighbor; ++l)
      {
        if(tmpNbr->neighbors[l]->indd == oldDelaunay[i].indd){index= l; break;}
      }
	
      int ia = -1;
      if(tmpNbr->has_mesh(no_indm[0])) ia  = 2 * i; 
      else if(tmpNbr->has_mesh(no_indm[1])) ia = 2 * i + 1;
      else {cout << "warning" << endl;getchar();}

      //delaunayNew[ia]->neighbors[delaunayNew[ia]->nNeighbor] = tmpNbr;
      indexofa2.push_back(ia); 
      indexofnbr2.push_back(delaunayNew[ia]->nNeighbor);
      index2.push_back(tmpNbr->indd);

      ++delaunayNew[ia]->nNeighbor;
      indexs.push_back(tmpNbr->indd);
      indexofnbr.push_back(index);
      indexofa.push_back(delaunayNew[ia]->indd);
    } 
  }

#ifndef _DEBUG

  for(int i = 0; i < index2.size(); ++i)
  {
    delaunayNew[indexofa2[i]]->neighbors[indexofnbr2[i]] = &Delaunay[index2[i]];
  }
 
  for(int i = 0; i < indexs.size(); ++i) 
  {
     Delaunay[indexs[i]].neighbors[indexofnbr[i]] = &Delaunay[indexofa[i]];
  }

#endif

#endif
  int sharedEdgeIdx[100][2];

  for(int ia = 0; ia < nsize; ++ia)
  {
    int ii = 0;
    for(int j = 0; j < 4; ++j)
    {
      if(oldDelaunay[ia].indm[j] == no_indm[0] || oldDelaunay[ia].indm[j] == no_indm[1])
      {
	if(oldDelaunay[ia].indm[j] == no_indm[0])
	  sharedEdgeIdx[ia][0] = j;
	else if(oldDelaunay[ia].indm[j] == no_indm[1])
	  sharedEdgeIdx[ia][1] = j;
	++ii;
	if(ii == 2) break;
      }	
    }
  }

#ifdef _DEBUG

  int nNeighbor[100];
  //cout <<"nto2n" << endl;

  for(int ia = 0; ia < 2*nsize; ++ia)
    nNeighbor[ia] = 3;

  for(int ia = 0; ia < 2*nsize; ++ia)
  {
    int idl = 0;
    faceClass *face = NULL; 

    if(delaunayNew[ia]->has_mesh(oldDelaunay[ia/2].indm[sharedEdgeIdx[ia/2][0]]))
    {
      face = oldDelaunay[ia/2].faces[sharedEdgeIdx[ia/2][1]];
    } else
    {
      face = oldDelaunay[ia/2].faces[sharedEdgeIdx[ia/2][0]];
    }
    
    if(face == NULL) continue;

    delaunayClass *interDel = face->interDel[0];

    if(interDel->indd == oldDelaunay[ia/2].indd)
    {interDel = face->interDel[1]; idl = 1;} 

    int inbr = face->inbrs[idl];

    interDel->neighbors[inbr] = delaunayNew[ia];

    idl = idl^1;
    delaunayNew[ia]->neighbors[3] = interDel;
    delaunayNew[ia]->faces[3] = face;
    face->interidx[idl] = 3;
    face->inbrs[idl] = 3;
    face->interDel[idl] = delaunayNew[ia];
    //++delaunayNew[ia]->nNeighbor;
    nNeighbor[ia]++;
  }

  for(int i = 0;i < 2*nsize; ++i)
  {
    delaunayNew[i]->nNeighbor = nNeighbor[i]; 
  }

#endif

  //tmptime = clock() - tmptime;
  //time4 += tmptime;

  //tmptime = clock();
  /*
  for(int i = 0; i < nsize; ++i)
  {
    oldDelaunay[i].eraseDelaunaryOfMesh(Mesh, params);
  }
  for(int i = 0; i < 2*nsize; ++i)
  {
    delaunayNew[i]->addDelaunayToMesh(Mesh, params);
  }
  */

  //tmptime = clock() - tmptime;
  //time4 += tmptime;
  //rval = check_flip(it, imesh, Mesh, Delaunay, delaunayNew, 2*nsize, params);
  //rval = check_flip(it, imesh, Mesh, Delaunay, NewIndex, 2*nsize, params);
  for(int i = 0; i < 2*nsize; ++i)
  {
    delaunayClass *tmpDelaunay = &Delaunay[NewIndex[i]];
//    stackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);
    //stackDelaunays[istack] = &Delaunay[tmpDelaunay->indd]; ++istack;
tmpStackDelaunays.push_back(&Delaunay[tmpDelaunay->indd]);

  }
    lastDel = NewIndex[2*nsize-1];



  return rval;
}


bool onEdge(int A, int B, int C, meshClass *Mesh, double &eps)
{
  bool flag = false;

  double r1[3] = {Mesh[A].getPosx() - Mesh[B].getPosx(), Mesh[A].getPosy() - Mesh[B].getPosy(), Mesh[A].getPosz() - Mesh[B].getPosz()};
  double r2[3] = {Mesh[A].getPosx() - Mesh[C].getPosx(), Mesh[A].getPosy() - Mesh[C].getPosy(), Mesh[A].getPosz() - Mesh[C].getPosz()};

  double dot = (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]);
  dot = dot * dot / (L(r1[0], r1[1], r1[2])) / (L(r2[0], r2[1], r2[2]));
  if(fabs(1.0 - dot) <= 1.0e-6)
  {
    flag = true;
  }
  eps = fabs(1.0-dot);

  return flag;
}

int getDelaunaySharingEdge(int *no_indm, int &nsize, int imesh, int it, delaunayClass *oldDelaunay, delaunayClass *Delaunay, meshClass *Mesh)
{
  int rval = 0;
  double epsmin = 1.0e+10;
  // 1.decide mesh index of sharing edge 
  /*
  for(int i = 0; i < 4; ++i)
  {
    for(int j = i + 1; j < 4; ++j)
    {
      double eps = 1.0e+10;
      //if(onEdge(Delaunay[it].indm[i], Delaunay[it].indm[j], imesh, Mesh, eps))
      onEdge(Delaunay[it].indm[i], Delaunay[it].indm[j], imesh, Mesh, eps);

      if(eps < epsmin)
      {
         no_indm[0] = Delaunay[it].indm[i]; 
	 no_indm[1] = Delaunay[it].indm[j];
	 epsmin = eps;
      }	
    } 
  }*/

  if(!(((no_indm[0] == men[0][0] || no_indm[0] == men[0][1] || no_indm[0] == men[0][2]) && 
     (no_indm[0] == men[1][0] || no_indm[0] == men[1][1] || no_indm[0] == men[1][2]))
     &&
     ((no_indm[1] == men[0][0] || no_indm[1] == men[0][1] || no_indm[1] == men[0][2]) && 
     (no_indm[1] == men[1][0] || no_indm[1] == men[1][1] || no_indm[1] == men[1][2]))))
  {
   int cc = 0;

   if(men[0][0] == men[1][0] || men[0][0] == men[1][1] || men[0][0] == men[1][2])
   { 
     no_indm[cc] = men[0][0];
     ++cc;
   }

   if(men[0][1] == men[1][0] || men[0][1] == men[1][1] || men[0][1] == men[1][2])
   { 
     no_indm[cc] = men[0][1];
     ++cc;
   }
   
   if(men[0][2] == men[1][0] || men[0][2] == men[1][1] || men[0][2] == men[1][2])
   { 
     no_indm[cc] = men[0][2];
     ++cc;
   } 
  }

  //2. get delaunay
  //for(int i = 0; i < Mesh[no_indm[0]].delaunay.size(); ++i)
  /*for(unordered_map<int, delaunayClass*>::iterator itr = Mesh[no_indm[0]].delaunaymap.begin(); itr != Mesh[no_indm[0]].delaunaymap.end();  ++itr)
  {
     //delaunayClass *tmpDelaunay = Mesh[no_indm[0]].delaunay[i];  
    delaunayClass *tmpDelaunay = itr->second;
     if(tmpDelaunay->has_mesh(no_indm[1]))
     {
       int index = tmpDelaunay->indd;
       oldDelaunay[nsize] = Delaunay[index];
       ++nsize;
     }
  }*/

  //int old = nsize;
  nsize = 0;
  delaunayClass *preDel = NULL, *Del = NULL;
  Del = &Delaunay[it];
  while(1)
  {
    bool flag = false;
     
    for(int inbr = 0; inbr < Del->nNeighbor; ++inbr)
    {
      delaunayClass *nbr = Del->neighbors[inbr];

      if(nbr->has_mesh(no_indm[0]) && nbr->has_mesh(no_indm[1]) && (preDel == NULL || nbr->indd != preDel->indd))
      {
	oldDelaunay[nsize].clone(Delaunay[Del->indd]);

	preDel = Del;
	Del = nbr;	

	++nsize;
	break;
      }

    }

    if(Del->indd == Delaunay[it].indd) 
      break;
     
  }

  if(nsize >= 50)
  {
    cout << "warning" << endl;
    getchar();
  }

  //cout << "end " << endl;

  return rval;

}
/*-------------------*/
int ch4to4(int it, int imesh, int jmesh, meshClass *Mesh, int in, delaunayClass *Delaunay, paramsClass &params)
{

  int rval = 0;
  //delaunayClass oldDelaunay[4];

  oldDelaunay[0] = Delaunay[it]; oldDelaunay[2] = Delaunay[in];

  int no_indm[2];
  bool flag = false;

  delaunayClass *Nbr;
  for(int inbr = 0; inbr < oldDelaunay[0].nNeighbor; ++inbr)
  {
    Nbr = oldDelaunay[0].neighbors[inbr]; 
    if(!Nbr->has_mesh(jmesh) && Nbr->has_mesh(imesh))
    {
      for(int i = 0; i < 3; ++i)
      {
         for(int j = i + 1; j < 3; ++j)
	 {
	   if(Nbr->has_mesh(oldDelaunay[0].indm[i]) && Nbr->has_mesh(oldDelaunay[0].indm[j]))
	   {
	     if(!(coplanar(oldDelaunay[0].indm[i], oldDelaunay[0].indm[j], jmesh, imesh, Mesh))) 
		 continue;
	     no_indm[0] = oldDelaunay[0].indm[i];
	     no_indm[1] = oldDelaunay[0].indm[j];
	     oldDelaunay[1] = Delaunay[Nbr->indd];
	     flag = true;
	     goto NEXT1;
	   }
	 }
      }	
    }
  } 

NEXT1:;

  if(!flag){ rval = 3; return rval;}
  // get oldDelaunay[3]
  flag = false;
  for(int inbr = 0; inbr < oldDelaunay[2].nNeighbor; ++inbr)
  {

    Nbr = oldDelaunay[2].neighbors[inbr]; 
    if(!Nbr->has_mesh(imesh) && Nbr->has_mesh(jmesh) && Nbr->has_mesh(no_indm[0]) && Nbr->has_mesh(no_indm[1]))
    {
      oldDelaunay[3] = Delaunay[Nbr->indd]; 
      flag = true;
    }
  }
  
  if(!flag) { rval = 3; return rval;}

  flag = false;
  for(int inbr = 0; inbr < oldDelaunay[1].nNeighbor; ++inbr)
  {
     if(oldDelaunay[3].indd == oldDelaunay[1].neighbors[inbr]->indd) 
     {flag = true; break;} 
  }

  if(!flag){ rval = 3;return rval;}

  rval = 4;
  return rval;
}

void settingNextMesh(meshClass* Mesh, Vec *hilbertCurves, int nGrid, paramsClass &params)
{
   const int nGrid3 = cube(nGrid);
   const double boxsize = params.periodic ? BOXSIZE + 2 * DX : BOXSIZE;
   const double dx = boxsize / nGrid;
   const double buf = params.periodic ? DX : 0;

   int *indC = new int[nGrid3];
   int ind;

   for(int i = 0; i < nGrid3; ++i)
   {
     int ix = hilbertCurves[i].x(); 
     int iy = hilbertCurves[i].y();
     int iz = hilbertCurves[i].z();
     ind = ix + nGrid * (iy + nGrid*iz);
     indC[ind] = i;
     hilbertCurves[i].first = NULL;
   }

   for(int imesh = 4; imesh < params.Ngas + 4 + params.Nghost; ++imesh)
   {
     Mesh[imesh].next = NULL; 
   }

   for(int imesh = 4; imesh < params.Ngas + 4 + params.Nghost; ++imesh)
   {
     meshClass *mesh = &Mesh[imesh];
     int ix = (mesh->getPosx() + buf) / dx;
     int iy = (mesh->getPosy() + buf) / dx;
     int iz = (mesh->getPosz() + buf) / dx;
     ind = ix + nGrid * (iy + nGrid * iz);

     int ih = indC[ind];

     meshClass *first = hilbertCurves[ih].first;
     hilbertCurves[ih].first = mesh;
     mesh->next = first;
   }
   
   delete [] indC;

}
    
faceClass::faceClass(){
  //facevec.push_back(this);
} 

splitType delaunayClass::getSplitType_WalkNeighbor(int &it, meshClass *mesh, paramsClass &params)
{
  
  int iterator = 0, rval;
  const int ndim = params.dim;

  it = -1;

  splitType splittype = _nosplit;
  delaunayClass *tmpDelaunay = this;

  splittype = tmpDelaunay->getSplitType(mesh->getPosx(), mesh->getPosy(), mesh->getPosz());
  
  if(splittype != _nosplit)
  {
    it = tmpDelaunay->indd;
    return splittype;
  }

  while(1)
  {
    double xg = 0.0, yg = 0.0, zg = 0.0, s = 0.0, t = 0.0, u = 0.0;

    if(ndim == 2)
    {
      xg = (tmpDelaunay->getPosx(0) + tmpDelaunay->getPosx(1) + tmpDelaunay->getPosx(2)) / 3.0;
      yg = (tmpDelaunay->getPosy(0) + tmpDelaunay->getPosy(1) + tmpDelaunay->getPosy(2)) / 3.0;
    } else if( ndim == 3)
    {
      xg = (tmpDelaunay->getPosx(0) + tmpDelaunay->getPosx(1) + tmpDelaunay->getPosx(2) + tmpDelaunay->getPosx(3)) / 4.0;
      yg = (tmpDelaunay->getPosy(0) + tmpDelaunay->getPosy(1) + tmpDelaunay->getPosy(2) + tmpDelaunay->getPosy(3)) / 4.0;
      zg = (tmpDelaunay->getPosz(0) + tmpDelaunay->getPosz(1) + tmpDelaunay->getPosz(2) + tmpDelaunay->getPosz(3)) / 4.0;
    }

    int id[3];
    if(ndim == 2)
    {

      for(int i = 0; i < 3; ++i)
      {
	for(int j = i + 1; j < 3; ++j)
	{
	  id[0] = tmpDelaunay->indm[i]; 
	  id[1] = tmpDelaunay->indm[j];
	  double p[2], a[2], b[2], det;
	  double a_inv[2], b_inv[2];

	  p[0] = mesh->getPosx() - xg;
	  p[1] = mesh->getPosy() - yg;

	  a[0] = tmpDelaunay->getPosx(i) - xg;
	  b[0] = tmpDelaunay->getPosy(i) - yg;

	  a[1] = tmpDelaunay->getPosx(j) - xg;
	  b[1] = tmpDelaunay->getPosy(j) - yg;

	  det = DET2D(a[0], a[1], b[0], b[1]);
	  a_inv[0] = b[1];
	  b_inv[0] = -a[1];
	  a_inv[1] = -b[0];
	  b_inv[1] = a[0];

	  s = (a_inv[0] * p[0] + b_inv[0] * p[1]) / det;
	  t = (a_inv[1] * p[0] + b_inv[1] * p[1]) / det;

	  if(s >= 0.0 && t >= 0.0) goto NEXT; 
	}
      }
    } else if (ndim == 3)
    {

      for(int i = 0; i < 4; ++i)
      {
	for(int j = i + 1; j < 4; ++j)
	{
	  for(int k = j + 1; k < 4; ++k)
	  {
	    id[0] = tmpDelaunay->indm[i]; id[1] =  tmpDelaunay->indm[j]; id[2] = tmpDelaunay->indm[k]; 
	    double p[3],a[3], b[3], c[3];
	    double a_inv[3], b_inv[3], c_inv[3];
	    p[0] = mesh->getPosx() - xg;
	    p[1] = mesh->getPosy() - yg;
	    p[2] = mesh->getPosz() - zg;
	    a[0] = tmpDelaunay->getPosx(i) - xg;
	    b[0] = tmpDelaunay->getPosy(i) - yg;
	    c[0] = tmpDelaunay->getPosz(i) - zg;
	    a[1] = tmpDelaunay->getPosx(j) - xg;
	    b[1] = tmpDelaunay->getPosy(j) - yg;
	    c[1] = tmpDelaunay->getPosz(j) - zg;
	    a[2] = tmpDelaunay->getPosx(k) - xg;
	    b[2] = tmpDelaunay->getPosy(k) - yg;
	    c[2] = tmpDelaunay->getPosz(k) - zg;

	    double det = DET3D(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
	    a_inv[0] = (b[1]*c[2] - b[2]*c[1]); a_inv[1] = (b[0]*c[2]-b[2]*c[0]) * (-1.0) ; a_inv[2] = (b[0]*c[1] - b[1]*c[0]);
	    b_inv[0] = -(a[1]*c[2] - a[2]*c[1]); b_inv[1] = (a[0]*c[2]-a[2]*c[0])         ; b_inv[2] = -(a[0]*c[1] - a[1]*c[0]);
	    c_inv[0] = (a[1]*b[2] - a[2]*b[1]); c_inv[1] = (a[0]*b[2]-a[2]*b[0]) * (-1.0) ; c_inv[2] = (a[0]*b[1] - a[1]*b[0]);

	    s = (a_inv[0] * p[0] + b_inv[0] * p[1] + c_inv[0] * p[2]) / det;
	    t = (a_inv[1] * p[0] + b_inv[1] * p[1] + c_inv[1] * p[2]) / det;
	    u = (a_inv[2] * p[0] + b_inv[2] * p[1] + c_inv[2] * p[2]) / det;
	    if(s >= 0 && t >= 0 && u >= 0) goto NEXT;
	  }
	} 
      }

    }
NEXT:;

     for(int inbr = 0; inbr < tmpDelaunay->nNeighbor; ++inbr)
     {
       if(ndim == 2)
       {
	 if(tmpDelaunay->neighbors[inbr]->has_mesh(id[0]) &&
	     tmpDelaunay->neighbors[inbr]->has_mesh(id[1]) ) 
	 {
	   tmpDelaunay = tmpDelaunay->neighbors[inbr];
	   break;
	 }
       } else if(ndim == 3)
       {
	 if(tmpDelaunay->neighbors[inbr]->has_mesh(id[0]) && tmpDelaunay->neighbors[inbr]->has_mesh(id[1]) 
	     && tmpDelaunay->neighbors[inbr]->has_mesh(id[2]))
	 {
	   tmpDelaunay = tmpDelaunay->neighbors[inbr];
	   break;
	 } 
       }
     }

     if(ndim == 2) rval = tmpDelaunay->crossdeterminant(mesh->getPosx(), mesh->getPosy());
     else if(ndim == 3)
     {
       splittype = tmpDelaunay->getSplitType(mesh->getPosx(), mesh->getPosy(), mesh->getPosz());
       rval = splittype == _nosplit ? 1 : 0;  
     } 

     if(rval == 0)
       break;
     ++iterator; 

     if(iterator == 1000) return _nosplit; 
  }

  it = tmpDelaunay->indd;
  return splittype;

}
