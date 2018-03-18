#ifndef _MESH_H
#define _MESH_H
#include "particle.h"
#include "params.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

enum splitType
{
  _1to4,
  _2to6,
  _nto2n,
  _nosplit,
};

enum flipType
{
  _2to3,
  _3to2,
  _4to4,
  _noflip,
};

class particleClass;
class meshClass;
class delaunayClass;
class faceClass;

class delaunayClass 
{
  private:
    double xp[4];
    double yp[4];
    double zp[4];

    double xc, yc, zc, radius, V;
    
  public:

    delaunayClass(){this->nNeighbor = 0; flag = true;checkFlag=false;}
    ~delaunayClass(){}
    //delaunayClass(int index, int imesh, int jmesh, int kmesh, meshClass *Mesh);

    string birthtype;
    int indd;
    int indm[4];
    bool flag; //trueの場合、violent faceのチャックを行う
    bool checkFlag;
    double orient;

    //vector<delaunayClass*> neighbors;
    delaunayClass *neighbors[4];
    faceClass *faces[4];
    int nNeighbor;
    double sss, ttt, uuu, vvv;
    //map<int, int> neighborsmap;

    void set(int index, int imesh, int jmesh, int kmesh, meshClass *Mesh);
    void set(int index, int imesh, int jmesh, int kmesh, int lmesh, meshClass *Mesh);

    inline void clone(delaunayClass &baseDel)
    {
       indm[0] = baseDel.indm[0];
       indm[1] = baseDel.indm[1];
       indm[2] = baseDel.indm[2];
       indm[3] = baseDel.indm[3];

       nNeighbor  = baseDel.nNeighbor; 

       indd = baseDel.indd;

       for(int i = 0; i < nNeighbor;++i)
       {
	 neighbors[i] = baseDel.neighbors[i];
       }

       for(int i = 0; i < 4; ++i)
       {
	 if(baseDel.faces[i] == NULL)
	   faces[i] = NULL;
	 else 
	   faces[i] = baseDel.faces[i];
       }

       birthtype = baseDel.birthtype;

    }

    inline void setPosx(int i, double value){xp[i] = value;}
    inline double getPosx(int i){return xp[i];}

    inline void setPosy(int i, double value){yp[i] = value;}
    inline double getPosy(int i){return yp[i];}

    inline void setPosz(int i, double value){zp[i] = value;}
    inline double getPosz(int i){return zp[i];}

    inline void setxc(double value){xc = value;}
    inline double getxc(){return xc;}

    inline void setyc(double value){yc = value;}
    inline double getyc(){return yc;}

    inline void setzc(double value){zc = value;}
    inline double getzc(){return zc;}

    inline void setRadius(double value){radius = value;}
    inline double getRadius(){return radius;}

    inline bool has_mesh(const int imesh)
    {
      return (imesh == this->indm[0] || imesh == this->indm[1] || imesh == this->indm[2] || imesh == this->indm[3]); 
    }

    int crossdeterminant(double xp, double yp);
    splitType getSplitType(double xp, double yp, double zp);
    flipType getFlipType(int imesh, int jmesh, meshClass *Mesh, bool flag, int *no_indm);
    int determinant2d(double xp, double yp);
    int determinant3d(double xp, double yp, double zp);
    int circumcenter(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double *xc, double *yc, double *radius);
    int circumcenter3d();
    
    int changeToNewDelaunay();
    int eraseDelaunaryOfMesh(meshClass *Mesh, paramsClass &params);
    int addDelaunayToMesh(meshClass *Mesh, paramsClass &params);

    bool skipflipFlag(double xp, double yp, double zp);

    splitType getSplitType_WalkNeighbor(int &it, meshClass *mesh, paramsClass &params);

    void dump()
    {
      cout << "n size = " << this->nNeighbor << " " << this->indm[0] << "  " << this->indm[1] << " " << this->indm[2]  << " " << this->indm[3]  << endl;
      //cout << this->getPosx(0) << " " << this->getPosy(0) <<  " " << getPosz(0 ) << endl;
      //cout << this->getPosx(1) << " " << this->getPosy(1) << " " << getPosz(1) << endl;
      //cout << this->getPosx(2) << " " << this->getPosy(2) << " " << getPosz(2) <<  endl;
      //cout << this->getPosx(3) << " " << this->getPosy(3) << " " << getPosz(3) <<  endl;
      
      cout << endl;

      for(int inbr = 0; inbr < this->nNeighbor; ++inbr)
      {
        cout << this->neighbors[inbr]->indd << " /   " << this->neighbors[inbr]->indm[0] << ", " 
	              << this->neighbors[inbr]->indm[1] << ", " 
		      << this->neighbors[inbr]->indm[2] << ", " 
		      << this->neighbors[inbr]->indm[3] << endl; 
      }
    }
    void check()
    {
      if(this->nNeighbor == 3) return;
      
      int count = 0;
      for(int inbr = 0; inbr < this->nNeighbor; ++inbr)
      {
	for(int i = 0; i < 4; ++i)
	  if(this->has_mesh(this->neighbors[inbr]->indm[i])) 
	    count++;
	if(count != 3){this->dump();getchar(); break;};
	count = 0;
      }

    }
    
    bool operator < (const delaunayClass &right) const 
    {
      return 1.0 / radius < 1.0 / right.radius;
    }


};

class faceClass
{
  private:
    double x;
    double y;
    double z;
  public:
    delaunayClass *interDel[2];
    int interidx[2];
    int inbrs[2];
    faceClass();
    faceClass(double _x, double _y, double _z)
    {
      x = _x;
      y = _y;
      z = _z;
    }
};

class meshClass
{
  
  private:
    particleClass *_particle;
    double xg, yg, zg;

  public:

    meshClass(){delaunay.reserve(100); neighbors.reserve(100);}
    int index;
    int npoly;
    double V;
    double S, f[3];
    double Aij;
    double dQi[5]; //dUi;

    meshClass *next;

    vector<delaunayClass*> delaunay; 
    std::unordered_map<int, delaunayClass*> delaunaymap;
    vector<meshClass> neighbors;
    
    inline void setParticle(particleClass *value){_particle = value;}
    inline particleClass* getParticle(){return _particle;}
    
    inline double getMass(){return _particle->getMass();}
    inline double getPosx(){return _particle->getPosx();}
    inline double getPosy(){return _particle->getPosy();}
    inline double getPosz(){return _particle->getPosz();}
    inline double getE(){return _particle->getE();}
    inline double getPx(){return _particle->getPx();}
    inline double getPy(){return _particle->getPy();}
    inline double getPz(){return _particle->getPz();}

    inline void setMass(double value){_particle->setMass(value);}
    inline void setPosx(double value){_particle->setPosx(value);}
    inline void setPx(double value){_particle->setPx(value);}
    inline void setPosy(double value){_particle->setPosy(value);}
    inline void setPy(double value){_particle->setPy(value);}
    inline void setPosz(double value){_particle->setPosz(value);}
    inline void setPz(double value){_particle->setPz(value);}

    inline void setE(double value){_particle->setE(value);}

    inline double pressure()
    {
      return ((5.0 /3.0 - 1.0) * (getE() - 0.5 * (getPx() * getPx() + getPy() * getPy() + getPz() * getPz()) / getMass()) / V); 
    }


    inline double getxg(){return this->xg;}
    inline double getyg(){return this->yg;}
    inline double getzg(){return this->zg;}

    void calcS2d();
    void calcV3d();
    void calcGravityCenter();
};

class data_tClass
{
  public: 
    int index;
    double theta;
    double l;
    double xp, yp, zp;
    
    bool operator < (const data_tClass &right) const 
    {
      return theta < right.theta;
    }

    void compute_theta(double *r1, double *r2, double *nvec)
    {
      double dot = r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
      double l1 = sqrt(sq(r1[0]) + sq(r1[1]) + sq(r1[2]));
      double l2 = sqrt(sq(r2[0]) + sq(r2[1]) + sq(r2[2])); 

      double lx = r1[1] * r2[2] - r1[2] * r2[1];
      double ly = r1[2] * r2[0] - r1[0] * r2[2];
      double lz = r1[0] * r2[1] - r1[1] * r2[0];
      
      double cos = dot / l1 / l2;

      if(cos > 1.0) cos = 1.0;
      else if(cos < -1.0) cos = -1.0;

      double sin = lx * nvec[0] + ly * nvec[1] + lz * nvec[2];
      if(fabs(sin) < 1.0e-10) sin = 0.0;
      theta = sin >= 0.0 ? acos(cos) : 2 * M_PI - acos(cos);
      l = l2;
      
    }

};


#endif
