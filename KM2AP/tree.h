#ifndef _TREE_H
#define _TREE_H
#include "mesh.h"
#include "particle.h"

template<typename type>
class tree
{
  private :
    int num;
    int dim;
    int nth;
    double xc, yc, zc;
    double xg, yg, zg;
    double mass;
    double dx;
    tree* nextlevel;
  public:
    tree(){}
    tree(int _num, double _xc, double _yc, double _dx, int _nth)
    {
      this->num = _num;
      this->nth = _nth;
      this->xc = _xc;
      this->yc = _yc;
      this->dx = _dx;
      this->xg = this->yg = 0.0;
      this->mass = 0.0;
      this->nextlevel = NULL;
      this->first = NULL;
    }
    tree(int _num, double _xc, double _yc, double _zc, double _dx, int _nth)
    {
      this->num = _num;
      this->nth = _nth;
      this->xc = _xc;
      this->yc = _yc;
      this->zc = _zc;
      this->dx = _dx;
      this->xg = this->yg = this->zg =0.0;
      this->mass = 0.0;
      this->nextlevel = NULL;
      this->first = NULL;
    }

    type *first;

    inline void setdx(double value){dx = value;}
    inline double getdx(){return dx;}

    inline void setxc(double value){xc = value;}
    inline double getxc(){return xc;}

    inline void setyc(double value){yc = value;}
    inline double getyc(){return yc;}

    inline void setzc(double value){zc = value;}
    inline double getzc(){return zc;}

    inline void setdim(double value){dim = value;}
    inline int getdim(){return dim;}

    int initialize(type *point);
    int construct();
    int build_tree();
    int delaunayfinder(int imesh, type *Mesh, splitType &splitype);
    int walk_tree(int imesh, meshClass *Mesh);
    int free();
    int gravity_force(particleClass *part);
    int gravity_shortrange_force(particleClass *part);
};
#endif
