#ifndef _PM_H_
#define _PM_H_
#include "particle.h"


class pmMeshClass
{
  private:
    int index;
    double dx;
    double xc, yc, zc;
    double rho, phix;
    double ax, ay, az;
   
  public:

    pmMeshClass(){}

    inline void setindex(int value){index = value;}
    inline int getindex(){return this->index;}

    inline void setdx(double value){dx = value;}
    inline double getdx(){return this->dx;}

    inline void setxc(double value){xc = value;}
    inline double getxc(){return this->xc;}
    
    inline void setyc(double value){yc = value;}
    inline double getyc(){return this->yc;}
    
    inline void setzc(double value){zc = value;}
    inline double getzc(){return this->zc;}

    inline void setrho(double value){rho = value;}
    inline double getrho(){return this->rho;}

    inline void setphix(double value){phix = value;}
    inline double getphix(){return this->phix;}
    
    inline void setax(double value){ax = value;}
    inline double getax(){return this->ax;}
    
    inline void setay(double value){ay = value;}
    inline double getay(){return this->ay;}
    
    inline void setaz(double value){az = value;}
    inline double getaz(){return this->az;}

};
#endif
