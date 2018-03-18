#ifndef _PARTICLE_H
#define _PARTICLE_H

//#include "mesh.h"
class meshClass;

enum particle_type 
{
  dark_matter,
  gas,
  star,
  dummy,
};

class particleClass
{
  private:

    particle_type ptype;
    int index;
    double mp, rho;
    double xp, yp, zp;
    double vx, vy, vz;
    double fx, fy, fz;
    double px, py, pz;
    double E, U;
    double p, T;
    meshClass *mgp;

  public:

    particleClass *next;

    inline void setParticleType(particle_type value){ ptype = value;}
    inline particle_type getParticleType(){return ptype;}

    inline void setIndex(int value){index = value;}
    inline int getIndex(){return index;}

    inline void setMass(double value){ mp = value;}
    inline double getMass(){return mp;}

    inline void setDensity(double value){rho = value;}
    inline double getDensity(){return rho;}

    inline void setPosx(double value){xp = value;}
    inline double getPosx(){return xp;}

    inline void setPosy(double value){yp = value;}
    inline double getPosy(){return yp;}

    inline void setPosz(double value){zp = value;}
    inline double getPosz(){return zp;}

    inline void setVelx(double value){vx = value;}
    inline double getVelx(){return vx;}

    inline void setVely(double value){vy = value;}
    inline double getVely(){return vy;}

    inline void setVelz(double value){vz = value;}
    inline double getVelz(){return vz;}

    inline void setFx(double value){fx = value;}
    inline double getFx(){return fx;}

    inline void setFy(double value){fy = value;}
    inline double getFy(){return fy;}

    inline void setFz(double value){fz = value;}
    inline double getFz(){return fz;}
 
    inline void setE(double value){E = value;}
    inline double getE(){return E;}

    inline void setU(double value){U = value;}
    inline double getU(double value){return U;}

    inline void setPressure(double value){p = value;}
    inline double getPressure(){return p;}

    inline void setTemperature(double value){T = value;}
    inline double getTemperature(){return T;}

    inline void setPx(double value){px = value;}
    inline double getPx(){return px;}

    inline void setPy(double value){py = value;}
    inline double getPy(){return py;}

    inline void setPz(double value){pz = value;}
    inline double getPz(){return pz;}

    inline void setMeshGenerationPoint(meshClass *value){mgp = value;}
    inline meshClass* getMeshGenerationPoint(){return mgp;}
};


#endif 

