#ifndef _PARAMS_H
#define _PARAMS_H

using namespace std;

extern double a;
extern double t;

enum var_type
{
  INT,
  DOUBLE,
  STRING,
};

struct adr_type
{
  void *var_adr;
  enum var_type type;
};

class paramsClass
{
  public:

    int periodic;
    int selfgravity;
    int dim;
    int Ndm;
    int Ngas;
    int Nghost;
    int Ndelaunay;
    int NPMGRID;
    int cosmological;
    int Ngrid;
    string dmfile;
    string gasfile; 
    double OmegaBaryon, OmegaCDM, OmegaMatter, OmegaLambda;
    double Hubble0;

    paramsClass()
    {
       periodic = 0;
       selfgravity = 0;
       dim = 0;
       Ndm = 0;
       Ngas = 0;
       Nghost = 0;
       Ndelaunay = 0;
       NPMGRID = 0;
       cosmological = 0;
       dmfile = "";
       OmegaLambda = OmegaCDM = OmegaMatter = OmegaLambda = 0.0;
       Hubble0 = 0.0;
    }
};

#endif
