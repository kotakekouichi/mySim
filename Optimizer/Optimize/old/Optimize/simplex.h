#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_
#include "opt.h"
#include "params.h"
#include "function.h"
#include <vector>

using namespace std;

/* --- class and struct --- */

namespace OptimizeName 
{
  class Simplex 
  {
    private:
      Function *fobj;
      vector<Function> g; 
    public:

      Simplex()
      {
	{this->fobj = new Function();g.reserve(10000000);}
      } 
      Function *getObjectiveFunction(){return fobj;}
      vector<Function> *getConstraintFunction(){return &g;} 
     
      void initialize();
      void free();
      void opt();
      void getFeasibleInitilize();
      void setArtificialVariable(Function *pfobj);
  };
};

#endif
