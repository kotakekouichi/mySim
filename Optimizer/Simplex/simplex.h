#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_
#include "params.h"

/* --- class and strunct --- */
namespace Simplex
{
  class coeff 
  {
    public:
      int index;
      int ix;
      coeff *nextCoeff; 
      coeff *preCoeff;
      double val;
      coeff()
      {
	this->val = 0.0;
      }
  };

  class Function
  {
    public:
      int ix;
      coeff *firstCoeff;
      coeff *Coeff;

      Function()
      {
	this->Coeff = new coeff[Nx];
	
      }

      int SwapVariable(int xpivod, int ix,Function Bi);
      void CheckFunc();
  };
  
  void run();
  void initialize(Function &fobj, Function *B);
  void free(Function &fobj, Function *B);
  void simplex(Function &fobj, Function *B);
  
};

#endif
