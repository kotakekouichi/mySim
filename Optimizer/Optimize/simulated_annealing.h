#ifndef _SA_H_
#define _SA_H_
#include "params.h"

namespace Simulated_Annealing
{

  class Function
  {
    public:
      double *Coeffval;
      Function()
      {
	this->Coeffval = new double[Nx];
      }
  };

  void run();
  void SA(Function &func); 
  void initialize(Function &func);
  void free(Function &func);
};

#endif
