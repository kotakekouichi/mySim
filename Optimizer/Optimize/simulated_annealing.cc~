#include <iostream>
#include "simulated_annealing.h"

using namespace std;

namespace Simulated_Annealing
{

  void run()
  {
    Function func;

    initialize(func);

    SA(func);

    free(func);
  }

  void initialize(Function &func)
  {
    for(int i = 0; i < Nx;++i)
    {
      func.Coeffval[i] = (double) i;
    }
  }

  void free(Function &func)
  {
    delete [] func.Coeffval;
  }

  void SA(Function &func)
  {

    cout << "start sa" << endl;

    cout << "end sa" << endl;
  }

};
