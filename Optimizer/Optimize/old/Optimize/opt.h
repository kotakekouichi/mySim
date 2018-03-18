#ifndef _OPT_H_
#define _OPT_H_
#include "params.h"
#include "function.h"
#include <vector>

using namespace std;

namespace OptimizeName
{
  const double eps = 1.0e-12;

  class Function : public AbstFunction
  {
    private:

    public:
      Function()
      {
	try
	{
	  this->Coeff = new coeff[Nx+1];
	  this->firstCoeff = NULL;

	  for(int i = 0; i < Nx; ++i)
	  {
	    this->Coeff[i].ix = i;
	    this->Coeff[i].index = i;
	    //this->Coeff[i].val = 0.0;
	    //this->Coeff[i].preCoeff = NULL;
	    //this->Coeff[i].nextCoeff = NULL;
	  }

	}
	catch(std::bad_alloc)
	{
	  cout << "a" << endl;
	}
      }

      void SwapVariable(int xpivod, int ix,Function &Bi);
      int getXpivodCol(vector<Function> &B);
      int getXpivodRow(int &xpivod, vector<Function> &B);
      void migration(int xpivodcol);

  };
};

#endif 
