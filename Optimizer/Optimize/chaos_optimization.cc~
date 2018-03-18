#include <iostream>
#include "chaos_optimization.h"

//namespace Chaos_Optimization
namespace OptimizeName 
{

  void Chaos_Optimization::initialize()
  {
  }

  void Chaos_Optimization::run()
  {
    //Function fobj(object, 0);
    initialize();
    //chaos_optimization(fobj);
    opt();
  }

  void Chaos_Optimization::opt()
  {

    Function *pfobj = this->getObjectiveFunction();
    AbstVector u(Nx), v(Nx), x(Nx);
    xConstraint *xc = new xConstraint[Nx];
    const double c = 1.0;
    const double a = 2.0;
    double dT = 0.01;
    
    for(int i = 0; i < Nx; ++i) xc[i] = xConstraint(-5.0, 5.0);	

    x.c[0] = 0.1;

    while(1)
    {

      for(int i = 0; i < Nx; ++i)
      {
	u.c[i] += dT * v.c[i]; 
	v.c[i] += -dT * (a * v.c[i] + c * pfobj->dxi(x.c,i));
	x.c[i] = sigmoid(u.c[i], xc[i]);
      }
      
      cout << u.c[0] << " " << v.c[0] << " " << x.c[0] << endl;
      cout << pfobj->f(x.c) << endl;
      getchar();
    }

    delete [] xc;
    u.clear();
    v.clear();
    x.clear();
  }

  double Chaos_Optimization::cooling_dT(const double dT)
  { 
    return dT;
  }
};
