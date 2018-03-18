#ifndef _CHAOS_H_
#define _CHAOS_H_
#include "params.h"
#include "matrix.h"
#include "function.h"
#include "utils.h"
#include "opt.h"

namespace OptimizeName  
{

  class Chaos_Optimization;
  class xConstraint;

  class xConstraint
  {
    private:
      double left;
      double right;
    public:
      xConstraint(){}
      xConstraint(double left, double right)
      {
	this->left = left; this->right = right;
      }
      inline double getleft(){return this->left;}
      inline double getright(){return this->right;}
  };

  class Chaos_Optimization : public OptClass 
  {
    public:
      void initialize();
      void run();
      void free();

      void opt();
      double cooling_dT(const double dT);
      inline double sigmoid(double ui, xConstraint &xc)
      {
	return (xc.getright() + xc.getleft() * exp(-ui)) / (1.0 + exp(-ui));
      }

  };
};

#endif 
