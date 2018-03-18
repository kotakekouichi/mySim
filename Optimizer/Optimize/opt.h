#ifndef _OPT_H_
#define _OPT_H_
#include "params.h"
#include "function.h"
#include "matrix.h"
#include "utils.h"

using namespace std;

namespace OptimizeName
{

  class OptClass;
  class Function;
  class Matrix;
  class Vector;
  const double eps = 1.0e-12;

  class Function : public AbstFunction
  {
    private:
      double dx;

    public:
      Function()
      {
	setdx(eps);
      }

      Function(efunc type, int id)
      {
	setdx(eps);
	setfunctype(type, id); // 関数の種類とIDを設定
	//AbstFunction(type, id);
      }


      int SwapVariable(int xpivod, int ix,Function &Bi);
      int getXpivodCol(vector<Function> &B);
      int getXpivodRow(int &xpivod, vector<Function> &B);
      void migration(int xpivodcol);
      void CheckFunc();

      inline void setdx(double val){dx = val;}
      inline double getdx(){return dx;}

      double f(double *x)
      {
	int id = this->getID();
	//const int constTerm = 1;
	efunc ftype = this->getfunctype();
	if(ftype == object)
	{
	  //return (x[0]*x[0]*x[0]);
	  //return sin(x[0]);
	  //return sq(x[0] / (sq(x[0]) + 1));
	  //return sq(x[0]) + sq(x[1]) + x[0] * x[1];
	  //return exp(x[0] * x[1]);
	  return -(x[0] + 2*x[1] + x[2] + x[3]); // + => min / - => max
	}
	else if(ftype == constraint && id == 0)
	{
	  /*for(idx = 0; idx < Nx;++idx)
	    {
	    val Coeff[idx + constterm].val * x[idx];
	    }
	   */
	  return -3.0 + 1.0 * x[0] - 2.0 * x[1] + 2.0 * x[2] + 3.0 * x[3];
	}
	else if(ftype == constraint && id == 1)
	{
	  return -6.0 + 3.0 * x[0] + 4.0 * x[1] - 5.0 * x[2] - 6.0 * x[3];
	}
	else if(ftype == constraint && id == 2)
	{
	  return -2.0 + 3.0 * x[0] + 1.0 * x[1] - 5.0 * x[2] - 7.0 * x[3];
	}
	else if(ftype == constraint && id == 3)
	{
	  return -2.0 - 2.0 * x[0] + 1.0 * x[1] + 2.0 * x[2] + 1.0 * x[3];
	}
	else if(ftype == constraint && id == 4)
	{
	  return -5.0 + 3.0 * x[0] - 3.0 * x[1] + 2.0 * x[2] + 2.0 * x[3];
	}
	else if(ftype == constraint && id == 5)
	{
	  return 0.0 - 1.0 * x[0];
	}
	else if(ftype == constraint && id == 6)
	{
	  return 0.0 - 1.0 * x[1];
	}
	else if(ftype == constraint && id == 7)
	{
	  return 0.0 - 1.0 * x[2];
	}
	else if(ftype == constraint && id == 8)
	{
	  return 0.0 - 1.0 * x[3];
	}

	return 1.0e30;
      }

      double dxi(double *x, const int i)
      {
	double dfdxi = 0.0;
	x[i] = x[i] + dx;
	dfdxi = f(x);
	x[i] = x[i] - 2 * dx;
	dfdxi = dfdxi - f(x);
	x[i] = x[i] + dx;

	return dfdxi / (2*dx);
      }

      double dxidxi(double *x, const int i)
      {
	double ddfd2xi = 0.0;
	x[i] = x[i] + dx;
	ddfd2xi = f(x);
	x[i] = x[i] - dx;
	ddfd2xi = ddfd2xi - 2 * f(x);
	x[i] = x[i] - dx;
	ddfd2xi = ddfd2xi + f(x);
	x[i] = x[i] + dx;
	ddfd2xi = fabs(ddfd2xi) < 1.0e-10 ? 0.0 : ddfd2xi / sq(dx);

	return ddfd2xi;
      }
      double dxidxj(double *x, const int i, const int j)
      {
	double dxdy = 0.0;
	double xi0 = x[i], xj0 = x[j];
	double f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0, f5 = 0.0, f6 = 0.0, f7 = 0.0;

	if(i == j) return this->dxidxi(x, i);

	x[i] = xi0 + dx;
	x[j] = xj0 + dx;
	f1 = f(x);	

	x[i] = xi0 - dx; 
	x[j] = xj0 - dx; 
	f2 = f(x);

	x[i] = xi0 + dx;
	x[j] = xj0;
	f3 = f(x);

	x[i] = xi0 - dx;
	x[j] = xj0;
	f4 = f(x);

	x[i] = xi0;
	x[j] = xj0 + dx;
	f5 = f(x);

	x[i] = xi0;
	x[j] = xj0 - dx;
	f6 = f(x);

	x[i] = xi0;
	x[j] = xj0;
	f7 = f(x);

	dxdy = fabs(f1 + f2 - f3 - f4 - f5 - f6 + 2.0 * f7 ) / (2.0 * sq(dx))< 1.0e-10 ?
	  0.0 : (f1 + f2 - f3 - f4 - f5 - f6 + 2.0 * f7) / (2.0 * sq(dx));

	return dxdy;
      }
  };

  class Matrix : public AbstMatrix
  {
    public :

      Matrix(){}
      Matrix(int Rank)
      {
	setrank(Rank);
	InitMatrix(Rank, Rank, NULL);
      }
      Matrix(int Rank, double **aij)
      {
	setrank(Rank);
	InitMatrix(Rank, Rank, aij);
      }
      Matrix(int nrow, int ncol, double **aij)
      {
	setrowrank(nrow);
	setcolrank(ncol);
	InitMatrix(nrow, ncol, aij);
      }
      Matrix(int nrow, int ncol)
      {
	setrowrank(nrow);
	setcolrank(ncol);
	InitMatrix(nrow, ncol, NULL);
      }

      void setMatrixIPM(Function &fobj, vector<Function> &g, Vector &xvec, Vector &zvec, Vector &uvec,Vector &lambdavec);
  };

  class Vector : public AbstVector
  {
    public:

      Vector(){}
      Vector(int rank)
      {
	this->setrank(rank);
	this->c  = new double[rank];
	for(int i = 0; i < rank; ++i)
	  this->c[i] = 0.0;
	this->L = NULL;
	this->U = NULL;
      }
      Vector(int rank, double *c)
      {
	this->setrank(rank);
	this->c = new double[rank];
	for(int i = 0; i < rank; ++i)
	  this->c[i] = c[i];
	this->L = NULL;
	this->U = NULL;
      }

      void setVecIPM(Function &fobj, vector<Function> &g, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec);
  };

  class OptClass
  {
    protected:
      Function *fobj;
      vector<Function> g; 
    public:

      OptClass(){this->fobj = new Function();g.reserve(10000000);}

      Function *getObjectiveFunction(){return fobj;}
      vector<Function> *getConstraintFunction(){return &g;} 
      
      virtual void initialize();
      virtual void opt();
      virtual void clear();
      virtual void SetObjective(string &str) ;
      virtual void AddConstraint(string &str) ;

      virtual ~OptClass(){}
  };

};

#endif 
