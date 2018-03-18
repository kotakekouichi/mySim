#ifndef _INTERIOR_POINT_METHOD_
#define _INTERIOR_POINT_METHOD_
#include "opt.h"
#include "params.h"
#include "function.h"
#include "utils.h"
#include "matrix.h"

/* --- class and struct --- */
namespace OptimizeName 
{

  //class Function;
  //class Matrix;
  //class Vector;
#if 0
  class Function : public AbstFunction 
  {
    private:  
      double dx;
    public:

      Function()
      {
	/*
	this->Coeff = new coeff[Nx + 1];

	this->firstCoeff = NULL;

	for(int i = 0; i < Nx + 1; ++i)
	{
	  this->Coeff[i].ix = i;
	  this->Coeff[i].index = i;
	}
	 */
	setdx(eps);
      }
      
      inline void setdx(double val){dx = val;}
      inline double getdx(){return this->dx;}

      double f(double *x)
      {
	int id = this->getID();
	const int constTerm = 1;
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
      Matrix(int rank)
      {
	this->setrank(rank);

	a = new double*[rank];
	for(int i = 0; i < rank; ++i)
	{
	  a[i] = new double[rank]; 
	  for(int j = 0; j < rank; ++j)
	  {
	    a[i][j] = 0.0;
	  }
	}

      }
      Matrix(int rank, double **aij)
      {
	this->setrank(rank);

	a = new double *[rank];
	for(int i = 0;i < rank; ++i) a[i] = new double[rank];

	for(int i = 0; i < rank; ++i)
	{
	  for(int j = 0; j < rank; ++j)
	  {
	    this->a[i][j] = aij[i][j];
	  }
	}

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
      }
      Vector(int rank, double *c)
      {
	this->setrank(rank);
	this->c = new double[rank];
	for(int i = 0; i < rank; ++i)
	  this->c[i] = c[i];
      }

      void setVecIPM(Function &fobj, vector<Function> &g, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec);
  };
#endif 

  class Interior_Point_Method : public OptClass
  {
    private:
      //Function fobj;
      //vector<Function> g;
    public:
      Interior_Point_Method()
      {
        OptClass(); 
      }
      void run();
      void initialize();
      void opt();
      void free();
      void free_Matrix_Vector(Matrix *Mhat, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec, Vector &drvec, Vector &gradvec);
      double getalpha(Vector &drvec, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec);

      void SetObjective(string &str);
      void AddContraint(string &str);
  };

};

#endif

