#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_
#include "opt.h"
#include "params.h"
#include "function.h"

/* --- class and struct --- */

namespace OptimizeName 
{
  class Simplex : public OptClass
  {
    private:
    public:

      Simplex()
      {
	OptClass();
      } 
      Simplex(Function *pfobj, vector<Function> *constfunc)
      {
	fobj = pfobj;
	g = *constfunc;
      }
      
      void run();
      void initialize();
      void free();
      void opt();
      void optRevised();
      void getFeasibleInitilize(int *basicVariableIdx, int *nonbasicVariableIdx);
      void setMatrixes(int *basicVariableIdx, int *nonbasicVariableIdx, Matrix *A, Matrix *B, Matrix *N);
      int setArtificialVariable(Function *pfobj);
      void setConstraintVector(Vector *tmpVec);
      void setConstraintMatrix(Matrix *A);

      int swapVariableRevised(int nonbasicidx, Vector *y, Vector *xB, Vector *xN, Vector *bvar, int *basicVariableIdx, int *nonbasicVariableIdx, Vector *cN, Vector *cB, Matrix *B, Matrix *N);  
      int swapVariableRevised_r2(int Bidx, int Nidx, int *basicVariableIdx, int *nonbasicVariableIdx, Vector *cB, Vector *cN, Matrix *B, Matrix *N, AbstMatrix *Bt0, AbstMatrix *Nt0);

      int getNonBasicIdx(Vector *cN, AbstMatrix *N, Vector *pi);

      void SetObjective(string &strFunc);
      void AddConstraint(string &strConstFunc);
      
      void FuncLog();

      void checkBasicParams(Function *fobj, int &artitermIdx);

      //==============
      int getzsNindex(AbstVector &vec);
      int getbasicindex(double &dbl, AbstVector &dx, AbstVector &b);
      //=============
  };
};

#endif
