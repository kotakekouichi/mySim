#ifndef _MATRIX_H
#define _MATRIX_H
#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#define IDX(i, j) ((j) + (i) * ncol)
using namespace std;

class AbstMatrix;
class AbstVector;
class ComponentClass;

enum  enmSimultaneousLinearEquations
{
  enmGaussJordan,
  enmLU,
  enmIterationMethod,
};

enum enmMatrixType
{
  normal,
  unit,
  sparce,
};

class ComponentClass 
{
  private:
  public:
    int i, j;
    double val;
    ComponentClass(){}
    ComponentClass(int _i, int _j, int _val)
    {
      this->i = _i;
      this->j = _j;
      this->val = _val; 
    }
};

class AbstMatrix
{
  private:
    int rank;
    int nrow;
    int ncol;

    enmMatrixType matrixType;
  
  public:
    double **a;
    double *a2;
    double det;

    int *preNonZeroColIndex;
    int *preNonZeroRowIndex;
    int *nextNonZeroColIndex;
    int *nextNonZeroRowIndex;

    int *firstNonZeroRowIndex;
    int *lastNonZeroRowIndex;

    map<int, double> compmap;
    map<int, ComponentClass> compmap2;

    AbstMatrix *Inv;
    AbstMatrix *_Trace;

    AbstMatrix(){}
    AbstMatrix(int Rank)
    {
      this->rank = Rank;
      //InitMatrix(rank, rank, NULL);
      InitMatrix(rank, rank);
    }

    AbstMatrix(int nRow, int nCol)
    {
      this->nrow = nRow;
      this->ncol = nCol; 
      this->rank = nCol > nRow ? nCol : nRow;
      //InitMatrix(nrow, ncol, NULL);
      InitMatrix(nrow, ncol);
    }

    AbstMatrix(int Rank, double **aij)
    {
      this->rank = Rank;
      InitMatrix(rank, rank, aij);
    }

    AbstMatrix(int nRow, int nCol, double **aij)
    {
      this->nrow = nRow;
      this->ncol = nCol;
      this->rank = nCol > nRow ? nCol : nRow;
      InitMatrix(nrow, ncol, aij);
    }


    ~AbstMatrix()
    {
    }

    inline void setMatrixType(enmMatrixType value){this->matrixType = value;}
    inline enmMatrixType getMatrixType(){return this->matrixType;}
   
    void InitMatrix(int nRow, int nCol)
    {
      this->nrow = nRow;
      this->ncol = nCol;
      
      a2 = new double [nRow * nCol];
      preNonZeroColIndex = new int[nrow*ncol];
      preNonZeroRowIndex = new int[nrow*ncol];
      nextNonZeroColIndex = new int[nrow*ncol];
      nextNonZeroRowIndex = new int[nrow*ncol];
     
      firstNonZeroRowIndex = new int[ncol];
      lastNonZeroRowIndex = new int[ncol];
/*
      for(int i = 0; i < nrow; ++i)
	for(int j = 0;j < ncol; ++j)
	  this->setComponent(i, j, 0.0);

      for(int i = 0; i < nrow * ncol; ++i)
      {
	preNonZeroRowIndex[i] = -1;
	preNonZeroColIndex[i] = -1;
	nextNonZeroColIndex[i] = -1;
	nextNonZeroRowIndex[i] = -1;
      }

*/

      this->Inv = NULL;
      this->_Trace = NULL;

    }

    void InitMatrix(int nRow, int nCol, double **aij)
    {
      this->nrow = nRow;
      this->ncol = nCol;
      //a = new double*[nrow];
      a2 = new double [nRow * nCol];
      preNonZeroColIndex = new int[nrow*ncol];
      preNonZeroRowIndex = new int[nrow*ncol];
      nextNonZeroColIndex = new int[nrow*ncol];
      nextNonZeroRowIndex = new int[nrow*ncol];

      firstNonZeroRowIndex = new int[ncol];
      lastNonZeroRowIndex = new int[ncol];

      //for(int i = 0; i < nrow; ++i) a[i] = new double[ncol];

      if(aij == NULL)
      {
//	for(int i = 0; i < nrow; ++i)
//	  for(int j = 0; j < ncol; ++j)
//	    this->a[i][j] = 0.0;

	for(int i = 0; i < nrow; ++i)
	  for(int j = 0;j < ncol; ++j)
	    this->setComponent(i, j, 0.0);
      }
      else 
      {

	//for(int i = 0; i < nrow; ++i)
	 // for(int j = 0; j < ncol; ++j)
	  //  this->a[i][j] = aij[i][j];
	
	for(int i = 0;i < nrow; ++i)
	  for(int j = 0; j < ncol; ++j)
	    this->setComponent(i, j, aij[i][j]);
      }
      
      for(int i = 0; i < nrow*ncol; ++i)
      {
	preNonZeroRowIndex[i] = -1;
	preNonZeroColIndex[i] = -1;
	nextNonZeroColIndex[i] = -1;
	nextNonZeroRowIndex[i] = -1;
      }

      this->Inv = NULL;
      this->_Trace = NULL;

    }

    void clear()
    {  
      /*
      const int nRow = this->getrowrank();
      for(int i = 0; i < nRow; ++i)
      {	
	delete [] this->a[i];
      }
      */
      //delete [] this->a;
      
      delete [] this->a2;
      delete [] preNonZeroRowIndex;
      delete [] preNonZeroColIndex;
      delete [] nextNonZeroColIndex;
      delete [] nextNonZeroRowIndex;
      delete [] firstNonZeroRowIndex;
      delete [] lastNonZeroRowIndex;

      if(this->Inv != NULL){ this->Inv->clear(); delete this->Inv;}
      if(this->_Trace != NULL){this->_Trace->clear(); delete this->_Trace;}
    }


    inline void setrank(int rank){this->rank = rank;}
    inline int getrank(){return this->rank;}

    inline void setrowrank(int rank){this->nrow = rank;}
    inline int getrowrank(){return this->nrow;}

    inline void setcolrank(int rank){this->ncol = rank;}
    inline int getcolrank(){return this->ncol;}

    //inline void setComponent(int i, int j, double val){this->a[i][j] = val;}
    //inline double getComponent(int i, int j){return this->a[i][j];}
    inline void setComponent(int i , int j, double val)
    {
      int index = IDX(i, j);
      this->a2[index] = val;
      //if(val != 0.0)
      if(fabs(val) > 1.0e-16)
      {
	//this->compmap[index] = val;
	ComponentClass val2(i, j, val);
	this->compmap2[index] = val2;
      }
      //else if(this->compmap.find(index) != this->compmap.end())
      //{
//	this->erase(i, j);
     // }
      else if(this->compmap2.find(index) != this->compmap2.end())
      {
        this->erase(i, j);
      }
    }
    inline double getComponent(int i, int j ) 
    {
      int index = IDX(i, j);
      return this->a2[index];
      //return this->compmap2.find(index) != this->compmap2.end() ? this->compmap2[index].val : 0.0;
      //return this->compmap2[index].val;
    }

    inline void erase(int i, int j)
    {
      int index = IDX(i, j);
      //map<int, double>::iterator ite = this->compmap.find(index);
      //if(ite != this->compmap.end()) 
//	this->compmap.erase(ite);

      map<int, ComponentClass>::iterator ite2 = this->compmap2.find(index);
      if(ite2 != this->compmap2.end())
	this->compmap2.erase(ite2);
      return;
    }

    inline double getdet(){return this->det;}
    void calcdet()
    {
      const double eps = 0.0;
      double buf = 0.0;
      int Rank = this->getrank();
      //AbstMatrix Ahat(Rank, this->a);
      AbstMatrix Ahat(Rank); 

      for(int i = 0; i < Rank; ++i)
      {
        for(int j = 0; j < Rank; ++j)
	{
	  Ahat.setComponent(i, j, this->getComponent(i, j));
	}	  
      }

      for(int i = 0; i < Rank; ++i)
	for(int j = 0; j < Rank; ++j)
	  if(i < j) 
	  {
	    if(Ahat.getComponent(i, i) == 0.0)
	    {
	      int pivodrow = Ahat.pivodOperator(i, i, eps);
	      if(pivodrow >= 0)
		Ahat.swapRow(i, pivodrow);
	      else 
		return;
	    }
	    buf = Ahat.getComponent(j, i) / Ahat.getComponent(i, i);

	    for(int k = 0; k < Rank; ++k) Ahat.setComponent(j, k, Ahat.getComponent(j, k) - Ahat.getComponent(i, k) * buf);
	  }

      for(int i = 0; i < Rank; ++i) det *= Ahat.getComponent(i, i);

      Ahat.clear();
    } 

    void calcInverseMatrix()
    {
      double buf = 0.0;
      const double eps = 1.0e-10;
      int Rank = this->getrank();

      //AbstMatrix Mhat(Rank, this->a);
      AbstMatrix Mhat(Rank);

      for(int i = 0; i < Rank; ++i)
      {
        for(int j = 0; j < Rank; ++j)
	{
	  Mhat.setComponent(i, j, this->getComponent(i, j));
	}
      }

      AbstMatrix *Ahat = &Mhat;
      AbstMatrix *Bhat = Inv;

      for(int i = 0; i < Rank; ++i)
      {
	for(int j = 0; j < Rank; ++j)
	{
	  Bhat->setComponent(i, j, i == j ? 1.0 : 0.0);
	}
      }

      for(int i = 0; i < Rank; ++i)
      {
	if(fabs(Ahat->getComponent(i, i)) < eps)
	{
	  int pivodrow = Ahat->pivodOperator(i, i, eps);
	  if(pivodrow >= 0)
	  {
	    Ahat->swapRow(i, pivodrow);
	    Bhat->swapRow(i, pivodrow);
	  }
	  else
	  {
	    this->Inv = NULL;
	    cout << i << " " << Rank << endl;
	    return;
	  }
	}

	buf = 1.0 / Ahat->getComponent(i, i);

	for(int j = 0; j < Rank;++j)
	{
	  Ahat->setComponent(i, j, Ahat->getComponent(i, j) * buf);
	  Bhat->setComponent(i, j, Bhat->getComponent(i, j) * buf);
	}

	for(int j = 0; j < Rank; ++j)
	{
	  if( i != j )
	  {
	    buf = Ahat->getComponent(j, i);
	    for(int k = 0; k < Rank; ++k)
	    {
	      Ahat->setComponent(j, k, Ahat->getComponent(j, k) - Ahat->getComponent(i, k) * buf);
	      Bhat->setComponent(j, k, Bhat->getComponent(j, k) - Bhat->getComponent(i, k) * buf); 
	    }
	  }
	}
      }

      Ahat->clear();

    }

    inline void AllocInvMatrix()
    {
      const int Rank = this->getrank();
      this->Inv = new AbstMatrix(Rank);
    }

    inline void AllocTraceMatrix()
    {
      const int nRow = this->getcolrank();
      const int nCol = this->getrowrank();
      this->_Trace = new AbstMatrix(nRow, nCol);
    } 

    AbstMatrix* getInverseMatrix()
    {
      return this->Inv;
    }

    AbstMatrix* getTraceMatrix()
    {
      const int nRow = this->getrowrank();
      const int nCol = this->getcolrank();

      for(int i = 0; i < nRow; ++i)
      {
	for(int j = 0; j < nCol; ++j)
	{
	  this->_Trace->setComponent(j, i, this->getComponent(i, j));
	}
      }

      return this->_Trace;
    }
    
    int pivodOperator(int i, int j, const double eps)
    {
      int Rank = this->getrank();
      int pivod = -1;
      double val = 0.0;
      double eps0 = eps;

      //for(int k = i + 1; k < Rank; ++k)
      for(int k = this->getFirstRowIndex(i + 1, j); k != -1; k = this->getNextNonZeroRowIndex(k, j))
      {
	if(fabs(this->getComponent(k, j)) > eps * 0)
	{
	  if(val < fabs(this->getComponent(k, j)))
	  {
	    pivod = k;
	    val = fabs(this->getComponent(k, j));
	  }
	}
      }	

      if(pivod >= 0) return pivod;

      while(eps0 > 1.0e-16)
      {
	eps0 /= 2;
	for(int k = i + 1; k < Rank; ++k)
	{
	  if(fabs(this->getComponent(k, j)) > eps0)
	  {
	    if(val < fabs(this->getComponent(k, j)))
	    {
	      pivod = k;
	      val = fabs(this->getComponent(k, j));
	    }
	  }
	}
      }

      return pivod;
    }
    
    int pivodOperatorLU(int i, int j, const double eps, AbstMatrix *L, AbstMatrix *U)
    {
#if 1
      int pivod = -1;
      int Rank = this->getrank();
      double dblval = 0.0;
      double sum = 0.0;

      for(int l = j; l < Rank; ++l)
      {
	sum = 0.0;
	for(int k = 0; k < l; ++k)
	  sum += L->getComponent(l, k) * U->getComponent(k, j);
	if(fabs(dblval) < fabs(this->getComponent(l, j) - sum))
	{
	  dblval = fabs(this->getComponent(l, j) - sum);
	  pivod = l;
	}
      }
    
#endif
#if 0
      int pivod = -1;
      int Rank = this->getrank();
      double *maxvals = new double[Rank];
      double *gamma = new double[Rank];

      for(int k = 0; k < Rank; ++k)
      {
	maxvals[k] = 0.0;
	for(int l = 0; l < j; ++l)
	{
	  maxvals[k] = max(maxvals[k], L->getComponent(k, l));
	} 

	for(int l = 0; l < j; ++l)
	{
	  //if(maxvals[k] == 0.0) cout << "warning in select pivod" << endl;
	  //L->a[k][l] = L->a[k][l] / maxvals[k];
	  //this->a[k][l] = this->a[k][l] / maxvals[k];
	}
      }

      double maxval = -1.0e+32;
      for(int k = j; k < Rank; ++k)
      {
	gamma[k] = this->a[k][j];
	for(int l = 0; l < j; ++l)
	{
	  gamma[k] -= L->a[k][l] * U->a[l][j];
	}

	if(maxval < gamma[k]) 
	{ 
	  maxval = gamma[k];
	  pivod = k;
	}
      }

      for(int k = 0; k < Rank; ++k)
      {
	for(int l = 0; l < j; ++l)
	{
	  //L->a[k][l] = L->a[k][l] * maxvals[k];
	  //this->a[k][l] = this->a[k][l] * maxvals[k];
	}
      }

      delete [] maxvals;
      delete [] gamma;
#endif
      return pivod;
    }

    void swapRow(int irow, int jrow)
    {
      //void *tmp = this->a[jrow];  
      //this->a[jrow] = this->a[irow];
      //this->a[irow] = static_cast<double*>(tmp);

      for(int i = 0; i < this->getrowrank(); ++i)
      {
	double tmp = this->getComponent(i, irow);
	this->setComponent(i, irow, this->getComponent(i, jrow));
	this->setComponent(i, jrow, tmp); 
      }
      return;
    }

    void swapCol(int icol, int jcol)
    {
      for(int i = 0; i < this->getrank(); ++i)
      { 
	double tmp = this->getComponent(i, icol);
	this->setComponent(i, icol, this->getComponent(i, jcol));
	this->setComponent(i, jcol, tmp);
      }
    }
    AbstMatrix  operator *(AbstMatrix &right)
    {
      int Rank = this->getrank();
      AbstMatrix matrix(Rank);

      for(int i = 0; i < Rank; ++i)
      {
	for(int j = 0; j < Rank; ++j)
	{
	  for(int k = 0; k < Rank; ++k)
	  {
	    matrix.setComponent(i, j, matrix.getComponent(i, j) + getComponent(i, k) * right.getComponent(k, j));
	  }
	}
      }
      return matrix;
    } 

    AbstVector operator *(AbstVector *rvec);

    void Cross(AbstMatrix *A,AbstMatrix *B) 
    {

      int Rank = this->getrank();

      for(int i = 0; i < Rank; ++i)
	for(int j = 0; j < Rank; ++j)
	  for(int k = 0; k < Rank; ++k)
	    //this->a[i][j] += A->getComponent(i, k) * B->getComponent(k, j);
	    this->setComponent(i, j, this->getComponent(i, j) + A->getComponent(i, k) * B->getComponent(k, j));
      return;
    } 

    void setPreNonZeroColIndex(int irow, int icol, int value)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      int index = icol + ncol * irow;
      this->preNonZeroColIndex[index] = value;
    }

    void setPreNonZeroRowIndex(int irow, int icol, int value)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      int index = icol + ncol * irow;
      this->preNonZeroRowIndex[index] = value;
    }

    void setNextNonZeroRowIndex(int irow, int icol, int value)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();

      int index = icol + ncol * irow;
      this->nextNonZeroRowIndex[index] = value;
    }

    void setNextNonZeroColIndex(int irow, int icol, int value)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();

      int index = icol + ncol * irow;
      this->nextNonZeroColIndex[index] = value;
    }

    int getPreNonZeroColIndex(int irow, int icol)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;
      int index = icol + ncol* irow;
      return this->preNonZeroColIndex[index];
    }

    int getPreNonZeroRowIndex(int irow, int icol)
    {
      //const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;
      int index = icol + ncol * irow;
      return this->preNonZeroRowIndex[index];
    }

    int getNextNonZeroColIndex(int irow, int icol)
    {
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;
      int index = icol + ncol * irow;
      return this->nextNonZeroColIndex[index];
    }

    int getNextNonZeroRowIndex(int irow, int icol)
    {
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;
      int index = icol + ncol * irow;
      return this->nextNonZeroRowIndex[index];
    }

    int getFirstColIndex(int irow, int icol)
    {
      const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;

      return getNextNonZeroColIndex(irow, icol);
      
      int index = icol + ncol * irow;
      double dbl = this->getComponent(irow, icol);
      while( dbl == 0.0)
      {
	++icol;
	index = icol + ncol * irow;
	if(index >= nrow * ncol) return -1;
	dbl = this->getComponent(irow, icol);
      }

      return icol;
    }

    int getFirstRowIndex(int irow, int icol)
    {
      const int nrow = this->getrowrank();
      const int ncol = this->getcolrank();
      if(!(irow < nrow && icol < ncol)) return -1;

      return getNextNonZeroRowIndex(irow, icol);

      int index = icol + ncol * irow;
      double dbl = this->getComponent(irow, icol);
      while( dbl == 0.0)
      {
        ++irow;
	index = icol + ncol * irow;
	if(index >= nrow * ncol) return -1;
	dbl = this->getComponent(irow, icol);
      }	
      
      return irow;
    }

     int getLastColIndex(int irow, int icol)
     {
       //const int nrow = this->getrowrank();
       const int ncol = this->getcolrank();

       if(!(irow < nrow && icol < ncol)) return -1;

       return this->getPreNonZeroRowIndex(irow, icol);

       int index = icol + ncol * irow;
       double dbl = this->getComponent(irow, icol);
       while( dbl == 0.0 )
       {
         --icol;
	 index = icol + ncol * irow;
	 if(index < 0 ) return -1;
	 dbl = this->getComponent(irow, icol);
       }

       return icol;
     }

     int getLastRowIndex(int irow, int icol)
     {
       const int nrow = this->getrowrank();
       const int ncol = this->getcolrank();

       if(!(irow < nrow && icol < ncol)) return -1;

       return this->getPreNonZeroColIndex(irow, icol);
       
       int index = icol + ncol * irow;
       double dbl = this->getComponent(irow, icol);
       while( dbl == 0.0 )
       {
         --irow;
	 index = icol + ncol * irow;
	 if( index < 0 ) return -1;
	 dbl = this->getComponent(irow, icol);
       }

       return irow;
     }

     void setNonZero()
     {
       //const int nrow = this->getrowrank();
       //const int ncol = this->getcolrank();

       for(int irow = 0; irow < nrow; ++irow)
       {
         int oldcol = -1;

	 for(int icol = 0; icol < ncol; ++icol)
	 {
	   if( this->getComponent(irow, icol) != 0.0)
	   {
	   //  cout << irow << " " << icol << "/ " << nrow << " " << ncol << "  " <<     icol + ncol * irow << endl; 
	     setPreNonZeroColIndex(irow, icol, oldcol);
	     if( 0 <= oldcol && oldcol < ncol)
	     {
	       setNextNonZeroColIndex(irow, oldcol, icol);
	     }
	     //setNextNonZeroColIndex(irow, icol, -1);
	     oldcol = icol;
	   } else 
	   {
	     //if(0 <= oldcol && oldcol < ncol)
	     {
	       setPreNonZeroColIndex(irow, icol, oldcol);
	     }
	   }
	 }
       }

       for(int icol = 0; icol < ncol; ++icol)
       {
         int oldrow = -1;
	 for(int irow = 0; irow < nrow; ++irow)
	 {
	   if( this->getComponent(irow, icol) != 0.0 )
	   {
	     setPreNonZeroRowIndex(irow,  icol,  oldrow);
	     if(0 <= oldrow && oldrow < nrow)
	     {
	       setNextNonZeroRowIndex(oldrow,  icol,  irow);
	       
	     }
	       setNextNonZeroRowIndex(irow, icol, -1);
	     oldrow = irow;
	   } else 
	   {
	     //if(0 <= oldrow && oldrow < nrow)
	     {
	       setPreNonZeroRowIndex(irow, icol, oldrow);
	     }
	   } 
	 }
       }

       for(int irow = 0; irow < nrow; ++irow)
       {
         int index = -1; 
	 for(int icol = ncol - 1; icol >= 0; --icol)
	 {
	    if( this->getComponent(irow, icol) == 0.0)
	    {
	      //if(0 <= index && index < ncol)
	      {
		setNextNonZeroColIndex(irow, icol, index);
	      }
	    } else 
	    {
	      index = icol;
	    } 
	 }
       }

       for(int icol = 0; icol < ncol; ++icol)
       {
         int index = -1;
	 for(int irow = nrow - 1; irow >= 0; --irow)
	 {
	   if( this->getComponent(irow, icol) == 0.0)
	   {
	     //if(0 <= index && index < nrow)
	     {
	       setNextNonZeroRowIndex(irow, icol, index);
	     } 
	   } else 
	   {
	     index = irow;
	   } 
	 }
       }

       return;
     }


     void setNonZeroNonZeroFirstRow()
     {
       int irow = -1;
       //const int ncol = this->getcolrank();

       //ffor(map<int, double>::iterator ite = compmap.begin(); ite != compmap.end(); ++ite)

       for(map<int, ComponentClass>::iterator ite = compmap2.begin(); ite != compmap2.end(); ++ite)
       {
	 //if(irow != (int) ite->first / ncol) 
	 if(irow != ite->second.i)
	 {
	   //irow = ite->first / ncol;
	   //firstNonZeroRowIndex[irow] = ite->first % ncol;
	   irow = ite->second.i;
	   firstNonZeroRowIndex[irow] = ite->second.j;
	 }  
       }

       return;
     }


};

class AbstVector
{
  private:
    int rank;
  public:

    double *c;
    map<int, double> cmap;
    AbstMatrix *L, *U;

    int *preNonZeroIndex;
    int *nextNonZeroIndex;
    int firstNonZeroIndex;
    int lastNonZeroIndex;

    vector<int> oldIndex;
    vector<AbstVector> oldx;

    AbstVector()
    {
    }
    AbstVector(int rank)
    {
      this->rank = rank;
      this->c = new double[rank];
      this->preNonZeroIndex = new int[rank];
      this->nextNonZeroIndex = new int[rank]; 
      this->firstNonZeroIndex = -1;
      this->lastNonZeroIndex = -1;
      
      for(int i = 0; i < rank; ++i)
      {
	this->c[i] = 0.0;
	this->preNonZeroIndex[i] = -1;
	this->nextNonZeroIndex[i] = -1;
      }

      L = NULL;
      U = NULL;
    }

    AbstVector(int rank, double *tmpc)
    {
      this->rank = rank;
      this->c = new double[rank];
      this->preNonZeroIndex = new int[rank];
      this->nextNonZeroIndex = new int[rank];
      this->firstNonZeroIndex = -1;
      this->lastNonZeroIndex = -1;
      for(int i = 0; i < rank; ++i)
      {
	this->c[i] = tmpc[i];
	this->preNonZeroIndex[i] = -1;
	this->nextNonZeroIndex[i] = -1;
      }

      L = NULL;
      U = NULL;
    }

    void clear()
    {
      delete [] this->c;

      if(L != NULL)
      {
	L->clear();
	delete L;
      }
      if(U != NULL)
      {	
	U->clear();
	delete U;
      }

      delete [] preNonZeroIndex;
      delete [] nextNonZeroIndex;
    }

    inline void setrank(int val){this->rank = val;}
    inline int getrank(){return this->rank;}
    //inline void setComponent(const int i, double val){this->c[i] = val;}
    //inline double getComponent(const int i){return this->c[i];}

    inline void setComponent(const int i, double val)
    {
      this->c[i] = val;
      //if(val != 0.0) 
      if(fabs(val) > 1.0e-16)
	this->cmap[i] = val;
      //else if( val == 0.0)
      else if(fabs(val) < 1.0e-16)
      {
        map<int, double>::iterator ite = this->cmap.find(i);
        if(ite != this->cmap.end())
          this->cmap.erase(ite);	  
      }
      return;
    }
    inline double getComponent(const int i)
    {
      //return this->cmap.find(i) != this->cmap.end() ? this->cmap[i] : 0.0;
      return this->c[i];
    }

    inline void setU(AbstMatrix *value)
    {
      const int nrank = this->getrank();
      this->U = new AbstMatrix(nrank);
      for(int i = 0; i < nrank; ++i)
      {
        for(int j = 0; j < nrank; ++j)
	{
   	  this->U->setComponent(i, j, value->getComponent(i, j));	  
	}
      }	
      //this->U = value;
    }
    inline AbstMatrix* getU(){return this->U;}
    inline void setL(AbstMatrix *value)
    {
      const int nrank = this->getrank();
      this->L = new AbstMatrix(nrank);
      for(int i = 0; i < nrank; ++i)
      {
	for(int j = 0; j < nrank; ++j)
	{
	  this->L->setComponent(i, j, value->getComponent(i, j));	  
	}
      }	
      //this->L = value;
    }
    inline AbstMatrix* getL(){return this->L;}


    double operator * (const AbstVector &right) 
    {
      double dot = 0;
      if(rank != right.rank) dot = 1.0e+10;

      for(int i = 0; i < rank; ++i)
	dot += c[i] * right.c[i];

      return dot;
    }

    inline double Scalar()
    {
      return *this * *this;
    }    

    void Add(AbstVector *x, AbstVector *y)
    {
      for(int i = 0; i < this->getrank(); ++i) this->c[i] = x->getComponent(i) + y->getComponent(i);
    }

    void Diff(AbstVector *x, AbstVector *y)
    {
      for(int i = 0; i < this->getrank(); ++i) this->c[i] = x->getComponent(i) - y->getComponent(i);
    }

    inline void Cross(AbstMatrix *A, AbstVector *vec)
    {
      const int Rank = this->getrank();
      /*
      for(int i = 0; i < Rank; ++i)
      {
	this->c[i] = 0.0;
	for(int j = 0; j < vec->getrank(); ++j)
	{
	  this->c[i] += A->getComponent(i, j) * vec->getComponent(j);	  
	}
      }
      */
/*
      for(int i = 0; i < Rank; ++i)
      {
        this->c[i] = 0.0;
	for(int j = vec->firstNonZeroIndex; j != -1; j = vec->getNextIndex(j))
	{
	  this->c[i] += A->getComponent(i, j) * vec->getComponent(j);
	}
      }
*/      
      int irow = -1, icol = -1;
      double dblvec = 0.0;

      this->cmap.clear();
      double *tmpc = this->c;
      for(int i = 0; i < Rank; ++i) this->c[i] = 0.0;

      //for(map<int, double>::iterator ite = A->compmap.begin(); ite != A->compmap.end(); ++ite)
      for(map<int, ComponentClass>::iterator ite= A->compmap2.begin(); ite != A->compmap2.end(); ++ite)
      {
	//int a = ite->first / ncol;
	//int a = ite->second.i;

	  //dblvec = this->getComponent(irow); 
	 // dblvec = this->c[irow];
	//icol = ite->first % ncol;
	icol = ite->second.j;
	double dd = vec->getComponent(icol);
	if(dd == 0.0) continue;
	
	//0112irow = ite->second.i;
	
	//dblvec += ite->second.val * vec->getComponent(icol);
	//this->setComponent(irow, dblvec);
	
	//0113dblvec = ite->second.val * dd;
	
	//this->c[irow] += dblvec;
	//if(dblvec == 0.0) continue;
	//this->setComponent(irow, this->getComponent(irow) + dblvec);
	
//	this->c[irow] += dblvec;
	tmpc[ite->second.i] += ite->second.val * dd;
      }

      for(int i = 0; i < this->getrank(); ++i)
      {
	if(this->c[i] != 0.0)
	  cmap[i] = this->c[i];
        //this->setComponent(i, this->c[i]);
      }

      return ;
    }

    inline void setPreIndex(int index, int value)
    {
      this->preNonZeroIndex[index] = value;
      return;
    }

    inline void setNextIndex(int index, int value)
    {
      this->nextNonZeroIndex[index] = value;
      return;
    }

    inline int getPreIndex(int index)
    {
      return this->preNonZeroIndex[index];
    }
    inline int getNextIndex(int index)
    {
      return this->nextNonZeroIndex[index]; 
    }

    int getFirstIndex(int index)
    {
      double ci = this->c[index];
      int rval = index;

      while(ci == 0.0)
      {
        ++rval; 
	if(rval >= rank) break;
	ci = this->c[rval];
      }
      if(ci == 0.0) rval = -1; 

      return rval;
    }

    void setNonZero()
    {
      const int rank = this->getrank();

      //set non zero col
      int old = -1;
      this->firstNonZeroIndex = rank;
      this->lastNonZeroIndex = -1;
      for(int i = 0; i < rank; ++i)
      {
	double ci = this->c[i];
	//if(ci != 0.0)
	if(fabs(ci) > 1.0e-16)
	{
	  if(0 <= old && old < rank)
	  {
	    setPreIndex(i, old);
	    setNextIndex(old, i);
	  }
	  old = i;
	  if(this->firstNonZeroIndex > i) 
	    this->firstNonZeroIndex = i;
	  if(this->lastNonZeroIndex < i ) 
            this->lastNonZeroIndex = i;
	}
      }

      return;
    }


    void updatevec(AbstVector drvec, const int startIdx, const int endIdx, double alpha)
    {
      for(int Idx = startIdx, cIdx = 0; Idx < endIdx; ++Idx, ++cIdx)
	//this->c[cIdx] += drvec.c[Idx] * alpha;
	this->setComponent(cIdx, this->getComponent(cIdx) + drvec.getComponent(Idx) * alpha);
    }

    void SimultaneousLinearEquations(AbstMatrix *A, AbstVector &bvec, enmSimultaneousLinearEquations enmMethod)
    {
      switch(enmMethod)
      {
	case enmGaussJordan:
	  GaussJordan(A, bvec);
	  break;
	case enmLU:
	  LUdecomposition(A, bvec);
	  break;
	case enmIterationMethod:
	  IterationMethod(A, bvec);
	  break;
      }
    }

    void LUdecomposition(AbstMatrix *A, AbstVector &bvec)
    {
      int Rank = this->getrank();
      const double eps = 1.0e-10;
      vector<int> swapLeft, swapRight;
      vector<int> swapLeft2, swapRight2;
      swapLeft.reserve(10); swapRight.reserve(10);
      swapLeft2.reserve(10); swapRight2.reserve(10);
      
      AbstMatrix *M = new AbstMatrix(Rank);
      
      for(int i = 0; i < Rank; ++i)
      {
        for(int j = 0; j < Rank; ++j)
	{
	  M->setComponent(i, j, A->getComponent(i, j));
	}
      }

      if(this->L == NULL) this->L = new AbstMatrix(Rank);
      if(this->U == NULL) this->U = new AbstMatrix(Rank);

      M->setNonZero();

      for(int i = 0; i < Rank; ++i)
      {
	int maxnum = 0; 
	int index = -1;
	for(int j = i; j < Rank; ++j)
	{
	  int num = 0;
	  /*
	  for(int k = i; k < Rank; ++k)
	  {
	    if(M->getComponent(j, k) == 0.0)
	    {
	      ++num; 
	    }	
	  } 
	  */
	  for(int k = M->getFirstColIndex(j, i); k != -1; k = M->getNextNonZeroColIndex(j, k))
	  {
	    ++num; 
	  }
	  
	  num = Rank - num;

	  if(num > maxnum)
	  {
	    maxnum = num;
	    index = j;
	  }
	}
	//swap row
	if(index != -1)
	{
	  M->swapRow(index, i);
	  swapLeft.push_back(index);
	  swapRight.push_back(i);
	  bvec.swap(i, index);
	}

	maxnum = 0;
	index = -1;
	for(int j = i; j < Rank; ++j)
	{
	  int num = 0;
	  /*for(int k = i; k < Rank; ++k)
	  {
	    if(M->getComponent(k, j) == 0.0)
	      ++num; 
	  } 
	  */
	  for(int k = M->getFirstRowIndex(i,j); k != -1; k = M->getNextNonZeroRowIndex(k, j))
	  {
	    ++num;
	  }
	  num = Rank - num;
	  if(num > maxnum) 
	  {
	    maxnum = num;
	    index = j;
	  }
	}

	if(index != -1)
	{
	  M->swapCol(index, i);
	  swapLeft2.push_back(index);
	  swapRight2.push_back(i);
	}

	double aii = M->getComponent(i, i);

	if(aii == 0.0) 
	{
	  int pivodrow = M->pivodOperator(i, i, eps);
	  if(pivodrow >= 0)
	  {
	    M->swapRow(i, pivodrow);
	    bvec.swap(i, pivodrow);
	    swapLeft.push_back(pivodrow);
	    swapRight.push_back(i);
	    aii = M->getComponent(i, i);
	  }
	}

	if(aii == 0.0) {cout << "wargnin ==================";}

	for(int j = i + 1; j < Rank; ++j)
	//for(int j = M->getFirstRowIndex(i + 1, i); j != -1; j = M->getNextNonZeroRowIndex(j, i))
	{
	  double c = -M->getComponent(j, i) / aii; 	   

	  if(c == 0.0) continue;
	  
	  for(int k = i+1; k < Rank; ++k)
	  //for(int k = M->getFirstColIndex(i, i + 1); k != -1; k = M->getNextNonZeroColIndex(i, k))  
	  {
	    M->setComponent(j, k, M->getComponent(j, k) + c * M->getComponent(i, k));
	  }
	}
      }
      
      for(int i = 0; i < Rank; ++i)
      {
	this->L->setComponent(i, i, 1.0);
	this->U->setComponent(i, i, M->getComponent(i, i));
	
	//for(int j = 0; j < Rank; ++j)
	for(int j = M->getFirstColIndex(i, 0); j != -1; j = M->getNextNonZeroColIndex(i ,j))
	{
	  if(j < i)
	  {
	    this->L->setComponent(i, j, M->getComponent(i, j) / M->getComponent(j,j));
	  } else
	  {
	    this->U->setComponent(i, j, M->getComponent(i, j));
	  }
	}	  
      }
      
      //L->setNonZero();
      //U->setNonZero();

      this->Lforwardsubs(bvec);
      this->Ubackwardsubs();

      for(int i = swapLeft.size() - 1 ; i >= 0; --i)      
      {
	bvec.swap(swapLeft[i], swapRight[i]);
      }
      for(int i = swapLeft2.size() - 1; i >= 0; --i)
      {
	this->swap(swapLeft2[i], swapRight2[i]);
      }

      M->clear();
      delete M;

      return;
    }

    void GaussJordan(AbstMatrix *A, AbstVector &bvec)
    {
      double buf = 0.0;
      const double eps = 1.0e-10;
      int Rank = this->getrank();

      //AbstMatrix Mhat(Rank, A->a);
      AbstMatrix Mhat(Rank);

      for(int i = 0; i < Rank; ++i)
      {
        for(int j = 0; j < Rank; ++j)
	{
	  Mhat.setComponent(i, j, A->getComponent(i, j));
	}
      }

      AbstMatrix *Ahat = &Mhat;

      this->cmap.clear();
      for(int i = 0; i < Rank; ++i)
      {
	//this->c[i] = bvec.getComponent(i);
	this->setComponent(i, bvec.getComponent(i));
      }

      for(int i = 0; i < Rank; ++i)
      {

	//if(fabs(Ahat->a[i][i]) < eps)
	if(fabs(Ahat->getComponent(i, i)) < eps)
	{
	  int pivodrow = Ahat->pivodOperator(i, i, eps);
	  if(pivodrow >= 0)
	  {
	    Ahat->swapRow(i, pivodrow);
	    this->swap(i, pivodrow);
	  }
	  else
	  {
	    return;
	  }
	}

	buf = 1.0 / Ahat->getComponent(i, i);

	for(int j = 0; j < Rank;++j)
	{
	  Ahat->setComponent(i, j, Ahat->getComponent(i,j) * buf);
	}

	this->c[i] *= buf;

	for(int j = 0; j < Rank; ++j)
	{
	  if( i != j )
	  {
	    buf = Ahat->getComponent(j, i);
	    for(int k = 0; k < Rank; ++k)
	    {
	      Ahat->setComponent(j, k, Ahat->getComponent(j, k) - Ahat->getComponent(i, k) * buf);
	    }

	    this->c[j]  -= this->c[i] * buf;
	  }
	}
      }

      Ahat->clear();
    }

    void Lforwardsubs(AbstVector &bvec)
    {
      const int Rank = this->getrank();

      this->cmap.clear();
      for(int i = 0; i < Rank; ++i) 
	this->setComponent(i, bvec.getComponent(i));

      if(this->L->getMatrixType() == unit) return;

      this->setNonZero();

      //for(int i = 0; i < Rank; ++i)
      for(int i = this->firstNonZeroIndex; i != -1; i = this->getNextIndex(i))
      {
	this->setComponent(i, this->getComponent(i) / this->L->getComponent(i, i));

	//for(int j = L->getFirstRowIndex(i+1, i); j != -1; j = L->getNextNonZeroRowIndex(j, i))
	for(int j = i + 1; j < Rank; ++j)
	{
	  this->setComponent(j , this->getComponent(j) - this->getComponent(i) * this->L->getComponent(j , i));
	}
	
      }
    }

    void Ubackwardsubs()
    {
      //const int Rank = this->getrank();
      if(this->U->getMatrixType() == unit) return;
      
      this->setNonZero();

      //for(int i = Rank - 1; i >= 0; --i)
      for(int i = this->lastNonZeroIndex; i != -1; i = this->getPreIndex(i))
      {
	this->setComponent(i, this->getComponent(i) / this->U->getComponent(i, i));

	//for(int j = U->getLastRowIndex(i - 1, i); j != -1; j = U->getPreNonZeroRowIndex(j, i))
	for(int j = i - 1; j >= 0; --j)
	{
	  this->setComponent(j, this->getComponent(j) - this->getComponent(i) * this->U->getComponent(j , i));
	}
      }
    }

    void IterationMethod(AbstMatrix *A, AbstVector &bvec)
    {
      int Rank = this->getrank();
      AbstVector dx(Rank);
      AbstVector bold(Rank);
      const double eps = 1.0e-10;
      double e = 1.0e+16;

      for(int i = 0; i < Rank; ++i) bold.c[i] = bvec.getComponent(i);

      this->GaussJordan(A, bvec);

      while(e > eps)
      {
	bvec.Cross(A, this);
	bvec.Diff(&bvec, &bold);

	e = bvec.Scalar();

	dx.GaussJordan(A, bvec);
	this->Diff(this, &dx);
      }
      dx.clear();
      bold.clear();
    }

    void swap(int i, int j)
    {
      double tmp = this->getComponent(i);
      this->setComponent(i, this->getComponent(j));
      this->setComponent(j, tmp);
      //this->c[i] = this->getComponent(j);
      //this->c[j] = tmp;
    }

    void LUdecompositionReusing(AbstMatrix *A, AbstVector &bvec)
    {
      const int Rank = this->getrank();
      const int count = this->oldx.size();


      this->cmap.clear();
      this->Lforwardsubs(bvec);
      this->Ubackwardsubs();

      for(int i = 0; i < count; ++i)
      {
	int icol = oldIndex[i];
	double dbl = this->getComponent(icol);
	if( dbl == 0.0) continue;

	for(map<int, double>::iterator ite = oldx[i].cmap.begin(); ite != oldx[i].cmap.end(); ++ite)
	{
  	   int j = ite->first; 
	   if(j == icol) continue;
	   //this->c[j] = this->c[j] - dbl / oldx[i].c[icol] * (oldx[i].c[j] - ( j == icol ? 1 : 0));
	   //this->c[j] = this->c[j] - dbl / oldx[i].c[icol] * oldx[i].c[j];

	   this->setComponent(j, this->getComponent(j) - dbl / oldx[i].getComponent(icol) * oldx[i].getComponent(j)); 
	}

	//this->c[icol] = this->c[icol] - dbl / oldx[i].c[icol] * (oldx[i].c[icol] - 1.0);
	this->setComponent(icol, this->getComponent(icol) - dbl / oldx[i].getComponent(icol) * (oldx[i].getComponent(icol) - 1.0));
	  
      }

      for(int j = 0; j < Rank; ++j)
      {
	if(this->c[j] != 0.0) 
	  this->cmap[j] = this->c[j];
      }

      return;
    }

    void LUdecompositionReusingForTrace(AbstMatrix *A, AbstVector &bvec)
    {
      const int Rank = this->getrank();
      const int count = this->oldx.size();

      bvec.cmap.clear();

      for(int i = 0; i < Rank; ++i)
      {
        if(bvec.c[i] != 0)
	{
	  bvec.cmap[i] = bvec.c[i];
	}
      }

      for(int i = count - 1; i >= 0; --i)
      {
	int icol = oldIndex[i];
	double dbl = 0.0;
	
	for(map<int, double>::iterator ite = bvec.cmap.begin(); ite != bvec.cmap.end(); ++ite)
	{
          int j = ite->first;
	  dbl += (j == icol ? oldx[i].c[j] - 1.0 : oldx[i].c[j] ) * bvec.c[j] / oldx[i].c[icol];
	}

	bvec.setComponent(icol, bvec.getComponent(icol) - dbl);
      }

      this->Lforwardsubs(bvec);
      this->Ubackwardsubs();

    }

    void calcLUMatrixForUnitMatrix()
    {
      const int nrank = this->getrank();
      if(this->L == NULL) this->L = new AbstMatrix(nrank); 
      if(this->U == NULL) this->U = new AbstMatrix(nrank);

      this->L->setMatrixType(unit);
      this->U->setMatrixType(unit);
      
      for(int i = 0; i < nrank; ++i)
      {
        this->L->setComponent(i, i, 1.0); 
	this->U->setComponent(i, i, 1.0);
      }
    
      return;
    }

};


#endif
