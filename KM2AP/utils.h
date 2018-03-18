#ifndef _UTILS_H
#define _UTILS_H
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

class AbstMatrix;
class AbstVector;

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

class AbstMatrix
{
  private:
    int rank;
  public:
    double **a;
    double det;
    double *a2;
    AbstMatrix *Inv;

    AbstMatrix(){}
    AbstMatrix(int rank)
    {
      this->rank = rank;

      a = new double*[rank];
      for(int i = 0; i < rank; ++i)
      {
	a[i] = new double[rank]; 
	for(int j = 0; j < rank; ++j)
	{
	  a[i][j] = 0.0;
	}
      }

      a2 = new double[rank*rank];
      for(int i = 0; i < rank*rank; ++i)
      {
	a2[i] = 0.0;
      }
    }
    AbstMatrix(int rank, double **aij)
    {
      this->rank = rank;

      a = new double *[rank];
      a2 = new double[rank*rank];

      for(int i = 0;i < rank; ++i) a[i] = new double[rank];

      int ii = 0;
      for(int i = 0; i < rank; ++i)
      {
	for(int j = 0; j < rank; ++j)
	{
	  this->a[i][j] = aij[i][j];
	  this->a2[ii] = aij[i][j];
	  ++ii;
	}
      }
    }

    ~AbstMatrix()
    {
    }

    void clear()
    {  
      const int Rank = this->getrank();
      for(int i = 0; i < Rank; ++i)
      {	
	delete [] this->a[i];
      }
      delete [] this->a;
      delete [] a2;
    }


    inline void setrank(int rank){this->rank = rank;}
    inline int getrank(){return this->rank;}
    inline void setComponent(int i, int j, double value){this->a[i][j] = value;}
    inline double getComponent(int i, int j){return this->a[i][j];}

    inline double getdet(){return this->det;}
    void calcdet()
    {
      const double eps = 0.0;
      double buf = 0.0;
      int Rank = this->getrank();
      AbstMatrix Ahat(Rank, this->a);

      for(int i = 0; i < Rank; ++i)
	for(int j = 0; j < Rank; ++j)
	  if(i < j) 
	  {
	    if(Ahat.a[i][i] == 0.0) 
	    {
	      int pivodrow = Ahat.pivodOperator(i, i, eps);
	      if(pivodrow >= 0)
		Ahat.swapRow(i, pivodrow);
	      else 
		return;
	    }
	    buf = Ahat.a[j][i] / Ahat.a[i][i];

	    for(int k = 0;k < Rank; ++k) Ahat.a[j][k] -= Ahat.a[i][k] * buf; 
	  }

      for(int i = 0; i < Rank; ++i) det *= Ahat.a[i][i];

      Ahat.clear();
    } 

    void calcInverseMatrix()
    {
      double buf = 0.0;
      const double eps = 1.0e-10;
      int Rank = this->getrank();

      AbstMatrix Mhat(Rank, this->a);
      AbstMatrix *Ahat = &Mhat;
      AbstMatrix *Bhat = Inv;

      for(int i = 0; i < Rank; ++i)
      {
	for(int j = 0; j < Rank; ++j)
	{
	  Bhat->a[i][j] = i == j ? 1.0 : 0.0;
	}
      }

      for(int i = 0; i < Rank; ++i)
      {

	if(fabs(Ahat->a[i][i]) < eps)
	  //if(i + 1 != Rank)
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

	buf = 1 / Ahat->a[i][i];

	for(int j = 0; j < Rank;++j)
	{
	  Ahat->a[i][j] *= buf;
	  Bhat->a[i][j] *= buf;
	}

	for(int j = 0; j < Rank; ++j)
	{
	  if( i != j )
	  {
	    buf = Ahat->a[j][i];
	    for(int k = 0; k < Rank; ++k)
	    {
	      Ahat->a[j][k] -= Ahat->a[i][k] * buf;
	      Bhat->a[j][k] -= Bhat->a[i][k] * buf;
	    }
	  }
	}
      }

      Ahat->clear();

    }

    AbstMatrix* getInverseMatrix()
    {
      return this->Inv;
    }

    int pivodOperator(int i, int j, const double eps)
    {
      int Rank = this->getrank();
      int pivod = -1;
      double val = 0.0;
      double eps0 = eps;

      for(int k = i + 1; k < Rank; ++k)
      {
	if(fabs(this->a[k][j]) > eps * 0)
	{
	  if(val < fabs(this->a[k][j]))
	  {
	    pivod = k;
	    val = fabs(this->a[k][j]);
	  }
	  //return k;
	}
      }	

      if(pivod >= 0) return pivod;

      while(eps0 > 1.0e-16)
      {
	eps0 /= 2;
	for(int k = i + 1; k < Rank; ++k)
	{
	  if(fabs(this->a[k][j]) >  eps0)
	  {
	    if(val < fabs(this->a[k][j]))
	    {
	      pivod = k;
	      val = fabs(this->a[k][j]);
	    }
	    //return k;
	  }
	}
      }

      return pivod;
    }

    int pivodOperatorLU(int i, int j, const double eps, AbstMatrix *L, AbstMatrix *U)
    {
      int pivod = -1;
      int Rank = this->getrank();
      double dblval = 0.0;
      double sum = 0.0;

      for(int l = j; l < Rank; ++l)
      {
	sum = 0.0;
	for(int k = 0; k < l; ++k)
	  sum += L->getComponent(l, k) * U->getComponent(k, j);
	if(fabs(dblval) < fabs(this->getComponent(l ,j) - sum))
	{
	  dblval = fabs(this->getComponent(l, j) - sum);
	  pivod = l;
	}
      }

      return pivod;
    }

    void swapRow(int irow, int jrow)
    {
      void *tmp = this->a[jrow];  

      this->a[jrow] = this->a[irow];
      this->a[irow] = static_cast<double*>(tmp);
#if 1 
      for(int icol = 0; icol < rank; ++icol)
      {
	double dbltmp = this->a2[icol + rank * irow];
	this->a2[icol + rank * irow] = this->a2[icol + rank *jrow];
	this->a2[icol + rank * jrow] = dbltmp;
      }	
#endif 
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
	    matrix.a[i][j] += getComponent(i, k) * right.getComponent(k, j);
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
	    this->a[i][j] += A->getComponent(i, k) * B->getComponent(k, j);
    } 

};

class AbstVector
{
  private:
    int rank;
  public:

    double *c;

    vector<int> oldIndex;
    vector<AbstVector> oldx;
    AbstMatrix *L, *U;

    AbstVector()
    {
    }
    AbstVector(int rank)
    {
      this->rank = rank;
      this->c = new double[rank];
      for(int i = 0; i < rank; ++i)
      {
	this->c[i] = 0.0; 
      }

      L = NULL;
      U = NULL;

    }
    AbstVector(int rank, double *c)
    {
      this->rank = rank;
      this->c = new double[rank];

      for(int i = 0; i < rank; ++i)
      {
	this->c[i] = c[i];
      }
      L =NULL;
      U = NULL;
    }

    void clear()
    {
      if(U != NULL) delete U;
      if(L != NULL) delete L;

      delete [] this->c;
    }

    inline void setrank(int val){this->rank = val;}
    inline int getrank(){return this->rank;}
    inline double getComponent(const int i){return this->c[i];}
    inline void setComponent(int i, double val){this->c[i] = val;}

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

    void Cross(AbstMatrix *A, AbstVector *vec)
    {
      int Rank = this->getrank();

      for(int i = 0; i < Rank; ++i)
      {
	this->c[i] = 0.0;
	for(int j = 0; j < Rank; ++j)
	{
	  this->c[i] += A->getComponent(i, j) * vec->getComponent(j);	  
	}
      }
    }

    void updatevec(AbstVector drvec, const int startIdx, const int endIdx, double alpha)
    {
      for(int Idx = startIdx, cIdx = 0; Idx < endIdx; ++Idx, ++cIdx)
	this->c[cIdx] += drvec.c[Idx] * alpha;
    }

    void LUdecomposition(AbstMatrix *A, AbstVector &bvec)
    {
      int Rank = this->getrank();
      double sum = 0.0;
      const double eps = 1.0e-10;

      if(L == NULL) L = new AbstMatrix(Rank);
      if(U == NULL) U = new AbstMatrix(Rank);

      for(int i = 0; i < Rank; ++i)
	for(int j = 0; j < Rank; ++j)
	{
	  L->a[i][j] = i == j ? 1.0 : 0.0;
	  U->a[i][j] = 0.0;
	}

      for(int i = 0; i < Rank; ++i)
	for(int j = 0; j < Rank; ++j)
	{
	  if(i > j)
	  {
	    sum = 0.0;
	    for(int k = 0; k < j; ++k)
	      sum += L->getComponent(i, k) * U->getComponent(k, j);
	    L->a[i][j] = (A->getComponent(i, j) - sum) / U->getComponent(j, j);
	    //L->a[i][j] = A->getComponent(i, j) / A->getComponent(j, j);
	  }
	  else 
	  {
	    if(i == j)
	    {
	      sum = 0.0;
	      for(int k = 0; k < i; ++k)
		sum += L->getComponent(i, k) * U->getComponent(k, j);
	      if(fabs(A->getComponent(i, j) - sum) < eps)
	      {
		int pivod = A->pivodOperatorLU(i, j, eps, L, U);
		if(pivod >= 0)
		{
		  A->swapRow(i, pivod);
		  //L->swapRow(i, pivod);
		  bvec.swap(i, pivod);
		}
	      }
	    }

	    sum = 0.0;
	    for(int k = 0; k < i; ++k)
	      sum += L->getComponent(i, k) * U->getComponent(k, j);
	    U->a[i][j] = A->getComponent(i, j) - sum;
	    //U->a[i][j] = A->getComponent(i, j);
	  }

	}

      AbstVector yvec(Rank);

      yvec.Lforwardsubs(L, U, bvec, *this);
      this->Ubackwardsubs(L, U, yvec);

      yvec.clear();
    }

    void GaussJordan(AbstMatrix *A, AbstVector &bvec)
    {
      double buf = 0.0;
      const double eps = 1.0e-10;
      int Rank = this->getrank();

      AbstMatrix Mhat(Rank, A->a);
      AbstMatrix *Ahat = &Mhat;

      for(int i = 0; i < Rank; ++i)
      {
	this->c[i] = bvec.getComponent(i);
      }

      for(int i = 0; i < Rank; ++i)
      {

	if(fabs(Ahat->a[i][i]) < eps)
	{
	  int pivodrow = Ahat->pivodOperator(i, i, eps);
	  if(pivodrow >= 0)
	  {
	    Ahat->swapRow(i, pivodrow);
	    bvec.swap(i, pivodrow);
	  }
	  else
	  {
	    return;
	  }
	}

	buf = 1 / Ahat->a[i][i];

	for(int j = 0; j < Rank;++j)
	{
	  Ahat->a[i][j] *= buf;
	}

	this->c[i] *= buf;

	for(int j = 0; j < Rank; ++j)
	{
	  if( i != j )
	  {
	    buf = Ahat->a[j][i];
	    for(int k = 0; k < Rank; ++k)
	    {
	      Ahat->a[j][k] -= Ahat->a[i][k] * buf;
	    }

	    this->c[j]  -= this->c[i] * buf;
	  }
	}
      }

      Ahat->clear();
    }

    void Lforwardsubs(AbstMatrix *L, AbstMatrix *U, AbstVector &bvec, AbstVector &xvec)
    {
      int Rank = this->getrank();
      //const double eps = 1.0e-16;
#if 0
      cout << "L = " << endl;
      for(int i = 0; i < Rank; ++i)
      {
	for(int j = 0; j < Rank; ++j)
	{
	  cout  << L->getComponent(i, j) << " " ;
	}
	cout << endl;
      }
#endif

      for(int i = 0; i < Rank; ++i) this->c[i] = bvec.getComponent(i);

      for(int i = 0; i < Rank; ++i)
      {
	this->c[i] /= L->getComponent(i, i);

	for(int j = i + 1; j < Rank; ++j)
	  //for(int j = L->getFirstRowIndex(i+1, i); j != -1; j = L->getNextNonZeroRowIndex(j, i))
	  this->c[j] -= this->getComponent(i) * L->getComponent(j, i);
      }

    }

    void Ubackwardsubs(AbstMatrix *L, AbstMatrix *U, AbstVector &yvec)
    {
      int Rank = this->getrank();

      for(int i = 0; i < Rank; ++i) this->c[i] = yvec.getComponent(i);

      for(int i = Rank - 1; i >= 0; --i)
      {
	this->c[i] /= U->getComponent(i, i);

	for(int j = i - 1; j >= 0; --j)
	  ///for(int j = U->getLastRowIndex(i - 1, i); j != -1; j = U->getPreNonZeroRowIndex(j, i))
	  this->c[j] -= this->getComponent(i) * U->getComponent(j, i);
      }
    }

    void swap(int i, int j)
    {
      double tmp = this->getComponent(i);
      this->c[i] = this->getComponent(j);
      this->c[j] = tmp;
    }

    void LUdecompositionReusing(AbstMatrix *A, AbstMatrix *L, AbstMatrix *U, AbstVector &bvec)
    {
      int Rank = this->getrank();

      AbstVector yvec(Rank);
      AbstVector xvec(Rank);
      AbstVector u(Rank);

      //U->setNonZero();
      //L->setNonZero();

      yvec.Lforwardsubs(L, U, bvec, u);
      u.Ubackwardsubs(L, U, yvec);

      const int count = this->oldx.size();

      for(int i = 0; i < count; ++i)
      {
	int icol = oldIndex[i];
	for(int j = 0; j < Rank; ++j)
	{
	  xvec.c[j] = u.c[j] - u.c[icol] / oldx[i].c[icol] * ( oldx[i].c[j] - ( j == icol ? 1 : 0));
	} 
	for(int j = 0; j < Rank; ++j)
	{
	u.c[j] = xvec.c[j];
	}
      }

      for(int i = 0; i < Rank; ++i)
	this->c[i] = xvec.c[i];

      u.clear();
      yvec.clear();

    }

    void LUdecompositionReusingForTrace(AbstMatrix *A, AbstMatrix *L, AbstMatrix *U, AbstVector &bvec)
    {
      int Rank = this->getrank();

      AbstVector yvec(Rank);
      AbstVector xvec(Rank);
      AbstVector u(Rank);

      //U->setNonZero();
      //L->setNonZero();

      //yvec.Lforwardsubs(L, U, bvec, u);
      //u.Ubackwardsubs(L, U, yvec);

      const int count = this->oldx.size();

      for(int i = count - 1; i >= 0; --i)
      //for(int i = 0; i < Rank; ++i)
      {
	int icol = oldIndex[i];
	for(int j = 0; j < Rank; ++j)
	{
	  double dbl = 0.0;
	  for(int k = 0; k < Rank; ++k) dbl += (k == icol ? oldx[i].c[k] - 1. : oldx[i].c[k] ) * bvec.c[k] / oldx[i].c[icol]; 

	  u.c[j] = bvec.c[j] - (j == icol ? 1 :0) * dbl; 
	} 
	for(int j = 0; j < Rank; ++j)
	{
	  bvec.c[j] = u.c[j];
	}
      }

      yvec.Lforwardsubs(L, U, u, u);
      this->Ubackwardsubs(L, U, yvec);

      //for(int i = 0; i < Rank; ++i)
	//this->c[i] = u.c[i];

      u.clear();
      yvec.clear();

      cout << "hoge" << endl;
    }

};



template<typename T1, typename T2>
class myPair
{
  public:
    myPair<T1, T2> *firstPair;
    myPair<T1, T2> *nextPair;
    myPair<T1, T2> *nextPair2;
    bool flag;
    T1 key;
    T2 value;
};


template<typename T1, typename T2>
class myMap
{
  public:

    int nsize;
    hash<T1> Hash;
    myPair<T1, T2> *Pair;
    vector< myPair<T1, T2>* > PairVec; 

    myMap()
    {
      nsize = 1000000;
      Pair = new myPair<T1, T2>[nsize];
      PairVec.reserve(nsize);

      for(int i = 0; i < nsize; ++i)
      {
	Pair[i].flag = false;
	Pair[i].firstPair = NULL; 
	Pair[i].nextPair = NULL;
      }
    }

    void add(T1 &key, T2 &value)
    {
      size_t sizet = Hash(key);
      sizet = sizet % nsize;

      for(myPair<T1, T2> *tmp = Pair[sizet].firstPair; tmp != NULL; tmp = tmp->nextPair)
      {
	if(value == tmp->value) return;
      }

      Pair[sizet].flag = true;
      myPair<T1, T2> *tmpPair = new myPair<T1, T2>;
      tmpPair->key = key;
      tmpPair->value = value;

      tmpPair->nextPair = Pair[sizet].firstPair;
      Pair[sizet].firstPair = tmpPair;

      if(PairVec.size() > 0) PairVec[PairVec.size() - 1]->nextPair2 = tmpPair;
      tmpPair->nextPair2 = NULL;
      PairVec.push_back(tmpPair);
    }

    myPair<T1, T2>* find(T1 &key)
    {
      size_t sizet = Hash(key) % nsize;

      myPair<T1, T2> *tmp = Pair[sizet].firstPair; 
      for(; tmp != NULL; tmp = tmp->nextPair)
      {
	if(key == tmp->key) break;
      }

      return tmp;
    }

    myPair<T1, T2>* begin()
    {
      return PairVec[0];
    }

    myPair<T1, T2>* end()
    {
      return NULL;
    }

    T2 operator [](T1 &key)
    {
      size_t sizet = Hash(key);
      sizet = sizet % nsize;

      for(myPair<T1, T2> *tmp = Pair[sizet].firstPair; tmp != NULL; tmp = tmp->nextPair)
      {
	if( key == tmp->key) break;
      }
      return Pair[sizet].firstPair->value;
    }

    int size()
    {
      return PairVec.size();
    }

    void clear()
    {

      for(int i = 0; i < PairVec.size(); ++i)
      {
	delete PairVec[i];
      }

      delete [] Pair;
    }
};

class UnionFinding
{
  public:

    UnionFinding(){}
    UnionFinding(int nsize)
    {
      par.reserve(nsize); 
      for(int i = 0; i < nsize; ++i)
      {
	par.push_back(i);
      }
    }

    vector<int> par;

    bool Union(int x, int y)
    {
      x = root(x);
      y = root(y);
      //y = par[y]; 
      par[y] = x;
      return true;
    }

    bool Find(int x, int y)
    {
      return root(x) == root(y);
    }

    int root(int x)
    {
      if(par[x] == x) return x;
      par[x] = root(par[x]);
      return par[x];
    }
};

inline std::vector<string> split(const string &str)
{
  vector<string> v;
  char delim[] = {'=', '\t', ',', ' ', ':'};
  string item;

  if(str.size() == 0) return v;

  for(int i = 0; i < (int) str.size(); ++i)
  {
    char ch = str[i];

    if(ch == delim[0] || ch == delim[1] || ch == delim[2] || ch == delim[3] || ch == delim[4])
    {
      if(!item.empty())
      {
	v.push_back(item);
      }
      item.clear();
    }
    else 
    {
      item += ch;
    }
  } 

  if(!item.empty())
  {
    v.push_back(item);
  }

  return v;
}

  template<typename var>
inline void vartostring(string &str, var &val)
{
  stringstream ss;
  ss << val;
  str = ss.str();
}

  template<typename var>
inline void stringtovar(string &str, var &val)
{
  istringstream (str) >> val;
}

  template<typename var>
inline var sq(var x)
{
  return x * x;
}

  template<typename var>
inline var cube(var x)
{
  return x * x * x;
}

  template<typename var>
inline void rnd1(var &val)
{

  val = static_cast<var> (rand()) / (static_cast<var>(RAND_MAX));
}

inline void init_rnd()
{
  srand((unsigned)time(NULL));
}

  template<typename var>
inline void stringtimetovartime(string &str, var &val)
{
  var hour;
  var minite;
  vector<string> time = split(str);

  val = hour * 60 +  minite;
} 

  template<typename var>
inline void vartimetostringtime(string &str, var &val)
{

  var hour = val / static_cast<var>(60.0);
  var minite = val % static_cast<var>(60.0);
  string strhour;
  string strminite;

  vartostring(strhour, hour);
  vartostring(strminite, minite);

  str = strhour + ":" + strminite;
}

inline bool fileexist(string &filepath)
{
  //bool flag = false;  
  FILE *fp;
  const char *filename = filepath.c_str();
  if((fp = fopen(filename, "r")) == NULL)
  {
    return false;
  }
  else 
  {
    fclose(fp);
    return true;
  }
}

inline bool mkdir(string &str)
{
  const char *dirpath = str.c_str(); 

  mode_t  mode = S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR | S_IWGRP;

  if(mkdir(dirpath, mode) != 0)
  {
    perror("error:mkdir");
    return false;
  }

  return true;
}

  template<class T> 
inline void shuffle(T ary[],int size, bool flag)
{
  for(int i=0;i<size;i++)
  {
    int j = rand()%size;
    //if(flag && (i == size - 1 || j == size - 1)) continue;
    if(flag && (i == 0 || j == 0)) continue;
    T t = ary[i];
    ary[i] = ary[j];
    ary[j] = t;
  }
}

inline double sign(double A)
{
  return (A > 0) - (A < 0);
}

#endif 

