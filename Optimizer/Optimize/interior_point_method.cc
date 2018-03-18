#include <iostream>
#include "matrix.h"
#include "params.h"
#include "interior_point_method.h"


//namespace Interior_Point_Methodname
namespace OptimizeName
{

  void checkfunc(Matrix *Mhat);

  static const int nRank = Nx + Neq + Nconst + Nconst;
  static const int xRank = Nx;
  static const int uRank = Neq;
  static const int lambdaRank = Nconst;
  static const int zRank = Nconst;
  static const int gRank = Nconst;
  static const int hRank = Neq;
  double rho = 100.0;

  void Interior_Point_Method::run()
  {
    //Function fobj;
    //Function *g = new Function[Nconst];

    initialize();

    opt();

    free();
  }

  void Interior_Point_Method::initialize()
  {
    Function *pobj = getObjectiveFunction();
    vector<Function> *pg = getConstraintFunction();
    pobj->setfunctype(object, 0);
    for(int i = 0; i < pg->size(); ++i)
      pg->at(i).setfunctype(constraint, i);
    /*this->fobj.setfunctype(object, 0);
    for(int i = 0; i < Nconst; ++i)
      this->g[i].setfunctype(constraint, i);
      */

  }
  
  void Interior_Point_Method::free()
  {
    //delete [] g; 
    vector<Function> *g = getConstraintFunction();
    //this->g.clear();
    g->clear();
  }

  void Interior_Point_Method::opt()
  {
    Matrix *Mhat = new Matrix(nRank);
    Mhat->Inv = new AbstMatrix(nRank);
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pg = this->getConstraintFunction();

    Vector xvec(xRank);
    Vector uvec(uRank);
    Vector lambdavec(lambdaRank);
    Vector zvec(zRank);

    for(int i = 0; i < zvec.getrank(); ++i) zvec.c[i] = 1.0;
    for(int i = 0; i < lambdavec.getrank(); ++i) lambdavec.c[i] = 1.0; 

    //Vector *drvec = new Vector(nRank);
    //Vector *gradvec = new Vector(nRank);

    Vector drvec(nRank);
    Vector gradvec(nRank);

    int loop  = 0;
    double alpha = 0.0;

    while(1)
    {
      cout << "step[" << loop++ << "]" << endl;
      gradvec.setVecIPM(*pfobj, *pg, xvec, zvec, uvec, lambdavec);

      Mhat->setMatrixIPM(*pfobj, *pg, xvec, zvec, uvec, lambdavec);
      Mhat->calcInverseMatrix();
      if(Mhat->Inv == NULL)
      {
	cout << "NULL" << endl;
	getchar();
      } 

      drvec.Cross(Mhat->Inv, &gradvec);
      //drvec.LUdecomposition(Mhat, gradvec);
      
      alpha = getalpha(drvec, xvec, zvec, uvec, lambdavec);

      // cout  
      cout << "alpha = " << alpha << endl;
      cout << xvec.c[0] << " " << xvec.c[1] << 
	" " << xvec.c[2] << " "  << xvec.c[3] << endl;
      cout << "f = " << pfobj->f(xvec.c) << endl;
      for(int i = 0; i < gRank; ++i)
	cout << "g" << i << "= " << pg->at(i).f(xvec.c) << endl;

      // update vector
      xvec.updatevec(drvec, 0, xRank, alpha);
      zvec.updatevec(drvec, xRank, xRank + zRank, alpha);
      uvec.updatevec(drvec, xRank + zRank, xRank + zRank + uRank, alpha);
      lambdavec.updatevec(drvec, xRank + zRank + uRank, xRank + zRank + uRank + lambdaRank, alpha);

      getchar();

    }

    free_Matrix_Vector(Mhat, xvec, zvec, uvec, lambdavec, drvec, gradvec);

  }

  void Matrix::setMatrixIPM(Function &fobj, vector<Function> &g, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec)
  {
    // 1 Columns Block
    for(int i = 0; i < nRank; ++i)
    {
      for(int j = 0; j < xRank; ++j)
      {
	if(i < xRank)
	{
	  this->a[i][j] = fobj.dxidxj(xvec.c, i, j);
#if 0
	  for(int hIdx = 0; hIdx < hRank; ++hIdx)
	  {
	    this->a[i][j] += h[hIdx].dxidxj(xvec.c, i, j) * uvec.c[hIdx];
	  }
#endif
	  //cout << "f = " << this->a[i][j] << endl;
	  for(int gIdx = 0; gIdx < gRank; ++gIdx)
	  {
	    this->a[i][j] += g[gIdx].dxidxj(xvec.c, i, j) * lambdavec.c[gIdx];
	    //cout << g[gIdx].dxidxj(xvec.c, i, j) << endl;
	  }

	}
	else if(i < xRank + lambdaRank)
	{
	  this->a[i][j] = 0.0;
	}
	else if(i < xRank + lambdaRank + hRank)
	{
#if 0
	  //int hIdx = hRank + lambdaRank + xRank - i;
	  int hIdx = i - xRank - lambdaRank;
	  this->a[i][j] = h[hIdx].dxi(xvec.c, j);
#endif
	}
	else if(i < xRank + lambdaRank + hRank + gRank)
	{
	  //int gIdx = gRank + hRank + lambdaRank + xRank - i;
	  int gIdx = i - xRank - lambdaRank - hRank;
	  this->a[i][j] = g[gIdx].dxi(xvec.c, j);
	}
      }
    }

    //2 Column Block
    for(int i = 0; i < nRank; ++i)
    {
      for(int j = xRank; j < xRank + lambdaRank; ++j)
      {
	if(i < xRank)
	{
	  this->a[i][j] = 0.0; 
	}
	else if(i < xRank + lambdaRank)
	{
	  //int lIdx = lambdaRank + xRank - i;
	  //int lJdx = lambdaRank + xRank - j;
	  int lIdx = i - xRank;
	  int lJdx = j - xRank;
	  this->a[i][j] = lIdx == lJdx ? lambdavec.c[lIdx] : 0.0;
	}
	else if(i < xRank + lambdaRank + hRank)
	{
	  this->a[i][j] = 0.0; 
	}
	else if(i < xRank + lambdaRank + hRank + gRank)
	{
	  //int Idx = gRank + hRank + lambdaRank + xRank - i;
	  //int Jdx = xRank + lambdaRank - j;
	  int Idx = i - xRank - lambdaRank - hRank;
	  int Jdx = j - xRank;
	  this->a[i][j] = Idx == Jdx  ? 1.0 : 0.0;
	}
      }
    }

    //3 Column Block
    for(int i = 0; i < nRank; ++i)
    {
      for(int j = xRank + lambdaRank; j < xRank + lambdaRank + hRank; ++j)
      {
	if(i < xRank)
	{
#if 0
	  //int hIdx = hRank + lambdaRank + xRank - j;
	  int hIdx = j - xRank - lambdaRank;
	  this->a[i][j] = h[hIdx].dx(xvec.c, i); 
#endif
	}
	else
	  this->a[i][j] = 0.0;
      } 
    }

    //4 Column Block
    for(int i = 0; i < nRank; ++i)
    {
      for(int j = xRank + lambdaRank + hRank; j < xRank + lambdaRank + hRank + gRank; ++j)
      {
	if(i < xRank) 
	{
	  int gIdx = j - xRank - lambdaRank - hRank;
	  this->a[i][j] = g[gIdx].dxi(xvec.c, i);
	}
	else if(i < xRank + lambdaRank)
	{
	  int zIdx = i - xRank;
	  int zJdx = j - xRank - lambdaRank - hRank;
	  this->a[i][j] = zIdx == zJdx ? zvec.c[zIdx] : 0.0; 

	}
	else 
	{
	  this->a[i][j] = 0.0;
	}
      }
    }

    return;
  }

  void Vector::setVecIPM(Function &fobj, vector<Function> &g, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec)
  {
    cout << "z * l = " << zvec * lambdavec << endl;
    cout << "rho = " << rho << endl; 

    for(int i = 0; i < xRank; ++i)
    {
      this->c[i] = -fobj.dxi(xvec.c, i);

      for(int gIdx = 0; gIdx < gRank; ++gIdx)
      {
	this->c[i] += -g[gIdx].dxi(xvec.c, i) * lambdavec.c[gIdx]; 
      }	
#if 0
      for(int hIdx = 0; hIdx < hRank; ++hIdx)
      {
	this->c[i] = -h[hIdx].dxi(xvec.c, v) * uvec.c[hIdx];
      }
#endif
    }

    for(int i = xRank, lIdx = 0; lIdx < lambdaRank; ++i, ++lIdx)
    {
      this->c[i] = -(zvec.c[lIdx] * lambdavec.c[lIdx] - rho); 
    }
#if 0
    for(int i = xRank + lambdaRank, hIdx = 0+ hIdx < hRank; ++i, ++hIdx)
    {
      this->c[i] = -h[hIdx].f(xvec.c);
    }
#endif

    for(int i = xRank + lambdaRank + hRank, gIdx = 0; gIdx < gRank; ++i, ++gIdx)
    {
      this->c[i] = -g[gIdx].f(xvec.c);
    }

    rho = rho / 1.2;

  }


  void Interior_Point_Method::free_Matrix_Vector(Matrix *Mhat, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec, Vector &drvec, Vector &gradvec)
  {

    xvec.clear();
    zvec.clear();
    uvec.clear();
    lambdavec.clear();
    drvec.clear();
    gradvec.clear();

    //delete drvec;
    //delete gradvec;

    if(Mhat->Inv == NULL) cout << "asf" << endl;

    //Mhat->Inv->clear();
    Mhat->clear();

    //delete Mhat->Inv;
    delete Mhat;

  }

  double Interior_Point_Method::getalpha(Vector &drvec, Vector &xvec, Vector &zvec, Vector &uvec, Vector &lambdavec)
  {
    double alpha = 1.0;
    double fact = 0.8;
    int ix = 0, iz = 0, iu = 0, il = 0;

    for(int i = 0; i < nRank; ++i)
    {
      if(i < xRank && xvec.c[ix] != 0.0)
      {
	//if(xvec.c[ix] + drvec.c[i] < 0.0) alpha = min(max(fabs(xvec.c[ix] / drvec.c[i]) *fact, 0.0), alpha);
	++ix;
      }
      else if(i < xRank + zRank )
      {
	if(zvec.c[iz] + drvec.c[i] < 0.0) alpha = min(fabs(zvec.c[iz] / drvec.c[i]) * fact, alpha);
	++iz;
      }
      else if(i < xRank + zRank + uRank)
      {
	if(uvec.c[iu] + drvec.c[i] < 0.0) alpha = min(fabs(uvec.c[iu] / drvec.c[i]) * fact, alpha);
	++iu;
      }
      else if(i < xRank + zRank + uRank + lambdaRank)
      {
	if(lambdavec.c[il] + drvec.c[i] < 0.0) alpha = min(fabs(lambdavec.c[il] / drvec.c[i]) * fact, alpha);
	++il;
      }
    }

    return alpha;
  }

  void checkfunc(Matrix *Mhat)
  {

    const int nrank = Mhat->getrank(); 
    Matrix I(nRank);
    I.Cross(Mhat, Mhat->Inv);
    bool flag = false;

    for(int i = 0; i < nrank; ++i)
    {
      for(int j = 0; j < nrank; ++j)
      {
	double tmp = I.getComponent(i, j) < 1.0e-8 ? 0.0 : 1.0;
	if(tmp > 0.0001) flag = true;

	if(i != j && tmp == 1.0)cout << "warning " << I.getComponent(i, j) << endl;
      }
    }
    I.clear();

    if(flag) return ;
    cout << "A = " << endl;
    for(int i = 0; i < nrank; ++i)
    {
      for(int j = 0; j < nrank; ++j)
      {
	cout  << I.getComponent(i, j) << " " ;
      }
      cout << endl;
    }

    cout << "Ainv = " << endl;

    for(int i = 0; i < nrank; ++i)
    {
      for(int j = 0; j < nrank; ++j)
      {
	cout << I.getComponent(i, j) << " " ;
      }
      cout << endl;
    }

  }

  void Interior_Point_Method::SetObjective(string &str)
  {
    Function *Fobj = getObjectiveFunction();
    Fobj->setfunctype(object, 0);
  }

  void Interior_Point_Method::AddContraint(string &str)
  {
    vector<Function> *g = getConstraintFunction();
    Function tmpFunction;
    tmpFunction.setfunctype(constraint, g->size());
    g->push_back(tmpFunction);
  }
};


