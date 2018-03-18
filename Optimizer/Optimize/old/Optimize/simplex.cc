#include <iostream>
#include <fstream>
#include <vector>
#include "simplex.h"
#include "params.h"
#include "function.h"
#include <math.h>

using namespace std;

double sq(double val){return val*val;}

namespace OptimizeName 
{
  /* --- Function --- */
  void Simplex::initialize()
  {
    //次数LP緩和
    for(int i = 0; i < Num; ++i)
    {
      Function tmp1, tmp2;
      g.push_back(tmp1);
      g.push_back(tmp2);
    } 

    for(int i = 0; i < Nx; ++i)
    {
      Function tmp;
      g.push_back(tmp);
    }
   
    int xpivod = 0;
    double time1 = 0.0; 
    clock_t t1;
    double dddx = 0.0, aaa = 0.0;

    int *cc = new int[Nx];
    for(int i = 0; i < Nx; ++i) cc[i] = 0;
    
    while(1)
    {

      xpivod = rand() % Nx;
      
      t1 = clock();
      for(int ic = 0; ic < g.size(); ++ic)
      {
	aaa = g.at(ic).Coeff[xpivod].val;
	if(aaa >= 0)
	{
	  continue;
	}
	dddx = -g.at(ic).constTermVal;
      }
      t1 = clock() - t1;
      time1 = (double) t1 / CLOCKS_PER_SEC;
      ++cc[xpivod];
      cout << cc[xpivod] << "  " << time1 << endl;
      getchar();
    }

 
    return;

  }

  void Simplex::free()
  {
    Function *pfobj = getObjectiveFunction();
    vector<Function> *B = getConstraintFunction();
    const int Nconstraint = B->size(); 

    delete [] pfobj->Coeff;
    delete pfobj;
    for(int i = 0; i < Nconstraint; ++i)
      delete [] B->at(i).Coeff;
    B->clear();
  }

  void Simplex::opt()
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pB = this->getConstraintFunction();
    const int Nconstraint = pB->size();
    int xpivod = 0;
    int xpivodrow = 0;
    double time1 = 0.0; 
    clock_t t1 = clock();
    //double dxpivodrow = 1.0e6, dddx, aaa;
    double dddx = 0.0, aaa = 0.0;
    
    //this->getFeasibleInitilize();
    
    while(1)
    {
      //xpivod = pfobj->getXpivodCol(*pB);
      cout << xpivod << endl;
      if(xpivod < 0) xpivod = 0;
      xpivodrow = 0;

      t1 = clock();
      for(int ic = 0; ic < Nconstraint; ++ic)
      {
	aaa = g.at(ic).Coeff[xpivod].val;
	if(aaa >= 0)
	{
	  continue;
	}
	dddx = -g.at(ic).constTermVal;
      }
      t1 = clock() - t1;
      time1 = (double) t1 / CLOCKS_PER_SEC;
      cout << time1 << endl;
      getchar();
    }

  }

  void Function::SwapVariable(int xpivod, int ix, Function &Bi)
  {
    double a = this->Coeff[xpivod].val; // apivod
    if(a == 0.0){return ;}

    this->Coeff[xpivod].ix = ix; //新しい基底変数に置き換える

    for(coeff *aCoeff = Bi.firstCoeff; aCoeff != NULL; aCoeff = aCoeff->nextCoeff)
    //for(int i = 0; i < Nx + 1; ++i)
    {
      //coeff *aCoeff = &(Bi.Coeff[i]);
      if(aCoeff->val == 0.0) continue;
      int index = aCoeff->index; // index of array 
      if(this->Coeff[index].val == 0.0)
      {
	coeff *CoeffNew = &this->Coeff[index];
	if(fabs(a * aCoeff->val) < eps){CoeffNew->val = 0.0;continue;}
	
	CoeffNew->index = aCoeff->index;
	CoeffNew->ix = aCoeff->ix;
	CoeffNew->val = a * aCoeff->val;
	CoeffNew->nextCoeff = this->firstCoeff;
	this->firstCoeff->preCoeff = CoeffNew;
	this->firstCoeff = CoeffNew;
	this->firstCoeff->preCoeff = NULL;
      }
      else 
      {
	this->Coeff[index].val = a * aCoeff->val 
	  + (index != xpivod ? this->Coeff[index].val : 0);  

	if(fabs(this->Coeff[index].val) < eps) this->Coeff[index].val = 0.0;

	if(this->Coeff[index].val == 0.0)
	{
	  if(this->firstCoeff != &this->Coeff[index])
	  {
	    this->Coeff[index].preCoeff->nextCoeff = this->Coeff[index].nextCoeff;
	    if(this->Coeff[index].nextCoeff != NULL) 
	      this->Coeff[index].nextCoeff->preCoeff = this->Coeff[index].preCoeff;
	  }
	  else 
	  {
	    this->firstCoeff = this->Coeff[index].nextCoeff;
	    this->firstCoeff->preCoeff = NULL;
	  }
	}

	  
      }
    }

    this->constTermVal = a * Bi.constTermVal + this->constTermVal;
    if(fabs(this->constTermVal) < eps) this->constTermVal = 0.0;

    return ;
  }

  int Function::getXpivodCol(vector<Function> &pB)
  {
    int xpivod = -1;
    double maxval = 0.0;
    //hoge //koko?
    //for(int i = 0; i < Nx + 1; ++i)
    for(coeff *aCoeff = this->firstCoeff; aCoeff != NULL; aCoeff = aCoeff->nextCoeff)
    {
      //coeff *aCoeff = &(this->Coeff[i]);
      if(aCoeff->val == 0.0) continue;
      if((this->getobject() == maximization && aCoeff->val > 0 && aCoeff->index >= 0) ||
	  (this->getobject() == minimization  && aCoeff->val < 0 && aCoeff->index >= 0))
      {
	if((maxval < aCoeff->val && this->getobject() == maximization) || 
	   (maxval > aCoeff->val && this->getobject() == minimization) ) 
	{
	  for(int ic = 0; ic < pB.size(); ++ic)
	  {
	    if(pB.at(ic).Coeff[aCoeff->index].val < 0 && fabs(pB.at(ic).Coeff[aCoeff->index].val) > eps)
	    {
	      maxval = aCoeff->val;
	      xpivod = aCoeff->index;
	      break;
	    }
	  }
	}
      }
    }

    return xpivod;
  }

  int Function::getXpivodRow(int &xpivod, vector<Function> &B)
  {
    int xpivodrow = 1000000;
    double dxpivodrow = 1.0e6;
    double b  = 0.0, dx = 0.0;
    const int Nconstraint = B.size();
    bool artificialFlag = false;

    for(int i = 0; i < Nconstraint; ++i)
    {
      b = B[i].constTermVal;
      if(b < 0.0) artificialFlag = true;
      if(artificialFlag) break;
    }

    dxpivodrow = artificialFlag ? -10000000 : 10000000;

    for(int i = 0; i < Nconstraint; ++i)
    {
      Function *gfunc = &B[i];
      double a = gfunc->Coeff[xpivod].val;
      if(a < 0 && fabs(a) > eps)
      {
	dx = -gfunc->constTermVal / a;
      }
      else 
      {
	dx = 1.0e+30;
      }

      if(!artificialFlag) 
      {
	if(dx >= 0 && dxpivodrow > dx)
	{
	  xpivodrow = i;
	  dxpivodrow = dx; 
	}
      }
      else 
      {
	if(dx >= 0 && dxpivodrow < dx)
	{
	  xpivodrow = i;
	  dxpivodrow = dx;
	}
      }
    }

    return xpivodrow;
  }

  void Function::migration(int xpivodcol)
  {
    int ix = this->ix;
    double a = -this->Coeff[xpivodcol].val;
    
    this->ix = this->Coeff[xpivodcol].ix;
    this->Coeff[xpivodcol].val = -1.0;
    this->Coeff[xpivodcol].ix = ix;
    
    int counter = 0;
    for(coeff *aCoeff = this->firstCoeff; aCoeff != NULL; aCoeff = aCoeff->nextCoeff)
    //for(int i = 0; i < Nx + 1; ++i)
    {
      //coeff *aCoeff = &(this->Coeff[i]);
      aCoeff->val /= a;
      ++counter;
    }

    this->constTermVal /= a;

  }

  void Simplex::getFeasibleInitilize()
  {
    cout << "getF" << endl;
    Function *fobj = new Function(); 
    Function *pfobj = getObjectiveFunction();
    setArtificialVariable(fobj);

    vector<Function> *pB = this->getConstraintFunction();
    int artitermIdx = Nx; // 人工変数が格納されている配列インデックス
    const int Nconstraint = pB->size();
    int xpivod = artitermIdx; 
    int xpivodrow = fobj->getXpivodRow(xpivod, *pB);
    int ix = pB->at(xpivodrow).ix;
    artitermIdx = xpivod;

#if 1
    pB->at(xpivodrow).migration(xpivod);

    fobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
    pfobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      if(ic == xpivodrow) continue;
      pB->at(ic).SwapVariable(xpivod, ix ,pB->at(xpivodrow));
    }
    
    cout << "start " << endl;
    while(xpivod >= 0)
    {
      xpivod = fobj->getXpivodCol(*pB);
    
      for(int i = 0 ; i < Nx; ++i)
      {
      double dbl = pfobj->Coeff[i].val;
      }
      if(xpivod == -1) break;

      xpivodrow = fobj->getXpivodRow(xpivod, *pB);
      ix = pB->at(xpivodrow).ix;
   
      pB->at(xpivodrow).migration(xpivod);
      if(pB->at(xpivodrow).Coeff[xpivod].val == 0.0){continue;}

      if(ix == Nx + Nconstraint) artitermIdx = xpivod; 
      // y = const + a1*x1 + ... + apivod * xpivod + .. 
      // Bxpivodrow = x_ix = bxpivodrow - a1,xpivodrow * x1 - ...   
      fobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
      pfobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
      for(int ic = 0; ic < Nconstraint; ++ic)
      {
	if(ic == xpivodrow) continue;
	pB->at(ic).SwapVariable(xpivod, ix ,pB->at(xpivodrow));
      }
    }
    
    pfobj->Coeff[artitermIdx].val = 0.0;
    if(pfobj->Coeff[artitermIdx].preCoeff == NULL)
    {  
      pfobj->firstCoeff = pfobj->Coeff[artitermIdx].nextCoeff;
      pfobj->firstCoeff->preCoeff = NULL;
    }
    else 
    {
      pfobj->Coeff[artitermIdx].preCoeff->nextCoeff = pfobj->Coeff[artitermIdx].nextCoeff;
      if(pfobj->Coeff[artitermIdx].nextCoeff != NULL)
	pfobj->Coeff[artitermIdx].nextCoeff->preCoeff = pfobj->Coeff[artitermIdx].preCoeff;
    }

    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      if(pB->at(ic).Coeff[artitermIdx].val == 0.0) continue;
      pB->at(ic).Coeff[artitermIdx].val = 0.0;
      if(pB->at(ic).Coeff[artitermIdx].preCoeff == NULL)
      {
	pB->at(ic).firstCoeff = pB->at(ic).Coeff[artitermIdx].nextCoeff;
	pB->at(ic).firstCoeff->preCoeff = NULL;
      }
      else 
      {
	pB->at(ic).Coeff[artitermIdx].preCoeff->nextCoeff = pB->at(ic).Coeff[artitermIdx].nextCoeff;
	if(pB->at(ic).Coeff[artitermIdx].nextCoeff != NULL )  
	  pB->at(ic).Coeff[artitermIdx].nextCoeff->preCoeff = pB->at(ic).Coeff[artitermIdx].preCoeff;
      }
    }
#endif
    delete [] fobj->Coeff;
    delete fobj;

    cout << "get Feasible Initialize " << endl;
  }

  void Simplex::setArtificialVariable(Function *pfobj)
  {
    vector<Function> *pB = this->getConstraintFunction();
    const int artitermIdx = Nx ;
    const int Nconstraint = pB->size();

    pfobj->Coeff[Nx].val = -1;
    pfobj->conectChain();
    pfobj->Coeff[artitermIdx].ix = Nx + pB->size();
    pfobj->Coeff[artitermIdx].index = Nx; 
    pfobj->setobject(maximization);

    for(int i = 0; i < Nconstraint; ++i)
    {
      Function *constFunc = &pB->at(i);
      constFunc->Coeff[Nx].ix = Nx + pB->size();
      constFunc->Coeff[artitermIdx].index = Nx;
      if(constFunc->constTermVal >= 0.0) continue;
      constFunc->Coeff[artitermIdx].val = 1.0; 
      
      coeff *tmpC = &(constFunc->Coeff[artitermIdx]);
      tmpC->nextCoeff = constFunc->firstCoeff;

      if(constFunc->firstCoeff != NULL)
      {
	constFunc->firstCoeff->preCoeff = tmpC;
      }
      constFunc->firstCoeff = tmpC;
     
    }
  }

};
