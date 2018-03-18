#include <iostream>
#include <fstream>
#include <vector>
#include "simplex.h"
#include "params.h"
#include "function.h"

using namespace std;

namespace OptimizeName 
{
  /* --- Function --- */
  void Simplex::run()
  {
    initialize();

    opt();

    free();
  }

  void Simplex::initialize()
  {

    string strFunc = "-4x0 - 3x1 ";
    SetObjective(strFunc);

    strFunc = "x0 - x1  < 1"; 
    AddConstraint(strFunc); 

    strFunc = "2x0 -x1 < 3";
    AddConstraint(strFunc); 

    strFunc = "x1 < 5";
    AddConstraint(strFunc);
    return;
#if 0
    string strFunc = "-x0 - 2x1 - x2 - x3";
    //string strFunc = "x0 + 2x1 + x2 + x3";
    SetObjective(strFunc);

    strFunc = "x0 - 2x1 + 2x2 + 3x3 < 3"; 
    AddConstraint(strFunc); 

    strFunc = "3x0 +4x1-5x2-6x3 < 6";
    AddConstraint(strFunc); 

    strFunc = "3x0+x1-5x2-7x3 < 2";
    AddConstraint(strFunc);

    strFunc = "-2x0 + x1 + 2x2 + x3 < 2";
    AddConstraint(strFunc); 

    strFunc = "3x0 - 3x1 + 2x2 + 2x3 < 5";
    AddConstraint(strFunc); 
    /*
       string strFunc = "2x0 + x1 -x2";
       SetObjective(strFunc);

       strFunc = "x0 - 2x1 +2x2 < -4";
       AddConstraint(strFunc);

       strFunc = "-2x1 - x2 < -1";
       AddConstraint(strFunc);

       strFunc = "2x0-3x1 +x2 < 1";
       AddConstraint(strFunc);
     */

#endif
  }

  void Simplex::free()
  {
    Function *pfobj = getObjectiveFunction();
    vector<Function> *B = getConstraintFunction();
    const int Nconstraint = B->size(); 

    //delete [] pfobj->Coeff;
    pfobj->clearCoeff();
    delete pfobj;
    for(int i = 0; i < Nconstraint; ++i)
      //delete [] B->at(i).Coeff;
      B->at(i).clearCoeff();
    B->clear();
  }

  void Simplex::SetObjective(string &str)
  {
    Function *Fobj = getObjectiveFunction();
    Fobj->Add(str);
    //Fobj->conectChain();
  }

  void Simplex::AddConstraint(string &str)
  {
    vector<Function> *B = getConstraintFunction();
    Function tmp;
    tmp.Add(str);
    //tmp.conectChain();
    tmp.ix = B->size() + Nx;
    B->push_back(tmp);
  }

  void Simplex::FuncLog()
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pB = this->getConstraintFunction();
    const int Nconstraint = pB->size();
    cout << "fobj = ";pfobj->CheckFunc();cout << endl;
    for(int ic = 0; ic < Nconstraint*0+2; ++ic)
    {
      cout << "B"  << pB->at(ic).ix << "= ";pB->at(ic).CheckFunc();cout << endl;
    }
  }

  void Simplex::opt()
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pB = this->getConstraintFunction();
    const int Nconstraint = pB->size();
    int xpivodrow = 0, xpivod = 0, *basicVariableIdx = new int[Nconstraint], *nonbasicVariableIdx = new int[Nx];

    this->getFeasibleInitilize(basicVariableIdx, nonbasicVariableIdx);

    clock_t time = clock();
    while(1)
    {
      xpivod = pfobj->getXpivodCol(*pB);
      if(xpivod < 0) break;

      xpivodrow  = pfobj->getXpivodRow(xpivod, *pB);
      int ix = pB->at(xpivodrow).ix;
      pB->at(xpivodrow).migration(xpivod);

      if(pB->at(xpivodrow).findCoeff(xpivod)){continue;}

      // y = const + a1*x1 + ... + apivod * xpivod + .. 
      // Bxpivodrow = x_ix = bxpivodrow - a1,xpivodrow * x1 - ...   
      pfobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
      for(int ic = 0; ic < Nconstraint; ++ic)
      {
	if(ic == xpivodrow) continue;
	pB->at(ic).SwapVariable(xpivod, ix ,pB->at(xpivodrow));
      }
      cout << pfobj->constTermVal << " " << xpivod << endl;
    }
    time = clock() - time;
    cout << (double) time / CLOCKS_PER_SEC << endl;
    getchar();

    this->FuncLog();
    cout << pfobj->constTermVal << endl;

    delete [] basicVariableIdx;
    delete [] nonbasicVariableIdx;
  }

  void Simplex::optRevised()
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pg = this->getConstraintFunction();
    const int Nconstraint = pg->size();
    int *basicVariableIdx = new int[Nconstraint], *nonbasicVariableIdx = new int[Nx];
    Vector *cB = new Vector(Nconstraint), *cN = new Vector(Nx);
    Vector *b = new Vector(Nconstraint);
    Vector *xB = new Vector(Nconstraint);
    Vector *xN = new Vector(Nx);
    Matrix *A = new Matrix(Nconstraint, Nx + Nconstraint);
    Matrix *B = new Matrix(Nconstraint);
    Matrix *N = new Matrix(Nconstraint, Nx);
    Vector oldVec(Nconstraint);

    B->AllocTraceMatrix();
    N->AllocTraceMatrix();

    //1.初期実行可能基底解のインデックスを取得
    this->getFeasibleInitilize(basicVariableIdx, nonbasicVariableIdx);
    this->setConstraintVector(b);
    this->setConstraintMatrix(A);

    for(int i = 0; i < Nconstraint; ++i)
      cB->setComponent(i, 0.0);

    int ii = 0;
    for(map<int, double>::iterator ite = pfobj->CoeffMap.begin(); ite != pfobj->CoeffMap.end() ; ++ite)
    {
      cN->setComponent(ii, ite->second);
      ++ii;
    }
    for(; ii < Nx; ++ii)
      cN->setComponent(ii, 0.0); 
    
    int count = 0;
    double dblinit = pfobj->constTermVal;
    clock_t time = clock();

    this->setMatrixes(basicVariableIdx, nonbasicVariableIdx, A, B, N);

#if 1
    /////////////////////////////////
    AbstVector zsN(Nx), xBs(Nconstraint);
    AbstVector dxB(Nconstraint), aj(Nconstraint);
    AbstVector dzN(Nx), uvec(Nconstraint), ei(Nconstraint);
    AbstMatrix *Bt0 = B->getTraceMatrix(), *Nt0 = N->getTraceMatrix(); 

    for(int i = 0; i < Nconstraint; ++i){xBs.setComponent(i, b->getComponent(i));}
    for(int i = 0; i < Nx; ++i){zsN.setComponent(i, -cN->getComponent(i));}

    dxB.calcLUMatrixForUnitMatrix();
    uvec.calcLUMatrixForUnitMatrix();
    
    B->setNonZeroNonZeroFirstRow();
    N->setNonZeroNonZeroFirstRow();
    Bt0->setNonZeroNonZeroFirstRow();
    Nt0->setNonZeroNonZeroFirstRow();

    time = clock();

    while(1)
    {
      int i0, j;
      double t;

      clock_t timeaLL  = clock();
      Vector oldVec2(Nconstraint);
      
      clock_t time1 = clock();
      //step1
      j = getzsNindex(zsN);
      if( j < 0 ){ break; }

      aj.cmap.clear();
      //for(int i = 0; i < Nconstraint; ++i){aj.setComponent(i, N->getComponent(i, j));}
      for(int i = 0; i < Nconstraint; ++i) aj.setComponent(i, 0.0);
#if 1  
      int in = Nt0->firstNonZeroRowIndex[j];
      map<int, ComponentClass>::iterator iteN = Nt0->compmap2.find(in + j * Nt0->getcolrank());
      while(1)
      {
	if(iteN == Nt0->compmap2.end()) break;
	if(j < iteN->second.i) break;
	aj.setComponent(iteN->second.j, iteN->second.val);
	++iteN;
      }
#endif
      
      time1 = clock() - time1;
      cout << "time1  = " << (double) time1 / CLOCKS_PER_SEC << endl; 

      //step2
      clock_t time2 = clock();
      if(count == 0)
      {
	//dxB.LUdecomposition(B, aj);
	dxB.cmap.clear();
	for(int i = 0; i < dxB.getrank(); ++i){ dxB.setComponent(i, aj.getComponent(i));}
      }
      else dxB.LUdecompositionReusing(B, aj);
      time2 = clock() - time2;
      cout << "time2  = " << (double) time2 / CLOCKS_PER_SEC << endl; 

      //step3
      clock_t time3 = clock();
      i0 = getbasicindex(t, dxB, xBs); 

      for(int i = 0; i < Nconstraint; ++i){ei.setComponent(i, i == i0 ? -1 : 0);}
      time3 = clock() - time3;
      cout << "time3  = " << (double) time3 / CLOCKS_PER_SEC << endl; 
      //step4
      clock_t time4 = clock();
      if(count == 0) 
      {
	//uvec.LUdecomposition(Bt0, ei);
	for(int i = 0; i < uvec.getrank(); ++i){ uvec.setComponent(i, ei.getComponent(i));}
      }	
      else uvec.LUdecompositionReusingForTrace(Bt0, ei);
      //uvec.setNonZero();
      dzN.Cross(Nt0, &uvec); 
      time4 = clock() - time4;
      cout << "time4  = " << (double) time4 / CLOCKS_PER_SEC << endl; 
      
      //step5
      clock_t time5 = clock();
      double s = zsN.getComponent(j) / dzN.getComponent(j);

      if(t != 0.0) 
      {
	for(map<int, double>::iterator ite = dxB.cmap.begin(); ite != dxB.cmap.end(); ++ite)
	  xBs.setComponent(ite->first, xBs.getComponent(ite->first) - t * ite->second);
      }
      if(s != 0.0)
      {
	for(map<int, double>::iterator ite = dzN.cmap.begin(); ite != dzN.cmap.end(); ++ite)
	  zsN.setComponent(ite->first, zsN.getComponent(ite->first) - s * ite->second);
      }

      xBs.setComponent(i0, t);
      zsN.setComponent(j, s);
      time5 = clock() - time5;
     cout << "time5  = " << (double) time5 / CLOCKS_PER_SEC << endl; 
      //step6
      clock_t time6 = clock();
      int rval = swapVariableRevised_r2(i0, j, basicVariableIdx, nonbasicVariableIdx, cB, cN, B, N, Bt0, Nt0);
      time6 = clock() - time6;

      cout << "time6  = " << (double) time6 / CLOCKS_PER_SEC << endl; 
      if(rval != 0) break;

      //step7 
      clock_t time7 = clock();
      oldVec2.cmap.clear();
      for(map<int, double>::iterator ite = dxB.cmap.begin(); ite != dxB.cmap.end(); ++ite)
      {
	oldVec2.setComponent(ite->first, ite->second); 
      }

      dxB.oldx.push_back(oldVec2); 
      dxB.oldIndex.push_back(i0);

      uvec.oldx.push_back(oldVec2);
      uvec.oldIndex.push_back(i0);
      
      time7 = clock() - time7;

      //cout << xBs * *cB << endl;
      timeaLL = clock() - timeaLL;
      
      cout << "time7  = " << (double) time7 / CLOCKS_PER_SEC << endl; 
      cout << "timeall  = " << (double) timeaLL / CLOCKS_PER_SEC << endl; 
      ++count;
    }
    cout << "opt val = " << xBs * *cB + dblinit << endl;
    time = clock()-time;
    cout << "time " << (double) time / CLOCKS_PER_SEC << endl;
    getchar();
    ///////////////////////////////////////////
#endif

    xB->clear();
    xN->clear();
    b->clear();
    A->clear();
    B->clear();
    N->clear();

    delete [] basicVariableIdx;
    delete [] nonbasicVariableIdx;
    delete b;
    delete xB;
    delete xN;
    delete A;
    delete B;
    delete N;
  }

  void Simplex::setConstraintVector(Vector *tmpVec)
  {
    vector<Function> *pg = this->getConstraintFunction();
    const int Nconstraint = pg->size();

    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      tmpVec->setComponent(ic, pg->at(ic).constTermVal); 
    }
  }

  void Simplex::setConstraintMatrix(Matrix *A)
  {
    vector<Function> *pg = this->getConstraintFunction();
    const int Nconstraint = pg->size();
    int i = 0, j = 0;

    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      Function *gi = &(pg->at(ic));

      //for(coeff *Coeff = gi->firstCoeff; Coeff != NULL; Coeff = Coeff->nextCoeff)
      for(map<int, double>::iterator ite = gi->CoeffMap.begin(); ite != gi->CoeffMap.end(); ++ite)
      {
	//	Coeff->val *= -1;
	ite->second *= -1;
	// xj = b - cixi -> b = xj + cixi;
      }
      
      gi->setCoeff(gi->ix, 1.0);
    }

    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      Function gi = pg->at(ic); 
      for(map<int, double>::iterator ite = gi.CoeffMap.begin(); ite != gi.CoeffMap.end(); ++ite)
      {
	i = ic;
	//j = Coeff->ix;
	//A->setComponent(i, j, Coeff->val);//移項した値を代入
	j = ite->first;
	if(ite->second == 0.0) continue;
	  A->setComponent(i, j, ite->second);
      }
    }

  }

  void Simplex::setMatrixes(int *basicVariableIdx, int *nonbasicVariableIdx, Matrix *A, Matrix *B, Matrix *N)
  {
    vector<Function> *pg = this->getConstraintFunction();
    const int Nconstraint = pg->size();
    int idx = 0;

    // set BasicMatrix
    for(int i = 0; i < Nconstraint; ++i)
    {
      for(int j  = 0; j < Nconstraint; ++j)
      {
	idx = basicVariableIdx[j];
	B->setComponent(i, j, A->getComponent(i, idx));
      }
    } 

    // set NonBasicMatrix
    for(int i = 0; i < Nconstraint; ++i)
    {
      for(int j = 0; j < Nx; ++j)
      {
	idx = nonbasicVariableIdx[j];
	N->setComponent(i, j, A->getComponent(i, idx));
      }
    }
  }

  int Function::SwapVariable(int xpivod, int ix, Function &Bi)
  {

    double a = this->getCoeff(xpivod);
    if(a == 0.0){return 0;}

    this->eraseCoeff(xpivod);

    for(map<int, double>::iterator ite = Bi.CoeffMap.begin(); ite != Bi.CoeffMap.end(); ++ite)
    {
      int index = ite->first; 

      if(!this->findCoeff(index))
      {
	this->setCoeff(index, a * ite->second);
      }
      else
      {
	if(index != xpivod)  
	{
	  this->CoeffMap[index] = a * ite->second + this->CoeffMap[index]; 
	}
	else
	{
	  this->CoeffMap[index] = a * ite->second;
	}

	if(fabs(this->CoeffMap[index]) < eps){this->setCoeff(index, 0.0); this->eraseCoeff(index);}

      }
    } 

    this->constTermVal = a * Bi.constTermVal + this->constTermVal;
    if(fabs(this->constTermVal) < eps) this->constTermVal = 0.0;

    return 0;
  }

  int Function::getXpivodCol(vector<Function> &pB)
  {

    int ix = -1;
    int xpivod = -1;
    double maxval = 0.0;
    //double minval = 1.0e+30;

#if 1
    for(map<int, double>::iterator ite = this->CoeffMap.begin(); ite != this->CoeffMap.end(); ++ite)
    {
      if((this->getobject() == maximization && ite->second > 0 && maxval < ite->second) ||
	  (this->getobject() == minimization && ite->second < 0 && maxval > ite->second))
      {

	maxval = ite->second; 
	ix = ite->first;
	xpivod = ite->first;
      }
    }	
#endif

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
      double a = 0.0;
      if(gfunc->findCoeff(xpivod))
	a = gfunc->getCoeff(xpivod);
      if((!artificialFlag && a < 0 && fabs(a) > eps) || (artificialFlag && a > 0))
	//if(-gfunc->constTermVal / a >= 0 && fabs(a) > eps)
      {
	dx = -gfunc->constTermVal / a;
	//if(artificialFlag) 
	//cout << i << " " << a <<  "  " << dxpivodrow << endl;
      }
      else 
      {
	dx = artificialFlag ? -1.0e+30 : 1.0e+30;
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

  void Function::CheckFunc()
  {
    int count = 0;
    cout << this->constTermVal << " + ";
    for(map<int, double>::iterator ite = this->CoeffMap.begin(); ite != this->CoeffMap.end(); ++ite)
    {
      cout << ite->second<<  "*x" << ite->first << " + ";
      ++count;
      if(count == 2) break;
    }
  }


  void Function::migration(int xpivodcol)
  {
    int ix = this->ix;
    double a = -this->getCoeff(xpivodcol);

    for(map<int, double>::iterator ite = this->CoeffMap.begin(); ite != this->CoeffMap.end(); ++ite)
    {
      if( ite->second / a > 1.0e+6)
      {
	this->eraseCoeff(xpivodcol);	
	return ;
      }
    }

    map<int, double>::iterator itecol = this->CoeffMap.find(xpivodcol);
    this->ix = itecol->first; 
    this->setCoeff(ix, -1.0);
    this->eraseCoeff(itecol->first);

    for(map<int, double>::iterator ite = this->CoeffMap.begin(); ite != this->CoeffMap.end(); ++ite)
    {
      ite->second /= a;
    }
    this->constTermVal /= a;

    return;
  }

  void Simplex::getFeasibleInitilize(int *basicVariableIdx, int *nonbasicVariableIdx)
  {

    cout << "start getFeasible Ini" << endl;
    Function *fobj = new Function(); 
    Function *pfobj = getObjectiveFunction();
    int rc = setArtificialVariable(fobj);

    //if(rc != 0){delete fobj;return;}

    vector<Function> *pB = this->getConstraintFunction();
    const int Nconstraint = pB->size();
    int artitermIdx = Nx + Nconstraint; // 人工変数が格納されている配列インデックス
    int xpivod = Nx + Nconstraint; // 人工変数の変数インデックス
    int xpivodrow = fobj->getXpivodRow(xpivod, *pB);

    this->FuncLog();

    if(xpivodrow > Nconstraint || rc != 0) 
    {
      for(int ic = 0; ic < Nconstraint; ++ic)
      {
	basicVariableIdx[ic] = pB->at(ic).ix;
      }

      int idx = 0;
      for(map<int, double>::iterator ite = pfobj->CoeffMap.begin(); ite != pfobj->CoeffMap.end(); ++ite)
      {
	nonbasicVariableIdx[idx] = ite->first; 
	++idx;
      }
      for(; idx < Nx ; ++idx) nonbasicVariableIdx[idx] = 0;
      //for(int idx = 0; idx < Nx; ++idx)
      //{
      //nonbasicVariableIdx[idx] = fobj->Coeff[idx].ix;
      //nonbasicVariableIdx[idx] = fobj->Coeff
      //}

      delete fobj;
      return ;
    } 

    cout << "get " << endl;
    int ix = pB->at(xpivodrow).ix;
    pB->at(xpivodrow).migration(xpivod);

    fobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
    pfobj->SwapVariable(xpivod, ix, pB->at(xpivodrow));
    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      if(ic == xpivodrow) continue;
      pB->at(ic).SwapVariable(xpivod, ix ,pB->at(xpivodrow));
    }

    while(xpivod >= 0)
    {
      xpivod = fobj->getXpivodCol(*pB);

      if(xpivod == -1) break;

      xpivodrow = fobj->getXpivodRow(xpivod, *pB);
      ix = pB->at(xpivodrow).ix;

      pB->at(xpivodrow).migration(xpivod);

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


    cout << "get " << endl;
    checkBasicParams(fobj, artitermIdx);

    pfobj->eraseCoeff(Nx + Nconstraint);
    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      if(!pB->at(ic).findCoeff(Nx + Nconstraint)) continue;
      pB->at(ic).eraseCoeff(Nx + Nconstraint); 
    }

    cout << "get Feasible Initialize " << endl;

    fobj->CheckFunc();

    delete fobj;
    return;
  }

  int Simplex::setArtificialVariable(Function *pfobj)
  {
    vector<Function> *pB = this->getConstraintFunction();
    const int artiterm = Nx + pB->size();
    const int Nconstraint = pB->size();

    string strNx, strFunc = "";
    vartostring(strNx, artiterm);
    strFunc = "-x" + strNx;

    pfobj->setCoeff(artiterm, -1.0);
    pfobj->setobject(maximization);

    int rc = -1;
    for(int i = 0; i < Nconstraint; ++i)
    {
      Function *constFunc = &pB->at(i);
      if(constFunc->constTermVal >= 0.0) continue;
      constFunc->setCoeff(artiterm, 1.0);
      rc = 0;
    }

    return rc;
  }

  int Simplex::getNonBasicIdx(Vector *cN, AbstMatrix *Nt, Vector *pi)
  {
    int irtn = -1;
    double dbl = 0.0;
    double min = 1.0e+10;
    const int Nconstraint = this->getConstraintFunction()->size();;

    // xi = cN - Nt * pi;
    for(int i = 0; i < Nx; ++i)
    {
      dbl = cN->getComponent(i);
      for(int j = 0; j < Nconstraint; ++j)
      {
	dbl -= Nt->getComponent(i, j) * pi->getComponent(j);  
      }
      if(dbl <= 0.0 && dbl < min)
	//if(dbl <= 0.0)
      {
	irtn = i;
	min = dbl;
	//break;
      }
    }
    return irtn;
  }

  int Simplex::swapVariableRevised(int nonbasicidx, Vector *y, Vector *xB, Vector *xN, Vector *bvar, int *basicVariableIdx, int *nonbasicVariableIdx, Vector *cN, Vector *cB, Matrix *B, Matrix *N)
  {
    int k = -1, tmpint, rval;
    const int nVec = y->getrank(); 
    double theta = INT_MAX;
    double tmpdbl = 0.0;

    //シータとkを計算する
    for(int i = 0; i < nVec; ++i)
    {
      tmpdbl = bvar->getComponent(i) / y->getComponent(i);
      if(theta > tmpdbl && y->getComponent(i) > 0.0)
      {
	theta = bvar->getComponent(i) / y->getComponent(i);
	k = i;
      }
    }

    for(int i = 0; i < nVec; ++i)
    {
      tmpdbl = bvar->getComponent(i) - theta * y->getComponent(i);
      if(i == k) tmpdbl = theta;
      xB->setComponent(i, tmpdbl);
    }

    for(int i = 0; i < xN->getrank(); ++i)
    {
      xN->setComponent(i, 0.0);
    }

    tmpint = nonbasicVariableIdx[nonbasicidx];
    nonbasicVariableIdx[nonbasicidx] = basicVariableIdx[k];
    basicVariableIdx[k] = tmpint;

    tmpdbl = cN->getComponent(nonbasicidx);
    cN->setComponent(nonbasicidx, cB->getComponent(k));
    cB->setComponent(k, tmpdbl);

    //BとNの列を交換
    for(int i = 0; i < nVec; ++i)
    {
      tmpdbl = N->getComponent(i, nonbasicidx);
      N->setComponent(i, nonbasicidx, B->getComponent(i, k));
      B->setComponent(i, k, tmpdbl);
    }

    rval = k;
    return rval;
  }

  //人工変数が非基底変数か
  //もしそうならば、基底変数へ
  void Simplex::checkBasicParams(Function *fobj, int &artitermIdx)
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *gfunc = this->getConstraintFunction();
    const int Nconstraint = gfunc->size();
    bool violationFlag = false; 
    int xpivodrow = -1;

    for(int i = 0; i < gfunc->size(); ++i)
    {
      if(gfunc->at(i).ix == artitermIdx)
      {
	violationFlag = true;
	xpivodrow = i;
	break;
      }
    }

    if(!violationFlag) return; 

    int xpivod = gfunc->at(xpivodrow).CoeffMap.begin()->first;
    int ix = gfunc->at(xpivodrow).CoeffMap.begin()->first;

    gfunc->at(xpivodrow).migration(xpivod);

    fobj->SwapVariable(xpivod, ix, gfunc->at(xpivodrow));
    pfobj->SwapVariable(xpivod, ix, gfunc->at(xpivodrow));
    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      if(ic == xpivodrow) continue;
      gfunc->at(ic).SwapVariable(xpivod, ix ,gfunc->at(xpivodrow));
    }

    artitermIdx = xpivod;
  }

  int Simplex::getzsNindex(AbstVector &vec)
  {
    int rval = -1;
    double minval = 1.0e+30, maxval = -1.0e+30;
    Function *pfobj = getObjectiveFunction();
    eObject object = pfobj->getobject();

    //for(int i = 0; i < nrank; ++i)
    for(map<int, double>::iterator ite = vec.cmap.begin(); ite != vec.cmap.end(); ++ite)
    {
      int i = ite->first;
      double tmpval = vec.getComponent(i);
      if(object == maximization && tmpval < 0.0)
      { 
	if(minval > tmpval)
	{
	  minval = tmpval;
	  rval = i;
	} 
      } else if(object == minimization && tmpval >= 0.0)
      {
	if(maxval < tmpval)
	{
	  maxval = tmpval;
	  rval = i;
	}
      }	

    }
    return rval;
  }

  int Simplex::getbasicindex(double &dbl, AbstVector &dx, AbstVector &b)
  {
    int rval = -1;
    dbl = 1.0e+30;

    for(map<int, double>::iterator ite = dx.cmap.begin(); ite != dx.cmap.end(); ++ite)
    {
      int i = ite->first;
      double tmpdbl = ite->second; 
      if(tmpdbl <= 0) continue;
      if(b.getComponent(i) / tmpdbl < dbl)
      {
        dbl = b.getComponent(i) / tmpdbl; 
	rval = i;
      }	
    }
    
    return rval;
  }

  int Simplex::swapVariableRevised_r2(int Bidx, int Nidx, int *basicVariableIdx, int *nonbasicVariableIdx, Vector *cB, Vector *cN, Matrix *B, Matrix *N, AbstMatrix *Bt0, AbstMatrix *Nt0)
  { 
    int rval = 0;

    int tmpint = nonbasicVariableIdx[Nidx];
    nonbasicVariableIdx[Nidx] = basicVariableIdx[Bidx];
    basicVariableIdx[Bidx] = tmpint;

    double tmpdbl = cN->getComponent(Nidx);
    cN->setComponent(Nidx, cB->getComponent(Bidx));
    cB->setComponent(Bidx, tmpdbl);

    int ib = Bt0->firstNonZeroRowIndex[Bidx];
    int in = Nt0->firstNonZeroRowIndex[Nidx];
    //map<int, double>::iterator iteB = Bt0->compmap.find(ib + Bidx * Bt0->getcolrank());
    //map<int, double>::iterator iteN = Nt0->compmap.find(in + Nidx * Nt0->getcolrank());
    map<int, ComponentClass>::iterator iteB = Bt0->compmap2.find(ib + Bidx * Bt0->getcolrank());
    map<int, ComponentClass>::iterator iteN = Nt0->compmap2.find(in + Nidx * Nt0->getcolrank());
    
    vector<int> bidxvec, nidxvec;
    vector<double> Bcompvec, Ncompvec;

    bidxvec.reserve(100); nidxvec.reserve(100);
    Bcompvec.reserve(100); Ncompvec.reserve(100);
#if 1
    while(1)
    {
      if(iteB == Bt0->compmap2.end()) break;
      //if(Bidx < iteB->first / Bt0->getcolrank()) break; 
      if(Bidx < iteB->second.i) break;
      //bidxvec.push_back(iteB->first % Bt0->getcolrank());
      bidxvec.push_back(iteB->second.j);
      //Bcompvec.push_back(iteB->second);
      Bcompvec.push_back(iteB->second.val);
      ++iteB;
    }

    while(1)
    {
      if(iteN == Nt0->compmap2.end()) break;
      //if(Nidx < iteN->first / Nt0->getcolrank()) break;
      if(Nidx < iteN->second.i) break;
      //nidxvec.push_back(iteN->first % Nt0->getcolrank());
      nidxvec.push_back(iteN->second.j);
      //Ncompvec.push_back(iteN->second);
      Ncompvec.push_back(iteN->second.val);
      ++iteN;
    }
    
    for(int i = 0; i < nidxvec.size(); ++i)
    {
      int irow = nidxvec[i];
      N->setComponent(irow, Nidx, 0.0); 
      Nt0->setComponent(Nidx, irow, 0.0);
    }

    for(int i = 0; i < bidxvec.size(); ++i)
    {
      int irow = bidxvec[i];
      B->setComponent(irow, Bidx, 0.0); 
      Bt0->setComponent(Bidx, irow, 0.0);
    }

    for(int i = 0; i < bidxvec.size(); ++i)
    {
      int irow = bidxvec[i];
      N->setComponent(irow, Nidx, Bcompvec[i]); 
      Nt0->setComponent(Nidx, irow, Bcompvec[i]);
      if(N->firstNonZeroRowIndex[irow] > Nidx ) 
	N->firstNonZeroRowIndex[irow] = Nidx;
    }

    for(int i = 0; i < nidxvec.size(); ++i)
    {
      int irow = nidxvec[i];
      B->setComponent(irow, Bidx, Ncompvec[i]); 
      Bt0->setComponent(Bidx, irow, Ncompvec[i]);
      if(B->firstNonZeroRowIndex[irow] > Bidx)
	B->firstNonZeroRowIndex[irow] = Bidx;
    }

    Bt0->firstNonZeroRowIndex[Bidx] = in;
    Nt0->firstNonZeroRowIndex[Nidx] = ib;
#endif 
    return rval;
  }
};
#if 0
void Simplex::optRevised()
  {
    Function *pfobj = this->getObjectiveFunction();
    vector<Function> *pg = this->getConstraintFunction();
    const int Nconstraint = pg->size();
    int *basicVariableIdx = new int[Nconstraint];
    int *nonbasicVariableIdx = new int[Nx];
    Vector *cB = new Vector(Nconstraint), *cN = new Vector(Nx);
    Vector *b = new Vector(Nconstraint);
    Vector *bvar = new Vector(Nconstraint);
    Vector *xB = new Vector(Nconstraint);
    Vector *xN = new Vector(Nx);
    Vector *pi = new Vector(Nconstraint);
    Matrix *A = new Matrix(Nconstraint, Nx + Nconstraint);
    Matrix *B = new Matrix(Nconstraint);
    Matrix *N = new Matrix(Nconstraint, Nx);
    Vector *y = new Vector(Nconstraint);
    Vector *ai = new Vector(Nconstraint);
    Vector oldVec(Nconstraint);

    B->AllocTraceMatrix();
    N->AllocTraceMatrix();

    //1.初期実行可能基底解のインデックスを取得
    this->getFeasibleInitilize(basicVariableIdx, nonbasicVariableIdx);
    this->setConstraintVector(b);
    this->setConstraintMatrix(A);

    for(int i = 0; i < Nconstraint; ++i)
      cB->setComponent(i, 0.0);

    int ii = 0;
    for(map<int, double>::iterator ite = pfobj->CoeffMap.begin(); ite != pfobj->CoeffMap.end() ; ++ite)
    {
      cN->setComponent(ii, ite->second);
      ++ii;
    }
    for(; ii < Nx; ++ii)
      cN->setComponent(ii, 0.0); 
    /*
       for(int i = 0; i < Nx; ++i)
       {
       if(pfobj->find(i) )
       cN->setComponent(i, p)
    //cN->setComponent(i, pfobj->Coeff[i].val);
    }
     */
    
    int count = 0;
    double dblinit = pfobj->constTermVal;
    clock_t time = clock();
    clock_t time1 = 0;
    clock_t time2 = clock();
    clock_t time3 = clock();


    this->setMatrixes(basicVariableIdx, nonbasicVariableIdx, A, B, N);

    /////////////////////////////////
    AbstVector zsN(Nx), xBs(Nconstraint);
    AbstVector dxB(Nconstraint), aj(Nconstraint);
    AbstVector dzN(Nx), uvec(Nconstraint), ei(Nconstraint);
    AbstMatrix *Bt0 = B->getTraceMatrix(), *Nt0 = N->getTraceMatrix(); 

    for(int i = 0; i < Nconstraint; ++i)
    {
      xBs.setComponent(i, b->getComponent(i));
      if(b->getComponent(i) != 0) xBs.cmap[i] = b->getComponent(i);
    }
    for(int i = 0; i < Nx; ++i)
    {
      zsN.setComponent(i, -cN->getComponent(i));
      if(-cN->getComponent(i) != 0) zsN.cmap[i] = -cN->getComponent(i);
    }

    dxB.calcLUMatrixForUnitMatrix();
    uvec.calcLUMatrixForUnitMatrix();
    
    //B->setNonZero();
    //N->setNonZero();
    time = clock();

    while(1)
    {
      time1 = clock();
      int i0, j;
      double t;

      clock_t time6 = clock();
      Vector oldVec2(Nconstraint);
      
      j = getzsNindex(zsN);
      if( j < 0 ){ break; }

      //step1
      for(int i = 0; i < Nconstraint; ++i){aj.setComponent(i, N->getComponent(i, j));}

      //step2
      time2 = clock();
      if(count == 0)
      {
	//dxB.LUdecomposition(B, aj);
	dxB.cmap.clear();
	for(int i = 0; i < dxB.getrank(); ++i)
	{ 
	  dxB.setComponent(i, aj.getComponent(i));
	  if(dxB.getComponent(i) != 0.0)
	  {
	    dxB.cmap[i] = dxB.getComponent(i);
	  }
	}
      }
      else dxB.LUdecompositionReusing(B, aj);

      time2 = clock() - time2;
      //step3
      i0 = getbasicindex(t, dxB, xBs); 

      for(int i = 0; i < Nconstraint; ++i){ei.setComponent(i, i == i0 ? -1 : 0);}

      //step4
      clock_t time4 = clock();
      if(count == 0) 
      {
	//uvec.LUdecomposition(Bt0, ei);
	for(int i = 0; i < uvec.getrank(); ++i){ uvec.setComponent(i, ei.getComponent(i));}
      }	
      else uvec.LUdecompositionReusingForTrace(Bt0, ei);
      uvec.setNonZero();
      dzN.Cross(Nt0, &uvec); 
      dzN.cmap.clear();
      for(int i = 0; i < Nx;++i){if(dzN.getComponent(i) != 0.0) dzN.cmap[i] = dzN.getComponent(i);}
      //for(int i = 0; i < Nx; ++i){dzN.setComponent(i, -dzN.getComponent(i));}
      time4 = clock() - time4;
      time6 = clock() - time6;
      //step5
      double s = zsN.getComponent(j) / dzN.getComponent(j);

      if(t != 0.0) 
      {
	//for(int i = 0; i < xBs.getrank(); ++i){xBs.setComponent(i, xBs.getComponent(i) - t * dxB.getComponent(i));}
	for(map<int, double>::iterator ite = dxB.cmap.begin(); ite != dxB.cmap.end(); ++ite)
	{
	  int i = ite->first;
	  xBs.setComponent(i, xBs.getComponent(i) - t * ite->second);
	  if(xBs.getComponent(i) == 0.0 && xBs.cmap.find(i) != xBs.cmap.end())
	    xBs.cmap.erase(xBs.cmap.find(i));
	  else if(xBs.getComponent(i) != 0.0)
	    xBs.cmap[i] = xBs.getComponent(i);
	}
      }
      if(s != 0.0)
      {
	//for(int i = 0; i < zsN.getrank(); ++i){zsN.setComponent(i, zsN.getComponent(i) - s * dzN.getComponent(i));}
	for(map<int, double>::iterator ite = dzN.cmap.begin(); ite != dzN.cmap.end(); ++ite)
	{
	  int i = ite->first;
	  zsN.setComponent(i, zsN.getComponent(i) - s * ite->second);
	  if(zsN.getComponent(i) == 0.0 && zsN.cmap.find(i) != zsN.cmap.end())
	    zsN.cmap.erase(zsN.cmap.find(i));
	  else if (zsN.getComponent(i) != 0.0)
	    zsN.cmap[i] = zsN.getComponent(i);
	}
      }

      xBs.setComponent(i0, t);
      zsN.setComponent(j, s);

      //step6
      clock_t time5 = clock();
      int rval = swapVariableRevised_r2(i0, j, basicVariableIdx, nonbasicVariableIdx, cB, cN, B, N, Bt0, Nt0);

      if(rval != 0) break;

      //step7 
      oldVec2.cmap.clear();
      for(map<int, double>::iterator ite = dxB.cmap.begin(); ite != dxB.cmap.end(); ++ite)
      {
	int i = ite->first;
	double dbl = dxB.getComponent(i);
	oldVec2.setComponent(i, dbl); 
	if(dbl != 0.0)
	  oldVec2.cmap[i] = dbl;
      }

      dxB.oldx.push_back(oldVec2); 
      dxB.oldIndex.push_back(i0);

      uvec.oldx.push_back(oldVec2);
      uvec.oldIndex.push_back(i0);

      //cout << xBs * *cB << endl;
      time5 = clock() - time5;
      time1 = clock() - time1;
      cout << (double) time1 / CLOCKS_PER_SEC << " " << 
	      (double) time2 / CLOCKS_PER_SEC << " " << 
	      (double) time4 / CLOCKS_PER_SEC << " " << 
	      (double) time5 / CLOCKS_PER_SEC << " " << 
	      (double) time6 / CLOCKS_PER_SEC << endl;
      
      ++count;
    }
    cout << "opt val = " << xBs * *cB + dblinit << endl;
    time = clock()-time;
    cout << "time " << (double) time / CLOCKS_PER_SEC << endl;
    getchar();

    ///////////////////////////////////////////

    //2.xBを計算する
    xB->LUdecomposition(B, *b);

    while(1)
    {
      double val = (*cB) * (*xB);
      cout << val << endl;
      bvar->LUdecomposition(B, *b);

      //転置行列を取得
      AbstMatrix *Bt = B->getTraceMatrix(), *Nt = N->getTraceMatrix(); 
      //3.シンプレックス乗数を計算する
      pi->LUdecomposition(Bt, *cB);
      //pi->SimultaneousLinearEquations(Bt, *cB , enmGaussJordan);

      int nonbasicidx = getNonBasicIdx(cN, Nt, pi);
      if(nonbasicidx < 0) break; 

      for(int i = 0; i < Nconstraint; ++i){ ai->setComponent(i, N->getComponent(i, nonbasicidx)); }

      //4.ベクトル y = B-1*aを計算する
      y->LUdecomposition(B,*ai);

      //5.非基底変数と基底変数を交換する
      swapVariableRevised(nonbasicidx, y, xB, xN, bvar, basicVariableIdx, nonbasicVariableIdx, cN, cB,B, N);

      ++count;
    }

    time = clock()-time;
    cout << (double) time / CLOCKS_PER_SEC << endl;
    cout << "done " << endl;
    getchar();

    ai->clear();
    bvar->clear();
    xB->clear();
    xN->clear();
    pi->clear();
    b->clear();
    A->clear();
    B->clear();
    N->clear();
    y->clear();

    delete [] basicVariableIdx;
    delete [] nonbasicVariableIdx;
    delete bvar;
    delete xB;
    delete xN;
    delete pi;
    delete b;
    delete A;
    delete B;
    delete N;
    delete y;
  }
#endif
