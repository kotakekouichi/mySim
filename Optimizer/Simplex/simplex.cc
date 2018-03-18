#include <iostream>
#include <fstream>
#include "simplex.h"
#include "params.h"

using namespace std;

namespace Simplex
{
  /* --- Function --- */
  void run()
  {
    Function fobj;
    Function *B = new Function[Nconst];

    initialize(fobj, B);

    simplex(fobj, B);

    free(fobj, B);

  }

  void initialize(Function &fobj, Function *B)
  {

    coeff *Coeff;

    // set value
    fobj.firstCoeff = NULL;
    fobj.Coeff[0].val = 0;
    fobj.Coeff[0].ix = 0;
    fobj.Coeff[0].index = 0;

    fobj.Coeff[1].val = 1;
    fobj.Coeff[1].ix = 1;
    fobj.Coeff[1].index = 1;

    fobj.Coeff[2].val = 1;
    fobj.Coeff[2].ix = 2;
    fobj.Coeff[2].index = 2;

    fobj.Coeff[3].val = 1;
    fobj.Coeff[3].ix = 3;
    fobj.Coeff[3].index = 3;

    B[0].firstCoeff = NULL;
    B[0].Coeff[0].val = 1;
    B[0].Coeff[0].ix = 0;
    B[0].Coeff[0].index = 0;

    B[0].Coeff[1].val = 3;
    B[0].Coeff[1].ix = 1;
    B[0].Coeff[1].index = 1;

    B[0].Coeff[2].val = -2;
    B[0].Coeff[2].ix = 2;
    B[0].Coeff[2].index = 2;

    B[0].Coeff[3].val = 2;
    B[0].Coeff[3].ix = 3;
    B[0].Coeff[3].index = 3;

    B[1].firstCoeff = NULL;
    B[1].Coeff[0].val = 3;
    B[1].Coeff[0].ix = 0;
    B[1].Coeff[0].index = 0;

    B[1].Coeff[1].val = 0;
    B[1].Coeff[1].ix = 1;
    B[1].Coeff[1].index = 1;

    B[1].Coeff[2].val = -1;
    B[1].Coeff[2].ix = 2;
    B[1].Coeff[2].index = 2;

    B[1].Coeff[3].val = -1;
    B[1].Coeff[3].ix = 3;
    B[1].Coeff[3].index = 3;

    B[2].firstCoeff = NULL;
    B[2].Coeff[0].val = 2;
    B[2].Coeff[0].ix = 0;
    B[2].Coeff[0].index = 0;

    B[2].Coeff[1].val = -2;
    B[2].Coeff[1].ix = 1;
    B[2].Coeff[1].index = 1;

    B[2].Coeff[2].val = 0;
    B[2].Coeff[2].ix = 2;
    B[2].Coeff[2].index = 2;

    B[2].Coeff[3].val = -1;
    B[2].Coeff[3].ix = 3;
    B[2].Coeff[3].index = 3;

    for(int i = 0; i < Nconst; ++i)
    {
      //B[i].firstCoeff = NULL;
      B[i].ix = i + Nx;

      for(int j = 0;j < Nx; ++j)
      {
	//B[i].Coeff[j].val = (j == 0 ? 1 : -1) * 2 * (j + i + 1);
	//B[i].Coeff[j].ix = j;
	//B[i].Coeff[j].index = j;
      }
    }

    // list of objective Function
    Coeff = fobj.firstCoeff;

    for(int i = 0; i < Nx; ++i)
    {
      if(fobj.Coeff[i].val == 0) continue;

      Coeff = &fobj.Coeff[i];
      Coeff->nextCoeff = fobj.firstCoeff;
      if(fobj.firstCoeff != NULL)
	fobj.firstCoeff->preCoeff = Coeff;
      fobj.firstCoeff = Coeff;
    }

    // list of constraint Function 
    for(int i = 0; i < Nconst; ++i)
    {
      B[i].firstCoeff = NULL;
      Coeff = B[i].firstCoeff;

      for(int j = 0;j < Nx; ++j)
      {
	if(B[i].Coeff[j].val == 0) continue;

	Coeff = &B[i].Coeff[j];
	Coeff->nextCoeff = B[i].firstCoeff;
	if(B[i].firstCoeff != NULL) 
	  B[i].firstCoeff->preCoeff = Coeff;
	B[i].firstCoeff = Coeff;
      }
    }

  }
  
  void free(Function &fobj, Function *B)
  {
    delete [] fobj.Coeff;
    for(int i = 0; i < Nconst; ++i)
      delete [] B[i].Coeff;
    delete [] B;
  }

  void simplex(Function &fobj, Function *B)
  {
    int xpivod; 
    double *xvec = new double[Nx + Nconst];
    double y;

    cout << "before" << endl;
    cout << "fobj = ";
    fobj.CheckFunc();
    cout << endl;
    for(int ic = 0; ic < Nconst; ++ic)
    {
      cout << "B = ";
      B[ic].CheckFunc();
      cout << endl;
    }

    xvec[0] = 1;
    for(int i = 1; i < Nx; ++i)xvec[i] = 0;
    for(int i = 0; i < Nconst; ++i) xvec[i + Nx] = B[i].Coeff[0].val;

    y = 0;
    // decide pivod

    while(xpivod != 0)
    {
      xpivod = 0; 
      for(coeff *Coeff = fobj.firstCoeff; Coeff != NULL; Coeff = Coeff->nextCoeff)
      {
	if(Coeff->val > 0 && Coeff->index >0)
	{
	  xpivod = Coeff->index;
	  break;
	}
      }


      if(xpivod == 0) break;

      int xmin = 1e6;
      double dxmin = 1.0e6;

      cout << xpivod << endl;

      for(int i = 0; i < Nconst; ++i)
      {
	double dx = -B[i].Coeff[0].val / B[i].Coeff[xpivod].val;
	if(B[i].Coeff[xpivod].val == 0.0) continue;
	if(dxmin > dx && dx > 0) 
	{
	  xmin = i;
	  dxmin = dx; 
	}
	//cout << "dx = " << dx << endl;
      }

      cout << "xpivod = " << xpivod <<  "xmin = " << xmin << endl;
      getchar();
      int ix = B[xmin].ix;
      double a = B[xmin].Coeff[xpivod].val * -1;

      B[xmin].ix = B[xmin].Coeff[xpivod].ix;
      B[xmin].Coeff[xpivod].val = -1.0;
      B[xmin].Coeff[xpivod].ix = ix;

      for(coeff *Coeff = B[xmin].firstCoeff; Coeff != NULL; Coeff = Coeff->nextCoeff)
      {
	Coeff->val /= a;
      }

      // y = const + a1*x1 + ... + apivod * xpivod + .. 
      // Bxmin = x_ix = bxmin - a1,xmin * x1 - ...   
      fobj.SwapVariable(xpivod, ix, B[xmin]);

      //cout << xpivod << " " << ix << endl;

      for(int ic = 0; ic < Nconst; ++ic)
      {
	if(ic == xmin) continue;
	B[ic].SwapVariable(xpivod, ix ,B[xmin]);
      }

#if 1
      cout << "after" << endl;
      cout << "fobj = ";
      fobj.CheckFunc(); cout << endl;
      for(int ic = 0; ic < Nconst; ++ic)
      {
	cout << "B" << B[ic].ix<< " = ";
	B[ic].CheckFunc();
	cout << endl;
      }
#endif

      getchar();

    }
    delete [] xvec;
  }


  int Function::SwapVariable(int xpivod, int ix, Function Bi)
  {

    double a = this->Coeff[xpivod].val; // apivod
    this->Coeff[xpivod].ix = ix;

    if(a == 0.0) return 0;

    //cout << "start" << endl;
    for(coeff *Coeff = Bi.firstCoeff; Coeff != NULL; Coeff = Coeff->nextCoeff)
    {
      int index = Coeff->index; // index of array 
      int ix = Coeff->ix; // index of x

      //cout << index << " " << ix << "->" << this->Coeff[index].val <<  endl;
      if(this->Coeff[index].val == 0.0)
      {
	coeff *CoeffNew = &this->Coeff[index];
	CoeffNew->index = Coeff->index;
	CoeffNew->ix = Coeff->ix;
	CoeffNew->val = a * Coeff->val;

	CoeffNew->nextCoeff = this->firstCoeff;
	this->firstCoeff->preCoeff = CoeffNew;
	this->firstCoeff = CoeffNew;

      }
      else 
      {
	this->Coeff[index].val = a * Coeff->val 
	  + (index != xpivod ? this->Coeff[index].val : 0);  

	if(this->Coeff[index].val == 0.0)
	{
	  //if(this->Coeff[index].nextCoeff != NULL)
	    this->Coeff[index].preCoeff->nextCoeff = this->Coeff[index].nextCoeff;
	}
      }
    }

    //cout << "endl" << endl;
    return 0;
  }

  void Function::CheckFunc()
  {
    for(coeff *Coeff = this->firstCoeff; Coeff != NULL; Coeff = Coeff->nextCoeff)
      cout << Coeff->val << "*x" << Coeff->ix << " + ";
  }

};
