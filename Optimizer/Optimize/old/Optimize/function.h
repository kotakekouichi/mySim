#ifndef _FUNC_H
#define _FUNC_H
#include <iostream>
#include "params.h"

//目的関数を最大化、最小化するか
enum eObject
{
  maximization,
  minimization,
};

//関数の係数クラス
class coeff 
{
  public:
    unsigned int index;
    unsigned int ix;
    //coeff *nextCoeff; 
    //coeff *preCoeff;
    //double val;
    coeff()
    {
    }
};

//関数クラス
class AbstFunction
{
  private:
      eObject objecttype;
  public :

    int ix;
    coeff *firstCoeff;
    coeff *Coeff;
    double constTermVal;

    AbstFunction()
    {
    }

    inline void setobject(eObject value){this->objecttype = value;}
    inline eObject getobject(){return this->objecttype;}
/*
    void conectChain()
    {
      this->firstCoeff = NULL;
      coeff *aCoeff = this->firstCoeff;

      for(unsigned int i = 0; i < Nx + 1; ++i)
      {
	if(this->Coeff[i].val == 0) continue;

	aCoeff = &this->Coeff[i];
	aCoeff->nextCoeff = this->firstCoeff;
	if(this->firstCoeff != NULL)
	  this->firstCoeff->preCoeff = aCoeff;
	this->firstCoeff = aCoeff;
      }
    }
    
    */
};
#endif
