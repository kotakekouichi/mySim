#ifndef _FUNC_H
#define _FUNC_H
#include "utils.h"
#include "params.h"

//関数の種類
enum efunc
{
  object,
  constraint,
  equation,
};

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
    int index;
    int ix;
    coeff *nextCoeff; 
    coeff *preCoeff;
    double val;
    coeff()
    {
      this->val = 0.0;
    }
};

//関数クラス
class AbstFunction
{
  private:
      efunc functype;
      eObject objecttype;
      int ID;
  public :

    int ix;
    coeff *firstCoeff;
    coeff *Coeff;
    double constTermVal;

    map<int, double> CoeffMap;

    AbstFunction()
    {
    }
    AbstFunction(efunc type, int id)
    {
      this->functype = type;
      this->ID = id;
    }

    inline void setfunctype(efunc functype, int id){this->functype = functype; this->ID = id;}
    //inline void setfunctype(efunc functype){this->functype = functype;}
    inline efunc getfunctype(){return this->functype;}
    inline void setID(int val){this->ID = val;}
    inline int getID(){return this->ID;}

    inline void setobject(eObject value){this->objecttype = value;}
    inline eObject getobject(){return this->objecttype;}

    void conectChain()
    {
      this->firstCoeff = NULL;
      coeff *Coeff = this->firstCoeff;

      for(int i = 0; i < Nx+1; ++i)
      {
	if(this->Coeff[i].val == 0) continue;

	Coeff = &this->Coeff[i];
	Coeff->nextCoeff = this->firstCoeff;
	if(this->firstCoeff != NULL)
	  this->firstCoeff->preCoeff = Coeff;
	this->firstCoeff = Coeff;
      }

    }
    
    void setCoeff(int key, double value)
    {
      this->CoeffMap[key] = value;
    }

    double getCoeff(int key)
    {
      if(this->CoeffMap.find(key) != this->CoeffMap.end())
	return this->CoeffMap[key];
      else 
	return 0.0;
    }

    void clearCoeff(){this->CoeffMap.clear();}

    void eraseCoeff(int key)
    {
      this->CoeffMap.erase(this->CoeffMap.find(key));
    }

    bool findCoeff(int key)
    {
      if(this->CoeffMap.find(key) != this->CoeffMap.end())
      {
        return true;
      }
      else 
      {
	return false;
      }
    }
    
    void Add(std::string &strConstFunc)
    {

      std::string item = "", strx = "";
      double dblval = 0.0;
      int i = 0, idx = 0;

      while(i < strConstFunc.size()) 
      {
	char ch = static_cast<char>(strConstFunc[i]);
	++i;
	if(ch == ' ' || ch == '\t') 
	{
	  cout << ch << endl;
	  continue; // 空白などは読みとばす
	}
	else if(ch == 'x')
	{
	  while(1)
	  {
	    ch = static_cast<char>(strConstFunc[i]); // xのインデックスを取得
	    strx.push_back(ch);

	    if(i + 1 < strConstFunc.size() ) 
	    {
	      char tmp = static_cast<char>(strConstFunc[i + 1]);
	      if(tmp != '\t' && tmp != ' '  && tmp != '+' && tmp != '-') 
	      {
		++i;
	      }
	      else break;
	    }
	    else break;
	  }
	  	  
	  stringtovar(strx, idx); // インデックスを数値に変換
	  strx.clear();
	  
	  if(item == "" || item == "+") item = "1";
	  else if(item == "-") item  = "-1";
	  //cout << "item = " << item << endl;
	  stringtovar(item, dblval); // 係数を数値に変換

	  if(idx > Nx) idx = Nx;

	  //this->Coeff[idx].val = dblval; 
	  this->setCoeff(idx, dblval);

	  item.clear(); 
	  ++i;
	  continue;
	}
	else if(ch == '<')
	{
	  //for(int i = 0; i <= Nx; ++i) this->Coeff[i].val *= -1; // 左辺から移行されるのでマイナスをかける
	  //for(int i = 0; i <= Nx ; ++i) this->Coeff[]
	  for(map<int , double>::iterator ite = this->CoeffMap.begin(); ite != this->CoeffMap.end(); ++ite){ite->second *= -1;}
	  while(i < strConstFunc.size())
	  {
	    ch = static_cast<char>(strConstFunc[i]); 
	    item.push_back(ch); 
	    ++i;
	  }

	  stringtovar(item, dblval);
	  //this->Coeff[0].val = dblval;
	  this->constTermVal = dblval;
	}
	item.push_back(ch);
      }

    }
    
};
#endif
