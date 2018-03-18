#ifndef _GA_H_
#define _GA_H_

#include "position.h"
#include "params.h"
#include <vector>
#include <map>

using namespace std;

enum Type
  {
    to_geno,
    to_pheno,
  };

namespace OptimizeName 
{
 
  class individual
  {
    public:
      int index;
      int Ngeneration;
      int generationtype;
      int parent[2];
      int cut;
      int BirthGeneration;
      double eval;
      double totaldistance;

      position *firstpos;
      position *poslist;
      //void **ad;
      vector<void*> ad;

      individual()
      {
	this->Ngeneration = 0;
	this->index = -1;
	this->firstpos = NULL;
	this->eval = 0.0;
	this->totaldistance = 0;
	//this->poslist = new position[Nc];
	this->BirthGeneration = 0;
	//this->ad = n <ew void*[Nc];
      }

      void calc_eval(position *pos);
      inline double geteval(){return (this->eval);}
      inline double getTotaldistance(){return (this->totaldistance);}
      void to_genotype();
      void to_phenotype();
      void outputRoute(position *pos);
      void mutation();
      void localsearch(position *pos);

      bool operator<(const individual & right) const
      {
	return eval < right.eval; 
      }     

      bool operator>(const individual & right) const
      {
	return eval > right.eval; 
      } 

      ~individual()
      {
	//delete [] this->poslist;
      }
  };

  class GAClass
  {
    public:
      GAClass(){}
    void run();
    void initialize(position *pos, individual *ind);
    void init_pos(position *pos);
    void init_indivual(position *pos, individual *ind);
    void GA(position *pos, individual *ind);
    void select_cross(individual *ind, position *pos);
    void free();
    void outputRoute(individual *ind, position *pos);
    void freefunction(individual *ind, position *pos);
    void select_parent(double p, int &iparent, int &jparent, individual *individual_parents);
    int select_parent_roulette(double value, individual *parents);
    void birth_individual_child(int cutIdx, individual *iparant, individual *jparent, individual *ichild , individual *jchild, position *pos);
    inline int getcutpoint()
    {
      int ir = -1;
      while(1)
      {
	ir = rand() % Nc; 
	if(ir != 0 && ir != Nc - 1) break;
      }
      return ir; 
    }
    
    map<int, int>mk_roulette(individual *ind);

  };

  inline double length(position A, position B)
  {
    return sqrt(sq(A.xp - B.xp) + sq(A.yp - B.yp));
  }

};

#endif 
