#include <iostream>
#include "genetic_algorithms.h"
#include "position.h"
#include "params.h"
#include "utils.h"
#include "tree.h"
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

namespace OptimizeName 
{

void GAClass:: run()
  {

    position *pos = new position[Nc];
    individual *ind = new individual[Nind]; 

    initialize(pos, ind);

    GA(pos, ind);

    for(int i = 0; i < Nind; ++i) delete [] ind[i].poslist;
    
    delete [] pos;
    delete [] ind;
  }

  void GAClass::initialize(position *pos, individual *ind)
  {
    init_pos(pos);
    init_indivual(pos, ind);
  }

  void GAClass::init_pos(position *pos)
  {
    //for(int i = 0; i < Nc; ++i)
    int i = 0;
    while(i != Nc)
    {
      double xp = 0.0;
      double yp = 0.0;
      rnd1(xp);rnd1(yp);
      
      //double l = sqrt(sq(xp - 0.5) + sq(yp - 0.5));
      //if(0.45 < l && l < 0.5)
      {
	position p = position(i, xp, yp);

	pos[i] = p;
	++i;
      }
    }

  }

  void GAClass::init_indivual(position *pos, individual *ind)
  {
    int *idx = new int[Nc];

    for(int i = 0; i < Nc; ++i) idx[i] = i;

    for(int i = 0; i < Nind; ++i)
    {
      ind[i].Ngeneration = 0;
      ind[i].poslist = new position[Nc];

      shuffle(idx, Nc, true);

      for(int j = Nc - 1; j >= 0; --j)
      {
	int ipos = idx[j];
	position *postmp = &(ind[i].poslist[j]);

	postmp->index = pos[ipos].index;
	postmp->xp = pos[ipos].xp;
	postmp->yp = pos[ipos].yp;
	postmp->pheno_num = pos[ipos].index;
	postmp->gene_num = pos[ipos].gene_num;

      }

      ind[i].calc_eval(pos);
    }

    sort(ind, ind + Nind, std::greater<individual>());

    delete [] idx;
  }

  void GAClass::GA(position *pos, individual *ind)
  {

    int Ng = 0;
    
    while(1)
    {

      outputRoute(ind, pos);

      select_cross(ind, pos);

      cout << "Ngeneration "  << Ng << endl;

      freefunction(ind, pos);

      Ng += 1;
      if(Ng == 20000)
	break;
    }

  }

  void GAClass::free()
  {

  }

  void individual::calc_eval(position *pos)
  {
    position *ipos, *jpos;
    this->eval = 0.0;	
    this->totaldistance = 0.0;

    for(int i = 0; i < Nc; ++i)
    {
      ipos = &pos[this->poslist[i].pheno_num];
      jpos = i + 1 != Nc ? &pos[this->poslist[i+1].pheno_num] 
	: &pos[this->poslist[0].pheno_num];

      this->totaldistance += sqrt(sq(ipos->xp - jpos->xp) + sq(ipos->yp - jpos->yp));
    }

    this->eval = 1.0 / this->getTotaldistance();
  }

  void individual::to_genotype()
  {
    int countofstore = 0;
    int dx = (Nc % 2) ? Nc + 1: Nc;
    int xcenter = dx / 2;
    int ipos = 0;
    tree Root(0, dx, xcenter);
    position *posarray = new position[Nc];

    //for(position *pos = this->firstpos; pos != NULL;pos = pos->nextpos)
    for(int i = 0; i < Nc; ++i)
    {
      position *pos = &this->poslist[i];
      posarray[ipos].pheno_num = pos->pheno_num;
      posarray[ipos].gene_num = pos->gene_num;
      ++ipos;
    }
    
    ipos = 0;
    //for(position *pos = this->firstpos; pos != NULL;pos = pos->nextpos)
    for(int i = 0; i < Nc; ++i)
    {
      position *pos = &this->poslist[i];
      countofstore = 0;

      Root.walkingTree(to_geno, &posarray[ipos], countofstore);

      pos->gene_num = pos->pheno_num - countofstore;

      ++ipos;
    }
    
    Root.deconstruct_tree(Root.nextlevel);
    
    delete [] posarray;
  }

  void individual::to_phenotype()
  {
    int countofstore = 0;
    int dx = (Nc % 2) ? Nc + 1: Nc;
    int xcenter = dx / 2;
    int ipos = 0;
    tree Root(0, dx, xcenter);

    position *posarray = new position[Nc];

    //for(position *pos = this->firstpos; pos != NULL; pos = pos->nextpos)
    for(int i = 0; i < Nc; ++i)
    {
      position *pos = &this->poslist[i];
      posarray[ipos].pheno_num = -1;
      posarray[ipos].gene_num = pos->gene_num;
      posarray[ipos].nextpos = NULL;
      ++ipos;
    }

    ipos = 0;
    //for(position *pos = this->firstpos; pos != NULL; pos = pos->nextpos)
    for(int i = 0; i < Nc; ++i)
    {
      position *pos = &this->poslist[i];
      countofstore = 0;

      Root.walkingTree(to_pheno, &posarray[ipos], countofstore);

      posarray[ipos].pheno_num = countofstore + pos->gene_num;
      pos->pheno_num = countofstore + pos->gene_num;

      ++ipos;
    }

    Root.deconstruct_tree(Root.nextlevel);
    
    delete [] posarray;
  }

  void individual::outputRoute(position *pos)
  {
    string Nstr;
    int n = 0;
    //vartostring(Nstr, this->Ngeneration);  
    vartostring(Nstr, n);
    ofstream routefile("./route" + Nstr + ".gnu");

    routefile << "#!/bin/gnuplot" << endl;
    routefile << "set xrange [-0.5:1.5]" << endl;
    routefile << "set yrange [-0.5:1.5]" << endl;
    routefile << "set size square" << endl;
    routefile << "plot '-' w l" << endl;

    //for(position *pos = this->firstpos; pos != NULL; pos = pos->nextpos)
    for(int i = 0; i < Nc; ++i)
    {
      position *Pos = &pos[this->poslist[i].pheno_num];
      routefile << Pos->xp << "\t" << Pos->yp << endl; 
    }
    //routefile << this->firstpos->xp << "\t" << this->firstpos->yp << endl;
    routefile << pos[this->poslist[0].pheno_num].xp << "\t" << pos[this->poslist[0].pheno_num].yp << endl;


    routefile << "e" << endl;
    routefile.close();
  }

  void individual::mutation()
  {
    //int mutate_point = 0;
    //int mutate_val = 0;
#if 0    
    while(mutate_point == 0 || mutate_point == Nc - 1)
    {
      mutate_point = rand() % Nc;
    }
    
    mutate_val = rand() % (Nc - mutate_point);

    
    this->poslist[mutate_point].gene_num = mutate_val;
#endif

    this->to_phenotype();

    int ipos = 0, jpos = 0;

    while(1)
    {
      ipos = rand() % Nc;
      jpos = rand() % Nc;

      if(ipos != 0 && jpos != 0 && ipos != jpos) 
       break;	
    }

    position iPos = this->poslist[ipos];
    position jPos = this->poslist[jpos];
    position tmp = jPos;

    jPos = iPos;
    iPos = tmp;
    
    this->poslist[ipos] = iPos;
    this->poslist[jpos] = jPos;

    this->to_genotype();
  
  }

  void individual::localsearch(position *pos)
  {
    for(int iloop = 0; iloop < 500; ++iloop)
    {
      int ipos = 0, jpos = 0;
      double oldval = this->geteval();

      ipos = rand() % (Nc - 2);
      jpos = rand() % (Nc - ipos - 2) + ipos + 2; 

      //ipos までは同じ経路
      double deltal = length(this->poslist[ipos], this->poslist[ipos + 1])
	            + length(this->poslist[jpos], this->poslist[jpos + 1])
		    - length(this->poslist[ipos], this->poslist[jpos]) 
		    - length(this->poslist[ipos + 1], this->poslist[jpos + 1]);
      
      if(deltal < 0.0) continue;

      individual indnew;
      indnew.poslist = new position[Nc];
      
      int npos = 0;
      for(int i = 0; i <= ipos; ++i, ++npos)
      {
	indnew.poslist[npos] = this->poslist[i];
      } 

      for(int i = jpos; i >= ipos + 1; --i, ++npos)
      {
	indnew.poslist[npos] = this->poslist[i];
      }

      for(int i = jpos + 1; i < Nc; ++i, ++npos)
      {
	indnew.poslist[npos] = this->poslist[i];
      }

      indnew.calc_eval(pos);

      if(indnew.geteval() > oldval)
      {
	for(int i = 0; i <Nc; ++i)
	{
	  this->poslist[i] = indnew.poslist[i];
	}
	this->eval = indnew.geteval();
	this->totaldistance = indnew.getTotaldistance();
      }

      delete [] indnew.poslist;
    }
  }
  void GAClass::select_cross(individual *ind, position *pos)
  {
    //const double eliterate = 0.1; 
    //const int Nelite = (int) (Nind * eliterate) % 2 == 0 ? Nind * eliterate : Nind * eliterate + 1;
    //const int Nelite = 20;
    const int Nelite = Nind / 2;
    const int Mind = Nind;
    double p = 0.0;
    int iM  = 0, count_new = 0;
    individual *individual_new = new individual[Nind + Mind];
    individual *individual_parents = new individual[Nind];
    individual *childI, *childJ;

    for(int i = 0; i < Nind; ++i) ind[i].to_genotype();

    /// エリート戦略 //
    for(int i = 0; i < Nelite; ++i)
    {
      individual_new[i] = ind[i];
      individual_new[i].Ngeneration += 1;
      ++count_new;
    }

    //エリート以外の個体を格納（親の候補）
    //for(int i = Nelite, j = 0; i < Nind; ++i, ++j)
    for(int i = 0, j = 0; i < Nind; ++i, ++j)
    {
      individual_parents[j] = ind[i];
      p += individual_parents[j].geteval();
    }
    
    //const int Num = Nind  / 2;
    //const int Num = Nelite;
    int count = 0;
    for(int i = 0; i < Nind; ++i)
    {
      //p += ind[i].geteval();
    }

    /*
    while(1)
    {
      double probability = 0.0;
      int kparent = 0;
      rnd1(probability);
      probability *= p;
      kparent = select_parent_roulette(probability, ind);

      individual_parents[count] = ind[kparent];
      ++count;
      if(count == Num) break;
    }
    */
    
    //親から子供を生む
    int iparent = 0, jparent = 0;
    
    while(1)
    {
      select_parent(p, iparent, jparent, individual_parents);
      /*
      while(1)
      {
	iparent = rand() % Num;
	jparent = rand() % Num;

	if(iparent != jparent)
	  break;
      }
      */

      int cutIdx = getcutpoint(); 

      childI = &individual_new[Nelite + iM]; ++iM; ++count_new; childI->poslist = new position[Nc];
      childJ = &individual_new[Nelite + iM]; ++iM; ++count_new; childJ->poslist = new position[Nc];
      
      birth_individual_child(cutIdx, &individual_parents[iparent], & individual_parents[jparent], childI, childJ, pos);
      childI->parent[0] = iparent;
      childJ->parent[1] = jparent;
      childI->cut = cutIdx;
      childJ->cut = cutIdx;

      if(count_new >= Nind + Mind) break;
      //if(count_new >= Nind) break;
    }

    sort(individual_new, individual_new + Nind + Mind, greater<individual>());
    
    for(int i = Nelite; i < Nind; ++i) delete [] ind[i].poslist;
    
    ind[0] = individual_new[0];
    count = 1;
    for(int i = 1; i < Nind + Mind; ++i)
    {
      if(count < Nind)
      {
	if(ind[count - 1].geteval() != individual_new[i].geteval())
	{
	  ind[count] = individual_new[i];
	  ++count;
	}
      }
      else 
      {
	delete [] individual_new[i].poslist;
      }
      //if(count == Nind) break;
      //else delete [] individual_new[i].poslist;
    }

    for(int i = 0; i < Nind; ++i)
    {
      ind[i].localsearch(pos);
    }

    delete [] individual_new; 
    delete [] individual_parents; 

  }

  void GAClass::outputRoute(individual *ind, position *pos)
  {
    int idx = 0; 
    ind[idx].outputRoute(pos);
  }

  void GAClass::freefunction(individual *ind, position *pos)
  {

    //cout << "Number of Generation " << ind[0].Ngeneration << endl;
    //return;
    double dist = 0.0; 
    
    for(int i = 0; i < 2; ++i)
    //for(int i = 0; i < Nind; ++i)
    {
      dist = ind[i].getTotaldistance();
      cout << "(" << i << ")" << "[" << ind[i].BirthGeneration << "] " << dist << " " << ind[i].geteval() << "->";
      //for(int j = 0; j < Nc; ++j)
      {
	//position *pos = &ind[i].poslist[j];
	//cout << pos->gene_num << "(" << pos->pheno_num << ") ";
      }
      cout << endl;
      //cout << ind[i].parent[0] << " " << ind[i].parent[1] << 
	//" " << ind[i].cut  <<  " " << ind[i].generationtype << endl;
    }

  }

  map<int, int> GAClass::mk_roulette(individual *ind)
  {
    map<int, int> Map;
    
#if 0
    int P = 0;
    
    for(int i = 0; i < Nind; ++i)
    {
      P = (int) ind[i].geteval();
    }

    int idx = 0;
    for(int i= 0; i < P; ++i)
    {

    }
#endif
    return Map;
  }

  int GAClass::select_parent_roulette(double value, individual *parents)
  {
    int indIdx =  0;
    double left = 0.0;
    double right = parents[indIdx].geteval(); 

    while(1)
    {

      if(left < value && value < right) 
	break;

      ++indIdx;
      left = right;
      right += parents[indIdx].geteval();
    }
    return indIdx;
  }

  void GAClass::select_parent(double p, int &iparent, int &jparent, individual *individual_parents)
  {
    while(1)
    {

      double probability = 0.0;
      rnd1(probability); 
      probability *= p; 
      iparent = select_parent_roulette(probability, individual_parents);

      rnd1(probability);
      probability *= p;
      jparent = select_parent_roulette(probability, individual_parents); 

      if(iparent != jparent) break;
    }
  }

  void GAClass::birth_individual_child(int cutIdx, individual *iparent, individual *jparent, individual *ichild, individual *jchild, position *pos)
  {

    const double probability_mutation = 0.2;
    const double probability_cross = 0.8;
    double probability = 0.0;
    ichild->BirthGeneration = iparent->Ngeneration + 1;
    ichild->Ngeneration = iparent->Ngeneration + 1;
    jchild->BirthGeneration = jparent->Ngeneration + 1;
    jchild->Ngeneration += jparent->Ngeneration + 1;
    position *itmp, *jtmp;

    rnd1(probability);
    if(probability < 1.0 - probability_cross) cutIdx = Nc; 
    
    for(int i = 0; i < Nc; ++i)
    {
      position ipos, jpos;

      itmp = &ichild->poslist[i];
      jtmp = &jchild->poslist[i];

      if(i <= cutIdx)
      {
	ipos = iparent->poslist[i];
	jpos = jparent->poslist[i];
      }
      else
      {
	ipos = jparent->poslist[i];
	jpos = iparent->poslist[i];
      }

      itmp->index = ipos.index;
      itmp->xp = ipos.xp;
      itmp->yp = ipos.yp;
      itmp->gene_num = ipos.gene_num;
      itmp->pheno_num = ipos.pheno_num; 

      jtmp->index = jpos.index;
      jtmp->xp = jpos.xp;
      jtmp->yp = jpos.yp;
      jtmp->gene_num = jpos.gene_num;
      jtmp->pheno_num = jpos.pheno_num;

    }	

    ichild->generationtype = 0;
    jchild->generationtype = 0;
    
    rnd1(probability); 
    if(probability < probability_mutation)
    {
      ichild->mutation();
      ichild->generationtype = 1;
    }
    else ichild->to_phenotype();
    ichild->calc_eval(pos);
    
    rnd1(probability);
    if(probability < probability_mutation)
    {   
      jchild->mutation();
      jchild->generationtype = 1;
    }

    else jchild->to_phenotype();
    jchild->calc_eval(pos);

  }
};

