#include <iostream>
#include "tree.h"
#include "params.h"
#include "position.h"
#include "genetic_algorithms.h"

using namespace std;

void tree::construct_tree(Type type, int num, int &count, tree *child[])
{
  int ichild = 0;
  int dx = this->dx / 2; 
  int xcenter = 0;

  if(dx % 2 && dx != 1)
    this->xc += 1;

  for(int i = 0; i < 2; ++i)
  {
    dx = this->dx * 0.5;
    xcenter = this->xc - dx * 0.5 + i * dx;
    if(dx % 2 && dx != 1) 
    {
      dx = (i == 0) ? dx + 1 : dx - 1; 
      xcenter = ( i == 0 ) ? this->xc - dx / 2: 
	this->xc + dx / 2;
    }

    child[i] = new tree(0, dx, xcenter);
  } 

  if(type == to_geno)
  {
    for(position *iposition = this->first; iposition != NULL;)
    {
      position *pnext = iposition->nextpos;
      if(iposition->pheno_num < this->xc)
	ichild = 0;
      else 
	ichild = 1;

      ++(child[ichild]->count);
      iposition->nextpos = child[ichild]->first;
      child[ichild]->first = iposition;
      iposition = pnext;
    }

    if(child[1]->first != NULL)
      if(child[1]->first->pheno_num == num)
	count += child[0]->count;
  }
  else 
  {
    position *jposition;

    for(position *iposition = this->first; iposition != NULL;)
    {
      position *pnext = iposition->nextpos;

      if(iposition->pheno_num == -1)
      {
	jposition = iposition;
	iposition = pnext;
	continue;
      }

      if(iposition->pheno_num < this->xc) ichild = 0;
      else ichild = 1; 

      ++(child[ichild]->count);
      iposition->nextpos = child[ichild]->first;
      child[ichild]->first = iposition;
      iposition = pnext;
    }

    if(num + 1 <= child[0]->dx - child[0]->count) 
      ichild = 0;
    else 
    {
      ichild = 1;
      num = num - (child[0]->dx - child[0]->count);
    }

    child[ichild]->count += 1; 

    jposition->nextpos = child[ichild]->first;
    child[ichild]->first = jposition;

    if(ichild == 1)
      count += child[0]->count;

  }

  for(int i = 0; i < 2; ++i)
  {
    if(child[i]->count > 1)
    {
      child[i]->construct_tree(type, num, count, child[i]->nextlevel);
    }
  }

}

void tree::deconstruct_tree(tree *child[])
{
  for(int i = 0; i < 2; ++i)
  {
    if(child[i] != NULL)
      child[i]->deconstruct_tree(child[i]->nextlevel);
    delete child[i];
  }
}

void tree::walkingTree(Type type, position *pos, int &countofstore)
{
  int ichild = 0;
  int num = type == to_geno ? pos->pheno_num : pos->gene_num;
  bool flag = false;
  tree *Tree = this, **child;

  while(Tree != NULL)
  {
    flag = false;
    ichild = 0;
    child = Tree->nextlevel;

    if(type == to_geno) 
      flag = num < Tree->xc; 
    else if(type == to_pheno && child[0] != NULL && child[1] != NULL) 
      flag = num + 1 <= child[0]->dx - child[0]->count; 
    
    ichild = flag ? 0 : 1; 

    ++Tree->count;
    pos->nextpos = Tree->first;
    Tree->first = pos;

    flag = Tree->count > 1 && child[ichild] == NULL; 
    if(flag)
    {
      Tree->construct_tree(type, num, countofstore, child);
      break;
    }

    flag = ichild == 1 && child[0] != NULL;
    if(flag) 
    { 
      countofstore += child[0]->count;
      if(type == to_pheno) num = num - (child[0]->dx - child[0]->count);  
    }

    Tree = child[ichild];
  }

}

