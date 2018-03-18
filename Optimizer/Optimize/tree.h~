#ifndef _TREE_H
#define _TREE_H

#include "params.h"
#include "genetic_algorithms.h"
#include "position.h"

class tree;

class tree
{

  public:
    int count;
    int dx;
    int xc;
    int num;

    tree *nextlevel[2];
    position *first;

    tree()
    {
      this->count = 0;
      this->dx = 0;
      this->xc = 0;
      this->num = -1;
      this->first = NULL;
      this->nextlevel[0] = NULL;
      this->nextlevel[1] = NULL;
    }

    tree(int count, int dx, int xc)
    {
      this->count = count;
      this->dx = dx;
      this->xc = xc;
      this->num = -1;
      this->first = NULL;
      this->nextlevel[0] = NULL;
      this->nextlevel[1] = NULL;

    }

    void construct_tree(Type type, int num, int &count, tree *child[]);
    void deconstruct_tree(tree *child[]);
    void construct_tree_v2(int num, int &count, tree *child[]);
    void walkingTree(Type type, position *pos, int &countofstore);

};

#endif
