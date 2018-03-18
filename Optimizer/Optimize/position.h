#ifndef _POSITION_H_
#define _POSITION_H_

#include "utils.h"
class position
{
  public:
    int index;
    int gene_num;
    int pheno_num;
    double xp;
    double yp;
    
    position *nextpos;

    position()
    {
      this->index = -1;
      this->gene_num = -1;
      this->pheno_num = -1; 
      this->xp = 0.0;
      this->yp = 0.0;
      this->nextpos = NULL;
    }
    
    position(int index, double xp, double yp)
    {
      this->index = index;
      this->gene_num = -1;
      this->pheno_num = index;
      this->xp = xp;
      this->yp = yp;
      this->nextpos = NULL;
    }

/*
    double operator * (const individual &right) const 
    {
      return sqrt(sq(this->xp - right.xp) + sq(this->yp - right.yp));
    }
    */
};

#endif
