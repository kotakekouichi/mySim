#ifndef _BRANCH_AND_CUT_H_
#define _BRANCH_AND_CUT_H_
#include "opt.h"
#include "params.h"
#include "function.h"
#include "simplex.h"
#include "node.h" 

using namespace std;

namespace OptimizeName
{

  class BranchAndCut : public OptClass
  {
    private:
      int **index2Nums; 
      int *xVariableIdx;
      NodeClass *nodes;
      vector<Function> partialRouteConstraintVec;
      vector<PathClass> initroutevec;
      map<int, int> initrouteIdxmap;

    public:

      BranchAndCut()
      {
  	{this->fobj = new Function();g.reserve(100000);this->fobj->constTermVal = 0.0;}
      }

      void initialize();
      void initializeFunction();
      void opt();
      void free();
      
      void partialRouteConstraint(vector<Function> &tmpvec, OptClass *opt);

      void outputRoute(OptClass *Opt);

      inline int** getIndex2Nums(){ return this->index2Nums;}
      
      void addPartialRouteConstraint(vector<Function> &value);
      inline vector<Function> *getPartialRouteConstraint(){return &(this->partialRouteConstraintVec);}

      void setInitRouteVec(vector<PathClass> &value){this->initroutevec = value;}
      inline vector<PathClass> *getInitRouteVec(){return &(this->initroutevec);}

      void setxVariableIdx(int *value){xVariableIdx = value;}
      inline int *getxVariableIdx(){return xVariableIdx;}
      
      void setInitRouteIdxMap(map<int, int> &value){this->initrouteIdxmap = value;}
      map<int, int> getInitRouteIdxMap(){return this->initrouteIdxmap;}
      
      ~BranchAndCut(){}
  };


};

#endif 

