#ifndef _MST_H_
#define _MST_H_
#include "opt.h"
#include "params.h"
#include "node.h"

/* --- class and strtuc --- */

namespace OptimizeName
{
  class MinimumSpanningTree : public OptClass
  {
    private:
      double T; // 最小全木の重み
      double distmin; 
      NodeClass *nodes;
      vector<PathClass*> AbsolutelyPath;
      vector<PathClass*> pathvec; //最小全域木の枝
      vector<PathClass> route; //2近似のルート
      map<double, PathClass> PathMap;

      bool addflag;

    public :

      MinimumSpanningTree(){T = 0; distmin = 1.0e+30;addflag = true;}

      void initialize();
      void opt();
      void free();

      inline void setNodes(NodeClass *value)
      {
	nodes = value;
      }

      inline NodeClass *getNodes()
      {
	return nodes;
      }

      inline void setAbsolutelyPath(vector<PathClass*> &value)
      {
	AbsolutelyPath = value;
      }
      inline vector<PathClass*> *getAbsolutelyPath()
      {
	return &AbsolutelyPath;
      }
      inline void setPathvec(vector<PathClass*> &value)
      {
	pathvec = value;
      }
      inline vector<PathClass*> *getPathvec()
      {
	return &pathvec;
      }

      inline void setPathMap(map<double, PathClass> &value)
      {
	PathMap = value;
      }
      inline map<double ,PathClass> *getPathMap()
      { 
	return &PathMap;
      }

      void outputNodes();
      
      void outputmst();
      void output2appRoute();

      inline double getWeightTree(){ return T;}
      inline vector<PathClass*> *getTree(){return &pathvec;}
      
      void calc2appRoute();
      inline double get2appDistMin(){ return distmin;}
      inline vector<PathClass> *get2appRoute(){ return &route;}

      inline void setaddflag(bool value)
      {
	addflag = value;
      }
      
  };
};

#endif 

