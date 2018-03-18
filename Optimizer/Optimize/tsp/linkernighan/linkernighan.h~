#ifndef _LINKERNIGHAN_H_
#define _LINKERNIGHAN_H_
#include "opt.h"
#include "node.h" 
#include "utils.h"
#include "params.h"

namespace OptimizeName
{

  class routeInfoClass;
  
  class routeInfoClass : public PathClass
  {
    public:
      routeInfoClass()
      {
	fowardFlag = false;
	backFlag = false; 
      } 
      bool fowardFlag;
      bool backFlag;
  };

  class Linkernighan : public OptClass
  {
    public:
      double distmin;
      NodeClass *nodes; 
      NodeClass *nodes2;

      NodeClass *forwardList;
      NodeClass *backList;
      vector<PathClass> bestroute;
      vector<PathClass> bestroute2;
      map<int, map<int, PathClass> > pathmap;
      map<int, map<int, int> > condidatemap;

      vector<routeInfoClass> routes_opt;

      Linkernighan(){distmin = 1.0e+10;}
      
      void initialize();
      void opt();
      void generationImplementation();
      void readyLinkerniham(
      int improvePath(int lambda, double distRoute, double &dist,  double &g, NodeClass *&forwardFirst, NodeClass *&forwardLast, NodeClass *&backFirst, NodeClass *&backLast);
      
      int improvePath_r2(int lambda, double &dist, double distRoute, double &g, vector<int> &optid, NodeClass *firstnode, NodeClass *lastnode);

      //int improvePath(int lambda, double &dist, double distRoute, double &g, NodeClass *route, map<int, int> &indexmap);
      void outputNode();
      void outputRoute();
      void setbestroute(NodeClass *route);
      void setbestroute2(NodeClass *pforward);
      void free(); 

      NodeClass* getnodes(){return nodes;}
      double getdistmin(){return distmin;}
  };

  };
#endif 
