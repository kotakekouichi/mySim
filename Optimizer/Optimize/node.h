#ifndef _NODE_H_
#define _NODE_H_

#include <math.h>
#include <vector>
#include <map>
#include "utils.h"

namespace OptimizeName
{
  class NodeClass;
  class PathClass;

  class NodeClass
  {
    private:
      int No;
      double xp;
      double yp;
      bool passingFlag;

      NodeClass *nextNode;
      NodeClass *firstNode;

      vector<NodeClass*> nextPathNode;

    public:
      NodeClass()
      {
	this->No = -1;
	this->xp = this->yp = -1.0;
	this->nextNode = NULL;
	this->passingFlag = false;
	this->firstNode = NULL;
	this->nextPathNode.reserve(100);
      }

      NodeClass(int _No, double _xp, double _yp)
      {
	this->No = _No;
	this->xp = _xp;
	this->yp = _yp;
	this->nextNode = NULL;
	this->passingFlag = false;
	this->firstNode = NULL;
	this->nextPathNode.reserve(100);
      }

      inline void setNo(int value){this->No = value;}
      inline int getNo(){return this->No;}

      inline void setPosx(double value){this->xp = value;}
      inline double getPosx(){return this->xp;}

      inline void setPosy(double value){this->yp = value;}
      inline double getPosy(){return this->yp;}

      inline void setNextNode(NodeClass *value){this->nextNode  = value;}
      inline NodeClass *getNextNode(){return this->nextNode;}

      inline void AddPathNode(NodeClass *value){this->nextPathNode.push_back(value);}
      inline vector<NodeClass*> *getNextPathNode(){return &(this->nextPathNode);} 

      inline void setPassingFlag(bool value){this->passingFlag = value;}
      inline bool getPassingFlag(){ return this->passingFlag;}

      inline void setFirstNode(NodeClass *value){this->firstNode = value;}
      inline NodeClass *getFirstNode(){return this->firstNode;}

      void goNextNode(NodeClass *preNode, NodeClass *&FromNode, NodeClass *&ToNode, vector<PathClass> &route);
  };

  class PathClass
  {
    private:
      NodeClass *FromNode;
      NodeClass *ToNode;
      double dblval;
      bool ApprovalFlag;

      PathClass *nextPath;
      PathClass *firstPath;
      PathClass *lastPath;

    public:

      PathClass()
      {
	this->FromNode = NULL;
	this->ToNode = NULL;
	this->dblval = 0.0;
	this->ApprovalFlag = true;
	this->nextPath = NULL;

	this->firstPath = NULL;
	this->lastPath = NULL;

      }
      PathClass(NodeClass *_FromNode, NodeClass *_ToNode)
      {
	double dx = _FromNode->getPosx() - _ToNode->getPosx();
	double dy = _FromNode->getPosy() - _ToNode->getPosy();
	this->FromNode = _FromNode;
	this->ToNode = _ToNode;
	this->dblval = sq(dx) + sq(dy);
	this->dblval = sqrt(this->dblval);
	this->ApprovalFlag = true;

	this->nextPath = NULL;
	this->firstPath = NULL;
	this->lastPath = NULL;
      }

      inline void setFromNode(NodeClass *value){this->FromNode = value;}
      inline NodeClass *getFromNode(){return this->FromNode;}

      inline void setToNode(NodeClass *value){this->ToNode = value;}
      inline NodeClass *getToNode(){return this->ToNode;}

      inline void setval(double value){this->dblval = value;}
      inline double getval(){return this->dblval;}

      inline void setApprovalFlag(bool value){this->ApprovalFlag = value;}
      inline bool getApprovalFlag(){return this->ApprovalFlag;}

      bool operator < (const PathClass& right) const 
      {
	return this->dblval < right.dblval;
      }

      inline void setFirstPath(PathClass *value){this->firstPath = value;}
      inline PathClass *getFirstPath(){ return this->firstPath;}

      inline void setNextPath(PathClass *value){this->nextPath = value;}
      inline PathClass *getNextPath(){ return this->nextPath;}

      inline void setLastPath(PathClass *value){this->lastPath = value;}
      inline PathClass *getLastPath(){ return this->lastPath;}

  };
}
#endif 


