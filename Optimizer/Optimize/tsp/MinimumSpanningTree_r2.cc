#include <iostream>
#include <math.h>
#include <map>
#include <vector>
#include <algorithm>
#include "utils.h"

#define N 30

using namespace std;

class PathClass;
class NodeClass;
class ForestClass;

class NodeClass
{
  private:
    int No;
    double xp;
    double yp;
    bool passingFlag;

    NodeClass *nextNode;
    ForestClass *pForest;
    ForestClass **ppForest;
    NodeClass *firstNode;

    vector<NodeClass*> nextPathNode;

  public:
    NodeClass()
    {
      this->No = -1;
      this->xp = this->yp = -1.0;
      this->pForest = NULL;
      this->nextNode = NULL;
      this->passingFlag = false;
      this->firstNode = NULL;
    }

    NodeClass(int _No, double _xp, double _yp)
    {
      this->No = _No;
      this->xp = _xp;
      this->yp = _yp;
      this->pForest = NULL;
      this->nextNode = NULL;
      this->passingFlag = false;
      this->firstNode = NULL;
    }

    inline void setNo(int value){this->No = value;}
    inline int getNo(){return this->No;}

    inline void setPosx(double value){this->xp = value;}
    inline double getPosx(){return this->xp;}

    inline void setPosy(double value){this->yp = value;}
    inline double getPosy(){return this->yp;}

    inline void setNextNode(NodeClass *value){this->nextNode  = value;}
    inline NodeClass *getNextNode(){return this->nextNode;}

    inline void setForest(ForestClass *value){this->pForest = value;}
    inline ForestClass *getForest(){return this->pForest;}

    inline void setpForest(ForestClass *value){this->ppForest = &value;}
    inline ForestClass **getpForest(){return this->ppForest;}

    inline void AddPathNode(NodeClass *value){this->nextPathNode.push_back(value);}
    inline vector<NodeClass*> getNextPathNode(){return this->nextPathNode;} 

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
  public:
    PathClass()
    {
      this->FromNode = NULL;
      this->ToNode = NULL;
      this->dblval = 0.0;
      this->ApprovalFlag = true;
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
};

class ForestClass
{
  private:
    NodeClass *firstNode;
    NodeClass *lastNode;
    int No;

  public:
    ForestClass()
    {
      firstNode = NULL;
      lastNode = NULL;
    }

    inline void setFirstNode (NodeClass *value){this->firstNode = value;}
    inline NodeClass *getFirstNode(){return this->firstNode;}

    inline void setLastNode (NodeClass *value){this->lastNode = value;}
    inline NodeClass *getLastNode(){return this->lastNode;}

    inline void setNo(int value){this->No = value;}
    inline int getNo(){return this->No;}
};

void MininumSpanningTree_r2(map<double, PathClass> &PathMap, vector<PathClass*>&AbsolutelyPath, double &distmin);
void MakeForest(PathClass *path, vector<ForestClass*> &forestvec, vector<PathClass*> &pathvec, int &nF);


int main()
{

  NodeClass *Nodes = new NodeClass[N];
  int No = 0;
  double distmin = 1.0e+10, rnd  = 0.0;


  for(int i = 0; i < N; ++i)
  {
    No = i;
    Nodes[i].setNo(No);

    rnd = (double) rand() / (double) RAND_MAX;
    Nodes[i].setPosx(rnd);

    rnd = (double) rand() / (double) RAND_MAX;
    Nodes[i].setPosy(rnd);
  }

  map<double, PathClass> PathMap;
  vector<PathClass*> AbsolutelyPath;

  for(int i = 0; i < N; ++i)
  {
    for(int j = i + 1; j < N; ++j)
    {
      PathClass tmpPath(&Nodes[i], &Nodes[j]);
      double dbl = tmpPath.getval(); // 距離を取得
      PathMap[dbl] = tmpPath;
    }
  }

  map<double, PathClass>::iterator tmpite = PathMap.end();
  AbsolutelyPath.push_back(&tmpite->second);
  --tmpite;
  AbsolutelyPath.push_back(&tmpite->second);

  tmpite = PathMap.begin();
  tmpite->second.setApprovalFlag(false);
  
  MininumSpanningTree_r2(PathMap, AbsolutelyPath, distmin);

  return 0;

}


  

 
