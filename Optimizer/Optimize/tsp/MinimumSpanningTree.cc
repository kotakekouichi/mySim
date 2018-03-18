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
    map<int, vector<NodeClass*> > nextPathNodeMap;

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
  public:
    PathClass()
    {
      this->FromNode = NULL;
      this->ToNode = NULL;
      this->dblval = 0.0;
    }
    PathClass(NodeClass *_FromNode, NodeClass *_ToNode)
    {
      this->FromNode = _FromNode;
      this->ToNode = _ToNode;
      this->dblval = sq(_FromNode->getPosx() - _ToNode->getPosx()) + sq(_FromNode->getPosy() - _ToNode->getPosy());
      this->dblval = sqrt(this->dblval);
    }

    inline void setFromNode(NodeClass *value){this->FromNode = value;}
    inline NodeClass *getFromNode(){return this->FromNode;}

    inline void setToNode(NodeClass *value){this->ToNode = value;}
    inline NodeClass *getToNode(){return this->ToNode;}

    inline void setval(double value){this->dblval = value;}
    inline double getval(){return this->dblval;}

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

void MininumSpanningTree();
double getdistance2app(vector<PathClass> &pathvec);

int main()
{

  MininumSpanningTree();

  return 0;
}

void MininumSpanningTree()
{
  NodeClass *Nodes = new NodeClass[N];
  int No = 0;
  double rnd = 0.0;

  //init_rnd(); 
  
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

  for(int i = 0; i < N; ++i)
  {
    for(int j = i + 1; j < N; ++j)
    {
      PathClass tmpPath(&Nodes[i], &Nodes[j]);
      double dbl = tmpPath.getval();
      PathMap[dbl] = tmpPath;
    }
  }

  int nF = 0;
  vector<PathClass*> pathvec;
  
  for(map<double, PathClass>::iterator ite = PathMap.begin(); ite != PathMap.end(); ++ite)
  {
    NodeClass *FromNode = ite->second.getFromNode();
    NodeClass *ToNode = ite->second.getToNode();

    ForestClass *Forest1 = FromNode->getForest();
    ForestClass *Forest2 = ToNode->getForest();

    if(Forest1 == NULL && Forest2 == NULL)
    {
      ForestClass *tmp = new ForestClass();
      tmp->setNo(nF); ++nF;

      tmp->setFirstNode(FromNode);
      FromNode->setNextNode(ToNode);
      tmp->setLastNode(ToNode);

      FromNode->setForest(tmp);
      ToNode->setForest(tmp);
    }
    else if(Forest1 != NULL && Forest2 == NULL)
    {
      ToNode->setForest(Forest1); 
      Forest1->getLastNode()->setNextNode(ToNode);
      Forest1->setLastNode(ToNode);
    }
    else if(Forest1 == NULL && Forest2 != NULL)
    {  
      FromNode->setForest(Forest2);
      Forest2->getLastNode()->setNextNode(FromNode);
      Forest2->setLastNode(FromNode);
    }
    else if(FromNode->getForest() == ToNode->getForest())
    { 
      continue;
    }
    else 
    {

      for(NodeClass *node = Forest2->getFirstNode(); node != NULL; node = node->getNextNode())
      {
	node->setForest(Forest1);
      }

      //Forest2->getFirstNode()->setForest(Forest1);
      Forest1->getLastNode()->setNextNode(Forest2->getFirstNode());
      Forest1->setLastNode(Forest2->getLastNode());
    }

    FromNode->AddPathNode(ToNode);
    ToNode->AddPathNode(FromNode);

    
    pathvec.push_back(&ite->second);
  }

  ofstream mst("mst.gnu");
  mst<< "#!/bin/gnuplot" << endl;
  mst << "plot '-' w l " << endl;
  for(vector<PathClass*>::iterator ite = pathvec.begin(); ite != pathvec.end(); ++ite)
  {
    NodeClass *FromNode = (*ite)->getFromNode();
    NodeClass *ToNode = (*ite)->getToNode();
    mst<< FromNode->getPosx() << " " << FromNode->getPosy() << endl;
    mst << ToNode->getPosx() << " " << ToNode->getPosy() << endl;
    mst << endl;
    mst << endl;
  }

  mst << "e" << endl;
  mst.close();

  int inode  =0;
  /*while(1) 
  {
    if(Nodes[inode].getNextPathNode().size() == 1 )break;
    ++inode;
  }
  */
  NodeClass *StartNode = &Nodes[inode];
  NodeClass *FromNode = &Nodes[inode];
  NodeClass *ToNode = NULL;
  vector<PathClass> route;

  FromNode->setPassingFlag(true);
  StartNode->goNextNode(NULL, FromNode, ToNode, route);

  PathClass onePath = PathClass(FromNode, &Nodes[inode]);
  route.push_back(onePath);

  ofstream ofile("route.gnu");

  ofile << "#!/bin/gnuplot" << endl;
  ofile << "plot '-' w l, 'nodes.dat'" << endl; 
  for(vector<PathClass>::iterator ite = route.begin(); ite != route.end(); ++ite)
  {
    ofile << ite->getFromNode()->getPosx() << " " << ite->getFromNode()->getPosy() << endl;
    ofile << ite->getToNode()->getPosx() << " " << ite->getToNode()->getPosy() << endl;
    ofile << endl;
  }
  ofile.close();

  ofstream nodefile("nodes.dat");
  for(int i = 0; i < N; ++i)
  {
    nodefile << Nodes[i].getPosx() << " " << Nodes[i].getPosy() << endl;
  }
  nodefile.close();

  delete [] Nodes;
}

double getdistance2app(vector<PathClass> pathvec)
{
  double dist = 0.0;


  return dist;

}

void NodeClass::goNextNode(NodeClass *preNode, NodeClass *&FromNode, NodeClass *&ToNode, vector<PathClass> &route)
{

  vector<NodeClass*> NextPathNodeVec = this->getNextPathNode();
#if 1

  for(vector<NodeClass*>::iterator iteNode = NextPathNodeVec.begin(); iteNode != NextPathNodeVec.end(); ++iteNode)
  {
    NodeClass *nextnode = *iteNode;
    if(nextnode == preNode) continue;
    if(nextnode->getPassingFlag()) continue;

    nextnode->setPassingFlag(true);
    ToNode = nextnode; 

    //cout << "from->" << FromNode->getNo() << " " << "to->" << ToNode->getNo() << endl;

    PathClass onePath = PathClass(FromNode, ToNode);
    route.push_back(onePath); 
    FromNode = ToNode;

    nextnode->goNextNode(this, FromNode, ToNode, route);
  }
#endif
#if 0
  while(1)
  {

    NodeClass *ToNode;
    double distmin = 1.0e30;
    for(vector<NodeClass*>::iterator iteNode = NextPathNodeVec.begin(); iteNode != NextPathNodeVec.end(); ++iteNode)
    {
      if( (*iteNode) == preNode) continue;
      if( (*iteNode)->getPassingFlag()) continue;
      //vector<PathClass>::iterator itend = route.end();
      NodeClass *tmpNode = NULL;
      if(route.size() == 0)
	tmpNode = FromNode;
      else 
	tmpNode = route[route.size() - 1].getFromNode();
      PathClass tmp = PathClass(tmpNode, (*iteNode));
      if(tmp.getval() < distmin)
      {
	distmin = tmp.getval();
	ToNode = *iteNode;
      }
    }

    if(distmin == 1.0e30) return;

    ToNode->setPassingFlag(true);
    PathClass onePath = PathClass(FromNode, ToNode);
    route.push_back(onePath);
    FromNode = ToNode;

    ToNode->goNextNode(this, FromNode, ToNode, route);
  }
#endif

  return;
}
