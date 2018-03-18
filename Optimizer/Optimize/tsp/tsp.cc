#include <iostream>
#include <math.h>
#include <map>
#include <vector>
#include <algorithm>
#include "utils.h"

#define N 26 

using namespace std;

class PathClass;
class NodeClass;
class ForestClass;

int nn;

class NodeClass
{
  private:
    int No;
    double xp;
    double yp;
    bool passingFlag;

    NodeClass *nextNode;
    ForestClass *pForest;
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

void MinimumSpanningTree(NodeClass *Nodes, double &distmin);

void BranchAndBound(NodeClass *Nodes, double &distmin);
void recursive(NodeClass *Nodes, int nodenum, int *inode, double &distmin, double &dist, NodeClass *FirstNode);

bool cross_check(double r1_x1, double r1_y1, double r1_x2, double r1_y2, double r2_x1, double r2_y1, double r2_x2, double r2_y2);

double getdistance2app(vector<PathClass> &pathvec);

int main()
{

  NodeClass *Nodes = new NodeClass[N];
  int No = 0;
  double rnd = 0.0;
  double distmin = 1.0e+30;

  //init_rnd(); 
  for(int i = 0; i < N; ++i)
  {
    No = i;
    Nodes[i].setNo(No);

    rnd = (double) rand() / (double) RAND_MAX;
    Nodes[i].setPosx(rnd);

    rnd = (double) rand() / (double) RAND_MAX;
    Nodes[i].setPosy(rnd);
    cout << Nodes[i].getPosx() << " " << Nodes[i].getPosy() << endl;
  } 

  nn = 0;
  MinimumSpanningTree(Nodes, distmin);
  BranchAndBound(Nodes, distmin);
  cout << nn << endl;

  delete [] Nodes;

  return 0;
}

void BranchAndBound(NodeClass *Nodes, double &distmin)
{
  int nodenum = 0;
  int *inode = new int[N];
  double dist = 0.0;

  for(int i = 0; i < N; ++i)
  {
    Nodes[i].setPassingFlag(false);
    Nodes[i].setNextNode(i != N - 1 ? &Nodes[i+1] : NULL);
  }

  NodeClass *FirstNode = &Nodes[0];
  NodeClass *PreNode = NULL;

  //for(int i = 0; i < N; ++i)
  for(NodeClass *Node = FirstNode; Node != NULL; Node = Node->getNextNode())
  {
    Node->setPassingFlag(true);

    if(PreNode != NULL)
    {
      PreNode->setNextNode(Node->getNextNode());
    }

    nodenum = 0;
    //inode[nodenum] = i;
    inode[nodenum] = Node->getNo();
    
    nodenum = 1;
    dist = 0.0;
    NodeClass *NextFirstNode = Node == FirstNode ? Node->getNextNode() : FirstNode;
    recursive(Nodes, nodenum, inode, distmin, dist, NextFirstNode);

    for(int j = 0; j < N; ++j)
    {
      Nodes[j].setPassingFlag(false);
    }

    if(PreNode != NULL)
    {
      PreNode->setNextNode(Node);
    }

    PreNode = Node;
    break;
  }  

  cout << "branch and bound / " << distmin << endl;

  delete [] inode;

}
void recursive(NodeClass *Nodes, int nodenum, int *inode, double &distmin, double &dist, NodeClass *FirstNode)
{

  double tmpdist;
  NodeClass *PreNode = NULL;

  //for(int i = 0; i < N; ++i)
  for(NodeClass *Node = FirstNode; Node != NULL; Node = Node->getNextNode())
  {
    //NodeClass *Node = &Nodes[i];
    //if(Node->getPassingFlag()){ cout << nodenum << "?"  << Node->getNo()  << endl; getchar();continue;}
    Node->setPassingFlag(true);

    inode[nodenum] = Node->getNo();
    //inode[nodenum] = i;
    ++nodenum; 

    if(nodenum != N)
    {
      int from = inode[nodenum - 2];
      int to = inode[nodenum - 1];

      tmpdist = sqrt(sq(Nodes[from].getPosx() - Nodes[to].getPosx()) + sq(Nodes[from].getPosy() - Nodes[to].getPosy())); 
      dist += tmpdist;

      bool crossFlag = false;

      if(PreNode != NULL)
      {
	PreNode->setNextNode(Node->getNextNode());
      }

      double distsum = 0.0;
#if 0
      vector<NodeClass*> vec;
      for(NodeClass *tmpNode = Node == FirstNode ? Node->getNextNode() : FirstNode; tmpNode != NULL; tmpNode = tmpNode->getNextNode())
      {
	vec.push_back(tmpNode);
      }

      map<double, double> distMap;
      for(int ii = 0; ii < vec.size(); ++ii)
      {

	for(int jj = ii + 1; jj < vec.size(); ++jj)
	{
	  double dbldist = sqrt(sq(vec[ii]->getPosx() - vec[jj]->getPosx()) + sq(vec[ii]->getPosy() - vec[jj]->getPosy()));
	  distMap[dbldist] = dbldist;
	}
      }
      
      int ncount = 0;
      
      for(map<double, double>::iterator ite = distMap.begin(); ite != distMap.end(); ++ite)
      {
	if(ncount + nodenum == N)
	  break;	  
	distsum += ite->first;
	++ncount;
      }
#endif
#ifdef CROSS
/*
      if( dist + distsum < distmin) 
      {
	for(int idx = 0; idx < nodenum - 1; ++idx)
	{
	  int from2 = inode[idx];
	  int to2 = inode[idx + 1];
	  if( cross_check(Nodes[from].getPosx(), Nodes[from].getPosy(), Nodes[to].getPosx(), Nodes[to].getPosy(), 
		Nodes[from2].getPosx(), Nodes[from2].getPosy(), Nodes[to2].getPosx(), Nodes[to2].getPosy()) )
	  { crossFlag = true; break;}
	}	
      }
      */
#endif
/*
      if(dist > distmin || crossFlag || dist + distsum > distmin)
      {
	dist -= tmpdist;
	--nodenum;

	if(PreNode != NULL)
	{
	  PreNode->setNextNode(Node);
	}

	PreNode = Node;

	Node->setPassingFlag(false);

	continue;
      }
     */ 
      NodeClass *NextFirstNode = Node == FirstNode ? Node->getNextNode() : FirstNode;
      recursive(Nodes, nodenum, inode, distmin, dist, NextFirstNode);

      if(PreNode != NULL)
      {
	PreNode->setNextNode(Node);
      }

      PreNode = Node;

      dist -= tmpdist;
	
      --nodenum;
    }
    if(nodenum == N)
    {
      ++nn;
      int from = inode[nodenum - 1];
      int to = inode[0];
      double tmpdist1 = sqrt(sq(Nodes[from].getPosx() - Nodes[to].getPosx()) + sq(Nodes[from].getPosy() - Nodes[to].getPosy()));
      dist += tmpdist1;
      
      int from2 = inode[nodenum -2];
      int to2 = inode[nodenum - 1];

      
      double tmpdist2 = sqrt(sq(Nodes[from2].getPosx() - Nodes[to2].getPosx()) + sq(Nodes[from2].getPosy() - Nodes[to2].getPosy()));
      dist += tmpdist2;
      
      if(dist < distmin) 
	distmin = dist;

      dist -= tmpdist1;
      dist -= tmpdist2;

    }

    Node->setPassingFlag(false);
    PreNode = Node;
  }
}

void MinimumSpanningTree(NodeClass *Nodes, double &distmin)
{
  int No = 0;

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
  vector<PathClass> pathvec;
  vector<ForestClass*> forestvec;
  
  for(map<double, PathClass>::iterator ite = PathMap.begin(); ite != PathMap.end(); ++ite)
  {
    NodeClass *FromNode = ite->second.getFromNode();
    NodeClass *ToNode = ite->second.getToNode();

    ForestClass *Forest1 = FromNode->getForest();
    ForestClass *Forest2 = ToNode->getForest();

    if(Forest1 == NULL && Forest2 == NULL)
    {
      ForestClass *tmp = new ForestClass();
      forestvec.push_back(tmp);
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

    pathvec.push_back(ite->second);
  }

  ofstream mst("mst.gnu");
  mst<< "#!/bin/gnuplot" << endl;
  mst << "plot '-' w l " << endl;
  for(vector<PathClass>::iterator ite = pathvec.begin(); ite != pathvec.end(); ++ite)
  {
    NodeClass *FromNode = ite->getFromNode();
    NodeClass *ToNode = ite->getToNode();
    mst<< FromNode->getPosx() << " " << FromNode->getPosy() << endl;
    mst << ToNode->getPosx() << " " << ToNode->getPosy() << endl;
    mst << endl;
    mst << endl;
  }

  mst << "e" << endl;
  mst.close();

  int inode  =0;
  while(1) 
  {
    if(Nodes[inode].getNextPathNode().size() == 1 )break;
    ++inode;
  }
  
  NodeClass *StartNode = &Nodes[inode];
  NodeClass *FromNode = &Nodes[inode];
  NodeClass *ToNode = NULL;
  vector<PathClass> route;

  FromNode->setPassingFlag(true);
  StartNode->goNextNode(NULL, FromNode, ToNode, route);

  PathClass onePath = PathClass(FromNode, &Nodes[inode]);
  route.push_back(onePath);

  ofstream ofile("route.gnu");

  double dist = 0;
  ofile << "#!/bin/gnuplot" << endl;
  ofile << "plot '-' w l, 'nodes.dat'" << endl; 
  for(vector<PathClass>::iterator ite = route.begin(); ite != route.end(); ++ite)
  {
    ofile << ite->getFromNode()->getPosx() << " " << ite->getFromNode()->getPosy() << endl;
    ofile << ite->getToNode()->getPosx() << " " << ite->getToNode()->getPosy() << endl;
    ofile << endl;
    dist += ite->getval();
  }
  distmin = dist;
  cout << "2.0 approximate method / " << dist << endl;
  ofile.close();

  ofstream nodefile("nodes.dat");
  for(int i = 0; i < N; ++i)
  {
    nodefile << Nodes[i].getPosx() << " " << Nodes[i].getPosy() << endl;
  }
  nodefile.close();

  for(int i = 0; i < forestvec.size(); ++i)
    delete forestvec[i];
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

  return;
}

bool cross_check(double r1_x1, double r1_y1, double r1_x2, double r1_y2, double r2_x1, double r2_y1, double r2_x2, double r2_y2)
{
  double s = 0.0, t = 0.0;
  double a11 = r1_x1 - r2_x1, a12 = r1_x2 - r2_x1;
  double a21 = r1_y1 - r2_y1, a22 = r1_y2 - r2_y1;
  double rx = r2_x2 - r2_x1, ry = r2_y2 - r2_y1; 
  double det = a11 * a22 - a12 * a21;

  if( fabs(det) < 1.0e-16 ) return false;

  s = ( a22 * rx - a12 * ry ) / det;
  t = ( -a21 * rx + a11 * ry ) / det;

  if( s > 0 && t > 0 && s + t > 1.0 )
  { 
    return true;
  }

  return false;
}

