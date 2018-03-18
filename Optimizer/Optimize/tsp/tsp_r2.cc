#include <iostream>
#include <math.h>
#include <map>
#include <vector> 
#include <algorithm>
#include "utils.h"

//#define N 67
#define N 1000 

using namespace std;

class PathClass;
class NodeClass;

int  nn;

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

// function //
void outputNodes(NodeClass* Nodes);
double MinimumSpanningTree_r3(NodeClass *Nodes, map<double, PathClass> &PathMap, vector<PathClass*>&AbsolutelyPath, vector<PathClass*> &pathvec, bool addflag, PathClass *FisrPath);
void outputmst(vector<PathClass*> &pathvec);
double get2appRoute(NodeClass *Nodes);
void BranchAndBound(map<double, PathClass> &PathMap, NodeClass *Nodes, double &distmin);

bool cross_check(NodeClass *Nodes, vector<PathClass*> &route);

bool cross_check_twopath(double r1_x1, double r1_y1, double r1_x2, double r1_y2, double r2_x1, double r2_y1, double r2_x2, double r2_y2);
double getWeightTree(NodeClass* Nodes, map<double, PathClass> &PathMap, vector<PathClass*> &AbsolutelyPath, PathClass *FirstPath);

void localsearch(vector<PathClass> &pathvec, double &distmin);

void recursive(map<double, PathClass>::iterator ite, map<double, PathClass> &PathMap, NodeClass *Nodes, vector<PathClass*> &route, double &distmin, double &dist, int *EdgeCounter, vector<PathClass*> &AbsolutelyPath);

// main //
int main()
{ 

  NodeClass *Nodes = new NodeClass[N];
  map<double, PathClass> PathMap;
  vector<PathClass*> AbsolutelyPath;
  vector<PathClass*> pathvec;
  int No = 0;
  double rnd = 0.0;
  double distmin = 1.0e+30;

  int i =  0;
  //for(int i = 0; i < N; ++i)
  while(i < N)
  {
    No = i;

    double xp =(double) rand() / (double) RAND_MAX;
    double yp =(double) rand() / (double) RAND_MAX;

    double l = sqrt(sq(xp - 0.5) + sq(yp - 0.5));
    if(0.42 < l && l < 0.48)
    {

    Nodes[i].setNo(No);
    Nodes[i].setPosx(xp);
    Nodes[i].setPosy(yp);
    ++i;
    }
  }

  for(int i = 0; i < N; ++i)
  {
    for(int j = i + 1; j < N; ++j)
    {
      PathClass tmpPath(&Nodes[i], &Nodes[j]);
      double dbl = tmpPath.getval();
      PathMap[dbl] = tmpPath;
    }
  }

  MinimumSpanningTree_r3(Nodes, PathMap, AbsolutelyPath, pathvec, true, NULL);
  distmin = get2appRoute(Nodes);
  
  //nodeと最小全域木を出力
  outputmst(pathvec);
  outputNodes(Nodes);

  nn = 0;
  BranchAndBound(PathMap, Nodes, distmin);
  cout << nn << endl;
  
  return 0;

}

void outputNodes(NodeClass *Nodes)
{
  ofstream nodefile("nodes.dat");

  for(int i = 0; i < N; ++i)
  {
    nodefile << Nodes[i].getPosx() << " " << Nodes[i].getPosy() << endl;
  }
  nodefile.close();
}

double get2appRoute(NodeClass *Nodes)
{
  int iStartNode = 0;
  double dist = 0.0;

  NodeClass *StartNode = &Nodes[iStartNode];
  NodeClass *ToNode = NULL;
  vector<PathClass> route;

  StartNode->setPassingFlag(true);
  StartNode->goNextNode(NULL, StartNode, ToNode, route);

  PathClass  lastPath = PathClass(StartNode, &Nodes[iStartNode]);
  route.push_back(lastPath);

  ofstream ofile("route.gnu");

  ofile << "#!/bin/gnuplot" << endl;
  ofile << "plot '-' w l,  'nodes.dat'" << endl;

  for(vector<PathClass>::iterator ite = route.begin(); ite != route.end(); ++ite)
  {
    ofile << ite->getFromNode()->getPosx() << " " << ite->getFromNode()->getPosy() << endl;
    ofile << ite->getToNode()->getPosx() << " " << ite->getToNode()->getPosy() << endl;
    ofile << endl;
    dist += ite->getval();
  }

  cout << "2.0 approximate method / " << dist << endl;
 
  ofile.close();
  
  localsearch(route, dist);
  cout << "localsearch : " << dist << endl;
  
  ofstream ofile1("route.gnu");

  ofile1 << "#!/bin/gnuplot" << endl;
  ofile1 << "plot '-' w l,  'nodes.dat'" << endl;

  dist = 0;
  for(vector<PathClass>::iterator ite = route.begin(); ite != route.end(); ++ite)
  {
    ofile1 << ite->getFromNode()->getPosx() << " " << ite->getFromNode()->getPosy() << endl;
    ofile1 << ite->getToNode()->getPosx() << " " << ite->getToNode()->getPosy() << endl;
    ofile1 << endl;
    dist += ite->getval();
  }
  ofile1.close();

  cout << dist << endl;
  getchar();
  return dist;
}

void outputmst(vector<PathClass*> &pathvec)
{

  ofstream mst("mst.gnu");

  mst << "#!/bin/gnuplot" << endl;
  mst << "plot '-' w l " << endl;

  for(vector<PathClass*>::iterator ite = pathvec.begin(); ite != pathvec.end(); ++ite)
  {
    NodeClass *FromNode = (*ite)->getFromNode();
    NodeClass *ToNode = (*ite)->getToNode();
    mst << FromNode->getPosx() << " " << FromNode->getPosy() << endl;
    mst << ToNode->getPosx() << " " << ToNode->getPosy() << endl;
    mst << endl;
    mst << endl;
  }

  mst << "e" << endl;
  mst.close();
  
  return;
}

//分枝限定法
void BranchAndBound(map<double, PathClass> &PathMap, NodeClass *Nodes, double &distmin)
{

  int *EdgeCounter = new int[N];
  double dist = 0.0;
  vector<PathClass*> AbsolutelyPath, route;

  // initialize //
  for(int i = 0; i < N; ++i) 
  {
    Nodes[i].setNextNode(i != N - 1 ? &Nodes[i+1] : NULL);
    EdgeCounter[i] = 0;
  }

  recursive(PathMap.begin(), PathMap, Nodes, route, distmin, dist, EdgeCounter, AbsolutelyPath);
  
  delete [] EdgeCounter;
  
  cout << "branch and bound / " << distmin << endl;

}

void recursive(map<double, PathClass>::iterator ite, map<double, PathClass> &PathMap, NodeClass *Nodes, vector<PathClass*> &route, double &distmin, double &dist, int *EdgeCounter, vector<PathClass*> &AbsolutelyPath)
{
  if(ite == PathMap.end()){return;}
  
  double distmsp = 0.0;
  double tmpdist = ite->second.getval();
  int ifrom = ite->second.getFromNode()->getNo();
  int ito = ite->second.getToNode()->getNo();
  bool flag = true;
  PathClass *FirstPath = NULL;
  PathClass *tmpPath = &(ite->second);

  if(route.size() >= N) flag = false;

  bool tmpbool  = true;
  if(EdgeCounter[ifrom] + 1 > 2 || EdgeCounter[ito] + 1 > 2){ flag = false; tmpbool = false;}
  
  if(dist + tmpdist > distmin) flag = false;

  AbsolutelyPath.push_back(tmpPath);
  
  route.push_back(tmpPath);

  if(flag) if(cross_check(Nodes, route)) flag = false;

  if(flag)
  {
    ite->second.setApprovalFlag(true);
    FirstPath = AbsolutelyPath.size() != 0 ? AbsolutelyPath[0] : NULL;
    distmsp = getWeightTree(Nodes, PathMap, AbsolutelyPath, FirstPath);
    if(distmin < distmsp) flag = false;
  }

  ++EdgeCounter[ifrom];
  ++EdgeCounter[ito];
  dist += tmpdist;
  ++ite;


  if(route.size() == N && tmpbool)
  {
    if(dist < distmin) distmin = dist;
  }

  if(flag) recursive(ite, PathMap, Nodes, route, distmin, dist, EdgeCounter, AbsolutelyPath); 

  if(route.size() > 0) route.pop_back();
  if(AbsolutelyPath.size() > 0) AbsolutelyPath.pop_back();
  dist -= tmpdist;
  --EdgeCounter[ifrom];
  --EdgeCounter[ito];

  flag = true;
  tmpPath->setApprovalFlag(false);
  FirstPath = AbsolutelyPath.size() != 0 ? AbsolutelyPath[0] : NULL;
  distmsp = getWeightTree(Nodes, PathMap, AbsolutelyPath, FirstPath);
  if(distmsp > distmin) flag = false;
  
  if(flag) recursive(ite, PathMap, Nodes, route, distmin, dist, EdgeCounter, AbsolutelyPath);

  tmpPath->setApprovalFlag(true);

  ++nn;

  return;
}

void NodeClass::goNextNode(NodeClass *preNode, NodeClass *&FromNode, NodeClass *&ToNode, vector<PathClass> &route)
{
  vector<NodeClass*> NextPathNodeVec = this->getNextPathNode();

  for(vector<NodeClass*>::iterator iteNode = NextPathNodeVec.begin(); iteNode != NextPathNodeVec.end(); ++iteNode)
  {
    NodeClass *nextnode = *iteNode;
    if(nextnode == preNode) continue;
    if(nextnode->getPassingFlag()) continue;

    nextnode->setPassingFlag(true);
    ToNode = nextnode;

    PathClass onepath = PathClass(FromNode, ToNode);
    route.push_back(onepath);
    FromNode = ToNode;

    nextnode->goNextNode(this, FromNode, ToNode, route);
  }

  return;
}

bool cross_check(NodeClass *Nodes, vector<PathClass*> &route)
{
  int nodenum = route.size();
  int from = route[nodenum - 1]->getFromNode()->getNo();
  int to = route[nodenum - 1]->getToNode()->getNo();

  for(int idx = 0; idx < nodenum - 1; ++idx)
  {
    int from2 = route[idx]->getFromNode()->getNo();
    int to2 = route[idx]->getToNode()->getNo();

    if( cross_check_twopath(Nodes[from].getPosx(), Nodes[from].getPosy(), Nodes[to].getPosx(), Nodes[to].getPosy(),
	                    Nodes[from2].getPosx(), Nodes[from2].getPosy(), Nodes[to2].getPosx(), Nodes[to2].getPosy()))
      return true;
  }
  return false;
   
}

bool cross_check_twopath(double r1_x1, double r1_y1, double r1_x2, double r1_y2, double r2_x1, double r2_y1, double r2_x2, double r2_y2)
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

double getWeightTree(NodeClass* Nodes, map<double, PathClass> &PathMap, vector<PathClass*> &AbsolutelyPath, PathClass *FirstPath)
{
  vector<PathClass*> pathvec; pathvec.reserve(100);
  return MinimumSpanningTree_r3(Nodes, PathMap, AbsolutelyPath, pathvec, false, FirstPath);
}

double MinimumSpanningTree_r3(NodeClass *Nodes, map<double, PathClass> &PathMap, vector<PathClass*>&AbsolutelyPath, vector<PathClass*> &pathvec, bool addflag, PathClass *FirstPath)
{
  int ncount = 0;
  double T = 0.0;
  UnionFinding uf(N);
  
  //初期化
  for(int i = 0; i < N; ++i) Nodes[i].getNextPathNode().clear();

  //必ず通る経路に対して木を作成
  for(int idxpath = 0; idxpath < AbsolutelyPath.size(); ++idxpath)
  {
    PathClass *AbsPath = AbsolutelyPath[idxpath];
    
    if(ncount == N -1){ T += AbsPath->getval(); break;}
    if(ncount >= N - 1) break;
    NodeClass *FromNode = AbsPath->getFromNode(), *ToNode = AbsPath->getToNode();

    if(uf.Find(FromNode->getNo(), ToNode->getNo())){ return 1.0e+10;}
    uf.Union(FromNode->getNo(), ToNode->getNo());

    if(addflag)
    {
      pathvec.push_back(AbsPath);
      FromNode->AddPathNode(ToNode);
      ToNode->AddPathNode(FromNode);
    }
    ++ncount;
    T += AbsPath->getval();
  }

  for(map<double, PathClass>::iterator ite = PathMap.begin(); ite != PathMap.end(); ++ite)
  {
    if(ncount == N -1){ T+= ite->second.getval();break;}
    if(ncount >= N - 1) break;
    NodeClass *FromNode = ite->second.getFromNode(), *ToNode = ite->second.getToNode();
    
    if(!ite->second.getApprovalFlag()) continue;
    if(uf.Find(FromNode->getNo(), ToNode->getNo())){continue;return 1.0e+10;}

    uf.Union(FromNode->getNo(), ToNode->getNo());

    if(addflag)
    {
      pathvec.push_back(&(ite->second));
      FromNode->AddPathNode(ToNode);
      ToNode->AddPathNode(FromNode);
    }

    T += ite->second.getval();
    ++ncount;
  }

  return T;
}

void localsearch(vector<PathClass> &pathvec, double &distmin)
{
  vector<NodeClass*> nodes;
  NodeClass *fromNode, *toNode;

  for(int i = 120; i < N; ++i)
  {
    fromNode = pathvec[i].getFromNode(); 
    nodes.push_back(fromNode);
  }

  for(int i = 0; i < 120; ++i)
  {
    fromNode = pathvec[i].getFromNode(); 
    nodes.push_back(fromNode);
  }

START:
  for(int ipos = 0; ipos < N-1; ++ipos)
  {
    for(int jpos = ipos + 2;jpos < N-1; ++jpos)
    {
    //int ipos = 0, jpos = 0;
    double oldval = distmin;

    //ipos = rand() % (N - 2);
    //jpos = rand() % (N - ipos - 2) + ipos + 1;

    double del = sqrt(sq(nodes[ipos]->getPosx() - nodes[ipos + 1]->getPosx()) + sq(nodes[ipos]->getPosy() - nodes[ipos + 1]->getPosy())) + 
                 sqrt(sq(nodes[jpos]->getPosx() - nodes[jpos + 1]->getPosx()) + sq(nodes[jpos]->getPosy() - nodes[jpos + 1]->getPosy())) -
		 sqrt(sq(nodes[ipos]->getPosx() - nodes[jpos]->getPosx()) + sq(nodes[ipos]->getPosy() - nodes[jpos]->getPosy())) -
		 sqrt(sq(nodes[ipos + 1]->getPosx() - nodes[jpos + 1]->getPosx()) + sq(nodes[ipos + 1]->getPosy() - nodes[jpos + 1]->getPosy()));

    if(del < 0.0) continue;

    vector<NodeClass*> tmpvec;

    int npos = 0;
    for(int i = 0; i < N ; ++i)
      tmpvec.push_back(nodes[i]); 

    nodes.clear();

    for(int i = 0; i <= ipos; ++i, ++npos)
      nodes.push_back(tmpvec[i]);

    for(int i = jpos; i >= ipos + 1; --i, ++npos)
      nodes.push_back(tmpvec[i]);

    for(int i = jpos + 1; i < N; ++i, ++npos)
      nodes.push_back(tmpvec[i]);

    distmin -= del;
goto START;
    }
  }

  pathvec.clear();
  for(int ipos = 0; ipos < N; ++ipos)
  {
    PathClass tmppath(nodes[ipos], nodes[ipos + 1 != N ? ipos + 1 : 0]);
    pathvec.push_back(tmppath);
  }

  return;
}

