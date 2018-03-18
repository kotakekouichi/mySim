#include "minimum_spanning_tree.h"
#include "params.h"
#include "utils.h"

using namespace std;

namespace OptimizeName
{

  void MinimumSpanningTree::initialize()
  {
    for(int i = 0; i < Num; ++i)
    {
      for(int j = i + 1; j < Num; ++j)
      {
	PathClass tmpPath = PathClass(&nodes[i], &nodes[j]);
	double dbl = tmpPath.getval();
	PathMap[dbl] = tmpPath;
      }
    }
  }

  void MinimumSpanningTree::free()
  {
    T = 0;
    distmin = 1.0e+30;
    vector<PathClass*> *vec = getAbsolutelyPath();
    vec->clear();
    vec = getPathvec();
    vec->clear();
    map<double, PathClass> *Map = getPathMap();
    Map->clear();
  }

  void MinimumSpanningTree::opt()
  {
    int ncount = 0;
    T = 0;
    NodeClass *Nodes = getNodes();

    vector<PathClass*> *abstPathvec = getAbsolutelyPath();
    vector<PathClass*> *PathVec = getPathvec();

    UnionFinding uf(Num);

    //初期化 
    for(int i = 0; i < Num; ++i)
      Nodes[i].getNextPathNode()->clear();

    for(int idxpath = 0; idxpath < abstPathvec->size(); ++idxpath)
    {
      PathClass *AbsPath = abstPathvec->at(idxpath);

      if(ncount == Num - 1){ T += AbsPath->getval(); break;}
      if(ncount >= Num) break;

      NodeClass *FromNode = AbsPath->getFromNode(), *ToNode = AbsPath->getToNode();
      if(uf.Find(FromNode->getNo(), ToNode->getNo())){ T = 1.0e+30;return;}
      uf.Union(FromNode->getNo(), ToNode->getNo());

      if(addflag)
      {
	PathVec->push_back(AbsPath);
	FromNode->AddPathNode(ToNode);
	ToNode->AddPathNode(FromNode);
      }

      ++ncount;
      T += AbsPath->getval();
    }

    map<double, PathClass> *pathmap = getPathMap();

    for(map<double, PathClass>::iterator ite = pathmap->begin(); ite != pathmap->end(); ++ite)
    {
      PathClass *tmpPath = &((ite)->second);

      if(ncount == Num - 1){ T += tmpPath->getval(); break;}
      if(ncount >= Num - 1) break;
      NodeClass *FromNode = tmpPath->getFromNode(), *ToNode = tmpPath->getToNode();


      if(!tmpPath->getApprovalFlag()) continue;

      if(uf.Find(FromNode->getNo(), ToNode->getNo())){continue;}

      uf.Union(FromNode->getNo(), ToNode->getNo());

      if(addflag)
      {
	PathVec->push_back(tmpPath);
	FromNode->AddPathNode(ToNode);
	ToNode->AddPathNode(FromNode);
      }

      T += tmpPath->getval();
      ++ncount;
    }
  }

  void MinimumSpanningTree::outputNodes()
  {
    NodeClass *Nodes = getNodes();
    ofstream nodefile("nodes.dat");

    for(int i = 0; i < Num; ++i)
      nodefile << Nodes[i].getPosx() << " " << Nodes[i].getPosy() << endl;

    nodefile.close();
  }

  void MinimumSpanningTree::outputmst()
  {
    vector<PathClass*> *PathVec = getPathvec();

    ofstream mst("mst.gnu");
    mst << "#!/bin/gnuplot" << endl;
    mst << "plot '-' w l," << endl;

    for(vector<PathClass*>::iterator ite = PathVec->begin(); ite != PathVec->end(); ++ite)
    {
      PathClass *tmpPath = *ite;
      NodeClass *FromNode = tmpPath->getFromNode();
      NodeClass *ToNode = tmpPath->getToNode();

      mst << FromNode->getPosx() << " " << FromNode->getPosy() << endl;
      mst << ToNode->getPosx() << " " << ToNode->getPosy() << endl;
      mst << endl;
      mst << endl;
    }

    mst << "e" << endl;
    mst.close();
  }

  void MinimumSpanningTree::calc2appRoute()
  {
    int iStartNode = 0;
    double dist = 0.0;
    NodeClass *Nodes = getNodes();
    NodeClass *StartNode = &Nodes[iStartNode];
    NodeClass *ToNode = NULL;
    vector<PathClass> *routevec = get2appRoute();

    StartNode->setPassingFlag(true);
    StartNode->goNextNode(NULL, StartNode, ToNode, *routevec);

    PathClass lastPath = PathClass(StartNode, &Nodes[iStartNode]);
    routevec->push_back(lastPath);

    for(int i = 0; i < routevec->size(); ++i)
      dist += routevec->at(i).getval();

    distmin = dist;
  }

  void MinimumSpanningTree::output2appRoute()
  {

    vector<PathClass> *routevec = get2appRoute();
    ofstream ofile("route.gnu");

    ofile << "#!/bin/gnuplot" << endl;
    ofile << "plot '-' w l,  'nodes.dat'" << endl;

    //for(vector<PathClass>::iterator ite = routevec->begin(); ite != routevec->end(); ++ite)
    for(int i = 0; i < routevec->size(); ++i)
    {
      //PathClass *tmpPath = *ite;
      PathClass *tmpPath = &(routevec->at(i));

      ofile << tmpPath->getFromNode()->getPosx() << " " << tmpPath->getFromNode()->getPosy() << endl;
      ofile << tmpPath->getToNode()->getPosx() << " " << tmpPath->getToNode()->getPosy() << endl;

      ofile << endl;

    }

    ofile.close();

  }

  void NodeClass::goNextNode(NodeClass *preNode, NodeClass *&FromNode, NodeClass *&ToNode, vector<PathClass> &routevec)
  {
    vector<NodeClass*> *NextPathNodeVec = this->getNextPathNode();

    //for(vector<NodeClass*>::iterator iteNode = NextPathNodeVec->begin(); iteNode != NextPathNodeVec->end(); ++iteNode)
    for(int i = 0; i < NextPathNodeVec->size(); ++i)
    {
      //NodeClass *nextnode = *iteNode;
      NodeClass *nextnode = NextPathNodeVec->at(i);
      if(nextnode == preNode) continue;
      if(nextnode->getPassingFlag()) continue;

      nextnode->setPassingFlag(true);
      ToNode = nextnode;

      PathClass onepath = PathClass(FromNode, ToNode);
      routevec.push_back(onepath);
      FromNode = ToNode;

      nextnode->goNextNode(this, FromNode, ToNode, routevec);
    } 
  }

};


