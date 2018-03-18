#include <iostream>
#include "node.h"
#include "linkernighan.h"
#include "params.h"
#include "minimum_spanning_tree.h"

using namespace std;

namespace OptimizeName
{

  void Linkernighan::initialize()
  {

    cout << "initialize " << endl;
    const int N = Num;
    double xp, yp;

    nodes = new NodeClass[N];

    int inode = 0;
    //for(int inode = 0; inode < N; ++inode)
    while(inode < N)
    {
      xp = (double) rand() / RAND_MAX; 
      yp = (double) rand() / RAND_MAX;

      if(sq(xp - 0.5) + sq(yp - 0.5) < 0.5 * 0.5  && 0.3 * 0.3 < sq(xp - 0.5) + sq(yp - 0.5))
      {
	nodes[inode].setNo(inode);
	nodes[inode].setPosx(xp);
	nodes[inode].setPosy(yp);
	nodes[inode].setPassingFlag(false);
	inode++;
      }
    }

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
	PathClass tmppath(&nodes[i], &nodes[j]);  
	pathmap[i][j] = tmppath;
      }
    }

    map<double, int> sortedmap;
    for(int i = 0; i < N; ++i)
    {
      sortedmap.clear();
      for(int j = 0; j < N; ++j)
      {
        if(i == j) continue; 

	PathClass tmppath(&nodes[i], &nodes[j]);
	double dbl = tmppath.getval();

	sortedmap[dbl] = j;

      }
 
      int c = 0;

      for(map<double, int>::iterator ite = sortedmap.begin(); ite != sortedmap.end(); ++ite, ++c)
      {

	int condidateIdx = ite->second;
	map<int, int> tmpmap;
	tmpmap[condidateIdx] = condidateIdx;
	condidatemap[i] = tmpmap;

	if(c == 10) break;

      }
    }
  }

  void Linkernighan::free()
  {
    cout << "free" << endl;
    delete [] nodes; 
    bestroute.clear();
  }

  void Linkernighan::opt()
  {
    const int N = Num;
    cout << "opt" << endl;
    generationImplementation();

  }

  void Linkernighan::generationImplementation()
  {
    int iteration = 0; 
    const int N = Num;
    double dist = 0.0;

    MinimumSpanningTree *mst = new MinimumSpanningTree();

    mst->setNodes(nodes);
    mst->initialize();
    mst->opt();
    mst->calc2appRoute();

    vector<PathClass> *route2app = mst->get2appRoute(); 

    NodeClass *route = new NodeClass[N];
    for(int i = 0; i < route2app->size(); ++i)
    {
      PathClass *tmpPath = &(route2app->at(i));
      route[i] = nodes[tmpPath->getFromNode()->getNo()];
    }

    for(int inode = 0; inode < N; ++inode) route[inode] = nodes[inode];
    
    for(int inode = 0; inode < N; ++inode)
    {
      nodes[inode].setPassingFlag(false);
      route[inode].setPassingFlag(false);
    }

    //setbestroute(route);

    // ===========
    NodeClass *forward = new NodeClass[N];
    NodeClass *back = new NodeClass[N];

    NodeClass *forwardFirst = NULL, *backFirst = NULL;
    NodeClass *forwardLast = NULL, *backLast = NULL;

    for(int inode = 0; inode < N; ++inode)
    {
      forward[inode] = route[inode];    
      back[inode] = route[N - 1 - inode];
    }

    forwardFirst = &forward[0];
    forwardFirst->setPreNode(NULL);
    for(int inode = 1; inode < N; ++inode)
    {
      forwardFirst->setNextNode(&forward[inode]);
      forwardFirst->getNextNode()->setPreNode(forwardFirst);
      forwardFirst = forwardFirst->getNextNode();
      forwardFirst->setNextNode(NULL);
    }
    forwardFirst = &forward[0];
    forwardLast = &forward[N - 1]; 

    backFirst = &back[0];
    backFirst->setPreNode(NULL);
    for(int inode = 1; inode < N; ++inode)
    {
      backFirst->setNextNode(&back[inode]);
      backFirst->getNextNode()->setPreNode(backFirst);
      backFirst = backFirst->getNextNode();
      backFirst->setNextNode(NULL);
    }
    backFirst = &back[0];
    backLast = &back[N - 1];

    setbestroute2(forwardFirst);
    // =============

    //distmin = mst->get2appDistMin();

    while( iteration < N )
    {
      dist = distmin;

//            for(int i = 0; i < N; ++i) route[i].setPassingFlag(false);

      double g = 0;
      improvePath(1, dist, g, forwardFirst, forwardLast, backFirst, backLast);

      if(dist < distmin) 
      {
	distmin = dist;
	cout << "distmin = " << distmin << endl;
	outputRoute();

	for(int i = 0; i < bestroute.size(); ++i)
	{
	  PathClass *tmpPath = &(bestroute.at(i));
	  route[i] = nodes[tmpPath->getFromNode()->getNo()];
	}
	iteration = 0;
      } else
      {
	++iteration;
      }

      /*NodeClass firstNode = route[0];
      for(int i = 0; i < N - 1; ++i)
      {
	route[i] = route[i + 1];
      }
      route[N - 1] = firstNode;
*/
      forwardLast->setNextNode(forwardFirst);
      forwardFirst->setPreNode(forwardLast);
      
      forwardLast = forwardFirst;
      forwardFirst = forwardFirst->getNextNode();

      forwardFirst->setPreNode(NULL);
      forwardLast->setNextNode(NULL);

      backLast->setNextNode(backFirst);
      backFirst->setPreNode(backLast);
      backFirst = backLast;
      backLast = backLast->getPreNode();
      
      backFirst->setPreNode(NULL);
      backLast->setNextNode(NULL);


    }

    outputNode();
    outputRoute();

    //delete [] nodes;
    mst->free();
    delete mst;
    delete [] route;
  }

  int Linkernighan::improvePath(int lambda, double &dist, double &g, NodeClass *&forwardFirst, NodeClass *&forwardLast, NodeClass *&backFirst, NodeClass *&backLast)
  {

    int rval = 0;
    const int N = Num;
    const int alpha = 3;
    const int endnode = N - 1;
    double g2 = 0.0, tmpdist = 0.0;
    double distRoute = 0.0;
    
    //for(int i = 0; i < N; ++i) oldroute[i] = route[i];
    
    for(NodeClass *p = forwardFirst; p != NULL; p = p->getNextNode())
    {
      int fromIdx = p->getNo(); 
      int toIdx = p->getNextNode() != NULL ? p->getNextNode()->getNo() : forwardFirst->getNo(); 
      distRoute += pathmap[fromIdx][toIdx].getval();
    }

    NodeClass *pforward = forwardFirst, *pback = backLast;

    if(lambda < alpha)
    {

      for(int path = 0; path < N - 1; pforward = pforward->getNextNode(), pback = pback->getPreNode(), ++path)
      //for(map<int, int>::iterator ite = condidatemap[pforward->getNo()].begin();
	//  ite != condidatemap[pforward->getNo()].end(); pforward = pforward->getNextNode(), pback = pback->getPreNode())
      {
	//int path = ite->first;

	//if(pforward->getPassingFlag()) continue;
	if((pforward->getPassingFlag() || condidatemap[pforward->getNo()].find(forwardLast->getNo()) == condidatemap[pforward->getNo()].end() ) && 
	    condidatemap[pforward->getNo()].find(pforward->getNextNode()->getNo()) != condidatemap[pforward->getNo()].end()) continue;

	//PathClass wxy(&route[path], &route[path + 1]);
	//2PathClass wex(&route[endnode], &route[path]);
	PathClass wxy(pforward, pforward->getNextNode());
	PathClass wex(forwardLast, pforward);
	g2 = wxy.getval() - wex.getval(); 
	
	if(g + g2 > 0.0)
	//if(g2 > 0)
	{
	  //PathClass wes(&route[endnode], &route[0]);
	  //PathClass wys(&route[path + 1], &route[0]);
	  PathClass wes(forwardLast, forwardFirst);
	  PathClass wys(pforward->getNextNode(), forwardFirst);
	  
	  tmpdist = distRoute - wxy.getval() + wex.getval() - wes.getval() + wys.getval();

	  if(tmpdist < dist && fabs(tmpdist - dist) > 1.0e-10)
	  {
	    //for(int ipath = path + 1; ipath < N; ++ipath)
	     // oldroute[ipath] = route[endnode - ipath + path + 1];
	    
	  //  for(int ipath = path + 1; ipath < N; ++ipath)
	    //  route[ipath] = oldroute[ipath]; 

	    ///////////////////////////
	    NodeClass *pforwardnext = pforward->getNextNode();
	    NodeClass *pbackpre = pback->getPreNode(); 
	    
	    pforward->getNextNode()->setPreNode(NULL);
	    pforward->setNextNode(backFirst);
	    backFirst->setPreNode(pforward);
	 
	    pback->getPreNode()->setNextNode(NULL);
	    pback->setPreNode(forwardLast);
	    forwardLast->setNextNode(pback);   
	    ///////////////////////////

	    //setbestroute(route);
	    setbestroute2(forwardFirst);
	    
	    dist = tmpdist;

	    g = g + g2;

	    rval = 1;

	    forwardLast = pbackpre;
	    backFirst = pforwardnext;
	    goto CLEANUP;

	  } else 
	  {
#if 0
	    for(int ipath = 0; ipath < path + 1; ++ipath)
	      oldroute[ipath] = route[ipath];
	    for(int ipath = path + 1; ipath < N; ++ipath)
	      oldroute[ipath] = route[endnode - ipath + path + 1];
#endif
	
	    ///////////////////////////

	    NodeClass *pforwardnext = pforward->getNextNode();
	    NodeClass *pbackpre = pback->getPreNode(); 

	    pforward->getNextNode()->setPreNode(NULL);
	    pforward->setNextNode(backFirst);
	    backFirst->setPreNode(pforward);
	 
	    pback->getPreNode()->setNextNode(NULL);
	    pback->setPreNode(forwardLast);
	    forwardLast->setNextNode(pback);   

	    ///////////////////////////
    
	    //oldroute[path].setPassingFlag(true);
	    pforward->setPassingFlag(true);

	    g = g + g2;

	    rval = improvePath(lambda + 1,  dist, g, forwardFirst, pbackpre, pforwardnext, backLast);

	    //oldroute[path].setPassingFlag(false);
	    pforward->setPassingFlag(false);
	    
	    if(rval) 
	    {
	      forwardLast = pbackpre;
	      backFirst = pforwardnext;

	      goto CLEANUP;
	    } else
	    {
	     forwardFirst->setPreNode(NULL);
	     forwardLast->setNextNode(NULL);
	     pforward->setNextNode(pforwardnext);
	     pforwardnext->setPreNode(pforward);

	     backFirst->setPreNode(NULL);
	     backLast->setNextNode(NULL);
	     pback->setPreNode(pbackpre);
	     pbackpre->setNextNode(pback);
	    
	    }

	  }

	} //if g
      }
    } else 
    {

      int idxpath = -1;
      double dbl = -1.0e+10;

      NodeClass *tmppforward = NULL, *tmppback = NULL;
      
      for(int path = 0; path < N - 1; ++path, pforward = pforward->getNextNode(), pback = pback->getPreNode())
      {

	//if(route[path].getPassingFlag()) continue;
	//if(pforward->getPassingFlag()) continue;
	if(pforward->getPassingFlag() && condidatemap[pforward->getNo()].find(forwardLast->getNo()) == condidatemap[pforward->getNo()].end()) continue;



	//PathClass wxy(&route[path], &route[path + 1]);
	//PathClass wex(&route[endnode], &route[path]);
	PathClass wxy(pforward, pforward->getNextNode());
	PathClass wex(forwardLast, pforward);
	g2 = wxy.getval() - wex.getval(); 

	if(g2 > dbl)
	{
	  //pforward = pforward->getNextNode();
	  //pback = pback->getPreNode();
	  tmppforward = pforward;
	  tmppback = pback;
	  dbl = g2;
	  idxpath = path;
	}
      }

      pforward = tmppforward;
      pback = tmppback;

      if(dbl > 0 && idxpath >= 0)
      {
	
	tmpdist = 0;
	//PathClass wes(&route[endnode], &route[0]);
	//PathClass wys(&route[idxpath + 1], &route[0]);
	PathClass wes(forwardLast, forwardFirst);
	PathClass wys(pforward->getNextNode(), forwardFirst);

	tmpdist = distRoute - dbl - wes.getval() + wys.getval();

	if(tmpdist < dist && fabs(tmpdist - dist) > 1.0e-10)
	{
#if 0 
	  for(int ipath = idxpath + 1; ipath < N; ++ipath)
	    oldroute[ipath] = route[endnode - ipath + idxpath + 1];

	  for(int ipath = idxpath + 1; ipath < N; ++ipath)
	    route[ipath] = oldroute[ipath]; 
#endif
	  ///////////////////////////

	  NodeClass *pforwardnext = pforward->getNextNode();
	  NodeClass *pbackpre = pback->getPreNode(); 

	  pforward->getNextNode()->setPreNode(NULL);
	  pforward->setNextNode(backFirst);
	  backFirst->setPreNode(pforward);

	  pback->getPreNode()->setNextNode(NULL);
	  pback->setPreNode(forwardLast);
	  forwardLast->setNextNode(pback);   
	  ///////////////////////////

	  dist = tmpdist;
	  rval = 1;

	  //setbestroute(route);
	  setbestroute2(forwardFirst);
	  
	  forwardLast = pbackpre;
	  backFirst = pforwardnext;

	  g = g  + dbl;
	  
	  goto CLEANUP;
	} else 
	{
#if 0
	  for(int ipath = 0; ipath < idxpath + 1; ++ipath)
	    oldroute[ipath] = route[ipath];
	  for(int ipath = idxpath + 1; ipath < N; ++ipath)
	    oldroute[ipath] = route[endnode - ipath + idxpath + 1];
#endif
	  
	  ///////////////////////////
	  NodeClass *pforwardnext = pforward->getNextNode();
	  NodeClass *pbackpre = pback->getPreNode(); 

	  pforward->getNextNode()->setPreNode(NULL);
	  pforward->setNextNode(backFirst);
	  backFirst->setPreNode(pforward);

	  pback->getPreNode()->setNextNode(NULL);
	  pback->setPreNode(forwardLast);
	  forwardLast->setNextNode(pback);   
	  ///////////////////////////

	  //oldroute[idxpath].setPassingFlag(true);
	  pforward->setPassingFlag(true);
	  g = g + dbl;
	  //rval = improvePath(lambda + 1, oldroute, dist, g, forwardFirst, forwardLast, backFirst, backLast); 
	  rval = improvePath(lambda + 1, dist, g, forwardFirst, pbackpre, pforwardnext, backLast);

	  if(rval)
	  {
	    forwardLast = pbackpre;
	    backFirst = pforwardnext;
	  }else
	  {
	  
	    forwardFirst->setPreNode(NULL);
	     forwardLast->setNextNode(NULL);
	     pforward->setNextNode(pforwardnext);
	     pforwardnext->setPreNode(pforward);

	     backFirst->setPreNode(NULL);
	     backLast->setNextNode(NULL);
	     pback->setPreNode(pbackpre);
	     pbackpre->setNextNode(pback);
	    
	  }

	  //oldroute[idxpath].setPassingFlag(false);
	  pforward->setPassingFlag(false);
	  //if(rval)
	  {
	    //for(int i = 0; i < N; ++i) route[i] = oldroute[i];
 	    goto CLEANUP;	    
	  }

	}
      }
    }

CLEANUP:


    return rval;
  }

  void Linkernighan::outputNode()
  {
    const int N = Num;
    NodeClass* Nodes = getnodes();
    ofstream nodefile("nodes.dat");

    for(int i = 0; i < N; ++i)
      nodefile << Nodes[i].getPosx() << "\t" << Nodes[i].getPosy() << endl;

    nodefile.close();
    return;
  }

  void Linkernighan::outputRoute()
  {
    ofstream ofile("route2.gnu");

    ofile << "#!/bin/gnuplot" << endl;
    ofile << "set size square" << endl;
    ofile << "set xrange [0.0:1.0]" << endl;
    ofile << "set yrange [0.0:1.0]" << endl;
    ofile << "plot '-' w l,  'nodes.dat'" << endl;

    for(int i = 0; i < bestroute2.size(); ++i)
    {
      PathClass *tmpPath = &(bestroute2[i]);

      ofile << tmpPath->getFromNode()->getPosx() << " " << tmpPath->getFromNode()->getPosy() << endl;
      ofile << tmpPath->getToNode()->getPosx() << " " << tmpPath->getToNode()->getPosy() << endl;

      ofile << endl;

    }

    ofile.close();

  }

  void Linkernighan::setbestroute(NodeClass *route)
  {
    //cout << "set best route"<< endl;
    const int N = Num;
    
    bestroute.clear();
    double dist = 0.0;
    for(int i = 0; i < N; ++i)
    {
      int fromIdx = i;
      int toIdx = i != N - 1 ? i + 1 : 0;
      //PathClass tmppath(&route[fromIdx], &route[toIdx]); 
      int inode = route[fromIdx].getNo();
      int jnode = route[toIdx].getNo();
      PathClass tmppath(&nodes[inode], &nodes[jnode]);

      dist += tmppath.getval();
      bestroute.push_back(tmppath);
    }
    cout << "testdist = " << dist << endl;
    return;
  }

  void Linkernighan::setbestroute2(NodeClass *pforward)
  {
    const int N = Num;

    bestroute2.clear();
    double dist = 0.0;

    for(NodeClass *tmpNode = pforward; tmpNode != NULL; tmpNode = tmpNode->getNextNode())
    {
      int inode = tmpNode->getNo();
      int jnode = tmpNode->getNextNode() != NULL ? tmpNode->getNextNode()->getNo() : pforward->getNo();    
      PathClass tmppath(&nodes[inode], &nodes[jnode]);

      dist += tmppath.getval();
      bestroute2.push_back(tmppath);
    }

    cout << "TEST DISS " << dist << endl;
    //getchar();

    return;
  }
};
