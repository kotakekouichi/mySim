#include <iostream>
#include "node.h"
#include "linkernighan.h"
#include "params.h"
#include "minimum_spanning_tree.h"

using namespace std;

namespace OptimizeName
{

  map<int, int> forwardIdxMap;
  map<int, int> backIdxMap;

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

      if(sq(xp - 0.5) + sq(yp - 0.5) < 0.5 * 0.5  && 0.4 * 0.4 < sq(xp - 0.5) + sq(yp - 0.5))
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

	sortedmap[dbl] = nodes[j].getNo();

      }

      int c = 0;

      map<int, int> tmpmap;
      for(map<double, int>::iterator ite = sortedmap.begin(); ite != sortedmap.end(); ++ite, ++c)
      {

	int condidateIdx = ite->second;
	tmpmap[condidateIdx] = condidateIdx;

	//if(c == N / 5) break;
	if(c == 20) break;

      }
      condidatemap[nodes[i].getNo()] = tmpmap;
    }

    return;
  }

  void Linkernighan::free()
  {
    //cout << "free" << endl;
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

    //nodes = new NodeClass[N];
    for(int inode = 0; inode < N; ++inode)
    {
      route[inode] = nodes[inode];
    }

    for(int inode = 0; inode < N; ++inode)
    {
      nodes[inode].setPassingFlag(false);
      route[inode].setPassingFlag(false);
    }

    //setbestroute(route);

    // ===========
    forwardList = new NodeClass[N];
    backList = new NodeClass[N];
    NodeClass *oldforward = new NodeClass[N];
    NodeClass *oldback = new NodeClass[N];

    NodeClass *forwardFirst = NULL, *backFirst = NULL;
    NodeClass *forwardLast = NULL, *backLast = NULL;

    for(int inode = 0; inode < N; ++inode)
    {
      forwardList[inode] = route[inode];    
      backList[inode] = route[N - 1 - inode];
    }

    readyLinkerniham(...);
    
    forwardFirst = &forwardList[0];
    forwardFirst->setPreNode(NULL);
    for(int inode = 1; inode < N; ++inode)
    {
      forwardFirst->setNextNode(&forwardList[inode]);
      forwardFirst->getNextNode()->setPreNode(forwardFirst);
      forwardFirst = forwardFirst->getNextNode();
      forwardFirst->setNextNode(NULL);
    }
    forwardFirst = &forwardList[0];
    forwardLast = &forwardList[N - 1]; 

    backFirst = &backList[0];
    backFirst->setPreNode(NULL);
    for(int inode = 1; inode < N; ++inode)
    {
      backFirst->setNextNode(&backList[inode]);
      backFirst->getNextNode()->setPreNode(backFirst);
      backFirst = backFirst->getNextNode();
      backFirst->setNextNode(NULL);
    }
    backFirst = &backList[0];
    backLast = &backList[N - 1];

    setbestroute2(forwardFirst);
    // =============

    //distmin = mst->get2appDistMin();
    clock_t time = clock();

    distmin = 0.0;
    for(NodeClass *p = forwardFirst; p != NULL; p = p->getNextNode())
    {
      int fromIdx = p->getNo(); 
      int toIdx = p->getNextNode() != NULL ? p->getNextNode()->getNo() : forwardFirst->getNo(); 
      distmin += pathmap[fromIdx][toIdx].getval();
    }
    
    while( iteration < N )
    {
      dist = distmin;

      //for(int i = 0; i < N; ++i) route[i].setPassingFlag(false);

      // set start-end path

      routes_opt.clear();
      routeInfoClass tmpinfo;
      const int firstIdx = 0;
      const int lastIdx = N - 1; 
      tmpinfo.fowardFlag = true;
      tmpinfo.backFlag = false;
      tmpinfo.setFromNode(&forwardList[firstIdx]);
      tmpinfo.setToNode(&forwardList[lastIdx]);
      routes_opt.push_back(tmpinfo);

      // set key index of array
      forwardIdxMap.clear();
      backIdxMap.clear();

      for(int i = 0; i < N; ++i)
      {
        forwardIdxMap[forwardList[i].getNo()] = i;
	backIdxMap[backList[i].getNo()] = i;
      }

      // improve path by k-opt
      double g = 0;
      //cout << "improvePath " << endl;
      improvePath(1, distmin, dist, g, forwardFirst, forwardLast, backFirst, backLast);

      // found improved plan ? 
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
	++iteration; // next iteration
      }

      /*NodeClass firstNode = route[0];
	for(int i = 0; i < N - 1; ++i)
	{
	route[i] = route[i + 1];
	}
	route[N - 1] = firstNode;
       */
      
      // change start node to second node
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

      int idx = 0;
      for(NodeClass *pfirst = forwardFirst; pfirst != NULL; pfirst = pfirst->getNextNode(), ++idx)
      {
        oldforward[idx].setNo(pfirst->getNo()); 
	oldforward[idx].setPosx(pfirst->getPosx());
	oldforward[idx].setPosy(pfirst->getPosy());
	oldforward[idx].setPassingFlag(false);
	oldforward[idx].setPreNode(NULL);
	oldforward[idx].setNextNode(NULL);
      }

      idx = 0;
      for(NodeClass *pfirst = backFirst; pfirst != NULL; pfirst = pfirst->getNextNode(), ++idx)
      {
      	oldback[idx].setNo(pfirst->getNo()); 
	oldback[idx].setPosx(pfirst->getPosx());
	oldback[idx].setPosy(pfirst->getPosy());
	oldback[idx].setPassingFlag(false);
	oldback[idx].setNextNode(NULL);
	oldback[idx].setPreNode(NULL);
      }

      for(idx = 0; idx < N; ++idx)
      {
        forwardList[idx] = oldforward[idx];
	backList[idx] = oldback[idx];
      }

      forwardFirst = &forwardList[0];
      forwardFirst->setPreNode(NULL);
      for(int inode = 1; inode < N; ++inode)
      {
	forwardFirst->setNextNode(&forwardList[inode]);
	forwardFirst->getNextNode()->setPreNode(forwardFirst);
	forwardFirst = forwardFirst->getNextNode();
	forwardFirst->setNextNode(NULL);
      }
      forwardFirst = &forwardList[0];
      forwardLast = &forwardList[N - 1]; 

      backFirst = &backList[0];
      backFirst->setPreNode(NULL);
      for(int inode = 1; inode < N; ++inode)
      {
	backFirst->setNextNode(&backList[inode]);
	backFirst->getNextNode()->setPreNode(backFirst);
	backFirst = backFirst->getNextNode();
	backFirst->setNextNode(NULL);
      }
      backFirst = &backList[0];
      backLast = &backList[N - 1];

      //getchar();
      //cout << "go next" << endl;
    }

    outputNode();
    outputRoute();

    time = clock() - time;
    cout << "time[sec] : " << (double) time / CLOCKS_PER_SEC << endl;

    //delete [] nodes;
    mst->free();
    delete mst;
    delete [] route;
  }

  int Linkernighan::improvePath(int lambda, double distRoute, double &dist, double &g, NodeClass *&forwardFirst, NodeClass *&forwardLast, 
      NodeClass *&backFirst, NodeClass *&backLast)
  {

    int rval = 0;
    const int N = Num;
    const int alpha = 5;
    const int endnode = N - 1;
    double g2 = 0.0, tmpdist = 0.0;

    NodeClass *pforward = forwardFirst, *pback = backLast;

    if(lambda < alpha)
    {

      for(map<int, int>::iterator ite = condidatemap[forwardLast->getNo()].begin(); ite != condidatemap[forwardLast->getNo()].end(); ++ite)
      {
	//int path = ite->first;

	int path = forwardIdxMap[ite->first];
	pforward = &forwardList[path];

	path = backIdxMap[ite->first];
	pback = &backList[path];

	if(pforward->getPassingFlag()) continue;
	//if((pforward->getPassingFlag() || condidatemap[pforward->getNo()].find(forwardLast->getNo()) == condidatemap[pforward->getNo()].end() ) && 
	  //  condidatemap[pforward->getNo()].find(pforward->getNextNode()->getNo()) != condidatemap[pforward->getNo()].end()) continue;

	bool f_flag = false, b_flag = false;

	int path2 = pforward->getNo();
	int counter = 0;
	routeInfoClass target;

	for(int ir = 0; ir < routes_opt.size(); ++ir, ++counter)
	{
	  routeInfoClass r = routes_opt[ir];
	  if(r.fowardFlag)
	  {
	    const int iii = forwardIdxMap[path2];
	    const int fromIdx = forwardIdxMap[r.getFromNode()->getNo()];
	    const int toIdx = forwardIdxMap[r.getToNode()->getNo()];

	    //if(fromIdx < path2 && path2 < toIdx) 
	    if(fromIdx < iii && iii < toIdx)
	    {
	      f_flag = true;
	      b_flag = false;
	      target = r;
	    }
	  }

	  if(r.backFlag)
	  {
	    const int iii = backIdxMap[path2];
	    const int fromIdx = backIdxMap[r.getFromNode()->getNo()];
	    const int toIdx = backIdxMap[r.getToNode()->getNo()];

	    //if(toIdx < path2 && path2 < fromIdx) 
	    //if(fromIdx < path2 && path2 < toIdx)
	    if(toIdx < iii && iii < fromIdx)
	    //if(fromIdx < iii && iii < toIdx)
	    {
	      f_flag = false;
	      b_flag = true; 
	      target = r;
	    }
	  }
	}

	if(!f_flag && !b_flag) continue;

	vector<routeInfoClass> tmpRoute = routes_opt;
	routes_opt.clear();

	for(int ii = 0; ii < counter - 1; ++ii)
	{
	  routes_opt.push_back(tmpRoute[ii]);
	}

	routeInfoClass r1; r1.setFromNode(target.getFromNode()); r1.setToNode(f_flag ? pforward : pback);

	r1.fowardFlag = false; r1.backFlag = false;
	if(f_flag == true) r1.fowardFlag = true;
	if(b_flag == true) r1.backFlag = true;

	routes_opt.push_back(r1);

	for(int ii = routes_opt.size() - 1; ii > counter + 1; --ii)
	{
	  routeInfoClass *pnode = &tmpRoute[ii];
	  routeInfoClass r2; r2.setFromNode(pnode->getToNode()); r2.setToNode(pnode->getFromNode());

	  if(pnode->fowardFlag){r2.fowardFlag = false; r2.backFlag = true;}
	  if(pnode->backFlag){r2.fowardFlag = true; r2.backFlag = false;}

	  routes_opt.push_back(r2);

	}

	routeInfoClass r3; r3.setFromNode(target.getToNode()); r3.setToNode(f_flag ? pforward->getNextNode() : pback->getNextNode());
	r3.fowardFlag = false; r3.backFlag = false;

	if(f_flag == true) r3.backFlag = true;
	if(b_flag == true) r3.fowardFlag = true;

	NodeClass *pnode = f_flag ? pforward : pback;
	NodeClass *pnodenext = pnode->getNextNode();
	NodeClass *plast = forwardLast;
	NodeClass *pnodeback = f_flag ? pback : pforward;
	PathClass wxy(pnode, pnodenext);
	PathClass wex(plast, pnode);

	g2 = wxy.getval() - wex.getval(); 

	if(g + g2 > 0.0)
	{
	  PathClass wes(forwardLast, forwardFirst);
	  PathClass wys(pnodenext, forwardFirst);

	  tmpdist = distRoute - wxy.getval() + wex.getval() - wes.getval() + wys.getval();

	  if(tmpdist < dist && fabs(tmpdist - dist) > 1.0e-10)
	  {
	    NodeClass *pforwardnext = pnode->getNextNode();
	    NodeClass *pbackpre = pnodeback->getPreNode();

	    pnodenext->setPreNode(NULL);
	    pnode->setNextNode(backFirst);
	    backFirst->setPreNode(pnode);
 
	    ///pforward->getNextNode()->setPreNode(NULL);
	    ///pforward->setNextNode(backFirst);
	    ///backFirst->setPreNode(pforward);

	    pbackpre->setNextNode(NULL);
	    pnodeback->setPreNode(forwardLast);
	    forwardLast->setNextNode(pnodeback);

	    //pback->getPreNode()->setNextNode(NULL);
	    //pback->setPreNode(forwardLast);
	    //forwardLast->setNextNode(pback);   
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

	    //NodeClass *pforwardnext = pnode->getNextNode();//pforward->getNextNode();
	    //NodeClass *pbackpre = pnodeback->getPreNode();//pback->getPreNode(); 

	    NodeClass *pforwardnext = f_flag  ? pnode->getNextNode() : pnodeback->getNextNode() ; 
	    NodeClass *pbackpre = f_flag ? pnodeback->getPreNode() : pnode->getPreNode(); 


	    //
	    pnodenext->setPreNode(NULL);
	    pnode->setNextNode(backFirst);
	    backFirst->setPreNode(pnode);
	    //pforward->getNextNode()->setPreNode(NULL);
	    //pforward->setNextNode(backFirst);
	    //backFirst->setPreNode(pforward);

	    pbackpre->setNextNode(NULL);
	    //pback->setPreNode(forwardLast);
	    pnodeback->setPreNode(forwardLast);
	    forwardLast->setNextNode(pnodeback);
	
	    //pback->getPreNode()->setNextNode(NULL);
	    //pback->setPreNode(forwardLast);
	    //forwardLast->setNextNode(pback);   
	    ///////////////////////////

	    pforward->setPassingFlag(true);
	    pback->setPassingFlag(true);

	    g = g + g2;

	    rval = improvePath(lambda + 1, tmpdist, dist, g, forwardFirst, pbackpre, pforwardnext, backLast);

	    pforward->setPassingFlag(false);
	    pback->setPassingFlag(false);

	    if(rval) 
	    {
	      forwardLast = pbackpre;
	      backFirst = pforwardnext;

	      goto CLEANUP;
	    } else
	    {

	      forwardFirst->setPreNode(NULL);
	      forwardLast->setNextNode(NULL);
	      //pforward->setNextNode(pforwardnext);
	      //pforwardnext->setPreNode(pforward);
	      pnode->setNextNode(pnodenext);
	      pnodenext->setPreNode(pnode);

	      backFirst->setPreNode(NULL);
	      backLast->setNextNode(NULL);
	      //pback->setPreNode(pbackpre);
	      //pbackpre->setNextNode(pback);
	      pnodeback->setPreNode(pbackpre);
	      pbackpre->setNextNode(pnodeback);
	    }

	  }

	} //if g

	routes_opt.clear();
	routes_opt = tmpRoute;
      }
    } else //if(1)
    {

      int idxpath = -1;
      double dbl = -1.0e+10;

      NodeClass *tmppforward = NULL, *tmppback = NULL;
      vector<routeInfoClass> tmpRoute = routes_opt;

      //for(int path = 0; path < N - 1; ++path, pforward = pforward->getNextNode(), pback = pback->getPreNode())

      int counter2 = 0;
      routeInfoClass target2;
      bool f_flag2 = false, b_flag2 = false;
      for(map<int, int>::iterator ite = condidatemap[forwardLast->getNo()].begin(); ite != condidatemap[forwardLast->getNo()].end(); ++ite)
      {

	int path = forwardIdxMap[ite->first];
	pforward = &forwardList[path];

	path = backIdxMap[ite->first];
	pback = &backList[path];
	if(pforward->getPassingFlag() && condidatemap[pforward->getNo()].find(forwardLast->getNo()) == condidatemap[pforward->getNo()].end()) continue;


	bool f_flag = false, b_flag = false;

	int path2 = pforward->getNo();
	int counter = 0;
	routeInfoClass target;

	for(int ir = 0; ir < routes_opt.size(); ++ir, ++counter)
	{
	  routeInfoClass r = routes_opt[ir];
	  if(r.fowardFlag)
	  {
	    const int iii = forwardIdxMap[path2];
	    const int fromIdx = forwardIdxMap[r.getFromNode()->getNo()];
	    const int toIdx = forwardIdxMap[r.getToNode()->getNo()];

	    if(fromIdx < iii && iii < toIdx)
	    {
	      f_flag = true;
	      b_flag = false;
	      target = r;
	    }
	  }

	  if(r.backFlag)
	  {
	    const int iii = backIdxMap[path2];
	    const int fromIdx = backIdxMap[r.getFromNode()->getNo()];
	    const int toIdx = backIdxMap[r.getToNode()->getNo()];

	    if(toIdx < iii && iii < fromIdx)
	    {
	      f_flag = false; 
	      b_flag = true;
	      target = r;
	    }
	  }
	}

	if(!f_flag && !b_flag) continue;
	
	NodeClass *pnode = f_flag ? pforward : pback;
	NodeClass *pnodenext = pnode->getNextNode();
	NodeClass *plast = forwardLast;
	NodeClass *pnodeback = f_flag ? pback : pforward;

	PathClass wxy(pforward, pnodenext);
	PathClass wex(forwardLast, pforward);
	g2 = wxy.getval() - wex.getval(); 

	if(g2 > dbl)
	{
	  //tmppforward = pforward;
	  //tmppback = pback;
	  tmppforward = pnode;
	  tmppback = pnodeback;
	  dbl = g2;
	  idxpath = path;
	  counter2 = counter;
	  target2 = target;
	  f_flag2 = f_flag;
	  b_flag2 = b_flag;
	}
      }

      pforward = tmppforward;
      pback = tmppback;

      //===============
      routes_opt.clear();

      for(int ii = 0; ii < counter2 - 1; ++ii)
      { routes_opt.push_back(tmpRoute[ii]); }

      routeInfoClass r1; r1.setFromNode(target2.getFromNode()); r1.setToNode(f_flag2 ? pforward : pback);

      r1.fowardFlag = false; r1.backFlag = false;
      if(f_flag2 == true) r1.fowardFlag = true;
      if(b_flag2 == true) r1.backFlag = true;

      routes_opt.push_back(r1);

      for(int ii = routes_opt.size() - 1; ii > counter2 + 1; --ii)
      {
	routeInfoClass *pnode = &tmpRoute[ii];	    
	routeInfoClass r2; r2.setFromNode(pnode->getToNode()); r2.setToNode(pnode->getFromNode());

	if(pnode->fowardFlag){r2.fowardFlag = false; r2.backFlag = true;}
	if(pnode->backFlag){r2.fowardFlag = true; r2.backFlag = false;}

	routes_opt.push_back(r2);
      }

      routeInfoClass r3; r3.setFromNode(target2.getToNode()); r3.setToNode(f_flag2 ? pforward->getNextNode() : pback->getNextNode());
      r3.fowardFlag = false; r3.backFlag = false;

      if(f_flag2 == true) r3.backFlag = true;
      if(b_flag2 == true) r3.fowardFlag = true;
      //===========

      if(dbl > 0 && idxpath >= 0)
      {

	tmpdist = 0;
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
	  rval = improvePath(lambda + 1, tmpdist, dist, g, forwardFirst, pbackpre, pforwardnext, backLast);

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
    int counter = 0;

    bestroute2.clear();
    double dist = 0.0;

    for(NodeClass *tmpNode = pforward; tmpNode != NULL; tmpNode = tmpNode->getNextNode(), ++counter)
    {
      int inode = tmpNode->getNo();
      int jnode = tmpNode->getNextNode() != NULL ? tmpNode->getNextNode()->getNo() : pforward->getNo();    
      PathClass tmppath(&nodes[inode], &nodes[jnode]);

      dist += tmppath.getval();
      bestroute2.push_back(tmppath);
    }

    cout << "TEST DISS " << counter  << " " << dist << endl;

    return;
  }
};
