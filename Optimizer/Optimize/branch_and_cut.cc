#include <iostream>
#include "params.h"
#include "utils.h"
#include "branch_and_cut.h"
#include "simplex.h"
#include "minimum_spanning_tree.h"

using namespace std;

namespace OptimizeName
{

  void BranchAndCut::initialize()
  {
    cout << "start initialize " << endl;

    Function *pfobj = getObjectiveFunction();
    pfobj->setobject(minimization);
    int idx = 0, ii = 0;

    xVariableIdx = new int[Num*Num];
    nodes = new NodeClass[Num];
    index2Nums = new int*[Nx];
    for(int i = 0; i < Nx; ++i) index2Nums[i] = new int[2];

    while(ii < Num)
    {
      double xp = (double) rand() / (double) RAND_MAX;
      double yp = (double) rand() / (double) RAND_MAX;
#if 1
      double l = sqrt(sq(xp - 0.5) + sq(yp - 0.5));
      if(!(0.42 < l && l < 0.48)) 
      {
	continue;
      }
#endif

      nodes[ii].setNo(ii);
      nodes[ii].setPosx(xp);
      nodes[ii].setPosy(yp);
      ++ii;
    }

    // x_idx = x_ij
    for(int i = 0; i < Num; ++i)
    {
      for(int j = i + 1; j < Num; ++j)
      {
	index2Nums[idx][0] = i;
	index2Nums[idx][1] = j;
	xVariableIdx[i + j*Num] = idx;
	++idx;
      }
    }

    //for TSP
    MinimumSpanningTree mst;
    mst.setNodes(nodes);
    mst.initialize();
    mst.opt();
    mst.calc2appRoute();
    vector<PathClass> *routevec = mst.get2appRoute();

    this->setInitRouteVec(*routevec);

    for(int i = 0; i < routevec->size(); ++i)
    {
      int fromIdx = routevec->at(i).getFromNode()->getNo();
      int toIdx = routevec->at(i).getToNode()->getNo(); 
      int index = fromIdx < toIdx ? fromIdx + Num * toIdx : toIdx + Num * fromIdx;
      index = xVariableIdx[index];
      initrouteIdxmap[index] = index;
    }

    cout << "end initialize " << endl;
    return;
  }


  void BranchAndCut::initializeFunction()
  {
    int idx = 0; 

    Function *pfobj = getObjectiveFunction();
    vector<Function> *gfunc = getConstraintFunction();
    vector<Function> *partialvec = getPartialRouteConstraint();

    pfobj->constTermVal = 0.0;
    for(int i = 0; i < Nx; ++i)
    {
      //pfobj->Coeff[i].val = val[i];
      int fromIdx = index2Nums[i][0], toIdx = index2Nums[i][1];
      double cij = sqrt(sq(nodes[fromIdx].getPosx() - nodes[toIdx].getPosx()) + sq(nodes[fromIdx].getPosy() - nodes[toIdx].getPosy()));
      if(initrouteIdxmap.find(i) == initrouteIdxmap.end())
      {
	pfobj->setCoeff(i, cij);
      }
      else
      {
	int index = Nx +  2 * Num + i; // x_index = 1 - xij >= 0  
	pfobj->constTermVal += cij;
	pfobj->setCoeff(index, -cij);
      }
    }

    //次数LP緩和
    for(int i = 0; i < Num; ++i)
    {
      Function tmp1, tmp2;
      tmp1.ix = gfunc->size() + Nx;
      tmp2.ix = gfunc->size() + Nx + 1;

      tmp1.constTermVal = 0;
      tmp2.constTermVal = 0;

      for(int k = 0; k < Num; ++k)
      {
	if(k < i)
	{
	  idx = xVariableIdx[k + i*Num];
	}
	else if(k > i)
	{ 
	  idx = xVariableIdx[i + k * Num]; 
	}
	else 
	  continue;

	if(initrouteIdxmap.find(idx) == initrouteIdxmap.end())
	{
	  tmp1.setCoeff(idx, 1);
	  tmp2.setCoeff(idx, -1);	  
	}
	else 
	{
	  int index = Nx +  2 * Num + idx; // x_index = 1 - xij >= 0  
	  tmp1.setCoeff(index, -1);
	  tmp2.setCoeff(index, 1);
	}
      }
      gfunc->push_back(tmp1);
      gfunc->push_back(tmp2);
    } 

    for(int i = 0; i < Nx; ++i)
    {
      Function tmp;

      if(initrouteIdxmap.find(i) == initrouteIdxmap.end())
      {
	tmp.ix = gfunc->size() + Nx; 
	tmp.constTermVal = 1; 
	tmp.setCoeff(i, -1.0);
      }
      else 
      {
	tmp.ix = i;
	tmp.constTermVal = 1;
	tmp.setCoeff(gfunc->size() + Nx, -1.0);       
      }

      //tmp.setCoeff(i, -1.0);
      gfunc->push_back(tmp);
    }

    for(int i = 0; i < partialvec->size(); ++i)
    {
      gfunc->push_back(partialvec->at(i));
    }

    return;

  }


  void BranchAndCut::opt()
  {
    Function *pfobj = getObjectiveFunction();
    vector<Function> *gfunc = getConstraintFunction();

    //initializeFunction();

    while(1)
    {

      vector<Function> newPartialVec;

      initializeFunction();
      
      OptClass *simplex = new Simplex(pfobj, gfunc);

      ((Simplex *)simplex)->optRevised();
      //simplex->opt();

      getchar();

      outputRoute(simplex);

      // 部分巡回路消去
      partialRouteConstraint(newPartialVec, simplex);
      addPartialRouteConstraint(newPartialVec);

      //for(int i = 0; i < newPartialVec.size(); ++i)
	//gfunc->push_back(newPartialVec[i]);
      
      if(newPartialVec.size() == 1) break;

      gfunc->clear();
      pfobj->clearCoeff();

      delete simplex;

    }

    return;
  }

  void BranchAndCut::partialRouteConstraint(vector<Function> &tmpvec, OptClass *Opt)
  {

    vector<Function> *gfunc = Opt->getConstraintFunction();
    const int Nconstraint = gfunc->size();
    int fromIdx = 0, toIdx = 0, ix = 0;
    int **Index2Nums = getIndex2Nums();

    UnionFinding uf(Num);

    for(int ic = 0; ic < Nconstraint; ++ic)
    {
      ix = gfunc->at(ic).ix;
      if(gfunc->at(ic).constTermVal == 0.0) continue;
      if(!(ix < Nx)) continue;

      fromIdx = Index2Nums[ix][0];
      toIdx = Index2Nums[ix][1];

      uf.Union(fromIdx, toIdx);
    }

    map<int, map<int, int> > counterMap;
    for(int i = 0;i < uf.par.size(); ++i)
    {
      counterMap[uf.root(i)][i] = i;
    }
    
    for(map<int, map<int, int> >::iterator oneloop = counterMap.begin(); oneloop != counterMap.end(); ++oneloop)
    {
      map<int, int> idxmap = oneloop->second;
      Function tmp1;
      tmp1.constTermVal = -2;
      tmp1.ix = gfunc->size() + Nx + tmpvec.size();

      int idx = 0;
      for(int i = 0; i < Num; ++i)
      {
	for(int j = i + 1; j < Num; ++j)
	{

	  if(idxmap.find(i) != idxmap.end() && idxmap.find(j) != idxmap.end())
	  {
	    ++idx;
	    continue; 
	  }
	  if(idxmap.find(i) == idxmap.end() && idxmap.find(j) == idxmap.end()) 
	  {
	    ++idx;
	    continue;
	  }

	  if(initrouteIdxmap.find(xVariableIdx[i+j*Num]) == initrouteIdxmap.end())
	  {
	    tmp1.setCoeff(idx, 1.0);
	  }
	  else
	  {
	    tmp1.constTermVal += 1.0;
	    tmp1.setCoeff(Nx + 2*Num + xVariableIdx[i + j * Num], -1.0); 
	  }

	  ++idx;
	} 
      }
      tmpvec.push_back(tmp1);
    }

#if 0
    for(int i = 0;i < tmpvec.size(); ++i)
    {
      for(int j = 0; j < gfunc->size(); ++j)
      {
	Function *constFunc = &gfunc->at(j);
	int ix = constFunc->ix;
        if(!tmpvec[i].findCoeff(ix)) continue;

	tmpvec[i].SwapVariable(ix, ix, *constFunc);
      
      }
    } 
#endif 
    return;
  }

  void BranchAndCut::addPartialRouteConstraint(vector<Function> &value)
  {
    vector<Function> *gfunc = this->getPartialRouteConstraint();
    for(int i = 0; i < value.size(); ++i)
    {
      gfunc->push_back(value[i]);
    }
  }

  //解を出力
  void BranchAndCut::outputRoute(OptClass *Opt)
  {
    vector<Function> *gfunc = Opt->getConstraintFunction();
    int fromIdx, toIdx, ix;
    int **Index2Nums = getIndex2Nums();
    const int Nconstraint = gfunc->size();
    ofstream ofile("outputRoute.gnu");

    ofile << "#!/bin/gnuplot" << endl;
    ofile << "set xrange[-0.05:1.05]" << endl;
    ofile << "set yrange[-0.05:1.05]" << endl;
    ofile << "set size square" << endl;
    ofile << "plot  'nodes.dat', '-' w l, " << endl;

    for(int i = 0; i < Nconstraint; ++i)
    {
      ix = gfunc->at(i).ix;
      //if(gfunc->at(i).constTermVal == 0.0) continue;
      if(fabs(gfunc->at(i).constTermVal) < eps) continue;
      if(!(ix< Nx)) continue;

      fromIdx = Index2Nums[ix][0]; 
      toIdx = Index2Nums[ix][1];

      ofile << nodes[fromIdx].getPosx() << " " << nodes[fromIdx].getPosy() << endl;
      ofile << nodes[toIdx].getPosx() << " " << nodes[toIdx].getPosy() << endl;
      ofile << endl;
      ofile << endl;
    }

    cout << "output ==> outputRoute.gnu" << endl;
    ofile.close();

    ofstream nodefile("nodes.dat");
    for(int i = 0; i < Num; ++i){nodefile << nodes[i].getPosx() << " " << nodes[i].getPosy() << endl;}
    nodefile.close();
    
    return;
  }

  //メモリの解放
  void BranchAndCut::free()
  {
    Function *pfobj = getObjectiveFunction();
    vector<Function> *gfunc = getConstraintFunction();
    int **Index2Nums = getIndex2Nums();

    for(int i = 0; i < Num; ++i)
      delete [] Index2Nums[i];
    delete [] Index2Nums;
    delete [] xVariableIdx;

    pfobj->clearCoeff();
    for(int i = 0; i < gfunc->size(); ++i)
      gfunc->at(i).clearCoeff();

    delete pfobj;
    gfunc->clear();
    delete nodes;
  }

};

