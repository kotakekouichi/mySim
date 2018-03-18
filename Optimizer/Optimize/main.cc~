#include <iostream>
#include "simplex.h"
#include "simulated_annealing.h"
#include "genetic_algorithms.h"
#include "interior_point_method.h"
#include "chaos_optimization.h"
#include "minimum_spanning_tree.h"
#include "branch_and_cut.h"
#include "params.h"

using namespace OptimizeName;

/* -- main --- */
int main(int argc, char *argv[])
{
#if 1  
  OptClass *branchAndCut = new BranchAndCut();
  branchAndCut->initialize();  
  branchAndCut->opt();
  delete branchAndCut;
#endif 
/*
  Simplex *simplex = new Simplex();
  simplex->initialize();
  simplex->opt();
*/
  //OptClass *simplex = new Simplex();
  /*
  OptClass *mst = new MinimumSpanningTree();
  NodeClass *nodes = new NodeClass[10];
  for(int i =0 ; i < 10; ++i)
  {
    int No = i;
    double xp = (double) rand() / (double) RAND_MAX;
    double yp = (double) rand() / (double) RAND_MAX;
    nodes[i].setNo(i);
    nodes[i].setPosx(xp);
    nodes[i].setPosy(yp);
  }
  ((MinimumSpanningTree*)mst)->setNodes(nodes);
  mst->initialize();
  mst->opt();
  ((MinimumSpanningTree*)mst)->outputNodes();
  ((MinimumSpanningTree*)mst)->outputmst();
  ((MinimumSpanningTree*)mst)->calc2appRoute();
  cout << ((MinimumSpanningTree*)mst)->get2appDistMin();
  ((MinimumSpanningTree*)mst)->output2appRoute();

  delete [] nodes;
  delete mst;
  */
/*
  simplex->SetObjective("2x1+x2-x3");
  simplex->AddConstraint("x1 - 2x2+2x3 < -4");
  simplex->AddConstraint("- 2x2-x3 < -1");
  simplex->AddConstraint("2x1 - 3x2+1x3 < 1");

  */
/*
  simplex->SetObjective("2x0+x1-x2");
  simplex->AddConstraint("x0 - 2x1+2x2 < -4");
  simplex->AddConstraint("- 2x1-x2 < -1");
  simplex->AddConstraint("2x0 - 3x1+1x2 < 1");
  */
  //simplex->run();


#if 1
  OptClass *simplex = new Simplex();
  simplex->initialize();
  ((Simplex*)simplex)->optRevised();
  delete simplex;
#endif

  //simplex->getFeasibleInitilize(x);
  
  //tmp.initialize();
  //tmp.run();
  //Simulated_Annealing::run();
  //Genetic_Algorithm::run();
  //Interior_Point_Method::run();
  //Chaos_Optimization::run();

  return 0;
}

