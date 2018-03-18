#include <iostream>
#include "opt.h"
#include "utils.h"
#include "params.h"
#include "minimum_spanning_tree.h"

using namespace std;
using namespace OptimizeName;

int main()
{

  const int N = Num;
  MinimumSpanningTree *mst = new MinimumSpanningTree();
  NodeClass *nodes = new NodeClass[N];

  for(int i = 0; i < N; ++i)
  {
    double xp = (double) rand() / RAND_MAX;
    double yp = (double) rand() / RAND_MAX;

    nodes[i].setNo(i);
    nodes[i].setPosx(xp);
    nodes[i].setPosy(yp);
  }

  mst->setNodes(nodes);

  mst->initialize();

  mst->opt();

  mst->calc2appRoute();

  cout << N << " " << mst->get2appDistMin()<< endl;

  mst->free();

  delete [] nodes;
  delete mst;

  return 0;
}
