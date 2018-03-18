#include <iostream>
#include <vector>
#include "function.h"
#include "opt.h"
#include "params2.h"

using namespace std;
using namespace OptimizeName;

int main()
{

  for(int i = 0; i < Nx; ++i)
  {
    Function tmp;
    for(int j = 0; j < (Nx+1);++j)
    {
//      tmp.Coeff[j].val = (double) rand() / RAND_MAX;
    }
  }

  cout << "hoge" << endl;getchar();
#if 0
  clock_t time1 = clock();
  double sum = 0.0;
  int cc=0;

  for(int i = 0; i < (Nx+1); ++i)
  {
    sum += tmp.Coeff[i].val;
    ++cc;
    if(cc == Nx) break;
  }
  time1 = clock() - time1;
  double t1 = (double) time1 / CLOCKS_PER_SEC;
  cout << t1 << endl;
  sum = 0.0;

  time1 = clock();

  //tmp.conectChain();
  tmp.firstCoeff = &tmp.Coeff[0];
  tmp.firstCoeff->nextCoeff =&tmp.Coeff[Nx];
  tmp.Coeff[Nx].nextCoeff = &tmp.Coeff[0];

  int counter = 0;
  for(coeff *Coeff = tmp.firstCoeff;  Coeff != NULL; Coeff = Coeff->nextCoeff)
  {
    sum += Coeff->val;
    ++counter;
    if(counter == Nx) break;
  }
  time1 = clock() - time1;
  t1 = (double) time1 / CLOCKS_PER_SEC;
  cout << t1 << endl;
  getchar();
#endif
  return 0;
}
