#include <iostream>
#include "matrix.h"

using namespace std;

int main()
{
  const int nrank = 6;
  AbstVector vec(nrank);

  AbstMatrix M(nrank, nrank);

  cout << M.getrowrank() << " " << M.getcolrank() << endl;
  for(int i = 0; i < nrank; ++i)
  {
    for(int j = 0; j < nrank; ++j)
    {
      double dbl = (double) rand() / RAND_MAX; 
      if(dbl < 0.4)
      {
        dbl = 0;
      }

      M.setComponent(i, j, dbl);
    }
  }

  for(int i = 0; i < nrank; ++i)
  {
    for(int j = 0; j < nrank; ++j)
    {
      int index = j + nrank * i;
      cout << M.getComponent(i, j) << "(" << index << "),";
    }
    cout << endl;
  }
#if 0
  M.setNonZero();
  cout << "-------" << endl;
  //for(int i = M.getFirstRowIndex(0, 0); i != -1; i = M.getNextNonZeroRowIndex(i, 0))
  for(int i = 0; i < nrank; ++i)
  {
    for(int j = M.getFirstColIndex(i, 0); j != -1; j = M.getNextNonZeroColIndex(i, j))
    {
      cout << M.getComponent(i, j) << "," ; 
    }
    cout << endl;
  }
#endif

  cout << "========" << endl;

  int irow = 0;
  for(map<int, double>::iterator ite = M.compmap.begin(); ite != M.compmap.end(); ++ite)
  {
    int index = ite->first;  
   if(index / M.getrank() > irow) {cout << endl; ++irow;}
    cout << "(" << ite->first << ") " << ite->second << ",";
  }
  cout << endl;
  return 0;
}
