#include <iostream>
#include <vector>

using namespace std;

class UnionFinding
{
  public:
    vector<int> par;

    bool Union(int x, int y)
    {
       y = par[y]; 
       par[y] = x;
       return true;
    }

    bool Find(int x, int y)
    {
      return root(x) == root(y);
    }

    int root(int x)
    {
      if(par[x] == x) return x;
      par[x] = root(par[x]);
      return par[x];
    }
};

int main()
{
  UnionFinding uf;
  for(int i = 0; i < 10; ++i)
  {
    uf.par.push_back(i);
  }


  bool flag  = uf.Union(0, 1);
  flag = uf.Union(0, 2);
  flag = uf.Union(3,4);
  flag = uf.Union(0, 3);
  /*flag = uf.Union(5, 6);
  flag = uf.Union(5, 7);
  flag = uf.Union(5, 8);
  flag = uf.Union(5, 9);
  flag = uf.Union(0,9);
  */

  for(int i = 0; i < 10; ++i)
    cout << uf.root(i) << endl;
  return 0;
}
