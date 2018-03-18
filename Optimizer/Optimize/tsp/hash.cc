#include <iostream>
#include <functional>
#include <vector> 
#include "utils.h"
#include <time.h>

using namespace std;

int main()
{

  cout << "hash" << endl;
  for(long int N = 1; N < 10000000000; N*=10)
  {
    myMap<string, int> Map;
    
    for(int i = 0; i < N; ++i)
    {
      string str = "";
      int rnd = rand() % N;
      vartostring(str, rnd);
      Map.add(str, rnd);

    }

    vector<string> strvec;
    for(int i = 0; i < N; ++i)
    {
      int rnd = rand() % (100*N);
      string str = "";
      vartostring(str, rnd);
      strvec.push_back(str);
    }
    
    clock_t time = clock();

    //for(int i = 0; i < N; ++i)
    for(int i = 0; i < strvec.size(); ++i)
    {
      if(Map.find(strvec[i]) != Map.end())
      {
      }
    }

    cout << N << " " << (double) (clock() - time) / CLOCKS_PER_SEC << endl;
  }


  return 0;
}
