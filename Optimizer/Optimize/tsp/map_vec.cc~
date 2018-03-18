#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include "utils.h"

using namespace std;

int main()
{

  cout << "map" << endl;
  for(long int N = 1; N < 10000000000; N *= 10)
  {
    map<string, int> Map;
    vector<string> Vec;

    for(int i = 0; i < N; ++i)
    {
      string str = "";
      int rnd = rand() % N;

      vartostring(str, rnd);
      Map[str] = rnd;
      Vec.push_back(str);
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

    for(int i = 0; i < strvec.size(); ++i)
    {
      if( Map.find(strvec[i]) != Map.end()){}
      //for(int j = 0; j < N; ++j){if(Vec[j] == str) break;}
    }
    cout << N << " " << (double) (clock() - time) / CLOCKS_PER_SEC << endl;

  }
  return 0;
}
