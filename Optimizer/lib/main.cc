#include <iostream>
#include "utils.h"

using namespace std;

int main()
{
  //read_params("./params,dat");
  
  init_rnd();
  double val = 0.0;
  rnd1(val);
  cout << val << endl;
  return 0;
}
