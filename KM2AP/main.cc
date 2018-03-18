#include <iostream>
#include "params.h"
#include "utils.h"
using namespace std;

void read_params(string &filename, paramsClass &params);
int km2ap(paramsClass &params);

int main(int argc, char *argv[])
{
  int rval = 0;
  paramsClass params;
  string str; 
  if(argc < 2) 
  {
    cout << "Set params file " << endl;
    goto CLEANUP;
  }

  str = argv[1];
  read_params(str, params);
  rval = km2ap(params);

CLEANUP:

  return 0;

}
