#include <iostream>
#include "simplex.h"
#include "params.h"

using namespace OptimizeName;

/* -- main --- */
int main(int argc, char *argv[])
{
  Simplex *simplex = new Simplex();
  simplex->initialize();
  simplex->opt();

  delete simplex;
  
  return 0;
}

