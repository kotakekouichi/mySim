#include <iostream>
#include <vector>
#include "../lib/utils.h"


using std::cout;
using std::endl;
using std::vector;
using std::string;

namespace VehicleAllocation
{

  class storeinfo;
  class trackinfo;

  class trackinfo 
  {
    public:
      int index;
      double estimateval;
      std::vector<storeinfo> StoreInfo;
      trackinfo()
      {}

      trackinfo(int index, double estimateval) : index(index), estimateval(estimateval)
      {
	this->StoreInfo.reserve(100);
      }
  };
  
  class storeinfo
  {
    public:
      int istore;
      int volume;
      std::string ordercode;
      storeinfo(){}
  };

  void hoge()
  {
    std::cout << "hoge" << std::endl;
  }

  int adaptive()
  {
    vector<trackinfo> Track;

    for(int itrack = 0; itrack < 19; itrack++)
    {
      //trackinfo tmp(itrack, rnd());

      double val = rnd();;
      //Track.push_back(tmp);
      //Track.push_back(tmp);
    }

    return 0;
  }

};

using namespace VehicleAllocation;

int main()
{
 
  return 0;
}


