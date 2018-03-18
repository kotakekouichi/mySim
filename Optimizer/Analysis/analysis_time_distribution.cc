#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

class store
{
  public:
    int No;
    int stime;
    int etime;
    int volume;
    int reduction;
    double listval;
    int rot[5];
};

void time_distribution(const char filename[]);
void inputfile(const char filename[], vector<store> &Store);
void outputfile(vector<store> &Store);
vector<string> split(const string &s, char delim[]);
inline int StringToInt(const string str);
inline double StringToDouble(const string str);

int main(int argv, char *argc[])
{

  time_distribution(argc[1]);

  return 0;
}

void time_distribution(const char filename[])
{

  vector<store> Store;

  inputfile(filename, Store);
  outputfile(Store);

}

void inputfile(const char filename[], vector<store> &Store)
{

  char delim[] = "\t:";
  string buf;
  vector<string> elems;
  ifstream ifs(filename);

  getline(ifs, buf);
  while(getline(ifs, buf))
  {
    store vec;
    elems = split(buf, delim);
    vec.No = StringToInt(elems[0]);
    vec.stime = StringToInt(elems[1]) * 60 + StringToInt(elems[2]);
    vec.etime = vec.stime + 120;
    vec.volume = StringToInt(elems[5]);
    vec.reduction = StringToInt(elems[6]);
    vec.listval = StringToDouble(elems[7]);

    for(int i = 0; i < 5; ++i)
      vec.rot[i] = StringToInt(elems[i + 8]);

    Store.push_back(vec);
  }

}

void outputfile(vector<store> &Store)
{
  int time;
  ofstream ofs("./time_distribution.dat");

  for(int i =0; i < (int) Store.size();++i)
  {
    time = (Store[i].stime + Store[i].etime) / 2;
    ofs << Store[i].No << "\t" << time << "\t" << Store[i].listval << "\t" << Store[i].volume  
      << "\t" << Store[i].reduction << "\t" << Store[i].rot[0] * Store[i].rot[1] << std::endl; 
  }

  ofs.clear();
}
  
vector<string>split(const string &s, char delim[])
{
  vector<string> v;
  //char delim[] = {'=', ' ', '\t'};
  string item;
  int i=0;
  while(i != s.size()){
    char ch = s[i];
    if (ch == delim[0] || ch == delim[1]) {
      if (!item.empty())
	v.push_back(item);
      item.clear();
    }
    else {
      item += ch;
    }
    i++;
  }
  if (!item.empty())
    v.push_back(item);

  return v;

}

inline int StringToInt(const string str){return atoi(str.c_str());}
inline double StringToDouble(const string str)
{
  istringstream is;
  is.str(str);
  double x;
  is >> x;
  return x;
}
