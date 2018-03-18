#ifndef _UTILS_H
#define _UTILS_H
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

enum var_type
{
  INT,
  DOUBLE,
  STRING,
};

std::vector<string> split(const string &str)
{
  vector<string> v;
  char delim[] = {'=', '\t', ',', ' ', ':'};
  string item;
  
  for(int i = 0; i < (int) str.size(); ++i)
  {
    char ch = str[i];

    if(ch == delim[0] || ch == delim[1] || ch == delim[2] || ch == delim[3] || ch == delim[4])
    {
      if(!item.empty())
      {
	v.push_back(item);
      }
      item.clear();
    }
    else 
    {
      item += ch;
    }
  } 
  
  if(!item.empty())
  {
    v.push_back(item);
  }

  return v;
}

template<typename var>
inline void vartostring(string &str, var val)
{
  stringstream ss;
  ss << val;
  str = ss.str();
}

template<typename var>
inline void stringtovar(string str, var &val)
{
  istringstream (str) >> val;
}

template<typename var>
inline var sq(var x)
{
  return x * x;
}

template<typename var>
inline var cube(var x)
{
  return x * x * x;
}

template<typename var>
inline void rnd1(var &val)
{

  val = static_cast<var> (rand()) / (static_cast<var>(RAND_MAX));
  //return static_cast<var> (rand());
}

inline void init_rnd()
{
  srand((unsigned)time(NULL));
}

template<typename var>
inline void stringtimetovartime(string str, var &val)
{
  var hour;
  var minite;
  vector<string> time = split(str);

  val = hour * 60 +  minite;
} 

template<typename var>
inline void vartimetostringtime(string &str, var val)
{

  var hour = val / static_cast<var>(60.0);
  var minite = val % static_cast<var>(60.0);
  string strhour;
  string strminite;

  vartostring(strhour, hour);
  vartostring(strminite, minite);

  str = strhour + ":" + strminite;
}

struct adr_type
{
  void *var_adr;
  enum var_type type;
};

class params
{
  public:
    int A;
    double B;
    string C;
};

void read_params(string paramsfile)
{
  ifstream ifs(paramsfile);
  string line;
  map<string, adr_type> variable_map;
  adr_type at;
  params Params;

  at.type = INT;
  at.var_adr = &Params.A;
  variable_map["A"] = at;

  at.type = DOUBLE;
  at.var_adr = &Params.B;
  variable_map["B"] = at;

  at.type = STRING;
  at.var_adr = &Params.C;
  variable_map["C"] = at;

  while(getline(ifs, line))
  {
    vector<string> v = split(line);

    if(v[0][0] == '#') continue;

    map<string, adr_type>::iterator itr = variable_map.find(v[0]);

    if(itr == variable_map.end()) continue;

    switch (itr->second.type)
    { 
      case(INT):
      {
	int *ivar = static_cast<int* >(itr->second.var_adr);
	stringtovar(v[1], *ivar);
	break;
      }
      case(DOUBLE):
      {
	double *dvar = static_cast<double* >(itr->second.var_adr);
	stringtovar(v[1], *dvar);
	break;
      }
      case(STRING):
      {
	string *svar = static_cast<string*>(itr->second.var_adr);
	*svar = v[1];
	break;
      }
    }
  }

  std::cout << Params.A << endl;
  cout << Params.B << endl;
  cout << Params.C << endl;

  ifs.close();
}

bool fileexist(string filepath)
{
  bool flag = false;  
  FILE *fp;
  const char *filename = filepath.c_str();
  if((fp = fopen(filename, "r")) == NULL)
  {
    return false;
  }
  else 
  {
    fclose(fp);
    return true;
  }
}

bool mkdir(string str)
{
  const char *dirpath = str.c_str(); 

  mode_t  mode = S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR | S_IWGRP;
 
  if(mkdir(dirpath, mode) != 0)
  {
    perror("error:mkdir");
    return false;
  }

  return true;
}

template<class T> void shuffle(T ary[],int size)
{
  for(int i=0;i<size;i++)
  {
    int j = rand()%size;
    T t = ary[i];
    ary[i] = ary[j];
    ary[j] = t;
  }
}

#endif 

