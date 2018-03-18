#ifndef _UTILS_H
#define _UTILS_H
#include <iostream>
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

template<typename T1, typename T2>
class myPair
{
  public:
    myPair<T1, T2> *firstPair;
    myPair<T1, T2> *nextPair;
    myPair<T1, T2> *nextPair2;
    bool flag;
    T1 key;
    T2 value;
};


template<typename T1, typename T2>
class myMap
{
  public:

    int nsize;
    hash<T1> Hash;
    myPair<T1, T2> *Pair;
    vector< myPair<T1, T2>* > PairVec; 
    
    myMap()
    {
      nsize = 10000000;
      Pair = new myPair<T1, T2>[nsize];
      PairVec.reserve(nsize);

      for(int i = 0; i < nsize; ++i)
      {
	Pair[i].flag = false;
	Pair[i].firstPair = NULL; 
	Pair[i].nextPair = NULL;
      }
    }

    void add(T1 &key, T2 &value)
    {
      size_t sizet = Hash(key);
      sizet = sizet % nsize;

      for(myPair<T1, T2> *tmp = Pair[sizet].firstPair; tmp != NULL; tmp = tmp->nextPair)
      {
	if(value == tmp->value) return;
      }

      Pair[sizet].flag = true;
      myPair<T1, T2> *tmpPair = new myPair<T1, T2>;
      tmpPair->key = key;
      tmpPair->value = value;

      tmpPair->nextPair = Pair[sizet].firstPair;
      Pair[sizet].firstPair = tmpPair;

      if(PairVec.size() > 0) PairVec[PairVec.size() - 1]->nextPair2 = tmpPair;
      tmpPair->nextPair2 = NULL;
      PairVec.push_back(tmpPair);
    }

    myPair<T1, T2>* find(T1 &key)
    {
      size_t sizet = Hash(key) % nsize;

      myPair<T1, T2> *tmp = Pair[sizet].firstPair; 
      for(; tmp != NULL; tmp = tmp->nextPair)
      {
	if(key == tmp->key) break;
      }

      return tmp;
    }

    myPair<T1, T2>* begin()
    {
      return PairVec[0];
    }

    myPair<T1, T2>* end()
    {
      return NULL;
    }
    
    T2 operator [](T1 &key)
    {
      size_t sizet = Hash(key);
      sizet = sizet % nsize;

      for(myPair<T1, T2> *tmp = Pair[sizet].firstPair; tmp != NULL; tmp = tmp->nextPair)
      {
	if( key == tmp->key) break;
      }
      return Pair[sizet].firstPair->value;
    }

    int size()
    {
      return PairVec.size();
    }

    void clear()
    {

      for(int i = 0; i < PairVec.size(); ++i)
      {
	delete PairVec[i];
      }

      delete [] Pair;
    }
};

class UnionFinding
{
  public:
    vector<int> par;

    UnionFinding(){}
    UnionFinding(int nsize){par.reserve(nsize);for(int i = 0; i < nsize; ++i) par.push_back(i);}
    bool Union(int x, int y)
    {
      x = root(x);
      y = root(y);
       //y = par[y]; 
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

inline std::vector<string> split(const string &str)
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
inline void vartostring(string &str, var &val)
{
  stringstream ss;
  ss << val;
  str = ss.str();
}

template<typename var>
inline void stringtovar(string &str, var &val)
{
  istringstream (str) >> val;
}

template<typename var>
inline var sq(var x)
{
  return x * x;
}

template<typename var>
inline var cube(var &x)
{
  return x * x * x;
}

template<typename var>
inline void rnd1(var &val)
{

  val = static_cast<var> (rand()) / (static_cast<var>(RAND_MAX));
}

inline void init_rnd()
{
  srand((unsigned)time(NULL));
}

template<typename var>
inline void stringtimetovartime(string &str, var &val)
{
  var hour;
  var minite;
  vector<string> time = split(str);

  val = hour * 60 +  minite;
} 

template<typename var>
inline void vartimetostringtime(string &str, var &val)
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

inline void read_params(string &paramsfile)
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

  cout << Params.A << endl;
  cout << Params.B << endl;
  cout << Params.C << endl;

  ifs.close();
}

inline bool fileexist(string &filepath)
{
  //bool flag = false;  
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

inline bool mkdir(string &str)
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

template<class T> 
inline void shuffle(T ary[],int size, bool flag)
{
  for(int i=0;i<size;i++)
  {
    int j = rand()%size;
    //if(flag && (i == size - 1 || j == size - 1)) continue;
    if(flag && (i == 0 || j == 0)) continue;
    T t = ary[i];
    ary[i] = ary[j];
    ary[j] = t;
  }
}
#endif 

