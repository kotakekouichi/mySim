#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <algorithm>

using namespace std;
class store
{
  public:
    int ordernum;
    int volume;
    int stime;
    int etime;
    int ctime;
    double ddepo;
    vector<double> dij;
    string name;
};

class track
{

  public:
    int index;
    int ctime;
    string name;
    string starttime;
    string endtime;
    vector<string> storename;
    vector<int> time;

    bool operator < (const track &right) const{
      return (ctime == right.ctime) ? storename.size() < right.storename.size() : ctime < right.ctime;
    }
};

inline int StringToInt(const string str);
inline double StringToDouble(const string str);
vector<string> split(const string &s, char delim[]);
void Function();
void ReadFile(string day_week, string day_week_JP, vector<store> &StoreVec);
void ReadFile_Track(string day_week, int istore, vector<track> &Track);
string GetFileName(string day_week, int num);
void Calcu_ctime(string day_week, int istore, vector<track> &Track);
void Analysis(int istore, string day_week, vector<track> &Track, vector<store> &StoreVec);

int main()
{

  Function();

  return 0;
}

void Function()
{
  string day_week[] = {"Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"};
  string day_week_JP[] = {"月","火","水","木","金","土","日"};

  for(int M = 0; M < 7;++M)
  {
    vector<store>StoreVec;
    vector<track>Track;

    ReadFile(day_week[M], day_week_JP[M], StoreVec);
    int istore = -1;
    ReadFile_Track(day_week_JP[M], istore, Track);
    Calcu_ctime(day_week_JP[M], istore, Track);
    Analysis(istore, day_week[M], Track, StoreVec);

    std::cout << day_week[M] << std::endl;

  }

}

void ReadFile(string day_week, string day_week_JP, vector<store> &StoreVec)
{
  int istore = 0, jstore = 0;
  char delim[] = {"\t :,"};
  string filename, buf;
  filename = "./input/yaoko/" + day_week + "/order.dat";
  ifstream ifs1(filename); 
  filename = "./input/yaoko/" + day_week + "/siteMaster.dat";
  ifstream ifs2(filename); 
  filename = "./input/yaoko/" + day_week + "/transportTime.dat";
  ifstream ifs3(filename); 

  // ORDER FILE
  getline(ifs1, buf);
  while(getline(ifs1, buf))
  {
    store vec;
    vector<string> elems = split(buf, delim);
    vec.ordernum = StringToInt(elems[1]);
    vec.volume = StringToInt(elems[2]);
    vec.stime = StringToInt(elems[3]) * 60 + StringToInt(elems[4]);
    vec.etime = StringToInt(elems[5]) * 60 + StringToInt(elems[6]);
    vec.ctime = 1000;//large time
    StoreVec.push_back(vec);
  }

  // SITEMASTER FILE
  
  getline(ifs2, buf);
  getline(ifs2, buf);
  while(getline(ifs2, buf))
  {
    vector<string> elems = split(buf, delim);
    StoreVec[istore].name = elems[1];
    ++istore;
  }

  for(int i = (int) StoreVec.size() - 1; i >= 0; --i)
  {
    int j = 0;
    while(StoreVec[i].ordernum != StoreVec[j].ordernum) ++j;
    StoreVec[i].name = StoreVec[j].name;
  }

  // TRANSPORTTIME FILE

  bool nextskip = false;
  getline(ifs3, buf);
  while(getline(ifs3, buf))
  {
    vector<string> elems = split(buf, delim);

    if(nextskip)
    {
      nextskip = false;
      continue;
    }

    if(StringToDouble(elems[4]) > 0.0)
    {
      nextskip = true;
    }

    if(elems[0] == "99999")
    {
      if(elems[1] == "99999") break;

      StoreVec[jstore].ddepo = StringToDouble(elems[3]) + StringToDouble(elems[4]) / 3.0;
      ++jstore;
      continue;
    }

    if(to_string(StoreVec[istore].ordernum) != elems[0]) ++istore;

    StoreVec[istore].dij.push_back(StringToDouble(elems[3]) + StringToDouble(elems[4]) / 3.0);

  }
}



vector<string>split(const string &s, char delim[])
{
  vector<string> v;
  //char delim[] = {'=', ' ', '\t'};
  string item;
  int i=0;
  while(i != s.size()){
    char ch = s[i];
    if (ch == delim[0] || ch == delim[1]
	|| ch == delim[2]) {
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

void ReadFile_Track(string day_week, int istore, vector<track> &Track)
{
  string buf;
  int it = 0;
  char delim[] = ",???";
  ifstream ifs("./input/yaoko/route/def/" + GetFileName(day_week, 0));
  //ifstream ifs("./input/yaoko/route/" + to_string(istore) + ...);
  getline(ifs, buf);
  while(getline(ifs, buf))
  {
    vector<string> elems = split(buf, delim);
    if(Track.size() == 0)
    {
      track Ti;
      Ti.name = elems[0];
      Ti.starttime = elems[12];
      Track.push_back(Ti);
    }
    else if(Track[(int)Track.size() - 1].name != elems[0])
    {
      track Ti;
      Ti.name = elems[0];
      Ti.starttime = elems[12];
      Track.push_back(Ti);
    }

    it = (int) Track.size() - 1;
    Track[it].index = it;
    Track[it].endtime = elems[13];
    Track[it].storename.push_back(elems[5]);

    char delim2[] = ":???";
    vector<string> strTime = split(elems[15], delim2);
    int Time = 60 * StringToInt(strTime[0]) + StringToInt(strTime[1]);
    Track[it].time.push_back(Time);
  }
}

string GetFileName(string day_week, int num)
{
  string str;
  if(num == 0)
    str = ".dat";
  else if(num == 1)
    str = ".dat";  
  return str;
}

void Calcu_ctime(string day_week, int istore, vector<track> &Track)
{
  int itrack = 0;
  char delim[] = "\t,??";
  string buf;
  ifstream ifs("./input/yaoko/route/def/" + GetFileName(day_week, 1));

  getline(ifs, buf);
  while(getline(ifs, buf))
  {
    vector<string> elems = split(buf, delim);
    Track[itrack].ctime = StringToInt(elems[5]);
    ++itrack;
  }
}

void Analysis(int istore, string day_week, vector<track>&Track, vector<store> &StoreVec)
{
  double vel = 20.0, eps = 17.0, d = 0.0;

  for(int i = 0;i < (int) StoreVec.size(); ++i)
  {
    double pot = 0.0;
    for(int j = 0; j < (int)StoreVec.size(); ++j)
    {
      if(i == j) continue;
      if(StoreVec[j].etime - StoreVec[i].etime < 0.0) continue;

      int k = StoreVec[i].ordernum - 1; 
      int l = StoreVec[j].ordernum - 1;
      double dti = (double) (StoreVec[i].etime + StoreVec[i].stime) / 2.0;
      double dtj = (double) (StoreVec[j].etime + StoreVec[j].stime) / 2.0;
      double dt = dtj - dti + 0.1;

      if(0 < StoreVec[i].volume + StoreVec[j].volume && StoreVec[i].volume + StoreVec[j].volume <= eps)
      {
        d = StoreVec[k].dij[l];	
	if(vel * abs(StoreVec[i].etime - StoreVec[j].etime) / 60.0 < d
	    && vel * (abs(StoreVec[i].etime - StoreVec[j].etime) + 120.0) / 60.0 > d
	    && d / vel * 60.0 + StoreVec[i].etime > StoreVec[j].stime)
	    pot += 1.0;
      }
      else 
      {
	int eps1 = StoreVec[i].volume > eps ? StoreVec[i].volume - eps : StoreVec[i].volume;
	int eps2 = StoreVec[j].volume > eps ? StoreVec[j].volume - eps : StoreVec[j].volume;

	if(eps1 + eps2 < eps)
	{
	  d = StoreVec[k].dij[l];
	  if(vel * abs(StoreVec[i].etime - StoreVec[j].etime) / 60.0 < d
	      && vel * (abs(StoreVec[i].etime - StoreVec[j].etime) + 120.0) / 60.0 > d
	      && d / vel * 60.0 + StoreVec[i].etime > StoreVec[j].stime)
	    pot += 1.0;
	}
	else 
	{
	  d = StoreVec[k].ddepo + StoreVec[l].ddepo;
	  if(vel * abs(StoreVec[i].etime - StoreVec[j].etime) / 60.0 < d
	      && vel * (abs(StoreVec[i].etime - StoreVec[j].etime) + 120.0) / 60.0 > d
	      && d / vel * 60.0 + StoreVec[i].etime > StoreVec[j].stime)
	    pot += 1.0;
	}
      }

    }// for j
  }//for i 

}
