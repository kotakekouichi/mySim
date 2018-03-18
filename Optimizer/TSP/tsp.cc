#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#define N 10

class _Node;
class _edge;

/*--- function ---*/
void TSP();
void Initialize(_Node *Node, _edge *Edge);
void Output(_Node *Node);
void twoopt(_Node *Node);
double distance(_Node *Node);
inline double SQ(double x){ return (x * x); }
inline double rnd(){ return ((double) rand() / (double) RAND_MAX); }
/*----------------*/
/*--- class ---*/
class _Node
{
  public:
    double x;
    double y;
};

class _edge
{
  public:
    _Node *Node[2];
    inline double dist(){return (sqrt(SQ(Node[0]->x - Node[1]->x) + SQ(Node[0]->y - Node[1]->y)));}
};

/*-------------*/

int main()
{

  TSP();

  return 0;
}

void TSP()
{

  _Node *Node = new _Node[N];
  _edge *Edge = new _edge[N];

  Initialize(Node, Edge);

  twoopt(Node);

  Output(Node);

  delete [] Node; 
}

void Initialize(_Node *Node, _edge *Edge)
{
  for(int i = 0;i < N; ++i)
  {
    Node[i].x = rnd();
    Node[i].y = rnd();
  }

  for(int i = 0; i < N; ++i)
  {
    int j = i + 1 != N ? i + 1 : 0;
    int k = i - 1 != -1 ? i - 1 : N - 1;
    Edge[i].Node[0] = &Node[j];
    Edge[i].Node[1] = &Node[k];
  }

}

void Output(_Node *Node)
{
  FILE *fp;

  fp = fopen("position.dat", "w");

  for(int i = 0; i < N; ++i)
    fprintf(fp, "%e	%e\n",Node[i].x, Node[i].y);

  fclose(fp);

  fp = fopen("tsp.gnu", "w");

  fprintf(fp, "plot '-' w l\n");
  for(int i = 0; i < N; ++i)
    fprintf(fp, "%e %e\n", Node[i].x, Node[i].y);
  //fprintf(fp, "%e %e\n", Node[0].x, Node[0].y);
  fprintf(fp,"e\n");

  fclose(fp);
}

void twoopt(_Node *Node)
{
#define Dist(A, B) (sqrt( SQ(Node[A].x - Node[B].x)  + SQ(Node[A].y - Node[B].y)))

  int A, B, C, D;
  int inode, jnode;
  double diff = 0.0;

  for(A = 0; A  < N; ++A)
  {

    for(C = B + 1; C < N ;++C)
    {
      double AB = Dist(A, B);
      double CD = Dist(C, D);
      double AC = Dist(A, C);
      double BD = Dist(B, D);
    }
  }
}

double distance(_Node *Node)
{
  int i, j;
  double dist = 0.0;

  for(i = 0; i < N; ++i)
  {
    int j = i != N - 1 ? i + 1 : 0;
    dist = sqrt(SQ(Node[j].x - Node[i].x) + SQ(Node[j].y - Node[i].y));
  }
 
  return dist;

} 
