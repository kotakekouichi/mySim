#include <iostream> 
#include "Vec.h"

void hilbertC(Vec *hilbertCurves, int &n, int s, double x, double y, double z, double dx, double dy, double dz, double dx2, 
    double dy2, double dz2, double dx3, double dy3, double dz3)
{
  if(s == 1)
  {
    hilbertCurves[n].set(x, y, z);
    if(n >= 16*16*16) getchar();
    ++n;
  }
  else
  {
    s/=2;
    if(dx<0) x-=s*dx;
    if(dy<0) y-=s*dy;
    if(dz<0) z-=s*dz;
    if(dx2<0) x-=s*dx2;
    if(dy2<0) y-=s*dy2;
    if(dz2<0) z-=s*dz2;
    if(dx3<0) x-=s*dx3;
    if(dy3<0) y-=s*dy3;
    if(dz3<0) z-=s*dz3;
    hilbertC(hilbertCurves, n, s, x, y, z, dx2, dy2, dz2, dx3, dy3, dz3, dx, dy, dz);
    hilbertC(hilbertCurves, n, s, x+s*dx, y+s*dy, z+s*dz, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2);
    hilbertC(hilbertCurves, n, s, x+s*dx+s*dx2, y+s*dy+s*dy2, z+s*dz+s*dz2, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2);
    hilbertC(hilbertCurves, n, s, x+s*dx2, y+s*dy2, z+s*dz2, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3);
    hilbertC(hilbertCurves, n, s, x+s*dx2+s*dx3, y+s*dy2+s*dy3, z+s*dz2+s*dz3, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3);
    hilbertC(hilbertCurves, n, s, x+s*dx+s*dx2+s*dx3, y+s*dy+s*dy2+s*dy3, z+s*dz+s*dz2+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2);
    hilbertC(hilbertCurves, n, s, x+s*dx+s*dx3, y+s*dy+s*dy3, z+s*dz+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2);
    hilbertC(hilbertCurves, n, s, x+s*dx3, y+s*dy3, z+s*dz3, dx2, dy2, dz2, -dx3, -dy3, -dz3, -dx, -dy, -dz);
  }

  return;
}

//m=0;
//hilbertC(256,0,0,0,1,0,0,0,1,0,0,0,1);}}}
