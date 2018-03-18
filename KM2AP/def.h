#ifndef _DEF_H
#define _DEF_H

#define BOXSIZE 1.0
#define LEFTSIDE 0.0
#define RIGHTSIDE 1.0
#define DX 0.12
#define DET2D(a,b,c,d) ((a)*(d) - (b)*(c))
#define DET3D(a,b,c,d,e,f,g,h,i) ((a)*(e)*(i) + (b)*(f)*(g) + (c)*(d)*(h) - (c)*(e)*(g)- (a)*(f)*(h) - (b)*(d)*(i))
#define IDX(a, b) ((a) + DIM * (b))
#define MINDOUBLE 1.0e-16

#endif
