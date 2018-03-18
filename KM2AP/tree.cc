#include <iostream>
#include "mesh.h"
#include "particle.h"
#include "tree.h"
#include "def.h"
#include "utils.h"

using namespace std;

template<typename ctype>
int tree<ctype>::initialize(ctype *point)
{
  int rval = 0;

  particleClass *pfirst = &point[4], *p = NULL;
  pfirst->next = NULL;
  this->first = pfirst;
  for(int ipart = 5; ipart < num + 4; ++ipart)
  {
    p = &point[ipart];
    p->next = first;
    pfirst = p;
    this->first = pfirst;
  }

  return rval;
}

template<typename ctype>
int tree<ctype>::construct()
{
   int rval = 0;

   this->build_tree();

   return rval;
}

template<typename ctype>
int tree<ctype>::delaunayfinder(int imesh, ctype *Mesh, splitType &splittype)
{
  int rval = 0;
  int ipart = 0;

  ipart = walk_tree(imesh, Mesh);

  if(ipart < 0) 
  {
    rval = -1; 
    return rval;
  }

  int id = 0;
  //delaunayClass *tmpDelaunay = Mesh[ipart].delaunay[0]; 
  unordered_map<int, delaunayClass*>::iterator itr = Mesh[ipart].delaunaymap.begin();
  if(itr == Mesh[ipart].delaunaymap.end())
  {
    cout << ipart << " " << Mesh[ipart].delaunaymap.size() << endl; 
    getchar();
    return -1;
  }
  delaunayClass *tmpDelaunay = itr->second;
  
  while(!tmpDelaunay->flag){tmpDelaunay=Mesh[ipart].delaunay[id];++id;}
  

  if(getdim() == 2)  rval = tmpDelaunay->crossdeterminant(Mesh[imesh].getPosx(), Mesh[imesh].getPosy());
  else if(getdim() == 3) 
  {
    splittype = tmpDelaunay->getSplitType(Mesh[imesh].getPosx(), Mesh[imesh].getPosy(), Mesh[imesh].getPosz());  
     
    if(splittype != _nosplit)  
      rval = 0;
    else rval = 1;
  }
  
  if(rval == 0) 
  {
    rval = tmpDelaunay->indd;
    return rval;
  }

  int iterator = 0;
  meshClass *mesh = &Mesh[imesh];

  while(1)
  {
    double xg = 0.0, yg = 0.0, zg = 0.0, s = 0.0, t = 0.0, u = 0.0;

    if(this->getdim() == 2)
    {
      xg = (tmpDelaunay->getPosx(0) + tmpDelaunay->getPosx(1) + tmpDelaunay->getPosx(2)) / 3.0;
      yg = (tmpDelaunay->getPosy(0) + tmpDelaunay->getPosy(1) + tmpDelaunay->getPosy(2)) / 3.0;
    } else if( this->getdim() == 3)
    {
      xg = (tmpDelaunay->getPosx(0) + tmpDelaunay->getPosx(1) + tmpDelaunay->getPosx(2) + tmpDelaunay->getPosx(3)) / 4.0;
      yg = (tmpDelaunay->getPosy(0) + tmpDelaunay->getPosy(1) + tmpDelaunay->getPosy(2) + tmpDelaunay->getPosy(3)) / 4.0;
      zg = (tmpDelaunay->getPosz(0) + tmpDelaunay->getPosz(1) + tmpDelaunay->getPosz(2) + tmpDelaunay->getPosz(3)) / 4.0;
    }

    int id[3];
    if(this->getdim() == 2)
    {

      for(int i = 0; i < 3; ++i)
      {
	for(int j = i + 1; j < 3; ++j)
	{
	  id[0] = tmpDelaunay->indm[i]; 
	  id[1] = tmpDelaunay->indm[j];
	  double p[2], a[2], b[2], det;
	  double a_inv[2], b_inv[2];

	  p[0] = Mesh[imesh].getPosx() - xg;
	  p[1] = Mesh[imesh].getPosy() - yg;

	  a[0] = tmpDelaunay->getPosx(i) - xg;
	  b[0] = tmpDelaunay->getPosy(i) - yg;

	  a[1] = tmpDelaunay->getPosx(j) - xg;
	  b[1] = tmpDelaunay->getPosy(j) - yg;

	  det = DET2D(a[0], a[1], b[0], b[1]);
	  a_inv[0] = b[1];
	  b_inv[0] = -a[1];
	  a_inv[1] = -b[0];
	  b_inv[1] = a[0];

	  s = (a_inv[0] * p[0] + b_inv[0] * p[1]) / det;
	  t = (a_inv[1] * p[0] + b_inv[1] * p[1]) / det;

	  if(s >= 0.0 && t >= 0.0) goto NEXT; 
	}
      }
    } else if (this->getdim() == 3)
    {
    
      for(int i = 0; i < 4; ++i)
      {
        for(int j = i + 1; j < 4; ++j)
	{
	  for(int k = j + 1; k < 4; ++k)
	  {
	     id[0] = tmpDelaunay->indm[i]; id[1] =  tmpDelaunay->indm[j]; id[2] = tmpDelaunay->indm[k]; 
	     double p[3],a[3], b[3], c[3];
	     double a_inv[3], b_inv[3], c_inv[3];
	     p[0] = Mesh[imesh].getPosx() - xg;
	     p[1] = Mesh[imesh].getPosy() - yg;
	     p[2] = Mesh[imesh].getPosz() - zg;
	     a[0] = tmpDelaunay->getPosx(i) - xg;
	     b[0] = tmpDelaunay->getPosy(i) - yg;
	     c[0] = tmpDelaunay->getPosz(i) - zg;
	     a[1] = tmpDelaunay->getPosx(j) - xg;
	     b[1] = tmpDelaunay->getPosy(j) - yg;
	     c[1] = tmpDelaunay->getPosz(j) - zg;
	     a[2] = tmpDelaunay->getPosx(k) - xg;
	     b[2] = tmpDelaunay->getPosy(k) - yg;
	     c[2] = tmpDelaunay->getPosz(k) - zg;

	     double det = DET3D(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
	     a_inv[0] = (b[1]*c[2] - b[2]*c[1]); a_inv[1] = (b[0]*c[2]-b[2]*c[0]) * (-1.0) ; a_inv[2] = (b[0]*c[1] - b[1]*c[0]);
	     b_inv[0] = -(a[1]*c[2] - a[2]*c[1]); b_inv[1] = (a[0]*c[2]-a[2]*c[0])         ; b_inv[2] = -(a[0]*c[1] - a[1]*c[0]);
	     c_inv[0] = (a[1]*b[2] - a[2]*b[1]); c_inv[1] = (a[0]*b[2]-a[2]*b[0]) * (-1.0) ; c_inv[2] = (a[0]*b[1] - a[1]*b[0]);

	     s = (a_inv[0] * p[0] + b_inv[0] * p[1] + c_inv[0] * p[2]) / det;
	     t = (a_inv[1] * p[0] + b_inv[1] * p[1] + c_inv[1] * p[2]) / det;
	     u = (a_inv[2] * p[0] + b_inv[2] * p[1] + c_inv[2] * p[2]) / det;
	     if(s >= 0 && t >= 0 && u >= 0) goto NEXT;
	  }
	} 
      }
    
    }
NEXT:;

     for(int inbr = 0; inbr < tmpDelaunay->nNeighbor; ++inbr)
     {
       if(getdim() == 2)
       {
	 if(tmpDelaunay->neighbors[inbr]->has_mesh(id[0]) &&
	     tmpDelaunay->neighbors[inbr]->has_mesh(id[1]) ) 
	 {
	   tmpDelaunay = tmpDelaunay->neighbors[inbr];
	   break;
	 }
       } else if(getdim() == 3)
       {
         if(tmpDelaunay->neighbors[inbr]->has_mesh(id[0]) && tmpDelaunay->neighbors[inbr]->has_mesh(id[1]) 
	    && tmpDelaunay->neighbors[inbr]->has_mesh(id[2]))
	 {
	   tmpDelaunay = tmpDelaunay->neighbors[inbr];
	   break;
	 } 
       }
    }

    if(getdim() == 2) rval = tmpDelaunay->crossdeterminant(mesh->getPosx(), mesh->getPosy());
    else if(getdim() == 3)
    {
      splittype = tmpDelaunay->getSplitType(mesh->getPosx(), mesh->getPosy(), mesh->getPosz());
      rval = splittype == _nosplit ? 1 : 0;  
    } 

    if(rval == 0)
      break;
    ++iterator; 

    if(iterator == 1000) return -1; 
  }

  rval = tmpDelaunay->indd;

  return rval;
}

template<typename ctype> 
int tree<ctype>::walk_tree(int imesh, meshClass *Mesh)
{
  int rval = 0;
  int i = 0, j = 0, k = 0;
  const int Ncount = 10;

  double xp = Mesh[imesh].getPosx();
  double yp = Mesh[imesh].getPosy();
  double zp = Mesh[imesh].getPosz();

  tree<meshClass> *Tree = this;
  meshClass *part = NULL;

  while(1)
  {
    i = j = k = 0;
    Tree->num += 1; 
    if(Tree->first != NULL && Tree->first->index != imesh) part = Tree->first;

    Mesh[imesh].next = Tree->first;
    Tree->first = &Mesh[imesh];

    if(Tree->num > Ncount)
    {
      if(Tree->nextlevel == NULL)
      {
	Tree->build_tree(); 
      }
    }

    if( xp > Tree->getxc() ) i = 1;
    if( yp > Tree->getyc() ) j = 1;
    if( zp > Tree->getzc() ) k = 1;

    if(Tree->nextlevel == NULL)
      break;
    int index = getdim() == 2 ? i + 2* j : i + 2 * (j + 2 *k );
    Tree = &(Tree->nextlevel[index]);
  }

  if(part != NULL && part->index != imesh)
    rval = part->index; 
  else
   rval = -1; 

  return rval;
}

template<typename ctype>
int tree<ctype>::build_tree()
{
  int rval = 0;
  int i = 0, j = 0, k = 0, ichild = 0;

  if(this->getdim() == 2) this->nextlevel = new tree[4];  
  else if (this->getdim() == 3) this->nextlevel = new tree[8];

  if( this->getdim() == 2) 
  {
    for(i = 0; i < 2; ++i)
    {
      for(j = 0; j < 2; ++j)
      {
	ichild = i + 2 * j;
	this->nextlevel[ichild].num = 0;
	this->nextlevel[ichild].setdx(this->dx / 2.0);
	this->nextlevel[ichild].setxc(this->xc - this->dx / 2.0 + i * this->dx);
	this->nextlevel[ichild].setyc(this->yc - this->dx / 2.0 + j * this->dx);      
	this->nextlevel[ichild].nextlevel = NULL;
	this->nextlevel[ichild].first = NULL;
	this->nextlevel[ichild].setdim(this->getdim());
	this->nextlevel[ichild].nth = this->nth;
      }
    }
  } else if(this->getdim() == 3)
  {

    for(i = 0; i < 2; ++i)
    {
      for(j = 0; j < 2; ++j)
      {
	for(k = 0; k < 2; ++k)
	{
	  ichild = i + 2 * (j + 2 *k);
	  this->nextlevel[ichild].num = 0;
	  this->nextlevel[ichild].setdx(this->dx / 2.0);
	  this->nextlevel[ichild].setxc(this->xc - this->dx / 2.0 + i * this->dx);
	  this->nextlevel[ichild].setyc(this->yc - this->dx / 2.0 + j * this->dx);
	  this->nextlevel[ichild].setzc(this->zc - this->dx / 2.0 + k * this->dx);
	  this->nextlevel[ichild].nextlevel = NULL;
	  this->nextlevel[ichild].first = NULL;
	  this->nextlevel[ichild].setdim(this->getdim());
	  this->nextlevel[ichild].nth = this->nth;
	}
      }
    }

  }

  for(ctype *pfirst = first; pfirst != NULL;)
  {
    ctype *pnext = pfirst->next;
    double xp = pfirst->getPosx();
    double yp = pfirst->getPosy();
    double zp = pfirst->getPosz();

    i = j = k = 0; 
    if(xp > this->xc) i = 1;   
    if(yp > this->yc) j = 1;
    if(zp > this->zc) k = 1;
    ichild = getdim() == 2 ? i + j * 2 : i + 2 * (j + 2*k);

    this->nextlevel[ichild].num += 1;
    pfirst->next = this->nextlevel[ichild].first;  
    this->nextlevel[ichild].first = pfirst; 

    this->nextlevel[ichild].mass += pfirst->getMass();
    this->nextlevel[ichild].xg += pfirst->getMass() * pfirst->getPosx();
    this->nextlevel[ichild].yg += pfirst->getMass() * pfirst->getPosy();
    this->nextlevel[ichild].zg += pfirst->getMass() * pfirst->getPosz();

    pfirst = pnext;
  } 

  if(this->getdim() == 2)
  {
     for(int ichild = 0; ichild < 4; ++ichild)
     {
       this->nextlevel[ichild].xg /= this->nextlevel[ichild].mass;
       this->nextlevel[ichild].yg /= this->nextlevel[ichild].mass;

     }
  } else if(this->getdim() == 3)
  {
    for(int ichild = 0; ichild < 8; ++ichild)
    {
      this->nextlevel[ichild].xg /= this->nextlevel[ichild].mass;
      this->nextlevel[ichild].yg /= this->nextlevel[ichild].mass;
      this->nextlevel[ichild].zg /= this->nextlevel[ichild].mass;
    }
  }
  

  if(this->getdim() == 2)
  {
    for(int i = 0; i < 2; ++i)
    {
      for(int j = 0; j < 2; ++j)
      {
	ichild = i + j * 2;    

	if(this->nextlevel[ichild].num > nth)
	  rval = this->nextlevel[ichild].build_tree();
      }
    }
  } else if(this->getdim() == 3)
  {
    for(int i = 0; i < 2; ++i)
    {
      for(int j = 0; j < 2; ++j)
      {
	for(int k = 0; k < 2; ++k)
	{
	  ichild = i + 2 * ( j + 2 * k);    

	  if(this->nextlevel[ichild].num > nth)
	    rval = this->nextlevel[ichild].build_tree();

	}
      }
    }

  }

  return rval;

}

template<typename ctype>
int tree<ctype>::gravity_force(particleClass *part)
{
  int rval = 0;
  const double theta = 0.7, softning = 0.000001;
  const double theta2 = theta * theta, softning2 = sq(softning);
  const double xp = part->getPosx(), yp = part->getPosy(), zp = part->getPosz();
  double length2 = 0.0, dx2 = 0.0, M = 0.0;

  if( this->getdim() == 2) 
  {
  } else if(this->getdim() == 3)
  {

    for(int ichild = 0; ichild < 8; ++ichild)
    {
      M  = this->nextlevel[ichild].mass;
      dx = sq(this->nextlevel[ichild].dx);
      length2 =  sq(this->nextlevel[ichild].xg - xp) + sq(this->nextlevel[ichild].yg - yp) + sq(this->nextlevel[ichild].zg - zp);

	if( theta2 < dx2 / length2)
	{
	  this->nextlevel[ichild].gravity_force(part);
	}else 
	{
	  part->setFx(part->getFx() + M / pow(length2 + softning2, 1.5) * (xp - this->xg));
	  part->setFy(part->getFy() + M / pow(length2 + softning2, 1.5) * (yp - this->yg));
	  part->setFz(part->getFz() + M / pow(length2 + softning2, 1.5) * (zp - this->zg));
	}   
    }

  }

  return rval = 0;
}

  template<typename ctype>
int tree<ctype>::gravity_shortrange_force(particleClass *part)
{
  int rval = 0;
  const double theta = 0.7, softning = 0.000001;
  const double theta2 = theta * theta, softning2 = sq(softning);
  const double xp = part->getPosx(), yp = part->getPosy(), zp = part->getPosz();
  double length = 0.0, length2 = 0.0, dx2 = 0.0, M = 0.0;
  double fx =0.0, fy = 0.0, fz = 0.0;
  double fac = 0.0;

  if( this->getdim() == 2) 
  {
  } else if(this->getdim() == 3)
  {

    for(int ichild = 0; ichild < 8; ++ichild)
    {
      M  = this->nextlevel[ichild].mass;
      dx = sq(this->nextlevel[ichild].dx);
      length2 =  sq(this->nextlevel[ichild].xg - xp) + sq(this->nextlevel[ichild].yg - yp) + sq(this->nextlevel[ichild].zg - zp);

	if( theta2 < dx2 / length2)
	{
	  this->nextlevel[ichild].gravity_force(part);
	}else 
	{
	  fx = part->getFx();
	  fy = part->getFy();
	  fz = part->getFz();

	  length = sqrt(length2);

	  fac = erfc(length / 2.0) + length / sqrt(M_PI) * exp(-length2 / 4.0); 
	  part->setFx(fx + M / pow(length2 + softning2, 1.5) * fac * (xp - this->xg));
	  part->setFy(fy + M / pow(length2 + softning2, 1.5) * fac * (yp - this->yg));
	  part->setFz(fz + M / pow(length2 + softning2, 1.5) * fac * (zp - this->zg));

	}   
    }

  }

  return rval = 0;
}


  template<typename ctype>
int tree<ctype>::free()
{
  int rval = 0;
  const int nchild = this->getdim() ? 4 : 8;

  if(this->nextlevel != NULL)
  {
    for(int i = 0; i < nchild; ++i)
      rval = this->nextlevel[i].free(); 
    delete [] this->nextlevel;
  }
  return 0;
}

template int tree<particleClass>::gravity_shortrange_force(particleClass *part);
template int tree<particleClass>::gravity_force(particleClass* part);
template int tree<particleClass>::construct();
template int tree<particleClass>::initialize(particleClass *point);
template int tree<meshClass>::delaunayfinder(int imesh, meshClass *Mesh, splitType &splittype); 
template int tree<meshClass>::walk_tree(int imes, meshClass *Mesh);
template int tree<meshClass>::free();
template int tree<particleClass>::free();
