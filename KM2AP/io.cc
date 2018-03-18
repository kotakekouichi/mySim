#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "utils.h"
#include "params.h"
#include "particle.h"
#include "mesh.h"

using namespace std;

void read_params(string &paramsfile, paramsClass &params)
{
  ifstream ifs(paramsfile);
  string line;
  map<string, adr_type> variable_map;
  adr_type at;

  at.type = INT;
  at.var_adr = &params.periodic;
  variable_map["periodic"] = at;

  at.type = INT;
  at.var_adr = &params.selfgravity;
  variable_map["selfgravity"] = at;

  at.type = INT;
  at.var_adr = &params.dim;
  variable_map["dim"] = at;

  at.type = INT;
  at.var_adr = &params.Ndm;
  variable_map["Ndm"] = at;

  at.type = INT;
  at.var_adr = &params.Ngas;
  variable_map["Ngas"] = at;

  at.type = INT;
  at.var_adr = &params.NPMGRID;
  variable_map["NPMGRID"] = at;

  at.type = DOUBLE;
  at.var_adr = &params.OmegaLambda;
  variable_map["OmegaLambda"] = at;

  at.type = DOUBLE;
  at.var_adr = &params.OmegaMatter;
  variable_map["OmegaMatter"] = at;

  at.type = DOUBLE;
  at.var_adr = &params.OmegaCDM;
  variable_map["OmegaCDM"] = at;

  at.type = DOUBLE;
  at.var_adr = &params.OmegaBaryon;
  variable_map["OmegaBaryon"] = at;

  at.type = DOUBLE;
  at.var_adr = &params.Hubble0;
  variable_map["Hubble0"] = at;

  at.type = STRING;
  at.var_adr = &params.dmfile;
  variable_map["dmfile"] = at;

  at.type = STRING;
  at.var_adr = &params.gasfile;
  variable_map["gasfile"] = at;

  at.type = INT;
  at.var_adr = &params.cosmological;
  variable_map["cosmological"] = at;

  at.type = INT;
  at.var_adr = &params.Ngrid;
  variable_map["Ngrid"] = at;


  while(getline(ifs, line))
  {
    if(line == "" ) continue;

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

  cout << "end read params file " <<endl;
  return;
}

void outputMesh(meshClass *Mesh, paramsClass &params)
{
  ofstream meshfile("mesh.dat");
  int buf = params.dim == 2 ? 3 : 4;
  const int Nmesh = params.Ngas + buf;

  for(int i = buf; i < Nmesh ; ++i)
  {
    meshfile << Mesh[i].getPosx() << " "  << Mesh[i].getPosy() << " " << Mesh[i].getPosz() << endl; 
  }
  cout << "output => mesh.dat" << endl;
  meshfile.close();
}

void outputDelaunay(delaunayClass *Delaunay, paramsClass &params)
{

  if(params.dim == 2)
  {
    //for(int imesh = 0; imesh < Delaunay.size(); ++imesh )
    {

      //string str = "";
      //vartostring(str, imesh);
      //string gnufile = str + ".gnu";
      ofstream delaunaygnufile("delaunay.gnu");
      ofstream delaunayfile("output/delaunay.dat");
      //ofstream delaunaygnufile(gnufile);
      delaunaygnufile << "#!/bin/gnuplot" << endl;
      delaunaygnufile << "set size square" << endl;
      delaunaygnufile << "set xrange[-1.0:2.0]" << endl;
      delaunaygnufile << "set yrange[-1.0:2.0]" << endl;
      delaunaygnufile << "set parametric" << endl;
      delaunaygnufile << "set trange[0:1]" << endl;
      delaunaygnufile << "set trange[0.0:2*pi] " << endl;
      delaunaygnufile << "const1 = 0" << endl;
      delaunaygnufile << "const2 = 1" << endl;
      delaunaygnufile << "plot const1, t/2.0/pi ls 1 , const2,t/2.0/pi ls 1,t/2.0/pi, const1 ls 1,t/2.0/pi ,const2 ls 1, 0.707107*cos(t)+0.5, 0.707107*sin(t) +0.5 ls 1, -0.914214, -0.724745 ls 1, 1.914214, -0.724745 ls 1, 0.500000, 1.914214  ls 1, 'mesh.dat' u 1:2, '-' w  l" << endl;
      //<< Delaunay[imesh].getxc() << " +" <<  Delaunay[imesh].getRadius() << " * cos(t),"
      //<< Delaunay[imesh].getyc() << " +" <<  Delaunay[imesh].getRadius() << " * sin(t) "<< endl;
      delaunaygnufile << "-0.724745       -0.207107" << endl;
      delaunaygnufile << "1.724745        -0.207107" << endl;
      delaunaygnufile << "0.500000        1.914214" << endl;
      delaunaygnufile << "-0.724745       -0.20710" << endl;

      for(int i = 0; i < params.Ndelaunay; ++i)
      {
	delaunaygnufile << endl;
	delaunaygnufile << Delaunay[i].getPosx(0) << "	" << Delaunay[i].getPosy(0) << endl;
	delaunaygnufile << Delaunay[i].getPosx(1) << "	" << Delaunay[i].getPosy(1) << endl; 
	delaunaygnufile << Delaunay[i].getPosx(2) << "	" << Delaunay[i].getPosy(2) << endl;
	delaunaygnufile << Delaunay[i].getPosx(0) << "	" << Delaunay[i].getPosy(0) << endl;

	delaunayfile << Delaunay[i].getPosx(0) << "\t" << Delaunay[i].getPosy(0) << endl;
	delaunayfile << Delaunay[i].getPosx(1) << "\t" << Delaunay[i].getPosy(1) << endl; 
	delaunayfile << Delaunay[i].getPosx(2) << "\t" << Delaunay[i].getPosy(2) << endl;
	delaunayfile << Delaunay[i].getPosx(0) << "\t" << Delaunay[i].getPosy(0) << endl;

	delaunayfile << endl;
      }

      delaunaygnufile << "e " << endl;
      delaunaygnufile.close();

      //cout << imesh << " " << Delaunay.size() << " " << Delaunay[imesh].neighbors.size() <<  "  "
      //<< Delaunay[imesh].getxc() << endl;
      for(int imesh = 0; imesh < params.Ndelaunay; ++imesh)
	Delaunay[imesh].determinant2d(0.0,0.0);
      //getchar();
    }
  } else if (params.dim == 3)
  {

    //const int id = 3;
    for(int id = 0; id < params.Ndelaunay; ++id)
    {

      string str = "";
      vartostring(str, id);

      string filename = "output/delaunay" + str + ".gnu";
      //ofstream delaunaygnufile("delaunay.gnu");
      ofstream delaunaygnufile(filename);
      ofstream delaunayfile("output/delaunay.dat");

#if 0 
      delaunaygnufile << "#!/bin/gnuplot" << endl;
      delaunaygnufile << "set size square" << endl;
      delaunaygnufile << "splot '-' w l" << endl;
      for(int i = 0; i < params.Ndelaunay; ++i)
      {
	//delaunaygnufile << Tetrahedron[i]->xc << "\t" << Tetrahedron[i]->yc << "\t" << Tetrahedron[i]->zc << endl; 
	delaunaygnufile <<  Delaunay[i].getPosx(0) << "\t" << Delaunay[i].getPosy(0) << "\t" << Delaunay[i].getPosz(0) << endl;
	delaunaygnufile <<  Delaunay[i].getPosx(1) << "\t" << Delaunay[i].getPosy(1) << "\t" << Delaunay[i].getPosz(1) << endl;
	delaunaygnufile <<  Delaunay[i].getPosx(2) << "\t" << Delaunay[i].getPosy(2) << "\t" << Delaunay[i].getPosz(2) << endl;
	delaunaygnufile <<  Delaunay[i].getPosx(3) << "\t" << Delaunay[i].getPosy(3) << "\t" << Delaunay[i].getPosz(3) << endl;
	delaunaygnufile <<  Delaunay[i].getPosx(0) << "\t" << Delaunay[i].getPosy(0) << "\t" << Delaunay[i].getPosz(0) << endl;

	delaunaygnufile << endl;
      } 
#endif
      delaunayClass *p = &Delaunay[id];
      for(int i= 0; i < 4; ++i)
      {
	for(int j = i + 1; j < 4; ++j)
	{
	  delaunaygnufile << p->getPosx(i) << "\t" << p->getPosy(i) << "\t" << p->getPosz(i) << endl;
	  delaunaygnufile << p->getPosx(j) << "\t" << p->getPosy(j) << "\t" << p->getPosz(j) << endl;
	  delaunaygnufile << endl;
	}
      }

      delaunaygnufile << "e" << endl;

      delaunaygnufile.close();
      delaunayfile.close();
    }

  }

  cout << "output => delaunay.gnu, output/delaunay.dat" << endl;
}

void outputVoronoi(meshClass *Mesh, paramsClass &params)
{
  ofstream voronoignufile("voronoi.gnu");
  ofstream voronoifile("output/voronoi.dat");

  cout << "output voronoi " << endl;
  switch (params.dim)
  {
    case(2):
      {

	const int Nmesh = params.Ngas + 3;

	voronoignufile << "#!/bin/gnuplot" << endl;
	voronoignufile << "set terminal png" << endl;
	voronoignufile << "set output 'voronoi.png' "<< endl;
	voronoignufile << "set size square" << endl;
	voronoignufile << "set xrange[0.0:1.0]" << endl;
	voronoignufile << "set yrange[0.0:1.0]" << endl;
	voronoignufile << "set parametric"  << endl;
	voronoignufile << "set trange[0:1]" << endl;
	voronoignufile << "set trange[0.0:2*pi]" << endl;
	voronoignufile << "const1 = 0" << endl; 
	voronoignufile << "const2 = 1" << endl;
	//voronoignufile << "plot const1, t/2.0/pi ls 1 , const2,t/2.0/pi ls 1,t/2.0/pi, const1 ls 1,t/2.0/pi ,const2 ls 1, 0.707107*cos(t)+0.5, 0.707107*sin(t)+0.5 ls 1, -0.914214, -0.724745 ls 1, 1.914214, -0.724745 ls 1, 0.500000, 1.914214  ls 1, 'mesh.dat' u 2:3 w d, '-' w  l" << endl;
	voronoignufile << "plot 'mesh.dat' u 2:3 w d , '-' w l" << endl;
	voronoignufile << "-0.724745       -0.207107" << endl;
	voronoignufile << "1.724745        -0.207107" << endl;
	voronoignufile << "0.500000        1.914214" << endl;
	voronoignufile << "-0.724745       -0.207107" << endl;
	voronoignufile << endl;

	for(int imesh = 3; imesh < Nmesh; ++imesh)
	{
	  for(int j = 0; j < Mesh[imesh].neighbors.size(); ++j)
	  {
	    for(int k = 0; k < Mesh[imesh].neighbors[j].delaunay.size(); ++k)
	    {
	      voronoignufile << Mesh[imesh].neighbors[j].delaunay[k]->getxc() << " " <<
		Mesh[imesh].neighbors[j].delaunay[k]->getyc() << endl;
	      voronoifile << Mesh[imesh].neighbors[j].delaunay[k]->getxc() << "\t" << 
		Mesh[imesh].neighbors[j].delaunay[k]->getyc() << endl;
	    }
	    voronoignufile << endl;
	    voronoifile << endl;
	  }
	}


	break;
      }

    case(3):
      {
	const int Nmesh = params.Ngas + 4 + params.Nghost;

	voronoignufile << "#!/bin/gnuplot" << endl;
	voronoignufile << "set term png" << endl;
	voronoignufile << "set output " << "'voronoi.png'" << endl;
	voronoignufile << "set xrange [0.0:1.0]" << endl;
	voronoignufile << "set yrange [0.0:1.0]" << endl;
	voronoignufile << "set zrange [0.0:1.0]"  << endl;
	voronoignufile << "splot '-' w l," << endl;
	
	for(int imesh = 4; imesh < Nmesh; ++imesh)
	{
	  string str ; vartostring(str, imesh);
	  string filename = "output/voronoi" + str + ".gnu";
	  ofstream voronoignufile2(filename);

	  voronoignufile2 << "#!/bin/gnuplot" << endl;
	  voronoignufile2 << "set term png" << endl;
	  voronoignufile2 << "set output " << "'" <<  str << ".png'" << endl;
	  voronoignufile2 << "set xrange [0.0:1.0]" << endl;
	  voronoignufile2 << "set yrange [0.0:1.0]" << endl;
	  voronoignufile2 << "set zrange [0.0:1.0]"  << endl;
	  voronoignufile2 << "splot '-' w l," << endl;
	  
	  for(int i = 0; i < Mesh[imesh].neighbors.size(); ++i)
	  {
	    string str2; vartostring(str2, i);
	    filename = "output/voronoi" + str + "_" + str2  + ".gnu"; 
	    //ofstream voronoignufile3(filename);
	    for(int id = 0; id < Mesh[imesh].neighbors[i].delaunay.size(); ++id)
	    {
	      delaunayClass *pdel = Mesh[imesh].neighbors[i].delaunay[id];
	      delaunayClass *pnbr = NULL;
	      int jmesh = Mesh[imesh].neighbors[i].index;

	      for(int inbr = 0; inbr < pdel->nNeighbor; ++inbr)
	      {
		pnbr = pdel->neighbors[inbr];
		//if(jmesh < 4) continue;
		if(pnbr->has_mesh(imesh) && pnbr->has_mesh(jmesh) && pdel->has_mesh(imesh) && pdel->has_mesh(jmesh))
		{
		  //cout << "write "  << imesh << " " << jmesh << endl;
		  /*
		  voronoignufile3 << pdel->getxc() << "\t" << pdel->getyc() << "\t" << pdel->getzc() << endl;
		  voronoignufile3 << pnbr->getxc() << "\t" << pnbr->getyc() << "\t" << pnbr->getzc() << endl;
		  voronoignufile3 << endl; 
		  voronoignufile3 << endl;
		   */

		  voronoignufile2 << pdel->getxc() << "\t" << pdel->getyc() << "\t" << pdel->getzc() << endl;
		  voronoignufile2 << pnbr->getxc() << "\t" << pnbr->getyc() << "\t" << pnbr->getzc() << endl;
		  voronoignufile2 << endl; 
		  voronoignufile2 << endl;

		  voronoignufile << pdel->getxc() << "\t" << pdel->getyc() << "\t" << pdel->getzc() << endl;
		  voronoignufile << pnbr->getxc() << "\t" << pnbr->getyc() << "\t" << pnbr->getzc() << endl;
		  voronoignufile << endl; 
		  voronoignufile << endl;



		}
	      }


	    }
	    //voronoignufile3.close();
	  }

	  voronoignufile2 << "e" << endl;
	  voronoignufile2.close();

	}

	voronoignufile << "e" << endl;
	voronoignufile.close();

	break; 

      }
  }

  cout << "output => voronoi.gnu, output/voronoi.dat" << endl;
  voronoignufile << "e" << endl;
  voronoifile.close();
  voronoignufile.close();
}

void outputParticle(string dmfilename, string gasfilename, particleClass *Particle, paramsClass &params)
{
  const int Ngas = params.Ngas;
  const int Ndm = params.Ndm;
  const int buf = params.dim == 2 ? 3 : 4;
  ofstream dmfile(dmfilename);
  ofstream gasfile(gasfilename);

  for(int i = buf; i < Ngas + Ndm + buf; ++i)
  {
    particleClass *part = &Particle[i]; 
    if(part->getParticleType() == dark_matter)
    {
      dmfile << dark_matter << "\t"<< i << "\t" << part->getMass() << "\t" << 
	part->getPosx() << "\t" << part->getPosy() << "\t" << part->getPosz() << endl; 
    } else if(part->getParticleType() == gas)
    {
      //gasfile << ... << endl; 
      //gasfile << gas << "\t" << i << "\t" << part->getPosx() << "\t" << part->getPosy() << "\t" << part->getPosz() << endl;	

      gasfile << i << "\t" 
	      << part->getMass() << "\t" 
	      << scientific << part->getPosx() << "\t" 
	      << scientific << part->getPosy() << "\t" 
	      << scientific << part->getPosz() << "\t" << "0\t0\t0\t0\t0" << endl;
    }
  } 

  cout << "output => dm and gas file" << endl;
  dmfile.close();
  gasfile.close();
}
