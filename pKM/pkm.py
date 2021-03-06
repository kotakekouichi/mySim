import os
import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class posClass:
  def __init__(self, index, m, x, y, z):
    self.index = index;
    self.m = m;
    self.x = x;
    self.y = y;
    self.z = z;

class pkm:
  def __init__(self, dirpath):
    self.dirpath = dirpath
    self.voronoiEdge = []
    self.delaunayEdge = []
    self.particleDic = {}

    self.dict_voronoiEdge = {}
    self.dict_delaunayEdge = {}
    
    self.readparticlefile(dirpath);
    self.readvoronoifile(dirpath);
    self.readdelaunayfile(dirpath);

  def load(self, filepath):
    return pkm(filepath);

  def readparticlefile(self, dirpath):
#read dark matter and gas particle file 
    self.particleDic["dark_matter"] = []
    self.particleDic["gas"] = []
    
    if(os.path.exists(dirpath + "/dm.dat") == True):
      print "file exits => " + dirpath + "/dm.dat"; 

      for line in open(dirpath + "/dm.dat"):
        itemList = line[:-1].split('\t'); 
	partiletype = int(itemList[0]);
	index = int(itemList[1]);
	mp = float(itemList[2]);
	xp = float(itemList[3]);
	yp = float(itemList[4]);
	zp = float(itemList[5]);
	self.particleDic["dark_matter"].append(posClass(index, mp, xp, yp, zp));
     
    else:
      print "not exits => " + dirpath + "/particle.dat"   
    
    if(os.path.exists(dirpath + "/gas.dat") == True):
      print "file exits => " + dirpath + "/gas.dat"

      for line in open(dirpath + "/gas.dat"):
        itemList = line[:-1].split('\t');
#partiletype = int(itemList[0]);
	index = int(itemList[0]);
	mp = float(itemList[1]);
	xp = float(itemList[2]);
	yp = float(itemList[3]);
	zp = float(itemList[4]);
	self.particleDic["gas"].append(posClass(index, mp, xp, yp, zp));

    else:
        print "not exits => " + dirpath + "/gas.dat"   

  def readvoronoifile(self, dirpath):

    if(os.path.exists(dirpath + "/voronoi") == True):
      print "read voronoi dir...=>" + dirpath + "/voronoi"

      for pos in self.particleDic["gas"]:
        index = pos.index
	filepath = dirpath + "/voronoi/voronoi" + str(index) + ".dat"

	if(os.path.exists(filepath) == True):

	  print "====> read voronoi file ...=>" + filepath 
	  self.dict_voronoiEdge[index] = []
	  posPair = []
	  for line in open(filepath, "r"):
	    itemList = line[:-1].split('\t');
	    if(len(itemList) == 1):
	      if(len(posPair) > 0):
		self.dict_voronoiEdge[index].append(posPair);
	      posPair = []
	      continue;
	    m = pos.m
	    x = float(itemList[0]);
	    y = float(itemList[1]);
	    z = 0.0;
	    if(len(itemList) == 3):
	      z = float(itemList[2]);
	    posPair.append(posClass(-1, m, x, y, z));

      filepath = dirpath + "/voronoi/voronoi.dat"

      if(os.path.exists(filepath) == True):
	self.dict_voronoiEdge["all"] = []
	posPair = []
	for line in open(filepath, "r"):
          itemList = line[:-1].split('\t');
	  if(len(itemList) == 1):
	    if(len(posPair) > 0):
              self.dict_voronoiEdge["all"].append(posPair);
	    posPair = []
	    continue;
	  m = pos.m
	  x = float(itemList[0]);
	  y = float(itemList[1]);
	  z = 0.0;
	  if(len(itemList) == 3):
	    z = float(itemList[2]);
	  posPair.append(posClass(-1, m, x, y, z));

    return;

  def readdelaunayfile(self, dirpath):
    
    if(os.path.exists(dirpath + "/delaunay") == True):
      print "read delaunay dir...=>" + dirpath + "/delaunay"

      for pos in self.particleDic["gas"]:
        index = pos.index
	filepath = dirpath + "/delaunay/delaunay" + str(index) + ".dat"

	if(os.path.exists(filepath) == True):

	  print "====> read delaunay file ...=>" + filepath 
	  self.dict_delaunayEdge[index] = []
	  posPair = []
	  for line in open(filepath, "r"):
	    itemList = line[:-1].split('\t');
	    if(len(itemList) == 1):
	      if(len(posPair) > 0):
		self.dict_delaunayEdge[index].append(posPair);
	      posPair = []
	      continue;
	    m = pos.m
	    x = float(itemList[0]);
	    y = float(itemList[1]);
	    z = 0.0;
	    if(len(itemList) == 3):
	      z = float(itemList[2]);
	    posPair.append(posClass(-1, m, x, y, z));

      filepath = dirpath + "/delaunay/delaunay.dat"

      if(os.path.exists(filepath) == True):
	self.dict_delaunayEdge["all"] = []
	posPair = []
	for line in open(filepath, "r"):
          itemList = line[:-1].split('\t');
	  if(len(itemList) == 1):
	    if(len(posPair) > 0):
              self.dict_delaunayEdge["all"].append(posPair);
	    posPair = []
	    continue;
	  m = pos.m
	  x = float(itemList[0]);
	  y = float(itemList[1]);
	  z = 0.0;
	  if(len(itemList) == 3):
	    z = float(itemList[2]);
	  posPair.append(posClass(-1, m, x, y, z));

    return;

class plotClass(plt.Axes):
  def __init__(self):
    self.name = "";

def load(dirpath):
  return pkm(dirpath);

def plotVoronoiMesh2d(Pkm, *arguments):

  plt.figure();

  if(not arguments):
    indexes = ["all"] 
  else:
    indexes = arguments

  for index in indexes:
    
    if (index in Pkm.dict_voronoiEdge ) == False:
      return plt;
    for edge in Pkm.dict_voronoiEdge[index]:
      plt.plot([edge[0].x, edge[1].x],[edge[0].y, edge[1].y],'r');
    plt.xlim(0,1);
    plt.ylim(0,1);

#plt.show();

  return plt;

def plotVoronoiMesh3d(Pkm, *arguments):

  if(not arguments):
    indexes = ["all"] 
  else:
    indexes = arguments

  fig = plt.figure();
  ax = Axes3D(fig)
  
  for index in indexes:
    if (index in Pkm.dict_voronoiEdge ) == False:
      return plt;
    for edge in Pkm.dict_voronoiEdge[index]:
      ax.plot([edge[0].x, edge[1].x], [edge[0].y, edge[1].y], [edge[0].z, edge[1].z], 'r-')

#plt.xlim(0, 1);
#plt.ylim(0, 1);
#plt.zlim(0, 1);

#plt.show();
  
  return plt;

def plotDelaunayMesh3d(Pkm, *arguments):

  if(not arguments):
    indexes = ["all"] 
  else:
    indexes = arguments

  plot.figure();
  ax = Axes3D(fig)

  for index in indexes:
    if(index in Pkm.dict_delaunayEdge) == False:
      return plt;

    for edge in Pkm.dict_delaunayEdge[index]:
      ax.plot([edge[0].x, edge[1],x, edge[2],x, edge[3].x, edge[4].x], [edge[0].y, edge[1].y, edge[2].y, edge[3].y, edge[4].y], [edge[0].z, edge[1].z, edge[2].z, edge[3].z, edge[4],z], 'r-')

#  plt.xlim(0, 1);
#  plt.ylim(0, 1);
#  plt.zlim(0, 1);

#plt.show();

  return plt;
 
def plotDelaunayMesh2d(Pkm, *arguments):

  if(not arguments):
    indexes = ["all"] 
  else:
    indexes = arguments

  plt.figure();

  for index in indexes:
    if(index in Pkm.dict_delaunayEdge) == False:
      return plt;

    for edge in Pkm.delaunayEdge:
      plt.plot([edge[0].x, edge[1].x, edge[2].x, edge[3].x], [edge[0].y, edge[1].y, edge[2].y, edge[3].y], 'r-');
    plt.xlim(0,1);
    plt.ylim(0,1);
#plt.show();
  return plt;

def projectionPlot(Pkm, xyz, key, NMESH):

  particleList = Pkm.particleDic[key];

  heatmap = [0.0000000001] * NMESH;
  for i in range(NMESH):
     heatmap[i] = [0.0000000001]*NMESH;

  dx = 1.0 / NMESH;
  for part in particleList:
    i = j = -1;
    if(xyz == "x"):
      i = int(part.y / dx);
      j = int(part.z / dx);
    elif(xyz == "y"):
      i = int(part.x / dx);
      j = int(part.z / dx);
    elif(xyz == "z"):
      i = int(part.x / dx);
      j = int(part.y / dx);
  
    index = int(i + j * NMESH);
    heatmap[i][j] += float(part.m);

  extent = [1.0,0.0,1.0,0.0]
  plt.imshow(heatmap, extent = extent);
  plt.show();

  return;
