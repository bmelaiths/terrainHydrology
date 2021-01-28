#!/usr/bin/env python
# coding: utf-8

# #### Boiler Plate

# In[ ]:


import random
import cv2 as cv
import networkx as nx
import numpy as np
import timeit
import datetime

from noise import pnoise2, snoise2
from matplotlib import pyplot as plt
from scipy.spatial import cKDTree
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy import interpolate
from poisson import PoissonGenerator

from PIL import Image

startTime = datetime.datetime.now()


# ### Data Structures

# In[ ]:


class ShoreModel:
    def __init__(self, inputFileName):
        self.img = cv.imread(inputFileName)
        
        self.imgray = cv.cvtColor(self.img, cv.COLOR_BGR2GRAY) # a black-and-white version of the input image
        self.shape = self.imgray.shape
        ret, thresh = cv.threshold(self.imgray, 127, 255, 0)
        contours, hierarchy = cv.findContours(thresh, cv.RETR_LIST, cv.CHAIN_APPROX_NONE)
        test = cv.cvtColor(thresh, cv.COLOR_GRAY2BGR) # not sure what this line does or why it exists
        if len(contours) > 1:
            print('WARNING: Multiple contours identified. The program may not have correctly')
            print('identified the land.')
        self.contour = contours[0]
        self.contour=self.contour.reshape(-1,2)
        self.contour=np.flip(self.contour,1)
        
        self.imgOutline = self.img.copy()
        cv.drawContours(self.imgOutline, contours, -1, (0,255,0), 2)
        
        # TODO raise exception if dimensions not square
        # TODO raise exception if multiple contours
    def distanceToShore(self, loc):
        #    for some reason this method is      y, x
        return cv.pointPolygonTest(self.contour,(loc[1],loc[0]),True)
    def isOnLand(self, loc):
        if 0 <= loc[0] < self.shape[0] and 0 <= loc[1] < self.shape[1]:
            return self.imgray[int(loc[1])][int(loc[0])] != 0
        else:
            return False
    def __getitem__(self, index):
        # TODO ensure that points returned are x,y
        return (self.contour[index][1],self.contour[index][0])
    def __len__(self):
        return len(self.contour)
    def tobytes(self):
        

class HydrologyNetwork:
    def __init__(self):
        self.nodeCounter = 0
        self.graph = nx.DiGraph()
        self.mouthNodes = []
    def addNode(self, loc, elevation, priority, contourIndex=None, parent=None):
        node = HydroPrimitive(self.nodeCounter, loc, elevation, priority, parent)
        if contourIndex is not None:
            node.contourIndex = contourIndex
            self.mouthNodes.append(self.nodeCounter)
        self.graph.add_node(
            self.nodeCounter,
            primitive=node
        )
        if parent is not None:
            self.graph.add_edge(parent.id,self.nodeCounter)
        self.nodeCounter += 1
        allpoints_list = np.array([self.graph.nodes[n]['primitive'].position for n in range(self.nodeCounter)])
        self.graphkd = cKDTree(allpoints_list)
        return node
    def query_ball_point(self, loc, radius):
        return self.graphkd.query_ball_point(loc,radius)
    def edgesWithinRadius(self, loc, radius):
        nodesToCheck = self.graphkd.query_ball_point(loc,radius)
        edges = [ self.graph.out_edges(n) for n in nodesToCheck ]
        edges = [item for edge in edges for item in edge]
        return [(self.graph.nodes[e[0]]['primitive'],self.graph.nodes[e[1]]['primitive']) for e in edges]
    def downstream(self, node):
        parent = list(self.graph.predecessors(node))
        if len(parent) > 0:
            return self.graph.nodes[parent[0]]['primitive']
        else:
            return None
    def upstream(self, node):
        return [self.graph.nodes[n]['primitive'] for n in list(self.graph.successors(node))]
    def adjacentNodes(self, node):
        downstream = self.downstream(node)
        upstream   = self.upstream(node)
        if downstream is None:
            downstream = [ ]
        return upstream + downstream
    def allUpstream(self, node):
        return [self.graph.nodes[n]['primitive'] for n in nx.descendants(self.graph, node)]
    def allNodes(self):
        return [self.graph.nodes[node]['primitive'] for node in list(self.graph.nodes)]
    def allMouthNodes(self):
        return [self.graph.nodes[id]['primitive'] for id in self.mouthNodes]
    def allLeaves(self, node):
        ids = [s for s in nx.descendants(self.graph,node) if len(self.graph.out_edges(s))==0]
        return [self.graph.nodes[id]['primitive'] for id in ids]
    def node(self, node):
        return self.graph.nodes[node]['primitive']
    def dfsPostorderNodes(self):
        ids = list(nx.dfs_postorder_nodes(self.graph))
        return [self.graph.nodes[id]['primitive'] for id in ids]
    def pathToNode(self, origin, destination):
        return [self.graph.nodes[n]['primitive'] for n in nx.shortest_path(self.graph,origin,destination)]
    def __len__(self):
        return len(self.graph.nodes)

class HydroPrimitive:
    def __init__(self, id, loc, elevation, priority, parent):
        self.id = id
        self.position = loc
        self.elevation = elevation
        self.priority = priority
        self.parent = parent
        self.rivers = [ ]
    def x(self):
        return self.position[0]
    def y(self):
        return self.position[1]

def openCVFillPolyArray(points):
    return [ np.array( [ [int(p[0]),int(p[1])] for p in points ] ) ]

class Q:
    def __init__(self, position, nodes, iv):
        self.position = position
        self.nodes = nodes
        self.vorIndex = iv # the index of the voronoi vertex this represents in vor.vertices

class TerrainHoneycomb:
    def __init__(self, shore, hydrology):
        self.shore     = shore
        self.hydrology = hydrology
        
        points = [node.position for node in hydrology.allNodes()]
        points.append((0,0))
        points.append((0,shore.shape[1]))
        points.append((shore.shape[0],0))
        points.append((shore.shape[0],shore.shape[1]))
        
        self.vor = Voronoi(points,qhull_options='Qbb Qc Qz Qx')
        
        self.imgvoronoi = np.zeros(shore.shape, dtype=np.uint16)
        for n in range(len(hydrology)):
            if self.vor_region_id(n) == -1:
                continue
            positions = self.ridgePositions(n)
            cv.fillPoly(
                self.imgvoronoi,
                openCVFillPolyArray(positions),
                np.int16(self.vor_region_id(n)+1).item()
            )
        ret, thresh = cv.threshold(shore.imgray, 127, 1, 0)
        mask = np.array(thresh, dtype=np.uint16) * (256*256-1)
        self.imgvoronoi = cv.bitwise_and(self.imgvoronoi,mask)
        
        self.qs = [ ]
        for iv in range(len(self.vor.vertices)):
            if not shore.isOnLand(self.vor.vertices[iv]):
                self.qs.append(None)
                continue
            regionIdxes = [self.vor.regions.index(bound) for bound in self.vor.regions if iv in bound]
            nodeIdxes = [list(self.vor.point_region).index(regionIndex) for regionIndex in regionIdxes]
            self.qs.append(Q(self.vor.vertices[iv],nodeIdxes, iv))
        
        self.cellsRidges = { }
        for n in range(len(hydrology)):
            connectedNodes = hydrology.upstream(n)
            if node.parent is not None:
                connectedNodes.append(node.parent.id)
            
            verts = self.vor.regions[self.vor_region_id(n)].copy()
            ridges = [ ]
            for ri in range(len(self.vor.ridge_points)):
                if n not in self.vor.ridge_points[ri]:
                    continue
                # ri points to a ridge of this node
                if self.vor.ridge_points[ri][0] in connectedNodes or self.vor.ridge_points[ri][1] in connectedNodes:
                    continue
                # ri does not point to a ridge that an edge goes through
                vertex0 = self.qs[self.vor.ridge_vertices[ri][0]]
                vertex1 = self.qs[self.vor.ridge_vertices[ri][1]]
                if vertex0 is None or vertex1 is None:
                    continue
                # both ends of the ridge are on land
                if vertex0.vorIndex in verts:
                    verts.remove(vertex0.vorIndex)
                if vertex1.vorIndex in verts:
                    verts.remove(vertex1.vorIndex)
                ridges.append((vertex0, vertex1))
            ridges += [(self.qs[vert],) for vert in verts if self.qs[vert] is not None]
            # ridges includes vertices that are not part of a ridge (but only vertices on land)
            self.cellsRidges[n] = ridges
    def vor_region_id(self, node):
        return self.vor.point_region[node]
    def ridgePositions(self, node):
        ridges = self.vor.regions[self.vor_region_id(node)] # the indices of the vertex boundaries
        return [self.vor.vertices[x] for x in ridges if x != -1] # positions of all the vertices
    def cellArea(self, loc):
        return np.count_nonzero(self.imgvoronoi == self.imgvoronoi[int(loc[1])][int(loc[0])])
    def cellQs(self, node):
        return [self.qs[vorIdx] for vorIdx in self.vor.regions[self.vor_region_id(node)] if self.qs[vorIdx] is not None]
    def allQs(self):
        return self.qs.copy()
    def boundingBox(self, n):
        idxes = np.where(self.imgvoronoi==self.vor.point_region[n]+1) # coordinates of all pixels in the voronoi region
        xllim = min(x for x in idxes[0]) # these lines get the bounding box of the voronoi region
        xulim = max(x for x in idxes[0])
        yllim = min(x for x in idxes[1])
        yulim = max(x for x in idxes[1])
        # I don't know why he creates another bounding box with opencv
        b = np.array([[xllim,yllim],[xllim,yulim],[xulim,yllim],[xulim,yulim]])
        b = cv.minAreaRect(b)
        pts = cv.boxPoints(b)
        xllim = int(min(x[0] for x in pts))
        xulim = int(max(x[0] for x in pts))
        yllim = int(min(x[1] for x in pts))
        yulim = int(max(x[1] for x in pts))
        return (xllim, xulim, yllim, yulim)
    def isInCell(self, p, n):
        return self.imgvoronoi[int(p[1])][int(p[0])]==self.vor.point_region[n]+1
    def cellRidges(self, n):
        return self.cellsRidges[n]
    def nodeID(self, point):
        id = list(self.vor.point_region).index(self.imgvoronoi[point[1]][point[0]]-1)
        return id if id != -1 else None

class T:
    def __init__(self, position, cell):
        self.position = position
        self.cell = cell

class Terrain:
    def __init__(self, hydrology, cells):
        self.cellTs = { }
        self.tList = [ ]
        
        disk = False                # this parameter defines if we look for Poisson-like distribution on a disk/sphere (center at 0, radius 1) or in a square/box (0-1 on x and y)
        repeatPattern = True        # this parameter defines if we look for "repeating" pattern so if we should maximize distances also with pattern repetitions
        num_iterations = 4          # number of iterations in which we take average minimum squared distances between points and try to maximize them
        first_point_zero = False    # should be first point zero (useful if we already have such sample) or random
        iterations_per_point = 128  # iterations per point trying to look for a new point with larger distance
        sorting_buckets = 0         # if this option is > 0, then sequence will be optimized for tiled cache locality in n x n tiles (x followed by y)
        num_dim = 2                 # 1, 2, 3 dimensional version
        num_rotations = 1           # number of rotations of pattern to check against
        
        poisson_generator = PoissonGenerator( repeatPattern, first_point_zero)
        points = poisson_generator.find_point_set(num_points, num_iterations, iterations_per_point, num_rotations)
        for n in range(len(hydrology)):
            xllim, xulim, yllim, yulim = cells.boundingBox(n)
            
            # I think this applies a mask to the poisson points, and adds those points as Tees for the cell
            points_projected = [ [p[0]*(yulim-yllim)+yllim,p[1]*(xulim-xllim)+xllim] for p in points ]
            points_filtered = [ (p[0],p[1]) for p in points_projected if cells.isInCell(p,n) ]
            cellTs = [T(p,n) for p in points_filtered]
            self.cellTs[n] = cellTs
            self.tList += cellTs
            
            clear_output(wait=True)
            print(n," out of ",len(nodes))  # use display(f) if you encounter performance issues
        
        allpoints_list = [[t.position[0],t.position[1]] for t in self.allTs()]
        allpoints_nd = np.array(allpoints_list)
        self.apkd = cKDTree(allpoints_nd)
    def allTs(self):
        return self.tList.copy()
    def cellTs(self, cell):
        return self.Ts[cell].copy()
    def query_ball_point(self, loc, radius):
        return [self.tList[i] for i in self.apkd.query_ball_point(loc,radius)]


# ### Global Parameters Of The Process

# In[ ]:


# Global Variables
Ps = 0.3 #0.05 ## probability of symetric branch
Pa = 0 #0.3 ## probability of asymetric branch
Pc = 1-(Ps+Pa) ## probability of continium (growth)
inputDomain='taiwan-outline-bigger.png'
inputTerrain='taiwan-terrain-bigger.png'
inputRiverSlope='taiwan-riverslope-bigger.png'
globalseed=4314
N_majorRivers=10
zeta = 10 # elevation range to include in candidate node selection
slopeRate = 30 # Maximum rate at which rivers climb
edgeLength = 30
eta = .75   #   eta * edgeLength is the minimum distance from a node to the coast
sigma = .75 # sigma * edgeLength is the minimum distance between two nodes
rwidth=6 # 
riverAngleDev = 1.7 # Used in picknewnodepos(). Standard Deviation of angle for new node. Increase for less straight rivers
maxTries = 15
outputResolution = 400
num_points = 50 # The (rough) number of terrain primitives for each cell

radius = edgeLength * 3

random.seed(globalseed)


# #### Load Base Image 

# In[ ]:


shore = ShoreModel(inputDomain)
plt.imshow(shore.imgOutline)

terrainSlope = Image.open(inputTerrain)
terrainSlope = terrainSlope.convert('L')
plt.show()
print('Terrain Slope Input')
plt.imshow(terrainSlope)
terrainSlope = terrainSlope.load()

riverSlope = Image.open(inputRiverSlope)
riverSlope = riverSlope.convert('L')
plt.show()
print('River Slope Input')
plt.imshow(riverSlope)
riverSlope = riverSlope.load()


# #### Find contour of the ROI, This will be Gamma

# In[ ]:





# In[ ]:





# #### Select first point at a ranndom offset, try to select the following points such as that the highest probability ( on a nomal distribution) is that they are furthest away from each other 
# 

# In[ ]:


hydrology = HydrologyNetwork()

firstIdx = random.randint(0,len(shore)-1)
point = shore[firstIdx]
hydrology.addNode(point, 0, random.randint(1,N_majorRivers), contourIndex=firstIdx)

dist = len(shore)/N_majorRivers
for i in range(1,N_majorRivers):
    idx = int((firstIdx+i*dist+random.gauss(0, dist/6))%len(shore))
    point = shore[idx]
    hydrology.addNode(point, 0, random.randint(1,N_majorRivers), contourIndex=idx)

imgMouthDots = shore.imgOutline.copy()
for node in hydrology.allMouthNodes():
    cv.circle(imgMouthDots, (node.x(),node.y()), int((shore.shape[0]/512)*10), (255,0,0), -1)
plt.imshow(imgMouthDots)


# # Network Generation

# In[ ]:



# Borrowed , all of it
def segments_distance(a,b,c,d):
  """ distance between two segments in the plane:
      one segment is a to b
      the other is   c to d
  """
  #print(a[0],a[1],b[0],b[1],c[0],c[1],d[1],d[0])
  #print(segments_distance_internal(a[0],a[1],b[0],b[1],c[0],c[1],d[1],d[0]))
  return segments_distance_internal(a[0],a[1],b[0],b[1],c[0],c[1],d[0],d[1])

def segments_distance_internal(x11, y11, x12, y12, x21, y21, x22, y22):
  """ distance between two segments in the plane:
      one segment is (x11, y11) to (x12, y12)
      the other is   (x21, y21) to (x22, y22)
  """
  if segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22): return 0
  # try each of the 4 vertices w/the other segment
  distances = []
  distances.append(point_segment_distance(x11, y11, x21, y21, x22, y22))
  distances.append(point_segment_distance(x12, y12, x21, y21, x22, y22))
  distances.append(point_segment_distance(x21, y21, x11, y11, x12, y12))
  distances.append(point_segment_distance(x22, y22, x11, y11, x12, y12))
  return min(distances)

def segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22):
  """ whether two segments in the plane intersect:
      one segment is (x11, y11) to (x12, y12)
      the other is   (x21, y21) to (x22, y22)
  """
  dx1 = x12 - x11
  dy1 = y12 - y11
  dx2 = x22 - x21
  dy2 = y22 - y21
  delta = dx2 * dy1 - dy2 * dx1
  if delta == 0: return False  # parallel segments
  s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta
  t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / (-delta)
  return (0 <= s <= 1) and (0 <= t <= 1)

import math
# I think this finds the distance between a point and a line segment
def point_segment_distance(px, py, x1, y1, x2, y2):
  dx = x2 - x1
  dy = y2 - y1
  if dx == dy == 0:  # the segment's just a point
    return math.hypot(px - x1, py - y1)

  # Calculate the t that minimizes the distance.
  t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  # See if this represents one of the segment's
  # end points or a point in the middle.
  if t < 0:
    dx = px - x1
    dy = py - y1
  elif t > 1:
    dx = px - x2
    dy = py - y2
  else:
    near_x = x1 + t * dx
    near_y = y1 + t * dy
    dx = px - near_x
    dy = py - near_y

  return math.hypot(dx, dy)

def point_segment_distance_is_endpoint(px, py, x1, y1, x2, y2):
  dx = x2 - x1
  dy = y2 - y1
  if dx == dy == 0:  # the segment's just a point
    return (math.hypot(px - x1, py - y1), True)

  # Calculate the t that minimizes the distance.
  t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  isEndpoint = True

  # See if this represents one of the segment's
  # end points or a point in the middle.
  if t < 0:
    dx = px - x1
    dy = py - y1
  elif t > 1:
    dx = px - x2
    dy = py - y2
  else:
    near_x = x1 + t * dx
    near_y = y1 + t * dy
    dx = px - near_x
    dy = py - near_y
    isEndpoint = False

  return (math.hypot(dx, dy), isEndpoint)

def point_segment_distance_tuple(p,a,b):
    return point_segment_distance(p[0],p[1],a[0],a[1],b[0],b[1])
    
def segments_intersect_tuple(a1,a2,b1,b2):
    return segments_intersect(a1[0],a1[1],a2[0],a2[1],b1[0],b1[1],b2[0],b2[1])

    


# In[ ]:





# In[ ]:





# In[ ]:


def selectNode(candidate_nodes,zeta):

    lowestCandidateZ = min([node.elevation for node in candidate_nodes]) # elevation of lowest candidate
    subselection = [n for n in candidate_nodes if n.elevation < lowestCandidateZ+zeta ] # 
    subselection.sort(key = lambda r : r.priority,reverse = True)
    subsubselection=[node for node in subselection if node.priority == subselection[0].priority]
    
    return subsubselection[0]


# In[ ]:


def alpha(node, candidates):         # alpha, as in the expansion rules in Table 1
    if node.priority==1:
        ruleBase(node, candidates)
    else:
        Pval = random.random();
        if Pval <= Pa:
            rulePa(node, candidates)
        elif Pval <= Pa+Pc:
            rulePc(node, candidates)
        else:
            rulePs(node, candidates)

            
def ruleBase(node, candidates): #filling
    #tao(priority,node)
    for i in range(random.randint(1,5)):
        beta(node, node.priority, candidates)
        
        
def rulePc(node, candidates): #rive growth
    #tao(priority,node)
    beta(node, node.priority,candidates)
    
    
def rulePs(node, candidates): #symetric junction
    #tao(priority,node)
    beta(node, node.priority-1, candidates)
    beta(node, node.priority-1, candidates)

    
def rulePa(priority,node): # asymetric junction
    #tao(priority,node)
    beta(node, node.priority, candidates)
    beta(random.randint(1,priority-1),node)
    
    
def beta(node, priority, candidates):
    point = picknewnodepos(node)
    if point is not None:
        slope = 2.0 * riverSlope[ int(node.x()) , int(node.y())] / 255
        newZ = node.elevation + random.random() * slopeRate * slope
        candidates.append(hydrology.addNode(point, priority=priority, elevation=newZ, parent=node))
    else:
        tao(node, candidates)

def tao(node, candidates):
    try:
        candidates.remove(node)
    except:
        None
    finally:
        None

def picknewnodepos(parentnode):
    parentsparent = hydrology.downstream(parentnode.id) # parent node of parentnode
    
    angle = None
    if parentsparent is None: # If there is no previous node, Go in a direction perpendicular to the coast
        angle = coastNormal(parentnode)
        if angle is None:
            return None
    else:
        # 'angle' is the 'direction' of the river
        angle = math.atan2( # y,x !
            parentnode.y() - parentsparent.y(), # y
            parentnode.x() - parentsparent.x()  # x
        )
    
    newNodePos = None
    # Try maxTries number of times to get a suitable point
    for i in range(maxTries):
        # Pick a random new point (generally in the same direction)
        newAngle = angle + random.gauss(0,riverAngleDev)
        newNodePos = (
            parentnode.x() + edgeLength*math.cos(newAngle),
            parentnode.y() + edgeLength*math.sin(newAngle)
        )
        if isAcceptablePosition(newNodePos):
            break
        else:
            newNodePos = None
    
    return newNodePos

def coastNormal(node): # Gets the angle that is approximately normal to the coast
    assert node.contourIndex is not None # assert that this is a mouth node
    p1 = shore[node.contourIndex+3]
    p2 = shore[node.contourIndex-3]
    theta = math.atan2(p2[1]-p1[1],p2[0]-p1[0])
    return theta + 0.5*math.pi

def isAcceptablePosition(point):
    if point is None:
        return False
    # is the point too close to the seeeeeeeeeeeee?
    if shore.distanceToShore(point) < eta*edgeLength:
        #print(f'Angle too close to the seeee (distance {cv.pointPolygonTest(contour,(point[1],point[0]),True)})')
        return False
    # is the point too close to other nodes?
    for node0,node1 in hydrology.edgesWithinRadius(point, 2*edgeLength): # Go through each edge
        dist = point_segment_distance( # Distance to the edge (edge is a line segment)
            point[0],point[1],
            node0.x(),node0.y(), # Line segment endpoint 1 (x,y)
            node1.x(),node1.y()  # Line segment endpoint 2 (x,y)
        )
        if dist < sigma*edgeLength:
            return False
    # otherwise return True
    return True


# In[ ]:





# In[ ]:


def calculateHorton_Strahler(selectedCandidate):
    #find the leaves from this node and calculate upstream
    leafs = hydrology.allLeaves(selectedCandidate.id)
    workingqueue=leafs
    nextQueue=set()
    while len(workingqueue)>0:
        nextQueue=set()
        for i in range(len(workingqueue)):
            priority=1
            children = hydrology.upstream(workingqueue[i].id)
            childrenPriorities = [child.priority for child in children]
            if len(children)>0:
                priority = max(childrenPriorities)
                if childrenPriorities.count(priority)>1:
                    priority=priority+1
            hydrology.node(workingqueue[i].id).priority = priority;
            parent = hydrology.downstream(workingqueue[i].id)
            if parent is not None:
                nextQueue.add(parent)
        workingqueue=list(nextQueue)

# I don't think this method is called anywhere
def calculateHorton_Strahler_():
    leafs = [x for x in G.nodes() if G.out_degree(x)==0]
    workingqueue=leafs
    nextQueue=set()
    while len(workingqueue)>0:
        nextQueue=set()
        for i in range(len(workingqueue)):
            priority=1
            children = G.successors(workingqueue[i])
            ChildrenPriorities = [G.nodes[x]['priority'] for x in children]
            if len(ChildrenPriorities)>0:
                priority = max(ChildrenPriorities)
                if ChildrenPriorities.count(priority)>1:
                    priority=priority+1
            G.nodes[workingqueue[i]]['priority']=priority;
            parent = G.predecessors(workingqueue[i])
            parent=[x for x in parent]
            if len(parent)==1:
                nextQueue.add(parent[0])
        workingqueue=list(nextQueue)


# In[ ]:


# I think this is where the rivers are built
candidates = hydrology.allMouthNodes()
from IPython.display import clear_output, display
while len(candidates)!=0:
    s= datetime.datetime.now()
    selectedCandidate = selectNode(candidates,zeta)
    a=datetime.datetime.now()
    alpha(selectedCandidate, candidates)
    b=datetime.datetime.now()
    calculateHorton_Strahler(selectedCandidate)
    c=datetime.datetime.now()
    clear_output(wait=True)
    #print("Select:   ",a-s) # time it takes to select a node
    #print("Expand:   ",b-a) # time it takes to expand the node
    #print("Classify: ",c-b) # time it takes to calculate the Horton-Strahler classification of the node
    print(len(hydrology))  # use display(f) if you encounter performance issues


# In[ ]:





# In[ ]:




# I don't think this segment runs

plt.figure(num=None, figsize=(16, 16), dpi=80)
pos=nx.get_node_attributes(G,'pos')
plt.imshow(img2)
#labels = list(map(lambda x: G.nodes[x]['priority'],nodes))
labels = list(map(lambda x: x,nodes))
labels = dict(zip(nodes,labels))
#nx.draw_networkx_labels(G,pos)
nx.draw(G,pos,node_size=60,labels=labels)
# In[ ]:


imgRiverHeights = np.zeros(shore.shape,dtype=np.uint16)

cells = TerrainHoneycomb(shore, hydrology)

pos = [node.position for node in hydrology.allNodes()]             # gets the 'pos' attributes of all nodes
labels = dict(zip(range(len(hydrology)),range(len(hydrology))))


for n in range(len(hydrology)):
    if cells.vor_region_id(n) == -1:
        continue
    positions = cells.ridgePositions(n)
    #print(f'Node ID: {n}, Positions: {positions}')
    cv.fillPoly(
        imgRiverHeights,
        openCVFillPolyArray(positions),
        np.int16(hydrology.node(n).elevation).item()
    )

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)
ax.imshow(shore.img)
ylim=ax.get_ylim();
xlim=ax.get_xlim();
nx.draw(hydrology.graph,pos,node_size=60,labels=labels,ax=ax)
voronoi_plot_2d(cells.vor, point_size=10, ax=ax,line_colors=['yellow']) # draws the voronoi cells?
ax.set_ylim(ylim);
ax.set_xlim(xlim);
kernel = cv.getStructuringElement(cv.MORPH_RECT,(2,2))#I have no idea what this is, and it isn't used anywhere else
plt.savefig('riverCellNetwork.png', dpi=100)
plt.show()

plt.show()
print("Same as last image but with the mask applied")
plt.imshow(cells.imgvoronoi)

plt.show()
print('River Elevations')
plt.imshow(imgRiverHeights, cmap=plt.get_cmap('terrain'))

### Breakdown of the image
# Giant red circles identify river mouths
# Blue dots identify river nodes
# Black arrows point upstream
# Black numbers identify the order of the nodes
# Green outline identifies the coast
# Yellow lines outline the voronoi cells around each river node
# Yellow dots identify the vertices of the voronoi cells


# In[ ]:





# In[ ]:





# ### voronoiimg has the watersheds where pixels of value  centernodeidx + 1 are pixels that belong to watershed of node centernodeidx, since there is only one s for each node ( that is the incoming edge from the parent ), we can calculate the watershed areas and store them into each node

# In[ ]:


# calculate watershed areas
for n in range(len(hydrology)):
    node = hydrology.node(n)
    node.localWatershed = cells.cellArea(node.position)
    node.inheritedWatershed = 0


# In[ ]:


# calculate the total area of the watershed behind the node
# also calculate the flow through the node
for node in hydrology.dfsPostorderNodes():  # search nodes in a depth-first post-ordering manner
    watershed = node.localWatershed + sum([n.inheritedWatershed for n in hydrology.allUpstream(node.id)])
    node.inheritedWatershed=watershed                         # calculate total watershed area
    node.flow = 0.42 * watershed**0.69 # calculate river flow


# In[ ]:


for q in cells.allQs():
    if q is None:
        continue
    nodes = [hydrology.node(n) for n in q.nodes]
    maxElevation = max([node.elevation for node in nodes])
    d = np.linalg.norm(q.position - nodes[0].position)
    slope = terrainSlope[int(q.position[0]),int(q.position[1])] / 255
    q.elevation = maxElevation + d * slope


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def classify(node):
    # Based on river slope and distance from Gamma
    # TODO: A real classification
    children = hydrology.upstream(node.id)
    for child in children:
        grade = (child.elevation - node.elevation) / edgeLength
        if grade > 0.1:
            child.rosgen = 'A+'
        elif grade >0.04:
            child.rosgen = 'A'
        elif grade > 0.02:
            child.rosgen = ['G','D','B'][random.randint(0,2)];
        elif grade > 0.005:
            child.rosgen = ['C','D','E','F'][random.randint(0,3)];
        else :
            child.rosgen = 'DA'
            
from itertools import islice

def window(seq, n=2): ##Borrowed as is
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
def clean_asin(sinAngle): ## Borrowed but modified
    return math.asin(min(1,max(sinAngle,-1)))


# In[ ]:


for n in range(len(hydrology)):
    classify(hydrology.node(n))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# This is a node graph, like the voronoi graph earlier. But the arrows are weighted by flow

plt.figure(num=None, figsize=(16, 16), dpi=80)
nodes = hydrology.allNodes()
ids = [node.id for node in nodes]
positions = [node.position for node in nodes]
normalizer = max([node.flow for node in nodes])
weights = [ 6*node.flow/normalizer for node in nodes if node.parent is not None ]

plt.imshow(shore.img)

nx.draw(hydrology.graph,positions,node_size=10,labels=labels,width=weights)


# In[ ]:


# Draw a map of rivers


# In[ ]:


# This is the least intuitive way to calculate Euclidean distance, but I respect it
def distance(a,b):
   return np.linalg.norm( np.subtract(a , b))


# In[ ]:


# This just gets the highest elevation in the entire map

maxNodeElevation = max([node.elevation for node in hydrology.allNodes()])
print(f'Highest Node: {maxNodeElevation}')
maxqElevation = max([q.elevation for q in cells.allQs() if q is not None])
print(f'Highest Ridge Elevation: {maxqElevation}')
maxZ = max([maxNodeElevation,maxqElevation])


# In[ ]:


# Same thing, but over imgvoronoi instead of the map

plt.figure(num=None, figsize=(16, 16), dpi=80)

plt.imshow(imgRiverHeights, plt.get_cmap('terrain'))

nx.draw(hydrology.graph,positions,node_size=60,labels=labels,width=weights)


# In[ ]:


Ts = Terrain(hydrology, cells)


# In[ ]:





# In[ ]:


fig = plt.figure(figsize=(16,16))

ax = fig.add_subplot(111)
ax.imshow(cells.imgvoronoi)
ax.scatter(*zip(*[t.position for t in Ts.allTs()]), color='r', alpha=0.6, lw=0)

ylim=ax.get_ylim();
xlim=ax.get_xlim();
nx.draw(hydrology.graph,positions,node_size=60,labels=labels,ax=ax)
voronoi_plot_2d(cells.vor, point_size=10, ax=ax,line_colors=['yellow'])

ax.set_ylim(ylim);
ax.set_xlim(xlim);
plt.show()


# In[ ]:


def projection(p,u,v):
    n = np.subtract(u,v)
    nnorm = np.linalg.norm(n, 2)
    n = n/nnorm
    ret = np.dot(np.subtract(p,v), n)
    proj = ret/nnorm
    if proj >1 :
        proj=1
    if proj <0 :
        proj =0
    return proj*n


# In[ ]:


import shapely.geometry as geom
import numpy as np
from shapely.geometry import asLineString

point = geom.Point(0.8, 10.5)

# Note that "line.distance(point)" would be identical

# This is a dictionary of nodes. I think it stores
# the path to the seeeee for each node
rivers = {}

for node in hydrology.allMouthNodes():
    # remember that out_edges gets upstream nodes
    leaves = hydrology.allLeaves(node.id)
    for leafNode in leaves: # essentially, this loops through all the highest nodes of a particular mouth
        # path to the leaf (there's only one, so it's the shortest)
        path = hydrology.pathToNode(node.id,leafNode.id)

        x = np.array([p.x() for p in path])
        y = np.array([p.y() for p in path])
        z = np.array([p.elevation for p in path])
        
        # it seems to me that, if the path is short, this block
        # adjusts the positions of the first three nodes
        if len(x)<4:
            x1 = (x[0]+x[1])/2
            x2 = (x[0]+x1)/2
            y1 = (y[0]+y[1])/2
            y2 = (y[0]+y1)/2
            z1 = (z[0]+z[1])/2
            z2 = (z[0]+z1)/2
            tmp = x[1:]
            x = [x[0],x2,x1]+list(tmp)
            x = np.array(x)
            tmp=y[1:]
            y = [y[0],y2,y1]+list(tmp)
            y = np.array(y)
            tmp=z[1:]
            z = [z[0],z2,z1]+list(tmp)
            z = np.array(z)
        
        # I think that this is where the river paths are smoothed
        tck, u = interpolate.splprep([x, y,z], s=0)
        unew = np.arange(0, 1.01, 0.05)
        out = interpolate.splev(unew, tck)
        
        lstr=[] # lstr is apparently "line string"
        dbg=[] # I think this is to verify that altitude increases continually
        for i in range(len(out[0])): # loops through each coordinate created in interpolation
            lstr.append((out[0][i],out[1][i],int(out[2][i])))
            dbg.append(int(out[2][i]))
        line = asLineString(lstr)
        
        for p in path: # for each node in the path to this particular leaf
            # I'm pretty sure this loop ensures that
            # the path to the sea is up to date
            p.rivers.append(line)


# In[ ]:


fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)
ax.imshow(shore.img)
ylim=ax.get_ylim();
xlim=ax.get_xlim();
voronoi_plot_2d(cells.vor, point_size=10, ax=ax,line_colors=['yellow']) # draws the voronoi cells?
ax.set_ylim(ylim);
ax.set_xlim(xlim);
for node in hydrology.allMouthNodes():
    for river in node.rivers:
        x = [coord[0] for coord in river.coords]
        y = [coord[1] for coord in river.coords]
        plt.plot(x,y)
for node in hydrology.allNodes():
    plt.text(node.x(),node.y(),node.id)
plt.savefig('interpolatedRiverCellNetwork.png', dpi=100)
plt.show()


# In[ ]:


n = 30
print(f'position: {hydrology.node(n).position}')
ridgeIdxes = [ ]
for ri in range(len(cells.vor.ridge_points)):
    if n in cells.vor.ridge_points[ri]:
        ridgeIdxes.append(ri)
print(f'ridges: {[cells.vor.ridge_points[ri] for ri in ridgeIdxes]}')
ridges = [cells.vor.ridge_vertices[ri] for ri in ridgeIdxes]
print(f'ridges: {ridges}')
ridges = [(cells.qs[ridge[0]],cells.qs[ridge[1]]) for ridge in ridges]
print(f'ridges: {ridges}')
ridges = [ridge for ridge in ridges if ridge[0] is not None and ridge[1] is not None]
print(f'ridges: {ridges}')
print(f'number of ridges: {len(ridges)}')


# In[ ]:


# TODO:Calculate Zees for Cees of Tees ( Elevations of center points of terrain primitives)

ba = [ ]
cb = [ ]
dc = [ ]

progressCounter = 0
numTs = len(Ts.allTs())
for t in Ts.allTs():
    a = datetime.datetime.now()
    ridges = cells.cellRidges(t.cell)
    b = datetime.datetime.now()
    ba.append((b-a).total_seconds())
    # find distance to closest sgment, and elevation at that point
    closestRdist = None
    ridgeElevation = None
    for ridge in ridges:
        if len(ridge) < 2:
            q0 = ridge[0]
            dist = distance(q0.position,t.position)
            if closestRdist is None or dist < closestRdist:
                closestRdist = dist
                ridgeElevation = q0.elevation
            continue
            
        
        q0 = ridge[0]
        q1 = ridge[1]
        dist, isToEndpoint = point_segment_distance_is_endpoint(
            t.position[0],t.position[1],
            q0.position[0],q0.position[1],
            q1.position[0],q1.position[1]
        )
        if closestRdist is not None and dist > closestRdist:
            continue
        if isToEndpoint:
            if distance(q0.position,t.position) < distance(q1.position,t.position):
                closestRdist = distance(q0.position,t.position)
                ridgeElevation = q0.elevation
            else:
                closestRdist = distance(q1.position,t.position)
                ridgeElevation = q1.elevation
        else:
            closestRdist = dist
            ridgeElevation =                 q0.elevation +                 (math.sqrt(distance(q0.position,t.position)**2 - dist**2) /                  distance(q0.position,q1.position)) *                 (q1.elevation - q0.elevation)
    
    c = datetime.datetime.now()
    cb.append((c-b).total_seconds())
    # see if the seeeeee is closer
    dist_gamma = shore.distanceToShore(t.position)
    if closestRdist is None or (dist_gamma < closestRdist):
        closestRdist = dist_gamma
        ridgeElevation = 0
    
    point = geom.Point(t.position[0],t.position[1])
    projected = None
    distancefromN = None
    node = hydrology.node(t.cell)
    if len(node.rivers) > 0:
        local_rivers = node.rivers # tries to get a line to the seeeee
        # index of the point on the interpolated river line that is closest to the Tee point
        rividx = [point.distance(x) for x in local_rivers].index(min([point.distance(x) for x in local_rivers]))
        # gets the point along the river that is the distance along the river to the point nearest to the Tee
        projected = local_rivers[rividx].interpolate(local_rivers[rividx].project(point))
        distancefromN = point.distance(local_rivers[rividx]) # distance to that point
    else: # handle cases of stub rivers
        node = hydrology.node(ridx)
        projected = geom.Point(node.x(),node.y(),node.elevation)
        distancefromN = point.distance(projected)
    
    if distancefromN==0 and closestRdist==0:
        distancefromN=1
    
    lerpedelevation = projected.z*(closestRdist/(closestRdist+distancefromN))+ridgeElevation*(distancefromN/(closestRdist+distancefromN))
    
    t.elevation = lerpedelevation
    
    d = datetime.datetime.now()
    dc.append((d-c).total_seconds())
    
    progressCounter = progressCounter + 1
    clear_output(wait=True)
    print(progressCounter," out of ",numTs)  # use display(f) if you encounter performance issues

print('printf() Profiling:')
print(f'cellRidges(): {sum(ba)/numTs}')
print(f'compute closest ridge: {sum(cb)/numTs}')
print(f'rest of the computation: {sum(dc)/numTs}')


# In[ ]:





# In[ ]:


def TerrainFunction(prePoint):
    point = [int(prePoint[0] * (shore.shape[0] / outputResolution)),int(prePoint[1] * (shore.shape[1] / outputResolution))]
    
    # if imgray[point[1]][point[0]]==0: This is why a new data model was implemented
    if not shore.isOnLand(point):
        return 0

    # Gets and computes influence and elevation values for nearby terrain primitives
    ts = Ts.query_ball_point(point,radius) # Gets all terrain primitives within a radius of the point
    if len(ts) < 1: # if there just aren't any T points around, just put it in the ocean
        return 0
    wts = [w(distance(point,t.position)) for t in ts] # "influence field" radii of those primitives
    # TODO: I think this end up getting different heights for
    hts = [ht(point,t) for t in ts]          # elevations of those primitives

    # Blends the terrain primitives
    ht_ = height_b(hts,wts) # Blends those terrain primitives
    wt_ = wts[0]            # I guess this is supposed to be the influence radius of the closest primitive?
    
    wi=wt_ # IDK why he converts these here
    hi=ht_
    
    nodeID = cells.nodeID(point)
    if nodeID is None:
        return hi
    node = hydrology.node(nodeID)
    geomp = geom.Point(point[0],point[1])     # Creates a Shapely point out of the input point
    rs = [ ]
    hrs = [ ]
    wrs = [ ]
    if len(node.rivers) > 0:
        rs  = [e for e in node.rivers if geomp.distance(e) < radius ]
        hrs = [hr(geomp,e) for e in rs]
        wrs = [w(geomp.distance(e)) for e in rs]
    else: # Sometimes there isn't a river, just a drainage point along the seeeee
        riverPoint = geom.Point(node.x(),node.y(),node.elevation)
        if geomp.distance(riverPoint) < radius:
            rs = [ geomp.distance(riverPoint) ]
            hrs = [ riverPoint.z ]
            wrs = [ w(geomp.distance(riverPoint)) ]
    
    # Height and "influence field" calculation per the last equation in Section 7
    # This is the so-called "replacement operator"
    for i in range(len(rs)): 
        hi=(1-wrs[i])*hi+wrs[i]*hrs[i] 
        wi = (1-wrs[i])*wi+wrs[i]**2

    if hi<0:
        pass
    
    return hi

def height_b(h,w): # height function of a blend node (section 7)
    try:
        ret = np.sum(np.multiply(h,w))/(np.sum(w))
        assert(ret>=0)
        assert(not np.isnan(ret)) # make sure ret is a number
        return ret
    except:
        return 0

scale = 100.0 # I think adjusting these values could be useful
octaves = 6
persistence = 0.5
lacunarity = 2.0
def ht(p,t): # Height of a terrain primitive
    return t.elevation# +pnoise2(p[0]/scale,p[1]/scale,octaves=octaves,persistence=persistence,lacunarity=lacunarity,repeatx=shore.shape[0],repeaty=shore.shape[1],base=0)*10

def hr(p,r): # Height of a river primitive?
    d=p.distance(r)
    # TODO profile based on Rosgen classification
    segma = 0.1 * min(rwidth**2,d**2) # I think this is the river profile (evidently the author can't read Greek)
    projected = r.interpolate(r.project(p))
    return projected.z+segma

def w(d): # This returns the "influence field" (section 7)
    if d <1:
        return 1;
    return (max(0,(radius+1)-d)/(((radius)+1)*d))**2


# In[ ]:


from multiprocessing import Process, Pipe, Queue
from tqdm import trange

def subroutine(conn, q):
    #print(f'Thread ID: {conn.recv()}')
    threadID = conn.recv()
    for i in range(threadID, outputResolution, 4):
        arr = np.zeros(outputResolution,dtype=np.double)
        for j in range(outputResolution):
            arr[j] = max(0,TerrainFunction((j,i)))
        try:
            q.put((i,arr.tobytes()))
        except:
            conn.close()
            exit()
        #print(f'row {i} complete')
    
    #conn.send(len(shore))
    conn.close()

imgTest = np.zeros((outputResolution,outputResolution),dtype=np.double)

#if __name__ == '__main__':
numProcs = 4
dataQueue = Queue()
pipes = []
processes = []
for p in range(numProcs):
    pipes.append(Pipe())
    processes.append(Process(target=subroutine, args=(pipes[p][1],dataQueue)))
    processes[p].start()
    pipes[p][0].send(p)
for i in trange(outputResolution):
    data = dataQueue.get()
    imgTest[data[0]] = np.frombuffer(data[1],dtype=np.double)
for p in range(numProcs):
    processes[p].join()
    pipes[p][0].close()

fig = plt.figure(figsize=(16,16))
plt.imshow(imgTest, cmap=plt.get_cmap('terrain'))


# In[ ]:





# In[ ]:


# this doesn't work because
immtt = np.array(imgTest)
normalizedImg = immtt.copy()
cv.normalize(immtt,  normalizedImg, 0, 255, cv.NORM_MINMAX)
normalizedImg = normalizedImg.astype('uint8')
cv.imwrite("taiwan-out.png",normalizedImg)


# In[ ]:


endTime = datetime.datetime.now()
print("Total time: ", endTime -startTime )


# In[ ]:





# In[ ]:


maxT.asdads


# In[ ]:




s= datetime.datetime.now()
print("start:",s)
#imgTest = np.zeros(imgray.shape,dtype='uint8')
 
print(TerrainFunction((326, 196)))
print(TerrainFunction((196, 326)))

e=datetime.datetime.now()
print("End:",e)
print("render time: ", e -s )


# In[ ]:


s= datetime.datetime.now()
print("start:",s)
#imgTest = np.zeros(imgray.shape,dtype='uint8')

imgTest = [[TerrainFunction((j,i)) for j in range(128,160)] for i in range(128,160)]

plt.imshow(imgTest)
e=datetime.datetime.now()
print("End:",e)
print("render time: ", e -s )


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


print(point_segment_distance_is_endpoint(
    6,0,
    0,0,
    5,0
))


# In[ ]:





# In[ ]:





# 
