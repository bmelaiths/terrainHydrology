import cv2 as cv
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from scipy.spatial import Voronoi
from poisson import PoissonGenerator

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
            print(f'\tCreating ridge primitive {iv} of {len(self.vor.vertices)}\r', end='')
        print()

        self.cellsRidges = { }
        for n in range(len(hydrology)):
            connectedNodes = hydrology.upstream(n)
            node = hydrology.node(n)
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
            print(f'\tCreating ridge primitive {iv} of {len(self.vor.vertices)}\r', end='')
        print()
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
    def __init__(self, hydrology, cells, num_points):
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
            
            print(f'\tPrimitives created: {n} of {len(hydrology)} \r', end='')  # use display(f) if you encounter performance issues
        print()

        allpoints_list = [[t.position[0],t.position[1]] for t in self.allTs()]
        allpoints_nd = np.array(allpoints_list)
        self.apkd = cKDTree(allpoints_nd)
    def allTs(self):
        return self.tList.copy()
    def cellTs(self, cell):
        return self.Ts[cell].copy()
    def query_ball_point(self, loc, radius):
        return [self.tList[i] for i in self.apkd.query_ball_point(loc,radius)]