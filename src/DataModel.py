import cv2 as cv
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from scipy.spatial import Voronoi
from poisson import PoissonGenerator
from PIL import Image
import struct

import typing

class RasterData:
    """A simple abstraction of raster data based on an image.

    Simply input a path to an image and a resolution, and this class will
    allow you to access the data for each point on the surface.

    :param inputFileName: This should be a path to the image
    :type inputFileName: str
    :param resolution: The resolution of the input image in meters per pixel
    :type resolution: float
    
    The image should be
    readable by PIL. The image should represent raster data that covers
    the entire area that will be accessed by the algorithm. The values
    are represented by the grayscale value of each pixel. The image's
    actual color model does not model, as it will be converted to
    grayscale in the constructor.
    """
    def __init__(self, inputFileName: str, resolution: float):
        self.raster = Image.open(inputFileName)
        self.xSize = self.raster.size[0]
        self.ySize = self.raster.size[1]
        self.raster = self.raster.convert('L')
        self.raster = self.raster.load()
        self.resolution = resolution
    def __getitem__(self, loc: typing.Tuple[float,float]) -> float:
        """Gets the value for a particular location

        :param loc: The location to interrogate
        :type loc: tuple[float,float]
        :return: The value of the data at loc
        :rtype: float
        """
        return self.raster[int(loc[0]/self.resolution),int(loc[1]/self.resolution)]
    def toBinary(self):
        binary = None
        for y in range(self.ySize):
            # print(f'Working on row {y}')
            row = None
            for x in range(self.xSize):
                if row is not None:
                    row = row + struct.pack('!f', self.raster[x,y])
                else:
                    row = struct.pack('!f', self.raster[x,y])
            if binary is not None:
                binary = binary + row
            else:
                binary = row
        return binary

class ShoreModel:
    """This class represents the land area.
    
    It is indexable as an array of points. These points
    represent the line that demarcates the boundary between
    land and sea. This class also has useful functions.

    :param inputFileName: The path to the image that defines the shoreline
    :type inputFileName: str
    :param resolution: The resolution of the input image in meters per pixel
    :type resolution: float

    :cvar realShape: The spatial dimensions of the area that the gamma image covers, in meters
    :vartype realShape: numpy.ndarray[float,float]

    .. note::
       Shape variables are all in order y,x.
    """
    def __init__(self, inputFileName: str, resolution: float):
        """Constructor
        """

        self.resolution = resolution

        self.img = cv.imread(inputFileName)
        
        self.imgray = cv.cvtColor(self.img, cv.COLOR_BGR2GRAY) # a black-and-white version of the input image
        self.rasterShape = self.imgray.shape
        self.realShape = (self.imgray.shape[0] * self.resolution, self.imgray.shape[1] * self.resolution)
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
    def distanceToShore(self, loc) -> float:
        """Gets the distance between a point and the shore

        The distance is in meters.

        :param loc: The location to test
        :type loc: a tuple of two floats
        :return: The distance between `loc` and the shore in meters
        :rtype: float
        """

        #    for some reason this method is      y, x
        return cv.pointPolygonTest(self.contour,(loc[1]/self.resolution,loc[0]/self.resolution),True) * self.resolution
    def isOnLand(self, loc) -> bool:
        """Determines whether or not a point is on land or not

        :param loc: The location to test
        :type loc: a tuple of two floats
        :return: True if the point is on land, False if otherwise
        :rtype: bool
        """
        
        if 0 <= loc[0] < self.realShape[1] and 0 <= loc[1] < self.realShape[0]:
            return self.imgray[int(loc[1]/self.resolution)][int(loc[0]/self.resolution)] == 255 # != 0
        else:
            return False
    def __getitem__(self, index: int):
        """Gets a point on the shore by index

        The shore can be thought of as an array of points. These points
        demarcate the boundary between land and sea, and can be indexed

        :param index: The index
        :type index: int
        :return: The index-th coordinate of the shoreline
        :rtype: tuple of two floats, representing the coordinates in meters
        """

        # TODO ensure that points returned are x,y
        return (self.contour[index][1]*self.resolution,self.contour[index][0]*self.resolution)
    def __len__(self):
        """The number of points that make up the shore

        :return: The number of points that make up the shore
        :rtype: int
        """
        return len(self.contour)

class HydroPrimitive:
    """Represents a certain stretch of river

    A HydroPrimitive is instantiated with the ``id``, ``position``,
    ``elevation``, ``priority``, and ``parent`` attributes, as
    applicable. The other attributes are computed later.

    :cvar id: The ID of this node. See :class:`HydrologyNetwork` for this value's significance
    :vartype id: int
    :cvar position: The location of the node in meters
    :vartype position: tuple[float,float]
    :cvar elevation: The node's elevation in meters
    :vartype elevation: float
    :cvar priority: The priority of the node. See :class:`HydrologyNetwork` for this value's significance
    :vartype priority: int
    :cvar parent: The **ID** of the parent node, or None if this node is a river mouth
    :vartype parent: int | None
    :cvar contourIndex: If this node is on the coast, this is the index in :class:`Shore` that is closest to this node
    :vartype contourIndex: int
    :cvar rivers: A :class:`LineString` representing the river's actual path. It only flows to a node where it joins with a larger river
    :vartype rivers: LineString
    :cvar localWatershed: The (rough) area of this cell
    :vartype localWatershed: float
    :cvar inheritedWatershed: The area of this cell, plus all the areas of all the cells that drain into it
    :vartype inheritedWatershed: float
    :cvar flow: The flow rate of water draining out of this cell (including flow from ancestor cells) in cubic meters per second
    :vartype flow: float
    """
    def __init__(self, id: int, loc: typing.Tuple[float,float], elevation: float, priority: int, parent: int):
        self.id = id
        self.position = loc
        self.elevation = elevation
        self.priority = priority
        self.parent = parent
        self.rivers = [ ]
    def x(self) -> float:
        """Gets the x location of this node

        :return: The x location of this node
        :rtype: float
        """
        return self.position[0]
    def y(self) -> float:
        """Gets the y position of this node

        :return: The y location of this node
        :rtype: float
        """
        return self.position[1]

class HydrologyNetwork:
    """This class represents the network of rivers that flow over the land

    A HydrologyNetwork is basically a forest of trees, with each tree
    representing a river that merges and drains into the ocean through a
    single mouth node.

    A HydrologyNetwork is empty when it is instantiated. The network is
    built incrementally using
    :func:`addNode()<DataModel.HydrologyNetwork.addNode>`. Edges connect
    all nodes.

    It should be noted that each node is associated with an integer ID.
    This ID is strictly the order that each node was added in, starting at
    0. Nodes cannot be removed, thus ``range(len(hydrology))`` will
    iterate over all the nodes in the network.

    .. note::
       Some methods return references to the actual HydroPrimitives that constitute
       the graph. Others return integers that refer to them. Always refer to the
       return type when querying an instance of this class.

    Internally, the data is held in a :class:`networkx DiGraph<networkx.DiGraph>`. A
    :class:`cKDTree<scipy.spatial.cKDTree>` is used for lookup by area.
    """
    def __init__(self, stream=None):
        self.nodeCounter = 0
        self.graph = nx.DiGraph()
        self.mouthNodes = []

        if stream is None:
            return

        buffer = stream.read(8)
        numberNodes = struct.unpack('!Q', buffer)[0]

        for i in range(numberNodes):
            buffer = stream.read(8)
            nodeID = struct.unpack('!Q', buffer)[0]

            buffer = stream.read(8)
            parentID = struct.unpack('!Q', buffer)[0]

            buffer = stream.read(1)
            numChildren = struct.unpack('!B', buffer)[0]

            for chiild in range(numChildren):
                buffer = stream.read(8)
                childID = struct.unpack('!Q', buffer)[0]

            buffer = stream.read(4)
            locX = struct.unpack('!f', buffer)[0]

            buffer = stream.read(4)
            locY = struct.unpack('!f', buffer)[0]

            buffer = stream.read(4)
            elevation = struct.unpack('!f', buffer)[0]

            self.addNode(
                (locX, locY),
                elevation,
                0, # Does priority actually matter elsewhere?
                parent = self.node(parentID) if parentID != nodeID else None # the parent ID of a mouth node will be its own ID
            )
    def addNode(self, loc: typing.Tuple[float,float], elevation: float, priority: int, contourIndex: int=None, parent: int=None) -> HydroPrimitive:
        """Creates and adds a HydrologyPrimitive to the network

        :param loc: The location of the new node
        :type loc: tuple[float,float]
        :param elevation: The elevation of the node in meters
        :type elevation: float
        :param priority: The priority of the node (for graph expansion)
        :type priority: int
        :param contourIndex: The index of the node's location on the shore (corresponds to :class:`ShoreModel[]<DataModel.ShoreModel>`)
        :type contourIndex: int, optional
        :param parent: The ID of the parent node, if applicable
        :type parent: int, optional

        The priority of the node affects its selection, as described in §4.2.1
        of Genevaux et al, within the graph expansion algorithm.

        The contourIndex should be used when the node represents the mouth of
        a river. In this case, contourIndex should be the index such that,
        when passed to ``ShoreModel[contourIndex]``, yields the position that
        most closely corresponds to the node. (Internally, this is used to
        figure out what the 'direction' of the river should be in the graph
        expansion.)

        :return: The node created
        :rtype: HydroPrimitive
        """
        node = HydroPrimitive(self.nodeCounter, loc, elevation, priority, parent)
        if parent is None or contourIndex is not None:
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
        
        # Classify the new leaf
        node.priority = 1

        # Classify priorities of affected nodes
        classifyNode = node.parent
        while True:
            if classifyNode is None:
                break
            children = self.upstream(classifyNode.id)
            maxNumber = max([child.priority for child in children])
            numMax = len([child.priority for child in children if child.priority == maxNumber])
            if numMax > 1 and classifyNode.priority < maxNumber + 1:
                # if there is more than one child with the maximum number,
                # and the parent isn't already set for it, then change it
                classifyNode.priority = maxNumber + 1
            elif classifyNode.priority < maxNumber:
                # if the parent isn't already set for the maximum number,
                # change it
                classifyNode.priority = maxNumber
            else:
                # if the parent does not need to be changed at all, then
                # none of its ancestors do, and the graph is fully adjsuted
                break
            classifyNode = classifyNode.parent

        return node
    def query_ball_point(self, loc: typing.Tuple[float,float], radius: float) -> typing.List[int]:
        """Gets all nodes that are within a certain distance of a location

        :param loc: The location to test
        :type loc: tuple[float,float]
        :param radius: The radius to search in
        :type radius: float
        :return: The IDs of all nodes that are within ``radius`` of ``loc``
        :rtype: list[int]
        """
        return self.graphkd.query_ball_point(loc,radius)
    def edgesWithinRadius(self, loc: typing.Tuple[float,float], radius: float) -> typing.List[typing.Tuple[HydroPrimitive,HydroPrimitive]]:
        """Gets all *edges* that are within a certain distance of a location

        :param loc: The location to test
        :type loc: tuple[float,float]
        :param radius: The radius to search in
        :type radius: float
        :return: Each tuple represents both ends of the edge
        :rtype: list[tuple[HydroPrimitive,HydroPrimitive]]
        """
        nodesToCheck = self.graphkd.query_ball_point(loc,radius)
        edges = [ self.graph.out_edges(n) for n in nodesToCheck ]
        edges = [item for edge in edges for item in edge]
        return [(self.graph.nodes[e[0]]['primitive'],self.graph.nodes[e[1]]['primitive']) for e in edges]
    def downstream(self, node: int) -> HydroPrimitive:
        """Gets the node that this node flows into

        :param node: The ID of the node whose parent you wish to identify
        :type node: int
        :return: The node that this node flows into
        :rtype: HydroPrimitive
        """
        parent = list(self.graph.predecessors(node))
        if len(parent) > 0:
            return self.graph.nodes[parent[0]]['primitive']
        else:
            return None
    def upstream(self, node: int) -> typing.List[HydroPrimitive]:
        """Gets all the nodes that flow into this node

        :param node: The ID of the node whose ancestors you wish to query
        :type node: int
        :return: A list of nodes that are upstream of this one
        :rtype: list[HydroPrimitive]
        """
        return [self.graph.nodes[n]['primitive'] for n in self.graph.successors(node)]
    def adjacentNodes(self, node: int) -> typing.List[HydroPrimitive]:
        """Basically just a concatenation of :func:`downstream()<DataModel.HydrologyNetwork.downstream>` and :func:`upstream()<DataModel.HydrologyNetwork.upstream>`

        :param node: The ID of the node whose adjacent nodes you wish to identify
        :type node: int
        :return: A list of all adjacent nodes, upstream and downstream
        :rtype: list[HydroPrimitive]
        """
        downstream = self.downstream(node)
        upstream   = self.upstream(node)
        if downstream is None:
            downstream = [ ]
        return upstream + downstream
    def allUpstream(self, node: int) -> typing.List[HydroPrimitive]:
        """Gets *all* nodes that are upstream of this one---all the way out to the leaf nodes

        :param node: The node whose upstream nodes you wish to identify
        :type node: int
        :return: All nodes that are upstream of this one
        :rtype: list[HydroPrimitive]
        """
        return [self.graph.nodes[n]['primitive'] for n in nx.descendants(self.graph, node)]
    def allNodes(self) -> typing.List[HydroPrimitive]:
        """All nodes in the graph

        This can be used if it is more convenient than ``range(len(hydrology))``

        :return: All nodes in the graph
        :rtype: list[HydroPrimitive]
        """
        return [self.graph.nodes[node]['primitive'] for node in list(self.graph.nodes)]
    def allEdges(self) -> typing.List[typing.Tuple[HydroPrimitive,HydroPrimitive]]:
        """Gets all edges in the graph

        :return: Every edge in the graph
        :rtype: list[tuple[HydroPrimitive,Hydroprimitive]]
        """
        return [(self.graph.nodes[u]['primitive'],self.graph.nodes[v]['primitive']) for u,v in self.graph.edges]
    def allMouthNodes(self) -> typing.List[HydroPrimitive]:
        """Gets all the mouth nodes (those that drain to the sea)

        :return: Every mouth node
        :rtype: list[HydroPrimitive]
        """
        return [self.graph.nodes[id]['primitive'] for id in self.mouthNodes]
    def allLeaves(self, node) -> typing.List[HydroPrimitive]:
        """Gets all leaf nodes that antecede this node

        :param node: The node whose leaf ancestors you wish to identify
        :type node: list[HydroPrimitive]
        """
        ids = [s for s in nx.descendants(self.graph,node) if len(self.graph.out_edges(s))==0]
        return [self.graph.nodes[id]['primitive'] for id in ids]
    def node(self, node: int) -> HydroPrimitive:
        """Gets a reference to the node that corresponds to a given ID

        :param node: The ID of the node you wish to query
        :type node: int
        :return: A reference to the HydroPrimitive that corresponds to the ID
        :rtype: HydroPrimitive
        """
        return self.graph.nodes[node]['primitive']
    def dfsPostorderNodes(self) -> typing.List[HydroPrimitive]:
        """Returns a list of all nodes in the network in a *depth-first, postorder* order

        See `Depth-first search <https://en.wikipedia.org/wiki/Depth-first_search/>`_.

        :return: All nodes in the network
        :rtype: list[HydroPrimitive]
        """
        ids = list(nx.dfs_postorder_nodes(self.graph))
        return [self.graph.nodes[id]['primitive'] for id in ids]
    def pathToNode(self, origin: int, destination: int) -> typing.List[HydroPrimitive]:
        """Returns the the path between any two nodes (but they should be in the same river system)

        :param origin: One end of the path
        :type origin: int
        :param destination: The other end of the path
        :type destination: int
        :return: References to the HydroPrimitives that make up the path
        :rtype: list[HydroPrimitive]
        """
        return [self.graph.nodes[n]['primitive'] for n in nx.shortest_path(self.graph,origin,destination)]
    def __len__(self) -> int:
        """Returns the number of nodes in the forest

        :return: The number of nodes in the forest
        :rtype: int
        """
        return len(self.graph.nodes)

def openCVFillPolyArray(points: typing.List[typing.Tuple[float,float]]) -> typing.List[np.ndarray]:
    """Formats a list of points into a format that OpenCV will understand

    :param points: The points for format
    :type points: list[tuple[float,float]]
    :return: Returns the points in a format that some OpenCV methods can use
    :rtype: list[np.array[[float,float]]]
    """
    return [ np.array( [ [int(p[0]),int(p[1])] for p in points ] ) ]

class Q:
    """Represents a ridge point

    Ridge points are created with a ``position``, ``nodes``, and ``vorIndex``.
    The elevation is computed later.

    :cvar position: Location
    :vartype position: tuple[float,float]
    :cvar nodes: A list of the IDs of each cell that this vertex borders
    :vartype nodes: list[int]
    :cvar vorIndex: The internal index of the voronoi vertex that this Q represents
    :vartype vorIndex: int
    :cvar elevation: The elevation of this point
    :vartype elevation: float
    """
    def __init__(self, position, nodes, iv):
        self.position = position
        self.nodes = nodes
        self.vorIndex = iv # the index of the voronoi vertex this represents in vor.vertices

class TerrainHoneycomb:
    """This class partitions the land into cells around the river nodes

    There is a cell around each river node. Every cell is a polygon. Each
    edge of any polygon is either transected by the flow of a river, or
    forms a mountainous ridge between two rivers.

    :param shore: The ShoreModel for the land area
    :type shore: ShoreModel
    :param hydrology: The filled-out HydrologyNetwork for the land area
    :type hydrology: HydrologyNetwork
    :param resolution: The resolution of the underlying rasters in meters per pixel
    :type resolution: float
    :param dryRun: Indicate that this object is being used for a dry run, thus ridge data will not be needed
    :type dryRun: bool

    .. note::
       ``resolution`` should be the same that was passed to the ShoreModel.

    Internally, this class encapsulates a :class:`scipy.spatial.Voronoi`
    instance and a couple of dictionaries to classify ridges and other
    edges.
    """
    def __init__(self, shore: ShoreModel, hydrology: HydrologyNetwork, resolution: float, dryRun: bool):
        self.resolution = resolution

        self.shore     = shore
        self.hydrology = hydrology
        
        points = [node.position for node in hydrology.allNodes()]
        points.append((0,0))
        points.append((0,shore.realShape[1]))
        points.append((shore.realShape[0],0))
        points.append((shore.realShape[0],shore.realShape[1]))
        
        self.vor = Voronoi(points,qhull_options='Qbb Qc Qz Qx')
        
        self.imgvoronoi = np.zeros(shore.rasterShape, dtype=np.uint16)
        for n in range(len(hydrology)):
            if self.vor_region_id(n) == -1:
                continue
            positions = [(int(p[0]/self.resolution),int(p[1]/self.resolution)) for p in self.ridgePositions(n)]
            cv.fillPoly(
                self.imgvoronoi,
                openCVFillPolyArray(positions),
                np.int16(self.vor_region_id(n)+1).item()
            )
        if not dryRun:
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
            self.cellsDownstreamRidges = { }
            # Classify all ridges
            for ri in range(len(self.vor.ridge_vertices)):
                for n in self.vor.ridge_points[ri]: # Each ridge separates exactly two nodes
                    if n >= len(self.hydrology):
                        continue # Apparently there are nodes that don't exist; skip these
                    node = hydrology.node(n)
                    otherNode = self.vor.ridge_points[ri][self.vor.ridge_points[ri] != n][0]
                    # if this ridge is the outflow ridge for this node, mark it as such and move on
                    if node.parent is not None and node.parent.id == otherNode:
                        # this ridge is the outflow ridge for this node
                        v1 = self.vor.vertices[self.vor.ridge_vertices[ri][0]]
                        v2 = self.vor.vertices[self.vor.ridge_vertices[ri][1]]
                        if not self.shore.isOnLand(v1) or not self.shore.isOnLand(v2):
                            # If one or both vertices is not on land, then don't bother
                            # trying to make the river flow through the ridge neatly
                            self.cellsDownstreamRidges[n] = None
                        else:
                            self.cellsDownstreamRidges[n] = (v1, v2)
                        break # the outflow ridge of this node is an inflow ridge for the other one
                    if otherNode in [nd.id for nd in hydrology.upstream(n)]:
                        continue # this is an inflow ridge, so it need not be considered further
                    # the ridge at ri is not transected by a river
                    if n not in self.cellsRidges:
                        self.cellsRidges[n] = [ ]
                    vertex0 = self.qs[self.vor.ridge_vertices[ri][0]]
                    vertex1 = self.qs[self.vor.ridge_vertices[ri][1]]
                    if vertex0 is None or vertex1 is None:
                        continue
                    self.cellsRidges[n].append((vertex0,vertex1))
                print(f'\tClassifying ridge {ri} of {len(self.vor.ridge_vertices)}\r', end='')
            print()

            # Add vertices that are not attached to ridges
            for n in range(len(self.hydrology)):
                verts = self.vor.regions[self.vor_region_id(n)].copy()
                if n not in self.cellsRidges:
                    self.cellsRidges[n] = [ ]
                # Eliminate vertices that are attached to ridges
                for ridge in self.cellsRidges[n]:
                    for q in ridge:
                        if q is not None and q.vorIndex in verts:
                            verts.remove(q.vorIndex)
                for v in verts:
                    if not self.shore.isOnLand(self.vor.vertices[v]):
                        continue
                    self.cellsRidges[n].append((self.qs[v],))
                print(f'\tFinding unaffiliated vertices for node {n} of {len(hydrology)}\r', end='')
            print()
    def vor_region_id(self, node: int) -> int:
        """Returns the index of the *voronoi region*

        This method is useful because the index of the voronoi region is not the same
        as the ID of the cell it is based on.

        :param node: The ID of the node/cell
        :type node: int

        :return: The ID of the *voronoi region* (or, cell) associated with that node
        :rtype: int

        .. note::
           This method is for internal use. If you are using it from outside this
           class, your approach is definitely breaking a key design principle.
        """
        return self.vor.point_region[node]
    def ridgePositions(self, node: int) -> typing.List[typing.Tuple[float,float]]:
        """Gets the position of each vertex that makes up the cell

        :param node: The ID of the node whose cell you wish to query
        :type node: int
        :return: A list of points that correspond to the vertices that make the cell
        :rtype: list[tuple[float,float]]
        """
        ridges = self.vor.regions[self.vor_region_id(node)] # the indices of the vertex boundaries
        return [self.vor.vertices[x] for x in ridges if x != -1] # positions of all the vertices
    def cellArea(self, loc: typing.Tuple[float,float]) -> float:
        """Calculates the (rough) area of the cell that a location is in

        This is a rough estimate. It is derived by counting the number of pixels
        that correspond to the given cell in an internal raster that is based on
        the rasters within ShoreModel. Thus, accuracy is proportional to the
        resolution of the user-supplied rasters. But it is good enough for what
        it is used for.

        :param loc: Any location within the cell you wish to query
        :type loc: tuple[float,float]
        """
        return np.count_nonzero(self.imgvoronoi == self.imgvoronoi[int(loc[1]/self.resolution)][int(loc[0]/self.resolution)]) * self.resolution**2
    def cellQs(self, node: int) -> typing.List[Q]:
        """Returns all the Qs binding the cell that corresponds to the given node

        :param node: The node you wish to query
        :type node: int
        :return: List of Q instances
        :rtype: list[Q]
        """
        return [self.qs[vorIdx] for vorIdx in self.vor.regions[self.vor_region_id(node)] if self.qs[vorIdx] is not None]
    def allQs(self) -> typing.List[Q]:
        """Simply returns all Qs of the land

        :return: All Qs of the land
        :rtype: list[Q]
        """
        return self.qs.copy()
    def boundingBox(self, n: int) -> typing.Tuple[float,float,float,float]:
        """Returns the measurements for a bounding box that would contain a cell

        This method had a rather specific application. You'll probably never use it.

        :param n: The ID of the node/cell you wish to get a bounding box for
        :type n: int
        :return: A tuple indicating the lower X, upper X, lower Y, and upper Y, respectively, in meters
        :rtype: tuple[float,float,float,float]
        """
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
        return (xllim * self.resolution, xulim * self.resolution, yllim * self.resolution, yulim * self.resolution)
    def isInCell(self, p: typing.Tuple[float,float], n: int) -> bool:
        """Determines if a point is within a given cell

        Like :func:`boundingBox`, it is a rough estimate based on a raster
        that is derived from the rasters in :class:`ShoreModel`

        :param p: The point you wish to test (measurements in meters)
        :type p: tuple[float,float]
        :param n: The ID of the cell you wish to test
        :type n: int
        :return: True of point ``p`` is in the cell that corresponds to ``n``
        :rtype: bool
        """
        return self.shore.isOnLand(p) and self.imgvoronoi[int(p[1]/self.resolution)][int(p[0]/self.resolution)]==self.vor.point_region[n]+1
    def cellRidges(self, n: int) -> typing.List[tuple]:
        """Returns the mountain ridges of a cell

        That is, this method returns the edges that are not transected by the
        flow of a river---in or out.

        .. note::
           This method returns a list of tuples, each of which contain exactly
           1 **OR** 2 :class:`Q` instances. If the tuple contains 2 Qs, then
           the Qs constitute a mountain ridge. If the tuple contains 1 Q, then
           it represents a mountain (or hill) that is not connected to a
           ridge.

        :param n: The ID of the cell that you wish to query
        :type n: int
        :return: A list of tuples, each contain exactly 1 or 2 Qs.
        :rtype: list[tuple]
        """
        return self.cellsRidges[n]
    def cellOutflowRidge(self, n: int) -> typing.Optional[typing.Tuple[int,int]]:
        """Returns the ridge through which the river flows out of the cell

        .. note::
           A ridge will only be returned if *both* vertices are on land. Otherwise
           ``None`` will be returned.

           Moreover, the method returns the *IDs* of the vertices, not the Qs
           themselves.

        .. todo::
           This method is implemented poorly, and should probably be reworked in
           the future.
        
        :return: The IDs of the vertices that define the outflow ridge, if both are on land
        :rtype: tuple[int,int] | None
        """
        return self.cellsDownstreamRidges[n]
    def nodeID(self, point: typing.Tuple[float,float]) -> int:
        """Returns the id of the node/cell in which the point is located

        Like :func:`cellArea`, this is not entirely precise, as it is derived
        from a raster based on those in ShoreModel.

        :param point: The point you wish to test
        :type point: tuple[float,float]
        :return: The ID of a node/cell
        :rtype: int
        """
        id = list(self.vor.point_region).index(self.imgvoronoi[int(point[1]/self.resolution)][int(point[0]/self.resolution)]-1)
        return id if id != -1 else None

class T:
    """Terrain primitive

    Terrain primitives are created with a ``position`` and ``cell``
    attribute. The elevation is computed separately.

    :cvar position: The position of the primitive
    :vartype position: tuple[float,float]
    :cvar cell: The ID of the cell in which this primitive is situated
    :vartype cell: int
    :cvar elevation: The elevation of the primitive
    :vartype elevation: float
    """
    def __init__(self, position, cell):
        self.position = position
        self.cell = cell

class Terrain:
    """Holds and organizes the terrain primitives (:class:`T`)

    :param hydrology: The HydrologyNetwork of the land
    :type hydrology: HydrologyNetwork
    :param cells: The TerrainHoneycomb of the land
    :type cells: TerrainHoneycomb
    :param num_points: (Roughly) the number of points in each cell
    :type num_points: int
    """
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
    def allTs(self) -> typing.List[T]:
        """Simply returns all the terrain primitives

        :return: A list of all the terrain primitives
        :rtype: list[T]
        """
        return self.tList.copy()
    def cellTs(self, cell: int) -> typing.List[T]:
        """Gets the terrain primitives within a given cell

        :param cell: The ID of the cell you wish to query
        :type cell: int
        :return: A list of the terrain primitives in the cell
        :rtype: list[T]
        """
        return self.Ts[cell].copy()
    def query_ball_point(self, loc: typing.Tuple[float,float], radius: float) -> typing.List[T]:
        """Gets all terrain primitives within a given radius of a given location

        :param loc: The location you wish to test
        :type loc: tuple[float,float]
        :param radius: The search radius
        :type radius: float
        :return: A list of the primitives within that area
        :rtype: list[T]
        """
        return [self.tList[i] for i in self.apkd.query_ball_point(loc,radius)]