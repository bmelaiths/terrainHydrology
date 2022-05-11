from re import I
import cv2 as cv
from networkx.algorithms.operators import binary
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from scipy.spatial import Voronoi
from poisson import PoissonGenerator
from PIL import Image
import shapely.geometry as geom
import struct
import math
from tqdm import trange
import datetime

import typing

import Math

def toImageCoordinates(loc: typing.Tuple[float,float], imgSize: typing.Tuple[float,float], resolution: float) -> typing.Tuple[float,float]:
    x = loc[0]
    x /= resolution
    x += imgSize[0] * 0.5

    y = loc[1]
    y /= resolution
    y = imgSize[1] * 0.5 - y

    return (x,y)

def fromImageCoordinates(loc: typing.Tuple[float,float], imgSize: typing.Tuple[float,float], resolution: float) -> typing.Tuple[float,float]:
    x = loc[0]
    x -= imgSize[0] * 0.5
    x *= resolution

    y = loc[1]
    y = imgSize[1] * 0.5 - y
    y *= resolution

    return (x,y)

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
        loc = toImageCoordinates(loc, (self.xSize,self.ySize), self.resolution)

        return self.raster[int(loc[0]), int(loc[1])]
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

    :param gammaFileName: The path to the image that defines the shoreline
    :type gammaFileName: str
    :param resolution: The resolution of the input image in meters per pixel
    :type resolution: float

    :cvar realShape: The spatial dimensions of the area that the gamma image covers, in meters
    :vartype realShape: numpy.ndarray[float,float]

    .. note::
       Shape variables are all in order y,x.
    """
    def __init__(self, resolution: float, gammaFileName: str=None, binaryFile: typing.IO=None):
        """Constructor
        """

        if gammaFileName is not None:
            self._initFromGammaImage(resolution, gammaFileName)
        elif binaryFile is not None:
            self._initFromBinaryFile(resolution, binaryFile)
        else:
            raise ValueError('You must either create a shore from an image, or reconstitute it from a binary file')
    def _initFromGammaImage(self, resolution, inputFileName):
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
    def _initFromBinaryFile(self, resolution, binaryFileName):
        self.resolution = resolution
        self.rasterShape = (
            struct.unpack('!Q', binaryFileName.read(struct.calcsize('!Q')))[0],
            struct.unpack('!Q', binaryFileName.read(struct.calcsize('!Q')))[0]
        )
        self.imgray = np.zeros(self.rasterShape)
        for d0 in range(self.rasterShape[0]):
            for d1 in range(self.rasterShape[1]):
                self.imgray[d0][d1] = struct.unpack('!B', binaryFileName.read(struct.calcsize('!B')))[0]
        self.realShape = (self.imgray.shape[0] * self.resolution, self.imgray.shape[1] * self.resolution)
        contourLength = struct.unpack('!Q', binaryFileName.read(struct.calcsize('!Q')))[0]
        self.contour = [ ]
        for i in range(contourLength):
            self.contour.append((
                struct.unpack('!Q', binaryFileName.read(struct.calcsize('!Q')))[0],
                struct.unpack('!Q', binaryFileName.read(struct.calcsize('!Q')))[0]
            ))
    def distanceToShore(self, loc) -> float:
        """Gets the distance between a point and the shore

        The distance is in meters.

        :param loc: The location to test
        :type loc: a tuple of two floats
        :return: The distance between `loc` and the shore in meters
        :rtype: float
        """
        loc = toImageCoordinates(loc, self.imgray.shape, self.resolution)

        #    for some reason this method is      y, x
        return cv.pointPolygonTest(self.contour,(loc[1],loc[0]),True) * self.resolution
    def isOnLand(self, loc) -> bool:
        """Determines whether or not a point is on land or not

        :param loc: The location to test
        :type loc: a tuple of two floats
        :return: True if the point is on land, False if otherwise
        :rtype: bool
        """
        loc = toImageCoordinates(loc, self.imgray.shape, self.resolution)
        
        if 0 <= loc[0] < self.imgray.shape[1] and 0 <= loc[1] < self.imgray.shape[0]:
            return self.imgray[int(loc[1])][int(loc[0])] == 255 # != 0
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
        # return (self.contour[index][1]*self.resolution,self.contour[index][0]*self.resolution)

        return fromImageCoordinates((self.contour[index][1],self.contour[index][0]), self.imgray.shape, self.resolution)
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
        self.inheritedWatershed = 0
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
    def __init__(self, stream=None, binaryFile=None):
        self.nodeCounter = 0
        self.graph = nx.DiGraph()
        self.mouthNodes = []

        if stream is not None:
            self._initFromStream(stream)
        elif binaryFile is not None:
            self._initFromBinary(binaryFile)
    def _initFromStream(self, pipe):
        buffer = pipe.read(8)
        numberNodes = struct.unpack('!Q', buffer)[0]

        allpoints_list = []

        for i in range(numberNodes):
            buffer = pipe.read(8)
            nodeID = struct.unpack('!Q', buffer)[0]

            buffer = pipe.read(8)
            parent = struct.unpack('!Q', buffer)[0]
            buffer = pipe.read(8)
            contourIndex = struct.unpack('!Q', buffer)[0]
            if parent == nodeID:
                parent = None
            else:
                parent = self.node(parent)
                contourIndex = None

            buffer = pipe.read(1)
            numChildren = struct.unpack('!B', buffer)[0]

            for chiild in range(numChildren):
                buffer = pipe.read(8)
                childID = struct.unpack('!Q', buffer)[0]

            buffer = pipe.read(4)
            locX = struct.unpack('!f', buffer)[0]

            buffer = pipe.read(4)
            locY = struct.unpack('!f', buffer)[0]

            buffer = pipe.read(4)
            elevation = struct.unpack('!f', buffer)[0]

            allpoints_list.append( (locX, locY) )

            node = HydroPrimitive(self.nodeCounter, (locX,locY), elevation, 0, parent)
            if contourIndex is not None:
                node.contourIndex = contourIndex

            self.graph.add_node(
                self.nodeCounter,
                primitive=node
            )
            if parent is None:
                self.mouthNodes.append(self.nodeCounter)
            else:
                self.graph.add_edge(parent.id, self.nodeCounter)

            self.nodeCounter += 1

        self.graphkd = cKDTree(allpoints_list)
    def _initFromBinary(self, file):
        buffer = file.read(8)
        numberNodes = struct.unpack('!Q', buffer)[0]

        allpoints_list = []

        for i in range(numberNodes):
            buffer = file.read(struct.calcsize('!I'))
            nodeID = struct.unpack('!I', buffer)[0]

            buffer = file.read(4)
            locX = struct.unpack('!f', buffer)[0]

            buffer = file.read(4)
            locY = struct.unpack('!f', buffer)[0]

            buffer = file.read(4)
            elevation = struct.unpack('!f', buffer)[0]

            buffer = file.read(struct.calcsize('!I'))
            parent = struct.unpack('!I', buffer)[0]
            if parent == nodeID:
                parent = None
            else:
                parent = self.node(parent)

            buffer = file.read(struct.calcsize('!I'))
            contourIndex = struct.unpack('!I', buffer)[0]

            rivers = [ ]
            buffer = file.read(struct.calcsize('!B'))
            numRivers = struct.unpack('!B', buffer)[0]
            for i in range(numRivers):
                points = [ ]
                buffer = file.read(struct.calcsize('!H'))
                numPoints = struct.unpack('!H', buffer)[0]
                for j in range(numPoints):
                    buffer = file.read(struct.calcsize('!f'))
                    riverX = struct.unpack('!f', buffer)[0]
                    buffer = file.read(struct.calcsize('!f'))
                    riverY = struct.unpack('!f', buffer)[0]
                    buffer = file.read(struct.calcsize('!f'))
                    riverZ = struct.unpack('!f', buffer)[0]
                    points.append((riverX,riverY,riverZ))
                rivers.append(geom.LineString(points))
            
            buffer = file.read(struct.calcsize('!f'))
            localWatershed = struct.unpack('!f', buffer)[0]

            buffer = file.read(struct.calcsize('!f'))
            inheritedWatershed = struct.unpack('!f', buffer)[0]

            buffer = file.read(struct.calcsize('!f'))
            flow = struct.unpack('!f', buffer)[0]

            allpoints_list.append( (locX, locY) )

            node = HydroPrimitive(self.nodeCounter, (locX,locY), elevation, 0, parent)

            if parent is not None:
                node.contourIndex = contourIndex
            node.rivers = rivers
            node.localWatershed = localWatershed
            node.inheritedWatershed = inheritedWatershed
            node.flow = flow

            self.graph.add_node(
                self.nodeCounter,
                primitive=node
            )
            if parent is None:
                self.mouthNodes.append(self.nodeCounter)
            else:
                self.graph.add_edge(parent.id, self.nodeCounter)

            self.nodeCounter += 1

        self.graphkd = cKDTree(allpoints_list)
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
        self.elevation = 0

class TerrainHoneycomb:
    """This class partitions the land into cells around the river nodes

    There is a cell around each river node. Every cell is a polygon. Each
    edge of any polygon is either transected by the flow of a river, or
    forms a mountainous ridge between two rivers.

    :param shore: The ShoreModel for the land area
    :type shore: ShoreModel
    :param hydrology: The filled-out HydrologyNetwork for the land area
    :type hydrology: HydrologyNetwork
    :param edgeLength: The edge length in the simulation
    :type edgeLength: float
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
    def __init__(self, shore: ShoreModel=None, hydrology: HydrologyNetwork=None, resolution: float=None, edgeLength: float=None, binaryFile: typing.IO=None):
        if binaryFile is not None and resolution is not None and edgeLength is not None and shore is not None and hydrology is not None:
            self._initFromBinaryFile(resolution, edgeLength, shore, hydrology, binaryFile)
        elif shore is not None and hydrology is not None and resolution is not None and edgeLength is not None:
            self._initFromModel(shore, hydrology, resolution, edgeLength)
    def _initFromModel(self, shore, hydrology, resolution, edgeLength):
        self.resolution = resolution
        self.edgeLength = edgeLength

        self.shore     = shore
        self.hydrology = hydrology
        
        points = [node.position for node in hydrology.allNodes()]

        # Add corners so that the entire area is covered
        points.append((-shore.realShape[0],-shore.realShape[1]))
        points.append((-shore.realShape[0],shore.realShape[1]))
        points.append((shore.realShape[0],shore.realShape[1]))
        points.append((shore.realShape[0],-shore.realShape[1]))
        
        self.vor = Voronoi(points,qhull_options='Qbb Qc Qz Qx')

        # Parallel to points. This indicates which voronoi region the point
        # is associated with. It is an index in vor.regions
        self.point_region = list(self.vor.point_region)

        # This is a list of lists. For each voronoi region, it is a list of
        # all the voronoi vertices that bind the region. Each number is an
        # index in vor.vertices
        self.regions = self.vor.regions

        # A list of coordinates, one for each voronoi vertex
        self.vertices = self.vor.vertices

        # This is just point_region, but for reverse lookup
        self.region_point = {self.point_region[i]: i for i in range(len(self.point_region))}

        # sort region vertices so that the polygons are convex
        print('\tOrganizing vertices into convex polygons...')
        for rid in trange(len(self.regions)):
            nodeID = self.id_vor_region(rid)
            if nodeID is None or nodeID >= len(self.hydrology):
                continue # if this region is not associated with a node, don't bother
            region = [iv for iv in self.regions[rid] if iv != -1]
            pivotPoint = self.hydrology.node(nodeID).position
            region.sort(key = lambda idx: math.atan2(
                    self.vertices[idx][1] - pivotPoint[1], # y
                    self.vertices[idx][0] - pivotPoint[0]  # x
            ))
            self.regions[rid] = region

        print('\tCreating ridge primitives...')
        self.qs = [ ]
        for iv in trange(len(self.vertices)):
            if not shore.isOnLand(self.vertices[iv]):
                self.qs.append(None)
                continue
            nearbyNodes = [ ]
            tryDistance = edgeLength * 2
            # Ensure that we at least find _some_ nodes
            while (len(nearbyNodes) < 1):
                nearbyNodes = self.hydrology.query_ball_point(self.vertices[iv], tryDistance)
                tryDistance *= 2
            borderedNodes = [ ]
            for nodeID in nearbyNodes:
                if iv in self.regions[self.vor_region_id(nodeID)]:
                    borderedNodes.append(nodeID)
            self.qs.append(Q(self.vertices[iv], borderedNodes, iv))

        print('\tClassifying ridges...')
        self.cellsRidges = { }
        self.cellsDownstreamRidges = { }
        # Classify all ridges
        for ri in trange(len(self.vor.ridge_vertices)):
            for n in self.vor.ridge_points[ri]: # Each ridge separates exactly two nodes
                if n >= len(self.hydrology):
                    continue # Apparently there are nodes that don't exist; skip these
                node = hydrology.node(n)
                otherNode = self.vor.ridge_points[ri][self.vor.ridge_points[ri] != n][0]
                # if this ridge is the outflow ridge for this node, mark it as such and move on
                if node.parent is not None and node.parent.id == otherNode:
                    # this ridge is the outflow ridge for this node
                    v1 = self.vertices[self.vor.ridge_vertices[ri][0]]
                    v2 = self.vertices[self.vor.ridge_vertices[ri][1]]
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

        print('\tFinding unaffiliated vertices...')
        # Add vertices that are not attached to ridges
        for n in trange(len(self.hydrology)):
            verts = self.regions[self.vor_region_id(n)].copy()
            if n not in self.cellsRidges:
                self.cellsRidges[n] = [ ]
            # Eliminate vertices that are attached to ridges
            for ridge in self.cellsRidges[n]:
                for q in ridge:
                    if q is not None and q.vorIndex in verts:
                        verts.remove(q.vorIndex)
            for v in verts:
                if not self.shore.isOnLand(self.vertices[v]):
                    continue
                self.cellsRidges[n].append((self.qs[v],))
    def _initFromBinaryFile(self, resolution, edgeLength, shore, hydrology, binaryFile):
        self.edgeLength = edgeLength
        self.shore = shore
        self.hydrology = hydrology

        points = [node.position for node in hydrology.allNodes()]
        points.append((0,0))
        points.append((0,shore.realShape[1]))
        points.append((shore.realShape[0],0))
        points.append((shore.realShape[0],shore.realShape[1]))
        
        self.vor = Voronoi(points,qhull_options='Qbb Qc Qz Qx')

        self.resolution = resolution
        self.point_region = [ ]
        numPoints = readValue('!Q', binaryFile)
        for i in range(numPoints):
            self.point_region.append(readValue('!Q', binaryFile))

        self.region_point = {self.point_region[i]: i for i in range(len(self.point_region))}

        self.regions = [ ]
        numRegions = readValue('!Q', binaryFile)
        for r in range(numRegions):
            regionArray = [ ]
            numVertices = readValue('!B', binaryFile)
            for v in range(numVertices):
                idx = readValue('!Q', binaryFile)
                if idx == 0xffffffffffffffff:
                    idx = -1
                regionArray.append(idx)
            self.regions.append(regionArray)
        
        self.vertices = [ ]
        numVertices = readValue('!Q', binaryFile)
        for i in range(numVertices):
            self.vertices.append((readValue('!f', binaryFile),readValue('!f', binaryFile)))

        self.qs = [ ]
        numQs = readValue('!Q', binaryFile)
        for i in range(numQs):
            if readValue('!B', binaryFile) == 0x00:
                self.qs.append(None)
                continue
            position = (readValue('!f', binaryFile),readValue('!f', binaryFile))
            nodeIdxes = [ ]
            numBorders = readValue('!B', binaryFile)
            for j in range(numBorders):
                nodeIdxes.append(readValue('!Q', binaryFile))
            voronoiIdx = readValue('!Q', binaryFile)
            q = Q(position, nodeIdxes, voronoiIdx)
            q.elevation = readValue('!f', binaryFile)
            self.qs.append(q)

        self.cellsRidges = { }
        numCells = readValue('!Q', binaryFile)
        for i in range(numCells):
            cellID = readValue('!Q', binaryFile)
            numRidges = readValue('!B', binaryFile)
            self.cellsRidges[cellID] = [ ]
            for j in range(numRidges):
                if readValue('!B', binaryFile) < 2:
                    self.cellsRidges[cellID].append( ( self.qs[readValue('!Q', binaryFile)],) )
                else:
                    self.cellsRidges[cellID].append((
                        self.qs[readValue('!Q', binaryFile)],
                        self.qs[readValue('!Q', binaryFile)]
                    ))
        
        self.cellsDownstreamRidges = { }
        numCells = readValue('!Q', binaryFile)
        for i in range(numCells):
            cellID = readValue('!Q', binaryFile)
            if readValue('!B', binaryFile) == 0xff:
                self.cellsDownstreamRidges[cellID] = None
            else:
                end0x = readValue('!f', binaryFile)
                end0y = readValue('!f', binaryFile)
                end1x = readValue('!f', binaryFile)
                end1y = readValue('!f', binaryFile)
                self.cellsDownstreamRidges[cellID] = (
                    (end0x, end0y), (end1x, end1y)
                )
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
        return self.point_region[node]
    def id_vor_region(self, regionID: int) -> int:
        """ Returns the index of the *node*

        This method is the opposite of :func:`TerrainHoneycomb.vor_region_id`

        :param regionID: The voronoi region id
        :type regionID: int

        :return: The ID of the hydrology node that this region corresponds to. If this region is not associated with an input point, None is returned
        :rtype: int
        .. note::
            This method is also for internal use
        """
        try:
            return self.region_point[regionID]
        except:
            return None # This region does not correspond to an input point
    def ridgePositions(self, node: int) -> typing.List[typing.Tuple[float,float]]:
        """Gets the position of each vertex that makes up the cell

        :param node: The ID of the node whose cell you wish to query
        :type node: int
        :return: A list of points that correspond to the vertices that make the cell
        :rtype: list[tuple[float,float]]
        """
        ridges = self.regions[self.vor_region_id(node)] # the indices of the vertex boundaries
        return [self.vertices[x] for x in ridges if x != -1] # positions of all the vertices
    def cellVertices(self, nodeID: int) -> typing.List[typing.Tuple[float,float]]:
        """Gets the coordinates of the vertices that define the shape of the node's cell

        .. todo::
            There are a number of degenerate cases where a cell may not have all of
            its vertices. This is most common with nodes adjacent to the seeee. In such
            cases, one or more vertices may not be on land, in which case those vertices
            will not be returned. This necessitates a number of workarounds in other
            methods, such as :func:`TerrainHoneycomb.cellArea`,
            :func:`TerrainHoneycomb.boundingBox`, :func:`TerrainHoneycomb.nodeID`, and
            anything that calls them. Moreover, :func:`TerrainHoneycomb.isInCell` also
            needs fixing.

        :param nodeID: The ID of the node whose shape you wish to query
        :type nodeID: int
        :return:
        :rtype: list[tuple[float,float]]
        """
        return [self.vertices[vi] for vi in self.regions[self.vor_region_id(nodeID)] if vi != -1 and self.shore.isOnLand(self.vertices[vi])]
    def cellArea(self, node: HydroPrimitive) -> float:
        """Calculates the (rough) area of the cell that a location is in

        This method derives the area based on the cell's shape. It is accurate.
        But in cases where the cell's shape is malformed, this method will not
        be accurate, and will simply return `resolution**2`, which is essentially
        the area of a single pixel.

        :param loc: Any location within the cell you wish to query
        :type loc: tuple[float,float]
        """
        try:
            return Math.convexPolygonArea(
                node.position,
                self.cellVertices(node.id)
            )
        except:
            return self.resolution**2
    def cellQs(self, node: int) -> typing.List[Q]:
        """Returns all the Qs binding the cell that corresponds to the given node

        :param node: The node you wish to query
        :type node: int
        :return: List of Q instances
        :rtype: list[Q]
        """
        return [self.qs[vorIdx] for vorIdx in self.regions[self.vor_region_id(node)] if self.qs[vorIdx] is not None]
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
        vertices = self.cellVertices(n) # vertices binding the region
        if len(vertices) < 1:
            # If this cell has a malformed shape, don't
            return (None, None, None, None)
        xllim = min([v[0] for v in vertices])
        xulim = max([v[0] for v in vertices])
        yllim = min([v[1] for v in vertices])
        yulim = max([v[1] for v in vertices])
        return (xllim, xulim, yllim, yulim)
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
        return self.shore.isOnLand(p) and Math.pointInConvexPolygon(p, self.cellVertices(n), self.hydrology.node(n).position)
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

        Like :func:`cellArea`, this is not entirely precise

        :param point: The point you wish to test
        :type point: tuple[float,float]
        :return: The ID of a node/cell (Returns None if it isn't in a valid cell)
        :rtype: int
        """
        # check hydrology nodes within a certain distance
        for id in self.hydrology.query_ball_point(point, self.edgeLength):
            # if this point is within the voronoi region of one of those nodes,
            # then that is the point's node
            vertices = self.cellVertices(id)
            if len(vertices) < 1:
                # Ignore cells with malformed shapes
                continue
            if Math.pointInConvexPolygon(point, vertices, self.hydrology.node(id).position):
                return id
        return None

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
    def __init__(self, hydrology=None, cells=None, num_points=None, binaryFile=None):
        if binaryFile is not None:
            self._initReconstitute(binaryFile)
        elif hydrology is not None and cells is not None and num_points is not None:
            self._initCreate(hydrology, cells, num_points)
        else:
            raise ValueError()
    def _initCreate(self, hydrology, cells, num_points):
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
        for n in trange(len(hydrology)):
            xllim, xulim, yllim, yulim = cells.boundingBox(n)
            if xllim is None:
                # Ignore cells that are too small
                continue
            
            # I think this applies a mask to the poisson points, and adds those points as Tees for the cell
            points_projected = [ [p[0]*(xulim-xllim)+xllim,p[1]*(yulim-yllim)+yllim] for p in points ]
            points_filtered = [ (p[0],p[1]) for p in points_projected if cells.isInCell(p,n) ]
            cellTs = [T(p,n) for p in points_filtered]
            self.cellTs[n] = cellTs
            self.tList += cellTs

        allpoints_list = [[t.position[0],t.position[1]] for t in self.allTs()]
        allpoints_nd = np.array(allpoints_list)
        self.apkd = cKDTree(allpoints_nd)
    def _initReconstitute(self, binaryFile):
        self.cellTs = { }
        self.tList = [ ]
        allpoints_list = [ ]

        numPrimitives = readValue('!Q', binaryFile)
        for i in range(numPrimitives):
            loc = (readValue('!f', binaryFile), readValue('!f', binaryFile))
            cellID = readValue('!I', binaryFile)
            elevation = readValue('!f', binaryFile)

            t = T(loc,cellID)
            t.elevation = elevation

            if cellID not in self.cellTs:
                self.cellTs[cellID] = [ ]
            self.cellTs[cellID].append(t)
            self.tList.append(t)
            allpoints_list.append(loc)

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
    def getT(self, tid: int) -> T:
        """Gets a terrain primitive identified by its index in the list of all primitives

        :param tid: The index of the primitive you wish to retrieve
        :type tid: int
        :return: The terrain primitive
        :rtype: T
        """
        return self.tList[tid]
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
    def __len__(self) -> int:
        """Returns the number of nodes in the forest

        :return: The number of primitives on the map
        :rtype: int
        """
        return len(self.tList)

def readValue(type, stream):
    buffer = stream.read(struct.calcsize(type))
    return struct.unpack(type, buffer)[0]