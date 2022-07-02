import random
import math
import struct

import Math

# imports and definitions for type hinting
import typing
import DataModel
HydroPrimitive = DataModel.HydroPrimitive
HydrologyNetwork = DataModel.HydrologyNetwork

def selectNode(candidate_nodes: typing.List[HydroPrimitive] , zeta: float) -> HydroPrimitive:
    """Given a list of candidate nodes, this function selects the next node to expand

    The selection is based on Genevaux et al §4.2.1

    :param candidate_nodes: The nodes that can still be expanded
    :type candidate_nodes: list[HydroPrimitive]
    :param zeta: Basically determines how much to prioritize elevation in the selection process (see Genevaux et al §4.2.1)
    :type zeta: int
    :return: The node to expand in the next step
    :rtype: HydroPrimitive
    """
    lowestCandidateZ = min([node.elevation for node in candidate_nodes]) # elevation of lowest candidate
    subselection = [n for n in candidate_nodes if n.elevation < lowestCandidateZ+zeta ] # 
    subselection.sort(key = lambda r : r.priority,reverse = True)
    subsubselection=[node for node in subselection if node.priority == subselection[0].priority]
    if len(subsubselection) > 1:
        subsubselection.sort(key = lambda n : n.elevation)
    
    return subsubselection[0]

class HydrologyParameters:
    """A simple struct that carries the paramaters relevant to expanding the river network

    It's pretty self-explanatory.

    .. todo::
       There must be a better way to do this. In the future it might make sense to
       implement the Factory Pattern to generate the :class:`DataModel.HydrologyNetwork`.
    """
    def __init__(
        self, shore, hydrology, Pa, Pc, maxTries, riverAngleDev, edgeLength, sigma, eta, zeta, riverSlope, slopeRate, candidates
    ):
        self.shore = shore
        self.hydrology = hydrology
        self.Pa = Pa
        self.Pc = Pc
        self.maxTries = maxTries
        self.riverAngleDev = riverAngleDev
        self.edgeLength = edgeLength
        self.sigma = sigma
        self.eta = eta
        self.zeta = zeta
        self.riverSlope = riverSlope
        self.slopeRate = slopeRate
        self.candidates = candidates
    
    def toBinary(self):
        binary =          struct.pack('!f', 0)
        binary = binary + struct.pack('!f', 0)
        binary = binary + struct.pack('!f', self.riverSlope.xSize * self.riverSlope.resolution)
        binary = binary + struct.pack('!f', self.riverSlope.ySize * self.riverSlope.resolution)
        binary = binary + struct.pack('!f', self.Pa)
        binary = binary + struct.pack('!f', self.Pc)
        binary = binary + struct.pack('!f', self.edgeLength)
        binary = binary + struct.pack('!f', self.sigma)
        binary = binary + struct.pack('!f', self.eta)
        binary = binary + struct.pack('!f', self.zeta)
        binary = binary + struct.pack('!f', self.slopeRate)
        binary = binary + struct.pack('!H', self.maxTries)
        binary = binary + struct.pack('!f', self.riverAngleDev)

        binary = binary + struct.pack('!I', self.riverSlope.xSize)
        binary = binary + struct.pack('!I', self.riverSlope.ySize)
        binary = binary + struct.pack('!f', self.shore.resolution)
        binary = binary + self.riverSlope.toBinary()

        # print('ToBinary complete!')

        binary = binary + struct.pack('!I', len(self.candidates))
        for candidate in self.candidates:
            binary = binary + struct.pack('!f', float(candidate.x()))
            binary = binary + struct.pack('!f', float(candidate.y()))
            binary = binary + struct.pack('!I', candidate.priority)
            binary = binary + struct.pack('!Q', candidate.contourIndex)
        
        binary = binary + struct.pack('!Q', len(self.shore.contour))
        for point in self.shore.contour:
            # point = point[0]
            # points in a contour array are y,x
            binary = binary + struct.pack('!Q', point[0])
            binary = binary + struct.pack('!Q', point[1])
        
        return binary

def alpha(node: HydroPrimitive, candidates: typing.List[HydroPrimitive], params: HydrologyParameters):         # alpha, as in the expansion rules in Table 1
    """Alpha node expansion rule

    :param node: The node to expand
    :type node: HydroPrimitive
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    if node.priority==1:
        ruleBase(node, candidates, params)
    else:
        Pval = random.random();
        if Pval <= params.Pa:
            rulePa(node, candidates, params)
        elif Pval <= params.Pa + params.Pc:
            rulePc(node, candidates, params)
        else:
            rulePs(node, candidates, params)

            
def ruleBase(node: HydroPrimitive, candidates: typing.List[HydroPrimitive], params: HydrologyParameters): #filling
    """Rule 1 from Table 1 in Genevaux et al

    :param node: The node to expand
    :type node: HydroPrimitive
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    #tao(priority,node)
    for i in range(random.randint(1,5)):
        beta(node, node.priority, candidates, params)
        
        
def rulePc(node: HydroPrimitive, candidates: typing.List[HydroPrimitive], params: HydrologyParameters): #rive growth
    """Rule 2.1 from Table 1 in Genevaux et al

    :param node: The node to expand
    :type node: HydroPrimitive
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    #tao(priority,node)
    beta(node, node.priority,candidates, params)
    
    
def rulePs(node: HydroPrimitive, candidates: typing.List[HydroPrimitive], params: HydrologyParameters): #symetric junction
    """Rule 2.2 from Table 1 in Genevaux et al

    :param node: The node to expand
    :type node: HydroPrimitive
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    #tao(priority,node)
    beta(node, node.priority-1, candidates, params)
    beta(node, node.priority-1, candidates, params)

    
def rulePa(node: HydroPrimitive, candidates: typing.List[HydroPrimitive], params: HydrologyParameters): # asymetric junction
    """Rule 2.3 from Table 1 in Genevaux et al

    :param node: The node to expand
    :type node: HydroPrimitive
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    #tao(priority,node)
    beta(node, node.priority, candidates, params)
    beta(random.randint(1,priority-1),node, params)
    
    
def beta(node: HydroPrimitive, priority: int, candidates: typing.List[HydroPrimitive], params: HydrologyParameters):
    """Instantiates a new node with a given priority

    :param node: The node to expand
    :type node: HydroPrimitive
    :param priority: The priority for the new node
    :type priority: int
    :param candidates: The set of candidate nodes
    :type candidates: list[HydroPrimitive]
    :param params: The parameter struct
    :type params: HydrologyParameters
    """
    point = picknewnodepos(node, params)
    if point is not None:
        slope = params.slopeRate * params.riverSlope[node.x() , node.y()]/255
        newZ = node.elevation + slope * params.edgeLength
        candidates.append(params.hydrology.addNode(point, priority=priority, elevation=newZ, parent=node))
    else:
        tao(node, candidates)

def tao(node: HydroPrimitive, candidates: typing.List[HydroPrimitive]):
    """Removes a node from the set of candidates (rule 3.2 from Table 1 of Genevaux et al)

    :param node: The node to remove
    :type node: HydroPrimitive
    :param candidates: The set (list) of candidates
    :type candidates: list[HydroPrimitive]
    """
    try:
        candidates.remove(node)
    except:
        None
    finally:
        None

def picknewnodepos(parentnode: HydroPrimitive, params: HydrologyParameters) -> typing.Optional[typing.Tuple[float,float]]:
    """Tries to pick a new position for a candidate node

    Tries to avoid the coast and crowding other nodes

    :param parentnode: The node to start from
    :type parentnode: HydroPrimitive
    :param params: The struct holding the relevant parameters
    :type params: HydrologyParameters
    :return: A new position, or None if one could not be found
    :rtype: tuple[float,float] | None
    """
    parentsparent = params.hydrology.downstream(parentnode.id) # parent node of parentnode
    
    angle = None
    if parentsparent is None: # If there is no previous node, Go in a direction perpendicular to the coast
        angle = coastNormal(parentnode, params)
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
    for i in range(params.maxTries):
        # Pick a random new point (generally in the same direction)
        newAngle = angle + random.gauss(0,params.riverAngleDev)
        newNodePos = (
            parentnode.x() + params.edgeLength*math.cos(newAngle),
            parentnode.y() + params.edgeLength*math.sin(newAngle)
        )
        if isAcceptablePosition(newNodePos, params):
            break
        else:
            newNodePos = None
    
    return newNodePos

def coastNormal(node: HydroPrimitive, params: HydrologyParameters) -> float: # Gets the angle that is approximately normal to the coast
    """Gets an angle that is (approximately) perpendicular to the coast

    This is why :class:`DataModel.HydroPrimitive` instances that represent
    a mouth node have a ``contourIndex`` attribute.

    :param node: The node to query
    :type node: HydroPrimitive:
    :param params: The struct holding the relevant parameters
    :type params: HydrologyParameters
    :return: The angle in radians
    :rtype: float
    """
    assert node.contourIndex is not None # assert that this is a mouth node
    p1 = params.shore[node.contourIndex+3]
    p2 = params.shore[node.contourIndex-3]
    theta = math.atan2(p2[1]-p1[1],p2[0]-p1[0])
    return theta - 0.5*math.pi

def isAcceptablePosition(point: typing.Tuple[float,float], params: HydrologyParameters) -> bool:
    """Determines whether or not a point is acceptable as a new node's location

    :param point: The point in question
    :type point: tuple[float,float]
    :param params: The struct holding the relevant parameters
    :type params: HydrologyParameters
    :return: True if the point meets the criteria, False if otherwise
    :rtype: bool
    """
    if point is None:
        return False
    # is the point too close to the seeeeeeeeeeeee?
    if params.shore.distanceToShore(point) < params.eta * params.edgeLength:
        #print(f'Angle too close to the seeee (distance {cv.pointPolygonTest(contour,(point[1],point[0]),True)})')
        return False
    # is the point too close to other nodes?
    for node0,node1 in params.hydrology.edgesWithinRadius(point, 2*params.edgeLength): # Go through each edge
        dist = Math.point_segment_distance( # Distance to the edge (edge is a line segment)
            point[0],point[1],
            node0.x(),node0.y(), # Line segment endpoint 1 (x,y)
            node1.x(),node1.y()  # Line segment endpoint 2 (x,y)
        )
        if dist < params.sigma * params.edgeLength:
            return False
    # otherwise return True
    return True

def classify(node: HydroPrimitive, hydrology: HydrologyNetwork, edgeLength: float):
    """Computes the Rosgen classification for a stretch of river

    :param node: The node to classify
    :type node: HydroPrimitive
    :param hydrology: The hydrological network
    :type hydrology: HydrologyNetwork
    :param edgeLength: The length of each edge between nodes
    :type edgeLength: float

    .. todo::
       This classification is just a slapdash placeholder. A real classification should
       be written in the future.
    """
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

def getLocalWatershed(node: HydroPrimitive, cells: DataModel.TerrainHoneycomb) -> float:
    """Gets the area of the watershed represented by this particular cell. This is just the cell's area.
    
    :param node: The node of the cell
    :type node: HydroPrimitive
    :param cells: The terrain honeycomb for the terrain
    :type cells: TerrainHoneycomb
    :return: The area of the local watershed in square meters
    :rtype: float
    """
    return cells.cellArea(node)

def getInheritedWatershed(node: HydroPrimitive, hydrology: HydrologyNetwork) -> float:
    """Gets the total area of the watershed.

    This is the area of all the cells that drain into this cell. This includes this cell's area.

    :param node: The node of the cell
    :type node: HydroPrimitive
    :param hydrology: The hydrology network for the terrain
    :type hydrology: HydrologyNetwork
    :return: The area of the inherited watershed in square meters
    :rtype: float
    
    .. note::
        To work properly, the nodes must all have their local watershed areas set. Use :func:`HydrologyFunctions.getLocalWatershed` for this.
    """
    return node.localWatershed + sum([n.inheritedWatershed for n in hydrology.upstream(node.id)])

def getFlow(inheritedWatershed: float) -> float:
    """This estimates the flow of water through a cell.

    This is based on a simple formula specified in §5.1 of Geneveaux et al

    :param inheritedWatershed: The total area of the watershed. Use :func:`HydrologyFunctions.getInheritedWatershed`
    :type inheritedWatershed: float
    :return: The approximate flow in cubic meters per second
    :rtype: float
    """
    return 0.42 * inheritedWatershed**0.69