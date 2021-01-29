import random
import math

import Math

def selectNode(candidate_nodes,zeta):

    lowestCandidateZ = min([node.elevation for node in candidate_nodes]) # elevation of lowest candidate
    subselection = [n for n in candidate_nodes if n.elevation < lowestCandidateZ+zeta ] # 
    subselection.sort(key = lambda r : r.priority,reverse = True)
    subsubselection=[node for node in subselection if node.priority == subselection[0].priority]
    
    return subsubselection[0]

class HydrologyParameters:
    def __init__(
        self, shore, hydrology, Pa, Pc, maxTries, riverAngleDev, edgeLength, sigma, eta, riverSlope, slopeRate
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
        self.riverSlope = riverSlope
        self.slopeRate = slopeRate

def alpha(node, candidates, params):         # alpha, as in the expansion rules in Table 1
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

            
def ruleBase(node, candidates, params): #filling
    #tao(priority,node)
    for i in range(random.randint(1,5)):
        beta(node, node.priority, candidates, params)
        
        
def rulePc(node, candidates, params): #rive growth
    #tao(priority,node)
    beta(node, node.priority,candidates, params)
    
    
def rulePs(node, candidates, params): #symetric junction
    #tao(priority,node)
    beta(node, node.priority-1, candidates, params)
    beta(node, node.priority-1, candidates, params)

    
def rulePa(priority,node, params): # asymetric junction
    #tao(priority,node)
    beta(node, node.priority, candidates, params)
    beta(random.randint(1,priority-1),node, params)
    
    
def beta(node, priority, candidates, params):
    point = picknewnodepos(node, params)
    if point is not None:
        slope = 2.0 * params.riverSlope[ int(node.x()) , int(node.y())] / 255
        newZ = node.elevation + random.random() * params.slopeRate * slope
        candidates.append(params.hydrology.addNode(point, priority=priority, elevation=newZ, parent=node))
    else:
        tao(node, candidates)

def tao(node, candidates):
    try:
        candidates.remove(node)
    except:
        None
    finally:
        None

def picknewnodepos(parentnode, params):
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

def coastNormal(node, params): # Gets the angle that is approximately normal to the coast
    assert node.contourIndex is not None # assert that this is a mouth node
    p1 = params.shore[node.contourIndex+3]
    p2 = params.shore[node.contourIndex-3]
    theta = math.atan2(p2[1]-p1[1],p2[0]-p1[0])
    return theta + 0.5*math.pi

def isAcceptablePosition(point, params):
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

def calculateHorton_Strahler(selectedCandidate, hydrology):
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

def classify(node, hydrology, edgeLength):
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