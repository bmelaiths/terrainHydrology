import numpy as np
from scipy.spatial import Voronoi
import math

from tqdm import trange

from DataModel import ShoreModel, TerrainHoneycomb, Q, HydrologyNetwork, RasterData

def getRidgeElevation(q: Q, hydrology: HydrologyNetwork, terrainSlope: RasterData, terrainSlopeRate: float) -> float:
    nodes = [hydrology.node(n) for n in q.nodes]
    maxElevation = max([node.elevation for node in nodes])
    d = np.linalg.norm(q.position - nodes[0].position)
    slope = terrainSlopeRate * terrainSlope[q.position[0],q.position[1]] / 255
    return maxElevation + d * slope

def id_vor_region(regionID: int, region_point) -> int:
    """ Returns the index of the *node*

    This function is the opposite of :func:`TerrainHoneycomb.vor_region_id`

    :param regionID: The voronoi region id
    :type regionID: int

    :return: The ID of the hydrology node that this region corresponds to. If this region is not associated with an input point, None is returned
    :rtype: int
    .. note::
        This method is also for internal use
    """
    try:
        return region_point[regionID]
    except:
        return None # This region does not correspond to an input point

def vor_region_id(node: int, point_region) -> int:
    """Returns the index of the *voronoi region*

    This function is useful because the index of the voronoi region is not the same
    as the ID of the cell it is based on.

    :param node: The ID of the node/cell
    :type node: int

    :return: The ID of the *voronoi region* (or, cell) associated with that node
    :rtype: int

    .. note::
        This method is for internal use. If you are using it from outside this
        class, your approach is definitely breaking a key design principle.
    """
    return point_region[node]

def initializeTerrainHoneycomb(shore: ShoreModel, hydrology: HydrologyNetwork, resolution: float, edgeLength: float) -> TerrainHoneycomb:
    cells = TerrainHoneycomb()
    
    points = [node.position for node in hydrology.allNodes()]

    # Add corners so that the entire area is covered
    points.append((-shore.realShape[0],-shore.realShape[1]))
    points.append((-shore.realShape[0],shore.realShape[1]))
    points.append((shore.realShape[0],shore.realShape[1]))
    points.append((shore.realShape[0],-shore.realShape[1]))
    
    vor = Voronoi(points,qhull_options='Qbb Qc Qz Qx')

    # Parallel to points. This indicates which voronoi region the point
    # is associated with. It is an index in vor.regions
    point_region = list(vor.point_region)

    # This is a list of lists. For each voronoi region, it is a list of
    # all the voronoi vertices that bind the region. Each number is an
    # index in vor.vertices
    regions = vor.regions

    # A list of coordinates, one for each voronoi vertex
    vertices = vor.vertices

    # This is just point_region, but for reverse lookup
    region_point = {point_region[i]: i for i in range(len(point_region))}

    num_regions = len(regions)
    num_vertices = len(vertices)
    num_ridges = len(vor.ridge_vertices)

    # sort region vertices so that the polygons are convex
    print('\tOrganizing vertices into convex polygons...')
    for rid in trange(num_regions):
        nodeID = id_vor_region(rid, region_point)
        if nodeID is None or nodeID >= len(hydrology):
            continue # if this region is not associated with a node, don't bother
        region = [iv for iv in regions[rid] if iv != -1]
        pivotPoint = hydrology.node(nodeID).position
        region.sort(key = lambda idx: math.atan2(
                vertices[idx][1] - pivotPoint[1], # y
                vertices[idx][0] - pivotPoint[0]  # x
        ))
        regions[rid] = region

    print('\tCreating ridge primitives...')
    qs = [ ]
    for iv in trange(num_vertices):
        if not shore.isOnLand(vertices[iv]):
            qs.append(None)
            continue
        nearbyNodes = [ ]
        tryDistance = edgeLength * 2
        # Ensure that we at least find _some_ nodes
        while (len(nearbyNodes) < 1):
            nearbyNodes = hydrology.query_ball_point(vertices[iv], tryDistance)
            tryDistance *= 2
        borderedNodes = [ ]
        for nodeID in nearbyNodes:
            if iv in regions[vor_region_id(nodeID, point_region)]:
                borderedNodes.append(nodeID)
        qs.append(Q(vertices[iv], borderedNodes, iv))

    cellsRidges = { }
    cellsDownstreamRidges = { }
    
    # Classify all ridges
    for ri in trange(num_ridges): # Remember that this is looping through the _indexes_ of ridge_vertices. This loop goes through the ridges
        for n in vor.ridge_points[ri]: # Each ridge separates exactly two nodes. This loop goes through them
            if n >= len(hydrology):
                continue # Apparently there are nodes that don't exist; skip these
            node = hydrology.node(n)
            otherNode = vor.ridge_points[ri][vor.ridge_points[ri] != n][0]
            
            # if this ridge is the outflow ridge for this node, mark it as such and move on
            if node.parent is not None and node.parent.id == otherNode:
                # this ridge is the outflow ridge for this node
                v1 = vertices[vor.ridge_vertices[ri][0]]
                v2 = vertices[vor.ridge_vertices[ri][1]]
                if not shore.isOnLand(v1) or not shore.isOnLand(v2):
                    # If one or both vertices is not on land, then don't bother
                    # trying to make the river flow through the ridge neatly
                    cellsDownstreamRidges[n] = None
                else:
                    cellsDownstreamRidges[n] = (v1, v2)
                break # the outflow ridge of this node is an inflow ridge for the other one
            if otherNode in [nd.id for nd in hydrology.upstream(n)]:
                continue # this is an inflow ridge, so it need not be considered further

            # the ridge at ri is not transected by a river
            if n not in cellsRidges:
                cellsRidges[n] = [ ]
            vertex0 = qs[vor.ridge_vertices[ri][0]]
            vertex1 = qs[vor.ridge_vertices[ri][1]]
            if vertex0 is None or vertex1 is None:
                continue
            cellsRidges[n].append((vertex0,vertex1))

    print('\tFinding unaffiliated vertices...')
    # Add vertices that are not attached to ridges
    for n in trange(len(hydrology)):
        verts = regions[vor_region_id(n, point_region)].copy()
        if n not in cellsRidges:
            cellsRidges[n] = [ ]
        # Eliminate vertices that are attached to ridges
        for ridge in cellsRidges[n]:
            for q in ridge:
                if q is not None and q.vorIndex in verts:
                    verts.remove(q.vorIndex)
        for v in verts:
            if not shore.isOnLand(vertices[v]):
                continue
            cellsRidges[n].append((qs[v],))

    cells.shore = shore
    cells.hydrology = hydrology

    cells.vertices = vertices
    cells.regions = regions
    cells.point_region = point_region

    cells.region_point = region_point

    cells.qs = qs

    cells.cellsRidges = cellsRidges
    cells.cellsDownstreamRidges = cellsDownstreamRidges

    return cells