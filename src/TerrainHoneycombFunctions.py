import numpy as np
from scipy.spatial import Voronoi
import math

from tqdm import trange

from DataModel import ShoreModel, TerrainHoneycomb, Q, HydrologyNetwork, RasterData

def getRidgeElevation(q: Q, hydrology: HydrologyNetwork, terrainSlope: RasterData, terrainSlopeRate: float) -> float:
    """Computes the elevation of a ridge/crest

    :param q: The ridge primitive to compute an elevation for
    :type q: Q
    :param hydrology: The hydrology network for the terrain
    :type hydrology: HydrologyNetwork
    :param terrainSlope: A raster that indicates how steep the terrain should climb in different areas
    :type terrainSlope: RasterData
    :param terrainSlopeRate: This is documented in hydrology.py under "Terrain Parameters"
    :type terrainSlopeRate: float
    """
    nodes = [hydrology.node(n) for n in q.nodes]
    maxElevation = max([node.elevation for node in nodes])
    d = np.linalg.norm(q.position - nodes[0].position)
    slope = terrainSlopeRate * terrainSlope[q.position[0],q.position[1]] / 255
    return maxElevation + d * slope

def initializeTerrainHoneycomb(shore: ShoreModel, hydrology: HydrologyNetwork, resolution: float, edgeLength: float) -> TerrainHoneycomb:
    """Creates the TerrainHoneycomb for the terrain

    This function does all the work of creating the Voronoi partition and classifying all cell edges.

    :param shore: The shore for the terrain
    :type shore: ShoreModel
    :param hydrology: The hydrology network for the terrain
    :type hydrology: HydrologyNetwork
    :param resolution: The resolution of the raster in meters per pixel
    :type resolution: float
    :param edgeLength:
    :type edgeLength: float
    """
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

    # Set these attributes earlier than the others so that we can use id_vor_region() and vor_region_id()
    cells.vertices = vertices
    cells.regions = regions
    cells.point_region = point_region

    # sort region vertices so that the polygons are convex
    print('\tOrganizing vertices into convex polygons...')
    for rid in trange(num_regions):
        nodeID = cells.id_vor_region(rid)
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
            if iv in regions[cells.vor_region_id(nodeID)]:
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
        verts = regions[cells.vor_region_id(n)].copy()
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

    cells.region_point = region_point

    cells.qs = qs

    cells.cellsRidges = cellsRidges
    cells.cellsDownstreamRidges = cellsDownstreamRidges

    return cells