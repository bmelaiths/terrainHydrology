import shapely.geometry as geom
import math
from scipy.spatial import cKDTree
import numpy as np
from tqdm import trange

from poisson import PoissonGenerator

from DataModel import T, ShoreModel, HydrologyNetwork, TerrainHoneycomb, Terrain
import Math

def initializeTerrain(hydrology: HydrologyNetwork, cells: TerrainHoneycomb, num_points: int) -> Terrain:
    terrain = Terrain()

    terrain.cellTs = { }
    terrain.tList = [ ]
    
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
        terrain.cellTs[n] = cellTs
        terrain.tList += cellTs

    allpoints_list = [[t.position[0],t.position[1]] for t in terrain.allTs()]
    allpoints_nd = np.array(allpoints_list)
    terrain.apkd = cKDTree(allpoints_nd)

    return terrain

def computePrimitiveElevation(t: T, shore: ShoreModel, hydrology: HydrologyNetwork, cells: TerrainHoneycomb) -> float:
    ridges = cells.cellRidges(t.cell)

    # find distance to closest sgment, and elevation at that point
    closestRdist = None
    ridgeElevation = None
    for ridge in ridges:
        if len(ridge) < 2:
            q0 = ridge[0]
            dist = Math.distance(q0.position,t.position)
            if closestRdist is None or dist < closestRdist:
                closestRdist = dist
                ridgeElevation = q0.elevation
            continue
        
        q0 = ridge[0]
        q1 = ridge[1]
        dist, isToEndpoint = Math.point_segment_distance_is_endpoint(
            t.position[0],t.position[1],
            q0.position[0],q0.position[1],
            q1.position[0],q1.position[1]
        )
        if closestRdist is not None and dist > closestRdist:
            continue
        if isToEndpoint:
            if Math.distance(q0.position,t.position) < Math.distance(q1.position,t.position):
                closestRdist = Math.distance(q0.position,t.position)
                ridgeElevation = q0.elevation
            else:
                closestRdist = Math.distance(q1.position,t.position)
                ridgeElevation = q1.elevation
        else:
            closestRdist = dist
            try:
                ridgeElevation = q0.elevation + (math.sqrt(Math.distance(q0.position,t.position)**2 - dist**2) / Math.distance(q0.position,q1.position)) * (q1.elevation - q0.elevation)
            except:
                print(f'That math domain error has occured')
                print(f'q0.elevation: {q0.elevation}, q0.position: {q0.position}, t.positon: {t.position}, dist: {dist}, q1.position: {q1.position}, q1.elevation: {q1.elevation}')
                exit()
    
    # breakpoint()
    
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
        # gets the river that is closest to the terrain primitive
        rividx = [point.distance(x) for x in local_rivers].index(min([point.distance(x) for x in local_rivers]))
        # gets the point along the river that is the distance along the river to the point nearest to the Tee
        projected = local_rivers[rividx].interpolate(local_rivers[rividx].project(point))
        distancefromN = point.distance(local_rivers[rividx]) # distance to that point
    else: # handle cases of stub rivers
        projected = geom.Point(node.x(),node.y(),node.elevation)
        distancefromN = point.distance(projected)
    
    if distancefromN==0 and closestRdist==0:
        distancefromN=1
    
    # this is the weighted average of the 2 elevations
    lerpedelevation = projected.z*(closestRdist/(closestRdist+distancefromN))+ridgeElevation*(distancefromN/(closestRdist+distancefromN))

    return lerpedelevation
