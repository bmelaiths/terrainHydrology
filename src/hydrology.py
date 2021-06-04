#!/usr/bin/env python

import datetime
import argparse
from os import read
import random
import matplotlib.pyplot as plt
import cv2 as cv
import numpy as np
from scipy.spatial import voronoi_plot_2d
import networkx as nx
from scipy import interpolate
import shapely.geometry as geom
import numpy as np
from shapely.geometry import asLineString
from multiprocessing import Process, Pipe, Queue
from tqdm import trange
import math
import rasterio
from rasterio.transform import Affine
import subprocess
import os.path
import struct

import DataModel
import HydrologyFunctions
import Math
import SaveFile

buildRiversExe = 'src/buildRivers'
computePrimitivesExe = 'src/terrainPrimitives'

# Get inputs

parser = argparse.ArgumentParser(
    description='Implementation of Genevaux et al., "Terrain Generation Using Procedural Models Based on Hydrology", ACM Transactions on Graphics, 2013'
)
parser.add_argument(
    '-g',
    '--gamma',
    help='An outline of the shore. Should be a grayscale image (but that doesn\'t have to be the actual color model)',
    dest='inputDomain',
    metavar='gamma.png',
    required=True
)
parser.add_argument(
    '-s',
    '--river-slope',
    help='Slope of the rivers (not terrain slope). Values in grayscale.',
    dest='inputRiverSlope',
    metavar='rivers.png',
    required=True
)
parser.add_argument(
    '-t',
    '--terrain-slope',
    help='Slope of the terrain (not river slope). Values in grayscale.',
    dest='inputTerrain',
    metavar='terrain.png',
    required=True
)
parser.add_argument(
    '-ri',
    '--input-resolution',
    help='The spatial resolution of the input images in meters per pixel',
    dest='resolution',
    metavar='87.5',
    required=True
)
parser.add_argument(
    '-p',
    '--num-points',
    help='The (rough) number of terrain primitives for each cell',
    dest='num_points',
    metavar='50',
    required=True
)
parser.add_argument(
    '--num-rivers',
    help='The number of drainages along the coast',
    dest='numRivers',
    metavar='10',
    default=10,
    required=False
)
parser.add_argument(
    '--accelerate',
    help='Accelerate Your Life™ using parallel processing with a natively-compiled module',
    action='store_true',
    dest='accelerate',
    required=False
)
parser.add_argument(
    '--num-procs',
    help='The number of processes/threads to use for calculating terrain primitives. This should be the number of cores you have on your system.',
    dest='num_procs',
    metavar='4',
    default=4,
    required=False
)
parser.add_argument(
    '-o',
    '--output',
    help='File that will contain the data model',
    dest='outputFile',
    metavar='outputFile',
    required=True
)
args = parser.parse_args()


## Global Variables

# Inputs
inputDomain = args.inputDomain
inputTerrain = args.inputTerrain
inputRiverSlope = args.inputRiverSlope
resolution = float(args.resolution) # meters per pixel length

# Random Number Generation
globalseed=4314
random.seed(globalseed)

## Hydrology Parameters

# Number of river mouths
N_majorRivers=int(args.numRivers)

# Branching Parameters
Ps = 0.3 #0.05 ## probability of symetric branch
Pa = 0 #0.3 ## probability of asymetric branch
Pc = 1-(Ps+Pa) ## probability of continium (growth)
zeta = 100 # elevation range to include in candidate node selection
riverAngleDev = 1.7 # Used in picknewnodepos(). Standard Deviation of angle for new node. Increase for less straight rivers
maxTries = 15

# Hydrological slope parameters
slopeRate = 0.1 # Maximum rate at which rivers climb in vertical meters per horizontal meter

# Hydrological resolution parameters
edgeLength = 2320.5 #4000 # in meters
eta = .75   #   eta * edgeLength is the minimum distance from a node to the coast
sigma = .75 # sigma * edgeLength is the minimum distance between two nodes

## Terrain Parameters
terrainSlopeRate = 1.0 # Maximum rate at which ridges climb in vertical meters per horizontal meter
num_points = int(args.num_points) # The (rough) number of terrain primitives for each cell
numProcs = int(args.num_procs) # The number of processes to use in calculating terrain primitives

## Output File
outputFile = args.outputFile


# Load input images

shore = DataModel.ShoreModel(resolution, gammaFileName=inputDomain)

terrainSlope = DataModel.RasterData(inputTerrain, resolution)
riverSlope = DataModel.RasterData(inputRiverSlope, resolution)


## Generate river mouths

hydrology = DataModel.HydrologyNetwork()

# generate first node
firstIdx = random.randint(0,len(shore)-1)
point = shore[firstIdx]
hydrology.addNode(point, 0, random.randint(1,N_majorRivers), contourIndex=firstIdx)

dist = len(shore)/N_majorRivers
for i in range(1,N_majorRivers):
    idx = int((firstIdx+i*dist+random.gauss(0, dist/6))%len(shore))
    point = shore[idx]
    hydrology.addNode(point, 0, 1, contourIndex=idx)


## Generate river nodes

print('Generating rivers...')

candidates = hydrology.allMouthNodes() # All mouth nodes are candidates
params = HydrologyFunctions.HydrologyParameters(
    # These parameters will be needed to generate the hydrology network
    shore, hydrology, Pa, Pc, maxTries, riverAngleDev, edgeLength,
    sigma, eta, zeta, riverSlope, slopeRate, candidates
)

start, end = None, None

if not args.accelerate: # Generate the hydrology in Python
    cyclesRun = 0
    start = datetime.datetime.now()
    while len(candidates)!=0:
        selectedCandidate = HydrologyFunctions.selectNode(candidates,zeta)
        HydrologyFunctions.alpha(selectedCandidate, candidates, params)
        print(f'\tCycles: {len(hydrology)}\t{cyclesRun/(datetime.datetime.now()-start).total_seconds()} cycles/sec\r', end='')
        cyclesRun = cyclesRun + 1
    end = datetime.datetime.now()
    print()
else: # Generate the hydrology using the native module
    if not os.path.exists(buildRiversExe):
        print('The executable does not exist. Run "make buildRivers" in the src/ directory to build it.')
        exit()
    proc = subprocess.Popen( # start the native module
        ['./' + buildRiversExe],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE
    )
    proc.stdin.write(params.toBinary()) # send the parameters to the native module
    print('\tData sent to native module...')

    # Display updates as native module builds the network
    cyclesRun = 0
    start = datetime.datetime.now()
    readByte = proc.stdout.read(1)
    cyclesRun = cyclesRun + 1
    while struct.unpack('B', readByte)[0] == 0x2e:
        print(f'\tCycles: {cyclesRun}\t{cyclesRun/(datetime.datetime.now()-start).total_seconds()} cycles/sec\r', end='')
        readByte = proc.stdout.read(1)
        cyclesRun = cyclesRun + 1
    end = datetime.datetime.now()
    print()

    # Recreate hydrology with data from the native module
    print('\tReading data...')
    hydrology = DataModel.HydrologyNetwork(stream=proc.stdout)

print(f'\tGenerated {len(hydrology)} nodes in {(end-start).total_seconds()} seconds')
print(f'\tRate: {len(hydrology)/(end-start).total_seconds()} node/sec')


## Create terrain partition (voronoi cells)
print('Generating terrain ridges...')
cells = DataModel.TerrainHoneycomb(shore, hydrology, resolution, edgeLength)


## Calculate watershed areas
print('Calculating watershed areas...')

# local watershed areas
for n in range(len(hydrology)):
    node = hydrology.node(n)
    node.localWatershed = cells.cellArea(node)
    node.inheritedWatershed = 0

# calculate inherited watershed areas and flow
for node in hydrology.dfsPostorderNodes():  # search nodes in a depth-first post-ordering manner
    watershed = node.localWatershed + sum([n.inheritedWatershed for n in hydrology.upstream(node.id)])
    node.inheritedWatershed=watershed  # calculate total watershed area
    node.flow = 0.42 * watershed**0.69 # calculate river flow


## Classify river nodes
print('Classifying river nodes...')
for n in range(len(hydrology)):
    HydrologyFunctions.classify(hydrology.node(n), hydrology, edgeLength)


## Calculate ridge elevations
print('Calculating ridge elevations...')

for q in cells.allQs():
    if q is None:
        continue
    nodes = [hydrology.node(n) for n in q.nodes]
    maxElevation = max([node.elevation for node in nodes])
    d = np.linalg.norm(q.position - nodes[0].position)
    slope = terrainSlopeRate * terrainSlope[q.position[0],q.position[1]] / 255
    q.elevation = maxElevation + d * slope


## Terrain pattern
print('Generating terrain primitives...')
Ts = DataModel.Terrain(hydrology, cells, num_points)


## Generate river paths
print('Interpolating river paths...')
for node in hydrology.allMouthNodes():
    # remember that out_edges gets upstream nodes
    leaves = hydrology.allLeaves(node.id)
    for leafNode in leaves: # essentially, this loops through all the highest nodes of a particular mouth
        # path to the leaf (there's only one, so it's the shortest)
        path = hydrology.pathToNode(node.id,leafNode.id)
        path.reverse()

        # terminates the path if there is another path with greater flow, and adds the downflow ridge as a point
        for ni in range(1,len(path)):
            #print(f'path: {[n.id for n in path]}')
            #print(f'ni: {ni}')
            #print(f'Upstream of node {path[ni].id} ({path[ni].position}): {[n.id for n in hydrology.upstream(path[ni].id)]}')
            upstreamFlow = max([n.flow for n in hydrology.upstream(path[ni].id)])
            if upstreamFlow > path[ni-1].flow:
                path = path[0:ni+1]
                break
        
        x = [ ]
        y = [ ]
        z = [ ]
        for pi in range(len(path)):
            p = path[pi]
            x.append(p.x())
            y.append(p.y())
            z.append(p.elevation)
            # makes the river flow through the cell's outflow ridge (so it doesn't transect a mountain)
            if p.parent is not None and pi < len(path)-1 and cells.cellOutflowRidge(p.id) is not None:
                ridge0, ridge1 = cells.cellOutflowRidge(p.id)
                x.append((ridge0[0] + ridge1[0])/2)
                y.append((ridge0[1] + ridge1[1])/2)
                z.append((p.elevation + p.parent.elevation)/2)

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


## Calculate elevations of terrain primitives
print('Calculating terrain primitive elevations...')
def subroutine(conn: Pipe, q: Queue):
    threadID = conn.recv()
    for ti in range(threadID, len(Ts), numProcs):
        t = Ts.getT(ti)
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
        
        lerpedelevation = projected.z*(closestRdist/(closestRdist+distancefromN))+ridgeElevation*(distancefromN/(closestRdist+distancefromN))
        t.elevation = lerpedelevation

        q.put(t)

# The terrain primitives will be calculated in parallel
if not args.accelerate: # Calculate the elevations in Python
    dataQueue = Queue()
    pipes = []
    processes = []
    for p in range(numProcs):
        pipes.append(Pipe())
        processes.append(Process(target=subroutine, args=(pipes[p][1],dataQueue)))
        processes[p].start()
        pipes[p][0].send(p)
    for ti in trange(len(Ts)):
        Ts.tList[ti] = dataQueue.get()
    for p in range(numProcs):
        processes[p].join()
        pipes[p][0].close()
else:
    # Write the binary data to a file
    with open('src/binaryFile', 'w+b') as file:
        SaveFile.writeToTerrainModule(file, shore, edgeLength, hydrology, cells, Ts)
        file.close()

    # Run the native module
    if not os.path.exists(computePrimitivesExe):
        print('The executable does not exist. Run "make buildRivers" in the src/ directory to build it.')
        exit()
    primitivesProc = subprocess.Popen( # start the native module
        ['./' + computePrimitivesExe],
        stdout=subprocess.PIPE
    )

    # Display updates as native module calculates the elevations
    for tid in trange(len(Ts)):
        readByte = primitivesProc.stdout.read(1)
    readByte = primitivesProc.stdout.read(1)
    assert struct.unpack('B',readByte)[0] == 0x21

    for t in Ts.allTs():
        readByte = primitivesProc.stdout.read(struct.calcsize('!f'))
        t.elevation = struct.unpack('!f', readByte)[0]

    # clean up
    os.remove('src/binaryFile')

## Save the data
print('Writing data model...')
SaveFile.writeDataModel(outputFile, edgeLength, shore, hydrology, cells, Ts)
print('Complete')