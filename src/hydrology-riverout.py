#! /usr/bin/env python

import argparse
import shapefile

import SaveFile

parser = argparse.ArgumentParser(
    description='Implementation of Genevaux et al., "Terrain Generation Using Procedural Models Based on Hydrology", ACM Transactions on Graphics, 2013'
)
parser.add_argument(
    '-i',
    '--input',
    help='The file that contains the data model you wish to render',
    dest='inputFile',
    metavar='output/data',
    required=True
)
parser.add_argument(
    '--lat',
    help='Center latitude for the output GeoTIFF',
    dest='latitude',
    metavar='-43.2',
    required=True
)
parser.add_argument(
    '--lon',
    help='Center longitude for the output GeoTIFF',
    dest='longitude',
    metavar='-103.8',
    required=True
)
parser.add_argument(
    '-o',
    '--output',
    help='Name of the file(s) to write to',
    dest='outputFile',
    metavar='output',
    required=True
)
args = parser.parse_args()

inputFile = args.inputFile
outputFile = args.outputFile

## Create the .prj file to be read by GIS software
with open(f'{outputFile}.prj', 'w') as prj:
    prjstr = f'PROJCS["unknown",GEOGCS["GCS_unknown",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Orthographic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Longitude_Of_Center",{args.longitude}],PARAMETER["Latitude_Of_Center",{args.latitude}],UNIT["Meter",1.0]]'
    prj.write(prjstr)
    prj.close()

## Read the data
resolution, edgeLength, shore, hydrology, cells, Ts = SaveFile.readDataModel(
    inputFile
)

realShape = shore.realShape

with shapefile.Writer(outputFile, shapeType=3) as w:
    # The only relevant field for rivers
    w.field('flow', 'F')

    # This loop adds rivers in the same way that they were created
    for node in hydrology.allMouthNodes():
        leaves = hydrology.allLeaves(node.id)
        for leafNode in leaves:
            # Path from the leaf to the sea
            path = hydrology.pathToNode(node.id, leafNode.id)
            path.reverse()

            # The flow of the river is the flow at the base of
            # its path to the sea, unless this stretch of the
            # river doesn't flow all the way to the sea
            riverFlow = path[len(path)-1].flow
            for ni in range(1,len(path)):
                upstreamFlow = max([n.flow for n in hydrology.upstream(path[ni].id)])
                # If this river is merely a tributary to a larger
                # stream, terminate the course and record the
                # flow at that point
                if upstreamFlow > path[ni-1].flow:
                    riverFlow = path[ni].flow
                    break
            w.record(riverFlow)

            coords = list(leafNode.rivers[0].coords)
            # Transform the coordinates
            coords = [((p[0])-(realShape[0]*0.5),(realShape[1]-p[1])-(realShape[1]*0.5)) for p in coords]
            w.line([list(coords)])
    w.close()
