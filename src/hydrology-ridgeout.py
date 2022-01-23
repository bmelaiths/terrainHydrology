#! /usr/bin/env python

import argparse
import shapefile
from tqdm.std import trange

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
    w.field('id', 'L')

    for ridge in cells.cellsRidges.values():
        coords = [ ]

        if type(ridge) != int:
            coords.append(cells.qs[ridge[0]].position)
            if len(ridge) > 1:
                coords.append(cells.qs[ridge[1]].position)
        else:
            continue
        
        coords = [((p[0])-(realShape[1]*0.5),(realShape[0]-p[1])-(realShape[0]*0.5)) for p in coords]

        w.record(True)

        w.line([list(coords)])

    w.close()
