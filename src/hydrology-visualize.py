#! /usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import numpy as np
from tqdm import trange
import networkx as nx
import cv2 as cv
from scipy.spatial import voronoi_plot_2d

import DataModel
import SaveFile
import Math

parser = argparse.ArgumentParser(
    description='Visualizes an already-generated data model. Useful for debugging or tinkering.'
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
    '-g',
    '--gamma',
    help='An outline of the shore. Should be a grayscale image (but that doesn\'t have to be the actual color model)',
    dest='inputDomain',
    metavar='gamma.png',
    required=True
)
parser.add_argument(
    '-xl', '--lower-x',
    help='The lower x bound of the area to display',
    dest='lowerX',
    metavar='123.5',
    required=False
)
parser.add_argument(
    '-yl', '--lower-y',
    help='The lower y bound of the area to display',
    dest='lowerY',
    metavar='123.5',
    required=False
)
parser.add_argument(
    '-xu', '--upper-x',
    help='The upper x bound of the area to display',
    dest='upperX',
    metavar='123.5',
    required=False
)
parser.add_argument(
    '-yu', '--upper-y',
    help='The upper y bound of the area to display',
    dest='upperY',
    metavar='123.5',
    required=False
)
parser.add_argument(
    '--hydrology-network',
    help='Show the hydrology graph',
    dest='showHydrology',
    action='store_true'
)
parser.add_argument(
    '--hydrology-network-flow',
    help='Show the hydrology graph with arrows adjusted for river flow',
    dest='showRiverFlow',
    action='store_true'
)
parser.add_argument(
    '--terrain-primitives',
    help='Show the terrain primitives in each cell',
    dest='showTerrainPrimitives',
    action='store_true'
)
parser.add_argument(
    '--river-paths',
    help='Show the interpolated river paths',
    dest='showRivers',
    action='store_true'
)
parser.add_argument(
    '--voronoi-cells',
    help='The background will be colored according to the voronoi ID of the cell',
    dest='showVoronoiCells',
    action='store_true'
)
parser.add_argument(
    '--river-heights',
    help='The background will be colored according to the height of the local river node',
    dest='showRiverHeights',
    action='store_true'
)
parser.add_argument(
    '-o',
    '--output-image',
    help='The image to write',
    dest='outputFile',
    metavar='output.png',
    required=True
)
args = parser.parse_args()

resolution, edgeLength, shore, hydrology, cells, Ts = SaveFile.readDataModel(
    args.inputFile # read the data model
)

imgOutline = cv.imread(args.inputDomain)

# boundaries of the debug image
lowerX = float(args.lowerX) if args.lowerX else 0
lowerY = float(args.lowerY) if args.lowerY else 0
upperX = float(args.upperX) if args.upperX else shore.realShape[0]
upperY = float(args.upperY) if args.upperY else shore.realShape[1]

plt.figure(figsize=(20,20)) # set the resolution of the debug image

# Generate a background image for context
if not args.showVoronoiCells and not args.showRiverHeights:
    # just show the shoreline
    imStretch = (0,int(shore.realShape[1]),int(shore.realShape[0]),0)
    plt.imshow(imgOutline, extent=imStretch)
elif args.showVoronoiCells:
    rastershape = (
        int(500 * (upperX-lowerX)/(upperY-lowerY) if upperY-lowerY > upperX-lowerX else 500),
        int(500 * (upperY-lowerY)/(upperX-lowerX) if upperX-lowerX > upperY-lowerY else 500)
    )
    imgvoronoi = np.zeros((rastershape[1], rastershape[0]),dtype=np.uint32)
    print('Rendering voronoi cell raster...')
    for y in trange(rastershape[0]):
        for x in range(rastershape[1]):
            nodeID = cells.nodeID((
                lowerX + (y/rastershape[0]) * (upperX-lowerX),
                lowerY + (x/rastershape[1]) * (upperY-lowerY)
            ))
            if nodeID is not None:
                imgvoronoi[x][y] = cells.vor_region_id(nodeID)
            else:
                imgvoronoi[x][y] = 0
    plt.imshow(imgvoronoi, extent=(lowerX, upperX, upperY, lowerY))
elif args.showRiverHeights:
    rastershape = (
        int(500 * (upperX-lowerX)/(upperY-lowerY) if upperY-lowerY > upperX-lowerX else 500),
        int(500 * (upperY-lowerY)/(upperX-lowerX) if upperX-lowerX > upperY-lowerY else 500)
    )
    imgRiverHeights = np.zeros((rastershape[1], rastershape[0]),dtype=np.float)
    print('Rendering river reight raster...')
    for y in trange(rastershape[0]):
        for x in range(rastershape[1]):
            point = (
                lowerX + (y/rastershape[0]) * (upperX-lowerX),
                lowerY + (x/rastershape[1]) * (upperY-lowerY)
            )
            if shore.isOnLand(point):
                nodeID = cells.nodeID(point)
                if nodeID is not None:
                    imgRiverHeights[x][y] = hydrology.node(nodeID).elevation
                else:
                    imgRiverHeights[x][y] = 0
            else:
                imgRiverHeights[x][y] = 0
    plt.imshow(imgRiverHeights, extent=(lowerX, upperX, upperY, lowerY))

# Show the terrain primitives, if applicable
if args.showTerrainPrimitives:
    print('Plotting the terrain primitives...')
    plt.gca().scatter(
        *zip(*[t.position for t in Ts.allTs()]),
        c=[t.elevation for t in Ts.allTs()],
        cmap=plt.get_cmap('terrain'),
        s=20, lw=0
    )
    for t in Ts.allTs():
        plt.annotate(
            f'{t.elevation:.1f}',
            xy=t.position,
            xytext=(7,0),
            textcoords='offset points',
            ha='left',
            va='center'
        )

# Show the hydrology network, if applicable
if args.showHydrology or args.showRiverFlow:
    print('Draw the hydrology network...')
    pos = [node.position for node in hydrology.allNodes()]
    labels = dict(zip(range(len(hydrology)),range(len(hydrology))))
    if args.showRiverFlow:
        normalizer = max([node.flow for node in hydrology.allNodes()])
        weights = [6*u.flow/normalizer for u,v in hydrology.allEdges()]
        nx.draw(hydrology.graph,pos,width=weights,node_size=1,labels=labels,ax=plt.gca())
    else:
        nx.draw(hydrology.graph,pos,node_size=1,labels=labels,ax=plt.gca())
    print('Drawing the voronoi cells...')
    voronoi_plot_2d(cells.vor, point_size=1, ax=plt.gca(),line_colors=['yellow'], show_vertices=False)

# Plot the rivers, if applicable
if args.showRivers:
    print('Plotting the river courses...')
    for mouth in hydrology.allMouthNodes():
        for leaf in hydrology.allLeaves(mouth.id):
            x = [coord[0] for coord in leaf.rivers[0].coords]
            y = [coord[1] for coord in leaf.rivers[0].coords]
            plt.plot(x,y)

plt.gca().set_ylim((upperY, lowerY))
plt.gca().set_xlim((lowerX, upperX))

plt.axis('on')
plt.tick_params(left=True,bottom=True,labelleft=True,labelbottom=True)
plt.tight_layout()
plt.show()
print('Saving the image...')
plt.savefig(args.outputFile)
print('All done')