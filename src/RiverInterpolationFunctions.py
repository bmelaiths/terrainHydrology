import numpy as np
from scipy import interpolate
from shapely.geometry import asLineString

from DataModel import HydroPrimitive, HydrologyNetwork, TerrainHoneycomb

def computeRivers(node: HydroPrimitive, hydrology: HydrologyNetwork, cells: TerrainHoneycomb):
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
