import typing

import numpy as np
from PIL import Image
from PIL import ImageDraw
from scipy import interpolate
from shapely.geometry import asLineString

from DataModel import ShoreModel, HydrologyNetwork, Q, TerrainHoneycomb, Terrain
import HydrologyFunctions

class RasterDataMock:
    def __init__(self, value: float) -> None:
        self.value = value
    def __getitem__(self, loc: typing.Tuple[float,float]) -> float:
        return self.value

def createTestImage():
    #create shore
    r = 100
    image = Image.new('L', (4*r,4*r))

    drawer = ImageDraw.Draw(image)

    #Keep in mind that these are _image_ coordinates
    drawer.polygon([         (150,134),  (100,200),(150,286),   (250,286),  (300,200),(250,134)], 255)
    #This should work out to (-100,132), (-200,0), (-100,-172), (100,-172), (200,0),  (100,132)

    image.save('imageFile.png')

def createTestObjects(): #-> typing.Tuple(ShoreModel, HydrologyFunctions.HydrologyParameters, HydrologyNetwork, TerrainHoneycomb):
    resolution = 2
    terrainSlopeRate = 1.0
    num_points = 20

    #create shore
    r = 100
    image = Image.new('L', (4*r,4*r))

    drawer = ImageDraw.Draw(image)

    #Keep in mind that these are _image_ coordinates
    drawer.polygon([         (150,134),  (100,200),(150,286),   (250,286),  (300,200),(250,134)], 255)
    #This should work out to (-100,132), (-200,0), (-100,-172), (100,-172), (200,0),  (100,132)

    image.save('imageFile.png')

    shore = ShoreModel(resolution, 'imageFile.png')

    #create candidate nodes
    hydrology = HydrologyNetwork()

    for i in range(0, 12):
        idx = int(42 * (i + 0.5))

        hydrology.addNode(shore[idx], 0, 1, contourIndex=idx)
    
    candidateNodes = hydrology.allMouthNodes()
    
    #create the parameters object
    Ps = 0.3
    Pa = 0
    Pc = 1 - (Ps + Pa)
    maxTries = 15
    riverAngleDev = 1.7
    edgelength = 50
    sigma = 0.75
    eta = 0.5
    zeta = 14
    riverSlope = RasterDataMock(100)
    slopeRate = 0.1
    params = HydrologyFunctions.HydrologyParameters(
        shore, hydrology, Pa, Pc, maxTries, riverAngleDev, edgelength,
        sigma, eta, zeta, riverSlope, slopeRate, candidateNodes
    )

    # This was used to generate the network for this test
    while len(candidateNodes) != 0:
        selectedCandiate = HydrologyFunctions.selectNode(candidateNodes, zeta)
        HydrologyFunctions.alpha(selectedCandiate, candidateNodes, params)

    cells = TerrainHoneycomb(shore, hydrology, resolution, edgelength)

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
        HydrologyFunctions.classify(hydrology.node(n), hydrology, edgelength)


    ## Calculate ridge elevations
    print('Calculating ridge elevations...')

    terrainSlope = RasterDataMock(128)

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
    Ts = Terrain(hydrology, cells, num_points)


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

    return shore, params, hydrology, cells

def hydrologyToCode(hydrology: HydrologyNetwork):
    code = ''
    for node in hydrology.allNodes():
        if node.parent is None:
            code += f'hydrology.addNode({node.position}, {node.elevation}, {node.priority}, contourIndex={node.contourIndex}) # ID: {node.id}\n'
        else:
            code += f'hydrology.addNode({node.position}, {node.elevation}, {node.priority}, parent=hydrology.node({node.parent.id})) # ID: {node.id}\n'
    return code

def terrainHoneycombToCode(cells: TerrainHoneycomb) -> str:
    code = ''

    # list
    i = 0
    code += 'cells.point_region = ['
    for datum in cells.point_region:
        code += f'{datum}'
        if i < len(cells.point_region) - 1:
            code += ','
        i += 1
    code += ']\n'

    # for a dictionary
    i = 0
    region_point_string = 'cells.region_point = {'
    for datum in cells.region_point:
        region_point_string += str(cells.region_point[datum]) + ' : ' + str(datum)
        if i < len(cells.region_point) - 1:
            region_point_string += ', '
        i += 1
    region_point_string += '}\n'
    code += region_point_string

    # for a list of lists
    i = 0
    region_string = 'cells.regions = ['
    for region in cells.regions:
        j = 0
        region_string += '['
        for point in region:
            region_string += str(point)
            if j < len(region) - 1:
                region_string += ','
                j += 1
        region_string += ']'
        if i < len(cells.regions) - 1:
            region_string += ','
        i += 1
    region_string += ']\n'
    code += region_string
    
    # list of coordinates
    i = 0
    code += 'cells.vertices = ['
    for coordinate in cells.vertices:
        code += f'[{coordinate[0]},{coordinate[1]}]'
        if i < len(cells.vertices) - 1:
            code += ','
        i += 1
    code += ']\n'

    # for a list of Qs
    i = 0
    q_string = 'cells.qs = [ ]\n'
    for q in cells.qs:
        if q is None:
            q_string += f'cells.qs.append(None) #{i}\n'
        else:
            # TODO: This needs to record positions as tuples, not lists
            q_string += f'cells.qs.append(Q([{q.position[0]},{q.position[1]}],{q.nodes},{q.vorIndex})) #{i}\n'
        i += 1
    code += q_string

    # for the cells Ridges
    cells_ridges_string = 'cells.cellsRidges = { }\n'
    for cellID in cells.cellsRidges:
        cells_ridges_string += f'cells.cellsRidges[{cellID}] = [ ]\n'
        for ridge in cells.cellsRidges[cellID]:
            if len(ridge) > 1:
                cells_ridges_string += f'cells.cellsRidges[{cellID}].append((cells.qs[{ridge[0].vorIndex}],cells.qs[{ridge[1].vorIndex}]))\n'
            else:
                cells_ridges_string += f'cells.cellsRidges[{cellID}].append((cells.qs[{ridge[0].vorIndex}],))\n'
    code += cells_ridges_string

    # for the cells Ridges
    cells_downstream_ridges_string = 'cells.cellsDownstreamRidges = { }\n'
    for cellID in cells.cellsDownstreamRidges:
        if cells.cellsDownstreamRidges[cellID] is not None:
            ridge = cells.cellsDownstreamRidges[cellID]
            #TODO this line needs to be fixed
            cells_downstream_ridges_string += f'cells.cellsDownstreamRidges[{cellID}] = ([{ridge[0][0]},{ridge[0][1]}],[{ridge[1][0]},{ridge[1][1]}])\n'
        else:
            cells_downstream_ridges_string += f'cells.cellsDownstreamRidges[{cellID}] = None\n'
    code += cells_downstream_ridges_string

    return code

def riversToCode(hydrology: HydrologyNetwork) -> str:
    code = ''

    for node in hydrology.allNodes():
        id = node.id
        for river in node.rivers:
            code += f'hydrology.node({id}).rivers.append(asLineString(' + str(list(river.coords)) + f'))\n'

    return code

def qElevationsToCode(cells: TerrainHoneycomb) -> str:
    code = ''
    
    i = 0
    for q in cells.allQs():
        if q is None:
            i += 1
            continue
        code += f'cells.qs[{i}].elevation = {q.elevation}\n'
        i += 1
    
    return code

def hydrologyAttributesToCode(hydrology: HydrologyNetwork) -> str:
    code = ''
    for node in hydrology.allNodes():
        id = node.id
        code += f'hydrology.node({id}).localWatershed = {node.localWatershed}\n'
        code += f'hydrology.node({id}).inheritedWatershed = {node.inheritedWatershed}\n'
        code += f'hydrology.node({id}).flow = {node.flow}\n'
    return code

def getPredefinedObjects0():
    edgelength = 2320.5
    resolution = 93.6
    r = 100
    image = Image.new('L', (4*r,4*r))

    drawer = ImageDraw.Draw(image)

    #Keep in mind that these are _image_ coordinates
    drawer.polygon([         (150,134),  (100,200),(150,286),   (250,286),  (300,200),(250,134)], 255)
    #This should work out to (-100,132), (-200,0), (-100,-172), (100,-172), (200,0),  (100,132)

    image.save('imageFile.png')

    shore = ShoreModel(resolution, '../in/test-example/in/gamma.png')

    hydrology = HydrologyNetwork()
    cells = TerrainHoneycomb()

    cells.shore = shore
    cells.resolution = resolution
    
    hydrology.addNode((-7768.799999999999, 2059.2), 0, 2, contourIndex=44) # ID: 0
    hydrology.addNode((-8049.599999999999, -2246.3999999999996), 0, 1, contourIndex=90) # ID: 1
    hydrology.addNode((-5054.4, -7394.4), 0, 1, contourIndex=145) # ID: 2
    hydrology.addNode((1123.1999999999998, -8049.599999999999), 0, 1, contourIndex=214) # ID: 3
    hydrology.addNode((4305.599999999999, -8049.599999999999), 0, 1, contourIndex=248) # ID: 4
    hydrology.addNode((6458.4, -5054.4), 0, 1, contourIndex=284) # ID: 5
    hydrology.addNode((8049.599999999999, 1684.8), 0, 1, contourIndex=356) # ID: 6
    hydrology.addNode((6832.799999999999, 3369.6), 0, 1, contourIndex=374) # ID: 7
    hydrology.addNode((280.79999999999995, 6177.599999999999), 0, 2, contourIndex=451) # ID: 8
    hydrology.addNode((-4867.2, 5990.4), 0, 1, contourIndex=2) # ID: 9
    hydrology.addNode((-6246.780372888135, 307.5788923724788), 173.81, 2, parent=hydrology.node(0)) # ID: 10
    hydrology.addNode((-5449.5362946522855, 2134.9371444985295), 173.81, 1, parent=hydrology.node(0)) # ID: 11
    hydrology.addNode((-5738.2285044404125, -2452.02601857411), 173.81, 1, parent=hydrology.node(1)) # ID: 12
    hydrology.addNode((-3779.8747892700185, -5455.249222671507), 173.81, 1, parent=hydrology.node(2)) # ID: 13
    hydrology.addNode((1735.4047436340666, -5811.313690817918), 173.81, 1, parent=hydrology.node(3)) # ID: 14
    hydrology.addNode((3913.3561082532797, -5762.491568073916), 173.81, 1, parent=hydrology.node(4)) # ID: 15
    hydrology.addNode((4575.801759157646, -3697.733455274554), 173.81, 1, parent=hydrology.node(5)) # ID: 16
    hydrology.addNode((6588.2148673450665, -117.71872224815615), 173.81, 1, parent=hydrology.node(6)) # ID: 17
    hydrology.addNode((4551.139609616782, 2946.8162338070397), 173.81, 1, parent=hydrology.node(7)) # ID: 18
    hydrology.addNode((1686.515368282502, 4331.337680237612), 173.81, 1, parent=hydrology.node(8)) # ID: 19
    hydrology.addNode((-267.90201010200553, 3922.9057071722514), 173.81, 1, parent=hydrology.node(8)) # ID: 20
    hydrology.addNode((-3628.3824225111284, 4028.245250826377), 173.81, 1, parent=hydrology.node(9)) # ID: 21
    hydrology.addNode((-3981.7104458665694, 811.7400528187368), 347.62, 1, parent=hydrology.node(10)) # ID: 22
    hydrology.addNode((-4397.990228017062, -1094.8102298855324), 347.62, 1, parent=hydrology.node(10)) # ID: 23
    hydrology.addNode((-3139.1312650010427, 2351.151965827316), 347.62, 1, parent=hydrology.node(11)) # ID: 24
    hydrology.addNode((-3652.2156918437145, -3468.52530321843), 347.62, 1, parent=hydrology.node(12)) # ID: 25
    hydrology.addNode((-1636.5946626095852, -4565.8277541525395), 347.62, 1, parent=hydrology.node(13)) # ID: 26
    hydrology.addNode((1544.4836554558808, -3498.6811242200897), 347.62, 1, parent=hydrology.node(14)) # ID: 27
    hydrology.addNode((4066.8172916595668, -1433.742496404423), 347.62, 1, parent=hydrology.node(16)) # ID: 28
    hydrology.addNode((4397.765121188957, 648.2121881900088), 347.62, 1, parent=hydrology.node(17)) # ID: 29
    hydrology.addNode((2306.434717613504, 2358.5814186043426), 347.62, 1, parent=hydrology.node(18)) # ID: 30
    hydrology.addNode((-1017.1741446275356, 1726.701818999854), 347.62, 1, parent=hydrology.node(20)) # ID: 31
    hydrology.addNode((-2307.913817105099, -795.4702929502098), 521.4300000000001, 1, parent=hydrology.node(22)) # ID: 32
    hydrology.addNode((-1496.566016258495, -2609.517313645959), 521.4300000000001, 1, parent=hydrology.node(25)) # ID: 33
    hydrology.addNode((1363.0351974719795, -1185.2860634716526), 521.4300000000001, 1, parent=hydrology.node(27)) # ID: 34
    hydrology.addNode((2670.0365109674985, 419.2884533342087), 521.4300000000001, 1, parent=hydrology.node(28)) # ID: 35
    hydrology.addNode((707.4833463640621, 676.8933493181478), 521.4300000000001, 1, parent=hydrology.node(30)) # ID: 36
    cells.point_region = [15,16,12,28,29,20,5,7,10,9,13,14,18,17,24,25,19,6,21,22,23,8,38,36,11,35,27,26,33,32,31,39,34,37,41,30,40,2,3,1,4]
    cells.region_point = {0 : 15, 1 : 16, 2 : 12, 3 : 28, 4 : 29, 5 : 20, 6 : 5, 7 : 7, 8 : 10, 9 : 9, 10 : 13, 11 : 14, 12 : 18, 13 : 17, 14 : 24, 15 : 25, 16 : 19, 17 : 6, 18 : 21, 19 : 22, 20 : 23, 21 : 8, 22 : 38, 23 : 36, 24 : 11, 25 : 35, 26 : 27, 27 : 26, 28 : 33, 29 : 32, 30 : 31, 31 : 39, 32 : 34, 33 : 37, 34 : 41, 35 : 30, 36 : 40, 37 : 2, 38 : 3, 39 : 1, 40 : 4}
    cells.regions = [[],[6,2,-1,1,5],[8,3,-1,0,7],[11,4,1,-1,0,10],[15,9,2,-1,3,14],[18,16,9,2,6],[17,25,26,16,18,19],[28,19,18,6,5,12,27],[31,32,34,33,30],[21,30,33,4,11],[34,20,12,5,1,4,33],[35,36,29,32,31],[7,8,13,38,37],[43,44,40,42,41,39],[39,41,35,31,30,21],[43,39,21,11,10],[7,37,44,43,10,0],[13,23,45,46,38],[37,38,46,47,40,44],[48,50,26,25,24],[49,15,9,16,26,50],[51,28,27,52],[53,52,27,12,20],[32,29,55,53,20,34],[60,62,63,59,61],[63,49,50,48,59],[58,61,59,48,24,22,57],[23,60,61,58,56,45],[8,3,14,62,60,23,13],[14,15,49,63,62],[66,68,67,65,64],[64,65,51,52,53,55,54],[67,17,19,28,51,65],[22,24,25,17,67,68],[69,70,75,74,73,72,71],[46,45,56,70,69,47],[47,69,71,42,40],[56,58,57,75,70],[42,71,72,36,35,41],[36,72,73,54,55,29],[73,74,66,64,54],[57,22,68,66,74,75]]
    cells.vertices = [[-46506.08407643312,-9.094947017729282e-13],[2.6147972675971687e-12,44226.7005988024],[46543.64331210191,-2.5011104298755527e-12],[5.9117155615240335e-12,-46570.47133757962],[-3525.521902377971,39972.85231539425],[13558.459960328704,28110.80657410031],[26619.009606328877,16377.840271237523],[-26698.104702750665,-16541.770008873114],[-5342.427753303965,-39560.66167400881],[35441.20183078683,-8340.111543380222],[-43498.34278442353,2227.431051158057],[-24953.50872483222,17779.580249280923],[5109.0255814912725,8395.459690146301],[-1996.7447595539256,-8015.650590079869],[2714.3999999999996,-43216.37197452229],[23612.9992884097,-19655.530738544476],[11641.06983801728,-2720.635933976302],[5208.635601476787,-547.9650058651541],[6239.082845541137,1659.004277335266],[5986.536369278422,1676.7167717523037],[543.7155555507849,4919.503757045385],[-6680.865690904744,4292.62943852493],[2628.517404889034,-2249.842814561697],[-1696.5804304393203,-7448.368722268813],[3110.058763294175,-2838.0481729036733],[6058.009213215691,-2175.2977391683517],[6801.766437182614,-2593.381713240389],[4914.451075025957,7354.28618418668],[5962.499774105598,1698.2241612392727],[-1737.0018899900447,3198.198620213414],[-5633.518041727106,4134.436097778672],[-4356.371628548453,2905.9619113121926],[-1959.7548342358916,3605.116544637308],[-2298.865028087399,6239.788272403492],[-1899.9120949673852,5514.184556143233],[-4268.698468751132,1969.1135115758668],[-2399.7332856176777,946.1571866549236],[-7144.421744184177,-5165.081742070794],[-6474.766974314268,-5072.428085853662],[-6591.153478267191,1545.413655907815],[-5811.891997043897,-1038.939525126176],[-5198.408619217473,937.7837124750154],[-4997.4197777103145,34.791145253456364],[-8437.747093365531,-59.129537389204245],[-6793.440632041436,-1219.8235314997867],[-2914.5426686927153,-4513.388678872556],[-5375.002315597173,-4355.289179137266],[-4425.902874448049,-2407.591207703676],[3011.7878956964087,-4334.590606933651],[5534.547354333574,-6661.643410927256],[5070.911181167803,-4995.22841079101],[3638.270932588358,1853.3084656616168],[3152.2050915292384,3708.1364069236884],[942.0084255460736,3013.6037341770443],[576.1703235384243,2402.732591397016],[574.891904979788,2409.4569155464214],[-2378.928797097079,-3529.5263151107347],[203.50001261420107,-2440.0462406337037],[-161.45820480964403,-3688.2482638015335],[2796.866017821202,-4559.486873042995],[-408.331712981741,-6427.83647054344],[202.65938098626384,-4773.653516381512],[2714.3999999999996,-7281.950244970516],[2854.2932737297983,-7121.312605780729],[1781.7422398560416,1256.4731251250837],[3397.988338982509,1559.506235989755],[1610.0097690595885,-51.86421349350661],[3644.2791292781067,-299.29466596039265],[2745.5022729159978,-976.7761938725419],[-3198.719658880232,-2022.0344053429003],[-2924.745798906439,-2159.8181666861387],[-3444.7439420729106,-304.2230231495334],[-2361.878510549001,823.5052256503409],[-854.634656375734,52.16237023542806],[-449.7788363341592,-776.981365945501],[-480.5358806262716,-1066.624705116723]]
    cells.qs = [ ]
    cells.qs.append(None) #0
    cells.qs.append(None) #1
    cells.qs.append(None) #2
    cells.qs.append(None) #3
    cells.qs.append(None) #4
    cells.qs.append(None) #5
    cells.qs.append(None) #6
    cells.qs.append(None) #7
    cells.qs.append(None) #8
    cells.qs.append(None) #9
    cells.qs.append(None) #10
    cells.qs.append(None) #11
    cells.qs.append(None) #12
    cells.qs.append(Q([-1996.7447595539256,-8015.650590079869],[13, 2, 3],13)) #13
    cells.qs.append(None) #14
    cells.qs.append(None) #15
    cells.qs.append(None) #16
    cells.qs.append(Q([5208.635601476787,-547.9650058651541],[28, 17, 29],17)) #17
    cells.qs.append(Q([6239.082845541137,1659.004277335266],[17, 6, 7],18)) #18
    cells.qs.append(Q([5986.536369278422,1676.7167717523037],[17, 29, 7],19)) #19
    cells.qs.append(Q([543.7155555507849,4919.503757045385],[20, 19, 8],20)) #20
    cells.qs.append(None) #21
    cells.qs.append(Q([2628.517404889034,-2249.842814561697],[27, 28, 34],22)) #22
    cells.qs.append(Q([-1696.5804304393203,-7448.368722268813],[13, 26, 3],23)) #23
    cells.qs.append(Q([3110.058763294175,-2838.0481729036733],[27, 16, 28],24)) #24
    cells.qs.append(Q([6058.009213215691,-2175.2977391683517],[16, 28, 17],25)) #25
    cells.qs.append(Q([6801.766437182614,-2593.381713240389],[5, 16, 17],26)) #26
    cells.qs.append(None) #27
    cells.qs.append(Q([5962.499774105598,1698.2241612392727],[29, 18, 7],28)) #28
    cells.qs.append(Q([-1737.0018899900447,3198.198620213414],[31, 24, 20],29)) #29
    cells.qs.append(Q([-5633.518041727106,4134.436097778672],[11, 9, 21],30)) #30
    cells.qs.append(Q([-4356.371628548453,2905.9619113121926],[11, 21, 24],31)) #31
    cells.qs.append(Q([-1959.7548342358916,3605.116544637308],[21, 24, 20],32)) #32
    cells.qs.append(None) #33
    cells.qs.append(Q([-1899.9120949673852,5514.184556143233],[21, 20, 8],34)) #34
    cells.qs.append(Q([-4268.698468751132,1969.1135115758668],[22, 11, 24],35)) #35
    cells.qs.append(Q([-2399.7332856176777,946.1571866549236],[22, 31, 24],36)) #36
    cells.qs.append(None) #37
    cells.qs.append(None) #38
    cells.qs.append(Q([-6591.153478267191,1545.413655907815],[10, 11, 0],39)) #39
    cells.qs.append(Q([-5811.891997043897,-1038.939525126176],[12, 23, 10],40)) #40
    cells.qs.append(Q([-5198.408619217473,937.7837124750154],[10, 22, 11],41)) #41
    cells.qs.append(Q([-4997.4197777103145,34.791145253456364],[23, 10, 22],42)) #42
    cells.qs.append(Q([-8437.747093365531,-59.129537389204245],[1, 10, 0],43)) #43
    cells.qs.append(Q([-6793.440632041436,-1219.8235314997867],[12, 1, 10],44)) #44
    cells.qs.append(Q([-2914.5426686927153,-4513.388678872556],[13, 26, 25],45)) #45
    cells.qs.append(Q([-5375.002315597173,-4355.289179137266],[12, 13, 25],46)) #46
    cells.qs.append(Q([-4425.902874448049,-2407.591207703676],[12, 25, 23],47)) #47
    cells.qs.append(Q([3011.7878956964087,-4334.590606933651],[27, 16, 15],48)) #48
    cells.qs.append(Q([5534.547354333574,-6661.643410927256],[4, 5, 15],49)) #49
    cells.qs.append(Q([5070.911181167803,-4995.22841079101],[5, 16, 15],50)) #50
    cells.qs.append(Q([3638.270932588358,1853.3084656616168],[29, 18, 30],51)) #51
    cells.qs.append(Q([3152.2050915292384,3708.1364069236884],[18, 19, 30],52)) #52
    cells.qs.append(Q([942.0084255460736,3013.6037341770443],[20, 19, 30],53)) #53
    cells.qs.append(Q([576.1703235384243,2402.732591397016],[31, 36, 30],54)) #54
    cells.qs.append(Q([574.891904979788,2409.4569155464214],[31, 20, 30],55)) #55
    cells.qs.append(Q([-2378.928797097079,-3529.5263151107347],[26, 33, 25],56)) #56
    cells.qs.append(Q([203.50001261420107,-2440.0462406337037],[33, 27, 34],57)) #57
    cells.qs.append(Q([-161.45820480964403,-3688.2482638015335],[26, 33, 27],58)) #58
    cells.qs.append(Q([2796.866017821202,-4559.486873042995],[27, 14, 15],59)) #59
    cells.qs.append(Q([-408.331712981741,-6427.83647054344],[26, 3, 14],60)) #60
    cells.qs.append(Q([202.65938098626384,-4773.653516381512],[26, 27, 14],61)) #61
    cells.qs.append(Q([2714.3999999999996,-7281.950244970516],[3, 4, 14],62)) #62
    cells.qs.append(Q([2854.2932737297983,-7121.312605780729],[4, 14, 15],63)) #63
    cells.qs.append(Q([1781.7422398560416,1256.4731251250837],[35, 36, 30],64)) #64
    cells.qs.append(Q([3397.988338982509,1559.506235989755],[35, 29, 30],65)) #65
    cells.qs.append(Q([1610.0097690595885,-51.86421349350661],[34, 35, 36],66)) #66
    cells.qs.append(Q([3644.2791292781067,-299.29466596039265],[28, 35, 29],67)) #67
    cells.qs.append(Q([2745.5022729159978,-976.7761938725419],[28, 34, 35],68)) #68
    cells.qs.append(Q([-3198.719658880232,-2022.0344053429003],[25, 23, 32],69)) #69
    cells.qs.append(Q([-2924.745798906439,-2159.8181666861387],[33, 25, 32],70)) #70
    cells.qs.append(Q([-3444.7439420729106,-304.2230231495334],[23, 32, 22],71)) #71
    cells.qs.append(Q([-2361.878510549001,823.5052256503409],[32, 22, 31],72)) #72
    cells.qs.append(Q([-854.634656375734,52.16237023542806],[32, 31, 36],73)) #73
    cells.qs.append(Q([-449.7788363341592,-776.981365945501],[32, 34, 36],74)) #74
    cells.qs.append(Q([-480.5358806262716,-1066.624705116723],[33, 32, 34],75)) #75
    cells.cellsRidges = { }
    cells.cellsRidges[8] = [ ]
    cells.cellsRidges[8].append((cells.qs[34],))
    cells.cellsRidges[8].append((cells.qs[20],))
    cells.cellsRidges[6] = [ ]
    cells.cellsRidges[6].append((cells.qs[18],))
    cells.cellsRidges[7] = [ ]
    cells.cellsRidges[7].append((cells.qs[18],cells.qs[19]))
    cells.cellsRidges[7].append((cells.qs[19],cells.qs[28]))
    cells.cellsRidges[1] = [ ]
    cells.cellsRidges[1].append((cells.qs[43],cells.qs[44]))
    cells.cellsRidges[3] = [ ]
    cells.cellsRidges[3].append((cells.qs[13],cells.qs[23]))
    cells.cellsRidges[3].append((cells.qs[23],cells.qs[60]))
    cells.cellsRidges[3].append((cells.qs[62],))
    cells.cellsRidges[2] = [ ]
    cells.cellsRidges[2].append((cells.qs[13],))
    cells.cellsRidges[9] = [ ]
    cells.cellsRidges[9].append((cells.qs[30],))
    cells.cellsRidges[0] = [ ]
    cells.cellsRidges[0].append((cells.qs[43],))
    cells.cellsRidges[0].append((cells.qs[39],))
    cells.cellsRidges[5] = [ ]
    cells.cellsRidges[5].append((cells.qs[49],cells.qs[50]))
    cells.cellsRidges[5].append((cells.qs[26],))
    cells.cellsRidges[4] = [ ]
    cells.cellsRidges[4].append((cells.qs[62],cells.qs[63]))
    cells.cellsRidges[4].append((cells.qs[49],))
    cells.cellsRidges[17] = [ ]
    cells.cellsRidges[17].append((cells.qs[18],cells.qs[19]))
    cells.cellsRidges[17].append((cells.qs[17],cells.qs[25]))
    cells.cellsRidges[17].append((cells.qs[25],cells.qs[26]))
    cells.cellsRidges[28] = [ ]
    cells.cellsRidges[28].append((cells.qs[17],cells.qs[25]))
    cells.cellsRidges[28].append((cells.qs[22],cells.qs[24]))
    cells.cellsRidges[28].append((cells.qs[17],cells.qs[67]))
    cells.cellsRidges[28].append((cells.qs[22],cells.qs[68]))
    cells.cellsRidges[16] = [ ]
    cells.cellsRidges[16].append((cells.qs[25],cells.qs[26]))
    cells.cellsRidges[16].append((cells.qs[24],cells.qs[48]))
    cells.cellsRidges[16].append((cells.qs[48],cells.qs[50]))
    cells.cellsRidges[19] = [ ]
    cells.cellsRidges[19].append((cells.qs[20],cells.qs[53]))
    cells.cellsRidges[19].append((cells.qs[52],cells.qs[53]))
    cells.cellsRidges[29] = [ ]
    cells.cellsRidges[29].append((cells.qs[19],cells.qs[28]))
    cells.cellsRidges[29].append((cells.qs[28],cells.qs[51]))
    cells.cellsRidges[29].append((cells.qs[65],cells.qs[67]))
    cells.cellsRidges[29].append((cells.qs[51],cells.qs[65]))
    cells.cellsRidges[29].append((cells.qs[17],cells.qs[67]))
    cells.cellsRidges[21] = [ ]
    cells.cellsRidges[21].append((cells.qs[30],cells.qs[31]))
    cells.cellsRidges[21].append((cells.qs[31],cells.qs[32]))
    cells.cellsRidges[21].append((cells.qs[32],cells.qs[34]))
    cells.cellsRidges[11] = [ ]
    cells.cellsRidges[11].append((cells.qs[30],cells.qs[31]))
    cells.cellsRidges[11].append((cells.qs[39],cells.qs[41]))
    cells.cellsRidges[11].append((cells.qs[35],cells.qs[41]))
    cells.cellsRidges[24] = [ ]
    cells.cellsRidges[24].append((cells.qs[31],cells.qs[32]))
    cells.cellsRidges[24].append((cells.qs[29],cells.qs[32]))
    cells.cellsRidges[24].append((cells.qs[29],cells.qs[36]))
    cells.cellsRidges[24].append((cells.qs[35],cells.qs[36]))
    cells.cellsRidges[20] = [ ]
    cells.cellsRidges[20].append((cells.qs[32],cells.qs[34]))
    cells.cellsRidges[20].append((cells.qs[29],cells.qs[32]))
    cells.cellsRidges[20].append((cells.qs[20],cells.qs[53]))
    cells.cellsRidges[20].append((cells.qs[53],cells.qs[55]))
    cells.cellsRidges[31] = [ ]
    cells.cellsRidges[31].append((cells.qs[29],cells.qs[36]))
    cells.cellsRidges[31].append((cells.qs[54],cells.qs[55]))
    cells.cellsRidges[31].append((cells.qs[72],cells.qs[73]))
    cells.cellsRidges[31].append((cells.qs[36],cells.qs[72]))
    cells.cellsRidges[31].append((cells.qs[54],cells.qs[73]))
    cells.cellsRidges[22] = [ ]
    cells.cellsRidges[22].append((cells.qs[35],cells.qs[36]))
    cells.cellsRidges[22].append((cells.qs[35],cells.qs[41]))
    cells.cellsRidges[22].append((cells.qs[42],cells.qs[71]))
    cells.cellsRidges[22].append((cells.qs[36],cells.qs[72]))
    cells.cellsRidges[12] = [ ]
    cells.cellsRidges[12].append((cells.qs[40],cells.qs[44]))
    cells.cellsRidges[12].append((cells.qs[40],cells.qs[47]))
    cells.cellsRidges[12].append((cells.qs[46],))
    cells.cellsRidges[10] = [ ]
    cells.cellsRidges[10].append((cells.qs[39],cells.qs[41]))
    cells.cellsRidges[10].append((cells.qs[40],cells.qs[44]))
    cells.cellsRidges[10].append((cells.qs[43],cells.qs[44]))
    cells.cellsRidges[10].append((cells.qs[42],))
    cells.cellsRidges[13] = [ ]
    cells.cellsRidges[13].append((cells.qs[13],cells.qs[23]))
    cells.cellsRidges[13].append((cells.qs[45],cells.qs[46]))
    cells.cellsRidges[25] = [ ]
    cells.cellsRidges[25].append((cells.qs[45],cells.qs[46]))
    cells.cellsRidges[25].append((cells.qs[45],cells.qs[56]))
    cells.cellsRidges[25].append((cells.qs[69],cells.qs[70]))
    cells.cellsRidges[25].append((cells.qs[47],cells.qs[69]))
    cells.cellsRidges[23] = [ ]
    cells.cellsRidges[23].append((cells.qs[40],cells.qs[47]))
    cells.cellsRidges[23].append((cells.qs[69],cells.qs[71]))
    cells.cellsRidges[23].append((cells.qs[47],cells.qs[69]))
    cells.cellsRidges[23].append((cells.qs[42],cells.qs[71]))
    cells.cellsRidges[27] = [ ]
    cells.cellsRidges[27].append((cells.qs[24],cells.qs[48]))
    cells.cellsRidges[27].append((cells.qs[48],cells.qs[59]))
    cells.cellsRidges[27].append((cells.qs[22],cells.qs[24]))
    cells.cellsRidges[27].append((cells.qs[57],cells.qs[58]))
    cells.cellsRidges[27].append((cells.qs[58],cells.qs[61]))
    cells.cellsRidges[15] = [ ]
    cells.cellsRidges[15].append((cells.qs[48],cells.qs[50]))
    cells.cellsRidges[15].append((cells.qs[49],cells.qs[50]))
    cells.cellsRidges[15].append((cells.qs[59],cells.qs[63]))
    cells.cellsRidges[15].append((cells.qs[48],cells.qs[59]))
    cells.cellsRidges[18] = [ ]
    cells.cellsRidges[18].append((cells.qs[28],cells.qs[51]))
    cells.cellsRidges[18].append((cells.qs[52],))
    cells.cellsRidges[30] = [ ]
    cells.cellsRidges[30].append((cells.qs[52],cells.qs[53]))
    cells.cellsRidges[30].append((cells.qs[53],cells.qs[55]))
    cells.cellsRidges[30].append((cells.qs[64],cells.qs[65]))
    cells.cellsRidges[30].append((cells.qs[54],cells.qs[55]))
    cells.cellsRidges[30].append((cells.qs[51],cells.qs[65]))
    cells.cellsRidges[14] = [ ]
    cells.cellsRidges[14].append((cells.qs[59],cells.qs[63]))
    cells.cellsRidges[14].append((cells.qs[60],cells.qs[61]))
    cells.cellsRidges[14].append((cells.qs[62],cells.qs[63]))
    cells.cellsRidges[26] = [ ]
    cells.cellsRidges[26].append((cells.qs[60],cells.qs[61]))
    cells.cellsRidges[26].append((cells.qs[58],cells.qs[61]))
    cells.cellsRidges[26].append((cells.qs[23],cells.qs[60]))
    cells.cellsRidges[26].append((cells.qs[56],cells.qs[58]))
    cells.cellsRidges[26].append((cells.qs[45],cells.qs[56]))
    cells.cellsRidges[33] = [ ]
    cells.cellsRidges[33].append((cells.qs[57],cells.qs[58]))
    cells.cellsRidges[33].append((cells.qs[56],cells.qs[58]))
    cells.cellsRidges[33].append((cells.qs[70],cells.qs[75]))
    cells.cellsRidges[33].append((cells.qs[57],cells.qs[75]))
    cells.cellsRidges[35] = [ ]
    cells.cellsRidges[35].append((cells.qs[64],cells.qs[65]))
    cells.cellsRidges[35].append((cells.qs[64],cells.qs[66]))
    cells.cellsRidges[35].append((cells.qs[65],cells.qs[67]))
    cells.cellsRidges[35].append((cells.qs[66],cells.qs[68]))
    cells.cellsRidges[36] = [ ]
    cells.cellsRidges[36].append((cells.qs[64],cells.qs[66]))
    cells.cellsRidges[36].append((cells.qs[73],cells.qs[74]))
    cells.cellsRidges[36].append((cells.qs[54],cells.qs[73]))
    cells.cellsRidges[36].append((cells.qs[66],cells.qs[74]))
    cells.cellsRidges[34] = [ ]
    cells.cellsRidges[34].append((cells.qs[66],cells.qs[68]))
    cells.cellsRidges[34].append((cells.qs[22],cells.qs[68]))
    cells.cellsRidges[34].append((cells.qs[74],cells.qs[75]))
    cells.cellsRidges[34].append((cells.qs[57],cells.qs[75]))
    cells.cellsRidges[34].append((cells.qs[66],cells.qs[74]))
    cells.cellsRidges[32] = [ ]
    cells.cellsRidges[32].append((cells.qs[69],cells.qs[70]))
    cells.cellsRidges[32].append((cells.qs[69],cells.qs[71]))
    cells.cellsRidges[32].append((cells.qs[70],cells.qs[75]))
    cells.cellsRidges[32].append((cells.qs[72],cells.qs[73]))
    cells.cellsRidges[32].append((cells.qs[73],cells.qs[74]))
    cells.cellsRidges[32].append((cells.qs[74],cells.qs[75]))
    cells.cellsDownstreamRidges = { }
    cells.cellsDownstreamRidges[17] = None
    cells.cellsDownstreamRidges[29] = ([5208.635601476787,-547.9650058651541],[5986.536369278422,1676.7167717523037])
    cells.cellsDownstreamRidges[18] = None
    cells.cellsDownstreamRidges[21] = None
    cells.cellsDownstreamRidges[19] = None
    cells.cellsDownstreamRidges[20] = ([543.7155555507849,4919.503757045385],[-1899.9120949673852,5514.184556143233])
    cells.cellsDownstreamRidges[24] = ([-4356.371628548453,2905.9619113121926],[-4268.698468751132,1969.1135115758668])
    cells.cellsDownstreamRidges[13] = None
    cells.cellsDownstreamRidges[10] = ([-6591.153478267191,1545.413655907815],[-8437.747093365531,-59.129537389204245])
    cells.cellsDownstreamRidges[23] = ([-5811.891997043897,-1038.939525126176],[-4997.4197777103145,34.791145253456364])
    cells.cellsDownstreamRidges[22] = ([-5198.408619217473,937.7837124750154],[-4997.4197777103145,34.791145253456364])
    cells.cellsDownstreamRidges[11] = None
    cells.cellsDownstreamRidges[12] = None
    cells.cellsDownstreamRidges[26] = ([-1696.5804304393203,-7448.368722268813],[-2914.5426686927153,-4513.388678872556])
    cells.cellsDownstreamRidges[25] = ([-5375.002315597173,-4355.289179137266],[-4425.902874448049,-2407.591207703676])
    cells.cellsDownstreamRidges[28] = ([3110.058763294175,-2838.0481729036733],[6058.009213215691,-2175.2977391683517])
    cells.cellsDownstreamRidges[16] = ([6801.766437182614,-2593.381713240389],[5070.911181167803,-4995.22841079101])
    cells.cellsDownstreamRidges[30] = ([3638.270932588358,1853.3084656616168],[3152.2050915292384,3708.1364069236884])
    cells.cellsDownstreamRidges[31] = ([-1737.0018899900447,3198.198620213414],[574.891904979788,2409.4569155464214])
    cells.cellsDownstreamRidges[27] = ([2796.866017821202,-4559.486873042995],[202.65938098626384,-4773.653516381512])
    cells.cellsDownstreamRidges[14] = ([-408.331712981741,-6427.83647054344],[2714.3999999999996,-7281.950244970516])
    cells.cellsDownstreamRidges[15] = ([5534.547354333574,-6661.643410927256],[2854.2932737297983,-7121.312605780729])
    cells.cellsDownstreamRidges[34] = ([2628.517404889034,-2249.842814561697],[203.50001261420107,-2440.0462406337037])
    cells.cellsDownstreamRidges[35] = ([3644.2791292781067,-299.29466596039265],[2745.5022729159978,-976.7761938725419])
    cells.cellsDownstreamRidges[36] = ([576.1703235384243,2402.732591397016],[1781.7422398560416,1256.4731251250837])
    cells.cellsDownstreamRidges[32] = ([-3444.7439420729106,-304.2230231495334],[-2361.878510549001,823.5052256503409])
    cells.cellsDownstreamRidges[33] = ([-2378.928797097079,-3529.5263151107347],[-2924.745798906439,-2159.8181666861387])
    cells.qs[13].elevation = 2510.848146702606
    cells.qs[17].elevation = 1430.0377736754585
    cells.qs[18].elevation = 1530.0604856168566
    cells.qs[19].elevation = 1765.230537522762
    cells.qs[20].elevation = 1136.505861851151
    cells.qs[22].elevation = 1760.0843151172155
    cells.qs[23].elevation = 2507.1671576529934
    cells.qs[24].elevation = 1620.3939383676675
    cells.qs[25].elevation = 1939.1323517252927
    cells.qs[26].elevation = 2035.0161948093946
    cells.qs[28].elevation = 1759.0638674535917
    cells.qs[29].elevation = 1574.6080430713005
    cells.qs[30].elevation = 1677.8005913522338
    cells.qs[31].elevation = 1349.5966564187677
    cells.qs[32].elevation = 1637.0123695839943
    cells.qs[34].elevation = 1881.1188752048213
    cells.qs[35].elevation = 1240.769098993515
    cells.qs[36].elevation = 1536.8215447968705
    cells.qs[39].elevation = 1136.1844080307794
    cells.qs[40].elevation = 1407.4866484907493
    cells.qs[41].elevation = 1263.8278902897393
    cells.qs[42].elevation = 1305.4620607560596
    cells.qs[43].elevation = 1837.7145738312977
    cells.qs[44].elevation = 1388.930269295892
    cells.qs[45].elevation = 1305.633636270446
    cells.qs[46].elevation = 1798.930087826579
    cells.qs[47].elevation = 1331.1409340287785
    cells.qs[48].elevation = 1612.4935604221987
    cells.qs[49].elevation = 1562.3750041782946
    cells.qs[50].elevation = 1214.0109633632885
    cells.qs[51].elevation = 1414.5690408196588
    cells.qs[52].elevation = 1540.5677517601375
    cells.qs[53].elevation = 1481.2688696412329
    cells.qs[54].elevation = 1817.853814207954
    cells.qs[55].elevation = 1645.13921626588
    cells.qs[56].elevation = 1476.2407112661158
    cells.qs[57].elevation = 1801.1240150880571
    cells.qs[58].elevation = 1807.0787555418087
    cells.qs[59].elevation = 1576.9641431819432
    cells.qs[60].elevation = 2018.4049591744592
    cells.qs[61].elevation = 1734.0241320286664
    cells.qs[62].elevation = 1497.0975772399393
    cells.qs[63].elevation = 1464.2142110213354
    cells.qs[64].elevation = 1435.7080405738143
    cells.qs[65].elevation = 1534.687839094508
    cells.qs[66].elevation = 1390.3061527744826
    cells.qs[67].elevation = 1428.1801498659543
    cells.qs[68].elevation = 1568.6364487605458
    cells.qs[69].elevation = 1656.8791403642763
    cells.qs[70].elevation = 1642.9419423078643
    cells.qs[71].elevation = 1449.0376719873552
    cells.qs[72].elevation = 1734.7478834797996
    cells.qs[73].elevation = 1781.5876427099447
    cells.qs[74].elevation = 1913.2784313448424
    cells.qs[75].elevation = 1905.1582977178052
    hydrology.node(0).localWatershed = 0.0
    hydrology.node(0).inheritedWatershed = 28131818.648611713
    hydrology.node(0).flow = 57968.99600441555
    hydrology.node(1).localWatershed = 0.0
    hydrology.node(1).inheritedWatershed = 14209117.884563802
    hydrology.node(1).flow = 36184.33548755921
    hydrology.node(2).localWatershed = 0.0
    hydrology.node(2).inheritedWatershed = 12969365.474605657
    hydrology.node(2).flow = 33975.29319657288
    hydrology.node(3).localWatershed = 3369737.8639595695
    hydrology.node(3).inheritedWatershed = 22453341.850410018
    hydrology.node(3).flow = 49617.30510437502
    hydrology.node(4).localWatershed = 183122.52971809474
    hydrology.node(4).inheritedWatershed = 5725769.280893651
    hydrology.node(4).flow = 19326.563717498964
    hydrology.node(5).localWatershed = 1998953.0865354785
    hydrology.node(5).inheritedWatershed = 18381498.51244879
    hydrology.node(5).flow = 43218.741926492454
    hydrology.node(6).localWatershed = 0.0
    hydrology.node(6).inheritedWatershed = 8173623.658922188
    hydrology.node(6).flow = 24706.736886128343
    hydrology.node(7).localWatershed = 2502.9336853702553
    hydrology.node(7).inheritedWatershed = 11859976.094863027
    hydrology.node(7).flow = 31942.368959858828
    hydrology.node(8).localWatershed = 0.0
    hydrology.node(8).inheritedWatershed = 14018798.502285188
    hydrology.node(8).flow = 35849.222236910085
    hydrology.node(9).localWatershed = 0.0
    hydrology.node(9).inheritedWatershed = 5441122.226075337
    hydrology.node(9).flow = 18658.39757381917
    hydrology.node(10).localWatershed = 5768213.687823514
    hydrology.node(10).inheritedWatershed = 19534470.370320205
    hydrology.node(10).flow = 45071.5426308281
    hydrology.node(11).localWatershed = 4348612.752602762
    hydrology.node(11).inheritedWatershed = 8597348.278291509
    hydrology.node(11).flow = 25583.550115525013
    hydrology.node(12).localWatershed = 3666329.790781794
    hydrology.node(12).inheritedWatershed = 14209117.884563802
    hydrology.node(12).flow = 36184.33548755921
    hydrology.node(13).localWatershed = 5021987.462126826
    hydrology.node(13).inheritedWatershed = 12969365.474605657
    hydrology.node(13).flow = 33975.29319657288
    hydrology.node(14).localWatershed = 6550001.862252988
    hydrology.node(14).inheritedWatershed = 19083603.98645045
    hydrology.node(14).flow = 44351.157977762596
    hydrology.node(15).localWatershed = 5542646.751175556
    hydrology.node(15).inheritedWatershed = 5542646.751175556
    hydrology.node(15).flow = 18897.926643135215
    hydrology.node(16).localWatershed = 6657668.371590553
    hydrology.node(16).inheritedWatershed = 16382545.425913312
    hydrology.node(16).flow = 39918.33726878896
    hydrology.node(17).localWatershed = 3527264.9577371967
    hydrology.node(17).inheritedWatershed = 8173623.658922188
    hydrology.node(17).flow = 24706.736886128343
    hydrology.node(18).localWatershed = 2117831.7071346184
    hydrology.node(18).inheritedWatershed = 11857473.161177658
    hydrology.node(18).flow = 31937.717428538777
    hydrology.node(19).localWatershed = 2244520.643887302
    hydrology.node(19).inheritedWatershed = 2244520.643887302
    hydrology.node(19).flow = 10128.020629321627
    hydrology.node(20).localWatershed = 6365407.845070189
    hydrology.node(20).inheritedWatershed = 11774277.858397886
    hydrology.node(20).flow = 31782.93091397958
    hydrology.node(21).localWatershed = 5441122.226075337
    hydrology.node(21).inheritedWatershed = 5441122.226075337
    hydrology.node(21).flow = 18658.39757381917
    hydrology.node(22).localWatershed = 3939326.1942312163
    hydrology.node(22).inheritedWatershed = 9564102.098520072
    hydrology.node(22).flow = 27535.552240048903
    hydrology.node(23).localWatershed = 4202154.583976618
    hydrology.node(23).inheritedWatershed = 4202154.583976618
    hydrology.node(23).flow = 15611.51272658518
    hydrology.node(24).localWatershed = 4248735.525688747
    hydrology.node(24).inheritedWatershed = 4248735.525688747
    hydrology.node(24).flow = 15730.715696602518
    hydrology.node(25).localWatershed = 5011626.696132958
    hydrology.node(25).inheritedWatershed = 10542788.093782008
    hydrology.node(25).flow = 29450.221791461456
    hydrology.node(26).localWatershed = 7947378.01247883
    hydrology.node(26).inheritedWatershed = 7947378.01247883
    hydrology.node(26).flow = 24232.809098167974
    hydrology.node(27).localWatershed = 6924975.677813589
    hydrology.node(27).inheritedWatershed = 12533602.12419746
    hydrology.node(27).flow = 33183.45975194562
    hydrology.node(28).localWatershed = 6033136.527569512
    hydrology.node(28).inheritedWatershed = 9724877.054322759
    hydrology.node(28).flow = 27854.112982808685
    hydrology.node(29).localWatershed = 4646358.701184991
    hydrology.node(29).inheritedWatershed = 4646358.701184991
    hydrology.node(29).flow = 16732.3552714103
    hydrology.node(30).localWatershed = 4880773.9598184135
    hydrology.node(30).inheritedWatershed = 9739641.454043038
    hydrology.node(30).flow = 27883.285100363544
    hydrology.node(31).localWatershed = 5408870.013327697
    hydrology.node(31).inheritedWatershed = 5408870.013327697
    hydrology.node(31).flow = 18582.01499218459
    hydrology.node(32).localWatershed = 5624775.904288855
    hydrology.node(32).inheritedWatershed = 5624775.904288855
    hydrology.node(32).flow = 19090.701859267556
    hydrology.node(33).localWatershed = 5531161.39764905
    hydrology.node(33).inheritedWatershed = 5531161.39764905
    hydrology.node(33).flow = 18870.89764112951
    hydrology.node(34).localWatershed = 5608626.446383872
    hydrology.node(34).inheritedWatershed = 5608626.446383872
    hydrology.node(34).flow = 19052.864816850353
    hydrology.node(35).localWatershed = 3691740.526753247
    hydrology.node(35).inheritedWatershed = 3691740.526753247
    hydrology.node(35).flow = 14277.060408616415
    hydrology.node(36).localWatershed = 4858867.494224624
    hydrology.node(36).inheritedWatershed = 4858867.494224624
    hydrology.node(36).flow = 17256.72899144959
    hydrology.node(0).rivers.append(asLineString([(-2307.9138171051, -795.4702929502099, 521.0), (-2410.2467404855824, -468.92229961709603, 494.0), (-2589.329360063512, -134.91271958570707, 467.0), (-2831.3001169335002, 182.3566023926544, 441.0), (-3122.2974521901565, 458.683821566686, 415.0), (-3448.459806928091, 669.867093185086, 388.0), (-3795.9256222419135, 791.7045724965521, 361.0), (-4151.091001563071, 801.3935929862843, 334.0), (-4506.991045729097, 712.1830927193814, 307.0), (-4863.721772535623, 575.6647549940313, 279.0), (-5221.7447113457165, 444.94402147870636, 250.0), (-5582.80921091596, 354.15517661449695, 222.0), (-5950.597772563143, 308.95479954551905, 194.0), (-6328.904136334305, 312.6569890399968, 168.0), (-6713.000142958087, 368.9172563554919, 143.0), (-7079.890150895739, 482.1229621134004, 119.0), (-7404.218439853203, 656.7560186697469, 95.0), (-7660.629289536417, 897.2983383805561, 72.0), (-7823.766979651329, 1208.2318336018511, 49.0), (-7868.275789903877, 1594.0384166896583, 25.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(0).rivers.append(asLineString([(-3139.1312650010436, 2351.1519658273155, 347.0), (-3377.346916022355, 2419.7214594739, 330.0), (-3612.6330740757653, 2459.4139283808986, 313.0), (-3845.3598789669527, 2473.5831382165215, 295.0), (-4075.897470501597, 2465.582854648985, 278.0), (-4304.615988485376, 2438.766843346499, 261.0), (-4531.885572723969, 2396.488869977277, 244.0), (-4758.076363023056, 2342.1027002095316, 226.0), (-4983.558499188315, 2278.9620997114757, 209.0), (-5208.702121025426, 2210.420834151323, 192.0), (-5433.877368340067, 2139.832669197284, 175.0), (-5659.454380937916, 2070.5513705175726, 157.0), (-5885.8032986246535, 2005.930703780401, 140.0), (-6113.294261205958, 1949.324434653983, 122.0), (-6342.297408487509, 1904.08632880653, 105.0), (-6573.182880274984, 1873.5701519062554, 88.0), (-6806.320816374065, 1861.1296696213717, 70.0), (-7042.081356590426, 1870.1186476200917, 53.0), (-7280.83464072975, 1903.890851570628, 35.0), (-7522.950808597715, 1965.800047141193, 17.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(1).rivers.append(asLineString([(-1496.566016258495, -2609.517313645959, 521.0), (-1902.0427869123275, -2564.1902269704024, 495.0), (-2252.993566029566, -2636.569132844976, 468.0), (-2567.0038080586255, -2790.8980761702683, 442.0), (-2861.6589674479214, -2991.421101846868, 415.0), (-3154.5444986458688, -3202.3822547753643, 388.0), (-3463.245856100884, -3388.0255798563458, 362.0), (-3804.701652299752, -3513.0509265168816, 336.0), (-4174.917221672106, -3556.909016048467, 311.0), (-4544.918186536708, -3516.652807268503, 286.0), (-4884.26610902004, -3390.366925618178, 261.0), (-5169.5825971566965, -3180.700487378298, 237.0), (-5410.488617232776, -2911.6374939582925, 212.0), (-5626.238695793359, -2613.390277025972, 187.0), (-5836.08851279935, -2316.171913960441, 162.0), (-6059.293748211665, -2050.1954821408135, 136.0), (-6315.110081991214, -1845.6740589461945, 110.0), (-6622.793194098913, -1732.820721755695, 83.0), (-7001.598764495668, -1741.8485479484234, 56.0), (-7470.782473142393, -1902.9706149034887, 28.0), (-8049.599999999999, -2246.3999999999996, 0.0)]))
    hydrology.node(2).rivers.append(asLineString([(-1636.5946626095858, -4565.827754152543, 347.0), (-1646.1730635493357, -5080.695788498145, 331.0), (-1708.2719073314167, -5462.95259193207, 316.0), (-1816.967819551049, -5728.148467464976, 300.0), (-1966.3374258034555, -5891.8337181075285, 286.0), (-2150.4573516838564, -5969.558646870382, 271.0), (-2363.4042227874747, -5976.873556764207, 256.0), (-2599.2546647095314, -5929.32875079966, 242.0), (-2852.085303045248, -5842.474531987407, 227.0), (-3115.9727633898465, -5731.861203338108, 212.0), (-3384.9936713385478, -5613.039067862422, 197.0), (-3653.224652486575, -5501.558428571012, 181.0), (-3914.7423324291485, -5412.969588474543, 165.0), (-4163.623336761489, -5362.822850583673, 148.0), (-4393.944291078819, -5366.668517909067, 130.0), (-4599.781820976363, -5440.056893461385, 111.0), (-4775.212552049337, -5598.538280251289, 91.0), (-4914.313109892967, -5857.662981289441, 70.0), (-5011.160120102473, -6232.981299586501, 48.0), (-5059.830208273076, -6740.043538153134, 25.0), (-5054.4, -7394.4, 0.0)]))
    hydrology.node(3).rivers.append(asLineString([(1363.0351974719797, -1185.286063471653, 521.0), (1336.143255634835, -1538.7001375219454, 495.0), (1353.7876112410288, -1890.806818069666, 468.0), (1399.585961103793, -2241.938110932066, 442.0), (1457.1560020363577, -2592.426021926398, 415.0), (1510.115430851954, -2942.6025568699124, 389.0), (1542.0819443638127, -3292.799721579862, 363.0), (1537.332270095805, -3643.3501728819615, 336.0), (1503.966575399799, -3994.6101010438947, 310.0), (1480.0316243673785, -4346.965278396289, 284.0), (1505.4616936351763, -4700.802289486349, 258.0), (1596.4084750411462, -5054.946165214905, 231.0), (1697.6362068321528, -5403.534668912686, 205.0), (1740.5704612390641, -5739.829753197196, 179.0), (1666.514067439315, -6058.414368534419, 153.0), (1494.5032524558872, -6364.2676362974225, 127.0), (1280.3573891727174, -6667.288106077629, 102.0), (1080.115304937674, -6977.403677606516, 76.0), (949.8158270986277, -7304.5422506155555, 51.0), (945.4977830034463, -7658.631724836225, 25.0), (1123.1999999999998, -8049.599999999999, 0.0)]))
    hydrology.node(4).rivers.append(asLineString([(3913.3561082532806, -5762.491568073917, 173.0), (3941.671617470887, -5875.354352961782, 165.0), (3969.6734937288447, -5988.2709265598605, 156.0), (3997.570825666919, -6101.205429728012, 147.0), (4025.572701924877, -6214.122003326092, 139.0), (4053.888211142483, -6326.984788213957, 130.0), (4082.621897639623, -6439.775854821535, 121.0), (4111.460128456644, -6552.548991859043, 112.0), (4139.984726314018, -6665.375917606766, 104.0), (4167.77751393221, -6778.328350344989, 95.0), (4194.420314031686, -6891.478008353995, 86.0), (4219.494949332916, -7004.896609914071, 78.0), (4242.583242556366, -7118.655873305503, 69.0), (4263.267016422504, -7232.827516808577, 60.0), (4281.128093651797, -7347.483258703574, 52.0), (4295.7482969647135, -7462.694817270785, 43.0), (4306.70944908172, -7578.533910790492, 34.0), (4313.593372723284, -7695.072257542983, 26.0), (4315.981890609875, -7812.381575808544, 17.0), (4313.456825461959, -7930.533583867455, 8.0), (4305.599999999999, -8049.599999999999, 0.0)]))
    hydrology.node(5).rivers.append(asLineString([(2670.0365109674985, 419.2884533342088, 521.0), (2758.264111976863, 21.011373584195212, 494.0), (2926.651794441099, -311.0764901513366, 466.0), (3154.1793564607824, -595.4881549065561, 438.0), (3419.8265961364887, -850.7366377156324, 411.0), (3702.5733115687926, -1095.334955612735, 383.0), (3981.399300858273, -1347.796125632034, 356.0), (4234.90454593771, -1626.0646861888808, 329.0), (4437.4543218884455, -1941.7470034297003, 302.0), (4560.7652406984325, -2302.4851353310232, 274.0), (4579.093414254613, -2714.2145129507085, 247.0), (4529.222011255771, -3142.9813144125806, 219.0), (4512.284236539314, -3513.7794918325594, 192.0), (4631.555599532296, -3750.4087423165456, 166.0), (4928.082793724015, -3829.7399568500273, 141.0), (5337.348556229627, -3818.67314092122, 117.0), (5784.599861521532, -3792.8377629774704, 94.0), (6195.083684072124, -3827.8632914661257, 72.0), (6494.046998353801, -3999.379194834534, 48.0), (6606.73677883896, -4383.014941530043, 25.0), (6458.4, -5054.4, 0.0)]))
    hydrology.node(6).rivers.append(asLineString([(4397.765121188958, 648.2121881900085, 347.0), (4651.83421662707, 740.611340037247, 330.0), (4896.645099683202, 769.6176426668867, 313.0), (5132.486817131576, 745.0848207510403, 296.0), (5359.648415746434, 676.8665989618216, 279.0), (5578.418942302003, 574.8167019713433, 262.0), (5789.087443572522, 448.78885445171846, 245.0), (5991.942966332221, 308.636781075061, 228.0), (6187.274557355335, 164.2142065134836, 211.0), (6375.371263416098, 25.374855439099605, 194.0), (6556.522131288741, -98.02754747597797, 176.0), (6731.0162077475, -196.13927755963596, 159.0), (6899.1425395666065, -259.10661013976136, 142.0), (7061.190173520294, -277.07582054424023, 125.0), (7217.448156382796, -240.19318410096025, 107.0), (7368.205534928347, -138.60497613780808, 89.0), (7513.751355931181, 37.54252801733003, 72.0), (7654.37466616553, 298.10305303656673, 54.0), (7790.364512405627, 652.9303235920152, 36.0), (7922.009941425705, 1111.8780643557889, 18.0), (8049.599999999999, 1684.8, 0.0)]))
    hydrology.node(7).rivers.append(asLineString([(707.4833463640622, 676.893349318148, 521.0), (728.42915148112, 1116.62928459072, 495.0), (850.9006452244694, 1459.65755506244, 470.0), (1055.781178299755, 1724.2255016738395, 446.0), (1323.9541014126216, 1928.580465365449, 421.0), (1636.3027652687138, 2090.9697870777986, 397.0), (1973.7105205736757, 2229.6408077514193, 372.0), (2317.0608532457245, 2362.8406987258013, 346.0), (2652.1008268062815, 2502.7161053142636, 320.0), (2982.293269900961, 2639.192276353892, 293.0), (3315.1335274120825, 2757.13635895912, 267.0), (3657.350771781155, 2843.220282797769, 240.0), (4009.286300147806, 2899.1631461433994, 213.0), (4368.107623518778, 2934.1601630049645, 187.0), (4730.959275855543, 2957.460671721965, 160.0), (5094.985791119567, 2978.314010633898, 134.0), (5457.331703272321, 3005.9695180802623, 107.0), (5815.141546275272, 3049.6765324005564, 81.0), (6165.55985408989, 3118.684391934279, 54.0), (6505.731160677642, 3222.242435020927, 27.0), (6832.799999999999, 3369.6, 0.0)]))
    hydrology.node(8).rivers.append(asLineString([(1686.5153682825016, 4331.337680237611, 173.0), (1616.2295998683765, 4423.650796225731, 165.0), (1545.9438314542517, 4515.9639122138515, 156.0), (1475.6580630401263, 4608.277028201971, 147.0), (1405.3722946260018, 4700.590144190091, 139.0), (1335.0865262118766, 4792.90326017821, 130.0), (1264.8007577977514, 4885.21637616633, 121.0), (1194.514989383626, 4977.529492154448, 112.0), (1124.2292209695013, 5069.842608142568, 104.0), (1053.9434525553763, 5162.155724130688, 95.0), (983.657684141251, 5254.468840118806, 86.0), (913.371915727126, 5346.781956106925, 78.0), (843.0861473130008, 5439.095072095044, 69.0), (772.8003788988756, 5531.408188083163, 60.0), (702.5146104847504, 5623.721304071281, 52.0), (632.2288420706254, 5716.034420059401, 43.0), (561.9430736565002, 5808.34753604752, 34.0), (491.6573052423751, 5900.66065203564, 26.0), (421.3715368282501, 5992.9737680237595, 17.0), (351.0857684141249, 6085.286884011879, 8.0), (280.79999999999995, 6177.599999999999, 0.0)]))
    hydrology.node(8).rivers.append(asLineString([(-1017.1741446275357, 1726.7018189998541, 347.0), (-964.6526348680636, 1956.5615396404905, 329.0), (-879.6740056147408, 2188.403951736881, 311.0), (-772.8167213108993, 2422.225661544413, 291.0), (-654.6592463998693, 2658.0232753184723, 272.0), (-535.7800453249822, 2895.7933993144434, 253.0), (-426.75758252956894, 3135.532639787715, 233.0), (-338.1703224569609, 3377.2376029936713, 214.0), (-280.59672955048904, 3620.9048951877, 195.0), (-264.61526825348403, 3866.5311226251847, 177.0), (-298.54746294559914, 4113.922071012762, 160.0), (-372.42535190600745, 4361.337180263369, 143.0), (-467.40858016649497, 4606.285744118389, 127.0), (-564.5982592182844, 4846.272107405172, 111.0), (-645.0955005525993, 5078.80061495108, 96.0), (-690.0014156606621, 5301.375611583464, 80.0), (-680.417116033696, 5511.501442129682, 65.0), (-597.4437131629237, 5706.682451417089, 49.0), (-422.1823185395683, 5884.42298427304, 33.0), (-135.73404365485237, 6042.227385524893, 17.0), (280.79999999999995, 6177.599999999999, 0.0)]))
    hydrology.node(9).rivers.append(asLineString([(-3628.3824225111284, 4028.2452508263764, 173.0), (-3690.323301385571, 4126.352988285057, 165.0), (-3752.264180260015, 4224.460725743738, 156.0), (-3814.2050591344578, 4322.5684632024195, 147.0), (-3876.1459380089027, 4420.676200661102, 139.0), (-3938.086816883345, 4518.783938119781, 130.0), (-4000.027695757788, 4616.891675578463, 121.0), (-4061.9685746322316, 4714.999413037142, 112.0), (-4123.909453506675, 4813.107150495825, 104.0), (-4185.85033238112, 4911.2148879545075, 95.0), (-4247.791211255563, 5009.322625413187, 86.0), (-4309.732090130007, 5107.430362871868, 78.0), (-4371.672969004452, 5205.538100330549, 69.0), (-4433.613847878894, 5303.64583778923, 60.0), (-4495.554726753338, 5401.753575247911, 52.0), (-4557.495605627782, 5499.861312706593, 43.0), (-4619.436484502226, 5597.969050165275, 34.0), (-4681.37736337667, 5696.076787623955, 26.0), (-4743.318242251114, 5794.184525082637, 17.0), (-4805.259121125557, 5892.292262541318, 8.0), (-4867.2, 5990.4, 0.0)]))
    hydrology.node(10).rivers.append(asLineString([(-2307.9138171051, -795.4702929502099, 521.0), (-2410.2467404855824, -468.92229961709603, 494.0), (-2589.329360063512, -134.91271958570707, 467.0), (-2831.3001169335002, 182.3566023926544, 441.0), (-3122.2974521901565, 458.683821566686, 415.0), (-3448.459806928091, 669.867093185086, 388.0), (-3795.9256222419135, 791.7045724965521, 361.0), (-4151.091001563071, 801.3935929862843, 334.0), (-4506.991045729097, 712.1830927193814, 307.0), (-4863.721772535623, 575.6647549940313, 279.0), (-5221.7447113457165, 444.94402147870636, 250.0), (-5582.80921091596, 354.15517661449695, 222.0), (-5950.597772563143, 308.95479954551905, 194.0), (-6328.904136334305, 312.6569890399968, 168.0), (-6713.000142958087, 368.9172563554919, 143.0), (-7079.890150895739, 482.1229621134004, 119.0), (-7404.218439853203, 656.7560186697469, 95.0), (-7660.629289536417, 897.2983383805561, 72.0), (-7823.766979651329, 1208.2318336018511, 49.0), (-7868.275789903877, 1594.0384166896583, 25.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(10).rivers.append(asLineString([(-4397.990228017062, -1094.8102298855326, 347.0), (-4498.859306167035, -1035.8036007119815, 338.0), (-4599.424615996055, -976.3965093063811, 330.0), (-4699.8886697180915, -916.8559304900978, 321.0), (-4800.453979547113, -857.4488390844976, 312.0), (-4901.323057697085, -798.4422099109462, 304.0), (-5002.597160274992, -739.969530380126, 295.0), (-5103.972518959884, -681.6303382599888, 286.0), (-5205.044109323825, -622.8906839078032, 278.0), (-5305.406906938878, -563.2166176808375, 269.0), (-5404.655887377105, -502.0741899363601, 260.0), (-5502.386026210574, -438.9294510316396, 252.0), (-5598.192299011345, -373.2484513239445, 243.0), (-5691.669681351484, -304.4972411705435, 234.0), (-5782.413148803056, -232.14187092870486, 225.0), (-5870.017676938121, -155.64839095569727, 217.0), (-5954.078241328751, -74.48285160878905, 208.0), (-6034.189817546999, 11.888696754751237, 199.0), (-6109.9473811649395, 104.00020377765492, 191.0), (-6180.945907754631, 202.38561910265378, 182.0), (-6246.780372888135, 307.5788923724788, 173.0)]))
    hydrology.node(11).rivers.append(asLineString([(-3139.1312650010436, 2351.1519658273155, 347.0), (-3377.346916022355, 2419.7214594739, 330.0), (-3612.6330740757653, 2459.4139283808986, 313.0), (-3845.3598789669527, 2473.5831382165215, 295.0), (-4075.897470501597, 2465.582854648985, 278.0), (-4304.615988485376, 2438.766843346499, 261.0), (-4531.885572723969, 2396.488869977277, 244.0), (-4758.076363023056, 2342.1027002095316, 226.0), (-4983.558499188315, 2278.9620997114757, 209.0), (-5208.702121025426, 2210.420834151323, 192.0), (-5433.877368340067, 2139.832669197284, 175.0), (-5659.454380937916, 2070.5513705175726, 157.0), (-5885.8032986246535, 2005.930703780401, 140.0), (-6113.294261205958, 1949.324434653983, 122.0), (-6342.297408487509, 1904.08632880653, 105.0), (-6573.182880274984, 1873.5701519062554, 88.0), (-6806.320816374065, 1861.1296696213717, 70.0), (-7042.081356590426, 1870.1186476200917, 53.0), (-7280.83464072975, 1903.890851570628, 35.0), (-7522.950808597715, 1965.800047141193, 17.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(12).rivers.append(asLineString([(-1496.566016258495, -2609.517313645959, 521.0), (-1902.0427869123275, -2564.1902269704024, 495.0), (-2252.993566029566, -2636.569132844976, 468.0), (-2567.0038080586255, -2790.8980761702683, 442.0), (-2861.6589674479214, -2991.421101846868, 415.0), (-3154.5444986458688, -3202.3822547753643, 388.0), (-3463.245856100884, -3388.0255798563458, 362.0), (-3804.701652299752, -3513.0509265168816, 336.0), (-4174.917221672106, -3556.909016048467, 311.0), (-4544.918186536708, -3516.652807268503, 286.0), (-4884.26610902004, -3390.366925618178, 261.0), (-5169.5825971566965, -3180.700487378298, 237.0), (-5410.488617232776, -2911.6374939582925, 212.0), (-5626.238695793359, -2613.390277025972, 187.0), (-5836.08851279935, -2316.171913960441, 162.0), (-6059.293748211665, -2050.1954821408135, 136.0), (-6315.110081991214, -1845.6740589461945, 110.0), (-6622.793194098913, -1732.820721755695, 83.0), (-7001.598764495668, -1741.8485479484234, 56.0), (-7470.782473142393, -1902.9706149034887, 28.0), (-8049.599999999999, -2246.3999999999996, 0.0)]))
    hydrology.node(13).rivers.append(asLineString([(-1636.5946626095858, -4565.827754152543, 347.0), (-1646.1730635493357, -5080.695788498145, 331.0), (-1708.2719073314167, -5462.95259193207, 316.0), (-1816.967819551049, -5728.148467464976, 300.0), (-1966.3374258034555, -5891.8337181075285, 286.0), (-2150.4573516838564, -5969.558646870382, 271.0), (-2363.4042227874747, -5976.873556764207, 256.0), (-2599.2546647095314, -5929.32875079966, 242.0), (-2852.085303045248, -5842.474531987407, 227.0), (-3115.9727633898465, -5731.861203338108, 212.0), (-3384.9936713385478, -5613.039067862422, 197.0), (-3653.224652486575, -5501.558428571012, 181.0), (-3914.7423324291485, -5412.969588474543, 165.0), (-4163.623336761489, -5362.822850583673, 148.0), (-4393.944291078819, -5366.668517909067, 130.0), (-4599.781820976363, -5440.056893461385, 111.0), (-4775.212552049337, -5598.538280251289, 91.0), (-4914.313109892967, -5857.662981289441, 70.0), (-5011.160120102473, -6232.981299586501, 48.0), (-5059.830208273076, -6740.043538153134, 25.0), (-5054.4, -7394.4, 0.0)]))
    hydrology.node(14).rivers.append(asLineString([(1363.0351974719797, -1185.286063471653, 521.0), (1336.143255634835, -1538.7001375219454, 495.0), (1353.7876112410288, -1890.806818069666, 468.0), (1399.585961103793, -2241.938110932066, 442.0), (1457.1560020363577, -2592.426021926398, 415.0), (1510.115430851954, -2942.6025568699124, 389.0), (1542.0819443638127, -3292.799721579862, 363.0), (1537.332270095805, -3643.3501728819615, 336.0), (1503.966575399799, -3994.6101010438947, 310.0), (1480.0316243673785, -4346.965278396289, 284.0), (1505.4616936351763, -4700.802289486349, 258.0), (1596.4084750411462, -5054.946165214905, 231.0), (1697.6362068321528, -5403.534668912686, 205.0), (1740.5704612390641, -5739.829753197196, 179.0), (1666.514067439315, -6058.414368534419, 153.0), (1494.5032524558872, -6364.2676362974225, 127.0), (1280.3573891727174, -6667.288106077629, 102.0), (1080.115304937674, -6977.403677606516, 76.0), (949.8158270986277, -7304.5422506155555, 51.0), (945.4977830034463, -7658.631724836225, 25.0), (1123.1999999999998, -8049.599999999999, 0.0)]))
    hydrology.node(15).rivers.append(asLineString([(3913.3561082532806, -5762.491568073917, 173.0), (3941.671617470887, -5875.354352961782, 165.0), (3969.6734937288447, -5988.2709265598605, 156.0), (3997.570825666919, -6101.205429728012, 147.0), (4025.572701924877, -6214.122003326092, 139.0), (4053.888211142483, -6326.984788213957, 130.0), (4082.621897639623, -6439.775854821535, 121.0), (4111.460128456644, -6552.548991859043, 112.0), (4139.984726314018, -6665.375917606766, 104.0), (4167.77751393221, -6778.328350344989, 95.0), (4194.420314031686, -6891.478008353995, 86.0), (4219.494949332916, -7004.896609914071, 78.0), (4242.583242556366, -7118.655873305503, 69.0), (4263.267016422504, -7232.827516808577, 60.0), (4281.128093651797, -7347.483258703574, 52.0), (4295.7482969647135, -7462.694817270785, 43.0), (4306.70944908172, -7578.533910790492, 34.0), (4313.593372723284, -7695.072257542983, 26.0), (4315.981890609875, -7812.381575808544, 17.0), (4313.456825461959, -7930.533583867455, 8.0), (4305.599999999999, -8049.599999999999, 0.0)]))
    hydrology.node(16).rivers.append(asLineString([(2670.0365109674985, 419.2884533342088, 521.0), (2758.264111976863, 21.011373584195212, 494.0), (2926.651794441099, -311.0764901513366, 466.0), (3154.1793564607824, -595.4881549065561, 438.0), (3419.8265961364887, -850.7366377156324, 411.0), (3702.5733115687926, -1095.334955612735, 383.0), (3981.399300858273, -1347.796125632034, 356.0), (4234.90454593771, -1626.0646861888808, 329.0), (4437.4543218884455, -1941.7470034297003, 302.0), (4560.7652406984325, -2302.4851353310232, 274.0), (4579.093414254613, -2714.2145129507085, 247.0), (4529.222011255771, -3142.9813144125806, 219.0), (4512.284236539314, -3513.7794918325594, 192.0), (4631.555599532296, -3750.4087423165456, 166.0), (4928.082793724015, -3829.7399568500273, 141.0), (5337.348556229627, -3818.67314092122, 117.0), (5784.599861521532, -3792.8377629774704, 94.0), (6195.083684072124, -3827.8632914661257, 72.0), (6494.046998353801, -3999.379194834534, 48.0), (6606.73677883896, -4383.014941530043, 25.0), (6458.4, -5054.4, 0.0)]))
    hydrology.node(17).rivers.append(asLineString([(4397.765121188958, 648.2121881900085, 347.0), (4651.83421662707, 740.611340037247, 330.0), (4896.645099683202, 769.6176426668867, 313.0), (5132.486817131576, 745.0848207510403, 296.0), (5359.648415746434, 676.8665989618216, 279.0), (5578.418942302003, 574.8167019713433, 262.0), (5789.087443572522, 448.78885445171846, 245.0), (5991.942966332221, 308.636781075061, 228.0), (6187.274557355335, 164.2142065134836, 211.0), (6375.371263416098, 25.374855439099605, 194.0), (6556.522131288741, -98.02754747597797, 176.0), (6731.0162077475, -196.13927755963596, 159.0), (6899.1425395666065, -259.10661013976136, 142.0), (7061.190173520294, -277.07582054424023, 125.0), (7217.448156382796, -240.19318410096025, 107.0), (7368.205534928347, -138.60497613780808, 89.0), (7513.751355931181, 37.54252801733003, 72.0), (7654.37466616553, 298.10305303656673, 54.0), (7790.364512405627, 652.9303235920152, 36.0), (7922.009941425705, 1111.8780643557889, 18.0), (8049.599999999999, 1684.8, 0.0)]))
    hydrology.node(18).rivers.append(asLineString([(707.4833463640622, 676.893349318148, 521.0), (728.42915148112, 1116.62928459072, 495.0), (850.9006452244694, 1459.65755506244, 470.0), (1055.781178299755, 1724.2255016738395, 446.0), (1323.9541014126216, 1928.580465365449, 421.0), (1636.3027652687138, 2090.9697870777986, 397.0), (1973.7105205736757, 2229.6408077514193, 372.0), (2317.0608532457245, 2362.8406987258013, 346.0), (2652.1008268062815, 2502.7161053142636, 320.0), (2982.293269900961, 2639.192276353892, 293.0), (3315.1335274120825, 2757.13635895912, 267.0), (3657.350771781155, 2843.220282797769, 240.0), (4009.286300147806, 2899.1631461433994, 213.0), (4368.107623518778, 2934.1601630049645, 187.0), (4730.959275855543, 2957.460671721965, 160.0), (5094.985791119567, 2978.314010633898, 134.0), (5457.331703272321, 3005.9695180802623, 107.0), (5815.141546275272, 3049.6765324005564, 81.0), (6165.55985408989, 3118.684391934279, 54.0), (6505.731160677642, 3222.242435020927, 27.0), (6832.799999999999, 3369.6, 0.0)]))
    hydrology.node(19).rivers.append(asLineString([(1686.5153682825016, 4331.337680237611, 173.0), (1616.2295998683765, 4423.650796225731, 165.0), (1545.9438314542517, 4515.9639122138515, 156.0), (1475.6580630401263, 4608.277028201971, 147.0), (1405.3722946260018, 4700.590144190091, 139.0), (1335.0865262118766, 4792.90326017821, 130.0), (1264.8007577977514, 4885.21637616633, 121.0), (1194.514989383626, 4977.529492154448, 112.0), (1124.2292209695013, 5069.842608142568, 104.0), (1053.9434525553763, 5162.155724130688, 95.0), (983.657684141251, 5254.468840118806, 86.0), (913.371915727126, 5346.781956106925, 78.0), (843.0861473130008, 5439.095072095044, 69.0), (772.8003788988756, 5531.408188083163, 60.0), (702.5146104847504, 5623.721304071281, 52.0), (632.2288420706254, 5716.034420059401, 43.0), (561.9430736565002, 5808.34753604752, 34.0), (491.6573052423751, 5900.66065203564, 26.0), (421.3715368282501, 5992.9737680237595, 17.0), (351.0857684141249, 6085.286884011879, 8.0), (280.79999999999995, 6177.599999999999, 0.0)]))
    hydrology.node(20).rivers.append(asLineString([(-1017.1741446275357, 1726.7018189998541, 347.0), (-964.6526348680636, 1956.5615396404905, 329.0), (-879.6740056147408, 2188.403951736881, 311.0), (-772.8167213108993, 2422.225661544413, 291.0), (-654.6592463998693, 2658.0232753184723, 272.0), (-535.7800453249822, 2895.7933993144434, 253.0), (-426.75758252956894, 3135.532639787715, 233.0), (-338.1703224569609, 3377.2376029936713, 214.0), (-280.59672955048904, 3620.9048951877, 195.0), (-264.61526825348403, 3866.5311226251847, 177.0), (-298.54746294559914, 4113.922071012762, 160.0), (-372.42535190600745, 4361.337180263369, 143.0), (-467.40858016649497, 4606.285744118389, 127.0), (-564.5982592182844, 4846.272107405172, 111.0), (-645.0955005525993, 5078.80061495108, 96.0), (-690.0014156606621, 5301.375611583464, 80.0), (-680.417116033696, 5511.501442129682, 65.0), (-597.4437131629237, 5706.682451417089, 49.0), (-422.1823185395683, 5884.42298427304, 33.0), (-135.73404365485237, 6042.227385524893, 17.0), (280.79999999999995, 6177.599999999999, 0.0)]))
    hydrology.node(21).rivers.append(asLineString([(-3628.3824225111284, 4028.2452508263764, 173.0), (-3690.323301385571, 4126.352988285057, 165.0), (-3752.264180260015, 4224.460725743738, 156.0), (-3814.2050591344578, 4322.5684632024195, 147.0), (-3876.1459380089027, 4420.676200661102, 139.0), (-3938.086816883345, 4518.783938119781, 130.0), (-4000.027695757788, 4616.891675578463, 121.0), (-4061.9685746322316, 4714.999413037142, 112.0), (-4123.909453506675, 4813.107150495825, 104.0), (-4185.85033238112, 4911.2148879545075, 95.0), (-4247.791211255563, 5009.322625413187, 86.0), (-4309.732090130007, 5107.430362871868, 78.0), (-4371.672969004452, 5205.538100330549, 69.0), (-4433.613847878894, 5303.64583778923, 60.0), (-4495.554726753338, 5401.753575247911, 52.0), (-4557.495605627782, 5499.861312706593, 43.0), (-4619.436484502226, 5597.969050165275, 34.0), (-4681.37736337667, 5696.076787623955, 26.0), (-4743.318242251114, 5794.184525082637, 17.0), (-4805.259121125557, 5892.292262541318, 8.0), (-4867.2, 5990.4, 0.0)]))
    hydrology.node(22).rivers.append(asLineString([(-2307.9138171051, -795.4702929502099, 521.0), (-2410.2467404855824, -468.92229961709603, 494.0), (-2589.329360063512, -134.91271958570707, 467.0), (-2831.3001169335002, 182.3566023926544, 441.0), (-3122.2974521901565, 458.683821566686, 415.0), (-3448.459806928091, 669.867093185086, 388.0), (-3795.9256222419135, 791.7045724965521, 361.0), (-4151.091001563071, 801.3935929862843, 334.0), (-4506.991045729097, 712.1830927193814, 307.0), (-4863.721772535623, 575.6647549940313, 279.0), (-5221.7447113457165, 444.94402147870636, 250.0), (-5582.80921091596, 354.15517661449695, 222.0), (-5950.597772563143, 308.95479954551905, 194.0), (-6328.904136334305, 312.6569890399968, 168.0), (-6713.000142958087, 368.9172563554919, 143.0), (-7079.890150895739, 482.1229621134004, 119.0), (-7404.218439853203, 656.7560186697469, 95.0), (-7660.629289536417, 897.2983383805561, 72.0), (-7823.766979651329, 1208.2318336018511, 49.0), (-7868.275789903877, 1594.0384166896583, 25.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(23).rivers.append(asLineString([(-4397.990228017062, -1094.8102298855326, 347.0), (-4498.859306167035, -1035.8036007119815, 338.0), (-4599.424615996055, -976.3965093063811, 330.0), (-4699.8886697180915, -916.8559304900978, 321.0), (-4800.453979547113, -857.4488390844976, 312.0), (-4901.323057697085, -798.4422099109462, 304.0), (-5002.597160274992, -739.969530380126, 295.0), (-5103.972518959884, -681.6303382599888, 286.0), (-5205.044109323825, -622.8906839078032, 278.0), (-5305.406906938878, -563.2166176808375, 269.0), (-5404.655887377105, -502.0741899363601, 260.0), (-5502.386026210574, -438.9294510316396, 252.0), (-5598.192299011345, -373.2484513239445, 243.0), (-5691.669681351484, -304.4972411705435, 234.0), (-5782.413148803056, -232.14187092870486, 225.0), (-5870.017676938121, -155.64839095569727, 217.0), (-5954.078241328751, -74.48285160878905, 208.0), (-6034.189817546999, 11.888696754751237, 199.0), (-6109.9473811649395, 104.00020377765492, 191.0), (-6180.945907754631, 202.38561910265378, 182.0), (-6246.780372888135, 307.5788923724788, 173.0)]))
    hydrology.node(24).rivers.append(asLineString([(-3139.1312650010436, 2351.1519658273155, 347.0), (-3377.346916022355, 2419.7214594739, 330.0), (-3612.6330740757653, 2459.4139283808986, 313.0), (-3845.3598789669527, 2473.5831382165215, 295.0), (-4075.897470501597, 2465.582854648985, 278.0), (-4304.615988485376, 2438.766843346499, 261.0), (-4531.885572723969, 2396.488869977277, 244.0), (-4758.076363023056, 2342.1027002095316, 226.0), (-4983.558499188315, 2278.9620997114757, 209.0), (-5208.702121025426, 2210.420834151323, 192.0), (-5433.877368340067, 2139.832669197284, 175.0), (-5659.454380937916, 2070.5513705175726, 157.0), (-5885.8032986246535, 2005.930703780401, 140.0), (-6113.294261205958, 1949.324434653983, 122.0), (-6342.297408487509, 1904.08632880653, 105.0), (-6573.182880274984, 1873.5701519062554, 88.0), (-6806.320816374065, 1861.1296696213717, 70.0), (-7042.081356590426, 1870.1186476200917, 53.0), (-7280.83464072975, 1903.890851570628, 35.0), (-7522.950808597715, 1965.800047141193, 17.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(25).rivers.append(asLineString([(-1496.566016258495, -2609.517313645959, 521.0), (-1902.0427869123275, -2564.1902269704024, 495.0), (-2252.993566029566, -2636.569132844976, 468.0), (-2567.0038080586255, -2790.8980761702683, 442.0), (-2861.6589674479214, -2991.421101846868, 415.0), (-3154.5444986458688, -3202.3822547753643, 388.0), (-3463.245856100884, -3388.0255798563458, 362.0), (-3804.701652299752, -3513.0509265168816, 336.0), (-4174.917221672106, -3556.909016048467, 311.0), (-4544.918186536708, -3516.652807268503, 286.0), (-4884.26610902004, -3390.366925618178, 261.0), (-5169.5825971566965, -3180.700487378298, 237.0), (-5410.488617232776, -2911.6374939582925, 212.0), (-5626.238695793359, -2613.390277025972, 187.0), (-5836.08851279935, -2316.171913960441, 162.0), (-6059.293748211665, -2050.1954821408135, 136.0), (-6315.110081991214, -1845.6740589461945, 110.0), (-6622.793194098913, -1732.820721755695, 83.0), (-7001.598764495668, -1741.8485479484234, 56.0), (-7470.782473142393, -1902.9706149034887, 28.0), (-8049.599999999999, -2246.3999999999996, 0.0)]))
    hydrology.node(26).rivers.append(asLineString([(-1636.5946626095858, -4565.827754152543, 347.0), (-1646.1730635493357, -5080.695788498145, 331.0), (-1708.2719073314167, -5462.95259193207, 316.0), (-1816.967819551049, -5728.148467464976, 300.0), (-1966.3374258034555, -5891.8337181075285, 286.0), (-2150.4573516838564, -5969.558646870382, 271.0), (-2363.4042227874747, -5976.873556764207, 256.0), (-2599.2546647095314, -5929.32875079966, 242.0), (-2852.085303045248, -5842.474531987407, 227.0), (-3115.9727633898465, -5731.861203338108, 212.0), (-3384.9936713385478, -5613.039067862422, 197.0), (-3653.224652486575, -5501.558428571012, 181.0), (-3914.7423324291485, -5412.969588474543, 165.0), (-4163.623336761489, -5362.822850583673, 148.0), (-4393.944291078819, -5366.668517909067, 130.0), (-4599.781820976363, -5440.056893461385, 111.0), (-4775.212552049337, -5598.538280251289, 91.0), (-4914.313109892967, -5857.662981289441, 70.0), (-5011.160120102473, -6232.981299586501, 48.0), (-5059.830208273076, -6740.043538153134, 25.0), (-5054.4, -7394.4, 0.0)]))
    hydrology.node(27).rivers.append(asLineString([(1363.0351974719797, -1185.286063471653, 521.0), (1336.143255634835, -1538.7001375219454, 495.0), (1353.7876112410288, -1890.806818069666, 468.0), (1399.585961103793, -2241.938110932066, 442.0), (1457.1560020363577, -2592.426021926398, 415.0), (1510.115430851954, -2942.6025568699124, 389.0), (1542.0819443638127, -3292.799721579862, 363.0), (1537.332270095805, -3643.3501728819615, 336.0), (1503.966575399799, -3994.6101010438947, 310.0), (1480.0316243673785, -4346.965278396289, 284.0), (1505.4616936351763, -4700.802289486349, 258.0), (1596.4084750411462, -5054.946165214905, 231.0), (1697.6362068321528, -5403.534668912686, 205.0), (1740.5704612390641, -5739.829753197196, 179.0), (1666.514067439315, -6058.414368534419, 153.0), (1494.5032524558872, -6364.2676362974225, 127.0), (1280.3573891727174, -6667.288106077629, 102.0), (1080.115304937674, -6977.403677606516, 76.0), (949.8158270986277, -7304.5422506155555, 51.0), (945.4977830034463, -7658.631724836225, 25.0), (1123.1999999999998, -8049.599999999999, 0.0)]))
    hydrology.node(28).rivers.append(asLineString([(2670.0365109674985, 419.2884533342088, 521.0), (2758.264111976863, 21.011373584195212, 494.0), (2926.651794441099, -311.0764901513366, 466.0), (3154.1793564607824, -595.4881549065561, 438.0), (3419.8265961364887, -850.7366377156324, 411.0), (3702.5733115687926, -1095.334955612735, 383.0), (3981.399300858273, -1347.796125632034, 356.0), (4234.90454593771, -1626.0646861888808, 329.0), (4437.4543218884455, -1941.7470034297003, 302.0), (4560.7652406984325, -2302.4851353310232, 274.0), (4579.093414254613, -2714.2145129507085, 247.0), (4529.222011255771, -3142.9813144125806, 219.0), (4512.284236539314, -3513.7794918325594, 192.0), (4631.555599532296, -3750.4087423165456, 166.0), (4928.082793724015, -3829.7399568500273, 141.0), (5337.348556229627, -3818.67314092122, 117.0), (5784.599861521532, -3792.8377629774704, 94.0), (6195.083684072124, -3827.8632914661257, 72.0), (6494.046998353801, -3999.379194834534, 48.0), (6606.73677883896, -4383.014941530043, 25.0), (6458.4, -5054.4, 0.0)]))
    hydrology.node(29).rivers.append(asLineString([(4397.765121188958, 648.2121881900085, 347.0), (4651.83421662707, 740.611340037247, 330.0), (4896.645099683202, 769.6176426668867, 313.0), (5132.486817131576, 745.0848207510403, 296.0), (5359.648415746434, 676.8665989618216, 279.0), (5578.418942302003, 574.8167019713433, 262.0), (5789.087443572522, 448.78885445171846, 245.0), (5991.942966332221, 308.636781075061, 228.0), (6187.274557355335, 164.2142065134836, 211.0), (6375.371263416098, 25.374855439099605, 194.0), (6556.522131288741, -98.02754747597797, 176.0), (6731.0162077475, -196.13927755963596, 159.0), (6899.1425395666065, -259.10661013976136, 142.0), (7061.190173520294, -277.07582054424023, 125.0), (7217.448156382796, -240.19318410096025, 107.0), (7368.205534928347, -138.60497613780808, 89.0), (7513.751355931181, 37.54252801733003, 72.0), (7654.37466616553, 298.10305303656673, 54.0), (7790.364512405627, 652.9303235920152, 36.0), (7922.009941425705, 1111.8780643557889, 18.0), (8049.599999999999, 1684.8, 0.0)]))
    hydrology.node(30).rivers.append(asLineString([(707.4833463640622, 676.893349318148, 521.0), (728.42915148112, 1116.62928459072, 495.0), (850.9006452244694, 1459.65755506244, 470.0), (1055.781178299755, 1724.2255016738395, 446.0), (1323.9541014126216, 1928.580465365449, 421.0), (1636.3027652687138, 2090.9697870777986, 397.0), (1973.7105205736757, 2229.6408077514193, 372.0), (2317.0608532457245, 2362.8406987258013, 346.0), (2652.1008268062815, 2502.7161053142636, 320.0), (2982.293269900961, 2639.192276353892, 293.0), (3315.1335274120825, 2757.13635895912, 267.0), (3657.350771781155, 2843.220282797769, 240.0), (4009.286300147806, 2899.1631461433994, 213.0), (4368.107623518778, 2934.1601630049645, 187.0), (4730.959275855543, 2957.460671721965, 160.0), (5094.985791119567, 2978.314010633898, 134.0), (5457.331703272321, 3005.9695180802623, 107.0), (5815.141546275272, 3049.6765324005564, 81.0), (6165.55985408989, 3118.684391934279, 54.0), (6505.731160677642, 3222.242435020927, 27.0), (6832.799999999999, 3369.6, 0.0)]))
    hydrology.node(31).rivers.append(asLineString([(-1017.1741446275357, 1726.7018189998541, 347.0), (-964.6526348680636, 1956.5615396404905, 329.0), (-879.6740056147408, 2188.403951736881, 311.0), (-772.8167213108993, 2422.225661544413, 291.0), (-654.6592463998693, 2658.0232753184723, 272.0), (-535.7800453249822, 2895.7933993144434, 253.0), (-426.75758252956894, 3135.532639787715, 233.0), (-338.1703224569609, 3377.2376029936713, 214.0), (-280.59672955048904, 3620.9048951877, 195.0), (-264.61526825348403, 3866.5311226251847, 177.0), (-298.54746294559914, 4113.922071012762, 160.0), (-372.42535190600745, 4361.337180263369, 143.0), (-467.40858016649497, 4606.285744118389, 127.0), (-564.5982592182844, 4846.272107405172, 111.0), (-645.0955005525993, 5078.80061495108, 96.0), (-690.0014156606621, 5301.375611583464, 80.0), (-680.417116033696, 5511.501442129682, 65.0), (-597.4437131629237, 5706.682451417089, 49.0), (-422.1823185395683, 5884.42298427304, 33.0), (-135.73404365485237, 6042.227385524893, 17.0), (280.79999999999995, 6177.599999999999, 0.0)]))
    hydrology.node(32).rivers.append(asLineString([(-2307.9138171051, -795.4702929502099, 521.0), (-2410.2467404855824, -468.92229961709603, 494.0), (-2589.329360063512, -134.91271958570707, 467.0), (-2831.3001169335002, 182.3566023926544, 441.0), (-3122.2974521901565, 458.683821566686, 415.0), (-3448.459806928091, 669.867093185086, 388.0), (-3795.9256222419135, 791.7045724965521, 361.0), (-4151.091001563071, 801.3935929862843, 334.0), (-4506.991045729097, 712.1830927193814, 307.0), (-4863.721772535623, 575.6647549940313, 279.0), (-5221.7447113457165, 444.94402147870636, 250.0), (-5582.80921091596, 354.15517661449695, 222.0), (-5950.597772563143, 308.95479954551905, 194.0), (-6328.904136334305, 312.6569890399968, 168.0), (-6713.000142958087, 368.9172563554919, 143.0), (-7079.890150895739, 482.1229621134004, 119.0), (-7404.218439853203, 656.7560186697469, 95.0), (-7660.629289536417, 897.2983383805561, 72.0), (-7823.766979651329, 1208.2318336018511, 49.0), (-7868.275789903877, 1594.0384166896583, 25.0), (-7768.799999999999, 2059.2, 0.0)]))
    hydrology.node(33).rivers.append(asLineString([(-1496.566016258495, -2609.517313645959, 521.0), (-1902.0427869123275, -2564.1902269704024, 495.0), (-2252.993566029566, -2636.569132844976, 468.0), (-2567.0038080586255, -2790.8980761702683, 442.0), (-2861.6589674479214, -2991.421101846868, 415.0), (-3154.5444986458688, -3202.3822547753643, 388.0), (-3463.245856100884, -3388.0255798563458, 362.0), (-3804.701652299752, -3513.0509265168816, 336.0), (-4174.917221672106, -3556.909016048467, 311.0), (-4544.918186536708, -3516.652807268503, 286.0), (-4884.26610902004, -3390.366925618178, 261.0), (-5169.5825971566965, -3180.700487378298, 237.0), (-5410.488617232776, -2911.6374939582925, 212.0), (-5626.238695793359, -2613.390277025972, 187.0), (-5836.08851279935, -2316.171913960441, 162.0), (-6059.293748211665, -2050.1954821408135, 136.0), (-6315.110081991214, -1845.6740589461945, 110.0), (-6622.793194098913, -1732.820721755695, 83.0), (-7001.598764495668, -1741.8485479484234, 56.0), (-7470.782473142393, -1902.9706149034887, 28.0), (-8049.599999999999, -2246.3999999999996, 0.0)]))
    hydrology.node(34).rivers.append(asLineString([(1363.0351974719797, -1185.286063471653, 521.0), (1336.143255634835, -1538.7001375219454, 495.0), (1353.7876112410288, -1890.806818069666, 468.0), (1399.585961103793, -2241.938110932066, 442.0), (1457.1560020363577, -2592.426021926398, 415.0), (1510.115430851954, -2942.6025568699124, 389.0), (1542.0819443638127, -3292.799721579862, 363.0), (1537.332270095805, -3643.3501728819615, 336.0), (1503.966575399799, -3994.6101010438947, 310.0), (1480.0316243673785, -4346.965278396289, 284.0), (1505.4616936351763, -4700.802289486349, 258.0), (1596.4084750411462, -5054.946165214905, 231.0), (1697.6362068321528, -5403.534668912686, 205.0), (1740.5704612390641, -5739.829753197196, 179.0), (1666.514067439315, -6058.414368534419, 153.0), (1494.5032524558872, -6364.2676362974225, 127.0), (1280.3573891727174, -6667.288106077629, 102.0), (1080.115304937674, -6977.403677606516, 76.0), (949.8158270986277, -7304.5422506155555, 51.0), (945.4977830034463, -7658.631724836225, 25.0), (1123.1999999999998, -8049.599999999999, 0.0)]))
    hydrology.node(35).rivers.append(asLineString([(2670.0365109674985, 419.2884533342088, 521.0), (2758.264111976863, 21.011373584195212, 494.0), (2926.651794441099, -311.0764901513366, 466.0), (3154.1793564607824, -595.4881549065561, 438.0), (3419.8265961364887, -850.7366377156324, 411.0), (3702.5733115687926, -1095.334955612735, 383.0), (3981.399300858273, -1347.796125632034, 356.0), (4234.90454593771, -1626.0646861888808, 329.0), (4437.4543218884455, -1941.7470034297003, 302.0), (4560.7652406984325, -2302.4851353310232, 274.0), (4579.093414254613, -2714.2145129507085, 247.0), (4529.222011255771, -3142.9813144125806, 219.0), (4512.284236539314, -3513.7794918325594, 192.0), (4631.555599532296, -3750.4087423165456, 166.0), (4928.082793724015, -3829.7399568500273, 141.0), (5337.348556229627, -3818.67314092122, 117.0), (5784.599861521532, -3792.8377629774704, 94.0), (6195.083684072124, -3827.8632914661257, 72.0), (6494.046998353801, -3999.379194834534, 48.0), (6606.73677883896, -4383.014941530043, 25.0), (6458.4, -5054.4, 0.0)]))
    hydrology.node(36).rivers.append(asLineString([(707.4833463640622, 676.893349318148, 521.0), (728.42915148112, 1116.62928459072, 495.0), (850.9006452244694, 1459.65755506244, 470.0), (1055.781178299755, 1724.2255016738395, 446.0), (1323.9541014126216, 1928.580465365449, 421.0), (1636.3027652687138, 2090.9697870777986, 397.0), (1973.7105205736757, 2229.6408077514193, 372.0), (2317.0608532457245, 2362.8406987258013, 346.0), (2652.1008268062815, 2502.7161053142636, 320.0), (2982.293269900961, 2639.192276353892, 293.0), (3315.1335274120825, 2757.13635895912, 267.0), (3657.350771781155, 2843.220282797769, 240.0), (4009.286300147806, 2899.1631461433994, 213.0), (4368.107623518778, 2934.1601630049645, 187.0), (4730.959275855543, 2957.460671721965, 160.0), (5094.985791119567, 2978.314010633898, 134.0), (5457.331703272321, 3005.9695180802623, 107.0), (5815.141546275272, 3049.6765324005564, 81.0), (6165.55985408989, 3118.684391934279, 54.0), (6505.731160677642, 3222.242435020927, 27.0), (6832.799999999999, 3369.6, 0.0)]))

    return edgelength, shore, hydrology, cells