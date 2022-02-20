import struct
import os

import DataModel

def writeDataModel(path: str, edgeLength: float, shore: DataModel.ShoreModel, hydrology: DataModel.HydrologyNetwork=None, cells: DataModel.TerrainHoneycomb=None, Ts: DataModel.Terrain=None):
    with open(path, 'wb') as file:
        file.write(struct.pack('!H', 2)) # version number. Increment every time a breaking change is made

        ## Contents ##
        tableOfContents = {
            'shore': 0,
            'hydrology': 0,
            'honeycomb': 0,
            'primitives': 0
        }
        file.write(struct.pack('!Q', 0))
        file.write(struct.pack('!Q', 0))
        file.write(struct.pack('!Q', 0))
        file.write(struct.pack('!Q', 0))

        tableOfContents['shore'] += struct.calcsize('!H') + struct.calcsize('!Q') + struct.calcsize('!Q') + struct.calcsize('!Q') + struct.calcsize('!Q')

        ## Parameters ##

        file.write(struct.pack('!B', 0)) # coordinate type
        file.write(struct.pack('!f', shore.resolution)) # resolution
        file.write(struct.pack('!f', edgeLength))

        sectionSize = struct.calcsize('!H') + struct.calcsize('!B') + struct.calcsize('!f')*2
        print(f"\tBasic parameters: {sectionSize} bytes")
        tableOfContents['shore'] += sectionSize

        ## Shore ##

        tableOfContents['hydrology'] = tableOfContents['shore']

        sectionSize = 0

        file.write(struct.pack('!Q', shore.rasterShape[0]))
        file.write(struct.pack('!Q', shore.rasterShape[1]))
        for d0 in range(shore.rasterShape[0]):
            for d1 in range(shore.rasterShape[1]):
                file.write(struct.pack('!B', shore.imgray[d0][d1]))

        sectionSize = struct.calcsize('!Q')*2 + struct.calcsize('!B')*shore.rasterShape[0] * shore.rasterShape[1]
        print(f"\tshore.imgray: {sectionSize} bytes")
        tableOfContents['hydrology'] += sectionSize

        sectionSize = 0

        # write the shore contour
        file.write(struct.pack('!Q', len(shore.contour)))
        for point in shore.contour:
            file.write(struct.pack('!Q', point[0]))
            file.write(struct.pack('!Q', point[1]))

        sectionSize = struct.calcsize('!Q') + struct.calcsize('!Q') * 2 * len(shore.contour)
        print(f"\tShore contour: {sectionSize} bytes")
        tableOfContents['hydrology'] += sectionSize

        ## Hydrology data structure ##

        if hydrology is None:
            # Abort writing the file; there's nothing left to write

            ## Write the table of contents ##
            file.seek(struct.calcsize('!H'))
            file.write(struct.pack('!Q', tableOfContents['shore']))
            file.write(struct.pack('!Q', 0))
            file.write(struct.pack('!Q', 0))
            file.write(struct.pack('!Q', 0))

            file.close()

            return

        tableOfContents['honeycomb'] = tableOfContents['hydrology']

        sectionSize = 0

        # write all hydrology primitives
        file.write(struct.pack('!Q', len(hydrology)))
        sectionSize += struct.calcsize('!Q')
        for node in hydrology.allNodes():
            file.write(struct.pack('!I', node.id))
            file.write(struct.pack('!f', node.x()))
            file.write(struct.pack('!f', node.y()))
            file.write(struct.pack('!f', node.elevation))
            file.write(struct.pack('!I', node.parent.id if node.parent is not None else node.id))
            file.write(struct.pack('!I', node.contourIndex if node.parent is None else 0))
            file.write(struct.pack('!B', len(node.rivers)))
            for river in node.rivers:
                file.write(struct.pack('!H', len(river.coords)))
                sectionSize += struct.calcsize('!H')
                for point in river.coords:
                    file.write(struct.pack('!f', point[0]))
                    file.write(struct.pack('!f', point[1]))
                    file.write(struct.pack('!f', point[2]))
                    sectionSize += struct.calcsize('!f')*3
            if cells is not None:
                file.write(struct.pack('!f', node.localWatershed))
                file.write(struct.pack('!f', node.inheritedWatershed))
                file.write(struct.pack('!f', node.flow))
            else:
                file.write(struct.pack('!f', 0.0))
                file.write(struct.pack('!f', 0.0))
                file.write(struct.pack('!f', 0.0))
            sectionSize += struct.calcsize('!I')*2 + struct.calcsize('!H') + struct.calcsize('!f')*6 + struct.calcsize('!B')

        print(f"\tHydrology nodes: {sectionSize} bytes")
        tableOfContents['honeycomb'] += sectionSize

        ## TerrainHoneycomb data structure ##

        if cells is None:
            # Abort writing the file; there's nothing left to write

            ## Write the table of contents ##
            file.seek(struct.calcsize('!H'))
            file.write(struct.pack('!Q', tableOfContents['shore']))
            file.write(struct.pack('!Q', tableOfContents['hydrology']))
            file.write(struct.pack('!Q', tableOfContents['honeycomb']))
            file.write(struct.pack('!Q', 0))

            file.close()

            return

        tableOfContents['primitives'] = tableOfContents['honeycomb']

        sectionSize = 0

        # write point_region
        file.write(struct.pack('!Q', len(cells.point_region)))
        sectionSize += struct.calcsize('!Q')
        for idx in cells.point_region:
            if idx != -1:
                file.write(struct.pack('!Q', idx))
                sectionSize += struct.calcsize('!Q')

        print(f"\tPoint_Region: {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        # write regions array

        sectionSize = 0

        file.write(struct.pack('!Q', len(cells.regions)))
        sectionSize += struct.calcsize('!Q')
        for region in cells.regions:
            file.write(struct.pack('!B', len(region)))
            sectionSize += struct.calcsize('!B')
            for point in region:
                if point != -1:
                    file.write(struct.pack('!Q', point))
                else:
                    file.write(struct.pack('!Q', 0xffffffffffffffff))
                sectionSize += struct.calcsize('!Q')

        print(f"\tRegions: {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        # write vertices
        sectionSize = 0
        file.write(struct.pack('!Q', len(cells.vertices)))
        sectionSize += struct.calcsize('!Q')
        for vertex in cells.vertices:
            file.write(struct.pack('!f', vertex[0]))
            file.write(struct.pack('!f', vertex[1]))
            sectionSize += struct.calcsize('!f')*2

        print(f"\tVertices: {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        # qs
        sectionSize = 0
        file.write(struct.pack('!Q', len(cells.qs)))
        sectionSize += struct.calcsize('!Q')
        for q in cells.qs:
            if q is None:
                file.write(struct.pack('!B', 0x00))
                sectionSize += struct.calcsize('!B')
            else:
                file.write(struct.pack('!B', 0xff))
                file.write(struct.pack('!f', q.position[0]))
                file.write(struct.pack('!f', q.position[1]))
                file.write(struct.pack('!B', len(q.nodes)))
                for node in q.nodes:
                    file.write(struct.pack('!Q', node))
                    sectionSize += struct.calcsize('!Q')
                file.write(struct.pack('!Q', q.vorIndex))
                file.write(struct.pack('!f', q.elevation))
                sectionSize += struct.calcsize('!Q') + struct.calcsize('!f')*3 + struct.calcsize('!B')*2

        print(f"\tQs: {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        # cellsRidges
        sectionSize = 0
        file.write(struct.pack('!Q', len(cells.cellsRidges)))
        sectionSize += struct.calcsize('!Q')
        for cellID in cells.cellsRidges:
            file.write(struct.pack('!Q', cellID))
            file.write(struct.pack('!B', len(cells.cellsRidges[cellID])))
            for ridge in cells.cellsRidges[cellID]:
                file.write(struct.pack('!B', len(ridge)))
                file.write(struct.pack('!Q', ridge[0].vorIndex))
                sectionSize += struct.calcsize('!Q') + struct.calcsize('!B')
                if len(ridge) > 1:
                    file.write(struct.pack('!Q', ridge[1].vorIndex))
                    sectionSize += struct.calcsize('!Q')
            sectionSize += struct.calcsize('!Q') + struct.calcsize('!B')

        print(f"\tcellsRidges: {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        # cellsDownstreamRidges
        sectionSize = 0
        file.write(struct.pack('!Q', len(cells.cellsDownstreamRidges)))
        sectionSize += struct.calcsize('!Q')
        for cellID in cells.cellsDownstreamRidges:
            file.write(struct.pack('!Q', cellID))
            sectionSize += struct.calcsize('!Q')
            if cells.cellsDownstreamRidges[cellID] is not None:
                file.write(struct.pack('!B', 0x00))
                file.write(struct.pack('!f', cells.cellsDownstreamRidges[cellID][0][0]))
                file.write(struct.pack('!f', cells.cellsDownstreamRidges[cellID][0][1]))
                file.write(struct.pack('!f', cells.cellsDownstreamRidges[cellID][1][0]))
                file.write(struct.pack('!f', cells.cellsDownstreamRidges[cellID][1][1]))
                sectionSize += struct.calcsize('!B') + struct.calcsize('!f')*4
            else:
                file.write(struct.pack('!B', 0xff))
                sectionSize += struct.calcsize('!B')

        print(f"\tcellsDownstreamRidges {sectionSize} bytes")
        tableOfContents['primitives'] += sectionSize

        ## Terrain primitives ##

        if Ts is None:
            # Abort writing the file; there's nothing left to write

            ## Write the table of contents ##
            file.seek(struct.calcsize('!H'))
            file.write(struct.pack('!Q', tableOfContents['shore']))
            file.write(struct.pack('!Q', tableOfContents['hydrology']))
            file.write(struct.pack('!Q', tableOfContents['honeycomb']))
            file.write(struct.pack('!Q', 0))

            file.close()

            return

        sectionSize = 0

        file.write(struct.pack('!Q', len(Ts)))
        sectionSize += struct.calcsize('!Q')
        for t in Ts.allTs():
            file.write(struct.pack('!f', t.position[0]))
            file.write(struct.pack('!f', t.position[1]))
            file.write(struct.pack('!I', t.cell))
            file.write(struct.pack('!f', t.elevation))
            sectionSize += struct.calcsize('!I') + struct.calcsize('!f')*3

        print(f"\tTerrain primitives: {sectionSize} bytes")

        ## Write the table of contents ##
        file.seek(struct.calcsize('!H'))
        file.write(struct.pack('!Q', tableOfContents['shore']))
        file.write(struct.pack('!Q', tableOfContents['hydrology']))
        file.write(struct.pack('!Q', tableOfContents['honeycomb']))
        file.write(struct.pack('!Q', tableOfContents['primitives']))

        file.close()

def writeToTerrainModule(pipe, shore, edgeLength, hydrology, cells, Ts):
    pipe.write(struct.pack('!f', 0))
    pipe.write(struct.pack('!f', 0))
    pipe.write(struct.pack('!f', shore.realShape[0]))
    pipe.write(struct.pack('!f', shore.realShape[1]))

    pipe.write(struct.pack('!f', edgeLength))
    pipe.write(struct.pack('!f', shore.resolution)) # resolution

    # write the shore contour
    pipe.write(struct.pack('!Q', len(shore.contour)))
    pipe.write(struct.pack('!I', shore.imgray.shape[0]))
    pipe.write(struct.pack('!I', shore.imgray.shape[1]))
    for point in shore.contour:
        pipe.write(struct.pack('!Q', point[0]))
        pipe.write(struct.pack('!Q', point[1]))

    # write all hydrology primitives
    pipe.write(struct.pack('!Q', len(hydrology)))
    for node in hydrology.allNodes():
        pipe.write(struct.pack('!I', node.id))
        pipe.write(struct.pack('!f', node.x()))
        pipe.write(struct.pack('!f', node.y()))
        pipe.write(struct.pack('!f', node.elevation))
        pipe.write(struct.pack('!I', node.parent.id if node.parent is not None else node.id))
        pipe.write(struct.pack('!I', node.contourIndex if node.parent is None else 0))
        pipe.write(struct.pack('!B', len(node.rivers)))
        for river in node.rivers:
            pipe.write(struct.pack('!H', len(river.coords)))
            for point in river.coords:
                pipe.write(struct.pack('!f', point[0]))
                pipe.write(struct.pack('!f', point[1]))
                pipe.write(struct.pack('!f', point[2]))
        pipe.write(struct.pack('!f', node.localWatershed))
        pipe.write(struct.pack('!f', node.inheritedWatershed))
        pipe.write(struct.pack('!f', node.flow if node.parent is not None else 0))

    # qs
    pipe.write(struct.pack('!Q', len(cells.qs)))
    for q in cells.qs:
        if q is None:
            pipe.write(struct.pack('!B', 0x00))
        else:
            pipe.write(struct.pack('!B', 0xff))
            pipe.write(struct.pack('!f', q.position[0]))
            pipe.write(struct.pack('!f', q.position[1]))
            pipe.write(struct.pack('!f', q.elevation))
            pipe.write(struct.pack('!Q', q.vorIndex))
            pipe.write(struct.pack('!B', len(q.nodes)))
            for node in q.nodes:
                pipe.write(struct.pack('!Q', node))

    # cellsRidges
    pipe.write(struct.pack('!Q', len(cells.cellsRidges)))
    for cellID in cells.cellsRidges:
        pipe.write(struct.pack('!Q', cellID))
        pipe.write(struct.pack('!B', len(cells.cellsRidges[cellID])))
        for ridge in cells.cellsRidges[cellID]:
            pipe.write(struct.pack('!B', len(ridge)))
            pipe.write(struct.pack('!Q', ridge[0].vorIndex))
            if len(ridge) > 1:
                pipe.write(struct.pack('!Q', ridge[1].vorIndex))

    pipe.write(struct.pack('!Q', len(Ts)))
    for t in Ts.allTs():
        pipe.write(struct.pack('!f', t.position[0]))
        pipe.write(struct.pack('!f', t.position[1]))
        pipe.write(struct.pack('!I', t.cell))

def readDataModel(path):
    with open(path, 'rb') as file:
        versionNumber = struct.unpack('!H', file.read(struct.calcsize('!H')))[0]
        if versionNumber != 2:
            raise ValueError('This file was created with a previous version of the hydrology script. It is not compatible')

        # Just ignore the TOC for now
        file.seek(struct.calcsize('!Q') * 4, os.SEEK_CUR)

        coordinateType = struct.unpack('!B', file.read(struct.calcsize('!B')))[0]
        rasterResolution = struct.unpack('!f', file.read(struct.calcsize('!f')))[0]
        edgeLength = struct.unpack('!f', file.read(struct.calcsize('!f')))[0]

        shore = DataModel.ShoreModel(rasterResolution, binaryFile=file)
        hydrology = DataModel.HydrologyNetwork(binaryFile=file)
        cells = DataModel.TerrainHoneycomb(resolution=rasterResolution, edgeLength=edgeLength, shore=shore, hydrology=hydrology, binaryFile=file)
        terrain = DataModel.Terrain(binaryFile=file)

        file.close()

        return (rasterResolution, edgeLength, shore, hydrology, cells, terrain)