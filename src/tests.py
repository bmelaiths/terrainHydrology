#! /bin/python

import unittest

from PIL import Image
from PIL import ImageDraw
from matplotlib.pyplot import draw
import os.path

import HydrologyFunctions
from DataModel import *
import Math

from TerrainPrimitiveFunctions import computePrimitiveElevation
from RiverInterpolationFunctions import computeRivers

import testcodegenerator
from testcodegenerator import RasterDataMock

class HydrologyFunctionTests(unittest.TestCase):
    def setUp(self):
        #create shore
        r = 100
        image = Image.new('L', (4*r,4*r))

        drawer = ImageDraw.Draw(image)

        #Keep in mind that these are _image_ coordinates
        drawer.polygon([(150,134),(100,200),(150,286),(250,286),(300,200),(250,134)], 255)
        #This should work out to (-500,660), (-1000,0), (-500,-860), (500,-860), (1000,0), (500,660)

        image.save('imageFile.png')

        shore = ShoreModel(10, 'imageFile.png')

        #create candidate nodes
        candidateNodes = [ ]

        candidateNodes.append(HydroPrimitive(0, None,  4.0, 1, None))
        candidateNodes.append(HydroPrimitive(1, None,  6.0, 2, None))
        candidateNodes.append(HydroPrimitive(2, None, 14.0, 3, None))
        candidateNodes.append(HydroPrimitive(3, None,  8.0, 3, None))
        candidateNodes.append(HydroPrimitive(4, None, 24.0, 1, None))
        candidateNodes.append(HydroPrimitive(5, None, 23.0, 4, None))

        #create hydrology
        hydrology = HydrologyNetwork()

        #create real mouth nodes
        self.node0 = hydrology.addNode(shore[215], 0, 0, 215) # This should be (130,-860)
        hydrology.addNode(shore[225], 0, 0, 225) # This should be (230,-860)
        hydrology.addNode(shore[235], 0, 0, 235) # This should be (330,-860)

        hydrology.addNode((130,-760), 10, 0, parent=self.node0)
        
        #create the parameters object
        edgelength = 100
        sigma = 0.75
        eta = 0.5
        zeta = 14
        self.params = HydrologyFunctions.HydrologyParameters(shore, hydrology, None, None, None, None, edgelength, sigma, eta, zeta, None, None, candidateNodes)

    def test_select_node(self):
        selectedNode = HydrologyFunctions.selectNode(self.params.candidates, self.params.zeta)

        self.assertEqual(self.params.zeta, 14.0)
        self.assertEqual(selectedNode.id, 3)
    
    def test_is_acceptable_position_not_on_land(self):
        acceptable0 = HydrologyFunctions.isAcceptablePosition((-100,-900), self.params)
        acceptable1 = HydrologyFunctions.isAcceptablePosition((-100,-700), self.params)
        
        self.assertFalse(acceptable0)
        self.assertTrue(acceptable1)
    
    def test_is_acceptable_position_too_close_to_seeee(self):
        acceptable = HydrologyFunctions.isAcceptablePosition((-100,-830), self.params)

        self.assertFalse(acceptable)
    
    def test_is_acceptable_position_too_close_to_nodes_or_edges(self):
        acceptable0 = HydrologyFunctions.isAcceptablePosition((80,-800), self.params)
        acceptable1 = HydrologyFunctions.isAcceptablePosition((100,-600), self.params)
        
        self.assertFalse(acceptable0)
        self.assertTrue(acceptable1)
    
    def test_coast_normal(self):
        angle = HydrologyFunctions.coastNormal(self.node0, self.params)
        
        self.assertAlmostEqual(angle, math.pi * 0.5, places=3)

    # def test_pick_new_node_position(self):
    #     pass
    
    def tearDown(self) -> None:
        os.remove('imageFile.png')

class HoneycombTests(unittest.TestCase):
    def setUp(self) -> None:
        resolution = 2

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
        self.hydrology = HydrologyNetwork()

        for i in range(0, 12):
            idx = int(42 * (i + 0.5))

            self.hydrology.addNode(shore[idx], 0, 1, contourIndex=idx)
        
        candidateNodes = self.hydrology.allMouthNodes()
        
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
        self.params = HydrologyFunctions.HydrologyParameters(
            shore, self.hydrology, Pa, Pc, maxTries, riverAngleDev, edgelength,
            sigma, eta, zeta, riverSlope, slopeRate, candidateNodes
        )
        
        self.hydrology.addNode((-132.0, 90.0), 0, 2, contourIndex=21) # ID: 0
        self.hydrology.addNode((-196.0, 6.0), 0, 1, contourIndex=63) # ID: 1
        self.hydrology.addNode((-154.0, -78.0), 0, 2, contourIndex=105) # ID: 2
        self.hydrology.addNode((-106.0, -162.0), 0, 2, contourIndex=147) # ID: 3
        self.hydrology.addNode((-26.0, -172.0), 0, 2, contourIndex=189) # ID: 4
        self.hydrology.addNode((58.0, -172.0), 0, 1, contourIndex=231) # ID: 5
        self.hydrology.addNode((124.0, -130.0), 0, 1, contourIndex=273) # ID: 6
        self.hydrology.addNode((174.0, -46.0), 0, 2, contourIndex=315) # ID: 7
        self.hydrology.addNode((172.0, 38.0), 0, 1, contourIndex=357) # ID: 8
        self.hydrology.addNode((108.0, 122.0), 0, 1, contourIndex=399) # ID: 9
        self.hydrology.addNode((26.0, 132.0), 0, 1, contourIndex=441) # ID: 10
        self.hydrology.addNode((-58.0, 132.0), 0, 1, contourIndex=483) # ID: 11
        self.hydrology.addNode((-107.56437944261481, 46.37775283212138), 1.9607843137254901, 2, parent=self.hydrology.node(0)) # ID: 12
        self.hydrology.addNode((-83.14370047559429, 100.63306149618042), 1.9607843137254901, 1, parent=self.hydrology.node(0)) # ID: 13
        self.hydrology.addNode((-146.86127194616188, -3.2404223524112297), 1.9607843137254901, 1, parent=self.hydrology.node(1)) # ID: 14
        self.hydrology.addNode((-107.16302325545817, -60.49863978339832), 1.9607843137254901, 2, parent=self.hydrology.node(2)) # ID: 15
        self.hydrology.addNode((-94.41042051990941, -113.36172651630545), 1.9607843137254901, 1, parent=self.hydrology.node(3)) # ID: 16
        self.hydrology.addNode((-60.439785486129125, -141.40226096266704), 1.9607843137254901, 1, parent=self.hydrology.node(3)) # ID: 17
        self.hydrology.addNode((12.605208355834506, -140.22520042860057), 1.9607843137254901, 1, parent=self.hydrology.node(4)) # ID: 18
        self.hydrology.addNode((-26.65373632854306, -122.00427389453432), 1.9607843137254901, 1, parent=self.hydrology.node(4)) # ID: 19
        self.hydrology.addNode((80.92201635535358, -127.56374041163114), 1.9607843137254901, 1, parent=self.hydrology.node(5)) # ID: 20
        self.hydrology.addNode((106.40526327284911, -83.19802098733166), 1.9607843137254901, 1, parent=self.hydrology.node(6)) # ID: 21
        self.hydrology.addNode((130.70103093067564, -20.996014766967807), 1.9607843137254901, 1, parent=self.hydrology.node(7)) # ID: 22
        self.hydrology.addNode((164.99414127031076, 3.1822580667142617), 1.9607843137254901, 1, parent=self.hydrology.node(7)) # ID: 23
        self.hydrology.addNode((126.90899179646235, 16.394422498148785), 1.9607843137254901, 1, parent=self.hydrology.node(8)) # ID: 24
        self.hydrology.addNode((74.36064195081241, 85.00819563689063), 1.9607843137254901, 1, parent=self.hydrology.node(9)) # ID: 25
        self.hydrology.addNode((16.250216208987762, 82.95979490225886), 1.9607843137254901, 1, parent=self.hydrology.node(10)) # ID: 26
        self.hydrology.addNode((-40.32118354618386, 85.22971617798027), 1.9607843137254901, 1, parent=self.hydrology.node(11)) # ID: 27
        self.hydrology.addNode((-65.39635263629756, 19.5103466442939), 3.9215686274509802, 1, parent=self.hydrology.node(12)) # ID: 28
        self.hydrology.addNode((-101.39458260246725, -3.2401230489691457), 3.9215686274509802, 1, parent=self.hydrology.node(12)) # ID: 29
        self.hydrology.addNode((-60.56633667743154, -78.6300710547887), 3.9215686274509802, 1, parent=self.hydrology.node(15)) # ID: 30
        self.hydrology.addNode((-69.18750506255246, -27.97387970772232), 3.9215686274509802, 1, parent=self.hydrology.node(15)) # ID: 31
        self.hydrology.addNode((14.452298520252938, -90.25932949728877), 3.9215686274509802, 1, parent=self.hydrology.node(18)) # ID: 32
        self.hydrology.addNode((-20.855300532352167, -72.34163056660168), 3.9215686274509802, 1, parent=self.hydrology.node(19)) # ID: 33
        self.hydrology.addNode((68.60994108873766, -79.10331529350837), 3.9215686274509802, 1, parent=self.hydrology.node(20)) # ID: 34
        self.hydrology.addNode((85.44485341294113, -37.80351450797149), 3.9215686274509802, 1, parent=self.hydrology.node(21)) # ID: 35
        self.hydrology.addNode((88.90246415194844, 6.442640725681283), 3.9215686274509802, 1, parent=self.hydrology.node(22)) # ID: 36
        self.hydrology.addNode((104.76698191106078, 61.22445036156617), 3.9215686274509802, 1, parent=self.hydrology.node(24)) # ID: 37
        self.hydrology.addNode((38.916306770475266, 49.74207728479644), 3.9215686274509802, 1, parent=self.hydrology.node(25)) # ID: 38
        self.hydrology.addNode((-11.264501327952445, 41.21125817505853), 3.9215686274509802, 1, parent=self.hydrology.node(26)) # ID: 39
        self.hydrology.addNode((-22.84914139209466, -6.752456240036139), 5.88235294117647, 1, parent=self.hydrology.node(28)) # ID: 40
        self.hydrology.addNode((29.099290915621435, -42.4527858680949), 5.88235294117647, 1, parent=self.hydrology.node(32)) # ID: 41
        self.hydrology.addNode((47.380499198782296, -5.382765931495825), 5.88235294117647, 1, parent=self.hydrology.node(35)) # ID: 42

        self.cells = TerrainHoneycomb(shore, self.hydrology, resolution, edgelength)
    
    def test_test(self) -> None:
        pass

    def tearDown(self) -> None:
        os.remove('imageFile.png')

class RiverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.edgeLength, self.shore, self.hydrology, self.cells = testcodegenerator.getPredefinedObjects0()
    
    def test_test(self) -> None:
        node = self.hydrology.node(3)
        node.rivers = [ ]

        computeRivers(node, self.hydrology, self.cells)

        # ensure that the river does not intersect any of the mountain ridges of any cells that it flows through
        allRidges = self.cells.cellRidges(3)
        allRidges += self.cells.cellRidges(14)
        allRidges += self.cells.cellRidges(27)
        allRidges += self.cells.cellRidges(34)

        self.assertEqual(1, len(node.rivers))

        river = list(node.rivers[0].coords)
        for i in range(len(river)-2):
            p0 = river[i]
            p1 = river[i+1]
            for ridge in allRidges:
                if len(ridge) < 2:
                    continue
                self.assertFalse(Math.segments_intersect_tuple(p0, p1, ridge[0].position, ridge[1].position))

    def tearDown(self) -> None:
        os.remove('imageFile.png')

class TerrainTests(unittest.TestCase):
    def setUp(self) -> None:
        self.edgeLength, self.shore, self.hydrology, self.cells = testcodegenerator.getPredefinedObjects0()

    def test_test(self) -> None:
        t = T((1519,-734), 34)

        z = computePrimitiveElevation(t, self.shore, self.hydrology, self.cells)
        
        self.assertAlmostEqual(z, 933.975, delta=10.0)

    def tearDown(self) -> None:
        os.remove('imageFile.png')