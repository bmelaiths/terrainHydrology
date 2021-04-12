#!/usr/bin/env python

import os.path
import subprocess
import sys
import struct
import cv2 as cv
import numpy as np

buildRiversExe = 'buildRivers' # path to the executable

if not os.path.exists(buildRiversExe):
    print('The executable does not exist. Run "make" to build it.')
    exit()

# i = int(2).to_bytes(1, sys.byteorder, signed=False)
# j = int(2).to_bytes(1, sys.byteorder, signed=False)

# input = i + j

# Create an image
r = 100
src = np.zeros((4*r, 4*r), dtype=np.uint8)

# Create a sequence of points to make a contour
vert = [None]*6
vert[0] = (3*r//2, int(1.34*r))
vert[1] = (1*r, 2*r)
vert[2] = (3*r//2, int(2.866*r))
vert[3] = (5*r//2, int(2.866*r))
vert[4] = (3*r, 2*r)
vert[5] = (5*r//2, int(1.34*r))

# Draw it in src
for i in range(6):
    cv.line(src, vert[i],  vert[(i+1)%6], ( 255 ), 3)

# Get the contours
contours, _ = cv.findContours(src, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

binary = struct.pack('l', len(contours[0]))

binary = binary + struct.pack('l', 4*r) + struct.pack('l', 4*r)

# points in a contour array are y,x
for point in contours[0]:
    point = point[0]
    binary = binary + struct.pack('f', float(point[0])) + struct.pack('f', float(point[1]))

# file = open('binaryFile', 'w+b')
# file.write(binary)
# file.close();

proc = subprocess.run(
    ['./' + buildRiversExe],
    input=binary,
    capture_output=True
)

result = proc.stdout

print(f'Result: {struct.unpack("B", result[0:1])[0]}')