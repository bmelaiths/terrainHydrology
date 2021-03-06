#!/usr/bin/env python

import os.path
import subprocess
import sys
import struct

buildRiversExe = 'buildRivers' # path to the executable

if not os.path.exists(buildRiversExe):
    print('The executable does not exist. Run "make" to build it.')
    exit()

i = int(2).to_bytes(1, sys.byteorder, signed=False)
j = int(2).to_bytes(1, sys.byteorder, signed=False)

input = i + j

proc = subprocess.run(
    ['./' + buildRiversExe],
    input=input,
    capture_output=True
)

result = proc.stdout

print(f'Result: {struct.unpack("B", result[0:1])[0]}')