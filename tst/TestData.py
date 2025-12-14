import math
import struct
from os.path import join, basename, splitext
import sys

import numpy as np

coresuffix = ".zdft"
updisuffix = ".zdfi"

print(f"The system's byte order is: {sys.byteorder}")
byteorder = '>' if (sys.byteorder == 'big') else '<'

if len(sys.argv) < 1 or len(sys.argv) > 3:
    print(f"Argument list Error. Must be {sys.argv[0]} {{/desired/zdf/dir encodingInteger, /path/to/output.zdfo}}", file=sys.stderr)
    sys.exit(1) 

if len(sys.argv) == 2:
    enc = int(splitext(basename(sys.argv[1]))[0])
    M = ((enc & 0xFFFF0000) >> 16).bit_count() * (((enc & 0x0000F00000000000) >> 44) - ((enc & 0x00000F0000000000) >> 40)) 
    derivFmt = byteorder+f'{M}f'
    recordSize = struct.calcsize(derivFmt) 
    with open(sys.argv[1], 'rb') as f:
        while True:
            derivx = f.read(recordSize)
            if not derivx:
                break
            print(struct.unpack(derivFmt, derivx))
    sys.exit(0) 


zdftfile = join(sys.argv[1], sys.argv[2]+coresuffix)  # Generate an initialization .zdft file
updifile = join(sys.argv[1], sys.argv[2]+updisuffix)  # Generate a file with observations to filter
N = int(sys.argv[2]) & 0x0000FFFF                # Length of filter

U = 256      # Number of new observations to generate
stdY = 0.2  # Observable random noise standard deviation 

T = 6*np.pi  # Time horizon for initialization series
rate = np.log(np.sqrt(2))/(2*np.pi)  # Attenuation for generated oscillatory series

C = np.exp(-1)/(1+rate*rate)
dt = T/N
tvec = np.linspace(-(N-1)*dt, T, num=U+N, endpoint=True)

y = -np.exp(-rate*tvec)*(np.cos(tvec) + rate*np.sin(tvec)) / (1+rate*rate) + C
yr = y + np.random.normal(loc=0, scale=stdY, size=U+N)

with open(zdftfile, 'wb') as f:
    yFmt = f'{N}f'
    f.write(struct.pack(byteorder+yFmt, *yr[0:N].flatten('F')))

with open(updifile, 'wb') as f:
    yFmt = f'{U}f'
    f.write(struct.pack(byteorder+yFmt, *yr[N:].flatten('F')))


