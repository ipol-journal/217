#! /usr/bin/python

import os
import sys
import numpy as np
import iio
import re

def load_pfm(filename):
    color = None
    width = None
    height = None
    scale = None
    endian = None
    
    file = open(filename, "r")
    header = file.readline().rstrip()
    
    if header == 'PF':
        color = True    
    elif header == 'Pf':
        color = False
    else:
        raise Exception('Not a PFM file.')

    dim_match = re.match(r'^(\d+)\s(\d+)\s$', file.readline())
    if dim_match:
        width, height = map(int, dim_match.groups())
    else:
        raise Exception('Malformed PFM header.')

    scale = float(file.readline().rstrip())
    if scale < 0: # little-endian
        endian = '<'
        scale = -scale
    else:
        endian = '>' # big-endian

    data = np.fromfile(file, endian + 'f')
    shape = (height, width, 3) if color else (height, width)
    return np.reshape(data, shape), scale


#-----------------------------------------------------------------------
if __name__ == '__main__':

    disp_map = sys.argv[1]

    image, scale = load_pfm(disp_map)

    vmin = image.min()
    vmax = image.min()                                                  
    for value in image:
        for v in value:
            if not np.isinf(v) and v > vmax:
                vmax = v

    image = 255 - (((image - vmin)/ (vmax-vmin)) * 255)
    image = np.flip(image,0)

    im = image.astype(np.uint8)
    iio.write(disp_map, im)

