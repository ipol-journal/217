#!/usr/bin/env python3


import os
import sys
import numpy as np
import iio

#-----------------------------------------------------------------------
if __name__ == '__main__':

    gt_name = sys.argv[1]
    disp_name = sys.argv[2]

    gt_tiff = iio.read(gt_name).squeeze()
    disp_tiff = iio.read(disp_name).squeeze()
    
    im_error = abs(gt_tiff-disp_tiff)
    error = 0
    for value in im_error:
        error += sum(value)

    error /= im_error.size
    f = open("error_value.txt","a")
    f.write("Error = "+str(error)+"\n")
    f.write("Calculated as abs(ground_truth - disparity_file)\n")
    f.close()

    a = im_error.flatten()
    b = np.sort(a)
    vmin, vmax = b[int(0.02*b.size)], b[int(0.98*b.size)]
    im_error =(((im_error - vmin)/ (vmax-vmin)) * 255)
    im = im_error.astype(np.uint8)
    iio.write("error_image.png", im)
    
