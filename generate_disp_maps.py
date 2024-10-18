#!/usr/bin/env python3

import os
import sys
import iio
#-----------------------------------------------------------------------
if __name__ == '__main__':
  
    """
    arg1 : disparity map estimated in TIF (input)
    arg2 : disparity map estimated in PNG (output)
    arg3 : disparity map truth in TIF (input)
    """
    in_disp_map = sys.argv[1]
    out_disp_png = sys.argv[2]

    if (len(sys.argv) == 3):
        in_disp_truth = sys.argv[1]
    elif (len(sys.argv) == 4):
        in_disp_truth = sys.argv[3]

    im = iio.read(in_disp_truth).squeeze()
    
    vmin = im.min()
    vmax = im.max()

    if not os.path.exists("max_and_min.txt"):
        f = open("max_and_min.txt","a")
        f.write("\nvmin= "+str(vmin)+"\n")
        f.write("vmax= "+str(vmax)+"\n")
        f.close()

    command="qeasy " + str(vmin) + " " + str(vmax) + " " + in_disp_map + " " + out_disp_png
    os.system(command)

    
