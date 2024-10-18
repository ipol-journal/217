#!/usr/bin/env python3

import os
import sys
import iio

#-----------------------------------------------------------------------
if __name__ == '__main__':
  
    """
    arg0 : the command
    arg1 : im1 (input)
    arg2 : im2 (input)
    arg3 : disparity map truth (input) optional
    """

    f = open("algo_info.txt","a")

    is_ok = True
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    file3 = sys.argv[3] if (len(sys.argv) == 4) else None
    
    # check the image format
    im1 = iio.read(file1)
    im2 = iio.read(file2)
    
    if im1 is None or im2 is None:
        f.write("message=Input pairs cannot be read in a good format.\n")
        is_ok = False
        
    if file3 != None:
        try:
            im3 = iio.read(file3)
            print("ici")
        except:
            print("la")
            f.write("message=Input disparity cannot be read in a good TIFF format.\n")
            is_ok = False            
    else:
        im3 = None

    if is_ok == False:
        f.write("message=Please, check your data.\n")
        f.close()
        exit(5)
        
    # check the triplets image size
    if (im1.shape != im2.shape):
        f.write("message=Input pairs have not the same size\n")
        is_ok = False
    if im3 is not None:
        if im3.shape[0:2] != im1.shape[0:2]:
            f.write("message=Input disparity have not the size of the pair.\n")
            is_ok = False
            
    # check the number of bands of the input pairs
    if im1.shape[2] > 3 or im2.shape[2] > 3:
        f.write("message=Input pairs have not 3 channels\n")
        is_ok = False
        
    if is_ok == False:
        f = open("algo_info.txt","a")
        f.write("message=Please check your data.\n")
        f.close()
        exit(5)
    
    exit(0)
    
