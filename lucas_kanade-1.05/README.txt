
/*
* README.txt
* 
* Copyright (C) 2016, Tristan Dagobert, CMLA, École Normale Supérieure de Cachan.
* 
* This software is a computer program.[describe
* functionalities and technical features of your software].
* 
* This software is governed by the CeCILL-C license under French law and
* abiding by the rules of distribution of free software.  You can  use, 
* modify and/ or redistribute the software under the terms of the CeCILL-C
* license as circulated by CEA, CNRS and INRIA at the following URL
* "http://www.cecill.info". 
* 
* As a counterpart to the access to the source code and  rights to copy,
* modify and redistribute granted by the license, users are provided only
* with a limited warranty  and the software's author,  the holder of the
* economic rights,  and the successive licensors  have only  limited
* liability. 
* 
* In this respect, the user's attention is drawn to the risks associated
* with loading,  using,  modifying and/or developing or reproducing the
* software by the user in light of its specific status of free software,
* that may mean  that it is complicated to manipulate,  and  that  also
* therefore means  that it is reserved for developers  and  experienced
* professionals having in-depth computer knowledge. Users are therefore
* encouraged to load and test the software's suitability as regards their
* requirements in conditions enabling the security of their systems and/or 
* data to be ensured and,  more generally, to use and operate it in the 
* same conditions as regards security. 
* 
* The fact that you are presently reading this means that you have had
* knowledge of the CeCILL-C license and that you accept its terms.
* 
*/
    
NAME
		lucas_kanade.exe

SYNOPSIS 

		./lucas_kanade.exe --version;
		./lucas_kanade.exe --help;
		./lucas_kanade.exe FIC_IM_L FIC_IM_R FIC_DISP NB_PROCS VX VY NB_SCALES SIGMA_DV \
		NB_ITER THRESHOLD TYPE_TEST;

OPTIONS
		--version            display the program version
		--help               display this help

DESCRIPTION
		This software is a multi-scale 1D version of the Lucas &
		Kanade's optical flow estimation algorithm. From two
		stereo-rectified images it produces a sparse disparity
		map. This map contains numerical values at pixels where a
		good estimation had been made, and NaN values at pixels where
		estimation wàs considered by the algorithm as erroneous.

		Create the disparity map FIC_DISP from left image FIC_IM_L
		and right image FIC_IM_R considered respectively as the
		reference and the shifted images. A neighborhood of VX by VY
		pixels is used around the central pixel to estimate its
		displacement, while a presmoothing by a Gaussian kernel with
		a SIGMA_DV standard deviation, is applied to both images at
		each scale. The inner gradient descent stops after NB_ITER
		iterations maximum. The left right test consistency is
		absolute TYPE_TEST=1 or relative TYPE_TEST=0 with a threshold
		THRESHOLD. In the absolute case, the threshold value
		corresponds to a distance in pixel. In the relative case, its
		value corresponds to a ratio. If the THRESHOLD is negative,
		the test consistency is not done : the dispary map will be dense. 

EXAMPLE
		./lucas_kanade.exe input_0.png input_1.png disp.tif 4 5 5 5 0.1 5 0.5 0;
		
REFERENCES
		Reference to the IPOL article :

		This algorithm was designed for the publication entitled :
		«Comparison of optical flow methods under  stereomatching with short baselines».

COPYRIGHT PATENT AND LICENSE INFORMATION
		THIS program is written by Tristan Dagobert. All the files
        are distributed under the terms of CeCILL-C License, except
	    for the following contributions :

		- bicubic_interpolation.c, bicubic_interpolation.h :
		  Javier Sánchez Pérez, Nelson Monzón López (initial coding)
		  These files are distributed under the terms of the simplified BSD
		  License.
		  
		- iio.c, iio.h by :
		  Enric Meinhardt-Llopis (initial coding)
		  Gabriele Facciolo (many patches, python bindings, cmake files)
		  Carlo de Franchis (patches)
		  Juan Cardelino (cmake configuration, bug reports)
		  These files are distributed under the terms of the GNU Affero General
		  Public License.


AUTHORS

		Tristan Dagobert	<tristan.dagobert@cmla.ens-cachan.fr>

VERSION
		2019/01/05      1.05
		2016/10/10      1.0

		Location of future releases and updates : http://www.ipol.im

CHANGES
		In the program since it was first published in IPOL
		No changes.

CONFIGURATION 
	    Libraries required to compile and use the program

		This software requires :
		- GNU GSL library >= 2.0
		- PNG, TIF, JPEG libraries

		Compilation instructions

		Before compiling, please check the Makefile in order to program the
		right paths for external libraries. 

		You can specify a mode with or without the Openmp.

		To compile, type :
		make (exe [with_openmp=yes] [with_bicubic=yes] [DEBOGUE=yes] | clean);

		This will produce the executable : lucas_kanade.exe
		
		Note that the option «with_bicubic» will to speed up the computation
		thanks to a numerical scheme different from the GSL.
		However results will vary from the original version.
				
NOTES
		The input and output images are all types supported by the iio library. 




