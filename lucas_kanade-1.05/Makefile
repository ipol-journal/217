
#
# Makefile
# 
# Copyright (C) 2016, Tristan Dagobert, CMLA, École Normale Supérieure de Cachan.
# 
# This software is a computer program.[describe
# functionalities and technical features of your software].
# 
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# 
#


# option of Make
ifeq ($(with_openmp), yes)
	FLAG_OPENMP = -DWITH_OPENMP -fopenmp
else 	
	FLAG_OPENMP = 
endif

ifeq ($(with_bicubic), yes)
	FLAG_BICUBIC = -DWITH_BICUBIC
else 	
	FLAG_BICUBIC = 
endif


# names of executables
EXE	= lucas_kanade.exe 

# search for C sources files
SRC_C := $(shell ls *.c)

# object files
OBJ_C   = $(SRC_C:.c=.o)

# compilation options

ifeq ($(DEBOGUE), yes)
	OPTS_C  = -O0 -ggdb -Wall -Wextra -Wundef -Wmissing-declarations \
	-Wmissing-prototypes -Wunused -Wformat \
	$(FLAG_BICUBIC) -fsanitize=address $(FLAG_OPENMP)

    # For debugging mode
	LIBS = -fsanitize=address
else
	OPTS_C  = -O3 -Wall -Wextra -Wundef -Wmissing-declarations\
	-Wmissing-prototypes -Wunused -Wformat $(FLAG_BICUBIC) $(FLAG_OPENMP)
endif

# GSL
DIR_GSL = -L/usr/local/lib
LIB_GSL = -lgsl -lgslcblas
INC_GSL =

# images
DIR_IMAGE =
LIB_IMAGE = -lpng -ltiff -ljpeg 
INC_IMAGE = 

# libraries
DIR_LIBS = $(DIR_GSL)
INC_LIBS = 
NOM_LIBS = -lm -lz $(LIB_IMAGE) $(LIB_GSL) $(FLAG_OPENMP) -lfftw3 $(LIBS)

default : 
	@echo "Type : make (exe [with_openmp=yes] [with_bicubic=yes] [DEBOGUE=yes] | clean);"
exe	: $(EXE)

# C source code compilation 
%.o	: %.c
	$(CC) $(WITH_OPENMP) $(OPTS_C) $(INC_DIRS) -c $< -o $@

# linkage of object code
$(EXE)	: $(OBJ_C)
	$(CC) $^ $(DIR_LIBS) $(NOM_LIBS) -o $@

# cleaning
clean	:
	$(RM) $(EXE) $(OBJ_C) *~

