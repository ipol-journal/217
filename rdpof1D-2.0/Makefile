
# executables
EXE2     = rdpof1D.exe

# looking after C++ source files
SRC_CPP := $(shell ls *.cpp)

# looking after C sources files
SRC_C := $(shell ls *.c)

# object files
OBJ_CPP	= $(SRC_CPP:.cpp=.o)
OBJ_C   = $(SRC_C:.c=.o)

# compilation options
OPTS     = -Wall -Wextra -O3 #-ggdb

CC 	= gcc
#CC	= gcc
#OPTS_C  = $(OPTS) -std=c99 -pedantic -DNDEBUG -D_GNU_SOURCE
OPTS_C  = $(OPTS) -std=c99 -DNDEBUG -D_GNU_SOURCE -DABSOLUTE_TEST=1

CPP	= g++
#CPP	= g++
OPTS_CPP= $(OPTS) -std=c++98 -fopenmp
#OPTS_CPP= $(OPTS) -fopenmp

# librairies
INC_LIBS = 
NOM_LIBS =  -lm -lpng -ljpeg -ltiff -fopenmp

#exempt = test_disp.o
exempt = 
# building by default 
default	 : $(EXE2)

# linking the object code
$(EXE2) : $(OBJ_CPP) $(filter-out $(exempt),$(OBJ_C))
	$(CPP) $^  -o $@  $(NOM_LIBS) $(DIR_LIBS)

# 1A) C++ code compilation
%.o	: %.cpp
	$(CPP) $(OPTS_CPP) $(OPTS_COUV) $(INC_LIBS) -c $< -o $@ 

# 1B) C code compilation
%.o	: %.c
	$(CC) $(OPTS_C) -c $< -o $@

clean   :
	$(RM) $(OBJ_CPP) $(OBJ_C) $(EXE1) $(EXE2)


