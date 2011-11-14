
# $Rev:: 319           $:  Revision of last commit
# $Author:: tom        $:  Author of last commit
# $Date:: 2009-04-28 2#$:  Date of last commit


# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



CC = g++
OMP_FLAG = -fopenmp
OBJDIR = obj

vpath %.cpp src
vpath %.hpp include
vpath %.o obj


platform = Linux

ifeq ($(OSTYPE), darwin9.0)
 platform = MacOS
endif

ifeq ($(OSTYPE), darwin10.0)
 platform = MacOS
endif


ifneq ($(CC),sunCC)
#	Set up some c++ compiler flags for gcc
ifeq ($(platform), Linux)
	OPT_FLAGS = -g -O3 -march=native -funroll-loops -ftree-vectorize -fprefetch-loop-arrays -funroll-all-loops  -fomit-frame-pointer -ffast-math
else
	OPT_FLAGS = -g -O3 -funroll-loops -ftree-vectorize -fprefetch-loop-arrays -funroll-all-loops  -fomit-frame-pointer -ffast-math
endif
	DEBUG_FLAGS = -O0 -g -Wall -Wextra -Wunused -Wuninitialized -Winit-self -Wshadow -W  -Wshadow -fno-common -g -ansi -pedantic -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums
	OMP_DEBUG_FLAG = $(OMP_FLAG)
	
else
#	Set up some c++ compiler flags for Sun
	OPT_FLAGS= -DFILE=__FILE -g0 -fast -fns -xipo -xtarget=generic -xarch=sse2 -xvector=simd -xalias_level=simple -xrestrict
	DEBUG_FLAGS = -g
	ifeq ($(OMP_FLAG),-fopenmp)
		OMP_FLAG = -xopenmp=parallel
		OMP_DEBUG_FLAG = -xopenmp=noopt
	endif
endif

SVNDEF := -D'SVN_REV=$(shell svnversion -n .)'


#	Common flags
CC_PNG_FLAGS = -DNO_FREETYPE
LD_PNG_FLAGS = -DNO_FREETYPE  -lpng -lpngwriter -lz
CC_FLAGS = -m64 -g -I/usr/local/include -I include
LD_COMMON_FLAGS = -m64 $(OMP_FLAG) -lmatio -lncurses -lm -lz
LD_DEBUG_COMMON_FLAGS = -m64 $(OMP_DEBUG_FLAG) -lncurses -lm

CC_DEBUG_FLAGS =  -c $(SVNDEF) $(CC_PNG_FLAGS) $(CC_FLAGS) $(OMP_DEBUG_FLAG) $(DEBUG_FLAGS)
CC_COMMON_FLAGS = -c $(SVNDEF) $(CC_PNG_FLAGS) $(CC_FLAGS) $(OMP_FLAG) $(OPT_FLAGS)


ifeq ($(platform), MacOS)
    LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -lcblas -framework Accelerate -lgfortran -lgsl
else
    ifneq ($(CC),sunCC)
	ifeq ($(OMP_FLAG),-fopenmp)
	        LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -L/usr/local/atlas/lib/ -llapack -lptf77blas -lptcblas -latlas -lgsl
	    else
		LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -L/usr/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas -lgsl
	endif
    else
	LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS)  $(OPT_FLAGS) -llapack -lgsl -lgslcblas
    endif
endif







CFLAGS = $(CC_COMMON_FLAGS)




CFLAGS = $(CC_COMMON_FLAGS)

SRC = types.cpp utils.cpp pgesv.cpp panel.cpp body.cpp system.cpp io.cpp tree.cpp node.cpp branch.cpp cell.cpp time_integrator.cpp
HEADERS = $(SRC:.cpp=.hpp)
SOURCES = $(SRC) main.cpp 


OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main

all: $(EXECUTABLE) $(SOURCES) $(HEADERS)

$(EXECUTABLE):$(addprefix $(OBJDIR)/, $(OBJECTS))
	$(CC) $^ $(LD_FLAGS) -o $@
	
sun: clean
	$(MAKE) all "CC= CC"

$(OBJDIR)/%.o : %.cpp
	$(CC) $(GSLLIBS) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f main *.o *~ $(OBJDIR)/*

valgrind: clean debug
	valgrind -v  --leak-check=full --show-reachable=yes ./main > valgrind.out

valgrind-openmp: clean debug-openmp
	valgrind -v  --leak-check=full --show-reachable=yes ./main > valgrind.out

rebuild: clean all

debug:
	$(MAKE) all "CFLAGS = $(CC_DEBUG_FLAGS)"

debug-openmp: clean 
	$(MAKE) all "CFLAGS= -c -Wall -g -m64 -O0 -fopenmp"

