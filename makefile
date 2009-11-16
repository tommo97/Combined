
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

#   -Werror
#      Consider warnings to be errors, so that compilation stops. This
#      prevents warnings from scrolling off the top of the screen and being
#      lost. You won't be able to compile the program until it is completely
#      warning-free.  You may be tempted to remove this.  DON'T! 
# 
#   -Wall
#      This turns on a set of warnings for common programming problems. You
#      need -Wall, but it is not enough on its own (for g++). 
# 
#   -W
#      This turns on some extra warnings in g++ not included in -Wall, such 
#      as missing return values and comparisons between signed and unsigned
#      integers. 
# 
#   -Wshadow
#      This warns whenever a local variable shadows another local variable.
#      If two variables have the same name then it is a potential source of
#      confusion. 
# 
#   -fno-common
#      This option prevents global variables being simultaneously defined in
#      different object files (you get an error at link time). Such a
#      variable should be defined in one file and referred to in other files
#      with an extern declaration. 
# 
#   -O3
#      Turn on optimization. The warnings for uninitialized variables in
#      -Wall rely on the optimizer to analyze the code. If there is no
#      optimization then the warnings aren't generated. 
# 
#   -g
#      It always makes sense to put debugging symbols in the executable so
#      that you can debug it using gdb. The only effect of debugging symbols
#      is to increase the size of the file, and you can use the "strip"
#      command to remove them later if necessary.
# 
# -ansi -pedantic
#      Use ISO C++, and reject any non-ANSI extensions. These flags help in
#      writing portable programs that will compile on other systems.
# 
#   -Wconversion
#      The main use of this option is to warn about conversions from signed
#      to unsigned integers. For example, unsigned int x = -1. If you need
#      to perform such a conversion you can use an explicit cast. 
# 
#   -Wpointer-arith -Wcast-qual -Wcast-align
#      These options warn if you try to do pointer arithmetic for types
#      which don't have a size, such as void, if you remove a const cast
#      from a pointer, or if you cast a pointer to a type which has a
#      different size, causing an invalid alignment. 
# 
#   -Wwrite-strings
#      This option gives string constants a const qualifier so that it will
#      be a compile-time error to attempt to overwrite them. 
# 
#   -fshort-enums
#      This option makes the type of enum as short as possible. Normally
#      this makes an enum different from an int. Consequently any attempts
#      to assign a pointer-to-int to a pointer-to-enum will generate a
#      cast-alignment warning. 


CC = g++
OMP_FLAG = -fopenmp 
OBJDIR = obj

vpath %.cpp src
vpath %.hpp include
vpath %.o obj


ifneq ($(CC),CC)
#	Set up some c++ compiler flags for gcc
	OPT_FLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -fprefetch-loop-arrays -funroll-all-loops  -fomit-frame-pointer -ffast-math
    DEBUG_FLAGS = -O0  -Wall -W  -Wshadow -fno-common -g -ansi -pedantic -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums
	OMP_DEBUG_FLAG = $(OMP_FLAG)

else
#	Set up some c++ compiler flags for Sun
	OPT_FLAGS = -xO5 -fast -xipo=2 -xvector=simd
	DEBUG_FLAGS = -g
	ifeq ($(OMP_FLAG),-fopenmp)
		OMP_FLAG = -xopenmp=parallel
		OMP_DEBUG_FLAG = -xopenmp=noopt
	endif

endif


#	Common flags
CC_PNG_FLAGS =
#-DNO_FREETYPE
LD_PNG_FLAGS =
#-DNO_FREETYPE  -L/usr/local/lib64 -lpng -lpngwriter -lz
CC_FLAGS = -m64 -g -I/usr/local/include -I include
LD_COMMON_FLAGS = -m64 $(OMP_FLAG) -lncurses -lm
LD_DEBUG_COMMON_FLAGS = -m64 $(OMP_DEBUG_FLAG) -lncurses -lm

CC_DEBUG_FLAGS =  -c $(CC_PNG_FLAGS) $(CC_FLAGS) $(OMP_DEBUG_FLAG) $(DEBUG_FLAGS)
CC_COMMON_FLAGS = -c $(CC_PNG_FLAGS) $(CC_FLAGS) $(OMP_FLAG) $(OPT_FLAGS)


platform = Linux
  ifeq ($(HOSTNAME), Chameleon)
   LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -lgsl -llapack -lcblas -latlas -lgfortran 
  endif
  ifeq ($(HOSTNAME), masternode)
   LD_FLAGS= $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -lgsl -lcblas -latlas -lgoto_opteron-r1.26 -lgfortran
endif

ifeq ($(OSTYPE), darwin9.0)
 platform = MacOS
 LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -lcblas -latlas -lgfortran -lgsl
endif

ifeq ($(OSTYPE), darwin10.0)
 platform = MacOS
 LD_FLAGS = $(LD_COMMON_FLAGS) $(LD_PNG_FLAGS) -lcblas -latlas -lgfortran -lgsl
endif





CFLAGS = $(CC_COMMON_FLAGS)

SRC = types.cpp system.cpp io.cpp tree.cpp node.cpp branch.cpp cell.cpp panel.cpp body.cpp time_integrator.cpp
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

