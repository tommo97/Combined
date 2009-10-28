/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev::                  $:  Revision of last commit
$Author::               $:  Author of last commit
$Date::                 $:  Date of last commit

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#ifndef INCLUDES
#define INCLUDES
#include <iostream> // std::if (WRITE_TO_SCREEN) cout etc
#include <stdio.h>
#include <math.h> // sqrt
#include <algorithm> // copy, min, max
#include <string> // copy
#include <sstream> // isstringstream
#include <memory> // auto_ptr
#include <vector> // vector
#include <iomanip> // setfill
#include <fstream> // ifstream
#include <queue> // queue
#include <utility> // pair
#include <list> // list
#include <set> // set
#include <queue> // queue
#include <sys/time.h> // for timer
#include <locale.h> // for greek symbols
#include <stdlib.h> // system
#include <unistd.h> // getcwd
#include <sys/param.h> // globalSystem->MaxPATHLEN
#include <memory.h> // memcpy
#include <math.h>
#include <time.h>
#ifndef TEST_MODE
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <ncurses.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _PNGWRITER
#include <pngwriter.h>
#include <png.h>
#endif

#endif

#endif
