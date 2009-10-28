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


#include "types.hpp"
#include "includes.hpp"
#include "array.hpp"

long unsigned int globalNum_NODES = 0, globalNum_BRANCHES = 0, globalNum_FVMCELLS = 0;
long unsigned int globalFactorial[24] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600, 6227020800, 87178291200,
    1307674368000, 2432902008176640000};

long int ticks() {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (now.tv_sec - ticks_global.start.tv_sec)*1000 + (now.tv_usec - ticks_global.start.tv_usec) / 1000;
}

long int globalTimer = ticks();
//  A small number of global variables, mainly for bookkeeping
SYSTEM *globalSystem;
IO *globalIO;
TIME_STEPPER *globalTimeStepper;
OCTREE *globalOctree;
/**************************************************************/
Array <REAL> globalLinspace(REAL start, REAL end, int n)
{
    Array <REAL> output(n);

    REAL d = (end - start)/(n-1);
    for (int i = 0; i < n; ++i)
        output[i] = start + i*d;

    return output;

}
/**************************************************************/
Array <Vect3> globalLinspace(Vect3 start, Vect3 end, int n)
{
    Array <Vect3> output(n);

    Vect3 d = (end - start)/(n-1);
    for (int i = 0; i < n; ++i)
        output[i] = start + i*d;

    return output;

}
/**************************************************************/
string globalGetStdoutFromCommand(string cmd)
{
  // setup
  string data;
  FILE *stream;
  char buffer[160];

  // do it
  stream = popen(cmd.c_str(), "r");
  while ( fgets(buffer, 160, stream) != NULL )
    data.append(buffer);
  pclose(stream);

  // exit
  return data;
}
