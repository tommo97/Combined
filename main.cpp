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


#include "includes.hpp"
#include "tree.hpp"
#include "node.hpp"
#include "system.hpp"


using namespace std;


int main(int argc, char *argv[]) {
    system("clear");
    SYSTEM System(0);


    string dir = "./neu_files/";

    //  Some default values
    globalSystem->GambitScale = 1;
    globalSystem->MaxP = 3;
    globalSystem->dtInit = 0.01;
    globalSystem->Del2 = 0.25;
    globalSystem->NeuFile = dir + "0012.neu";


    switch (argc) {

        case 2:
        {
            globalIO->read_input(argv[1]);
            break;
        }

        case 3:
        {
            globalSystem->NeuFile = dir + argv[1];
            globalSystem->GambitScale = atof(argv[2]);
            break;
        }
        case 4:
        {
            globalSystem->NeuFile = dir + argv[1];
            globalSystem->GambitScale = atof(argv[2]);
            globalSystem->MaxP = atoi(argv[3]);
            break;
        }
                case 5:
        {
            globalSystem->NeuFile = dir + argv[1];
            globalSystem->GambitScale = atof(argv[2]);
            globalSystem->MaxP = atoi(argv[3]);
            globalSystem->dtInit = atof(argv[4]);
            break;
        }
        default:
        {
            break;
        }

    }

    globalIO->PrepOutputDir();
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "globalSystem->MaxP set to " << globalSystem->MaxP << "; dtInit " << globalSystem->dtInit << endl;
#endif
    globalSystem->uinf = 0.0;
    globalSystem->vinf = 0.0;
    globalSystem->winf = 0.0;
    globalSystem->Initialise();

    globalSystem->TimeStep();
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000 << " seconds" << endl;
#endif
}
/**************************************************************/
void globalDirectVel(Vect3 diff, Vect3 omega, Vect3 & vel) {

    REAL mult, nrm;
    nrm = sqrt(globalSystem->Del2 + diff.Dot(diff));
    mult = -1 / (four_pi * nrm * nrm * nrm);
    vel += mult * diff.Cross(omega);
}
