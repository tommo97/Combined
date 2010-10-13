/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 35               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2009-11-16 00:1#$:  Date of last commit

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

        if (argc < 2){
        cout << "Incorrect number of input arguments. Command expected is:" << endl;
        cout << "\t./main case_name" << endl;
        cout << "where case_name.tar.gz is a bundle in the tarballs directory" << endl;
        cout << "containing the directory structure:" << endl;
        cout << "\tcase_name/" << endl << "\tcase_name/case_name.cas" << endl;
        cout << "\tcase_name/case_name.neu" << endl << "\tcase_name/case_name.mat" << endl;
        cout << "Aborting." << endl;
        return 1;
     }

    
    SYSTEM System(0);


    //  Some default values
    globalSystem->GambitScale = 1;
    globalSystem->MaxP = 3;
    globalSystem->dtInit = 0.01;
    globalSystem->Del2 = .001;
    globalSystem->DS = .3;
    globalSystem->NeuFile = "neu_files/0012.neu";








    string input = argv[1];



    //  Make a scratch directory
    string command = "mkdir -p scratch";
    cout << command << endl;
    string output = globalGetStdoutFromCommand(command);
    //  Copy the input file into scratch and untar
    command =  "cp tarballs/" + input + ".tar.gz scratch/";
    cout << command << endl;
    output = globalGetStdoutFromCommand(command);
    command =  "tar -xvf scratch/" + input + ".tar.gz -C scratch/";
    cout << command << endl;
    output = globalGetStdoutFromCommand(command);
    //  Copy the files from that tar.gz into the new directory - in
    globalSystem->WorkingDir = "scratch/" + input + "/";


 

    //command = "cp tarballs/" + input + ".t* scratch/; tar -xvf scratch/" + input + ".tar.gz -C scratch/";
    cout << command << endl;
    output = globalGetStdoutFromCommand(command);

    globalIO->read_input("scratch/" + input + "/" + input + ".cas");

    globalIO->PrepOutputDir();
    globalIO->WriteBinary();
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "globalSystem->MaxP set to " << globalSystem->MaxP << "; dtInit " << globalSystem->dtInit << endl;
#endif

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

/**************************************************************/
Vect3 globalDirectVel(Vect3 diff, Vect3 omega) {

    REAL mult, nrm;
    nrm = sqrt(globalSystem->Del2 + diff.Dot(diff));
    mult = -1 / (four_pi * nrm * nrm * nrm);
    return mult * diff.Cross(omega);
}
