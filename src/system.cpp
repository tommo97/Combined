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



#include <deque>

#include "system.hpp"
#include "types.hpp"
#ifndef USEGSL
#include "pgesv.hpp"
#endif

/**************************************************************/
SYSTEM::~SYSTEM() {

    delete globalIO;
    delete globalOctree;
    delete globalTimeStepper;

}

/**************************************************************/
SYSTEM::SYSTEM(int NT) {
    LiftingLineMode =  false;
    num_out = 0;
    SysDumpInterval = 0; //  Used to dump the result of system calls
    scaledVinf.x = scaledVinf.y = scaledVinf.z = 0;
    unscaledVinf = scaledVinf;
    MaxP = 3;
    DS = .3;
    Temp = 288.15; //  Kelvin
    Rho = 1.226; //1027; //  Kg/m3
    Mu = (sqrt(pow(Temp, 3)) * 1.458e-6) / (Temp + 110.4); //   Dynamic Viscocity kg/ms
    GambitScale = 1;
    NumThreads = 1;
#ifdef _OPENMP
    if (NT != 0) omp_set_num_threads(NT);
#pragma omp parallel
    NumThreads = omp_get_num_threads();
#else
    if (NT != 0) NumThreads = 1;
#endif
    srand((unsigned) time(NULL));



    globalSystem = this;
    globalIO = new IO();

}

/**************************************************************/
void SYSTEM::Initialise() {

    scaledVinf = unscaledVinf*GambitScale;
    Nu = GambitScale * GambitScale * Mu / Rho;
    globalTimeStepper = new TIME_STEPPER();

    NumTransVars = BODY::Bodies.size();

    globalOctree = new OCTREE();

    if (WRITE_TO_SCREEN) cout << "rho: " << Rho << "; mu: " << Mu << "; nu: " << Nu << endl;
    if (WRITE_TO_SCREEN) cout << "FVM mesh scale Factor: " << GambitScale << endl;

    int nss = globalSystem->NumSubSteps;

    globalTimeStepper->time_step();

    globalSystem->NumSubSteps = nss;



}

/**************************************************************/
void SYSTEM::TimeStep() {

    while (TIME_STEPPER::SimTime < TIME_STEPPER::MaxTime)
        globalTimeStepper->time_loop();

    if (WRITE_TO_SCREEN) cout << "Finished at sim time: " << TIME_STEPPER::SimTime << endl;
}

/**************************************************************/
void SYSTEM::SetupGlobalInfluenceMatrices() {
}

/**************************************************************/
void SYSTEM::UpdateGlobalInfluenceMatrices() {
}

/**************************************************************/
void SYSTEM::GetPressures(REAL dt) {
}

/**************************************************************/
void SYSTEM::GetGlobalRHS() {
}

/**************************************************************/
void SYSTEM::LinAlg() {

}

/**************************************************************/
void SYSTEM::BodySubStep(REAL delta_t, int n_steps) {

    unsigned long int t0 = ticks();
    REAL dt = delta_t / n_steps;

    for (int SubStep = 1; SubStep <= n_steps; ++SubStep) {



        globalTimeStepper->substep_time = SubStep*dt;
        globalTimeStepper->sim_time += dt;

    }

}

/**************************************************************/
void SYSTEM::ReadNeuGetBodies() {

}

/**************************************************************/
void SYSTEM::PrintInfluenceMatrices() {

}

/**************************************************************/
void SYSTEM::PrintBodiesVels() {

}

/**************************************************************/
void SYSTEM::PrintBodiesAndWakes() {

}

/**************************************************************/
void SYSTEM::WriteBodiesAndWakes(ostream& out_stream) {
}

/**************************************************************/
void SYSTEM::PutWakesInTree() {
    for (int J = 0; J < BODY::Bodies.size(); ++J){
        for (int i = 0; i < BODY::Bodies[J]->VortonX.size(); ++i)
            for (int j = 0; j < BODY::Bodies[J]->VortonX[i].size(); ++j){
                for (int k = 0; k < BODY::Bodies[J]->VortonX[i][j].size(); ++k) {
                    Vect3 Om = BODY::Bodies[J]->VortonOM[i][j][k]*(4*globalSystem->GambitScale*globalSystem->GambitScale);
                    Vect3 P = BODY::Bodies[J]->VortonX[i][j][k]*globalSystem->GambitScale;
                    OctreeCapsule C(P,Om, true);
                    C.AssociatedBody = J;
                    globalOctree->Root->EvalCapsule(C);
                }
                BODY::Bodies[J]->VortonX[i].pop_back();
                BODY::Bodies[J]->VortonOM[i].pop_back();
                BODY::Bodies[J]->VortonVel[i].pop_back();
            }

    }
    
}
/**************************************************************/
void SYSTEM::GetFaceVels() {
//    for (int i = 0; i < NumBodies; ++i)
//        for (int j = 0; j < Bodies[i]->NumFaces; ++j)
//            globalOctree->Root->RecursivePanelVel(*(Bodies[i]->Faces[i]));

//    globalOctree->Root->RecursivePassPanelVelsDown();
   
    
    
    
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    //    for (int j = 0; j < globalOctree->AllCells.size(); ++j)
    //        for (int i = 0; i < NumBodies; ++i)
    //#ifdef COLLAPSE_TO_FACES
    //            for (int k = 0; k < 6; ++k)
    //                if ((k == 0) || (k == 2) || (k == 4) || !globalOctree->AllCells[j]->Neighb[k])
    //                    globalOctree->AllCells[j]->FaceVels[k] +=
    //                        Bodies[i]->GetVel(globalOctree->AllCells[j]->Position + Node::NeighbOffset[k]);
    //#else
    //            globalOctree->AllCells[j]->Velocity += Bodies[i]->GetVel(globalOctree->AllCells[j]->Position);
    //#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();

}

/**************************************************************/
void SYSTEM::GetPanelFMMVelocities() {


    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->Vfmm = globalOctree->TreeVel(globalSystem->GambitScale*BODY::AllBodyFaces[i]->Centroid)/globalSystem->GambitScale;
    }
     
     
     
//    for (int i = 0; i < NumBodies; ++i)
//        for (int j = 0; j < Bodies[i]->Faces.size(); ++j)
//            Bodies[i]->Faces[j]->Vfmm =
//                globalOctree->TreeVel(Bodies[i]->Faces[j]->CollocationPoint->vP);

}

/**************************************************************/
void SYSTEM::MoveBodies(REAL dt, bool update) {
//    for (int j = 0; j < NumBodies; ++j)
//        Bodies[j]->MoveBody(dt);
//
//    if (update && NumBodies > 1)
//        UpdateGlobalInfluenceMatrices();
}

/**************************************************************/
void SYSTEM::WriteDomain() {
    //  Here we want to find the extents of the domain
    Vect3 Max(.5, .5, .5), Min(.5, .5, .5);
    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
        Max = max(Max, globalOctree->AllCells[i]->Position);
        Min = min(Min, globalOctree->AllCells[i]->Position);
    }

    Vect3 Sz = Max - Min;
    cout << Max << "\t " << Min << "\t " << Sz << endl;


    int buffer = 10;
    stringstream outstream;
    //  Make the domain slightly larger and get velocities

    for (int i = 0; i < Sz.x + 2 * buffer; ++i)
        for (int j = 0; j < Sz.y + 2 * buffer; ++j)
            for (int k = 0; k < Sz.z + 2 * buffer; ++k) {
                Vect3 X;
                X = Min - buffer + Vect3(i, j, k);
                Vect3 V = globalOctree->TreeVel(X);
                outstream << i + 1 << "\t" << j + 1 << "\t" << k + 1 << "\t" << X << "\t" << V + scaledVinf << endl;
            }

    globalIO->write_file(outstream, string("data"), string("dat"), true);
}

/**************************************************************/
void SYSTEM::WriteVorticity() {
    stringstream outstream;
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        if (globalOctree->AllCells[i]->Omega.Mag() > VORTICITY_CUTOFF) {
            outstream << globalOctree->AllCells[i]->Position + .5 << " " << globalOctree->AllCells[i]->Omega;
            for (int q = 0; q < NumTransVars; ++q)
                outstream << " " << globalOctree->AllCells[i]->TransVars[q];
            outstream << endl;
        }


    globalIO->write_file(outstream, string("f"), string("dat"), true);
}

/**************************************************************/
void SYSTEM::WriteBodies() {
//   
}

/**************************************************************/
void SYSTEM::WritePanelVels() {
//  
}
