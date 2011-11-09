/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 35               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-16 00:1#$:  Date of last commit

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

    cout << "----------------------- " << TIME_STEPPER::SimTime  << " " << TIME_STEPPER::MaxTime << endl;
    while (TIME_STEPPER::SimTime < (TIME_STEPPER::MaxTime - 1e-3))
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

    int n = (int) GambitScale;
    REAL h = 5, d = 1.0, u = 0;
    REAL M4x = 0.0, M4y = 0.0, M4z = 0.0, M4 = 0;;
  

    for (int j = 0; j < BODY::VortexPositions.size(); ++j) {
        Vect3 X = BODY::VortexPositions[j] * GambitScale;
        Vect3 OM = -1.0 * BODY::VortexOmegas[j] * GambitScale * GambitScale;

        for (REAL a = -h; a <= h; a+=1.0)
            for (REAL b = -h; b <= h; b+=1.0)
                for (REAL c = -h; c <= h; c+=1.0) {

                    M4x = 0.0, M4y = 0.0, M4z = 0.0;

                    u = abs(a / (h * d));
                    if (u <= 2)
                        M4x = .5 * (1 - u * (2 - u)*(2 - u));
                    if (u <= 1)
                        M4x = 1 - 2.5 * u * u + 1.5 * u * u * u;

                    u = abs(b / (h * d));
                    if (u <= 2)
                        M4y = .5 * (1 - u * (2 - u)*(2 - u));
                    if (u <= 1)
                        M4y = 1 - 2.5 * u * u + 1.5 * u * u * u;

                    u = abs(c / (h * d));
                    if (u <= 2)
                        M4z = .5 * (1 - u * (2 - u)*(2 - u));
                    if (u <= 1)
                        M4z = 1 - 2.5 * u * u + 1.5 * u * u * u;


                    M4x = M4x / h;
                    M4y = M4y / h;
                    M4z = M4z / h;
                    M4 = M4x*M4y*M4z;

                    Vect3 G = M4*OM;

                    OctreeCapsule C(X + Vect3(a, b, c), G, true);
                    C.AssociatedBody = BODY::VortexOwnerID[j] - 1;
                    globalOctree->Root->EvalCapsule(C);
                }
    }

    BODY::VortexPositions.clear();
    BODY::VortexOmegas.clear();
    BODY::VortexOwnerID.clear();


}
/**************************************************************/
void SYSTEM::GetFaceVels() {
    for (int i = 0; i < Node::NumNodes; ++i)
        Node::AllNodes[i]->PanelVel = Vect3(0.0);

//    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
//        globalOctree->Root->RecursivePanelVel(BODY::AllBodyFaces[i]);
//
//    globalOctree->Root->RecursivePassPanelVelsDown();
//
//    for (int i = 0; i < globalOctree->AllCells.size(); ++i){
////        cout << globalOctree->AllCells[i]->Velocity << " " << globalOctree->AllCells[i]->PanelVel << endl;
//        globalOctree->AllCells[i]->Velocity += globalOctree->AllCells[i]->PanelVel;
//    }
    
    
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
void SYSTEM::GetPanelFMMVelocities(REAL dt) {


    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->Vfmm0 = globalOctree->TreeVel(globalSystem->GambitScale * BODY::AllBodyFaces[i]->Centroid);
    }

    BODY::Time += dt;
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        BODY::Bodies[i]->MoveBody();

    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->Vfmm1 = globalOctree->TreeVel(globalSystem->GambitScale * BODY::AllBodyFaces[i]->Centroid);
    }
    BODY::Time -= dt;
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        BODY::Bodies[i]->MoveBody();
    
    
//    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) 
//        cout << "D " << BODY::AllBodyFaces[i]->VWake << endl << "F " << BODY::AllBodyFaces[i]->Vfmm0 << endl << "R " << BODY::AllBodyFaces[i]->Vfmm0/BODY::AllBodyFaces[i]->VWake << endl << "---" << endl;
        
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


   
    Array < Array < REAL > > DataOut = UTIL::zeros(globalOctree->AllCells.size(), 6 + (NumTransVars * 3));
    int count = 0;
    for (int i = 0; i < globalOctree->AllCells.size(); ++i){
            DataOut[count][0] = globalOctree->AllCells[i]->Position.x + 0.5;
            DataOut[count][1] = globalOctree->AllCells[i]->Position.y + 0.5;
            DataOut[count][2] = globalOctree->AllCells[i]->Position.z + 0.5;
            DataOut[count][3] = globalOctree->AllCells[i]->Omega.x;
            DataOut[count][4] = globalOctree->AllCells[i]->Omega.y;
            DataOut[count][5] = globalOctree->AllCells[i]->Omega.z;
            int cnt = 6;
            for (int q = 0; q < NumTransVars; ++q) {
                DataOut[count][cnt] = globalOctree->AllCells[i]->TransVars[q].x;
                cnt++;
                DataOut[count][cnt] = globalOctree->AllCells[i]->TransVars[q].y;
                cnt++;
                DataOut[count][cnt] = globalOctree->AllCells[i]->TransVars[q].z;
                cnt++;
            }
            count++;
        }

    globalIO->write_2D_mat(DataOut, string("Domain"), string("Domain"), true);


}

/**************************************************************/
void SYSTEM::WriteBodies() {
//   
}

/**************************************************************/
void SYSTEM::WritePanelVels() {
//  
}
