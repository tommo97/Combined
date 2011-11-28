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

    while ((TIME_STEPPER::SimTime < (TIME_STEPPER::MaxTime)) && (!globalTimeStepper->last_step))
        globalTimeStepper->time_loop();

    //  Final loop since last_step will exit the previous loop
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



    unsigned long int t0 = ticks();



    int Num2Insert = 0, Num2Keep = 0;

    Array <bool> toInsert(BODY::VortexPositions.size(), false);
    for (int i = 0; i < BODY::VortexPositions.size(); ++i)
        if ((BODY::VortexPositions[i] - *BODY::VortexOrigins[i]).Mag() > (GambitScale * .5)) {
            Num2Insert++;
            toInsert[i] = true;
        } else
            Num2Keep++;


    if (Num2Insert == 0)
        return;


    Array <Vect3> XtoInsert(Num2Insert), XtoKeep(Num2Keep), OMtoInsert(Num2Insert), OMtoKeep(Num2Keep);
    Array <Vect3*> Origins2Keep(Num2Keep);
    Array <int> IDtoInsert(Num2Insert), IDtoKeep(Num2Keep);

    int InsertCount = 0, KeepCount = 0;
    for (int i = 0; i < BODY::VortexPositions.size(); ++i)
        if (toInsert[i]) {
            XtoInsert[InsertCount] = BODY::VortexPositions[i];
            OMtoInsert[InsertCount] = BODY::VortexOmegas[i];
            IDtoInsert[InsertCount] = BODY::VortexOwnerID[i];
            InsertCount++;
        } else {
            XtoKeep[KeepCount] = BODY::VortexPositions[i];
            OMtoKeep[KeepCount] = BODY::VortexOmegas[i];
            IDtoKeep[KeepCount] = BODY::VortexOwnerID[i];
            Origins2Keep[KeepCount] = BODY::VortexOrigins[i];
            KeepCount++;
        }



    

    
    REAL d = 1.0, u = 0;
    REAL M4x = 0.0, M4y = 0.0, M4z = 0.0, M4 = 0;
    int count = 0;
    for (REAL a = -h; a <= h; a += 1.0)
        count++;


    Array <Vect3> Xs(Num2Insert * (count * count * count)), Oms(Num2Insert * (count * count * count));
    Array <int> Owners(Num2Insert * (count * count * count));

    count = 0;

    for (int j = 0; j < XtoInsert.size(); ++j) {
        Vect3 X = XtoInsert[j];
        Vect3 OM = OMtoInsert[j];
        int ID = IDtoInsert[j];
        //        OctreeCapsule C(X, OM, true);
        //        C.AssociatedBody = BODY::VortexOwnerID[j] - 1;
        //        globalOctree->Root->EvalCapsule(C);

        //
        for (REAL a = -h; a <= h; a += 1.0) {
            M4x = 0.0;
            u = abs(a / (h * d));
            if (u <= 2)
                M4x = .5 * (1 - u * (2 - u)*(2 - u));
            if (u <= 1)
                M4x = 1 - 2.5 * u * u + 1.5 * u * u * u;
            M4x = M4x / h;
            for (REAL b = -h; b <= h; b += 1.0) {
                M4y = 0.0;
                u = abs(b / (h * d));
                if (u <= 2)
                    M4y = .5 * (1 - u * (2 - u)*(2 - u));
                if (u <= 1)
                    M4y = 1 - 2.5 * u * u + 1.5 * u * u * u;
                M4y = M4y / h;
                for (REAL c = -h; c <= h; c += 1.0) {

                    u = abs(c / (h * d));
                    if (u <= 2)
                        M4z = .5 * (1 - u * (2 - u)*(2 - u));
                    if (u <= 1)
                        M4z = 1 - 2.5 * u * u + 1.5 * u * u * u;

                    M4z = M4z / h;
                    M4 = M4x * M4y*M4z;

                    Vect3 G = M4*OM;
                    Xs[count] = Vect3(floor(X.x + a) + 0.5, floor(X.y + b) + 0.5, floor(X.z + c) + 0.5);
                    Oms[count] = (G);
                    Owners[count] = (ID - 1);
                    count++;
                }
            }
        }
    }





    Array <OctreeCapsule> Test(Xs.size());

    for (int i = 0; i < Test.size(); ++i) {
        Test[i] = OctreeCapsule(Xs[i], Oms[i], true);
        Test[i].AssociatedBody = Owners[i];
    }


    Array<OctreeCapsule>::QuickSortB(Test);
    Vect3 min_diff = 1e16;
    int num_unique = 1;
    for (int i = 1; i < Test.size(); ++i) {
        if (Test[i] < Test[i - 1])
            cout << "Error in sort..." << endl;

        min_diff = min(min_diff, Test[i].Position - Test[i - 1].Position);

        if (Test[i] != Test[i - 1])
            num_unique++;
    }


    Array <OctreeCapsule> unique_data(num_unique);
    num_unique = 1;
    unique_data[num_unique - 1] = Test[0];

    for (int i = 1; i < Test.size(); ++i) {
        if (Test[i] != Test[i - 1]) {
            num_unique++;
            unique_data[num_unique - 1] = Test[i];
        } else
            unique_data[num_unique - 1].Omega += Test[i].Omega;
    }

    for (int i = 0; i < unique_data.size(); ++i) {
        globalOctree->Root->EvalCapsule(unique_data[i]);
    }

    unsigned long int t1 = ticks();
//    cout << "----- " << Num2Insert << " " << Num2Keep << " " << Test.size() << " " << num_unique << " " << t1-t0 << endl;


    BODY::VortexPositions = XtoKeep;
    BODY::VortexOmegas = OMtoKeep;
    BODY::VortexOwnerID = IDtoKeep;
    BODY::VortexOrigins = Origins2Keep;


}
/**************************************************************/
void SYSTEM::GetFaceVels() {

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
        globalOctree->AllCells[i]->PanelVel = Vect3(0.0);

        for (int j = 0; j < BODY::AllBodyFaces.size(); ++j)
            globalOctree->AllCells[i]->PanelVel += BODY::AllBodyFaces[j]->BodyPanelVelocity(globalOctree->AllCells[i]->Position);


        for (int j = 0; j < BODY::AllProtoWakes.size(); ++j)
            globalOctree->AllCells[i]->PanelVel += BODY::AllProtoWakes[j]->VortexPanelVelocity(globalOctree->AllCells[i]->Position);

        globalOctree->AllCells[i]->Velocity += globalOctree->AllCells[i]->PanelVel;
        
        globalOctree->AllCells[i]->SetVelsEqual();
    }
}

/**************************************************************/
void SYSTEM::GetPanelFMMVelocities(REAL dt) {
    
    int sz = BODY::AllBodyFaces.size();
    Array <Vect3> P1(sz), P2(sz), V2(sz), V1(sz);

#pragma omp parallel for
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
        P1[i] = BODY::AllBodyFaces[i]->CollocationPoint;


    BODY::Time += dt;
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        BODY::Bodies[i]->MoveBody();

#pragma omp parallel for
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
        P2[i] = BODY::AllBodyFaces[i]->CollocationPoint;

    
    BODY::Time -= dt;
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        BODY::Bodies[i]->MoveBody();
    

    V1 = Vect3(0.0);
    V2 = V1;
//#pragma omp parallel for
//    for (int i = 0; i < sz; ++i) {
//        for (int j = 0; j < globalOctree->AllCells.size(); ++j) {
//            V1[i] += globalDirectVel(P1[i] - globalOctree->AllCells[j]->Position,
//                    globalOctree->AllCells[j]->Omega);
//            V2[i] += globalDirectVel(P2[i] - globalOctree->AllCells[j]->Position,
//                    globalOctree->AllCells[j]->Omega);
//        }
//    }

    for (int i = 0; i < sz; ++i)
        V1[i] = globalOctree->TreeVel(P1[i]);

    for (int i = 0; i < sz; ++i)
        V2[i] = globalOctree->TreeVel(P2[i]);
    
    
        
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i){
            BODY::AllBodyFaces[i]->dVFMM_dt = (1/dt) * (V2[i]-V1[i]);
            BODY::AllBodyFaces[i]->Vfmm = V1[i];
        }
        

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
void SYSTEM::WriteData() {

   MATLABOutputStruct Output;
   
   
    Array < Array < REAL > > DataOut = UTIL::zeros(globalOctree->AllCells.size(), 6 + (NumTransVars * 3));
    int count = 0;
    Vect3 Mins(1e32,1e32,1e32), Maxs(-1e32,-1e32,-1e32);
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
            
            Mins = min(globalOctree->AllCells[i]->Position,Mins);
            Maxs = max(globalOctree->AllCells[i]->Position,Maxs);
            
        }

    cout << "Domain Size: " << Maxs - Mins << endl;
    
    Array < Array <Vect3> > BodyVelSliceY, TreeVelSliceY,  Posns;
    int imax = 1, jmax = 1;
    for (int i = (int) Mins.x - 48; i < (int) Maxs.x + 48; ++i)
        imax++;
    for (int j = (int) Mins.z - 48; j < (int) Maxs.z + 48; ++j)
        jmax++;

    
    Posns = TreeVelSliceY = BodyVelSliceY = UTIL::zerosv(imax,jmax);
    for (int i = 0; i < imax; ++i)
        for (int j = 0; j < jmax; ++j)
        {
            Posns[i][j] = Mins + Vect3(i-48.0,0.0,j-48.0);
            TreeVelSliceY[i][j] = globalOctree->TreeVel(Posns[i][j]);
            
            for (int k = 0; k < BODY::AllBodyFaces.size(); ++k)
            {
                BodyVelSliceY[i][j] += BODY::AllBodyFaces[k]->SourceVel(Posns[i][j]);
                BodyVelSliceY[i][j] += BODY::AllBodyFaces[k]->BodyPanelVelocity(Posns[i][j]);
            }
        }
    
    


    Output.Double2DArrays.push_back(DataOut);
    Output.Double2DArrayStrings.push_back(string("Domain"));

    Output.Double2DArrays.push_back(BODY::CpHistoryAll);
    Output.Double2DArrayStrings.push_back(string("CpHistoryAll"));

    Output.Double2DArrays.push_back(BODY::CpHistoryAllD);
    Output.Double2DArrayStrings.push_back(string("CpHistoryAllD"));


    Output.Double1DArrays.push_back(BODY::SubTIMES);
    Output.Double1DArrayStrings.push_back(string("CpHistoryAllD"));

    Output.Vect2DArrays.push_back(TreeVelSliceY);
    Output.Vect2DArrayStrings.push_back(string("TreeVelSliceY"));

    Output.Vect2DArrays.push_back(BodyVelSliceY);
    Output.Vect2DArrayStrings.push_back(string("BodyVelSliceY"));  
    
    Output.Vect2DArrays.push_back(Posns);
    Output.Vect2DArrayStrings.push_back(string("SlicePositions"));  
    
    globalIO->writeMATLABOutputStruct(Output, string("RunData"), true);
    BODY::CpHistoryAll.clear();
    BODY::CpHistoryAllD.clear();
    BODY::SubTIMES.clear();
}

/**************************************************************/
void SYSTEM::WriteBodies() {
//   
}

/**************************************************************/
void SYSTEM::WritePanelVels() {
//  
}
