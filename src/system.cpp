/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2013
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velocity vorticity form

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

int SYSTEM::MaxP, SYSTEM::NumTransVars, SYSTEM::M4Radius_in_cells;
REAL SYSTEM::GambitScale;
/**************************************************************/
SYSTEM::~SYSTEM() {

    delete globalIO;
    delete globalOctree;
    delete globalTimeStepper;

}

/**************************************************************/
SYSTEM::SYSTEM(int NT) {
    LiftingLineMode = false;
    useBodies = true;
    num_out = 0;
    SysDumpInterval = 0; //  Used to dump the result of system calls
    scaledVinf.x = scaledVinf.y = scaledVinf.z = 0;
    unscaledVinf = scaledVinf;
    SYSTEM::MaxP = 3;
    DS = .3;
    dtInit = 0.0;
    NumSubSteps = 0;
    g = 9.80665; // m/s/s
    Temp = 288.15; //  Kelvin
    Rho = 1027; //1.226; //1027; //  Kg/m3
    Mu = 1.5e-3; //(sqrt(pow(Temp, 3)) * 1.458e-6) / (Temp + 110.4); //   Dynamic Viscocity
    SYSTEM::GambitScale = 1;
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

    useFMM = true;
    scaledVinf = unscaledVinf*SYSTEM::GambitScale;
    Nu = SYSTEM::GambitScale * SYSTEM::GambitScale * Mu / Rho; // <- this is Reynolds scaling, for constant rho
    globalTimeStepper = new TIME_STEPPER();

    if (globalSystem->useBodies)
        SYSTEM::NumTransVars = max(1, BODY::Bodies.size());


    globalOctree = new OCTREE();

    if (WRITE_TO_SCREEN) cout << "rho: " << Rho << "; mu: " << Mu << "; nu: " << Nu << endl;
    if (WRITE_TO_SCREEN) cout << "FVM mesh scale Factor: " << SYSTEM::GambitScale << endl;

    int nss = globalSystem->NumSubSteps;
    globalTimeStepper->dt_prev = globalSystem->dtInit;
    globalTimeStepper->time_step();

    globalSystem->NumSubSteps = nss;


    cout << "nss" << nss << endl;
}

/**************************************************************/
void SYSTEM::TimeStep() {

    while ((TIME_STEPPER::SimTime < (TIME_STEPPER::MaxTime)) && (!globalTimeStepper->last_step)) {
        globalTimeStepper->time_loop();
    }
    //  Final loop since last_step will exit the previous loop
    globalTimeStepper->time_loop();

    if (WRITE_TO_SCREEN) cout << "Finished at sim time: " << TIME_STEPPER::SimTime << endl;
}

/**************************************************************/
void SYSTEM::AddVortonsToTree(Array <Vect3> &XtoInsert, Array <Vect3> &OMtoInsert, Array <int> &IDtoInsert) {

    int Num2Insert = XtoInsert.size();
    REAL d = 1.0, u = 0;
    REAL M4x = 0.0, M4y = 0.0, M4z = 0.0, M4 = 0;
    int count = 0, h = SYSTEM::M4Radius_in_cells;

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
                    Xs[count] = Vect3(round(X.x + a), round(X.y + b), round(X.z + c));
                    Oms[count] = (G);
                    Owners[count] = (ID);
                    count++;
                }
            }
        }
    }





    Array <OctreeCapsule> Test(Xs.size());

    for (int i = 0; i < Test.size(); ++i) {
        Test[i] = OctreeCapsule(Xs[i], -1.0 * Oms[i], true);
        Test[i].AssociatedBody = Owners[i];
    }


    for (int i = 0; i < Test.size(); ++i) {
        globalOctree->Root->EvalCapsule(Test[i]);
    }




}

/**************************************************************/
void SYSTEM::PutWakesInTree() {
#ifdef TIME_STEPS
    long unsigned int t1 = ticks();
#endif
    int Num2Insert = 0, Num2Keep = 0;

    Array <bool> toInsert(BODY::VortexPositions.size(), false);
    for (int i = 0; i < BODY::VortexPositions.size(); ++i)
        if ((BODY::VortexPositions[i] - *BODY::VortexOrigins[i]).Mag() > (0.0625 * SYSTEM::GambitScale)) {
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
    int count = 0, h = SYSTEM::M4Radius_in_cells;
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
                    Xs[count] = Vect3(round(X.x + a), round(X.y + b), round(X.z + c));
                    Oms[count] = (G);
                    Owners[count] = (ID - 1);
                    count++;
                }
            }
        }
    }





    Array <OctreeCapsule> Test(Xs.size());

    for (int i = 0; i < Test.size(); ++i) {
        Test[i] = OctreeCapsule(Xs[i], 1.0 * Oms[i], true);
        Test[i].AssociatedBody = Owners[i];
    }


    for (int i = 0; i < Test.size(); ++i) {
        globalOctree->Root->EvalCapsule(Test[i]);
    }
    //    Array<OctreeCapsule>::QuickSortB(Test);
    //    Vect3 min_diff = 1e16;
    //    int num_unique = 1;
    //    for (int i = 1; i < Test.size(); ++i) {
    //        if (Test[i] < Test[i - 1])
    //            cout << "Error in sort..." << endl;
    //
    //        min_diff = min(min_diff, Test[i].Position - Test[i - 1].Position);
    //
    //        if (Test[i] != Test[i - 1])
    //            num_unique++;
    //    }
    //
    //
    //    Array <OctreeCapsule> unique_data(num_unique);
    //    num_unique = 1;
    //    unique_data[num_unique - 1] = Test[0];
    //
    //    for (int i = 1; i < Test.size(); ++i) {
    //        if (Test[i] != Test[i - 1]) {
    //            num_unique++;
    //            unique_data[num_unique - 1] = Test[i];
    //        } else
    //            unique_data[num_unique - 1].Omega += Test[i].Omega;
    //    }
    //    cout << "I1 " << globalIO->ReturnMemPercent() << endl;
    //    for (int i = 0; i < unique_data.size(); ++i) {
    //        globalOctree->Root->EvalCapsule(unique_data[i]);
    //    }
    //    cout << "I2 " << globalIO->ReturnMemPercent() << endl;



    BODY::VortexPositions = XtoKeep;
    BODY::VortexOmegas = OMtoKeep;
    BODY::VortexOwnerID = IDtoKeep;
    BODY::VortexOrigins = Origins2Keep;

#ifdef TIME_STEPS
    long unsigned int t2 = ticks();

    stringstream tmp;
    tmp << "Bin panel wake into tree : " << double(t2 - t1) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif

}

/**************************************************************/
void SYSTEM::GetFaceVels() {
#ifdef TIME_STEPS
    long unsigned int t6 = ticks();
#endif
    //    if (globalSystem->useBodies) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    //    for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
    //        for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
    //            Vect3 Target = FVMCell::AllCells[i]->Position;
    //            //      The following are divided by 2 since they are the free space rather than the surface velocities
    //            Vect3 VelS = BODY::AllBodyFaces[j]->Sigma * BODY::AllBodyFaces[j]->SourcePanelVelocity(Target) / 2.0;
    //            Vect3 VelD = BODY::AllBodyFaces[j]->Mu * BODY::AllBodyFaces[j]->DoubletPanelVelocity(Target) / 2.0;
    //            FVMCell::AllCells[i]->Velocity += (VelS - VelD);
    //        }
    //    }
    //    }



    for (int i = 0; i < FVMCell::AllCells.size(); ++i)
        FVMCell::AllCells[i]->SetVelsEqual();


    //    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
    //        Vect3 Vel(0.,0.,0.);
    //
    //        for (int j = 0; j < FVMCell::AllCells.size(); ++j)
    //            globalDirectVel(FVMCell::AllCells[j]->Position - FVMCell::AllCells[i]->Position, FVMCell::AllCells[j]->Omega, Vel);
    //
    //        cout << Vel << " " << FVMCell::AllCells[i]->Velocity << endl;
    //
    //    }
#ifdef TIME_STEPS
    long unsigned int t7 = ticks();
    stringstream tmp;
    tmp << "GetFaceVels()            : " << double(t7 - t6) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void SYSTEM::GetPanelFMMVelocities(REAL dt) {
#ifdef TIME_STEPS
    long unsigned int t9 = ticks();
#endif


#pragma omp parallel for
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->Xfmm0 = BODY::AllBodyFaces[i]->CollocationPoint;
    }
    if (dt > 0) {
        BODY::Time += dt;
        for (int i = 0; i < BODY::Bodies.size(); ++i)
            BODY::Bodies[i]->MoveBody();

#pragma omp parallel for
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
            BODY::AllBodyFaces[i]->Xfmm1 = BODY::AllBodyFaces[i]->CollocationPoint;
        }


        BODY::Time -= dt;
        for (int i = 0; i < BODY::Bodies.size(); ++i)
            BODY::Bodies[i]->MoveBody();
    }

#pragma omp parallel for
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        if (dt > 0)
            BODY::AllBodyFaces[i]->Vfmm1 = UTIL::interp3Pointer <Vect3 > (BODY::AllBodyFaces[i]->Xp, BODY::AllBodyFaces[i]->Vp, BODY::AllBodyFaces[i]->Xfmm1);
        else
            BODY::AllBodyFaces[i]->Vfmm0 = UTIL::interp3Pointer <Vect3 > (BODY::AllBodyFaces[i]->Xp, BODY::AllBodyFaces[i]->Vp, BODY::AllBodyFaces[i]->Xfmm0);
    }

    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        if (dt > 0)
            BODY::AllBodyFaces[i]->dVFMM_dt = (1.0 / dt) * (BODY::AllBodyFaces[i]->Vfmm1 - BODY::AllBodyFaces[i]->Vfmm0);
        else
            BODY::AllBodyFaces[i]->Vfmm = BODY::AllBodyFaces[i]->Vfmm0;
    }
#ifdef TIME_STEPS
    long unsigned int t10 = ticks();
    stringstream tmp;
    tmp << "GetPanelFMMVelocities(.) : " << double(t10 - t9) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif


    


    //    cout << "--------------- Calcing Inerp Vels --------------" << endl;
    //    cout << FVMCell::AllCells.size()*BODY::AllBodyFaces.size() << " Interactions " << endl;
    //    Array < Array < Array < Vect3 > > > Posns, Vels;
    //            
    //    Posns = Array < Array < Array <Vect3> > > (DX, Array < Array < Vect3 > > (DY, Array < Vect3 > (DZ, Vect3(0.0))));
    //    Vels = Posns;
    //    for (int i = 0; i < DX; ++i)
    //        for (int j = 0; j < DY; ++j)
    //            for (int k = 0; k < DZ; ++k) {
    //                Vect3 XP(MinX + i, MinY + j, MinZ + k);
    //                Posns[i][j][k] = XP;
    //                        #pragma omp parallel for
    //                                for (int l = 0; l < FVMCell::AllCells.size(); ++l)
    //                                    Vels[i][j][k] += UTIL::globalDirectVel(FVMCell::AllCells[l]->Position - XP, FVMCell::AllCells[l]->Omega, globalSystem->Del2);
    ////                Vels[i][j][k] = globalOctree->TreeVel(XP);
    //            }

    //#ifdef TIME_STEPS
    //        long unsigned int t11 = ticks();
    //        cout << "GetPanelFMMVelocities(.) : " << double(t10 - t9) / 1000.0 << " Direct/Interp : " << double(t11 - t10) / 1000.0 << endl;
    //#endif
}


/**************************************************************/
void SYSTEM::WriteDomain() {
    //  Here we want to find the extents of the domain
    Vect3 Max(.5, .5, .5), Min(.5, .5, .5);
    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
        Max = max(Max, FVMCell::AllCells[i]->Position);
        Min = min(Min, FVMCell::AllCells[i]->Position);
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
#ifdef TIME_STEPS
    unsigned long int t1 = ticks();
#endif
    MATLABOutputStruct Output;

    Vect3 Maxs, Mins;
    Maxs = -1e32;
    Mins = 1e32;
    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
        Maxs.x = max(FVMCell::AllCells[i]->Position.x, Maxs.x);
        Maxs.y = max(FVMCell::AllCells[i]->Position.y, Maxs.y);
        Maxs.z = max(FVMCell::AllCells[i]->Position.z, Maxs.z);

        Mins.x = min(FVMCell::AllCells[i]->Position.x, Mins.x);
        Mins.y = min(FVMCell::AllCells[i]->Position.y, Mins.y);
        Mins.z = min(FVMCell::AllCells[i]->Position.z, Mins.z);

    }

//    int sx = (Maxs.x - Mins.x) + 1;
//    int sy = (Maxs.y - Mins.y) + 1;
//    int sz = (Maxs.z - Mins.z) + 1;


    //    cout << sx << " " << sy << " " << sz << " " << Maxs << " " << Mins << endl;
    //    
    //    int cnt = 0;
    //    REAL tmpx = Mins.x;
    //    while (tmpx <= Maxs.x) {
    //        REAL tmpy = Mins.y;
    //        while (tmpy <= Maxs.y) {
    //            REAL tmpz = Mins.z;
    //            while (tmpz <= Maxs.z) {
    //                tmpz += 1;
    //                cnt++;
    //            }
    //            tmpy += 1;
    //        }
    //        tmpx += 1;
    //    }
    //    cout << cnt << endl;
    //    Array <Vect3> Positions(cnt), Velocities(cnt);
    //    cnt = 0;
    //    tmpx = Mins.x;
    //    while (tmpx <= Maxs.x) {
    //        REAL tmpy = Mins.y;
    //        while (tmpy <= Maxs.y) {
    //            REAL tmpz = Mins.z;
    //            while (tmpz <= Maxs.z) {
    //                Positions[cnt] = Vect3(tmpx,tmpy,tmpz);
    //                tmpz += 1;
    //                cnt++;
    //            }
    //            tmpy += 1;
    //        }
    //        tmpx += 1;
    //    }
    //    
    //    
    //    for (int i = 0; i < Positions.size(); ++i)
    //        Velocities[i] = globalOctree->TreeVel(Positions[i]);
    //    
    //    
    //
    //    Output.Vect1DArrays.push_back(Positions);
    //    Output.Vect1DArrayStrings.push_back(string("CellPositions"));
    //
    //    Output.Vect1DArrays.push_back(Velocities);
    //    Output.Vect1DArrayStrings.push_back(string("CellVelocities"));   



    //    REAL DistTravelled = (BODY::Bodies[0]->CG - BODY::Bodies[0]->CGo).Mag();
    //    REAL tmpX = BODY::Bodies[0]->CG.x - (2.0 * SYSTEM::GambitScale * 0.4);
    //    int nSlices = 2 + floor(DistTravelled / (SYSTEM::GambitScale * 0.4));
    //    int I = 0;
    //    Array <REAL> SlicePlanesX(7);
    //
    //    SlicePlanesX[0] = 0.1 * SYSTEM::GambitScale;
    //    SlicePlanesX[1] = 0.2 * SYSTEM::GambitScale;
    //    SlicePlanesX[2] = 0.5 * SYSTEM::GambitScale;
    //    SlicePlanesX[3] = 1.0 * SYSTEM::GambitScale;
    //    SlicePlanesX[4] = 2.0 * SYSTEM::GambitScale;
    //    SlicePlanesX[5] = 4.0 * SYSTEM::GambitScale;
    //    SlicePlanesX[6] = 6.0 * SYSTEM::GambitScale;
    //
    //    Array <string> Names(7);
    //    Names[0] = "0pt1";
    //    Names[1] = "0pt2";
    //    Names[2] = "0pt5";
    //    Names[3] = "1pt0";
    //    Names[4] = "2pt0";
    //    Names[5] = "4pt0";
    //    Names[6] = "6pt0";
    //
    //    int NC = 0;
    //
    //    int JKmin = floor(-5 * SYSTEM::GambitScale), JKmax = ceil(5 * SYSTEM::GambitScale);
    //    Array <REAL> PS;
    //    for (int J = JKmin; J < JKmax; ++J) {
    //        NC++;
    //        PS.push_back(J * 1.0);
    //    }
    //
    //    for (int I = 0; I < SlicePlanesX.size(); ++I) {
    //        cout << I << endl;
    //        Array < Array < Vect3 > > SliceWakeVel = UTIL::zerosv(NC, NC);
    //        Array < Array < Vect3 > > SliceVortonWakeVel = UTIL::zerosv(NC, NC);
    //        Array < Array < Vect3 > > SliceBodyVel = UTIL::zerosv(NC, NC);
    //        Array < Array < Vect3 > > SliceProtoWakeVel = UTIL::zerosv(NC, NC);
    //
    //        Array < Array < Vect3 > > SlicePosn = UTIL::zerosv(NC, NC);
    //        Array < Array < Vect3 > > SlicePosnPlus = UTIL::zerosv(NC, NC);
    //        Array < Array < Vect3 > > SlicePosnMinus = UTIL::zerosv(NC, NC);
    //
    //        Array < Array < REAL > > SliceBodyPhi = UTIL::zeros(NC, NC);
    //        Array < Array < REAL > > SliceBodyPhiPlus = UTIL::zeros(NC, NC);
    //        Array < Array < REAL > > SliceBodyPhiMinus = UTIL::zeros(NC, NC);
    //
    //        Array < Array < REAL > > SliceProtoWakePhi = UTIL::zeros(NC, NC);
    //        Array < Array < REAL > > SliceProtoWakePhiPlus = UTIL::zeros(NC, NC);
    //        Array < Array < REAL > > SliceProtoWakePhiMinus = UTIL::zeros(NC, NC);
    //
    //        for (int J = 0; J < PS.size(); ++J)
    //            for (int K = 0; K < PS.size(); ++K) {
    //                Vect3 Pos(BODY::Bodies[0]->CG.x + SlicePlanesX[I], PS[J], PS[K]);
    //
    //                SlicePosn[J][K] = Pos;
    //                SlicePosnPlus[J][K] = Pos + Vect3(1, 0, 0);
    //                SlicePosnMinus[J][K] = Pos - Vect3(1, 0, 0);
    //
    //                SliceWakeVel[J][K] = globalOctree->TreeVel(Pos);
    //
    //            }
    //
    //#pragma omp parallel for
    //        for (int J = 0; J < PS.size(); ++J)
    //            for (int K = 0; K < PS.size(); ++K) {
    //                Vect3 Pos = SlicePosn[J][K];
    //                for (int k = 0; k < BODY::VortexPositions.size(); ++k)
    //                    SliceVortonWakeVel[J][K] += globalDirectVel(Pos - BODY::VortexPositions[k], BODY::VortexOmegas[k]);
    //
    //                for (int k = 0; k < BODY::AllProtoWakes.size(); ++k)
    //                    SliceProtoWakeVel[J][K] += BODY::AllProtoWakes[k]->VortexPanelVelocity(Pos);
    //                for (int k = 0; k < BODY::AllBodyFaces.size(); ++k)
    //                    SliceBodyVel[J][K] += BODY::AllBodyFaces[k]->BodyPanelVelocity(Pos);
    //
    //
    //
    //                for (int k = 0; k < BODY::AllProtoWakes.size(); ++k) {
    //                    Pos = SlicePosn[J][K];
    //                    SliceProtoWakePhi[J][K] += BODY::AllProtoWakes[k]->WakePanelPotential(Pos);
    //                    Pos = SlicePosnPlus[J][K];
    //                    SliceProtoWakePhiPlus[J][K] += BODY::AllProtoWakes[k]->WakePanelPotential(Pos);
    //                    Pos = SlicePosnMinus[J][K];
    //                    SliceProtoWakePhiMinus[J][K] += BODY::AllProtoWakes[k]->WakePanelPotential(Pos);
    //                }
    //
    //                for (int k = 0; k < BODY::AllBodyFaces.size(); ++k) {
    //                    Pos = SlicePosn[J][K];
    //                    SliceBodyPhi[J][K] += BODY::AllBodyFaces[k]->BodyPanelPotential(Pos);
    //                    Pos = SlicePosnPlus[J][K];
    //                    SliceBodyPhiPlus[J][K] += BODY::AllBodyFaces[k]->BodyPanelPotential(Pos);
    //                    Pos = SlicePosnMinus[J][K];
    //                    SliceBodyPhiMinus[J][K] += BODY::AllBodyFaces[k]->BodyPanelPotential(Pos);
    //                }
    //            }
    //
    //        string SliceName = Names[I];
    //        Output.Vect2DArrays.push_back(SliceWakeVel);
    //        Output.Vect2DArrayStrings.push_back(string("SliceWakeVel" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SliceVortonWakeVel);
    //        Output.Vect2DArrayStrings.push_back(string("SliceVortonWakeVel" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SliceBodyVel);
    //        Output.Vect2DArrayStrings.push_back(string("SliceBodyVel" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SliceProtoWakeVel);
    //        Output.Vect2DArrayStrings.push_back(string("SliceProtoWakeVel" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SlicePosn);
    //        Output.Vect2DArrayStrings.push_back(string("SlicePosn" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceProtoWakePhi);
    //        Output.Double2DArrayStrings.push_back(string("SliceProtoWakePhi" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceProtoWakePhiPlus);
    //        Output.Double2DArrayStrings.push_back(string("SliceProtoWakePhiPlus" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceProtoWakePhiMinus);
    //        Output.Double2DArrayStrings.push_back(string("SliceProtoWakePhiMinus" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceBodyPhi);
    //        Output.Double2DArrayStrings.push_back(string("SliceBodyPhi" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceBodyPhiPlus);
    //        Output.Double2DArrayStrings.push_back(string("SliceBodyPhiPlus" + SliceName));
    //
    //        Output.Double2DArrays.push_back(SliceBodyPhiMinus);
    //        Output.Double2DArrayStrings.push_back(string("SliceBodyPhiMinus" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SlicePosnMinus);
    //        Output.Vect2DArrayStrings.push_back(string("SlicePosnMinus" + SliceName));
    //
    //        Output.Vect2DArrays.push_back(SlicePosnPlus);
    //        Output.Vect2DArrayStrings.push_back(string("SlicePosnPlus" + SliceName));
    //    }
    //

    //  DataOut is vectors of: Position Omega [Transvars1 ... TransvarsN] CFL DerivConv DerivVisc DerivStretch DerivArt01


    Array < Array < Vect3 > > TransVars(FVMCell::AllCells.size(), Array <Vect3 > (SYSTEM::NumTransVars, Vect3(0.0)));
    Array < Vect3 > CellPos(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellOms(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellVel(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellCFL(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellConvDeriv(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellTiltDeriv(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellViscDeriv(FVMCell::AllCells.size(), Vect3(0.0));
    Array < Vect3 > CellArtDeriv(FVMCell::AllCells.size(), Vect3(0.0));
    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
        CellPos[i] = FVMCell::AllCells[i]->Position + Vect3(0.5, 0.5, 0.5);
        CellOms[i] = FVMCell::AllCells[i]->Omega;
        CellVel[i] = FVMCell::AllCells[i]->Velocity;
        for (int j = 0; j < SYSTEM::NumTransVars; ++j)
            TransVars[i][j] = FVMCell::AllCells[i]->TransVars[j];
    }

    Output.Vect1DArrays.push_back(CellPos);
    Output.Vect1DArrayStrings.push_back(string("CellPos"));

    Output.Vect1DArrays.push_back(CellOms);
    Output.Vect1DArrayStrings.push_back(string("CellOms"));

    Output.Vect1DArrays.push_back(CellVel);
    Output.Vect1DArrayStrings.push_back(string("CellVel"));

    Output.Vect1DArrays.push_back(CellConvDeriv);
    Output.Vect1DArrayStrings.push_back(string("CellConvDeriv"));

    Output.Vect1DArrays.push_back(CellTiltDeriv);
    Output.Vect1DArrayStrings.push_back(string("CellTiltDeriv"));

    Output.Vect1DArrays.push_back(CellViscDeriv);
    Output.Vect1DArrayStrings.push_back(string("CellViscDeriv"));

    Output.Vect1DArrays.push_back(CellArtDeriv);
    Output.Vect1DArrayStrings.push_back(string("CellArtDeriv"));

    Output.Vect2DArrays.push_back(TransVars);
    Output.Vect2DArrayStrings.push_back(string("TransVars"));


    Array < Array < REAL > > DataOut = UTIL::zeros(FVMCell::AllCells.size(), 21 + (SYSTEM::NumTransVars * 3));
    Mins = Vect3(1e32, 1e32, 1e32);
    Maxs = Vect3(-1e32, -1e32, -1e32);
    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {

        Mins = min(FVMCell::AllCells[i]->Position, Mins);
        Maxs = max(FVMCell::AllCells[i]->Position, Maxs);

    }
    if (globalSystem->useBodies) {
        Output.Double2DArrays.push_back(BODY::CpHistoryAll);
        Output.Double2DArrayStrings.push_back(string("CpHistoryAll"));

        Output.Double2DArrays.push_back(BODY::CpHistoryAllD);
        Output.Double2DArrayStrings.push_back(string("CpHistoryAllD"));

        Output.Double1DArrays.push_back(BODY::SubTIMES);
        Output.Double1DArrayStrings.push_back(string("SubTimes"));


    };
    globalIO->writeMATLABOutputStruct(Output, string("RunData"));

    if (globalSystem->useBodies)
        BODY::SubTIMES.clear();

#ifdef TIME_STEPS
    unsigned long int t2 = ticks();
    stringstream tmp;
    tmp << "Writedata()              : " << double(t2 - t1) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void SYSTEM::WriteBodies() {
    //   
}

/**************************************************************/
void SYSTEM::WritePanelVels() {
    //  
}
