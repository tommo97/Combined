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



#include <deque>

#include "system.hpp"
#include "types.hpp"
#ifndef USEGSL
#include "pgesv.hpp"
#endif

/**************************************************************/
SYSTEM::~SYSTEM() {


#ifdef USEGSL
    gsl_matrix_free(globalA);
    gsl_matrix_free(globalB);
    gsl_permutation_free(globalP);
    gsl_vector_free(globalMu);
    gsl_vector_free(globalSigma);
    gsl_vector_free(globalRHS);
#else
    for (int i = 0; i < NumBodyPanels; ++i) {
        delete [] A[i];
        delete [] B[i];
    }
    delete [] A;
    delete [] B;
    delete rhs;
    delete mu;
    delete sigma;
    delete ipiv;
#endif


    for (int i = 0; i < BodyPoints.size(); ++i) delete BodyPoints[i];
    for (int i = 0; i < Bodies.size(); ++i) delete Bodies[i];
    //    for (int i = 0; i < AllBodyPanels.size(); ++i) delete AllBodyPanels[i];

    delete globalIO;
    delete globalOctree;
    delete globalTimeStepper;

}

/**************************************************************/
SYSTEM::SYSTEM(int NT) {
    OS = getenv("OSTYPE");
    LiftingLineMode = false;
    SysDumpInterval = 0; //  Used to dump the result of system calls
    uinf = vinf = winf = 0;
    MaxP = 3;
    Temp = 288.15; //  Kelvin
    Rho = 1027; //  Kg/m3
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

#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "Platform identified itself as: " << getenv("OSTYPE");
    if (WRITE_TO_SCREEN) cout << " " << VERS << "; Number of threads: " << NumThreads << endl;
#endif

    globalSystem = this;

}

/**************************************************************/
void SYSTEM::Initialise() {
    uinf *= GambitScale;
    vinf *= GambitScale;
    winf *= GambitScale;
    Vinf = Vect3(uinf, vinf, winf);
    Rho /= GambitScale * GambitScale*GambitScale; //  Density falls as volume increases
    Mu /= GambitScale; //   Dynamic Viscocity kg/ms (should this be times scale?)
    Nu = Mu / Rho;
    Array <Vect3> OMEGAS;
    OMEGAS.push_back(Vect3(0, 0, 0));

    globalIO = new IO();
    globalTimeStepper = new TIME_STEPPER();


    ReadNeuGetBodies();

    NumTransVars = NumBodies;

    globalOctree = new OCTREE();


    for (int i = 0; i < NumBodies; ++i)
        Bodies[i]->InitNascentWake(dtInit / NumSubSteps);

    globalTimeStepper->time_step();


    cout << Rho << " " << Mu << " " << Nu << endl;
}

/**************************************************************/
void SYSTEM::TimeStep() {

    while (globalTimeStepper->t < globalTimeStepper->max_t)
        globalTimeStepper->time_loop();

    if (WRITE_TO_SCREEN) cout << "Finished at sim time: " << globalTimeStepper->t << endl;
}

/**************************************************************/
void SYSTEM::SetupGlobalInfluenceMatrices() {
    NumBodyPanels = (int) AllBodyPanels.size();

#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "Allocating space..." << NumBodyPanels * NumBodyPanels * sizeof (REAL) << "... " << NumBodyPanels << endl;
#endif
#ifdef USEGSL
    globalA = gsl_matrix_alloc(NumBodyPanels, NumBodyPanels); // doublet matrix
    globalB = gsl_matrix_alloc(NumBodyPanels, NumBodyPanels); // source matrix
    globalSigma = gsl_vector_alloc(NumBodyPanels);
    globalMu = gsl_vector_alloc(NumBodyPanels);
    globalRHS = gsl_vector_alloc(NumBodyPanels);
#else
    A = new REAL*[NumBodyPanels];
    B = new REAL*[NumBodyPanels];
    rhs = new REAL[NumBodyPanels];
    mu = new REAL[NumBodyPanels];
    sigma = new REAL[NumBodyPanels];
    ipiv = new int[NumBodyPanels];
#endif
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "done" << endl;
#endif


    for (int i = 0; i < NumBodyPanels; ++i) {

#ifdef USEGSL
        gsl_vector_set(globalRHS, i, 0.0);
        gsl_vector_set(globalSigma, i, 0.0);
        gsl_vector_set(globalMu, i, 0.0);
#else
        mu[i] = rhs[i] = sigma[i] = 0.0;
        A[i] = new REAL[NumBodyPanels];
        B[i] = new REAL[NumBodyPanels];
#endif


    }


    UpdateGlobalInfluenceMatrices();


    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
#ifdef USEGSL
        AllBodyPanels[i]->sigma = gsl_vector_ptr(globalSigma, i);
        AllBodyPanels[i]->mu = gsl_vector_ptr(globalMu, i);
#else
        AllBodyPanels[i]->sigma = &sigma[i];
        AllBodyPanels[i]->mu = &mu[i];
#endif
    }



}

/**************************************************************/
void SYSTEM::UpdateGlobalInfluenceMatrices() {
    long int t = ticks();
#ifndef use_NCURSES
    //    if (WRITE_TO_SCREEN) cout << "Calculating influence coefficients..." << endl;
#endif
    //  Temporarily change globalSystem->Del2

    REAL hold = globalSystem->Del2;
    globalSystem->Del2 = 0.;
    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
        for (int j = 0; j < (int) AllBodyPanels.size(); ++j) {
            REAL a = 0, b = 0;
            if (!LiftingLineMode) {
                REAL PHI[2] = {0.};
                AllBodyPanels[j]->SourceDoubletPotential(AllBodyPanels[i]->CollocationPoint, 1, 1, PHI);

                a = PHI[0];
                b = PHI[1];
                if (AllBodyPanels[j]->isBound) {
                    PHI[0] = 0.0;
                    PHI[1] = 0.0;
                    AllBodyPanels[j]->Wake->SourceDoubletPotential(AllBodyPanels[i]->CollocationPoint, 1, 0, PHI);
                    a += (AllBodyPanels[j]->isTop) ? PHI[0] : - PHI[0];
                }
            } else {
                Vect3 V = AllBodyPanels[j]->BodyPanelVelocity(AllBodyPanels[i]->CollocationPoint->vP, 1.0);
                a = V.Dot(AllBodyPanels[i]->TRANS[2]);
                b = 0.0;
                if (i == j)
                    b = 1.0;

            }

#ifdef USEGSL
            gsl_matrix_set(globalA, i, j, a);
            gsl_matrix_set(globalB, i, j, b);
#else
            A[i][j] = a;
            B[i][j] = b;
#endif
        }
    }
#ifndef use_NCURSES
    //    if (WRITE_TO_SCREEN) cout << " time elapsed: " << ticks() - t << endl;
#endif
    t = ticks();

#ifdef USEGSL
#ifndef use_NCURSES
    //    if (WRITE_TO_SCREEN) cout << "Calculating permutation..." << endl;
#endif
    int s;
    globalP = gsl_permutation_alloc(NumBodyPanels);
    gsl_linalg_LU_decomp(globalA, globalP, &s);
#ifndef use_NCURSES
    //    if (WRITE_TO_SCREEN) cout << " time elapsed: " << ticks() - t << endl;
#endif
#endif


    globalSystem->Del2 = hold;
}

/**************************************************************/
void SYSTEM::GetPressures(REAL dt) {
    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
        AllBodyPanels[i]->GetCollocationPoint();
        Vect3 Pos = AllBodyPanels[i]->CollocationPoint->vP;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = AllBodyPanels[i]->Owner->EulerRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = AllBodyPanels[i]->Owner->CG.vV + Vrot;
        // 	Include freestream
        Vect3 V = Vinf - Vkin; // + AllBodyPanels[i]->CollocationPoint->vVfmm;
        for (int J = 0; J < (int) NumBodies; ++J)
            V += Bodies[J]->GetVel(Pos);


        REAL Vref2 = V.Dot(V);



        Vect3 Vlocal = VectMultMatrix(AllBodyPanels[i]->TRANS, V);
        //        if (WRITE_TO_SCREEN) cout << V.Dot(AllBodyPanels[i]->TRANS[2]) << " " << Vlocal << endl;

        if (Pos.y == 0.) {
            // if (WRITE_TO_SCREEN) cout << AllBodyPanels[i]->CollocationPoint->vP << " " << V << " " << AllBodyPanels[i]->TRANS[2] << endl;

        }

        REAL dmu_dt = (*(AllBodyPanels[i]->mu) - (AllBodyPanels[i]->mu_prev)) / dt; //	Get dphi_by_dt

        REAL CP = 1.0 - (V.Dot(V) / Vref2) - (2.0 * dmu_dt / Vref2);


        REAL mult = -CP * 0.5 * _RHO * Vref2 * AllBodyPanels[i]->Area;

        Vect3 deltaF = mult * AllBodyPanels[i]->TRANS[2];
        //        AllBodyPanels[i]->dFORCE_HISTORY.push_back(deltaF);

        AllBodyPanels[i]->Owner->vFORCE += deltaF;
        AllBodyPanels[i]->Owner->vTORQUE += deltaF.Cross(Pos - AllBodyPanels[i]->Owner->CG.vP);
    }
    for (int I = 0; I < NumBodies; ++I) {
        Bodies[I]->ForceRecords.push_back(Bodies[I]->vFORCE);
        Bodies[I]->TorqueRecords.push_back(Bodies[I]->vTORQUE);
        if (WRITE_TO_SCREEN) cout << Bodies[I]->vFORCE << " ";
        if (WRITE_TO_SCREEN) cout << Bodies[I]->vTORQUE << endl;
    }
}

/**************************************************************/
void SYSTEM::GetGlobalRHS() {
#ifdef _OPENMP
    //#pragma omp parallel for
#endif
    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {

        AllBodyPanels[i]->GetCollocationPoint();
        Vect3 Pos = AllBodyPanels[i]->CollocationPoint->vO;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = AllBodyPanels[i]->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = AllBodyPanels[i]->Owner->CG.vV + Vrot;
        // 	Include freestream and FMM wake interpolation
        Vect3 V = Vinf - Vkin + AllBodyPanels[i]->Vfmm; //AllBodyPanels[i]->VelInterp[SubStep];       //  Substep counting starts at 1
        //  Iterate over all wake panels
        for (int J = 0; J < (int) NumBodies; ++J)
            V += Bodies[J]->GetWakeVel(Pos);

        AllBodyPanels[i]->CollocationPoint->vV = V;

        REAL Vn = V.Dot(AllBodyPanels[i]->TRANS[2]);

        REAL Vt1 = V.Dot(AllBodyPanels[i]->TRANS[0]), Vt2 = V.Dot(AllBodyPanels[i]->TRANS[1]);

        AllBodyPanels[i]->alpha = atan(Vn / Vt1)*180 / pi;




        *(AllBodyPanels[i]->sigma) = -Vn;
#ifdef USEGSL
        gsl_vector_set(globalSigma, i, -Vn);
#else
        sigma[i] = -Vn;
        rhs[i] = 0.0;
#endif
    }

#ifndef USEGSL
    for (int i = 0; i < NumBodyPanels; ++i) for (int j = 0; j < NumBodyPanels; ++j) rhs[i] += B[i][j] * sigma[j];
#else
    REAL alpha = 1.0, beta = 0.0;
#ifdef DOUBLE_PRECISION
    if (gsl_blas_dgemv(CblasNoTrans, alpha, globalB, globalSigma, beta, globalRHS))
#else
    if (gsl_blas_sgemv(CblasNoTrans, alpha, globalB, globalSigma, beta, globalRHS))
#endif
        if (WRITE_TO_SCREEN) cout << "GSL Matrix/Vector product error." << endl;
#endif
}

/**************************************************************/
void SYSTEM::LinAlg() {

    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
        AllBodyPanels[i]->mu_prev = AllBodyPanels[i]->gamma_prev = *(AllBodyPanels[i]->mu);
    }

#ifndef use_NCURSES
#ifdef DEBUG
    long int t = ticks();
#endif
#endif

#ifdef USEGSL
    gsl_linalg_LU_solve(globalA, globalP, globalRHS, globalMu);
#else
    pgesv(A, rhs, NumBodyPanels, mu);
#endif
#ifndef use_NCURSES
#ifdef DEBUG
    if (WRITE_TO_SCREEN) cout << "Linear algebra took: " << ticks() - t << endl;
#endif
#endif
    for (int i = 0; i < (int) AllBodyPanels.size(); ++i)
        AllBodyPanels[i]->gamma = *(AllBodyPanels[i]->mu);


    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < (int) Bodies[i]->ProtoWake.size(); ++j) {
            Bodies[i]->ProtoWake[j]->gamma_prev = Bodies[i]->ProtoWake[j]->gamma;
            if (!LiftingLineMode)
                Bodies[i]->ProtoWake[j]->gamma = Bodies[i]->ProtoWake[j]->Shedder->OtherBoundarySurface->gamma -
                    Bodies[i]->ProtoWake[j]->Shedder->gamma;
            else
                Bodies[i]->ProtoWake[j]->gamma = Bodies[i]->ProtoWake[j]->Shedder->gamma;
        }

    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < (int) Bodies[i]->ProtoWake.size(); ++j)
            Bodies[i]->ProtoWake[j]->WakeNeighbSet();
}

/**************************************************************/
void SYSTEM::BodySubStep(REAL delta_t, int n_steps) {

    unsigned long int t0 = ticks();
    REAL dt = delta_t / n_steps;
    //    cout << "dt " << dt << " " << n_steps << endl;


    //    for (int i = 0; i < (int) AllBodyPanels.size(); ++i)
    //        AllBodyPanels[i]->VelInterp = globalLinspace(AllBodyPanels[i]->CollocationPoint->vVfmm[0],AllBodyPanels[i]->CollocationPoint->vVfmm[1],n_steps);




    for (int SubStep = 1; SubStep <= n_steps; ++SubStep) {

        globalTimeStepper->substep_time = SubStep*dt;

        MoveBodies(dt, true);


        //	Sort wake convection
        Vect3 Vwake, Vbody;
        for (int i = 0; i < NumBodies; ++i) {

            for (int j = 0; j < Bodies[i]->WakePoints.size(); ++j) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(i,j,dt) private(Vwake, Vbody)
#endif
                for (int k = 0; k < Bodies[i]->WakePoints[j].size(); ++k) {
//                    for (int l = 0; l < NumBodies; ++l)
//                        Vbody = Bodies[l]->GetVel(Bodies[i]->WakePoints[j][k]->vP);

                    Bodies[i]->WakePoints[j][k]->vV = (Vinf + Vbody);
//                    Bodies[i]->WakePoints[j][k]->vV += globalOctree->TreeVel(Bodies[i]->WakePoints[j][k]->vP);
                }
            }
        }



        for (int i = 0; i < NumBodies; ++i)
            for (int j = 0; j < Bodies[i]->WakePoints.size(); ++j)
                for (int k = 0; k < Bodies[i]->WakePoints[j].size(); ++k)
                    Bodies[i]->WakePoints[j][k]->vP += Bodies[i]->WakePoints[j][k]->vV * dt;


        GetGlobalRHS();
        LinAlg();



        for (int j = 0; j < NumBodies; ++j)
            Bodies[j]->DissolveWake(dt);




        //        if (SubStep % 1 == 0) globalIO->write_m();
#ifndef USE_NCURSES
        // #ifdef DEBUG
        //        if (WRITE_TO_SCREEN) {
        //            cout << "Substep " << SubStep << " of " << n_steps;
        //            cout << " Phi " << RAD2DEG(Bodies[0]->EulerAngles.x) << " ";
        //            cout << (ticks() - t0) / 1000 << " dt " << dt << endl;
        //        }
        // #endif
#endif




    }

}

/**************************************************************/
void SYSTEM::ReadNeuGetBodies() {
    Vect3 ORIGIN[2], ATTITUDE, VELOCITY[2], RATES[2];

    ORIGIN[0] = Vect3(0, 0, 0);
    ORIGIN[1] = Vect3(2.5, 0, 0);
    RATES[0] = Vect3(-7.5,0.,0.);
    RATES[1] = Vect3(7.5, 0, 0);
    VELOCITY[0] = VELOCITY[1] = Vect3(-10., 0., 0.);
    // These are the BODY rates about BODY axis

    Array <Vect3> X;
    ARRAY2(int) PNLS, GROUPS, BCS;
    Array <string> NAMES;
    globalIO->read_neu(NeuFile.data(), X, PNLS, GROUPS, BCS, NAMES);


    for (int i = 0; i < (int) X.size(); ++i) {
        POINT *P = new POINT(GambitScale * X[i]);
        BodyPoints.push_back(P);
    }

    for (int i = 0; i < (int) PNLS.size(); ++i) {
        POINT *P1, *P2, *P3, *P4;
        P1 = BodyPoints[PNLS[i][0]];
        P2 = BodyPoints[PNLS[i][1]];
        P3 = BodyPoints[PNLS[i][2]];
        P4 = BodyPoints[PNLS[i][3]];

        P1->vO = P1->vP;
        P2->vO = P2->vP;
        P3->vO = P3->vP;
        P4->vO = P4->vP;

        PANEL *temp = new PANEL(P1, P2, P3, P4);

        AllBodyPanels.push_back(temp);

    }


    // Get BCs sorted: here we want to create the proto-wake panels
    for (int i = 0; i < (int) BCS.size(); ++i) {
        int j = BCS[i][0];
        AllBodyPanels[j]->isBound = true;

        if (BCS[i][1] == 0) {
            AllBodyPanels[j]->edgeX1 = AllBodyPanels[j]->C1;
            AllBodyPanels[j]->edgeX2 = AllBodyPanels[j]->C2;
        }
        if (BCS[i][1] == 1) {
            AllBodyPanels[j]->edgeX1 = AllBodyPanels[j]->C2;
            AllBodyPanels[j]->edgeX2 = AllBodyPanels[j]->C3;
        }
        if (BCS[i][1] == 2) {
            AllBodyPanels[j]->edgeX1 = AllBodyPanels[j]->C3;
            AllBodyPanels[j]->edgeX2 = AllBodyPanels[j]->C4;
        }
        if (BCS[i][1] == 3) {
            AllBodyPanels[j]->edgeX1 = AllBodyPanels[j]->C4;
            AllBodyPanels[j]->edgeX2 = AllBodyPanels[j]->C1;
        }
    }

    NumBodies = (int) GROUPS.size();


    for (int i = 0; i < NumBodies; ++i) {

        //        ORIGIN *= GambitScale;
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << NAMES[i] << " " << RATES[i] << " " << ORIGIN[i] << endl;
        if (WRITE_TO_SCREEN) cout << "making body " << endl;
#endif

        BODY *Btemp = new BODY(ORIGIN[i], ATTITUDE, VELOCITY[i], RATES[i], NAMES[i], this);
        Btemp->BodyID = i;
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "done. Getting panels:" << endl;
#endif
        Btemp->GetPanels(GROUPS[i], AllBodyPanels);
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "done" << endl;
#endif
        Bodies.push_back(Btemp);
    }

    for (int i = 0; i < BodyPoints.size(); ++i)
        BodyPoints[i]->vO -= BodyPoints[i]->Owner->CG.vP;
    
    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
        AllBodyPanels[i]->GetNewGlobalPosition();
        AllBodyPanels[i]->GetCollocationPoint();
        AllBodyPanels[i]->Owner->Rmax = max(AllBodyPanels[i]->Owner->Rmax, AllBodyPanels[i]->C1->vO.Mag());
        AllBodyPanels[i]->Owner->Rmax = max(AllBodyPanels[i]->Owner->Rmax, AllBodyPanels[i]->C2->vO.Mag());
        AllBodyPanels[i]->Owner->Rmax = max(AllBodyPanels[i]->Owner->Rmax, AllBodyPanels[i]->C3->vO.Mag());
        AllBodyPanels[i]->Owner->Rmax = max(AllBodyPanels[i]->Owner->Rmax, AllBodyPanels[i]->C4->vO.Mag());
    }



        for (int i = 0; i < NumBodies; ++i) {
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << Bodies[i]->Name << " " << Bodies[i]->NumFaces;
        if (WRITE_TO_SCREEN) cout << " " << Bodies[i]->BoundaryFaces.size() << endl;
#endif
    }

    for (int i = 0; i < (int) NumBodies; ++i)
        for (int j = 0; j < (int) Bodies[i]->Faces.size(); ++j) {
            Bodies[i]->Faces[j]->Owner = Bodies[i];
        }




    for (int i = 0; i < (int) AllBodyPanels.size(); ++i)
        for (int j = 0; j < (int) AllBodyPanels.size(); ++j)
            AllBodyPanels[i]->CheckNeighb(AllBodyPanels[j]);

    for (int i = 0; i < (int) AllBodyPanels.size(); ++i) {
        if (AllBodyPanels[i]->Neighb.L) {
            REAL N1N2 = AllBodyPanels[i]->TRANS[2].Dot(AllBodyPanels[i]->Neighb.L->TRANS[2]);
            REAL costheta = N1N2 / (AllBodyPanels[i]->TRANS[2].Mag() * AllBodyPanels[i]->Neighb.L->TRANS[2].Mag());
            costheta = min(fabs(costheta), 1.0);
            //                cout << "L " << costheta << " " << acos(costheta) << " " << RAD2DEG(acos(costheta)) << endl;
        }
        if (AllBodyPanels[i]->Neighb.R) {
            REAL N1N2 = AllBodyPanels[i]->TRANS[2].Dot(AllBodyPanels[i]->Neighb.R->TRANS[2]);
            REAL costheta = N1N2 / (AllBodyPanels[i]->TRANS[2].Mag() * AllBodyPanels[i]->Neighb.R->TRANS[2].Mag());
            costheta = min(fabs(costheta), 1.0);
            //                cout << "R " << costheta << " " << acos(costheta) << " " << RAD2DEG(acos(costheta)) << endl;
        }
        if (AllBodyPanels[i]->Neighb.T) {
            REAL N1N2 = AllBodyPanels[i]->TRANS[2].Dot(AllBodyPanels[i]->Neighb.T->TRANS[2]);
            REAL costheta = N1N2 / (AllBodyPanels[i]->TRANS[2].Mag() * AllBodyPanels[i]->Neighb.T->TRANS[2].Mag());
            costheta = min(fabs(costheta), 1.0);
            //                cout << "T " << costheta << " " << acos(costheta) << " " << RAD2DEG(acos(costheta)) << endl;
        }
        if (AllBodyPanels[i]->Neighb.B) {
            REAL N1N2 = AllBodyPanels[i]->TRANS[2].Dot(AllBodyPanels[i]->Neighb.B->TRANS[2]);
            REAL costheta = N1N2 / (AllBodyPanels[i]->TRANS[2].Mag() * AllBodyPanels[i]->Neighb.B->TRANS[2].Mag());
            costheta = min(fabs(costheta), 1.0);
            //                cout << "B " << costheta << " " << acos(costheta) << " " << RAD2DEG(acos(costheta)) << endl;
        }

    }
}

/**************************************************************/
void SYSTEM::PrintInfluenceMatrices() {
    //if (WRITE_TO_SCREEN) cout.setf(ios::fixed, ios::floatfield);
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "A = zeros(" << NumBodyPanels << "); B = A;" << endl;
    if (WRITE_TO_SCREEN) cout << "%Doublet (A) Matrix: " << endl;
#endif
    REAL a, b;
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "A_global = [";
#endif
    for (int i = 0; i < NumBodyPanels; i++) {
        for (int j = 0; j < NumBodyPanels; j++) {
#ifdef USEGSL
            a = gsl_matrix_get(globalA, i, j);
#else
            a = A[i][j];
#endif
#ifndef use_NCURSES
            if (WRITE_TO_SCREEN) cout << a << ", ";
#endif
        }
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "; " << endl;
#endif
    }
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "]; " << endl;


    if (WRITE_TO_SCREEN) cout << "%Source (B) Matrix: " << endl;
    if (WRITE_TO_SCREEN) cout << "B_global = [";
#endif
    for (int i = 0; i < NumBodyPanels; i++) {
        for (int j = 0; j < NumBodyPanels; j++) {
#ifdef USEGSL
            b = gsl_matrix_get(globalB, i, j);
#else
            b = B[i][j];
#endif
#ifndef use_NCURSES
            if (WRITE_TO_SCREEN) cout << b << ", ";
#endif
        }
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "; " << endl;
#endif
    }
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "]; " << endl;

    if (WRITE_TO_SCREEN) cout << "%RHS (B.Sigma) Array: " << endl;
    if (WRITE_TO_SCREEN) cout << "rhs_global = [";
#endif
    for (int i = 0; i < NumBodyPanels; i++) {
#ifdef USEGSL
        b = gsl_vector_get(globalRHS, i);
#else
        b = sigma[i];
#endif
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << b << "; ";
#endif
    }
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "]; " << endl;

    if (WRITE_TO_SCREEN) cout << "%Mu Vector: " << endl;
    if (WRITE_TO_SCREEN) cout << "moo_global = [";
#endif
    for (int i = 0; i < NumBodyPanels; i++) {
#ifdef USEGSL
        b = gsl_vector_get(globalMu, i);
#else
        b = mu[i];
#endif
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << b << "; ";
#endif
    }
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "]; " << endl;
#endif
}

/**************************************************************/
void SYSTEM::PrintBodiesVels() {
#ifndef use_NCURSES
    for (int i = 0; i < NumBodies; ++i) {
        // 	Bodies[i]->PrintVels();
    }
#endif
}

/**************************************************************/
void SYSTEM::PrintBodiesAndWakes() {
#ifndef use_NCURSES
    for (int i = 0; i < NumBodies; ++i) {
        Bodies[i]->PrintSurface();
        Bodies[i]->PrintWake();
    }
    if (WRITE_TO_SCREEN) cout << "lighting phong" << endl << "set(gcf,'Renderer','zbuffer')" << endl;
#endif
}

/**************************************************************/
void SYSTEM::WriteBodiesAndWakes(ostream& out_stream) {
    for (int i = 0; i < NumBodies; ++i) {
        Bodies[i]->WriteSurface(out_stream);
        Bodies[i]->WriteWake(out_stream);
    }
    if (WRITE_TO_FILE) out_stream << "lighting phong" << endl << "set(gcf,'Renderer','zbuffer')" << endl;
}

/**************************************************************/
void SYSTEM::PutWakesInTree() {
    for (int i = 0; i < NumBodies; ++i) {
        for (int j = 0; j < globalSystem->NumSubSteps; ++j) {
            for (int k = 0; k < Bodies[i]->WakePoints.back().size(); ++k) {
                OctreeCapsule C(Bodies[i]->WakePoints.back()[k]->vP, Bodies[i]->WakePoints.back()[k]->vO, true);
                C.AssociatedBody = i;
                globalOctree->Root->EvalCapsule(C);
                delete Bodies[i]->WakePoints.back()[k];
            }
            Bodies[i]->WakePoints.pop_back();
        }
    }
}

/**************************************************************/
void SYSTEM::GetFaceVels() {
    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < globalOctree->AllCells.size(); ++j)
#ifdef COLLAPSE_TO_FACES
            for (int k = 0; k < 6; ++k)
                if ((k == 0) || (k == 2) || (k == 4) || !globalOctree->AllCells[j]->Neighb[k])
                    globalOctree->AllCells[j]->FaceVels[k] +=
                        Bodies[i]->GetVel(globalOctree->AllCells[j]->Position + Node::NeighbOffset[k]);
#else
            globalOctree->AllCells[j]->Velocity += Bodies[i]->GetVel(globalOctree->AllCells[j]->Position);
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();

}

/**************************************************************/
void SYSTEM::GetPanelFMMVelocities() {
    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < Bodies[i]->Faces.size(); ++j)
            Bodies[i]->Faces[j]->Vfmm =
                globalOctree->TreeVel(Bodies[i]->Faces[j]->CollocationPoint->vP);

}

/**************************************************************/
void SYSTEM::MoveBodies(REAL dt, bool update) {
    for (int j = 0; j < NumBodies; ++j)
        Bodies[j]->MoveBody(dt);

    if (update && NumBodies > 1)
        UpdateGlobalInfluenceMatrices();
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
                outstream << i + 1 << "\t" << j + 1 << "\t" << k + 1 << "\t" << X << "\t" << V + Vinf << endl;
            }

    globalIO->write_file(outstream, string("data"), string("dat"));
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


    globalIO->write_file(outstream, string("f"), string("dat"));
}

/**************************************************************/
void SYSTEM::WritePanelVels() {
    //  Here we want to find the extents of the domain
    Vect3 Max(.5, .5, .5), Min(.5, .5, .5);
    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
        Max = max(Max, globalOctree->AllCells[i]->Position);
        Min = min(Min, globalOctree->AllCells[i]->Position);
    }


    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < Bodies[i]->WakePoints.size(); ++j)
            for (int k = 0; k < Bodies[i]->WakePoints[j].size(); ++k) {
                Max = max(Max, Bodies[i]->WakePoints[j][k]->vP);
                Min = min(Min, Bodies[i]->WakePoints[j][k]->vP);
            }


    Max = floor(Max) + .5;
    Min = floor(Min) + .5;
    Vect3 Sz = Max - Min;
    cout << Max << "\t " << Min << "\t " << Sz << endl;


    int buffer = 10;
    stringstream outstream;
    //  Make the domain slightly larger and get velocities

    for (int i = 0; i < Sz.x + 2 * buffer; ++i)
        for (int j = 0; j < Sz.y + 2 * buffer; ++j)
            for (int k = 0; k < Sz.z + 2 * buffer; ++k) {
                Vect3 X, V;
                X = Min - buffer + Vect3(i, j, k);
                for (int l = 0; l < NumBodies; ++l)
                    V += Bodies[l]->GetVel(X);

                outstream << i + 1 << "\t" << j + 1 << "\t" << k + 1 << "\t" << X << "\t" << V + Vinf << endl;
            }

    globalIO->write_file(outstream, string("panel_vels"), string("dat"));
}
