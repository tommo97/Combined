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
#define PANEL_MODE

#include "body.hpp"
#include "waves.hpp"

/**************************************************************/
#include "body.hpp"
#include "panel.hpp"
#include "pgesv.hpp"

Array <Vect3> BODY::VortexPositions, BODY::VortexOmegas, BODY::VortexVelocities;
Array <Vect3*> BODY::VortexOrigins;
Array <int> BODY::VortexOwnerID;
Array <string> BODY::NAMES;
Array <REAL> BODY::Times;
REAL BODY::WaitLenghts;
Array <Vect3> BODY::CGS, BODY::VX, BODY::VO;
Array <Vect3> BODY::RATES;
Array <Vect3> BODY::VELOCITY;
Array <Vect3> BODY::ATTITUDE;
Array <REAL> BODY::TimePrev = Array <REAL> (4, 0.0);
REAL BODY::Time = 0;
REAL BODY::RHO = 1;
int BODY::BodyPanelIDCounter = 0, BODY::BodyPointIDCounter = 0, BODY::SubStep = 0;
Array <Vect3> BODY::TorqueHist;
Array <Vect3> BODY::ForceHist;
Array < BODY*> BODY::Bodies;
Array <PANEL*> BODY::AllBodyFaces;
Array < Array < Array < int > > > BODY::Surfaces;
Array < Array < Array < PANEL* > > > BODY::ptSurfaces;
Array < Array <REAL> > BODY::A;
Array < Array <REAL> > BODY::B;
Array < Array <REAL> > BODY::C;
Array < Array <REAL> > BODY::D;
Array < Array <REAL> > BODY::LocalChordRadius;
Array <REAL> BODY::AlphaHistory;
Array <REAL> BODY::AlphaDotHistory;
Array < Array <REAL> > BODY::CpHistory, BODY::CpHistoryAll;
Array < Array <REAL> > BODY::CpHistoryD, BODY::CpHistoryAllD;
Array <REAL> BODY::SubTIMES;
Array <REAL> BODY::RHS;
Array <REAL> BODY::Mu;
Array <REAL> BODY::Sigma;
Array <REAL> BODY::DPhi;
Array <REAL> BODY::DCp;
Array <REAL> BODY::DphiDCp;
Array <REAL> BODY::LiftHist;
Array <REAL> BODY::dDeltaCp_dDeltaPhi;
Array <PANEL*> BODY::AllProtoWakes;
Array < Array <Array <int> > > BODY::PtIDS;
Array < Array <Array <int> > > BODY::TipInboardUSIDS;
Array < Array <Array <int> > > BODY::TipInboardLSIDS;
Array < Array <Array <int> > > BODY::TipOutboardUSIDS;
Array < Array <Array <int> > > BODY::TipOutboardLSIDS;
Array < Array < Array < int > > > BODY::InnerTipUSPanelIDS, BODY::OuterTipUSPanelIDS;
Array < Array < Array < int > > > BODY::InnerTipLSPanelIDS, BODY::OuterTipLSPanelIDS;
Array <Vect3> BODY::PointsAsRead;
Array <Array < int > > BODY::PanelsAsRead;
Array <Vect3> BODY::AllBodyPoints;
bool BODY::IPKC = false;
int BODY::NumFaces = 0;
int BODY::NumBodies = 0;
REAL BODY::Radius, BODY::TSR;
bool BODY::LiftingLineMode;

/**************************************************************/
/*  Body functions      */

/**************************************************************/
void BODY::MakeWake() {

    if (WakePanels[0].size() == 0)
        cout << "!!!!!!!!!!!!! No Wake Panels!" << endl;
    Vect3 C1, C2, C3, C4;
    for (int i = 0; i < ProtoWakes.size(); ++i) {
        Array <PANEL*> tmp;
        for (int j = 0; j < ProtoWakes[i].size(); ++j) {
            if (WakePanels[i].size() == 0) {
                C1 = ProtoWakeLastC2[i][j];
                C2 = ProtoWakes[i][j].C2;
                C3 = ProtoWakes[i][j].C3;
                C4 = ProtoWakeLastC3[i][j];
            } else {
                C1 = WakePanels[i].back()[j]->C2;
                C2 = ProtoWakes[i][j].C2;
                C3 = ProtoWakes[i][j].C3;
                C4 = WakePanels[i].back()[j]->C3;
            }
            tmp.push_back(new PANEL(C1, C2, C3, C4));
            tmp.back()->Gamma = ProtoWakes[i][j].Gamma;
            tmp.back()->SheddingSurface = ProtoWakes[i][j].SheddingSurface;
        }

        for (int j = 0; j < tmp.size(); ++j)
            for (int k = 0; k < tmp.size(); ++k)
                tmp[j]->CheckNeighb(tmp[k]);

        WakePanels[i].push_back(tmp);


    }

    int n = 2;  // also try 2...
    Array <Vect3> Pts, tPts;
    Array <Vect3> Oms, tOms;
    Array <Vect3*> Origins, tOrigins;
    for (int i = 0; i < WakePanels.size(); ++i)
        if (WakePanels[i].size() > 15) {
            for (int j = 0; j < WakePanels[i][0].size(); ++j) {
                PANEL *tmp = (WakePanels[i][0][j]);
                for (int k = 0; k < n; ++k) {
                    tPts.push_back(tmp->C1 + (k * 1.0 / (n - 1)) * (tmp->C2 - tmp->C1));
                    tPts.push_back(tmp->C2 + (k * 1.0 / (n - 1)) * (tmp->C3 - tmp->C2));
                    tPts.push_back(tmp->C3 + (k * 1.0 / (n - 1)) * (tmp->C4 - tmp->C3));
                    tPts.push_back(tmp->C4 + (k * 1.0 / (n - 1)) * (tmp->C1 - tmp->C4));
                    tOms.push_back((1.0 / n) * tmp->Gamma * (tmp->C2 - tmp->C1));
                    tOms.push_back((1.0 / n) * tmp->Gamma * (tmp->C3 - tmp->C2));
                    tOms.push_back((1.0 / n) * tmp->Gamma * (tmp->C4 - tmp->C3));
                    tOms.push_back((1.0 / n) * tmp->Gamma * (tmp->C1 - tmp->C4));
                    tOrigins.push_back(&(tmp->SheddingSurface->CollocationPoint));
                    tOrigins.push_back(&(tmp->SheddingSurface->CollocationPoint));
                    tOrigins.push_back(&(tmp->SheddingSurface->CollocationPoint));
                    tOrigins.push_back(&(tmp->SheddingSurface->CollocationPoint));
                }

                Pts.push_back(tPts);
                Oms.push_back(tOms);
                Origins.push_back(tOrigins);
                tPts.clear();
                tOms.clear();
                tOrigins.clear();
                delete tmp;

            }
            WakePanels[i].pop_front();
        }




    //      Find unique vortices - this is a total bodge, but I cannot be bothered writing proper code

    Array <Vect3> Posns, Vorts;
    Array <Vect3*> Orgns;
    bool test;

    for (int i = 0; i < Pts.size(); ++i) {
        test = false;
        for (int j = 0; j < Posns.size(); ++j) {
            if (Posns[j] == Pts[i]) {
                Vorts[j] += Oms[i];
                test = true;
                break;
            }
        }
        if (test == false) {
            Posns.push_back(Pts[i]);
            Vorts.push_back(Oms[i]);
            Orgns.push_back(Origins[i]);
        }
    }


    BODY::VortexPositions.push_back(Posns); // = Array <Vect3 > (BODY::VortexPositions.size() + Pts.size());
    BODY::VortexOmegas.push_back(Vorts); // = Array <Vect3 > (BODY::VortexOmegas.size() + Pts.size());
    BODY::VortexOwnerID.push_back(Array <int> (Pts.size(), this->ID)); // = Array <int> (BODY::VortexOwnerID.size() + Pts.size());
    BODY::VortexOrigins.push_back(Orgns);


    BODY::VortexVelocities = Array <Vect3 > (BODY::VortexOwnerID.size(), Vect3(0.0));

}

/**************************************************************/
void BODY::GetLinearRHS() {




#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = BODY::AllBodyFaces[i]->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = BODY::AllBodyFaces[i]->Owner->Velocity + Vrot;
        // 	Include freestream and FMM wake interpolation
        Vect3 V = globalSystem->unscaledVinf - Vkin; // + BODY::AllBodyFaces[i]->Vfmm; //BODY::AllBodyFaces[i]->VelInterp[SubStep];       //  Substep counting starts at 1

        Vect3 VWake = BODY::AllBodyFaces[i]->Vfmm;

        REAL PhiWake = 0.0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < BODY::VortexPositions.size(); ++j)
            VWake = VWake + UTIL::globalDirectVel(BODY::AllBodyFaces[i]->CollocationPoint - BODY::VortexPositions[j], BODY::VortexOmegas[j]);




        //        Iterate over all wake panels : using the potential formulation on the RHS seems to work better than using the velocity via sigma
        //        for (int J = 0; J < (int) NumBodies; ++J)
        //            VWake = Bodies[J]->GetWakeVel(BODY::AllBodyFaces[i]->CollocationPoint);


        //        for (int J = 0; J < (int) NumBodies; ++J)
        //            PhiWake += Bodies[J]->GetWakePhi(BODY::AllBodyFaces[i]->CollocationPoint);


        BODY::AllBodyFaces[i]->Phi = PhiWake;
        BODY::AllBodyFaces[i]->Vkin = Vkin;
        BODY::AllBodyFaces[i]->VWake = VWake;
        BODY::AllBodyFaces[i]->VCentroid = globalSystem->unscaledVinf - Vkin + VWake;

        V = BODY::AllBodyFaces[i]->VCentroid;

        REAL Vn = BODY::AllBodyFaces[i]->Vn = V.Dot(BODY::AllBodyFaces[i]->TRANS[2]);

        BODY::AllBodyFaces[i]->Sigma = Vn;

        BODY::Sigma[i] = Vn;
        BODY::RHS[i] = -PhiWake;


        //                    REAL PhiE = 0, PhiW = 0, PhiN = 0, PhiS = 0, PhiT = 0, PhiB = 0, dlta = 1e-4;
        //            for (int J = 0; J < (int) NumBodies; ++J)
        //                VWake = Bodies[J]->GetWakeVel(Pos);
        //
        //
        //            for (int J = 0; J < (int) NumBodies; ++J) {
        //                PhiE = Bodies[J]->GetWakePhi(Pos + Vect3(dlta, 0, 0));
        //                PhiW = Bodies[J]->GetWakePhi(Pos - Vect3(dlta, 0, 0));
        //                PhiN = Bodies[J]->GetWakePhi(Pos + Vect3(0, dlta, 0));
        //                PhiS = Bodies[J]->GetWakePhi(Pos - Vect3(0, dlta, 0));
        //                PhiT = Bodies[J]->GetWakePhi(Pos + Vect3(0, 0, dlta));
        //                PhiB = Bodies[J]->GetWakePhi(Pos - Vect3(0, 0, dlta));
        //            }
        //
        //
        //            cout << "Direct  : " << VWake << endl;
        //            cout << "Gradient: " << (PhiE - PhiW) / (2 * dlta);
        //            cout << " " << (PhiN - PhiS) / (2 * dlta);
        //            cout << " " << (PhiT - PhiB) / (2 * dlta) << endl;
    }





    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
            BODY::RHS[i] += BODY::B[i][j] * BODY::Sigma[j];
            //            BODY::RHS[i] -= BODY::C[i][j] * BODY::Mu[j];
        }


    }
}

/**************************************************************/
void BODY::GetNonLinearRHS() {
#ifdef _OPENMP
#pragma omp parallel for  
#endif 
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;

        Vect3 U;

        //  Iterate over all protowake panels
        REAL PhiWake = 0.0;
        for (int j = 0; j < BODY::AllProtoWakes.size(); ++j) {
            U += BODY::AllProtoWakes[j]->VortexPanelVelocity(Pos);
            //            REAL fid = 0;
            //            PANEL::DoubletPotential(BODY::AllProtoWakes[j], Pos, fid, 1, 2);
            //            PhiWake += BODY::AllProtoWakes[j]->Gamma * fid;
        }

        Vect3 VWake = BODY::AllBodyFaces[i]->Vfmm;
        BODY::AllBodyFaces[i]->VCentroid = globalSystem->unscaledVinf + VWake - BODY::AllBodyFaces[i]->Vkin + BODY::AllBodyFaces[i]->VWake + U;

        Vect3 V = BODY::AllBodyFaces[i]->VCentroid;

        REAL Vn = BODY::AllBodyFaces[i]->Vn = V.Dot(BODY::AllBodyFaces[i]->TRANS[2]);

        BODY::AllBodyFaces[i]->Sigma = Vn;

        BODY::Sigma[i] = Vn;
        BODY::RHS[i] = -(BODY::AllBodyFaces[i]->Phi + PhiWake);
    }

    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
        for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
            BODY::RHS[i] += BODY::B[i][j] * BODY::Sigma[j];
        }

}

/**************************************************************/
void BODY::GetPanelVels() {
    //    for (int i = 0; i < Faces.size(); ++i) {
    //        Faces[i].Vd = Vect3(0.0);
    //        for (int j = 0; j < Faces.size(); ++j) {
    //            Faces[i].Vd += -Faces[j].Mu * localVD[i][j];
    //            Faces[i].Vd += Faces[j].Sigma * localVS[i][j];
    //        }
    //    }
}

/**************************************************************/
void BODY::SplitUpLinearAlgebra() {



    /*      This function calculates the new body surface circulation strengths in the case
     *      that it is problematic (read time consuming) to solve a large contiguous
     *      system (A) matrix. The process is:
     *      1)      for all body faces determine wake/kinematic velocity
     *      2)      Use this to calculate each bodies independent circulation (mu) vector
     *      3)      If |mu(k) - mu(k-1)| < tol, return, else
     *      4)      If first time round the loop generate influence matrices for inter-body
     *              interactions
     *      5)      Calculate new local RHSs using body mus and influent matrices 
     */


    //  Calculate velocity at all faces
    Array <PANEL*> srcs;
    for (int I = 0; I < BODY::Bodies.size(); ++I)
        for (int J = 0; J < BODY::Bodies[I]->WakePanels.size(); ++J)
            for (int j = 0; j < BODY::Bodies[I]->WakePanels[J].size(); ++j)
                for (int k = 0; k < BODY::Bodies[I]->WakePanels[J][j].size(); ++k) {
                    PANEL *src = BODY::Bodies[I]->WakePanels[J][j][k];
                    srcs.push_back(src);
                }


#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        PANEL *trg = BODY::AllBodyFaces[i];

        Vect3 Pos = trg->CollocationPoint - trg->Owner->CG;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = trg->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = trg->Owner->Velocity + Vrot;
        // 	Include freestream and FMM wake interpolation
        Vect3 V = globalSystem->unscaledVinf - Vkin; // + trg->Vfmm; //trg->VelInterp[SubStep];       //  Substep counting starts at 1

        Vect3 VWake = 0.0;

        REAL PhiWake = 0.0;

        for (int j = 0; j < BODY::VortexPositions.size(); ++j)
            VWake += 2.0 * UTIL::globalDirectVel(trg->CollocationPoint - BODY::VortexPositions[j], BODY::VortexOmegas[j]);

        /*        for (int j = 0; j < FVMCell::AllCells.size(); ++j)
                    VWake += 2.0*globalDirectVel(trg->CollocationPoint - FVMCell::AllCells[j]->Position,
                        FVMCell::AllCells[j]->Omega);
         */
        VWake += 1.0 * trg->Vfmm;
        
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
//        REAL a = 0.150/2.0, om = 2.0*pi, d = 2.4, g = 9.80665, k = 4.0257;
        REAL t = TIME_STEPPER::SimTime;
//        REAL PhiWave = (a * g / om) * cosh(k * (d + (trg->CollocationPoint.z - 0.8))) * cos(om*t - k*trg->CollocationPoint.x)/cosh(k*d);
//        REAL UWave = -(a*g*k*sin(trg->CollocationPoint.x*k - om*t)*cosh(k*((trg->CollocationPoint.z - 0.8) + d)))/(om*cosh(d*k));
//        REAL VWave = (a*g*k*cos(trg->CollocationPoint.x*k - om*t)*sinh(k*((trg->CollocationPoint.z - 0.8) + d)))/(om*cosh(d*k));
//        Vect3 WaveVel = WaveField::Cnoidal.CnoidalVelocity(trg->CollocationPoint  - Vect3(0.,0.,0.275),t );
//        VWake += WaveVel*SYSTEM::GambitScale;
//        VWake.x += UWave*SYSTEM::GambitScale;
//        VWake.z += VWave*SYSTEM::GambitScale;
        for (int j = 0; j < srcs.size(); ++j){
            srcs[j]->Mu = srcs[j]->Gamma;
            PhiWake -= srcs[j]->GetTriTesselatedDoubletPhi(trg->CollocationPoint);
        }
//            PhiWake += srcs[j]->WakePanelPotential(trg->CollocationPoint);


        trg->Phi = PhiWake + WaveField::Cnoidal.CnoidalPerturbationPotential(trg->CollocationPoint  - Vect3(0.,0.,0.200),t );// + 0.0*PhiWave;// + trg->PhiWakePrev;
        trg->Vkin = Vkin;
        trg->VWake = VWake;
        trg->VCentroid = globalSystem->unscaledVinf - Vkin + VWake;

        V = trg->VCentroid;

        REAL Vn = trg->Vn = V.Dot(trg->TRANS[2]);

        trg->Sigma = Vn;


    }


    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->MuPrev = BODY::AllBodyFaces[i]->GammaPrev = BODY::AllBodyFaces[i]->Mu;
    }



    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        int n = BODY::Bodies[I]->Faces.size();

        BODY::Bodies[I]->localSigma = Array <REAL > (n, 0.0);

        for (int i = 0; i < n; ++i)
            BODY::Bodies[I]->localSigma[i] = BODY::Bodies[I]->Faces[i].Sigma;

        BODY::Bodies[I]->localRHS = Array <REAL > (n, 0.0);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                BODY::Bodies[I]->localRHS[i] += BODY::Bodies[I]->localB[i][j] * BODY::Bodies[I]->localSigma[j];

        for (int i = 0; i < n; ++i)
            BODY::Bodies[I]->localRHS[i] += BODY::Bodies[I]->Faces[i].Phi;

        BODY::Bodies[I]->localMu = Array <REAL > (n, 0.0);


#ifdef USE_MATRIX_INVERSE     
        Array <REAL> t(n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                t[i] += BODY::Bodies[I]->localA[i][j] * BODY::Bodies[I]->localRHS[j];
            }
        BODY::Bodies[I]->localMu = t;
#else

        BODY::Bodies[I]->localMu = pgesv(BODY::Bodies[I]->localA, BODY::Bodies[I]->localRHS);


#endif

        for (int i = 0; i < n; ++i)
            BODY::Bodies[I]->Faces[i].Mu = BODY::Bodies[I]->Faces[i].Gamma = BODY::Bodies[I]->localMu[i];
    }




    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->PhiPrev[3] = BODY::AllBodyFaces[i]->PhiPrev[2];
        BODY::AllBodyFaces[i]->PhiPrev[2] = BODY::AllBodyFaces[i]->PhiPrev[1];
        BODY::AllBodyFaces[i]->PhiPrev[1] = BODY::AllBodyFaces[i]->PhiPrev[0];
        BODY::AllBodyFaces[i]->PhiPrev[0] = BODY::AllBodyFaces[i]->Mu;
    }

    /*  This is when we use the direct face vel calc
    for (int I = 0; I < BODY::Bodies.size(); ++I)
        BODY::Bodies[I]->GetPanelVels();
     */
    REAL Lift = 0.0, incrementalCp = 0.0;
    Vect3 Torque(0.), Force(0.);
    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        REAL Cp = BODY::AllBodyFaces[i]->GetCp();
        BODY::CpHistory[BODY::SubStep - 1][i] = Cp;
        Cp = 0; //BODY::AllBodyFaces[i]->GetCpD();
        BODY::CpHistoryD[BODY::SubStep - 1][i] = Cp;
        Vect3 F = Cp * BODY::AllBodyFaces[i]->Area * BODY::AllBodyFaces[i]->TRANS[2];
        Lift += F.z;
        Vect3 dQ = BODY::AllBodyFaces[i]->dF.Cross(BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG);

        Force += Vect3(BODY::AllBodyFaces[i]->dF.Dot(Vect3(1, 0, 0)), BODY::AllBodyFaces[i]->dF.Dot(Vect3(0, 1, 0)), BODY::AllBodyFaces[i]->dF.Dot(Vect3(0, 0, 1)));
        Torque += dQ;
        incrementalCp += dQ.x * BODY::Bodies[0]->BodyRates.x / (0.5 * BODY::RHO * globalSystem->unscaledVinf.Mag() * globalSystem->unscaledVinf.Mag() * globalSystem->unscaledVinf.Mag() * pi * BODY::Radius * BODY::Radius);
    }

    BODY::ForceHist.push_back(Force);
    BODY::TorqueHist.push_back(Torque);
    //    cout << "Thrust: " << Force[0] << " Torque: " << Torque << " Ct: " << Force[0] / (0.5 * 1027 * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * pi * BODY::Radius * BODY::Radius);
    //    cout << " Cp: " << Torque.x * BODY::Bodies[0]->BodyRates.x / (0.5 * 1027 * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * pi * BODY::Radius * BODY::Radius) << endl;
    BODY::LiftHist.push_back(Lift);
    for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
        BODY::AllProtoWakes[i]->GammaPrev = BODY::AllProtoWakes[i]->Gamma;
        BODY::AllProtoWakes[i]->Gamma =
                (BODY::AllProtoWakes[i]->SheddingSurface->Mu - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->Mu);
    }


    if (BODY::dDeltaCp_dDeltaPhi.size() == 0) {
        BODY::dDeltaCp_dDeltaPhi = Array <REAL > (BODY::AllProtoWakes.size(), 0.0);
        for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
            REAL dCp = BODY::AllProtoWakes[i]->SheddingSurface->GetCp() - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->GetCp();
            REAL dPhi = BODY::AllProtoWakes[i]->Gamma;
            BODY::dDeltaCp_dDeltaPhi[i] = dCp / dPhi;
        }
    }


}

/**************************************************************/
void BODY::LinAlg() {

    Array <REAL> dCP0, dCP1, dPHITE0, dPHITE1;

    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->MuPrev = BODY::AllBodyFaces[i]->GammaPrev = BODY::AllBodyFaces[i]->Mu;
    }

#ifdef USE_MATRIX_INVERSE     
    Array <REAL> t(BODY::Mu.size(), 0.0);

    if (!BODY::IPKC) {
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
            for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
                t[i] += BODY::A[i][j] * BODY::RHS[j];
            }
        //    } else {
        //        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
        //            for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
        //                t[i] += BODY::C[i][j] * BODY::RHS[j];
        //            }
    }

    BODY::Mu = t;
#else

    BODY::Mu = pgesv(BODY::A, BODY::RHS);

#endif


    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->PhiPrev[3] = BODY::AllBodyFaces[i]->PhiPrev[2];
        BODY::AllBodyFaces[i]->PhiPrev[2] = BODY::AllBodyFaces[i]->PhiPrev[1];
        BODY::AllBodyFaces[i]->PhiPrev[1] = BODY::AllBodyFaces[i]->PhiPrev[0];
        BODY::AllBodyFaces[i]->PhiPrev[0] = BODY::AllBodyFaces[i]->Mu;

        BODY::AllBodyFaces[i]->Gamma = BODY::AllBodyFaces[i]->Mu = BODY::Mu[i];

    }


    REAL Lift = 0.0;
    Vect3 Torque(0.), Force(0.);
    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        REAL Cp = BODY::AllBodyFaces[i]->GetCp();
        Vect3 F = Cp * BODY::AllBodyFaces[i]->Area * BODY::AllBodyFaces[i]->TRANS[2];
        Lift += F.z;


        Force += Vect3(BODY::AllBodyFaces[i]->dF.Dot(Vect3(1, 0, 0)), BODY::AllBodyFaces[i]->dF.Dot(Vect3(0, 1, 0)), BODY::AllBodyFaces[i]->dF.Dot(Vect3(0, 0, 1)));
        Torque += BODY::AllBodyFaces[i]->dF.Cross(BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG);

    }

    BODY::ForceHist.push_back(Force);
    BODY::TorqueHist.push_back(Torque);
    //    cout << "Thrust: " << Force[0] << " Torque: " << Torque << " Ct: " << Force[0] / (0.5 * 1027 * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * pi * BODY::Radius * BODY::Radius);
    //    cout << " Cp: " << Torque.x * BODY::Bodies[0]->BodyRates.x / (0.5 * 1027 * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * globalSystem->Vinf.Mag() * pi * BODY::Radius * BODY::Radius) << endl;
    BODY::LiftHist.push_back(Lift);
    for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
        BODY::AllProtoWakes[i]->GammaPrev = BODY::AllProtoWakes[i]->Gamma;
        BODY::AllProtoWakes[i]->Gamma =
                (BODY::AllProtoWakes[i]->SheddingSurface->Mu - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->Mu);
    }


    if (BODY::dDeltaCp_dDeltaPhi.size() == 0) {
        BODY::dDeltaCp_dDeltaPhi = Array <REAL > (BODY::AllProtoWakes.size(), 0.0);
        for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
            REAL dCp = BODY::AllProtoWakes[i]->SheddingSurface->GetCp() - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->GetCp();
            REAL dPhi = BODY::AllProtoWakes[i]->Gamma;
            BODY::dDeltaCp_dDeltaPhi[i] = dCp / dPhi;
        }
    }



    if (BODY::IPKC) {

        REAL NRM = 0;
        dCP1 = Array <REAL > (BODY::AllProtoWakes.size());

        for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
            dCP1[i] = BODY::AllProtoWakes[i]->SheddingSurface->GetCp() - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->GetCp();
            NRM = max(NRM, dCP1[i] * dCP1[i]);
        }



        int cnt = 0, maxits = 50;
        REAL alpha = 0.5; //  underrelaxation factor
        REAL dlta = 1e-9;
        REAL tol = 1e-3;
        while (NRM > tol) {


            //  Get new CPs
            //  Get TE pressure difference.
            dCP0 = Array <REAL > (BODY::AllProtoWakes.size());
            dPHITE0 = Array <REAL > (BODY::AllProtoWakes.size());
            for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
                dCP0[i] = BODY::AllProtoWakes[i]->SheddingSurface->GetCp() - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->GetCp();
                dPHITE0[i] = BODY::AllProtoWakes[i]->Gamma;
            }


            //  Modify wake strength

            for (int i = 0; i < BODY::AllProtoWakes.size(); ++i)
                BODY::AllProtoWakes[i]->Gamma += dlta;

            BODY::GetNonLinearRHS();
#ifdef USE_MATRIX_INVERSE
            Array <REAL> Soln(BODY::Mu.size(), 0.0);

            for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
                for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
                    Soln[i] += BODY::C[i][j] * BODY::RHS[j];
                }

            BODY::Mu = Soln;
#else
            BODY::Mu = pgesv(BODY::C, BODY::RHS);
#endif
            for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
                BODY::AllBodyFaces[i]->Gamma = BODY::AllBodyFaces[i]->Mu = Mu[i];

            //  Get new CPs
            //  Get TE pressure difference.
            NRM = 0;
            dCP1 = Array <REAL > (BODY::AllProtoWakes.size());
            dPHITE1 = Array <REAL > (BODY::AllProtoWakes.size());
            for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
                dCP1[i] = BODY::AllProtoWakes[i]->SheddingSurface->GetCp() - BODY::AllProtoWakes[i]->SheddingSurface->OtherBoundarySurface->GetCp();
                dPHITE1[i] = BODY::AllProtoWakes[i]->Gamma;
                NRM = max(NRM, dCP1[i] * dCP1[i]);
            }
            NRM = sqrt(NRM);
            //  Modify wake strength

            if (NRM > tol) {


                for (int i = 0; i < BODY::AllProtoWakes.size(); ++i)
                    BODY::AllProtoWakes[i]->Gamma -= alpha * (dCP0[i] / ((dCP1[i] - dCP0[i]) / dlta)) + dlta;

                BODY::GetNonLinearRHS();
#ifdef USE_MATRIX_INVERSE
                Soln = 0.0;
                for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
                    for (int j = 0; j < BODY::AllBodyFaces.size(); ++j) {
                        Soln[i] += BODY::C[i][j] * BODY::RHS[j];
                    }

                BODY::Mu = Soln;
#else
                BODY::Mu = pgesv(BODY::C, BODY::RHS);
#endif
                for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
                    BODY::AllBodyFaces[i]->Gamma = BODY::AllBodyFaces[i]->Mu = Mu[i];


                cnt++;

                cout << "%\tIteration: " << cnt << " residual " << NRM << endl;


                if (cnt > maxits) {

                    //                    BODY::Mu = pgesv(BODY::A, BODY::RHS);
                    //
                    //                    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i)
                    //                        BODY::AllBodyFaces[i]->Gamma = BODY::AllBodyFaces[i]->Mu = Mu[i];
                    //                    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i)
                    //                        if (BODY::AllBodyFaces[i]->isBound && !BODY::AllBodyFaces[i]->isTop)
                    //                            BODY::AllBodyFaces[i]->AttachedProtoWake->Gamma -= BODY::AllBodyFaces[i]->Gamma;
                    //

                    NRM = 1e-6;
                }

            }

        }
    }

}

/**************************************************************/
void BODY::BodySubStep(REAL delta_t, int n_steps) {
#ifdef TIME_STEPS
    long unsigned int t11 = ticks();
#endif
    REAL dt = delta_t / n_steps;

    
    
    BODY::CpHistoryAll.push_back(BODY::CpHistory);
    BODY::CpHistory = UTIL::zeros(n_steps, BODY::AllBodyFaces.size());

    BODY::CpHistoryAllD.push_back(BODY::CpHistoryD);
    BODY::CpHistoryD = UTIL::zeros(n_steps, BODY::AllBodyFaces.size());
    BODY::SubTIMES = BODY::Times;

    for (BODY::SubStep = 1; BODY::SubStep <= n_steps; ++BODY::SubStep) {
        cout << BODY::SubStep << " " << n_steps << endl;
        BODY::TimePrev[3] = BODY::TimePrev[2];
        BODY::TimePrev[2] = BODY::TimePrev[1];
        BODY::TimePrev[1] = BODY::TimePrev[0];
        BODY::TimePrev[0] = TIME_STEPPER::SimTime = BODY::Time;

        TIME_STEPPER::SubStepTime += dt;
        TIME_STEPPER::SimTime += dt;
        BODY::Time = TIME_STEPPER::SimTime;
        BODY::Times.push_back(BODY::Time);

        SplitUpLinearAlgebra();


        //        //      Interpolate face vels for subtimestep...
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
            BODY::AllBodyFaces[i]->Vfmm += dt * BODY::AllBodyFaces[i]->dVFMM_dt;
        }
        //  Check the order of the next few lines
        for (int i = 0; i < BODY::Bodies.size(); ++i) {
            BODY::Bodies[i]->MoveBody();
            BODY::Bodies[i]->AngleHist.push_back(BODY::Bodies[i]->EulerAngles);
        }

        //  Sort Wake Convection



        //      Pretend that panels are just vortons for the sake of speed...
        Array <PANEL*> WakePans;
        for (int I = 0; I < BODY::Bodies.size(); ++I)
            for (int i = 0; i < BODY::Bodies[I]->WakePanels.size(); ++i)
                for (int j = 0; j < BODY::Bodies[I]->WakePanels[i].size(); ++j)
                    for (int k = 0; k < BODY::Bodies[I]->WakePanels[i][j].size(); ++k) {
                        WakePans.push_back((BODY::Bodies[I]->WakePanels[i][j][k]));
                    }

        //
        //
        //        Array <Vect3> PanelEdgeMidpoints(4 * WakePans.size()), PanelEdgeVorticities(4 * WakePans.size());
        //
        //        int count = 0;
        //        for (int i = 0; i < WakePans.size(); ++i) {
        //            PanelEdgeMidpoints[count + 0] = WakePans[i]->C1;
        //            PanelEdgeMidpoints[count + 1] = WakePans[i]->C2;
        //            PanelEdgeMidpoints[count + 2] = WakePans[i]->C3;
        //            PanelEdgeMidpoints[count + 3] = WakePans[i]->C4;
        //
        //            PanelEdgeVorticities[count + 0] = 0.5 * WakePans[i]->Gamma * (WakePans[i]->C2 - WakePans[i]->C4);
        //            PanelEdgeVorticities[count + 1] = 0.5 * WakePans[i]->Gamma * (WakePans[i]->C3 - WakePans[i]->C1);
        //            PanelEdgeVorticities[count + 2] = 0.5 * WakePans[i]->Gamma * (WakePans[i]->C4 - WakePans[i]->C2);
        //            PanelEdgeVorticities[count + 3] = 0.5 * WakePans[i]->Gamma * (WakePans[i]->C1 - WakePans[i]->C3);
        //            count += 4;
        //        }
        //
        //
        //        //      Find unique vortices - this is a total bodge, but I cannot be bothered writing proper code
        //
        //        Array <Vect3> Posns, Vorts;
        //        bool test;
        //
        //        for (int i = 0; i < PanelEdgeMidpoints.size(); ++i) {
        //            test = false;
        //            for (int j = 0; j < Posns.size(); ++j) {
        //                if (Posns[j] == PanelEdgeMidpoints[i]) {
        //                    Vorts[j] += PanelEdgeVorticities[i];
        //                    test = true;
        //                    break;
        //                }
        //            }
        //            if (test == false) {
        //                Posns.push_back(PanelEdgeMidpoints[i]);
        //                Vorts.push_back(PanelEdgeVorticities[i]);
        //            }
        //        }
        //
        //        PanelEdgeMidpoints = Posns;
        //        PanelEdgeVorticities = Vorts;
        //
        //
        //
        //
        //        Array <Vect3> V1(WakePans.size()), V2(WakePans.size()), V3(WakePans.size()), V4(WakePans.size());
        //        V1 = V2 = V3 = V4 = globalSystem->unscaledVinf;

        //#ifdef _OPENMP
        //#pragma omp parallel for
        //#endif
        //        for (int i = 0; i < WakePans.size(); ++i) {
        //            for (int j = 0; j < BODY::VortonPositions.size(); ++j) {
        //                V1[i] += globalDirectVel(WakePans[i]->C1 - BODY::VortonPositions[j], BODY::VortonStrengths[j]);
        //                V2[i] += globalDirectVel(WakePans[i]->C2 - BODY::VortonPositions[j], BODY::VortonStrengths[j]);
        //                V3[i] += globalDirectVel(WakePans[i]->C3 - BODY::VortonPositions[j], BODY::VortonStrengths[j]);
        //                V4[i] += globalDirectVel(WakePans[i]->C4 - BODY::VortonPositions[j], BODY::VortonStrengths[j]);
        //            }
        //
        //            for (int j = 0; j < PanelEdgeMidpoints.size(); ++j) {
        //                V1[i] += globalDirectVel(WakePans[i]->C1 - PanelEdgeMidpoints[j], PanelEdgeVorticities[j]);
        //                V2[i] += globalDirectVel(WakePans[i]->C2 - PanelEdgeMidpoints[j], PanelEdgeVorticities[j]);
        //                V3[i] += globalDirectVel(WakePans[i]->C3 - PanelEdgeMidpoints[j], PanelEdgeVorticities[j]);
        //                V4[i] += globalDirectVel(WakePans[i]->C4 - PanelEdgeMidpoints[j], PanelEdgeVorticities[j]);
        //            }
        //        }


        Vect3 dx = dt * globalSystem->unscaledVinf;


        for (int i = 0; i < WakePans.size(); ++i) {
            WakePans[i]->C1 += dx;
            WakePans[i]->C2 += dx;
            WakePans[i]->C3 += dx;
            WakePans[i]->C4 += dx;
            WakePans[i]->GetNormal();
        }



        //        Array <Vect3> Vels(BODY::VortexPositions.size(), Vect3(globalSystem->unscaledVinf));


        //#ifdef _OPENMP
        //#pragma omp parallel for
        //#endif
        //        for (int i = 0; i < BODY::VortexPositions.size(); ++i) {
        //            for (int j = 0; j < BODY::VortexPositions.size(); ++j)
        //                Vels[i] += globalDirectVel(BODY::VortexPositions[i] - BODY::VortexPositions[j], BODY::VortexOmegas[j]);
        //
        //
        //        for (int i = 0; i < BODY::Bodies.size(); ++i)
        //            for (int j = 0; j < PosPtr.size(); ++j)
        //                Vels[i] += BODY::Bodies[i]->GetVel(*(PosPtr[j]));



        for (int i = 0; i < BODY::VortexPositions.size(); ++i)
            BODY::VortexPositions[i] += dx;



        for (int i = 0; i < NumBodies; ++i)
            BODY::Bodies[i]->MakeWake();


    }


    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->PhiWakePrev = 0.0; //BODY::AllBodyFaces[i]->Phi;
        //        BODY::AllBodyFaces[i]->VWakePrev = BODY::AllBodyFaces[i]->VWake;
    }
    
#ifdef TIME_STEPS
    long unsigned int t12 = ticks();
    stringstream tmp;
    tmp << "BodySubStep              : " << double(t12 - t11) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void BODY::SetUpInfluenceMatrices() {
    //BODY::A = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
    //BODY::B = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
    //if (BODY::IPKC)
    //    BODY::C = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
    BODY::RHS = Array <REAL > (BODY::NumFaces, 0.0);
    BODY::Mu = Array <REAL > (BODY::NumFaces, 0.0);
    BODY::Sigma = Array <REAL > (BODY::NumFaces, 0.0);



    BODY::UpdateGlobalInfluenceMatrices();


    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
#ifdef USEGSL
        BODY::AllBodyFaces[i]->sigma = gsl_vector_ptr(globalSigma, i);
        BODY::AllBodyFaces[i]->mu = gsl_vector_ptr(globalMu, i);
#else
        AllBodyFaces[i]->Sigma = BODY::Sigma[i];
        AllBodyFaces[i]->Mu = BODY::Mu[i];
#endif
    }


}

/**************************************************************/
void BODY::UpdateGlobalInfluenceMatrices() {
    long int t = ticks();
    cout << "%\tCalculating influence coefficients..." << endl;
    /*    Array <Array <Vect3> > tmp(1);
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
            BODY::AllBodyFaces[i]->GetNormal();
            if (BODY::AllBodyFaces[i]->isBound)
                BODY::AllBodyFaces[i]->AttachedProtoWake->GetNormal();
        }

            for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
                tmp[0].push_back(BODY::AllBodyFaces[i]->TRANS[2]);
                BODY::AllBodyFaces[i]->GetNormal();
                PANEL *trg = BODY::AllBodyFaces[i];
        #pragma omp parallel for
                for (int j = 0; j < (int) BODY::AllBodyFaces.size(); ++j) {
                    REAL a = 0, b = 0, c = 0;
                    REAL PhiS = 0.0, PhiD = 0.0;
                    PANEL *src = BODY::AllBodyFaces[j];
                    PANEL::SourceDoubletPotential(src, trg, PhiD, PhiS, i, j);

                    //            Vect3 V = 0.0;
                    //            PhiS = 0.0, PhiD = 0.0;
                    //            src->SubPan(100, trg->CollocationPoint, 1, 0, PhiD, PhiS, V);
                    a = PhiD;


                    //            V = 0.0;
                    //            PhiS = 0.0, PhiD = 0.0;
                    //            src->SubPan(100, trg->CollocationPoint, 0, 1, PhiD, PhiS, V);
                    b = PhiS;



                    c = PhiD;
                    if (BODY::AllBodyFaces[j]->isBound) {
                        PhiS = 0.0, PhiD = 0.0;
                        src = BODY::AllBodyFaces[j]->AttachedProtoWake;

                        PANEL::SourceDoubletPotential(src, trg, PhiD, PhiS, 1, 2);
                        a += (BODY::AllBodyFaces[j]->isTop) ? PhiD : -PhiD;
                        //                V = 0.0;
                        //                PhiS = 0.0, PhiD = 0.0;
                        //                src->SubPan(100, trg->CollocationPoint, 1, 0, PhiD, PhiS, V);
                        //                a += (BODY::AllBodyFaces[j]->isTop) ? PhiD : -PhiD;
                    }


                    A[i][j] = a;
                    B[i][j] = b;
                    if (BODY::IPKC)
                        C[i][j] = c;

                }
            }
            cout << "%\tTime elapsed: " << ticks() - t << ". Inverting " << A.size() << "x" << A[0].size() << " matrices: " << endl;
            t = ticks();
        #ifdef USE_MATRIX_INVERSE
            inverse(A);
            if (BODY::IPKC)
                inverse(C);
        #endif
            cout << "%\tTime elapsed: " << ticks() - t << endl;
 
     */

    //        UTIL::WriteMATLABMatrix2D("A","Output.mat", A);
    //        UTIL::WriteMATLABMatrix2D("B","Output.mat", B);
    //        UTIL::WriteMATLABMatrix2D("C","Output.mat", C);


    cout << "\tCalculating local A and B matrices... " << endl;
    t = ticks();
    REAL TempFarField = PANEL::FarField;
    PANEL::FarField = 1e32;
    
    
    int npts = 10;
    Array < Array < Array < Array <REAL> > > > Areas(BODY::Bodies.size());
    Array < Array < Array < Array <Vect3> > > > CP(BODY::Bodies.size()), N(BODY::Bodies.size());
    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        Areas[I].allocate(BODY::Bodies[I]->Faces.size());
        CP[I].allocate(BODY::Bodies[I]->Faces.size());
        N[I].allocate(BODY::Bodies[I]->Faces.size());
        int n = BODY::Bodies[I]->Faces.size();
        for (int j = 0; j < n; ++j) {
            PANEL *src = &BODY::Bodies[I]->Faces[j];
            src->DivPanel(npts, CP[I][j], N[I][j], Areas[I][j]);
        }
    }
    
    
    
    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        cout << "\t Body " << I + 1 << endl;
        int n = BODY::Bodies[I]->Faces.size();
        BODY::Bodies[I]->localA = UTIL::zeros(n, n);
        BODY::Bodies[I]->localB = UTIL::zeros(n, n);
        //BODY::Bodies[I]->localVS = UTIL::zerosv(n, n);
        //BODY::Bodies[I]->localVD = UTIL::zerosv(n, n);
        for (int i = 0; i < n; ++i) {
            PANEL *trg = &BODY::Bodies[I]->Faces[i];
            cout << i << " " << n << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int j = 0; j < n; ++j) {
                REAL PhiS = 0.0, PhiD = 0.0;
                PANEL *src = &BODY::Bodies[I]->Faces[j];
//                PANEL::SourceDoubletPotential(src, trg->CollocationPoint, PhiD, PhiS, i, j);

//                if (i == j) {
                    PANEL::SourceDoubletPotential(src, trg->CollocationPoint, PhiD, PhiS, i, j);
//                    {
//                        for (int k = 0; k < npts; ++k)
//                            for (int l = 0; l < npts; ++l) {
//
//                                Vect3 DX = (CP[I][j][k][l] - trg->CollocationPoint);
//                                REAL DXMag = DX.Mag();
//
//                                PhiS += Areas[I][j][k][l] / (two_pi * DXMag);
//                            }
//                    }
//                                } else
//                PhiS = -src->CurvedSourcePhi(trg->CollocationPoint);
                //                PhiD = -src->CurvedDoubletPhi(trg->CollocationPoint);
                //                
                src->Mu = 1.0;
                src->Gamma = 1.0;
                //cout << PhiD << " " << -src->CurvedDoubletPhi(trg->CollocationPoint) << " " << -src->GetTriTesselatedDoubletPhi(trg->CollocationPoint) << endl;
                if (i!=j)
                    PhiD = -src->GetTriTesselatedDoubletPhi(trg->CollocationPoint);
                else
                    PhiD = 1.0;
                
                src->Mu = .0;
               
                src->Gamma = .0;
//              
                
                REAL a = PhiD; 
                REAL b = PhiS;

 
                if (BODY::Bodies[I]->Faces[j].isBound) {
                    PhiS = 0.0, PhiD = 0.0;
                    src = BODY::Bodies[I]->Faces[j].AttachedProtoWake;

                    src->Mu = src->Gamma = 1.0;
                    PhiD = -src->GetTriTesselatedDoubletPhi(trg->CollocationPoint);
                    src->Mu = src->Gamma = 0.0;
                    a += (BODY::Bodies[I]->Faces[j].isTop) ? PhiD : -PhiD;
                }

                BODY::Bodies[I]->localA[i][j] = a;
                BODY::Bodies[I]->localB[i][j] = b;

            }

        }
#ifdef USE_MATRIX_INVERSE
        UTIL::WriteMATLABMatrix2D("Amatrix", "Output.mat", BODY::Bodies[I]->localA);
        UTIL::WriteMATLABMatrix2D("Bmatrix", "Output.mat", BODY::Bodies[I]->localB);
        inverse(BODY::Bodies[I]->localA);
#endif

    }
    PANEL::FarField = TempFarField;
    cout << "%\tTime elapsed: " << ticks() - t << endl;
}

/**************************************************************/
void BODY::PollFaces() {
    BODY::AllBodyFaces.clear();
    BODY::NumFaces = 0;
    for (int i = 0; i < BODY::Bodies.size(); ++i) {

        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j) {
            BODY::AllBodyFaces.push_back(&(BODY::Bodies[i]->Faces[j]));
            BODY::Bodies[i]->Faces[j].GetNormal();
            BODY::Bodies[i]->Faces[j].GetEdgeInfo();
        }
    }
    NumFaces = BODY::AllBodyFaces.size();
}

/**************************************************************/
void BODY::SetUpProtoWakes(REAL dt) {



    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        Array <PANEL> PWtmp;
        Array <Vect3> ProtoWakePoint1, ProtoWakePoint2;
        Array <REAL> LL, LR;

        BODY::Time = -dt * globalSystem->DS;
        BODY::Bodies[I]->MoveBody();
        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
            P1 = P1 + Vect3(globalSystem->unscaledVinf) * dt * globalSystem->DS;
            P2 = P2 + Vect3(globalSystem->unscaledVinf) * dt * globalSystem->DS;

            ProtoWakePoint1.push_back(P1);
            ProtoWakePoint2.push_back(P2);

        }

        BODY::Time = 0.0;
        BODY::Bodies[I]->MoveBody();



        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
            if (BODY::Bodies[I]->BoundaryFaces[i]->isTop) {


                PWtmp.push_back(PANEL(P1, ProtoWakePoint1[i], ProtoWakePoint2[i], P2));
                PWtmp.back().SheddingSurface = BODY::Bodies[I]->BoundaryFaces[i];
                LL.push_back((P1 - ProtoWakePoint1[i]).Mag());
                LR.push_back((P2 - ProtoWakePoint2[i]).Mag());
            }
        }


        //        PWtmp.clear();






        //        Array <int> NeighbLeft(4), NeighbRight(4);
        //        NeighbLeft[0] = 3;
        //        NeighbRight[0] = 1;
        //        NeighbLeft[1] = 0;
        //        NeighbRight[1] = 2;
        //        NeighbLeft[2] = 1;
        //        NeighbRight[2] = 3;
        //        NeighbLeft[3] = 2;
        //        NeighbRight[3] = 0;
        //
        //        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
        //            PANEL *tmpPan = BODY::Bodies[I]->BoundaryFaces[i];
        //
        //            Vect3 P1 = tmpPan->edgeX1;
        //            Vect3 P2 = tmpPan->edgeX2;
        //
        //            Vect3 N1 = tmpPan->TRANS[2];
        //            Vect3 N2 = tmpPan->OtherBoundarySurface->TRANS[2];
        //            Vect3 N = 0.5 * (N1 + N2);
        //            N = N / N.Mag();
        //
        //            Vect3 LeftNormal = N;
        //            Vect3 RightNormal = N;
        //
        //
        //            REAL LeftL = 0.0;
        //            REAL RightL = 0.0;
        //
        //            if (tmpPan->BoundBC == 0) {
        //                LeftL = (tmpPan->C4 - tmpPan->C1).Mag();
        //                RightL = (tmpPan->C3 - tmpPan->C2).Mag();
        //            }
        //
        //            if (tmpPan->BoundBC == 1) {
        //                LeftL = (tmpPan->C1 - tmpPan->C2).Mag();
        //                RightL = (tmpPan->C3 - tmpPan->C4).Mag();
        //            }
        //
        //            if (tmpPan->BoundBC == 2) {
        //                LeftL = (tmpPan->C3 - tmpPan->C2).Mag();
        //                RightL = (tmpPan->C1 - tmpPan->C4).Mag();
        //            }
        //
        //            if (tmpPan->BoundBC == 3) {
        //                LeftL = (tmpPan->C4 - tmpPan->C3).Mag();
        //                RightL = (tmpPan->C2 - tmpPan->C1).Mag();
        //            }
        //
        //
        //
        //
        //
        //
        //            for (int j = 0; j < 4; ++j) {
        //                if (tmpPan->BoundBC == j) {
        //                    if (tmpPan->Neighb[NeighbLeft[j]] && tmpPan->Neighb[NeighbLeft[j]]->isBound) {
        //                        Vect3 NNormal = tmpPan->Neighb[NeighbLeft[j]]->TRANS[2] + tmpPan->Neighb[NeighbLeft[j]]->OtherBoundarySurface->TRANS[2];
        //                        NNormal = NNormal / NNormal.Mag();
        //                        LeftNormal = 0.5 * (LeftNormal + NNormal);
        //                        LeftNormal = LeftNormal / LeftNormal.Mag();
        //                    }
        //
        //                    if (tmpPan->Neighb[NeighbRight[j]] && tmpPan->Neighb[NeighbRight[j]]->isBound) {
        //                        Vect3 NNormal = tmpPan->Neighb[NeighbRight[j]]->TRANS[2] + tmpPan->Neighb[NeighbRight[j]]->OtherBoundarySurface->TRANS[2];
        //                        NNormal = NNormal / NNormal.Mag();
        //                        RightNormal = 0.5 * (RightNormal + NNormal);
        //                        RightNormal = RightNormal / RightNormal.Mag();
        //                    }
        //                }
        //            }
        //            
        //            if (BODY::Bodies[I]->BoundaryFaces[i]->isTop) {
        //                PWtmp.push_back(PANEL(P1, P1 + LL[count]*LeftNormal, P2 + LR[count]*RightNormal, P2));
        //                PWtmp.back().SheddingSurface = BODY::Bodies[I]->BoundaryFaces[i];
        //                count++;
        //            }
        //        }

        for (int i = 0; i < PWtmp.size(); ++i)
            for (int j = 0; j < PWtmp.size(); ++j)
                PWtmp[i].CheckNeighb(&PWtmp[j]);


        //        PWtmp.clear();
        //
        //        BODY::Bodies[I]->MoveBody(-dt * DS);
        //        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
        //            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
        //            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
        //            P1 = P1 + 0.0 * DS * dt * Vect3(BODY::Vinf);
        //            P2 = P2 + 0.0 * DS * dt * Vect3(BODY::Vinf);
        //
        //            ProtoWakePoint1.push_back(P1);
        //            ProtoWakePoint2.push_back(P2);
        //
        //        }
        //
        //
        //        BODY::Bodies[I]->MoveBody(dt * DS);
        //        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
        //            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
        //            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
        //            if (BODY::Bodies[I]->BoundaryFaces[i]->isTop) {
        //                PWtmp.push_back(PANEL(P1, ProtoWakePoint1[i], ProtoWakePoint2[i], P2));
        //                PWtmp.back().SheddingSurface = BODY::Bodies[I]->BoundaryFaces[i];
        //            }
        //        }
        //
        //
        //        for (int i = 0; i < PWtmp.size(); ++i)
        //            for (int j = 0; j < PWtmp.size(); ++j)
        //                PWtmp[i].CheckNeighb(&PWtmp[j]);


        int PWcnt = 0;
        for (int i = 0; i < PWtmp.size(); ++i)
            if (!PWtmp[i].Neighb.T)
                PWcnt++;

        BODY::Bodies[I]->ProtoWakes.allocate(PWcnt);
        BODY::Bodies[I]->WakePanels.allocate(PWcnt);

        PWcnt = 0;
        for (int i = 0; i < PWtmp.size(); ++i)
            if (!PWtmp[i].Neighb.T) {
                
                PANEL *tmp = &PWtmp[i];
                while (tmp) {
                    BODY::Bodies[I]->ProtoWakes[PWcnt].push_back(PANEL(*tmp));
                    tmp = tmp->Neighb.B;
                }
                PWcnt++;
            }






        for (int i = 0; i < BODY::Bodies[I]->ProtoWakes.size(); ++i) {
            for (int j = 0; j < BODY::Bodies[I]->ProtoWakes[i].size(); ++j) {
                for (int k = 0; k < BODY::Bodies[I]->ProtoWakes[i].size(); ++k)
                    BODY::Bodies[I]->ProtoWakes[i][j].CheckNeighb(&BODY::Bodies[I]->ProtoWakes[i][k]);

                BODY::Bodies[I]->ProtoWakes[i][j].GetEdgeInfo();
                BODY::Bodies[I]->ProtoWakes[i][j].SheddingSurface->AttachedProtoWake = &(BODY::Bodies[I]->ProtoWakes[i][j]);
                BODY::Bodies[I]->ProtoWakes[i][j].SheddingSurface->OtherBoundarySurface->AttachedProtoWake = &(BODY::Bodies[I]->ProtoWakes[i][j]);
                BODY::Bodies[I]->ProtoWakes[i][j].Owner = BODY::Bodies[I];
                BODY::Bodies[I]->ProtoWakes[i][j].isWake = true;
                BODY::Bodies[I]->ProtoWakes[i][j].C1o = VectMultMatrixTranspose(BODY::Bodies[I]->ProtoWakes[i][j].Owner->TRANS, BODY::Bodies[I]->ProtoWakes[i][j].C1o);
                BODY::Bodies[I]->ProtoWakes[i][j].C2o = VectMultMatrixTranspose(BODY::Bodies[I]->ProtoWakes[i][j].Owner->TRANS, BODY::Bodies[I]->ProtoWakes[i][j].C2o);
                BODY::Bodies[I]->ProtoWakes[i][j].C3o = VectMultMatrixTranspose(BODY::Bodies[I]->ProtoWakes[i][j].Owner->TRANS, BODY::Bodies[I]->ProtoWakes[i][j].C3o);
                BODY::Bodies[I]->ProtoWakes[i][j].C4o = VectMultMatrixTranspose(BODY::Bodies[I]->ProtoWakes[i][j].Owner->TRANS, BODY::Bodies[I]->ProtoWakes[i][j].C4o);
            }

            Array <Vect3> Pts;
            for (int j = 0; j < BODY::Bodies[I]->ProtoWakes[i].size(); ++j) {

                BODY::AllProtoWakes.push_back(&BODY::Bodies[I]->ProtoWakes[i][j]);
                BODY::AllProtoWakes.back()->GetNormal();
                if (!BODY::Bodies[I]->ProtoWakes[i][j].Neighb.B) {
                    PANEL *tmp = &(BODY::Bodies[I]->ProtoWakes[i][j]);
                    Pts.push_back(tmp->C2);
                    while (tmp) {
                        Pts.push_back(tmp->C3);
                        tmp = tmp->Neighb.T;
                    }
                }
            }

            Pts.clear();

        }
    }

    //    BODY::D = zeros(BODY::AllBodyFaces.size(),BODY::AllProtoWakes.size());


}

/**************************************************************/
void BODY::GetEulerRates() {
    REAL cosphi = cos(EulerAngles.x), tanthe = tan(EulerAngles.y + 1e-16);
    REAL sinphi = sin(EulerAngles.x), secthe = 1. / (cos(EulerAngles.y) + 1e-16);
    Array <Vect3> vTransform(3);
    vTransform[0] = Vect3(1., tanthe*sinphi, tanthe * cosphi);
    vTransform[1] = Vect3(0., cosphi, -sinphi);
    vTransform[2] = Vect3(0., secthe*sinphi, secthe * cosphi);
    EulerRates = VectMultMatrix(vTransform, BodyRates);
}

/**************************************************************/
void BODY::GetBodyRates() {
    REAL cosphi = cos(EulerAngles.x), costhe = cos(EulerAngles.y);
    REAL sinphi = sin(EulerAngles.x), sinthe = sin(EulerAngles.y);
    Array <Vect3> vTransform(3);
    vTransform[0] = Vect3(1., 0., -sinthe);
    vTransform[1] = Vect3(0., costhe, costhe * sinphi);
    vTransform[2] = Vect3(0., -sinphi, costhe * cosphi);
    BodyRates = VectMultMatrix(vTransform, EulerRates);
}

/**************************************************************/
void BODY::SetEulerTrans() {
    REAL cosphi = cos(EulerAngles.x), costhe = cos(EulerAngles.y), cospsi = cos(EulerAngles.z);
    REAL sinphi = sin(EulerAngles.x), sinthe = sin(EulerAngles.y), sinpsi = sin(EulerAngles.z);
    REAL a1 = costhe*cospsi;
    REAL a2 = costhe*sinpsi;
    REAL a3 = -sinthe;
    REAL b1 = sinphi * sinthe * cospsi - cosphi*sinpsi;
    REAL b2 = sinphi * sinthe * sinpsi + cosphi*cospsi;
    REAL b3 = sinphi*costhe;
    REAL c1 = cosphi * sinthe * cospsi + sinphi*sinpsi;
    REAL c2 = cosphi * sinthe * sinpsi - sinphi*cospsi;
    REAL c3 = cosphi*costhe;
    TRANS[0] = Vect3(a1, b1, c1);
    TRANS[1] = Vect3(a2, b2, c2);
    TRANS[2] = Vect3(a3, b3, c3);
}

/**************************************************************/
void BODY::MoveBody() {

    //    First off calculate new Euler angles


    //    if (BODY::Time > 0.0) {
    //        EulerRates = Vect3(0.0, BODY::AlphaDotHistory.back(), 0.0);
    //        EulerAngles = Vect3(0.0, BODY::AlphaHistory.back(), 0.0);
    //    }
    //    Now get new cg position in global reference frame
    EulerAngles = EulerRates * BODY::Time;
    CG = Velocity * BODY::Time;
    //    Now set appropriate body rates, and transforms etc.
    SetEulerTrans();


    for (int i = 0; i < Faces.size(); ++i) 
        Faces[i].GetNewGlobalPosition();

    ProtoWakeLastC1 = ProtoWakeLastC2 = ProtoWakeLastC3 = ProtoWakeLastC4 = Array < Array < Vect3 > > (ProtoWakes.size());
    for (int i = 0; i < ProtoWakes.size(); ++i) {
        ProtoWakeLastC1[i] = ProtoWakeLastC2[i] = ProtoWakeLastC3[i] = ProtoWakeLastC4[i] = Array <Vect3 > (ProtoWakes[i].size());
        for (int j = 0; j < ProtoWakes[i].size(); ++j) {
            ProtoWakeLastC1[i][j] = ProtoWakes[i][j].C1;
            ProtoWakeLastC2[i][j] = ProtoWakes[i][j].C2;
            ProtoWakeLastC3[i][j] = ProtoWakes[i][j].C3;
            ProtoWakeLastC4[i][j] = ProtoWakes[i][j].C4;
            ProtoWakes[i][j].GetNewGlobalPosition();
        }
    }
}

/**************************************************************/
void BODY::GetPanels(Array <PANEL> &PANELS) {
    //  For all panels in the group
    Faces = PANELS;
    Rmax = 0;
    for (int i = 0; i < PANELS.size(); ++i) {
        PANELS[i].Owner = this;
        Rmax = max(Rmax, PANELS[i].Centroid);
    }

    for (int i = 0; i < (int) Faces.size(); ++i)
        if (Faces[i].isBound)
            BoundaryFaces.push_back(&Faces[i]);



    NumFaces = (int) Faces.size();

    BODY::Mu = Array <REAL > (NumFaces, 0.0);


    for (int i = 0; i < (int) BoundaryFaces.size(); ++i) {
        Vect3 xs1 = (BoundaryFaces[i]->edgeX1);
        Vect3 xs2 = (BoundaryFaces[i]->edgeX2);
        for (int j = 0; j < (int) BoundaryFaces.size(); ++j) {
            if (i != j) {
                Vect3 xt1 = (BoundaryFaces[j]->edgeX1);
                Vect3 xt2 = (BoundaryFaces[j]->edgeX2);
                //  Check the end-points of each boundary edge against all the others

                if (((xs1 == xt1) && (xs2 == xt2)) || ((xs1 == xt2) && (xt1 == xs2))) {
                    BoundaryFaces[i]->OtherBoundarySurface = BoundaryFaces[j];
                    BoundaryFaces[i]->isTop = true;
                    BoundaryFaces[j]->OtherBoundarySurface = BoundaryFaces[i];
                    BoundaryFaces[j]->isTop = false;
                    //  Theres no real point optimising this bit - it's already very fast
                }
            }
        }
    }

    if (NumFaces == BoundaryFaces.size()) {
        cout << "┌\t\t\t\t\t\t\t\t┐" << endl;
        cout << "|\tNumber of wake shedding elements = number of elements\t|" << endl;
        cout << "|\tUsing Lifting Line Mode\t\t\t\t\t|" << endl;
        cout << "└\t\t\t\t\t\t\t\t┘" << endl;

        BODY::LiftingLineMode = true;
        for (int i = 0; i < (int) BoundaryFaces.size(); ++i) {
            BoundaryFaces[i]->isTop = true;
            BoundaryFaces[i]->OtherBoundarySurface = NULL;
        }
    }
}

/**************************************************************/
void BODY::ReadNeuGetBodies(string neu_file, string name, Vect3 dpos, Vect3 cg, Vect3 vel, Vect3 att, Vect3 rates, bool mirror, int plane) {
    //  The trick here (and what must be done) is to set it that all data read from files is combined AT THE END of the function
    //  call into the global (BODY::etc) variables, without prematurly overwriting already existant data. The ID numbering system
    //  for points and panels can remain locally independant, but when it is combined they point and panel ID numbers must include
    //  all previously numbered panels. This should be easy.....


    //  It might also be an idea to make each BODY::Body() have it's own list of panels with indendant IDs, thus the requirement for a
    //  global ID list is reduced or eliminated.



    {
        cout << "%\t" << "Reading " << neu_file << "..." << endl;
    }
    //  Some local variables.
    Array <Vect3> X;
    Array <Array <int> > PNLS, GROUPS, BCS;
    Array < Array < Array < int> > > Surfs;
    Array < Array <REAL> > LCrds;
    Array < Array <Array <int> > > PtIDS_local;
    Array < Array <Array <int> > > TipInboardUSIDS_local;
    Array < Array <Array <int> > > TipInboardLSIDS_local;
    Array < Array <Array <int> > > TipOutboardUSIDS_local;
    Array < Array <Array <int> > > TipOutboardLSIDS_local;
    Array < Array < Array < int > > > InnerTipUSPanelIDS_local, OuterTipUSPanelIDS_local;
    Array < Array < Array < int > > > InnerTipLSPanelIDS_local, OuterTipLSPanelIDS_local;
    Array <string> Names;

    int ngrps = UTIL::read_neu(neu_file, X, PNLS, GROUPS, BCS, Names, Surfs, LCrds,
            PtIDS_local,
            TipInboardUSIDS_local,
            TipInboardLSIDS_local,
            TipOutboardUSIDS_local,
            TipOutboardLSIDS_local,
            InnerTipUSPanelIDS_local,
            InnerTipLSPanelIDS_local,
            OuterTipUSPanelIDS_local,
            OuterTipLSPanelIDS_local);


    for (int i = 0; i < X.size(); ++i)
        X[i] = SYSTEM::GambitScale * X[i];

    BODY::PointsAsRead = X;
    BODY::PanelsAsRead = PNLS;

    //  now make sure surfaces actually match ID numbers.....

    for (int i = 0; i < Surfs.size(); ++i)
        for (int j = 0; j < Surfs[i].size(); ++j)
            for (int k = 0; k < Surfs[i][j].size(); ++k)
                Surfs[i][j][k] += BODY::BodyPanelIDCounter;

    BODY::Surfaces.push_back(Surfs);
    BODY::LocalChordRadius.push_back(LCrds);
    BODY::PtIDS.push_back(PtIDS_local + BODY::AllBodyPoints.size());


    //  Not sure how important it is to get the following correct ....
    BODY::TipInboardUSIDS.push_back(TipInboardUSIDS_local);
    BODY::TipInboardLSIDS.push_back(TipInboardLSIDS_local);
    BODY::TipOutboardUSIDS.push_back(TipOutboardUSIDS_local);
    BODY::TipOutboardLSIDS.push_back(TipOutboardLSIDS_local);
    BODY::InnerTipUSPanelIDS.push_back(InnerTipUSPanelIDS_local);
    BODY::InnerTipLSPanelIDS.push_back(InnerTipLSPanelIDS_local);
    BODY::OuterTipUSPanelIDS.push_back(OuterTipUSPanelIDS_local);
    BODY::OuterTipLSPanelIDS.push_back(OuterTipLSPanelIDS_local);

    //  .... to here.

    Array <Vect3> BodyPoints;
    Array < Array <PANEL> > BodyPanels(ngrps);


    for (int i = 0; i < (int) X.size(); ++i) {
        if (mirror) {
            //  [yz=1, xz=2, xy=3]
            Vect3 a(0.);
            if (plane == 1)
                a = Vect3(1, 0, 0);
            if (plane == 2)
                a = Vect3(0, 1, 0);
            if (plane == 3)
                a = Vect3(0, 0, 1);

            // performing any flip mucks up the normal direction

            X[i] = X[i] - (2 * (X[i].Dot(a)) / (a.Dot(a))) * a;

        }
        Vect3 P = Vect3((X[i] + dpos));
        BodyPoints.push_back(P);
    }

    BODY::AllBodyPoints.push_back(BodyPoints);

    {
        Array <PANEL> tmp(PNLS.size());
        for (int i = 0; i < (int) PNLS.size(); ++i) {
            Vect3 P1, P2, P3, P4;
            P1 = BodyPoints[PNLS[i][0]];
            P2 = BodyPoints[PNLS[i][1]];
            P3 = BodyPoints[PNLS[i][2]];
            P4 = BodyPoints[PNLS[i][3]];

            if (!mirror)
                tmp[i] = PANEL(P1, P2, P3, P4);
            else
                tmp[i] = PANEL(P1, P4, P3, P2);

            tmp[i].ID = BODY::BodyPanelIDCounter + 1;
            BODY::BodyPanelIDCounter++;
            if (LCrds.size() > 0) {
                REAL R1 = LCrds[PNLS[i][0]][1];
                REAL R2 = LCrds[PNLS[i][1]][1];
                REAL R3 = LCrds[PNLS[i][2]][1];
                REAL R4 = LCrds[PNLS[i][3]][1];
                REAL C1 = LCrds[PNLS[i][0]][0];
                REAL C2 = LCrds[PNLS[i][1]][0];
                REAL C3 = LCrds[PNLS[i][2]][0];
                REAL C4 = LCrds[PNLS[i][3]][0];
                tmp[i].Rloc = 0.25 * (R1 + R2 + R3 + R4);
                tmp[i].Cloc = 0.25 * (C1 + C2 + C3 + C4);
            }

        }
        // Get BCs sorted: here we want to create the proto-wake panels
        for (int i = 0; i < (int) BCS.size(); ++i) {
            int j = BCS[i][0];
            tmp[j].isBound = true;
            tmp[j].BoundBC = BCS[i][1];
            if (mirror) {

                if (BCS[i][1] == 0)
                    tmp[j].BoundBC = 3;

                if (BCS[i][1] == 1)
                    tmp[j].BoundBC = 2;

                if (BCS[i][1] == 2)
                    tmp[j].BoundBC = 1;

                if (BCS[i][1] == 3)
                    tmp[j].BoundBC = 0;


            }
            if (tmp[j].BoundBC == 0) {
                tmp[j].edgeX1 = tmp[j].C1;
                tmp[j].edgeX2 = tmp[j].C2;

            }
            if (tmp[j].BoundBC == 1) {
                tmp[j].edgeX1 = tmp[j].C2;
                tmp[j].edgeX2 = tmp[j].C3;
            }
            if (tmp[j].BoundBC == 2) {
                tmp[j].edgeX1 = tmp[j].C3;
                tmp[j].edgeX2 = tmp[j].C4;
            }
            if (tmp[j].BoundBC == 3) {
                tmp[j].edgeX1 = tmp[j].C4;
                tmp[j].edgeX2 = tmp[j].C1;
            }
        }

        //  Sort into groups

        for (int i = 0; i < GROUPS.size(); ++i) {
            for (int j = 0; j < GROUPS[i].size(); ++j)
                tmp[GROUPS[i][j]].Group = i;

            //            
        }


        BodyPanels.allocate(GROUPS.size());
        Array <int> cnt(GROUPS.size(), 0);
        for (int i = 0; i < tmp.size(); ++i)
            cnt[tmp[i].Group]++;

        for (int i = 0; i < GROUPS.size(); ++i)
            BodyPanels[i].allocate(cnt[i]);

        cnt = 0;
        for (int i = 0; i < tmp.size(); ++i) {
            BodyPanels[tmp[i].Group][cnt[tmp[i].Group]] = tmp[i];
            cnt[tmp[i].Group]++;
        }

    }


    cout << "%\t------------|" << endl;
    for (int i = 0; i < GROUPS.size(); ++i) {

        //        ORIGIN *= SYSTEM::GambitScale;
        cout << "%\t" << name << " " << rates << " " << cg << endl;
        cout << "%\tmaking body " << endl;

        BODY::Bodies.push_back(new BODY(cg, att, vel, rates, name));
        BODY::Bodies.back()->ID = BODY::Bodies.size();
        cout << "%\tdone. Getting panels: " << BodyPanels[i].size();

        BODY::Bodies.back()->GetPanels(BodyPanels[i]);
        BODY::NumBodies++;
        cout << "\tdone" << endl;
    }



    /*
        for (int i = 0; i < BodyPoints.size(); ++i)
            BodyPoints[i]->Om -= BodyPoints[i]->Owner->CG.P;

        for (int i = 0; i < (int) BodyPanels.size(); ++i) {
            BodyPanels[i][j].GetNewGlobalPosition();
            BodyPanels[i][j].GetCentroid();
            BodyPanels[i][j].Owner->Rmax = max(BodyPanels[i][j].Owner->Rmax, BodyPanels[i][j].C1->Om.Mag());
            BodyPanels[i][j].Owner->Rmax = max(BodyPanels[i][j].Owner->Rmax, BodyPanels[i][j].C2->Om.Mag());
            BodyPanels[i][j].Owner->Rmax = max(BodyPanels[i][j].Owner->Rmax, BodyPanels[i][j].C3->Om.Mag());
            BodyPanels[i][j].Owner->Rmax = max(BodyPanels[i][j].Owner->Rmax, BodyPanels[i][j].C4->Om.Mag());
        }
     */


    for (int i = 0; i < NumBodies; ++i) {
        cout << "%\t" << BODY::Bodies[i]->Name << " " << BODY::Bodies[i]->NumFaces;
        cout << " " << BODY::Bodies[i]->BoundaryFaces.size() << endl;

        for (int j = 0; j < (int) BODY::Bodies[i]->Faces.size(); ++j) {
            BODY::Bodies[i]->Faces[j].Owner = BODY::Bodies[i];
            BODY::Bodies[i]->Faces[j].GetNewGlobalPosition();
            BODY::Bodies[i]->Faces[j].GetNormal();
        }

    }

    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j)
            for (int k = 0; k < BODY::Bodies[i]->Faces.size(); ++k)
                BODY::Bodies[i]->Faces[j].CheckNeighb(&(BODY::Bodies[i]->Faces[k]));

    for (int i = 0; i < NumBodies; ++i) {
        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j)
            for (int k = 0; k < 4; ++k)
                if (BODY::Bodies[i]->Faces[j].Neighb[k]) {
                    {
                        Vect3 n1 = BODY::Bodies[i]->Faces[j].Normal;
                        Vect3 n2 = BODY::Bodies[i]->Faces[j].Neighb[k]->Normal;
                        Vect3 cp = n1.Cross(n2);
                        BODY::Bodies[i]->Faces[j].Theta[k] = atan2(cp.Mag(), n1.Dot(n2));
                    }
                }
    }

    BODY::WaitLenghts = 0.0;
    for (int i = 0; i < NumBodies; ++i)
        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j)
            BODY::WaitLenghts = max(BODY::WaitLenghts, BODY::Bodies[i]->Faces[j].MaxDiagonal);

    cout << "WaitLength -----------------> " << BODY::WaitLenghts << endl;
}


