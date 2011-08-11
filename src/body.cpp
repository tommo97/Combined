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
#define PANEL_MODE

#include "body.hpp"

/**************************************************************/
#include "body.hpp"
#include "panel.hpp"
#include "pgesv.hpp"

Array <string> BODY::NAMES;
Array <REAL> BODY::Times;
Array <Vect3> BODY::CGS;
Array <Vect3> BODY::RATES;
Array <Vect3> BODY::VELOCITY;
Array <Vect3> BODY::ATTITUDE;
Array <Vect3> BODY::VortonPositions;
Array <Vect3> BODY::VortonStrengths;
Array <REAL> BODY::TimePrev = Array <REAL> (4, 0.0);
REAL BODY::Time = 0;
REAL BODY::RHO = 1;
int BODY::BodyPanelIDCounter = 0, BODY::BodyPointIDCounter = 0;
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
Array <Vect3> BODY::AllBodyPoints;
bool BODY::IPKC = false;
int BODY::NumFaces = 0;
int BODY::NumBodies = 0;
REAL BODY::Radius, BODY::TSR;
bool BODY::LiftingLineMode;

/**************************************************************/
/*  Body functions      */
/**************************************************************/

/**************************************************************/
Vect3 BODY::GetVel(Vect3 Target) {

    //    //    // Why can this section not be multi-threaded????
    Vect3 U(0.0), V(0.0), W(0.0);
    //        for (int i = 0; i < ProtoWakes.size(); ++i)
    //            for (int j = 0; j < ProtoWakes[i].size(); ++j)
    //                U += ProtoWakes[i][j].WakePanelVelocity(Target);


    //    V = GetWakeVel(Target);

    //#pragma omp parallel for
    //        for (int l = 0; l < (int) Faces.size(); ++l) {
    //            W += Faces[l].SourceVel(Target);
    //            W += Faces[l].BodyPanelVelocity(Target); //, globalSystem->Del2);
    //        }

    return U + V + W;
}

/**************************************************************/
void BODY::MakeWake() {
    REAL Gamma;
    Array <Vect3> Pts;
    Array <REAL> Gms;
    for (int i = 0; i < FirstProtoWakes.size(); ++i) {
        PANEL *tmp = FirstProtoWakes[i];
        
        VortonVel[i].push_front(Array < Vect3 > (VortonVel[i].front().size(), Vect3(0.0)));
        VortonX[i].push_front(Array < Vect3 > (VortonX[i].front().size(), Vect3(0.0)));
        VortonOM[i].push_front(Array < Vect3 > (VortonOM[i].front().size(), Vect3(0.0)));


        VortonX[i][0][0] = tmp->C1;

        

        for (int j = 0; j < VortonX[i].front().size() - 1; ++j) {
            VortonX[i][0][j + 1] = tmp->C4;


            Gamma = tmp->Gamma;

            VortonOM[i][0][j] += 0.5 * Gamma * (tmp->C2 - tmp->C1);
            VortonOM[i][0][j] += 0.5 * Gamma * (tmp->C1 - tmp->C4);

            VortonOM[i][1][j] += 0.5 * Gamma * (tmp->C3 - tmp->C2);
            VortonOM[i][1][j] += 0.5 * Gamma * (tmp->C2 - tmp->C1);

            VortonOM[i][1][j + 1] += 0.5 * Gamma * (tmp->C4 - tmp->C3);
            VortonOM[i][1][j + 1] += 0.5 * Gamma * (tmp->C3 - tmp->C2);

            VortonOM[i][0][j + 1] += 0.5 * Gamma * (tmp->C1 - tmp->C4);
            VortonOM[i][0][j + 1] += 0.5 * Gamma * (tmp->C4 - tmp->C3);

            tmp = tmp->Neighb[2];
        }

        if (VortonX[i].size() > 2)
        {
            Array <Vect3> X = VortonX[i][2];
            Array <Vect3> Om = VortonOM[i][2];
            Array <Vect3> Vel = VortonVel[i][2];
            
            VortonX[i][2].clear();
            VortonOM[i][2].clear();
            VortonVel[i][2].clear();
            
            int n = 5;
            for (int j = 0; j < X.size() - 1; ++j)
            {
                Array <Vect3> Xint = UTIL::globalLinspace(X[j],X[j+1],n);
                Array <Vect3> OMint = UTIL::globalLinspace(Om[j],Om[j+1],n);
                
                for (int k = 0; k < n-1; ++k)
                {
                    VortonX[i][2].push_back(Xint[k]);
                    VortonOM[i][2].push_back(OMint[k]/n);
                    VortonVel[i][2].push_back(Vect3(0.));
                    
                    
                }
                
            }
            
            
        }
        
        
        
//
//            VortonX[i][2].clear();
//            VortonOM[i][2].clear();
//            VortonVel[i][2].clear();

        tmp = FirstProtoWakes[i];
        Pts.push_back(tmp->C1);
        while (tmp) {
            Pts.push_back(tmp->C4);
            Gms.push_back(tmp->Gamma);
            tmp = tmp->Neighb[2];
        }



        WakePoints[i].push_front(Pts);
        WakeGamma[i].push_front(Gms);
        Gms.clear();
        Pts.clear();
    }



}

/**************************************************************/
void BODY::GetLinearRHS() {




#pragma omp parallel for
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

#pragma omp parallel for
        for (int j = 0; j < BODY::VortonPositions.size(); ++j)
            VWake = VWake - UTIL::globalDirectVel(BODY::AllBodyFaces[i]->CollocationPoint - BODY::VortonPositions[j], BODY::VortonStrengths[j], globalSystem->Del2);




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
#pragma omp parallel for  
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;

        Vect3 U;

        //  Iterate over all protowake panels
        REAL PhiWake = 0.0;
        for (int j = 0; j < BODY::AllProtoWakes.size(); ++j) {
            U += BODY::AllProtoWakes[j]->WakePanelVelocity(Pos);
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
//    cout << globalSystem->Vinf << endl;
#pragma omp parallel for
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = BODY::AllBodyFaces[i]->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = BODY::AllBodyFaces[i]->Owner->Velocity + Vrot;
        // 	Include freestream and FMM wake interpolation


        Vect3 VWake = BODY::AllBodyFaces[i]->Vfmm;

        REAL PhiWake = 0.0;

#pragma omp parallel for
        for (int j = 0; j < BODY::VortonPositions.size(); ++j)
            VWake = VWake - UTIL::globalDirectVel(BODY::AllBodyFaces[i]->CollocationPoint - BODY::VortonPositions[j], BODY::VortonStrengths[j], globalSystem->Del2);

        BODY::AllBodyFaces[i]->Phi = PhiWake;
        BODY::AllBodyFaces[i]->Vkin = Vkin;
        BODY::AllBodyFaces[i]->VWake = VWake;
        BODY::AllBodyFaces[i]->VCentroid = globalSystem->unscaledVinf - Vkin + VWake;

        Vect3 V = BODY::AllBodyFaces[i]->VCentroid;

        REAL Vn = BODY::AllBodyFaces[i]->Vn = V.Dot(BODY::AllBodyFaces[i]->TRANS[2]);

        BODY::AllBodyFaces[i]->Sigma = Vn;


    }

    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        BODY::AllBodyFaces[i]->MuPrev = BODY::AllBodyFaces[i]->GammaPrev = BODY::AllBodyFaces[i]->Mu;
    }


//    cout << BODY::AllBodyFaces[0]->Vn << endl;

    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        int n = BODY::Bodies[I]->Faces.size();

        BODY::Bodies[I]->localSigma = Array <REAL > (n, 0.0);

        for (int i = 0; i < n; ++i)
            BODY::Bodies[I]->localSigma[i] = BODY::Bodies[I]->Faces[i].Sigma;

        BODY::Bodies[I]->localRHS = Array <REAL > (n, 0.0);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                BODY::Bodies[I]->localRHS[i] += BODY::Bodies[I]->localB[i][j] * BODY::Bodies[I]->localSigma[j];

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

    REAL Lift = 0.0, incrementalCp = 0.0;
    Vect3 Torque(0.), Force(0.);
    for (int i = 0; i < (int) BODY::AllBodyFaces.size(); ++i) {
        REAL Cp = BODY::AllBodyFaces[i]->GetCp();
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

    REAL dt = delta_t / n_steps;
    for (int SubStep = 1; SubStep <= n_steps; ++SubStep) {




        BODY::TimePrev[3] = BODY::TimePrev[2];
        BODY::TimePrev[2] = BODY::TimePrev[1];
        BODY::TimePrev[1] = BODY::TimePrev[0];
        BODY::TimePrev[0] = TIME_STEPPER::SimTime;
        
        TIME_STEPPER::SubStepTime += dt;
        TIME_STEPPER::SimTime += dt;
        BODY::Times.push_back(BODY::Time);


        //cout << "%\t" << SubStep << " " << globalSystem->Vinf << " " << dt << endl;
        for (int i = 0; i < NumBodies; ++i)
            BODY::Bodies[i]->MakeWake();


        ////
        ////
        ////        //  Sort Wake Convection
        ////
        int NumVortons = 0;
        for (int J = 0; J < (int) NumBodies; ++J)
            for (int i = 0; i < BODY::Bodies[J]->VortonX.size(); ++i)
                for (int j = 0; j < BODY::Bodies[J]->VortonX[i].size(); ++j)
                    for (int k = 0; k < BODY::Bodies[J]->VortonX[i][j].size(); ++k)
                        NumVortons++;

        Array <Vect3> Vels(NumVortons, Vect3(globalSystem->unscaledVinf));
        Array <Vect3*> PosPtr(NumVortons);
        int count = 0;
        for (int J = 0; J < (int) NumBodies; ++J)
            for (int i = 0; i < BODY::Bodies[J]->VortonX.size(); ++i)
                for (int j = 0; j < BODY::Bodies[J]->VortonX[i].size(); ++j)
                    for (int k = 0; k < BODY::Bodies[J]->VortonX[i][j].size(); ++k) {
                        PosPtr[count] = (&BODY::Bodies[J]->VortonX[i][j][k]);
                        count++;
                    }












#ifdef InterpPhi

        Vect3 DMin(1e32, 1e32, 1e32), DMax(-1e32, -1e32, -1e32);

        for (int i = 0; i < PosPtr.size(); ++i) {
            DMin = min(DMin, *(PosPtr[i]));
            DMax = max(DMax, *(PosPtr[i]));
        }

        DMin = floor(DMin);
        DMax = ceil(DMax);


        int n = 2;

        int ni = n * ((DMax.x - DMin.x)) + 1;
        int nj = n * ((DMax.y - DMin.y)) + 1;
        int nk = n * ((DMax.z - DMin.z)) + 1;

        Array <REAL> DomainX = UTIL::globalLinspace(DMin.x, DMax.x, ni);
        Array <REAL> DomainY = UTIL::globalLinspace(DMin.y, DMax.y, nj);
        Array <REAL> DomainZ = UTIL::globalLinspace(DMin.z, DMax.z, nk);

        Array < Array < Array <REAL> > > DomainPhi = UTIL::zeros(ni, nj, nk);


        //#pragma omp parallel for
        for (int l = 0; l < (int) BODY::AllBodyFaces.size(); ++l) {
            PANEL *src = BODY::AllBodyFaces[l];
            Vect3 MuT = src->Mu * src->Area * src->TRANS[2];
            REAL SigmaT = src->Sigma * src->Area;
            REAL PhiS = 0.0, PhiD = 0.0;
            Vect3 VT, CP = src->CollocationPoint;
            for (int i = 0; i < ni; ++i)
                for (int j = 0; j < nj; ++j)
                    for (int k = 0; k < nk; ++k) {
                        PhiS = 0.0, PhiD = 0.0;
                        PANEL::PointDoublet(CP, Vect3(DomainX[i], DomainY[j], DomainZ[k]), VT, MuT, PhiD);
                        PANEL::PointSource(CP, Vect3(DomainX[i], DomainY[j], DomainZ[k]), VT, SigmaT, PhiS);
                        DomainPhi[i][j][k] += PhiD + PhiS;
                    }
        }



        cout << ni << " " << nj << " " << nk << " " << DomainPhi[ni - 1][nj - 1][nk - 2] << endl;
#endif





//#pragma omp parallel for
//        for (int i = 0; i < PosPtr.size(); ++i)
//            for (int j = 0; j < BODY::VortonPositions.size(); ++j)
//                Vels[i] += UTIL::globalDirectVel(BODY::VortonPositions[j] - *(PosPtr[i]), BODY::VortonStrengths[j], globalSystem->Del2);
//
//
//
//        for (int i = 0; i < BODY::Bodies.size(); ++i)
//            for (int j = 0; j < PosPtr.size(); ++j)
//                Vels[i] += BODY::Bodies[i]->GetVel(*(PosPtr[j]));



#pragma omp parallel for
        for (int i = 0; i < PosPtr.size(); ++i)
            *(PosPtr[i]) += dt * Vels[i];




        //#pragma omp parallel for
        //        for (int i = 0; i < BODY::VortonPositions.size(); ++i) {
        //            for (int j = 0; j < BODY::VortonPositions.size(); ++j)
        //                WV[i] += UTIL::globalDirectVel(BODY::VortonPositions[i] - BODY::VortonPositions[j], BODY::VortonStrengths[j]);
        //
        //            for (int m = 0; m < BODY::Bodies.size(); ++m)
        //                WV[i] += BODY::Bodies[m]->GetVel(BODY::VortonPositions[i]);
        //        }


        //        for (int i = 0; i < BODY::Bodies.size(); ++i) {
        //            for (int j = 0; j < BODY::Bodies[i]->WakePoints.size(); ++j) {
        //                Array < Array < Vect3 > > WakeVels = BODY::Bodies[i]->WakePoints[j];
        //#pragma omp parallel for
        //                for (int k = 0; k < BODY::Bodies[i]->WakePoints[j].size(); ++k)
        //                    for (int l = 0; l < BODY::Bodies[i]->WakePoints[j][k].size(); ++l) {
        //                        WakeVels[k][l] = Vect3(uinf, vinf, winf);
        ////                        for (int m = 0; m < BODY::Bodies.size(); ++m)
        ////                            WakeVels[k][l] += BODY::Bodies[m]->GetVel(BODY::Bodies[i]->WakePoints[j][k][l]);
        //                    }

        //#pragma omp parallel for
        //                for (int k = 0; k < BODY::Bodies[i]->WakePoints[j].size(); ++k)
        //                    for (int l = 0; l < BODY::Bodies[i]->WakePoints[j][k].size(); ++l)
        //                        BODY::Bodies[i]->WakePoints[j][k][l] += dt * WakeVels[k][l];
        //
        //            }
        //
        //        }


        //  Check the order of the next few lines
        for (int i = 0; i < BODY::Bodies.size(); ++i) {
            BODY::Bodies[i]->MoveBody(dt);
            BODY::Bodies[i]->AngleHist.push_back(BODY::Bodies[i]->EulerAngles);
        }

        SplitUpLinearAlgebra();
//        GetLinearRHS();

        //        if (BODY::Bodies,size() > 1)
        //            BODY::UpdateGlobalInfluenceMatrices();

//        LinAlg();


        //  Get the trailing edge pressure coefficients...

        //        for (int i = 0; i < BODY::Bodies.size(); ++i)
        //            for (int j = 0; j < BODY::Bodies[i]->BoundaryFaces.size(); ++j)
        //                cout << BODY::Bodies[i]->BoundaryFaces[j]->GetCp() - BODY::Bodies[i]->BoundaryFaces[j]->OtherBoundarySurface->GetCp() << endl;






        BODY::VortonPositions.allocate(NumVortons);
        BODY::VortonStrengths.allocate(NumVortons);
        count = 0;
        for (int J = 0; J < (int) NumBodies; ++J)
            for (int i = 0; i < BODY::Bodies[J]->VortonX.size(); ++i)
                for (int j = 0; j < BODY::Bodies[J]->VortonX[i].size(); ++j)
                    for (int k = 0; k < BODY::Bodies[J]->VortonX[i][j].size(); ++k) {
                        BODY::VortonPositions[count] = (BODY::Bodies[J]->VortonX[i][j][k]);
                        BODY::VortonStrengths[count] = (BODY::Bodies[J]->VortonOM[i][j][k]);
                        count++;
                    }
        //        for (int i = 0; i < BODY::VortonPositions.size(); ++i)
        //            BODY::VortonPositions[i] += WV[i] * dt;

        //        for (int i = 0; i < BODY::Bodies.size(); ++i)
        //            BODY::Bodies[i]->WakeToVortons();
    }
    //  Check the order of the next few lines
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        BODY::Bodies[i]->MoveBody(-dt);
}

/**************************************************************/
void BODY::SetUpInfluenceMatrices() {
    BODY::A = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
    BODY::B = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
    if (BODY::IPKC)
        BODY::C = UTIL::zeros(BODY::NumFaces, BODY::NumFaces);
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

    for (int I = 0; I < BODY::Bodies.size(); ++I) {
        cout << "\t Body " << I + 1 << endl;
        int n = BODY::Bodies[I]->Faces.size();
        BODY::Bodies[I]->localA = UTIL::zeros(n, n);
        BODY::Bodies[I]->localB = UTIL::zeros(n, n);

        for (int i = 0; i < n; ++i) {
            PANEL *trg = &BODY::Bodies[I]->Faces[i];
#pragma omp parallel for
            for (int j = 0; j < n; ++j) {
                REAL a = 0, b = 0, c = 0;
                REAL PhiS = 0.0, PhiD = 0.0;
                PANEL *src = &BODY::Bodies[I]->Faces[j];
                PANEL::SourceDoubletPotential(src, trg, PhiD, PhiS, i, j);

                a = PhiD;


                b = PhiS;


                c = PhiD;
                if (BODY::Bodies[I]->Faces[j].isBound) {
                    PhiS = 0.0, PhiD = 0.0;
                    src = BODY::Bodies[I]->Faces[j].AttachedProtoWake;

                    PANEL::SourceDoubletPotential(src, trg, PhiD, PhiS, 1, 2);
                    a += (BODY::Bodies[I]->Faces[j].isTop) ? PhiD : -PhiD;
                }


                BODY::Bodies[I]->localA[i][j] = a;
                BODY::Bodies[I]->localB[i][j] = b;
                //            if (BODY::IPKC)
                //                C[i][j] = c;


            }

        }
#ifdef USE_MATRIX_INVERSE
        inverse(BODY::Bodies[I]->localA);
#endif

    }
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

        BODY::Bodies[I]->MoveBody(-dt * globalSystem->DS);


        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
            P1 = P1 + Vect3(globalSystem->unscaledVinf) * dt*globalSystem->DS;
            P2 = P2 + Vect3(globalSystem->unscaledVinf) * dt*globalSystem->DS;

            ProtoWakePoint1.push_back(P1);
            ProtoWakePoint2.push_back(P2);

        }


        BODY::Bodies[I]->MoveBody(dt * globalSystem->DS);
        for (int i = 0; i < BODY::Bodies[I]->BoundaryFaces.size(); ++i) {
            Vect3 P1 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX1);
            Vect3 P2 = (BODY::Bodies[I]->BoundaryFaces[i]->edgeX2);
            if (BODY::Bodies[I]->BoundaryFaces[i]->isTop) {
                PWtmp.push_back(PANEL(P1, ProtoWakePoint1[i], ProtoWakePoint2[i], P2));
                PWtmp.back().SheddingSurface = BODY::Bodies[I]->BoundaryFaces[i];
            }
        }


        for (int i = 0; i < PWtmp.size(); ++i)
            for (int j = 0; j < PWtmp.size(); ++j)
                PWtmp[i].CheckNeighb(&PWtmp[j]);


        int PWcnt = 0;
        for (int i = 0; i < PWtmp.size(); ++i)
            if (!PWtmp[i].Neighb[2])
                PWcnt++;

        BODY::Bodies[I]->ProtoWakes.allocate(PWcnt);

        PWcnt = 0;
        for (int i = 0; i < PWtmp.size(); ++i)
            if (!PWtmp[i].Neighb[2]) {

                PANEL *tmp = &PWtmp[i];
                while (tmp) {
                    BODY::Bodies[I]->ProtoWakes[PWcnt].push_back(PANEL(*tmp));
                    tmp = tmp->Neighb[0];
                }
                PWcnt++;
            }


        BODY::Bodies[I]->WakeGamma.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->WakePoints.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VortonX.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VortonOM.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VortonVel.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VXSizes.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VOMSizes.allocate(BODY::Bodies[I]->ProtoWakes.size());
        BODY::Bodies[I]->VVelSizes.allocate(BODY::Bodies[I]->ProtoWakes.size());

        BODY::Bodies[I]->FirstProtoWakes.allocate(BODY::Bodies[I]->ProtoWakes.size());
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
                if (!BODY::Bodies[I]->ProtoWakes[i][j].Neighb[0]) {
                    PANEL *tmp = &(BODY::Bodies[I]->ProtoWakes[i][j]);
                    BODY::Bodies[I]->FirstProtoWakes[i] = tmp;
                    Pts.push_back(tmp->C2);
                    while (tmp) {
                        Pts.push_back(tmp->C3);
                        tmp = tmp->Neighb[2];
                    }
                }
            }

            BODY::Bodies[I]->WakePoints[i].push_back(Pts);
            BODY::Bodies[I]->VortonX[i].push_back(Pts);
            BODY::Bodies[I]->VortonOM[i].push_back(Array < Vect3 > (Pts.size(), Vect3(0.0)));
            BODY::Bodies[I]->VortonVel[i].push_back(Array < Vect3 > (Pts.size(), Vect3(0.0)));
            BODY::Bodies[I]->VXSizes[i] = Pts.size();
            BODY::Bodies[I]->VOMSizes[i] = Pts.size();
            BODY::Bodies[I]->VVelSizes[i] = Pts.size();
            Pts.clear();

        }
    }

    //    BODY::D = zeros(BODY::AllBodyFaces.size(),BODY::AllProtoWakes.size());


}

/**************************************************************/
void BODY::GetEulerRates() {
    REAL cosphi = cos(EulerAngles.x), tanthe = tan(EulerAngles.y + 1e-16);
    REAL sinphi = sin(EulerAngles.x), secthe = 1 / (cos(EulerAngles.y) + 1e-16);
    Vect3 vTransform[3] = {Vect3(1, tanthe*sinphi, tanthe * cosphi),
        Vect3(0, cosphi, -sinphi),
        Vect3(0, secthe*sinphi, secthe * cosphi)};
    EulerRates = VectMultMatrix(vTransform, BodyRates);
}

/**************************************************************/
void BODY::GetBodyRates() {
    REAL cosphi = cos(EulerAngles.x), costhe = cos(EulerAngles.y);
    REAL sinphi = sin(EulerAngles.x), sinthe = sin(EulerAngles.y);
    Vect3 vTransform[3] = {Vect3(1, 0, -sinthe),
        Vect3(0, costhe, costhe * sinphi),
        Vect3(0, -sinphi, costhe * cosphi)};
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
void BODY::MoveBody(REAL dt) {

    //    First off calculate new Euler angles
    EulerAngles += dt*EulerRates;
    //    Now get new cg position in global reference frame

    CG += Velocity*dt;
//    cout << CG << endl;
    //    Now set appropriate body rates, and transforms etc.
    SetEulerTrans();


    for (int i = 0; i < Faces.size(); ++i) {
        Faces[i].GetNewGlobalPosition();
    }

    for (int i = 0; i < ProtoWakes.size(); ++i)
        for (int j = 0; j < ProtoWakes[i].size(); ++j) {
            ProtoWakes[i][j].GetNewGlobalPosition();
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

        //        ORIGIN *= GambitScale;
        cout << "%\t" << name << " " << rates << " " << cg << endl;
        cout << "%\tmaking body " << endl;

        BODY::Bodies.push_back(new BODY(cg, att, vel, rates, name));
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
}


