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



#include "includes.hpp"
#include "tree.hpp"
#include "node.hpp"
#include "system.hpp"
#include "waves.hpp"
#include "unit_tests.hpp"
//      maximum memory use before bombing out

using namespace std;
/**************************************************************/

void sheet(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL radius, REAL scale, REAL THETA);

inline bool srt_vect3(Vect3 A, Vect3 B) {
    if (A.x == B.x) {
        if (A.y == B.y)
            return (A.z < B.z);
        else
            return (A.y < B.y);
    } else
        return (A.x < B.x);
}




















/**************************************************************/
int main(int argc, char *argv[]) {

    system("clear");
    system("rm Output.mat");

    BODY::OutputSubStepCollocationPoints = false;
    SYSTEM::PanelOnly = false;
    PANEL::Initialise();
    UTIL::GetCellPans();

//    TEST::SimpleTestPanel();
//    return 0;


    //TEST::SmallTestFMM();
    //return 0;
    //TEST::SolveMatfileVels();
    //return 0;
    //    TEST::TestBulkLoader(100000);
    //        TEST::TestFMM(argc, argv);
    //        TEST::SimpleTestPanel();
    //        TEST::TestBiotSavart();
    //        return 0;
    //    TEST::TestBEMFVM();
    //        return 0;



    TIME_STEPPER::SpinUpTime = 2.5; // time over which starting vortices are linearly damped
    SYSTEM System(0);

    //    globalSystem->dtInit = 0.01; //     This *maybe* gets changed according to the maximum kinematic velocity of the body(s)

    UTIL::cpu_t = ticks();
    //    WaveField::Period = 1.0;//7.72;
    //    WaveField::Height = 0.15;//18.6666*0.3048;
    //    WaveField::Depth = 2.245;//30*0.3048;
    //    WaveField::Cnoidal.WaveFieldCnoidal();
    //    cout << WaveField::Cnoidal.CnoidalVelocity(Vect3(1.2, 0.0, -1.2) ,2.5) << endl;
    //    cout << WaveField::Cnoidal.CnoidalPerturbationPotential(Vect3(1.2, 0.0, -1.2) ,2.5) << endl;
    //    return 0;
    UTIL::PreAmble();



    if (!SYSTEM::PanelOnly) {
        globalSystem->Initialise();

        globalSystem->VortonsXs.clear();
        globalSystem->VortonOmegas.clear();


        globalSystem->TimeStep();

    }
    UTIL::PostAmble(string("Output.mat"));


#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - UTIL::cpu_t) / 1000 << " seconds" << endl;
#endif
}

/**************************************************************/
void UTIL::PostAmble(string fname) {

    
    //  Minimal PostAmble as seemed to take too much time


    {
        int n = BODY::AllBodyFaces.size();

        Array < Array <REAL> > TmpPtsx, TmpPtsy, TmpPtsz;
        TmpPtsx = UTIL::zeros(n, 4);
        TmpPtsy = UTIL::zeros(n, 4);
        TmpPtsz = UTIL::zeros(n, 4);
        Array <REAL> CPs(n), Areas(n), Sigmas(n), Mus(n);

        for (int i = 0; i < n; ++i) {
            TmpPtsx[i][0] = BODY::AllBodyFaces[i]->C1.x;
            TmpPtsy[i][0] = BODY::AllBodyFaces[i]->C1.y;
            TmpPtsz[i][0] = BODY::AllBodyFaces[i]->C1.z;

            TmpPtsx[i][1] = BODY::AllBodyFaces[i]->C2.x;
            TmpPtsy[i][1] = BODY::AllBodyFaces[i]->C2.y;
            TmpPtsz[i][1] = BODY::AllBodyFaces[i]->C2.z;

            TmpPtsx[i][2] = BODY::AllBodyFaces[i]->C3.x;
            TmpPtsy[i][2] = BODY::AllBodyFaces[i]->C3.y;
            TmpPtsz[i][2] = BODY::AllBodyFaces[i]->C3.z;

            TmpPtsx[i][3] = BODY::AllBodyFaces[i]->C4.x;
            TmpPtsy[i][3] = BODY::AllBodyFaces[i]->C4.y;
            TmpPtsz[i][3] = BODY::AllBodyFaces[i]->C4.z;

            CPs[i] = BODY::AllBodyFaces[i]->GetCp();
            Mus[i] = BODY::AllBodyFaces[i]->Mu;
            Areas[i] = BODY::AllBodyFaces[i]->Area;
            Sigmas[i] = BODY::AllBodyFaces[i]->Sigma;

        }

        n = BODY::AllWakePanels.size();
        Array <Vect3> C1(n), C2(n), C3(n), C4(n);
        Array <REAL> GMA(n);

        for (int i = 0; i < n; ++i) {
            C1[i] = BODY::AllWakePanels[i]->C1;
            C2[i] = BODY::AllWakePanels[i]->C2;
            C3[i] = BODY::AllWakePanels[i]->C3;
            C4[i] = BODY::AllWakePanels[i]->C4;
            GMA[i] = BODY::AllWakePanels[i]->Gamma;
        }



        WriteMATLABMatrix2D("BodyPointsX", fname, TmpPtsx);
        WriteMATLABMatrix2D("BodyPointsY", fname, TmpPtsy);
        WriteMATLABMatrix2D("BodyPointsZ", fname, TmpPtsz);
        WriteMATLABMatrix1D("Cp", fname, CPs);
        WriteMATLABMatrix1D("Mu", fname, Mus);
        WriteMATLABMatrix1D("Area", fname, Areas);
        WriteMATLABMatrix1D("Sigma", fname, Sigmas);
        WriteMATLABMatrix1DVect3("WakePanC1", fname, C1);
        WriteMATLABMatrix1DVect3("WakePanC2", fname, C2);
        WriteMATLABMatrix1DVect3("WakePanC3", fname, C3);
        WriteMATLABMatrix1DVect3("WakePanC4", fname, C4);
        WriteMATLABMatrix1D("WakePanGamma", fname, GMA);


        for (int i = 0; i < BODY::Bodies.size(); ++i) {
            UTIL::WriteMATLABMatrix1DVect3("BodyCG" + UTIL::toString(i), fname, BODY::Bodies[i]->CG);
        }

        for (int i = 0; i < BODY::Surfaces.size(); ++i) {
            UTIL::WriteMATLABMatrix2D("BodySurface" + UTIL::toString(i), fname, BODY::Surfaces[i]);
            UTIL::WriteMATLABMatrix2D("BodyMainPointIDS" + UTIL::toString(i), fname, BODY::PtIDS[i]);
        }

        Array <Vect3> AllBodyPoints;
        for (int i = 0; i < BODY::AllBodyPoints0.size(); ++i) {
            Vect3 X0 = BODY::AllBodyPoints0[i]; // nb: already relative to cg0
            Vect3 CG0 = BODY::AllBodyPointCG0[i];
            Vect3 Rates = BODY::AllBodyPointRates0[i];
            Vect3 Angles0 = BODY::AllBodyPointEulerAngles0[i];
            Vect3 Vels0 = BODY::AllBodyPointCGVels0[i];

            Vect3 EulerAngles = UTIL::ODE4Final(BODY::EulerDot, 0.0, BODY::Time, 0.001, Angles0, Rates);

            Vect3 CG = CG0 + BODY::Time * Vels0;
            Array <Vect3> TRANS = BODY::ReturnTrans(EulerAngles);
            Vect3 X = CG + VectMultMatrix(TRANS, X0);
            AllBodyPoints.push_back(X);
        }
        UTIL::WriteMATLABMatrix1DVect3("AllBodyPoints", fname, AllBodyPoints);
        WriteMATLABMatrix1D("Time", fname, BODY::Time);
        WriteMATLABMatrix1D("NBodies", fname, BODY::Bodies.size());
        WriteMATLABMatrix1D("NParts", fname, BODY::Surfaces.size());
        WriteMATLABMatrix1D("GambitScale", fname, SYSTEM::GambitScale);

    }
    

    
    
    return;
    
    
    
    
    
    
    
    WriteMATLABMatrix1DVect3("PointsAsRead", fname, BODY::PointsAsRead);
    WriteMATLABMatrix1D("AlphaHistory", fname, BODY::AlphaHistory);
    WriteMATLABMatrix1D("AlphaDotHistory", fname, BODY::AlphaDotHistory);
    WriteMATLABMatrix2D("CpHistory", fname, BODY::CpHistory);
    WriteMATLABMatrix2D("PanelsAsRead", fname, BODY::PanelsAsRead);
    for (int i = 0; i < BODY::Surfaces.size(); ++i) {
        UTIL::WriteMATLABMatrix2D("BodySurface" + UTIL::toString(i), fname, BODY::Surfaces[i]);
        UTIL::WriteMATLABMatrix2D("BodyMainPointIDS" + UTIL::toString(i), fname, BODY::PtIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipInboardUSPointIDS" + UTIL::toString(i), fname, BODY::TipInboardUSIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipInboardLSPointIDS" + UTIL::toString(i), fname, BODY::TipInboardLSIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipOutboardUSPointIDS" + UTIL::toString(i), fname, BODY::TipOutboardUSIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipOutboardLSPointIDS" + UTIL::toString(i), fname, BODY::TipOutboardLSIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipInboardUSPanelIDS" + UTIL::toString(i), fname, BODY::InnerTipUSPanelIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipInboardLSPanelIDS" + UTIL::toString(i), fname, BODY::InnerTipLSPanelIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipOutboardUSPanelIDS" + UTIL::toString(i), fname, BODY::OuterTipUSPanelIDS[i]);
        UTIL::WriteMATLABMatrix2D("BodyTipOutboardLSPanelIDS" + UTIL::toString(i), fname, BODY::OuterTipLSPanelIDS[i]);
    }

    Array < Array <REAL> > TRANS1, TRANS2, TRANS3, TmpPtsx, TmpPtsy, TmpPtsz, TmpPWPtsx, TmpPWPtsy, TmpPWPtsz, PhiPrev;
    Array <REAL> Mus, Sigmas, PWGammas, CPs, TPrev, Cloc, Rloc, Areas;
    Array <Vect3> Centroids, VCentroids;
    for (int i = 0; i < BODY::Bodies.size(); ++i) {
        UTIL::WriteMATLABMatrix1DVect3("BodyCG" + UTIL::toString(i), fname, BODY::Bodies[i]->CG);
        UTIL::WriteMATLABMatrix1DVect3("BodyDisplacement" + UTIL::toString(i), fname, BODY::Bodies[i]->Disp);
        UTIL::WriteMATLABMatrix1DVect3("BodyVel" + UTIL::toString(i), fname, BODY::Bodies[i]->Velocity);
        UTIL::WriteMATLABMatrix1DVect3("EulerAngles" + UTIL::toString(i), fname, BODY::Bodies[i]->EulerAngles);
        UTIL::WriteMATLABMatrix1DVect3("TRANS1_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[0]);
        UTIL::WriteMATLABMatrix1DVect3("TRANS2_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[1]);
        UTIL::WriteMATLABMatrix1DVect3("TRANS3_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[2]);
        if (BODY::OutputSubStepCollocationPoints) {
            UTIL::WriteMATLABMatrix2DVect3("CollocPtHistory" + UTIL::toString(i), fname, BODY::SubStepCollocPts[i]);
            UTIL::WriteMATLABMatrix2DVect3("PanelNormalHistory" + UTIL::toString(i), fname, BODY::SubStepCollocPtNrms[i]);
        }
        for (int j = 0; j < BODY::Bodies[i]->VortonX.size(); ++j) {
            UTIL::WriteMATLABMatrix2DVect3("VortonPoints" + UTIL::toString(i) + "_" + UTIL::toString(j), fname, BODY::Bodies[i]->VortonX[j]);
            UTIL::WriteMATLABMatrix2DVect3("VortonStrength" + UTIL::toString(i) + "_" + UTIL::toString(j), fname, BODY::Bodies[i]->VortonOM[j]);
        }


        for (int j = 0; j < BODY::Bodies[i]->WakePoints.size(); ++j) {
            UTIL::WriteMATLABMatrix2DVect3("WakePoints" + UTIL::toString(j), fname, BODY::Bodies[i]->WakePoints[j]);
            UTIL::WriteMATLABMatrix2D("WakeGamma" + UTIL::toString(j), fname, BODY::Bodies[i]->WakeGamma[j]);
        }








        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j) {
            Array <REAL> t;
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[0][0]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[0][1]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[0][2]);
            TRANS1.push_back(t);
            t.clear();
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[1][0]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[1][1]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[1][2]);
            TRANS2.push_back(t);
            t.clear();
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[2][0]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[2][1]);
            t.push_back(BODY::Bodies[i]->Faces[j].TRANS[2][2]);
            TRANS3.push_back(t);
            t.clear();
            t.push_back(BODY::Bodies[i]->Faces[j].C1.x);
            t.push_back(BODY::Bodies[i]->Faces[j].C2.x);
            t.push_back(BODY::Bodies[i]->Faces[j].C3.x);
            t.push_back(BODY::Bodies[i]->Faces[j].C4.x);
            TmpPtsx.push_back(t);
            t.clear();
            t.push_back(BODY::Bodies[i]->Faces[j].C1.y);
            t.push_back(BODY::Bodies[i]->Faces[j].C2.y);
            t.push_back(BODY::Bodies[i]->Faces[j].C3.y);
            t.push_back(BODY::Bodies[i]->Faces[j].C4.y);
            TmpPtsy.push_back(t);
            t.clear();
            t.push_back(BODY::Bodies[i]->Faces[j].C1.z);
            t.push_back(BODY::Bodies[i]->Faces[j].C2.z);
            t.push_back(BODY::Bodies[i]->Faces[j].C3.z);
            t.push_back(BODY::Bodies[i]->Faces[j].C4.z);
            TmpPtsz.push_back(t);
            t.clear();
            t = BODY::Bodies[i]->Faces[j].PhiPrev;
            PhiPrev.push_back(t);
            t.clear();
            Mus.push_back(BODY::Bodies[i]->Faces[j].Mu);
            CPs.push_back(BODY::Bodies[i]->Faces[j].GetCp());
            Areas.push_back(BODY::Bodies[i]->Faces[j].Area);
            Sigmas.push_back(BODY::Bodies[i]->Faces[j].Sigma);
            Centroids.push_back(BODY::Bodies[i]->Faces[j].CollocationPoint);
            VCentroids.push_back(BODY::Bodies[i]->Faces[j].VCentroid);
            Rloc.push_back(BODY::Bodies[i]->Faces[j].Rloc);
            Cloc.push_back(BODY::Bodies[i]->Faces[j].Cloc);

        }
    }


    Array <Vect3> C1, C2, C3, C4;
    Array <REAL> GMA;
    for (int I = 0; I < BODY::Bodies.size(); ++I)
        for (int i = 0; i < BODY::Bodies[I]->WakePanels.size(); ++i)
            for (int j = 0; j < BODY::Bodies[I]->WakePanels[i].size(); ++j)
                for (int k = 0; k < BODY::Bodies[I]->WakePanels[i][j].size(); ++k) {
                    C1.push_back(BODY::Bodies[I]->WakePanels[i][j][k]->C1);
                    C2.push_back(BODY::Bodies[I]->WakePanels[i][j][k]->C2);
                    C3.push_back(BODY::Bodies[I]->WakePanels[i][j][k]->C3);
                    C4.push_back(BODY::Bodies[I]->WakePanels[i][j][k]->C4);
                    GMA.push_back(BODY::Bodies[I]->WakePanels[i][j][k]->Gamma);
                }


    for (int i = 0; i < BODY::AllProtoWakes.size(); ++i) {
        C1.push_back(BODY::AllProtoWakes[i]->C1);
        C2.push_back(BODY::AllProtoWakes[i]->C2);
        C3.push_back(BODY::AllProtoWakes[i]->C3);
        C4.push_back(BODY::AllProtoWakes[i]->C4);
        GMA.push_back(BODY::AllProtoWakes[i]->Gamma);
    }
    WriteMATLABMatrix1DVect3("WakePanC1", fname, C1);
    WriteMATLABMatrix1DVect3("WakePanC2", fname, C2);
    WriteMATLABMatrix1DVect3("WakePanC3", fname, C3);
    WriteMATLABMatrix1DVect3("WakePanC4", fname, C4);
    WriteMATLABMatrix1D("WakePanGamma", fname, GMA);
    //    WriteMATLABMatrix2D("Amatrix", fname, BODY::Bodies[0]->localA);
    //    WriteMATLABMatrix2D("Bmatrix", fname, BODY::Bodies[0]->localB);
    WriteMATLABMatrix2D("PhiPrev", fname, PhiPrev);
    WriteMATLABMatrix2D("TRANS1", fname, TRANS1);
    WriteMATLABMatrix2D("TRANS2", fname, TRANS2);
    WriteMATLABMatrix2D("TRANS3", fname, TRANS3);
    WriteMATLABMatrix2D("BodyPointsX", fname, TmpPtsx);
    WriteMATLABMatrix2D("BodyPointsY", fname, TmpPtsy);
    WriteMATLABMatrix2D("BodyPointsZ", fname, TmpPtsz);
    WriteMATLABMatrix1D("Alpha", fname, BODY::ATTITUDE[0].y);
    WriteMATLABMatrix1D("IPKC", fname, (int) BODY::IPKC);
    WriteMATLABMatrix1D("Mu", fname, Mus);
    WriteMATLABMatrix1D("Cp", fname, CPs);
    WriteMATLABMatrix1D("Area", fname, Areas);
    WriteMATLABMatrix1D("Sigma", fname, Sigmas);
    WriteMATLABMatrix1D("Rloc", fname, Rloc);
    WriteMATLABMatrix1D("Cloc", fname, Cloc);
    WriteMATLABMatrix1D("TPrev", fname, BODY::TimePrev);
    WriteMATLABMatrix1D("Time", fname, BODY::Time);
    WriteMATLABMatrix1D("Lift", fname, BODY::LiftHist);
    WriteMATLABMatrix1DVect3("Force", fname, BODY::ForceHist);
    WriteMATLABMatrix1DVect3("Torque", fname, BODY::TorqueHist);
    WriteMATLABMatrix1D("NBodies", fname, BODY::Bodies.size());
    WriteMATLABMatrix1D("NParts", fname, BODY::Surfaces.size());
    WriteMATLABMatrix1DVect3("CollocPts", fname, Centroids);
    WriteMATLABMatrix1DVect3("VCollocPts", fname, VCentroids);
    WriteMATLABMatrix1DVect3("VortonX", fname, BODY::VortexPositions);
    WriteMATLABMatrix1DVect3("VortonO", fname, BODY::VortexOmegas);
    WriteMATLABMatrix1D("VortonOwnerID", fname, BODY::VortexOwnerID);

    {
        Array <Vect3> AllBodyPoints;
        for (int i = 0; i < BODY::AllBodyPoints0.size(); ++i) {
            Vect3 X0 = BODY::AllBodyPoints0[i];         // nb: already relative to cg0
            Vect3 CG0 = BODY::AllBodyPointCG0[i];
            Vect3 Rates = BODY::AllBodyPointRates0[i];
            Vect3 Angles0 = BODY::AllBodyPointEulerAngles0[i];
            Vect3 Vels0 = BODY::AllBodyPointCGVels0[i];

            Vect3 EulerAngles = UTIL::ODE4Final(BODY::EulerDot, 0.0, BODY::Time, 0.001, Angles0, Rates);

            Vect3 CG = CG0 + BODY::Time * Vels0;
            Array <Vect3> TRANS = BODY::ReturnTrans(EulerAngles);
            Vect3 X = CG + VectMultMatrix(TRANS, X0);
            AllBodyPoints.push_back(X);
        }
        UTIL::WriteMATLABMatrix1DVect3("AllBodyPoints", fname, AllBodyPoints);
    }

    
    
    
    
    
    
    
    
    
    //    WriteMATLABMatrix1D("VortonOwnerID", fname, BODY::VortexOwnerID);
    Array <Vect3> AllBodyPointsTransformed(BODY::AllBodyPoints0.size(), Vect3(0.0));

    //        for (int i = 0; i < BODY::AllBodyPoints.size(); ++i)
    //            AllBodyPointsTransformed[i] = VectMultMatrix(BODY::AllBodyFaces->, vCorner_g);

//    UTIL::WriteMATLABMatrix1DVect3("AllBodyPoints", fname, BODY::AllBodyPoints0);
    UTIL::WriteMATLABMatrix2D("ProtoWakePointsX", fname, TmpPWPtsx);
    UTIL::WriteMATLABMatrix2D("ProtoWakePointsY", fname, TmpPWPtsy);
    UTIL::WriteMATLABMatrix2D("ProtoWakePointsZ", fname, TmpPWPtsz);
    UTIL::WriteMATLABMatrix1D("ProtoWakeGamma", fname, PWGammas);

    //        cout << "for i=1:" << TmpPtsx.size() << endl;
    //        SURF(string("BodyPointsX(i,1)"), string("BodyPointsY(i,1)"), string("BodyPointsZ(i,1)"), string("BodyPointsX(i,2)"), string("BodyPointsY(i,2)"), string("BodyPointsZ(i,2)"), string("BodyPointsX(i,3)"), string("BodyPointsY(i,3)"), string("BodyPointsZ(i,3)"), string("BodyPointsX(i,4)"), string("BodyPointsY(i,4)"), string("BodyPointsZ(i,4)"), "Cp(i)", cout)
    //        cout << "end" << endl;
    //        cout << "for i=1:" << TmpPWPtsx.size() << endl;
    //        SURF(string("ProtoWakePointsX(i,1)"), string("ProtoWakePointsY(i,1)"), string("ProtoWakePointsZ(i,1)"), string("ProtoWakePointsX(i,2)"), string("ProtoWakePointsY(i,2)"), string("ProtoWakePointsZ(i,2)"), string("ProtoWakePointsX(i,3)"), string("ProtoWakePointsY(i,3)"), string("ProtoWakePointsZ(i,3)"), string("ProtoWakePointsX(i,4)"), string("ProtoWakePointsY(i,4)"), string("ProtoWakePointsZ(i,4)"), "ProtoWakeGamma(i)", cout)
    //        cout << "end" << endl;


    Array <int> OtherBCIDs(BODY::AllBodyFaces.size(), 0), IDs(BODY::AllBodyFaces.size(), 0), isBound(BODY::AllBodyFaces.size(), 0), BoundEdge(BODY::AllBodyFaces.size(), 0), NeighbIDs4(BODY::AllBodyFaces.size(), 0), NeighbIDs2(BODY::AllBodyFaces.size(), 0), NeighbIDs3(BODY::AllBodyFaces.size(), 0), NeighbIDs1(BODY::AllBodyFaces.size(), 0);
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        IDs[i] = BODY::AllBodyFaces[i]->ID;

        if (BODY::AllBodyFaces[i]->Neighb[3])
            NeighbIDs4[i] = BODY::AllBodyFaces[i]->Neighb[3]->ID;
        else
            NeighbIDs4[i] = -1;
        if (BODY::AllBodyFaces[i]->Neighb[1])
            NeighbIDs2[i] = BODY::AllBodyFaces[i]->Neighb[1]->ID;
        else
            NeighbIDs2[i] = -1;
        if (BODY::AllBodyFaces[i]->Neighb[2])
            NeighbIDs3[i] = BODY::AllBodyFaces[i]->Neighb[2]->ID;
        else
            NeighbIDs3[i] = -1;
        if (BODY::AllBodyFaces[i]->Neighb[0])
            NeighbIDs1[i] = BODY::AllBodyFaces[i]->Neighb[0]->ID;
        else
            NeighbIDs1[i] = -1;

        isBound[i] = BODY::AllBodyFaces[i]->isBound;
        BoundEdge[i] = BODY::AllBodyFaces[i]->BoundBC;
        if (isBound[i])
            OtherBCIDs[i] = BODY::AllBodyFaces[i]->OtherBoundarySurface->ID;
    }

    UTIL::WriteMATLABMatrix1D("IDs", fname, IDs);
    UTIL::WriteMATLABMatrix1D("NeighbIDs4", fname, NeighbIDs4);
    UTIL::WriteMATLABMatrix1D("NeighbIDs2", fname, NeighbIDs2);
    UTIL::WriteMATLABMatrix1D("NeighbIDs3", fname, NeighbIDs3);
    UTIL::WriteMATLABMatrix1D("NeighbIDs1", fname, NeighbIDs1);

    UTIL::WriteMATLABMatrix1D("isBound", fname, isBound);
    UTIL::WriteMATLABMatrix1D("BoundEdge", fname, BoundEdge);
    UTIL::WriteMATLABMatrix1D("OtherBCIDs", fname, OtherBCIDs);
    UTIL::WriteMATLABMatrix1D("Times", fname, BODY::Times);


    for (int i = 0; i < BODY::Bodies.size(); ++i) {
        UTIL::WriteMATLABMatrix1DVect3("BodyRates" + UTIL::toString(i), fname, BODY::Bodies[i]->BodyRates);
        UTIL::WriteMATLABMatrix1DVect3("EulerHist" + UTIL::toString(i), fname, BODY::Bodies[i]->AngleHist);
    }



}

/**************************************************************/
void UTIL::PreAmble() {

    stringstream outstream;
    int nSteps, nBodies, BodyDatum, IPKC, makeLog;
    REAL nRevs, maxT, rho;
    bool useRevs = false;


    cout << setfill('=') << setw(80) << "=" << endl;

    cout << "\t Engine Initialisation & Setup: " << endl;
    IO::VerboseMode = true;
    IO::FormattedQuery("Verbose Menu Mode?", "quiet=0, verbose=1, suggest quiet for script","Verbose Menu? Quiet=0, Verbose=1", outstream, IO::VerboseMode);
    if (!IO::VerboseMode)
        cout << endl;
    
    IO::FormattedQuery("Enter FMM maximum expansion order, pₘₐₓ ", "integer, minimum 3, suggest 5","FMM maximum expansion order", outstream, SYSTEM::MaxP);
    IO::FormattedQuery("Enter FMM smoothing parameter radius, δ", "real, suggest 0.0316","FMM smoothing parameter radius", outstream, SYSTEM::Del2);
    SYSTEM::Del2 *= SYSTEM::Del2;
    IO::FormattedQuery("Enter scale factor for model","real","Scale factor", outstream, SYSTEM::GambitScale);
    IO::FormattedQuery("Enter number of quadrature points to use during multipole expansions","integer, minimum 2, suggest 24","Number of quadrature points", outstream, UTIL::NumCellGaussPts);
    UTIL::lgwt(UTIL::NumCellGaussPts, UTIL::QuadPts, UTIL::QuadWts);
    //  This is to get the right multipliers for a cell, e.g. a range [-0.5,0.5] rather than [-1,1]
    for (int i = 0; i < UTIL::QuadPts.size(); ++i) {
        UTIL::QuadPts[i] = UTIL::QuadPts[i] / 2.0;
        UTIL::QuadWts[i] = UTIL::QuadWts[i] / 2.0;
    }
    IO::FormattedQuery("Enter M'₄ interpolation radius, measured in cells", "integer, suggest 2","M'4 interpolation radius", outstream, SYSTEM::M4Radius_in_cells);
    IO::FormattedQuery("Enter number of sub-time-steps to be approximated per FVM step, nₛₛ", "integer, suggest 5-10","Initial number of sub-time-steps", outstream, globalSystem->NumSubSteps);
    IO::FormattedQuery("Enter period of data output", "real","Data sample period", outstream, TIME_STEPPER::dt_out);
    IO::FormattedQuery("Enter fraction of timestep for protowake length, Δₛ", "real, suggest 0.3","Fraction of timestep for protowake length", outstream, globalSystem->DS);

        
    if (IO::VerboseMode){
    cout << "Select flux limiter:" << endl;
    cout << "         0) First order upwind " << endl;
    cout << "         1) minmod " << endl;
    cout << "         2) superbee " << endl;
    cout << "         3) van Leer " << endl;
    cout << "         4) Koren (3rd order)" << endl;
    cout << "         5) MUSCL (TVD LUD/Fromm) " << endl;
    cout << "         6) UMIST (TVD QUICK) " << endl;
    cout << "         7) CADA" << endl;
    }
    int FLIM;
    IO::FormattedQuery("         8) MC","","Flux limiter", outstream, FLIM);


    switch (FLIM) {
        case 0:
            FVMCell::limiter = FVMCell::O1UW;
            break;
        case 1:
            FVMCell::limiter = FVMCell::minmod;
            break;
        case 2:
            FVMCell::limiter = FVMCell::superbee;
            break;
        case 3:
            FVMCell::limiter = FVMCell::vanLeer;
            break;
        case 4:
            FVMCell::limiter = FVMCell::Koren;
            break;
        case 5:
            FVMCell::limiter = FVMCell::MUSCL;
            break;
        case 6:
            FVMCell::limiter = FVMCell::UMIST;
            break;
        case 7:
            FVMCell::limiter = FVMCell::CADA;
            break;
        case 8:
            FVMCell::limiter = FVMCell::MC;
            break;
        default:
            FVMCell::limiter = FVMCell::O1UW;
    }





    cout << setfill('=') << setw(80) << "=" << endl;
    cout << "\t Simulation Setup: " << endl;
    bool WaveModel = false;
    bool cnoidal = false;

    IO::FormattedQuery("Use wave model?", "no=0, yes=1","Wave model?", outstream, WaveModel);

    if (WaveModel) {
        IO::FormattedQuery("          Enter wave period in seconds","real","Period", outstream, WaveField::Period);
        IO::FormattedQuery("          Enter wave height in metres","real","Height", outstream, WaveField::Height);
        IO::FormattedQuery("          Enter water depth in metres","real","Water depth", outstream, WaveField::Depth);
        IO::FormattedQuery("Preferred wave theory","linear=0, cnoidal=1","Wave theory? linear=0, cnoidal=1", outstream, cnoidal);
    }

    IO::FormattedQuery("Enter number of bodies","integer","Number of bodies", outstream, nBodies);

    Array <string> infname(nBodies);

    IO::FormattedQuery("Enter body to use for specifying # revs [enter +ve integer], or use time [enter 0]","","Body for # revs [+ve], or use time [0]", outstream, BodyDatum);
    //cout << "Enter body to use for specifying # revs [enter +ve integer], or use time [enter 0]:\t";

    if (BodyDatum > 0) {
        useRevs = true;
        IO::FormattedQuery("Enter number of revs for body " + toString(BodyDatum) + " ","real","# revs for body " + toString(BodyDatum), outstream, nRevs);
    } else {
        IO::FormattedQuery("Enter simulation duration in seconds","real","Max sim time", outstream, maxT);
    }
    bool useTSR, useRadians;

    IO::FormattedQuery("Prefer to specify degrees or radians for attitudes?"," degrees=0, radians=1","Degrees[0] or radians[1] for attitudes?", outstream, useRadians);
    IO::FormattedQuery("Prefer to specify TSR plus axis or rates?", " rates=0, TSR=1","Rates or TSR? rates=0, TSR=1", outstream, useTSR);

    int defRates = 2;
    if (!useTSR) {
        IO::FormattedQuery("Prefer to specify rates RPM, Hz, or rad/s?", "RPM=0, Hz=1, rad/s=2","RPM, Hz or rad/s? RPM=0, Hz=1, rad/s=2", outstream, defRates);
    }
    IO::FormattedQuery("Enter freestream velocity U∞ V∞ W∞","3 x real","Freestream velocity", outstream, globalSystem->unscaledVinf);
    IO::FormattedQuery("Enter fluid density in kg/m³","real","Fluid density", outstream, globalSystem->Rho);
    IO::FormattedQuery("Use iterative pressure Kutta condition?","no=0, yes=1","Use IPKC? no=0, yes=1", outstream, IPKC);
    if (IPKC == 0)
        BODY::IPKC = false;
    else
        BODY::IPKC = true;



    IO::FormattedQuery("Make an input log file?","no=0, yes=1","Make log? no=0, yes=1", outstream, makeLog);

    Array <REAL> tsr(nBodies, 0.0), rads(nBodies, 0.0);
    BODY::NAMES = Array <string > (nBodies);
    BODY::CGS = Array <Vect3 > (nBodies, Vect3(0.));
    BODY::RATES = Array <Vect3 > (nBodies, Vect3(0.));
    BODY::VELOCITY = Array <Vect3 > (nBodies, Vect3(0.));
    BODY::ATTITUDE = Array <Vect3 > (nBodies, Vect3(0.));
    Array <Vect3> Disp(nBodies), Ax(nBodies);
    Array <bool> flip(nBodies, false);
    Array <int> plane(nBodies, 0);

    for (int i = 0; i < nBodies; ++i) {
        cout << setfill('=') << setw(80) << "=" << endl;
        cout << "\t Body " << i + 1 << " Setup:" << endl;

        IO::FormattedQuery("Enter input neutral file","string","Neutral file", outstream, infname[i]);
        
        REAL scale = 1.;
        IO::FormattedQuery("Enter geometry scale factor","real, normally suggest 1.0","Geometry scale", outstream, scale);
        
        IO::FormattedQuery("Mirror/flip input geometry?","no=0, yes=1","Mirror/flip geometry? no=0, yes=1", outstream, flip[i]);
       
        if (flip[i])
            IO::FormattedQuery("Mirror using yz, xz or xy plane?","yz=1, xz=2, xy=3","Flip plane yz=1, xz=2, xy=3", outstream, plane[i]);


         int tmp_int = 0;
        IO::FormattedQuery("Use either translation[0]/rotation[1]/both[2] for protowake?","int, suggest 2 only if confirmed working","Translation[0]/rotation[1]/both[2] for protowake", outstream, tmp_int);
        BODY::ProtowakeRotTrans.push_back(tmp_int);
        
        tmp_int = 0;
        IO::FormattedQuery("Number of panels to produce before conversion to vortex particles?","int, suggest 5-25","# panels before conversion to particles", outstream, tmp_int);
        BODY::nPanelsBeforeVortons.push_back(tmp_int);
        
        REAL tmp_real = 0;
                IO::FormattedQuery("Minimum distance before binning particle into octree?","real, suggest 1-5 chord lengths","Min distance before binning particle into tree", outstream, tmp_real);
        BODY::MinDistance.push_back(tmp_real);
        IO::FormattedQuery("Enter displacement of neutral file origin in global frame","3 x real","Displacement of origin [dx,dy,dz]", outstream, Disp[i]);
        Disp[i] = SYSTEM::GambitScale * Disp[i];
        IO::FormattedQuery("Enter body CG position (i.e. centre of rotation) xcg ycg zcg","3 x real","CG position [x,y,z]", outstream, BODY::CGS[i]);
        BODY::CGS[i] = SYSTEM::GambitScale * BODY::CGS[i];
        IO::FormattedQuery("Enter body CG translational velocity Vx Vy Vz", "3 x real","CG velocity [u,v,w]", outstream, BODY::VELOCITY[i]);

        if (useRadians)
            IO::FormattedQuery("Enter body attitude ψ, θ, and φ in radians", "3 x real","Body Euler angles (radians) [psi, theta, phi]", outstream, BODY::ATTITUDE[i]);
        else
        {
            IO::FormattedQuery("Enter body attitude ψ, θ, and φ in degrees", "3 x real","Body Euler angles (degrees) [psi, theta, phi]", outstream, BODY::ATTITUDE[i]);
            BODY::ATTITUDE[i] = deg2rad(BODY::ATTITUDE[i]);
        }

        if (useTSR) {
            IO::FormattedQuery("Enter tip-speed ratio, λ","real","TSR", outstream, tsr[i]);

            IO::FormattedQuery("Enter rotation axis unit vector", "3 x real","Rotation axis vector [vx,vy,vz]", outstream, Ax[i]);

            IO::FormattedQuery("Enter radius in metres","real","Radius", outstream, rads[i]);

            BODY::Radius = rads[0];
            BODY::RATES[i] = Ax[i]*((globalSystem->unscaledVinf - BODY::VELOCITY[i]).Mag() * tsr[i] / rads[i]);



        } else {
            Vect3 rt;


            if (defRates == 2) {
                IO::FormattedQuery("Enter body rotational velocity p q r in rads/s", "3 x real", "Body rates [p,q,r] in rads/s", outstream, rt);
                BODY::RATES[i] = rt;
            }

            if (defRates == 1) {
                IO::FormattedQuery("Enter body rotational velocity p q r in Hz", "3 x real", "Body rates [p,q,r] in Hz", outstream, rt);
                BODY::RATES[i] = rt * two_pi;
            }

            if (defRates == 0) {
                IO::FormattedQuery("Enter body rotational velocity p q r in RPM", "3 x real", "Body rates [p,q,r] in RPM", outstream, rt);
                BODY::RATES[i] = rt * two_pi / 60.;
            }

        }

        BODY::VELOCITY[i] = SYSTEM::GambitScale * BODY::VELOCITY[i];

        if (useRevs) {
            Vect3 Om(fabs(BODY::RATES[BodyDatum - 1].x), fabs(BODY::RATES[BodyDatum - 1].y), fabs(BODY::RATES[BodyDatum - 1].z));
            maxT = two_pi * nRevs / (max(Om));
        }
        string str_tmp = "tmp";
        BODY::ReadNeuGetBodies(infname[i], str_tmp, Disp[i], BODY::CGS[i], BODY::VELOCITY[i], BODY::ATTITUDE[i], BODY::RATES[i], flip[i], plane[i],scale);
    }

    cout << setfill('=') << setw(80) << "=" << endl;
    if (makeLog > 0) {

        cout << "\tInput file start:" << endl;
        cout << setfill('-') << setw(80) << "-" << endl;
        cout << outstream.str() << endl;
        cout << setfill('-') << setw(80) << "-" << endl;
        cout << "\tInput file end." << endl;
        cout << setfill('=') << setw(80) << "=" << endl;
    }
    UTIL::WriteMATLABString("Input", "Output.mat", outstream.str());
    globalSystem->InputStr = outstream.str();
    TIME_STEPPER::MaxTime = maxT;
    nSteps = globalSystem->NumSubSteps;



    cout << "\tSimulation Summary:" << endl;
    cout << "\tRuntime / Sample Number: \t" << maxT << "/" << nSteps << endl;
    cout << "\tInflow Velocity: \t\t" << globalSystem->unscaledVinf << endl << endl;


    for (int i = 0; i < BODY::Bodies.size(); ++i) {
        cout << setfill('.') << setw(80) << "." << endl;
        cout << "\tBody " << i + 1 << " Summary:" << endl;
        cout << "\tCG position: \t\t\t" << BODY::Bodies[i]->CG << endl;
        cout << "\tCG translational velocity: \t" << BODY::Bodies[i]->Velocity << endl;
        cout << "\tAttitude: \t\t\t" << BODY::Bodies[i]->EulerAngles << endl;
        cout << "\tBody rates: \t\t\t" << BODY::Bodies[i]->BodyRates << endl;
        cout << "\tNeutral file: \t\t\t" << infname[i] << endl;



    }
    //            cout << i+1 << " \t " << BODY::Bodies[i]->CG << " \t " << BODY::Bodies[i]->Velocity << " \t " << BODY::Bodies[i]->BodyAngles << " \t " << BODY::Bodies[i]->BodyRates << " \t " << BODY::Bodies[i]->Name << endl;


    cout << setfill('=') << setw(80) << "=" << endl;

    for (int i = 0; i < BODY::Bodies.size(); ++i)
        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j)
            BODY::Bodies[i]->Faces[j].GetNewGlobalPosition();

    BODY::PollFaces();
    globalSystem->dtInit = TIME_STEPPER::MaxTime / globalSystem->NumSubSteps;
    if (!SYSTEM::PanelOnly)
        globalSystem->dtInit /= 100.;

    BODY::SetUpProtoWakes(globalSystem->dtInit);
    Array < Array <REAL> > TRANS1, TRANS2, TRANS3, TmpPtsx, TmpPtsy, TmpPtsz, TmpPWPtsx, TmpPWPtsy, TmpPWPtsz, PhiPrev;
    for (int i = 0; i < BODY::Bodies.size(); ++i)
        for (int j = 0; j < BODY::Bodies[i]->ProtoWakes.size(); ++j)
            for (int k = 0; k < BODY::Bodies[i]->ProtoWakes[j].size(); ++k) {
                Array <REAL> t;
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C1.x);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C2.x);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C3.x);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C4.x);
                TmpPWPtsx.push_back(t);
                t.clear();
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C1.y);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C2.y);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C3.y);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C4.y);
                TmpPWPtsy.push_back(t);
                t.clear();
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C1.z);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C2.z);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C3.z);
                t.push_back(BODY::Bodies[i]->ProtoWakes[j][k].C4.z);
                TmpPWPtsz.push_back(t);
                t.clear();
            }

    WriteMATLABMatrix2D("ProtoWakePointsXInit", "Output.mat", TmpPWPtsx);
    WriteMATLABMatrix2D("ProtoWakePointsYInit", "Output.mat", TmpPWPtsy);
    WriteMATLABMatrix2D("ProtoWakePointsZInit", "Output.mat", TmpPWPtsz);

    BODY::PollFaces();

    BODY::SetUpInfluenceMatrices();

    if (WaveModel) {
        cout << setfill('=') << setw(80) << "=" << endl;
        if (cnoidal)
            WaveField::Cnoidal.WaveFieldCnoidal();
        else
            WaveField::Linear.WaveFieldLinear();
        cout << setfill('=') << setw(80) << "=" << endl;
    }

    IO::VerboseMode = true;
    if (SYSTEM::PanelOnly) {
        BODY::BodySubStep(maxT, nSteps);
    }
}

/**************************************************************/

void sheet(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL radius, REAL scale, REAL THETA) {
    int count = 0, nr = 1000, nt = 360;

    Array <Array <REAL> > R(3, Array <REAL > (3, 0.0));

    R[0][0] = cos(THETA);
    R[1][0] = 0.;
    R[2][0] = sin(THETA);
    R[0][1] = 0.;
    R[1][1] = 1.;
    R[2][1] = 0.;
    R[0][2] = -sin(THETA);
    R[1][2] = (THETA);
    R[2][2] = cos(THETA);


    REAL r = 0;
    X.allocate(nr * nt);
    Omega.allocate(nr * nt);
    for (int i = 0; i < nr; ++i) {
        REAL gamma = amplitude * sqrt(1 - ((r / radius) * (r / radius)));
        REAL theta = 0;
        for (int j = 0; j < nt; ++j) {
            Vect3 x = scale * Vect3(r * cos(theta), r * sin(theta), 0.0);

            Vect3 xt = scale * centre + Vect3(x.x * R[0][0] + x.z * R[2][0], x.y, x.x * R[0][2] + x.z * R[2][2]);
            X[count] = xt;

            Vect3 om = scale * scale * Vect3(-gamma * sin(theta), gamma * cos(theta), 0.0);
            Vect3 omt = Vect3(om.x * R[0][0] + om.z * R[2][0], om.y, om.x * R[0][2] + om.z * R[2][2]);
            Omega[count] = omt;
            theta += two_pi / nt;
            count++;
        }
        r += radius / nr;
    }

    Array <Vect3> Xn, Omn;
    Xn.push_back(X[0]);
    Omn.push_back(Omega[0]);
    for (int i = 1; i < X.size(); ++i) {
        X[i] = Vect3(round(X[i].x), round(X[i].y), round(X[i].z));
        bool put_in = true;
        for (int j = 0; j < Xn.size(); ++j)
            if (X[i].x == Xn[j].x && X[i].y == Xn[j].y && X[i].z == Xn[j].z) put_in = false;

        if (put_in) {
            Xn.push_back(X[i]);
            Omn.push_back(Omega[i]);
        }

    }

    X = Xn;
    Omega = Omn;
}
