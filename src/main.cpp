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
#define MAXMEM 90
//      maximum memory use before bombing out

using namespace std;
/**************************************************************/
void TestFMM(int argc, char *argv[]);
void WeeAmble();
void SolveMatfileVels(string fname, int pmax, REAL del2);
void sheet(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL radius, REAL scale, REAL THETA);
void PanelRecursiveDivide(PANEL &, Array <PANEL> &);
void PanelTriangleDivide(PANEL &, Array <PANEL> &);
REAL  PanelMaxTheta(PANEL &);

/**************************************************************/
class OutOfMemory {
};



class RunningStat {
public:

    RunningStat() : m_n(0), mmax(-1e32), mmin(1e32) {
    }

    void Clear() {
        m_n = 0;
    }

    void Push(REAL x) {
        m_n++;

        if (x > mmax) mmax = x;
        if (x < mmin) mmin = x;
        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1) {
            m_oldM = m_newM = x;
            m_oldS = 0.0;
        } else {
            m_newM = m_oldM + (x - m_oldM) / m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

            // set up for next iteration
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
    }

    int NumDataValues() const {
        return m_n;
    }

    REAL Mean() const {
        return (m_n > 0) ? m_newM : 0.0;
    }

    REAL Variance() const {
        return ( (m_n > 1) ? m_newS / (m_n - 1) : 0.0);
    }

    REAL StandardDeviation() const {
        return sqrt(Variance());
    }

    REAL Max() const {
        return mmax;
    }

    REAL Min() const {
        return mmin;
    }
private:
    int m_n;
    REAL mmax, mmin;
    REAL m_oldM, m_newM, m_oldS, m_newS;
};

/**************************************************************/
int main(int argc, char *argv[]) {
    system("clear");
    PANEL::Initialise();
    UTIL::GetCellPans();

    //  Get quadrature points and weights
    UTIL::NumCellGaussPts = 8;
    UTIL::lgwt(UTIL::NumCellGaussPts, UTIL::QuadPts, UTIL::QuadWts);
    for (int i = 0; i < UTIL::QuadPts.size(); ++i)
    {
        UTIL::QuadPts[i] = UTIL::QuadPts[i] / 2.0;
        UTIL::QuadWts[i] = UTIL::QuadWts[i] / 2.0;
    }

    //TestFMM(argc, argv);
   // WeeAmble();

//    return 0;
//    {
//        SYSTEM System(0);
//        globalSystem->Del2 = 0.001;//sqrt(2.)*0.5015;// * globalSystem->GambitScale*globalSystem->GambitScale;
//
//        WaveField::Depth = 25.0;
//        WaveField Field0, Field1;
//        Field0.SetPeakFreq(0.1/(2.*pi));
//        Field0.SetHSig(1.0);
//
//        Field1.SetModalFreq(0.1/(2.*pi)); 
//        Field1.SetHSig(1.0);
//        
//        Field0.getWaveFeld(&WaveField::JONSWAP,2.5,20.0,150);
//        Field1.getWaveFeld(&WaveField::Bretschneider,2.5,20.0,150);
//        Array <REAL> S0 = Field0.Spectrum(), S1 = Field1.Spectrum();
//        
////        for (int i = 0; i < S1.size(); ++i)
////            cout << S0[i] << " " << S1[i] << endl;
////        
////        cout << TIME_STEPPER::SimTime << endl;
//        
//        
//        
//      
//        int N = 1000;
//        Array <REAL> xs = UTIL::globalLinspace(-3,3,N);
//        
//        Array <Vect3> Vels(N), Velps(N), Velqs(N), VelGrad0(N), VelCubicGrad0(N), VelGrad1(N), VelCubicGrad1(N), VelGrad2(N), VelCubicGrad2(N);
//        for (int I = 0; I < N; ++I) {
//
//            Vect3 Target(xs[I], 0., 0.), Omega(1., 1.0, 1.0);
//
//            Vect3 VelP = UTIL::globalCubicDirectVel(Target,Omega);
//            Array <Vect3> GradsC(3, Vect3(0.0,0.0,0.0));
//            UTIL::globalCubicDirectVelGrads(Target, Omega, GradsC);
//            VelCubicGrad0[I] = GradsC[0];
//            VelCubicGrad1[I] = GradsC[1];
//            VelCubicGrad2[I] = GradsC[2];
//
//            int n = 100;
//            Array <REAL> pts = UTIL::globalLinspace(-0.5, 0.5, n);
//            Vect3 Vel(0., 0., 0.);
//            Array <Vect3> Grads(3, Vect3(0.0,0.0,0.0));
//            for (int i = 0; i < n; ++i)
//                for (int j = 0; j < n; ++j)
//                    for (int k = 0; k < n; ++k){
//                        Vel += UTIL::globalDirectVel(Target - Vect3(pts[i], pts[j], pts[k]), Omega / (REAL(n * n * n)));
//                        UTIL::globalDirectVelGrads(Target - Vect3(pts[i], pts[j], pts[k]), Omega / (REAL(n * n * n)), Grads);
//                    }
//            Vels[I] = Vel;
//            Velps[I] = VelP;
//            VelGrad0[I] = Grads[0];
//            VelGrad1[I] = Grads[1];
//            VelGrad2[I] = Grads[2];
//
//            Vect3 VelQ(0., 0., 0.);
//            for (int i = 0; i < UTIL::QuadPts.size(); ++i)
//                for (int j = 0; j < UTIL::QuadPts.size(); ++j)
//                    for (int k = 0; k < UTIL::QuadPts.size(); ++k)
//                        VelQ += UTIL::QuadWts[i] * UTIL::QuadWts[j] * UTIL::QuadWts[k] *
//                            UTIL::globalDirectVel(Target - Vect3(UTIL::QuadPts[i], UTIL::QuadPts[j], UTIL::QuadPts[k]), Omega);
//            Velqs[I] = VelQ;      // divided by 8 since the range of points is only -0.5->0.5 rather than -1->1
//        }
//        
//        UTIL::WriteMATLABMatrix1DVect3("Vels","Vels.mat",Vels);
//        UTIL::WriteMATLABMatrix1DVect3("VelPs","Vels.mat",Velps);
//        UTIL::WriteMATLABMatrix1DVect3("VelQs","Vels.mat",Velqs);
//        UTIL::WriteMATLABMatrix1D("xs","Vels.mat",xs);
//        UTIL::WriteMATLABMatrix1DVect3("Grads0","Vels.mat",VelGrad0);
//        UTIL::WriteMATLABMatrix1DVect3("Grads1","Vels.mat",VelGrad1);
//        UTIL::WriteMATLABMatrix1DVect3("Grads2","Vels.mat",VelGrad2);
//        UTIL::WriteMATLABMatrix1DVect3("GradsC0","Vels.mat",VelCubicGrad0);
//        UTIL::WriteMATLABMatrix1DVect3("GradsC1","Vels.mat",VelCubicGrad1);
//        UTIL::WriteMATLABMatrix1DVect3("GradsC2","Vels.mat",VelCubicGrad2);
//        
////        return 0;
//        
//        
//        
//    }
//    return 0;
//    
    
    
    /*
     *  This code has several "modes:"
     *  1) It can be used as a BEM simulator with or without a BEM wake. If a BEM 
     *     wake is included, it can take one of 2 forms:
     *     a) Panel representation using constant vortex panels;
     *     b) Vortex blob representation using point vortices.
     *     Either way, a direct NxN calculation must be done for the wake rollup,
     *     and this is not accelerated.
     *  2) It can be used with a BEM source for a FVM vortex code. This is
     *     accelerated using a FMM method for the Biot-Savart calculation.
     *  3) It can be used to model vortex flows in the absence of a body, e.g.
     *     a) The normal or 
     *     b) oblique collision of two vortex rings.
     *  4) It can be used to calculate the velocity due to a vorticity field, as
     *     provided at runtime.
     */

//    {
//        SYSTEM System(0);
//        globalOctree->Root->SetUpISBIndices() ;
//    }
//    return 0;

//
//    TestFMM(argc, argv);
//
//    return 0;
//


//
//    cout << "Enter filename:" << endl;
//    string fname;
//
//
//     SolveMatfileVels(fname, 8, 0.25);
// 
// 
// 
//     return 0;


    //    WeeAmble();
    //    return 1;

    SYSTEM System(0);

    //  Some default values
    globalSystem->GambitScale = 1.0;
    globalSystem->MaxP = 5;
    globalSystem->Del2 = 0.001;// * globalSystem->GambitScale*globalSystem->GambitScale;
    globalSystem->DS = .3;
    globalSystem->dtInit = 0.05;
    globalSystem->h = 2;
    globalSystem->unscaledVinf = Vect3(0.0);
    globalSystem->NumSubSteps = 0;


    UTIL::cpu_t = ticks();

    //UTIL::PreAmble();


        TIME_STEPPER::MaxTime = 100.0;
        globalSystem->useBodies = false;

        globalSystem->NumTransVars = 2;


    globalSystem->Initialise();
    globalSystem->VortonsXs.clear();
    globalSystem->VortonOmegas.clear();
        
    
    
    {
        
        
        REAL Radius = 20.0;
        REAL gamma = 5.0;
        int n = 10000;
        Array <REAL> thetas = globalLinspace(0.0,2*pi,n);
        
        Vect3 Centre(0.0,0.0,5);
        
        REAL arc_length = Radius * (thetas[1]-thetas[0]);
        
        Array <Vect3> Xs(n), Oms(n);
        Array <int> IDs(n,0);
        for (int i = 0; i < n; ++i)
        {
            Xs[i] = Centre + Radius * Vect3(cos(thetas[i]), sin(thetas[i]),0.0);
            Oms[i] = gamma * arc_length * Radius * Vect3(-sin(thetas[i]), cos(thetas[i]), 0.0);
        }
        
        
        
        globalSystem->AddVortonsToTree(Xs,Oms, IDs);
        
         for (int i = 0; i < n; ++i)
        {
            Xs[i].z *= -1;
            Oms[i] = -Oms[i];
            IDs[i] = 1;
        }
        
        
        globalSystem->AddVortonsToTree(Xs,Oms, IDs);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
    
           
//        REAL radius = 1.0, amplitude = 1.0;
//        Vect3 centre(1.,0,0);
//        Array <Vect3> X, Omega;
//        
//    
//        REAL THETA = 3*0.523598775598299;
//    
//        sheet(centre, X, Omega, -amplitude, radius, globalSystem->GambitScale, THETA);
//        
//        cout << X.size() << " " << Omega.size() << endl;
//        
//        Array <int> IDs(X.size(),0);
//        
//        globalSystem->AddVortonsToTree(X,Omega, IDs);
//    
//        centre = Vect3(-1.,0,0);
//     
//        sheet(centre, X, Omega, -amplitude, radius, globalSystem->GambitScale, -THETA);
//        IDs = 1;
//        cout << X.size() << " " << Omega.size() << endl;    
//        globalSystem->AddVortonsToTree(X,Omega, IDs);
//    
    #ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "globalSystem->MaxP set to " << globalSystem->MaxP << "; dtInit " << globalSystem->dtInit << endl;
    #endif
    
        cout << "Number of cells: " << FVMCell::NumCells << endl;
    globalSystem->TimeStep();


    UTIL::PostAmble(string("Output.mat"));


#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000 << " seconds" << endl;
#endif
}




/**************************************************************/
void UTIL::PostAmble(string fname) {

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


        UTIL::WriteMATLABMatrix1DVect3("EulerAngles" + UTIL::toString(i), fname, BODY::Bodies[i]->EulerAngles);
        UTIL::WriteMATLABMatrix1DVect3("TRANS1_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[0]);
        UTIL::WriteMATLABMatrix1DVect3("TRANS2_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[1]);
        UTIL::WriteMATLABMatrix1DVect3("TRANS3_" + UTIL::toString(i), fname, BODY::Bodies[i]->TRANS[2]);


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
    //WriteMATLABMatrix2D("Amatrix", fname, BODY::Bodies[0]->localA);
    //WriteMATLABMatrix2D("Bmatrix", fname, BODY::Bodies[0]->localB);
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
    Array <Vect3> AllBodyPointsTransformed(BODY::AllBodyPoints.size(), Vect3(0.0));

    //        for (int i = 0; i < BODY::AllBodyPoints.size(); ++i)
    //            AllBodyPointsTransformed[i] = VectMultMatrix(BODY::AllBodyFaces->, vCorner_g);

    UTIL::WriteMATLABMatrix1DVect3("AllBodyPoints", fname, BODY::AllBodyPoints);
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
    REAL nRevs, maxT, uinf, vinf, winf, rho;
    bool useRevs = false;


    cout << setfill('=') << setw(80) << "=" << endl;
    cout << "\t Simulation Setup: " << endl;
    cout << "Enter number of bodies [int]:";
    cin >> nBodies;
    outstream << nBodies << endl;

    Array <string> infname(nBodies);

    cout << "Enter body to use for specifying # revs [enter +ve integer], or use time [enter 0]:" << endl;
    cin >> BodyDatum;
    outstream << BodyDatum << endl;

    if (BodyDatum > 0) {
        useRevs = true;
        cout << "Enter number of revs for body " << BodyDatum << " [real]:" << endl;
        cin >> nRevs;
        outstream << nRevs << endl;
    } else {
        cout << "Enter simulation duration in seconds [real]:" << endl;
        cin >> maxT;
        outstream << maxT << endl;
    }
    bool useTSR, useRadians;

    cout << "Prefer to specify degrees or radians for attitudes [degrees=1, radians=0]?:" << endl;
    cin >> useRadians;
    outstream << useRadians << endl;
    cout << "Prefer to specify TSR plus axis or rates [TSR=1, rates=0]:" << endl;
    cin >> useTSR;
    outstream << useTSR << endl;

    int defRates = 2;
    if (!useTSR) {
        cout << "Prefer to specify rates RPM, Hz, or rad/s [RPM=0, Hz=1, rad/s=2]:" << endl;
        cin >> defRates;
        outstream << defRates << endl;
    }





    cout << "Enter number of timestep samples [integer]:" << endl;
    cin >> nSteps;
    outstream << nSteps << endl;

    cout << "Enter freestream velocity Uinf Vinf Winf as [3 x real]:" << endl;

    cin >> globalSystem->unscaledVinf.x >> globalSystem->unscaledVinf.y >> globalSystem->unscaledVinf.z;
    outstream << globalSystem->unscaledVinf << endl;
    cout << "Enter fluid density in kg/m3 [real]:" << endl;
    cin >> rho;
    outstream << BODY::RHO << endl;
    cout << "Use iterative pressure kutta condition? [bool 1/0]:" << endl;
    cin >> IPKC;
    outstream << IPKC << endl;
    if (IPKC == 0)
        BODY::IPKC = false;
    else
        BODY::IPKC = true;



    cout << "Make an input log file? [bool 1/0]:" << endl;
    cin >> makeLog;
    outstream << makeLog << endl;

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

        cout << "Enter input neutral file [string]:" << endl;
        cin >> infname[i];
        outstream << infname[i] << endl;


        cout << "Mirror/flip input geometry? [bool 1/0]" << endl;
        cin >> flip[i];
        outstream << flip[i] << endl;

        if (flip[i]) {
            cout << "Mirror using yz, xz or xy plane? [yz=1, xz=2, xy=3]" << endl;
            cin >> plane[i];
            outstream << plane[i] << endl;

        }

        cout << "Enter displacement of neutral file origin in global frame as [3 x real]:" << endl;
        cin >> Disp[i].x >> Disp[i].y >> Disp[i].z;
        outstream << Disp[i] << endl;
        Disp[i] = globalSystem->GambitScale * Disp[i];
        cout << "Enter body CG position (i.e. centre of rotation) xcg ycg zcg as [3 x real]:" << endl;
        cin >> BODY::CGS[i].x >> BODY::CGS[i].y >> BODY::CGS[i].z;
        outstream << BODY::CGS[i] << endl;
        BODY::CGS[i] = globalSystem->GambitScale * BODY::CGS[i];
        cout << "Enter body CG translational velocity Vx Vy Vz as [3 x real]:" << endl;
        cin >> BODY::VELOCITY[i].x >> BODY::VELOCITY[i].y >> BODY::VELOCITY[i].z;

        outstream << BODY::VELOCITY[i] << endl;


        cout << "Enter body attitude psi theta phi in degrees as [3 x real]:" << endl;
        cin >> BODY::ATTITUDE[i].x >> BODY::ATTITUDE[i].y >> BODY::ATTITUDE[i].z;
        outstream << BODY::ATTITUDE[i] << endl;

        if (useTSR) {
            cout << "Enter TSR (lambda) [real]:" << endl;
            cin >> tsr[i];
            outstream << tsr[i] << endl;

            cout << "Enter rotation axis as [3 x real]:" << endl;
            cin >> Ax[i].x >> Ax[i].y >> Ax[i].z;
            outstream << Ax[i] << endl;

            cout << "Enter radius as [real]:" << endl;
            cin >> rads[i];
            outstream << rads[i] << endl;

            BODY::Radius = rads[0];
            BODY::RATES[i] = Ax[i]*((globalSystem->unscaledVinf - BODY::VELOCITY[i]).Mag() * tsr[i] / rads[i]);



        } else {
            cout << "Enter body rotational velocity p q r as [3 x real]:" << endl;
            Vect3 rt;
            cin >> rt.x >> rt.y >> rt.z;
            outstream << rt << endl;

            if (defRates == 2)
                BODY::RATES[i] = rt;

            if (defRates == 1)
                BODY::RATES[i] = rt * 2 * pi;

            if (defRates == 0)
                BODY::RATES[i] = rt * 2 * pi / 60;
        }

        BODY::VELOCITY[i] = globalSystem->GambitScale * BODY::VELOCITY[i];

        if (useRevs) {
            Vect3 Om(fabs(BODY::RATES[BodyDatum - 1].x), fabs(BODY::RATES[BodyDatum - 1].y), fabs(BODY::RATES[BodyDatum - 1].z));
            maxT = 2 * pi * nRevs / (max(Om));
        }
        BODY::ReadNeuGetBodies(infname[i], "tmp", Disp[i], BODY::CGS[i], BODY::VELOCITY[i], BODY::ATTITUDE[i], BODY::RATES[i], flip[i], plane[i]);
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
    globalSystem->NumSubSteps = nSteps;



    cout << "\tSimulation Summary:" << endl;
    cout << "\tRuntime / Sample Number: \t" << maxT << "/" << nSteps << endl;
    cout << "\tInflow Velocity: \t\t" << globalSystem->unscaledVinf << endl << endl;


    for (int i = 0; i < BODY::Bodies.size(); ++i) {
        cout << setfill('.') << setw(80) << "." << endl;
        cout << "\tBody " << i + 1 << " Summary:" << endl;
        cout << "\tCG position: \t\t\t" << BODY::Bodies[i]->CG << endl;
        cout << "\tCG translational velocity: \t" << BODY::Bodies[i]->Velocity << endl;
        cout << "\tAttitude: \t\t\t" << BODY::Bodies[i]->BodyAngles << endl;
        cout << "\tBody rates: \t\t\t" << BODY::Bodies[i]->BodyRates << endl;
        cout << "\tNeutral file: \t\t\t" << infname[i] << endl;



    }
    //            cout << i+1 << " \t " << BODY::Bodies[i]->CG << " \t " << BODY::Bodies[i]->Velocity << " \t " << BODY::Bodies[i]->BodyAngles << " \t " << BODY::Bodies[i]->BodyRates << " \t " << BODY::Bodies[i]->Name << endl;


    cout << setfill('=') << setw(80) << "=" << endl;

    for (int i = 0; i < BODY::Bodies.size(); ++i)
        for (int j = 0; j < BODY::Bodies[i]->Faces.size(); ++j)
            BODY::Bodies[i]->Faces[j].GetNewGlobalPosition();

    BODY::PollFaces();

    BODY::SetUpProtoWakes(0.01);
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


    //    BODY::BodySubStep(maxT, nSteps);

}
/**************************************************************/
void UTIL::GetCellPans() {
    Array <Vect3> Verts(8);
    Verts[0] = 0.5 * Vect3(-1, -1, -1);
    Verts[1] = 0.5 * Vect3(-1, -1, 1);
    Verts[2] = 0.5 * Vect3(1, -1, 1);
    Verts[3] = 0.5 * Vect3(1, -1, -1);
    Verts[4] = 0.5 * Vect3(-1, 1, -1);
    Verts[5] = 0.5 * Vect3(-1, 1, 1);
    Verts[6] = 0.5 * Vect3(1, 1, 1);
    Verts[7] = 0.5 * Vect3(1, 1, -1);
    

    UTIL::Pans.push_back(new PANEL(Verts[0], Verts[1], Verts[2], Verts[3]));
    UTIL::Pans.push_back(new PANEL(Verts[3], Verts[2], Verts[6], Verts[7]));
    UTIL::Pans.push_back(new PANEL(Verts[7], Verts[6], Verts[5], Verts[4]));
    UTIL::Pans.push_back(new PANEL(Verts[4], Verts[5], Verts[1], Verts[0]));
    UTIL::Pans.push_back(new PANEL(Verts[1], Verts[5], Verts[6], Verts[2]));
    UTIL::Pans.push_back(new PANEL(Verts[4], Verts[0], Verts[3], Verts[7]));

    for (int i = 0; i < 6; ++i) {
        UTIL::Pans[i]->GetNormal();
        UTIL::Pans[i]->Sigma = UTIL::Pans[i]->Gamma = UTIL::Pans[i]->Mu = 1.0;
    }
}
/**************************************************************/
Vect3 UTIL::globalDirectVel(Vect3 diff, Vect3 omega) {

    REAL mult, nrm;
    nrm = sqrt(globalSystem->Del2 + diff.Dot(diff));
    mult = -1. / (four_pi * nrm * nrm * nrm);
    return mult * diff.Cross(omega);
}
/**************************************************************/
Vect3 UTIL::globalCubicDirectVel(Vect3 diff, Vect3 omega) {
   if (diff.Dot(diff) > 4) {
        Vect3 VelQ(0., 0., 0.);
        for (int i = 0; i < UTIL::QuadPts.size(); ++i)
            for (int j = 0; j < UTIL::QuadPts.size(); ++j)
                for (int k = 0; k < UTIL::QuadPts.size(); ++k)
                    VelQ += UTIL::QuadWts[i] * UTIL::QuadWts[j] * UTIL::QuadWts[k] *
                        UTIL::globalDirectVel(diff - Vect3(UTIL::QuadPts[i], UTIL::QuadPts[j], UTIL::QuadPts[k]), omega);
        return VelQ; 
    } else {
        Vect3 VelP(0., 0., 0.);
        for (int i = 0; i < 6; ++i) {
            REAL Phi = UTIL::Pans[i]->HyperboloidSourcePhi(diff);
            VelP -= UTIL::Pans[i]->TRANS[2].Cross(omega) * Phi / 2.0;
        }
        return VelP;
    }
}

/**************************************************************/
void UTIL::globalDirectVelGrads(Vect3 diff, Vect3 omega, Array <Vect3> &Grads) {
    //  Diff here is Xtarget - Xsource (hence -ve (-=) below)
    REAL R, R3, R5;
    
    R = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z + globalSystem->Del2);
    R3 = R*R*R;
    R5 = R3*R*R;
    
    Grads[0].x  -= -(3.0*diff.x*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5); 
    Grads[0].y  -= (3.0*diff.x*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5) - omega.z/(four_pi*R3);
    Grads[0].z  -= omega.y/(four_pi*R3) - (3.0*diff.x*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    Grads[1].x  -= omega.z/(four_pi*R3) - (3.0*diff.y*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5);
    Grads[1].y  -= (3.0*diff.y*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5);
    Grads[1].z  -= -omega.x/(four_pi*R3) - (3.0*diff.y*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    Grads[2].x  -= -omega.y/(four_pi*R3) - (3.0*diff.z*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5);
    Grads[2].y  -= omega.x/(four_pi*R3) + (3.0*diff.z*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5);
    Grads[2].z  -= -(3.0*diff.z*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    

}
/**************************************************************/
void UTIL::globalCubicDirectVelGrads(Vect3 diff, Vect3 Omega, Array <Vect3> &Grads) {
    //  Diff here is Xtarget - Xsource (hence -ve (-=) below)


    if (diff.Dot(diff) > 4) {
        for (int i = 0; i < UTIL::QuadPts.size(); ++i)
            for (int j = 0; j < UTIL::QuadPts.size(); ++j)
                for (int k = 0; k < UTIL::QuadPts.size(); ++k)
                    UTIL::globalDirectVelGrads(diff - Vect3(UTIL::QuadPts[i], UTIL::QuadPts[j], UTIL::QuadPts[k]), UTIL::QuadWts[i] * UTIL::QuadWts[j] * UTIL::QuadWts[k] * Omega, Grads);
    } else {
        for (int i = 0; i < 6; ++i) {
            Vect3 V = UTIL::Pans[i]->SourcePanelVelocity(diff) * two_pi / four_pi;
            Vect3 C = UTIL::Pans[i]->TRANS[2].Cross(Omega);
            Grads[0] += V.x*C;
            Grads[1] += V.y*C;  // NB this is -ve because the return value for V.y in SourceVel is -ve (why?)
            Grads[2] += V.z*C;
        }
    }
}
/**************************************************************/
void SolveMatfileVels(string fname, int pmax, REAL del2) {

    system("clear");

    cout << "TUI driven file based specification of vortex points and target points." << endl;
    Array <Vect3> Posns, Omegas;



    Array <REAL> data;

    Array <int> dims;
    cout << "Enter name of .mat file containing [Nx3] list of vortex locations" << endl;
    getline(cin, fname);
    string varname = "Posns";
    int err = UTIL::readmat(fname, varname, data, dims, true);
    int numel = dims[0], count = 0;
    Posns.allocate(numel);
    Vect3 Mins(1e32), Maxs(-1e32);
    for (int i = 0; i < dims[0]; i++) {
        Posns[i].x = data[dims[0] * 0 + i] - 0.5;
        Posns[i].y = data[dims[0] * 1 + i] - 0.5;
        Posns[i].z = data[dims[0] * 2 + i] - 0.5;
        Mins = min(Mins, Posns[i]);
        Maxs = max(Maxs, Posns[i]);
        //cout << Posns[i] << endl;
    }
    cout << "Domain bounds: [" << Mins.x << " " << Maxs.x << "][" << Mins.y << " " << Maxs.y << "][" << Mins.z << " " << Maxs.z << "]" << endl;
    cout << "Enter name of .mat file containing [" << dims[0] << "x3] list of vortex strengths" << endl;
    getline(cin, fname);
    varname = "Omegas";
    err = UTIL::readmat(fname, varname, data, dims, true);
    numel = dims[0];
    count = 0;
    Omegas.allocate(numel);


    for (int i = 0; i < dims[0]; i++) {
        Omegas[i].x = data[dims[0] * 0 + i];
        Omegas[i].y = data[dims[0] * 1 + i];
        Omegas[i].z = data[dims[0] * 2 + i];
        //cout << Omegas[i] << endl;
    }

    data.clear();

    cout << "done reading. Finding active cells";

    Array <bool> isActive(Omegas.size(), false);
    unsigned long int activeCount = 0;

    for (int i = 0; i < Posns.size(); ++i) {
        if (Omegas[i].Mag() > 1e-6) {
            isActive[i] = true;
            activeCount++;
        } else {
            isActive[i] = false;
        }
    }



    Array <Vect3> ActiveP(activeCount), ActiveO(activeCount);
    count = 0;
    for (int i = 0; i < Posns.size(); ++i)
        if (isActive[i]) {
            ActiveP[count] = Posns[i];
            ActiveO[count] = Omegas[i];
            count++;
        }
    //        else
    //        {
    //            Vect3 Pos = globalSystem->GambitScale*Posns[i];
    //            Vect3 ParentPos = floor(Pos) + 1.0;
    //        }

    //    for (int i = 0; i < Posns.size(); ++i)

    cout << " done. " << activeCount << " active cells." << endl;




//        cout << "Not using tree Doing direct calc. This might take some time.." << endl;
//        
//        
//        
//        Array < Array < Vect3 > > PosI, VelI;
//
//    PosI = UTIL::zerosv(200, 50);
//    VelI = PosI;
//
//    Vect3 Pos(-10, -100, -50);
//    for (int i = 0; i < 200; ++i)
//        for (int j = 0; j < 50; ++j)
//            PosI[i][j] = (Pos + Vect3(0.0, REAL(i), REAL(j)))/7.5;
//
//
//    count = 0;
//    Array <Vect3> Vels(numel);
//#pragma omp parallel for
//    for (int i = 0; i < PosI.size(); ++i)
//        for (int j = 0; j < PosI[i].size(); ++j) 
//            for (int k = 0; k < activeCount; ++k) 
//                VelI[i][j] += UTIL::globalDirectVel(ActiveP[k] - PosI[i][j], ActiveO[k], 0.01);
//
//
//
//        UTIL::WriteMATLABMatrix2DVect3(string("Posns"), string("OutArrayA.mat"), PosI);
//        UTIL::WriteMATLABMatrix2DVect3(string("Vels"), string("OutArrayA.mat"), VelI);
//
//        return;
//    
//

        SYSTEM System(0);

    cout << "Enter GambitScale" << endl;
    cin >> globalSystem->GambitScale;
    cout << "Enter MaxP" << endl;
    cin >> globalSystem->MaxP;
    cout << "Enter Del2" << endl;
    cin >> globalSystem->Del2;


    cout << "Beginning insertion into tree" << endl;

    globalSystem->Initialise();

    REAL t0 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;

    for (int i = 0; i < Posns.size(); ++i) {
//        if (isActive[i]) {
            OctreeCapsule C(globalSystem->GambitScale * Posns[i], globalSystem->GambitScale * globalSystem->GambitScale * Omegas[i], true);
            C.AssociatedBody = 0;
            globalOctree->Root->EvalCapsule(C);




            if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {

                string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
                REAL MEM_PERCENT, temp;
                stringstream psdata;
                psdata << top_data;
                psdata >> temp >> MEM_PERCENT;
                cout << FVMCell::NumCells << " " << Node::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
                if (MEM_PERCENT > 90) {
                    cout << setfill('!') << setw(80) << "!" << endl;

                    cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
                    cout << "Memory used: " << MEM_PERCENT << "%" << endl;

                    cout << setfill('!') << setw(80) << "!" << endl;
                    throw OutOfMemory();
                }
            }
//        }
    }

    globalTimeStepper->PruneNow = false;
    globalOctree->Reset();
    Mins = 1e32;
    Maxs = -1e32;

    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
        Mins = min(globalOctree->AllCells[i]->Position, Mins);
        Maxs = max(globalOctree->AllCells[i]->Position, Maxs);
    }
    cout << "Scaled domain bounds: [" << Mins.x << " " << Maxs.x << "][" << Mins.y << " " << Maxs.y << "][" << Mins.z << " " << Maxs.z << "]" << endl;
    cout << "Calculating velocities" << endl;
    globalOctree->ResetAllVelsAndFields();
    globalOctree->GetVels();
    REAL t1 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    cout << "Done. Time elapsed: " << t1 - t0 << endl;
//
    int n = globalOctree->AllCells.size();
    Posns.allocate(n);
    Omegas.allocate(n);
    Array <Vect3> FMMVels(n);
    Array <Vect3> DirVels(n,Vect3(0.0,0.0,0.0));





    for (int i = 0; i < n; ++i) {
        Posns[i] = (globalOctree->AllCells[i]->Position);
        Omegas[i] = (globalOctree->AllCells[i]->Omega);
        FMMVels[i] = (globalOctree->AllCells[i]->Velocity);
    }
    
    
//    
#pragma omp parallel for
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j)
            DirVels[i] += UTIL::globalDirectVel(Posns[j] - Posns[i], Omegas[j]) ;
        
        //cout << DirVels[i] << " " << FMMVels[i] << (DirVels[i] - FMMVels[i]).Mag() << endl;
    }
    
    

    cout << "Clearing tree ";
    globalOctree->ClearNodes();
    cout << "done." << endl;
    string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
    REAL MEM_PERCENT, temp;
    stringstream psdata;
    psdata << top_data;
    psdata >> temp >> MEM_PERCENT;
    cout << FVMCell::NumCells << " " << Node::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
    cin.ignore();
    cout << "Enter name of .mat file to write list of [" << numel << "x9] vortex positions, strengths and velocities" << endl;
    getline(cin, fname);
    string vname = "DomainData";

    UTIL::WriteMATLABMatrix1DVect3(string("Posns"), fname, Posns);
    UTIL::WriteMATLABMatrix1DVect3(string("Omegas"), fname, Omegas);
    UTIL::WriteMATLABMatrix1DVect3(string("Vels"), fname, FMMVels);
    UTIL::WriteMATLABMatrix1DVect3(string("DVels"), fname, DirVels);

    return;
}
//
//
//    
//
//
//
//
////    Array <REAL> L2;
////
////
////    REAL t0 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
////    globalOctree->Reset();
////    globalOctree->InitVelsGetLaplacian();
////    globalOctree->GetVels();
////    REAL t1 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
////    cout << "Done. Time elapsed: " << t1 - t0 << endl << "Performing Direct Calculation" << endl;
////    Posns.clear();
////    Omegas.clear();
////    int n = globalOctree->AllCells.size();
////    Posns.allocate(n);
////    Omegas.allocate(n);
////    Array <Vect3> FMMVels(n);
////
////
////
////
////
////    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
////        Posns[i] = (globalOctree->AllCells[i]->Position);
////        Omegas[i] = (globalOctree->AllCells[i]->Omega);
////        FMMVels[i] = (globalOctree->AllCells[i]->Velocity);
////    }
//
//
////
////    Array <int> Indices;
////    REAL Mult = 1.0;
////    if (Posns.size() > 10000) {
////        Mult = Posns.size() / 10000.0;
////        Indices = Array <int> (10000, 0);
////        int count = 0;
////        //while (RandIndices.size() < 10000)
////        while (count < 10000) {
////            int randint = int (Posns.size()*(REAL(rand()) / RAND_MAX));
////            bool isin = false;
////            for (int i = 0; i < 10000; ++i)
////                if (Indices[i] == randint)
////                    isin = true;
////
////            if (!isin) {
////                Indices[count] = randint;
////                count++;
////            }
////        }
////
////    } else {
////        Indices = Array <int> (Posns.size(), 0);
////        for (int i = 0; i < Posns.size(); ++i)
////            Indices[i] = i;
////    }
////
////
////    Array <Vect3> DirectVels(Indices.size());
////
////    REAL t2 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
////#pragma omp parallel for
////    for (int i = 0; i < Indices.size(); ++i) {
////        Vect3 V(0, 0, 0);
////        for (int j = 0; j < Posns.size(); ++j) {
////            Vect3 D = Posns[j] - Posns[Indices[i]];
////            V += globalDirectVel(D, Omegas[j]);
////
////        }
////        DirectVels[i] = (V);
////    }
////    REAL t3 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
////    cout << "Done. Time elapsed: " << t3 - t2 << endl << "Calculating error L2 norm";
////    REAL l2 = 0, Vmean = 0;
////
////
////    Array <REAL> AbsErrs(Indices.size()), RelErrs(Indices.size());
////
////    REAL AbsErrsSum2 = 0., RelErrsSum2 = 0., AbsErrsSum = 0., RelErrsSum = 0., AbsErrsMax = 0., RelErrsMax = 0.;
////
////    RunningStat AbsErrsStats, RelErrsStats;
////
////
////    for (int i = 0; i < Indices.size(); ++i) {
////        cout << "Cell " << i << "\tError L2 Norm " << (DirectVels[i] - FMMVels[Indices[i]]).Mag() << " \t " << DirectVels[i] << endl << "\t\t\t\t\t\t " << FMMVels[Indices[i]] << endl;
////        l2 += (DirectVels[i] - FMMVels[Indices[i]]).Mag()*(DirectVels[i] - FMMVels[Indices[i]]).Mag();
////        Vmean += DirectVels[i].Mag();
////
////        AbsErrs[i] = (DirectVels[i] - FMMVels[Indices[i]]).Mag();
////        RelErrs[i] = AbsErrs[i] / DirectVels[i].Mag();
////
////        AbsErrsStats.Push(AbsErrs[i]);
////        RelErrsStats.Push(RelErrs[i]);
////
////        AbsErrsSum += AbsErrs[i];
////        RelErrsSum += RelErrs[i];
////        AbsErrsSum2 += AbsErrs[i] * AbsErrs[i];
////        RelErrsSum2 += RelErrs[i] * RelErrs[i];
////
////        AbsErrsMax = max(AbsErrsMax, AbsErrs[i]);
////        RelErrsMax = max(RelErrsMax, RelErrs[i]);
////
////
////    }
////
////
////    REAL AbsErrsStd = sqrt(AbsErrsSum2 / (Posns.size() - 1));
////    REAL RelErrsStd = sqrt(RelErrsSum2 / (Posns.size() - 1));
////    REAL AbsErrsMean = AbsErrsSum / Posns.size();
////    REAL RelErrsMean = RelErrsSum / Posns.size();
////
////    int TreeSize = OCTREE_SIZE;
////    int TreeLevs = OCTREE_LEVS;
////    Vmean = Vmean / n;
////    cout << "Done." << endl << "Standard deviation: " << sqrt(l2 / Posns.size()) << endl;
////
////    ofstream myfile;
////    myfile.open("output.dat", ios::out | ios::app);
////    myfile << globalSystem->MaxP << " " << FVMCell::NumCells << " " << FVMCell::NumNodes << " " << globalSystem->NumThreads << " " << (t1 - t0)*1000 << " " << Mult * (t3 - t2)*1000 << " " << TreeLevs << " " << TreeSize << " " << AbsErrsStats.Max() << " " << AbsErrsStats.Mean() << " " << AbsErrsStats.StandardDeviation() << " " << RelErrsStats.Max() << " " << RelErrsStats.Mean() << " " << RelErrsStats.StandardDeviation() << endl;
////    //    myfile << n << " " << globalSystem->MaxP << " " << Vmean << " " << sqrt(l2 / Posns.size()) << " " << sqrt(l2 / Posns.size())/Vmean << " " << (t1 - t0)*1000 << " " << (t3 - t2)*1000 << endl;
////    myfile.close();
////
////#ifndef use_NCURSES
////    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000 << " seconds" << endl;
////    cout << "Output written to output.dat:" << endl;
////    cout << "p_max FVMCell::NumCells FVMCell::NumNodes nthreads t_fmm (ms)  t_dir (ms) TreeLevs TreeSize AbsErrsMax AbsErrsMean AbsErrsStd RelErrsMax RelErrsMean RelErrsStd sparse/dense sphere/cube;" << endl;
////#endif
//
//
//}

/**************************************************************/
void TestFMM(int argc, char *argv[]) {

    system("clear");

    if (argc < 5) {
        cout << "FMM Test Mode. Expects argument: Nvortices maxP numThreads sparse/dense cube/sphere" << endl;
        //cout << "Returns velocities at these points as calculated via FMM and directly." << endl;
        return;
    }

    string sparsity = argv[4];
    string sparse("sparse"), dense("dense");

    string geometry = argv[5];
    string cube("cube"), sphere("sphere");
#ifdef _OPENMP
    omp_set_num_threads(atoi(argv[3]));
#endif

    SYSTEM System(0);

    //  Some default values
    globalSystem->GambitScale = 1;
    globalSystem->MaxP = atoi(argv[2]);
    globalSystem->Del2 = 0.001;


    //    globalIO->read_input(dir2 + argv[1]);
    //
    //    globalIO->PrepOutputDir();
    //    globalIO->WriteBinary();
#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "globalSystem->MaxP set to " << globalSystem->MaxP << "; dtInit " << globalSystem->dtInit << endl;
#endif

    globalSystem->Initialise();



    Array <Vect3> Posns, Omegas;
    Array <REAL> L2;



    srand((unsigned) time(NULL));


    int n = atoi(argv[1]);

    Posns.allocate(n);
    Omegas.allocate(n);
    cout << "Generating " << n << " points Preparing ";

    if (geometry.compare(cube) == 0) {
        REAL rho = 4095 * 2;
        cout << "cube with a ";
        if (sparsity.compare(dense) == 0) {
            cout << "dense cell arrangement" << endl;
            
            //  Get the side length of the cube
            REAL l = pow(REAL(n),1./3.);
            int nl = ceil(l);
            cout << "Number of cells in cube " << nl << "^3 = " << nl*nl*nl << endl;
            
            if (round((REAL)nl/2.0)!=(REAL)nl/2.0){
                cout << "Odd number - increasing by 1" << endl;
                nl+=1;
            }
            else
                cout << "Even number" << endl;
            
            REAL mins = -(REAL)nl/2;
            
            for (REAL x = mins; x < -mins; x+=1.)
                for (REAL y = mins; y < -mins; y+=1.)
                    for (REAL z = mins; z < -mins; z+=1.) {

                        //                        cout << x << " " << y << " " << z << endl;
                        Vect3 OM = Vect3((0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX));
                        Vect3 PX = Vect3(x, y, z);
                        OctreeCapsule C(PX, 1 * OM, true);
                        C.AssociatedBody = 0;
                        globalOctree->Root->EvalCapsule(C);


                        if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {

                            string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
                            REAL MEM_PERCENT, temp;
                            stringstream psdata;
                            psdata << top_data;
                            psdata >> temp >> MEM_PERCENT;
                            cout << FVMCell::NumCells << " " << Node::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
                            if (MEM_PERCENT > MAXMEM) {
                                cout << setfill('!') << setw(80) << "!" << endl;

                                cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
                                cout << "Memory used: " << MEM_PERCENT << "%" << endl;

                                cout << setfill('!') << setw(80) << "!" << endl;
                                throw OutOfMemory();
                            }
                        }




                    }
        }
        
            
//            
//            while (FVMCell::NumCells < n) {
//                REAL x = rho * (0.5 - REAL(rand()) / RAND_MAX);
//                REAL y = rho * (0.5 - REAL(rand()) / RAND_MAX);
//                REAL z = rho * (0.5 - REAL(rand()) / RAND_MAX);
//                Vect3 OM = Vect3((0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX));
//                Vect3 PX = Vect3(x, y, z);
//                OctreeCapsule C(PX, 10000 * OM, true);
//                C.AssociatedBody = 0;
//                globalOctree->Root->EvalCapsule(C);
//
//
//                if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {
//
//                    string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
//                    REAL MEM_PERCENT, temp;
//                    stringstream psdata;
//                    psdata << top_data;
//                    psdata >> temp >> MEM_PERCENT;
//                    cout << FVMCell::NumCells << " " << Node::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
//                    if (MEM_PERCENT > MAXMEM) {
//                        cout << setfill('!') << setw(80) << "!" << endl;
//
//                        cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
//                        cout << "Memory used: " << MEM_PERCENT << "%" << endl;
//
//                        cout << setfill('!') << setw(80) << "!" << endl;
//                        throw OutOfMemory();
//                    }
//                }
            
        if (sparsity.compare(sparse) == 0) {
            cout << "sparse ";
            int S = 1;
            while (FVMCell::NumCells < n) {
                REAL phi = asin(2 * REAL(rand()) / RAND_MAX - 1);
                REAL theta = 2 * pi * REAL(rand()) / RAND_MAX;
                REAL x = rho * (0.5 - REAL(rand()) / RAND_MAX);
                REAL y = rho * (0.5 - REAL(rand()) / RAND_MAX);
                REAL z = rho * (0.5 - REAL(rand()) / RAND_MAX);
                if (S == 3) {
                    x = rho * (0.5 - round(REAL(rand()) / RAND_MAX));
                    S = 2;
                } else
                    if (S == 2) {
                    y = rho * (0.5 - round(REAL(rand()) / RAND_MAX));
                    S = 3;
                } else
                    if (S == 3) {
                    y = rho * (0.5 - round(REAL(rand()) / RAND_MAX));
                    S = 1;
                }
                Vect3 OM = Vect3((0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX));
                Vect3 PX = Vect3(x, y, z);
                OctreeCapsule C(PX, 10000 * OM, true);
                C.AssociatedBody = 0;
                globalOctree->Root->EvalCapsule(C);
                if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {

                    string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
                    REAL MEM_PERCENT, temp;
                    stringstream psdata;
                    psdata << top_data;
                    psdata >> temp >> MEM_PERCENT;
                    cout << FVMCell::NumCells << " " << FVMCell::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
                    if (MEM_PERCENT > MAXMEM) {
                        cout << setfill('!') << setw(80) << "!" << endl;

                        cout << "Out of memory. Terminating to avoid swapping to disk." << endl;
                        cout << "Memory used: " << MEM_PERCENT << "%" << endl;

                        cout << setfill('!') << setw(80) << "!" << endl;
                        throw OutOfMemory();
                    }
                }
            }
        }
    }

    if (geometry.compare(sphere) == 0) {
        cout << "sphere with a ";
        if (sparsity.compare(dense) == 0) {
            cout << "dense ";
            while (FVMCell::NumCells < n) {
                REAL phi = asin(2 * REAL(rand()) / RAND_MAX - 1);
                REAL theta = 2 * pi * REAL(rand()) / RAND_MAX;
                REAL rho = REAL(rand()) / RAND_MAX;
                rho *= 4095;
                REAL x = rho * cos(phi) * cos(theta);
                REAL y = rho * cos(phi) * sin(theta);
                REAL z = rho * sin(phi);
                Vect3 OM = Vect3((0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX));
                Vect3 PX = Vect3(x, y, z);

                OctreeCapsule C(PX, 10000 * OM, true);
                C.AssociatedBody = 0;
                globalOctree->Root->EvalCapsule(C);
                if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {

                    string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
                    REAL MEM_PERCENT, temp;
                    stringstream psdata;
                    psdata << top_data;
                    psdata >> temp >> MEM_PERCENT;
                    cout << FVMCell::NumCells << " " << FVMCell::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
                    if (MEM_PERCENT > MAXMEM) {
                        cout << setfill('!') << setw(80) << "!" << endl;

                        cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
                        cout << "Memory used: " << MEM_PERCENT << "%" << endl;

                        cout << setfill('!') << setw(80) << "!" << endl;
                        throw OutOfMemory();
                    }
                }
            }
        }
        if (sparsity.compare(sparse) == 0) {
            cout << "sparse ";
            while (FVMCell::NumCells < n) {
                REAL phi = asin(2 * REAL(rand()) / RAND_MAX - 1);
                REAL theta = 2 * pi * REAL(rand()) / RAND_MAX;
                REAL rho = 4095;
                REAL x = rho * cos(phi) * cos(theta);
                REAL y = rho * cos(phi) * sin(theta);
                REAL z = rho * sin(phi);
                Vect3 OM = Vect3((0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX), (0.5 - REAL(rand()) / RAND_MAX));
                Vect3 PX = Vect3(x, y, z);
                OctreeCapsule C(PX, 10000 * OM, true);
                C.AssociatedBody = 0;
                globalOctree->Root->EvalCapsule(C);
                if (!fmod((REAL) FVMCell::NumCells, 50000.0)) {

                    string top_data = "\t\t" + globalGetStdoutFromCommand(globalIO->top_command);
                    REAL MEM_PERCENT, temp;
                    stringstream psdata;
                    psdata << top_data;
                    psdata >> temp >> MEM_PERCENT;
                    cout << FVMCell::NumCells << " " << FVMCell::NumNodes << " mem used: " << MEM_PERCENT << " percent" << endl;
                    if (MEM_PERCENT > MAXMEM) {
                        cout << setfill('!') << setw(80) << "!" << endl;

                        cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
                        cout << "Memory used: " << MEM_PERCENT << "%" << endl;

                        cout << setfill('!') << setw(80) << "!" << endl;
                        throw OutOfMemory();
                    }
                }
            }
        }
    }

    cout << "points arrangement" << endl;




    //    
    //    string line;
    //    ifstream myfile("infile.dat");
    //    cout << "Reading file" << endl;
    //    if (myfile.is_open()) {
    //        while (myfile.good()) {
    //            getline(myfile, line);
    //            
    //            Vect3 X, O;
    //            istringstream strm(line);
    //            strm >> X.x >> X.y >> X.z >> O.x >> O.y >> O.z;
    //            
    //           // cout << line << endl;
    //           // cout << X << " " << O << endl;
    //            Posns.push_back(X);
    //            Omegas.push_back(O);
    //            
    //        }
    //        myfile.close();
    //    }
    //    
    //    cout << "Done." << endl << "Putting Cells In Tree" << endl;

    //    for (int i = 0; i < Posns.size(); ++i) {
    //           
    //    }
	cout << "Done." << endl << "Calculating FMM" << endl;


    REAL t0 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    
    cout << "Resetting" << endl;
    
    
    //  This is a version of the RESET() function, without the neighbour checking
       //  Everything which changes the shape of the tree must be done recursively
    globalOctree->Root->ApplyRecursively(&Node::MarkWithoutLoad, &Node::MarkWithoutLoad, &Node::DoNothing);
    globalOctree->Root->ApplyRecursively(&Node::DoNothing, &Node::CheckLoad, &Node::DoNothing);

     


    globalOctree->AllCells.clear();
    globalOctree->AllBranches.allocate(OCTREE_LEVS);
    globalOctree->BranchCount.assign(OCTREE_LEVS - 1, 0);
    globalOctree->Root->ApplyRecursively(&Branch::BranchCount, &Node::DoNothing, &Node::DoNothing);

    for (int i = 0; i < globalOctree->BranchCount.size(); ++i) {
        if (globalOctree->BranchCount[i] > 0)
            globalOctree->AllBranches[i].assign(globalOctree->BranchCount[i], NULL);
        globalOctree->BranchCount[i] = 0; //  Recycle for use as a pseudo iterator
    }

    globalOctree->CellCount = 0; //  This is used as a pseudo iterator for assigning into AllCells
    globalOctree->AllCells.assign(FVMCell::NumCells, NULL);
    Node::AllNodes.allocate(Node::NumNodes);
    Node::NodeCount = 0;
    Node::UpList.clear();
    Node::DownList.clear();
    globalOctree->Root->ApplyRecursively(&Node::DoNothing, &Node::ReList, &Node::ReList);
    
    
    globalOctree->ResetAllVelsAndFields();
    cout << "done" << endl << "Setting Vels to 0" << endl;
    globalOctree->Root->ApplyRecursively(&Branch::SetVelsZero, &FVMCell::SetVelsZero, &Node::DoNothing);
    cout << "done" << endl << "Getting Vels" << endl;
    globalOctree->GetVels();
//    globalOctree->Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::SetVelsZero, &Node::DoNothing);
//    globalOctree->Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::PassMmnts2Prnt, &Node::DoNothing);
//    globalOctree->Root->ApplyRecursivelyP(&Branch::GetVelField, &Node::DoNothing, &Node::DoNothing);
//    globalOctree->Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::CollapseVField, &Node::DoNothing);
//    //globalOctree->Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::SetVelsEqual, &Node::DoNothing);
    
    REAL t1 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    cout << "Done. Time elapsed: " << t1 - t0 << endl << "Performing Direct Calculation (including velocity gradients)" << endl;
    Posns.clear();
    Omegas.clear();
    n = globalOctree->AllCells.size();
    Posns.allocate(n);
    Omegas.allocate(n);
    Array <Vect3> FMMVels(n);





    for (int i = 0; i < globalOctree->AllCells.size(); ++i) {
        Posns[i] = (globalOctree->AllCells[i]->Position);
        Omegas[i] = (globalOctree->AllCells[i]->Omega);
        FMMVels[i] = (globalOctree->AllCells[i]->Velocity);
    }


    int nIndices = 1000;

    Array <int> Indices;
    REAL Mult = 1.0;
    if (Posns.size() > nIndices) {
        Mult = Posns.size() / (REAL) nIndices;
        Indices = Array <int> (nIndices, 0);
        int count = 0;
        //while (RandIndices.size() < 10000)
        while (count < nIndices) {
            int randint = int (Posns.size()*(REAL(rand()) / RAND_MAX));
            bool isin = false;
            for (int i = 0; i < nIndices; ++i)
                if (Indices[i] == randint)
                    isin = true;

            if (!isin) {
                Indices[count] = randint;
                count++;
            }
        }

    } else {
        Indices = Array <int> (Posns.size(), 0);
        for (int i = 0; i < Posns.size(); ++i)
            Indices[i] = i;
    }


    Array <Vect3> DirectVels(Indices.size());
    
    Array < Array <Vect3> > DirectVelGrads(Indices.size(), Array <Vect3> (3, Vect3(0.0))), NumVelGradients(Indices.size(), Array <Vect3> (3, Vect3(0.0)));;

    REAL t2 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
#pragma omp parallel for
    for (int i = 0; i < Indices.size(); ++i) {
        Vect3 V(0, 0, 0);

        for (int j = 0; j < Posns.size(); ++j) {
            Vect3 D = Posns[Indices[i]] - Posns[j];
            V += UTIL::globalCubicDirectVel(D, Omegas[j]);

            UTIL::globalCubicDirectVelGrads(D, Omegas[j], DirectVelGrads[i]);
        }


        REAL h = 0.001;

        Vect3 VelN(0.), VelS(0.), VelE(0.), VelW(0.), VelT(0.), VelB(0.);

        for (int j = 0; j < Posns.size(); ++j) {
            VelE += UTIL::globalCubicDirectVel(Posns[Indices[i]] + Vect3(h, 0., 0.) - Posns[j], Omegas[j]);
            VelW += UTIL::globalCubicDirectVel(Posns[Indices[i]] - Vect3(h, 0., 0.) - Posns[j], Omegas[j]);
            VelN += UTIL::globalCubicDirectVel(Posns[Indices[i]] + Vect3(0., h, 0.) - Posns[j], Omegas[j]);
            VelS += UTIL::globalCubicDirectVel(Posns[Indices[i]] - Vect3(0., h, 0.) - Posns[j], Omegas[j]);
            VelT += UTIL::globalCubicDirectVel(Posns[Indices[i]] + Vect3(0., 0., h) - Posns[j], Omegas[j]);
            VelB += UTIL::globalCubicDirectVel(Posns[Indices[i]] - Vect3(0., 0., h) - Posns[j], Omegas[j]);
        }
        NumVelGradients[i][0] = (VelE - VelW) / (2. * h);
        NumVelGradients[i][1] = (VelN - VelS) / (2. * h);
        NumVelGradients[i][2] = (VelT - VelB) / (2. * h);
        DirectVels[i] = V;
    }
    REAL t3 = (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    cout << "Done. Time elapsed: " << t3 - t2 << endl << "Calculating error L2 norm";
    REAL l2 = 0, Vmean = 0;


    Array <REAL> AbsErrs(Indices.size()), RelErrs(Indices.size());

    REAL AbsErrsSum2 = 0., RelErrsSum2 = 0., AbsErrsSum = 0., RelErrsSum = 0., AbsErrsMax = 0., RelErrsMax = 0.;

    RunningStat AbsErrsStats, RelErrsStats;

    REAL ErrorPercent = 0.0, MaxErr = 0.0;
    for (int i = 0; i < Indices.size(); ++i) {
        cout << std::scientific << std::setprecision(3) << "---" << endl << "Vel error percent: " << 100 * (DirectVels[i] - FMMVels[Indices[i]]).Mag() / FMMVels[Indices[i]].Mag() << "%. " << DirectVels[i] << " " << FMMVels[Indices[i]] << endl;
        cout << "xGrads: err " << 100 * (DirectVelGrads[i][0] - globalOctree->AllCells[Indices[i]]->VelGrads[0]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[0].Mag() << " " << globalOctree->AllCells[Indices[i]]->VelGrads[0] << " \t" << DirectVelGrads[i][0] << "\t" << NumVelGradients[i][0] << endl;
        cout << "yGrads: err " << 100 * (DirectVelGrads[i][1] - globalOctree->AllCells[Indices[i]]->VelGrads[1]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[1].Mag() << " " << globalOctree->AllCells[Indices[i]]->VelGrads[1] << " \t" << DirectVelGrads[i][1] << "\t" << NumVelGradients[i][1] << endl;
        cout << "zGrads: err " << 100 * (DirectVelGrads[i][2] - globalOctree->AllCells[Indices[i]]->VelGrads[2]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[2].Mag() << " " << globalOctree->AllCells[Indices[i]]->VelGrads[2] << " \t" << DirectVelGrads[i][2] << "\t" << NumVelGradients[i][2] << endl;
        //cout << "Cell " << i << "\tError L2 Norm " << (DirectVels[i] - FMMVels[Indices[i]]).Mag() << " \t " << DirectVels[i] << endl << "\t\t\t\t\t\t " << FMMVels[Indices[i]] << endl;
        l2 += (DirectVels[i] - FMMVels[Indices[i]]).Mag()*(DirectVels[i] - FMMVels[Indices[i]]).Mag();
        Vmean += DirectVels[i].Mag();

        AbsErrs[i] = (DirectVels[i] - FMMVels[Indices[i]]).Mag();
        RelErrs[i] = AbsErrs[i] / FMMVels[Indices[i]].Mag();

        AbsErrsStats.Push(AbsErrs[i]);
        RelErrsStats.Push(RelErrs[i]);

        AbsErrsSum += AbsErrs[i];
        RelErrsSum += RelErrs[i];
        AbsErrsSum2 += AbsErrs[i] * AbsErrs[i];
        RelErrsSum2 += RelErrs[i] * RelErrs[i];

        AbsErrsMax = max(AbsErrsMax, AbsErrs[i]);
        RelErrsMax = max(RelErrsMax, RelErrs[i]);

        ErrorPercent += 100 * (DirectVels[i] - FMMVels[Indices[i]]).Mag() / FMMVels[Indices[i]].Mag();

        MaxErr = max(MaxErr, 100 * (DirectVels[i] - FMMVels[Indices[i]]).Mag() / FMMVels[Indices[i]].Mag());


    }


    REAL AbsErrsStd = sqrt(AbsErrsSum2 / (Posns.size() - 1));
    REAL RelErrsStd = sqrt(RelErrsSum2 / (Posns.size() - 1));
    REAL AbsErrsMean = AbsErrsSum / Posns.size();
    REAL RelErrsMean = RelErrsSum / Posns.size();

    int TreeSize = OCTREE_SIZE;
    int TreeLevs = OCTREE_LEVS;
    Vmean = Vmean / n;
    cout << "Done." << endl << "Mean Error: " << ErrorPercent/Indices.size() << "%" << endl << "Max Error: " << MaxErr << "%" << endl <<  "Standard deviation: " << 100*RelErrsStats.StandardDeviation() << "%" << endl;

    ofstream myfile;
    myfile.open("output.dat", ios::out | ios::app);
    myfile << globalSystem->MaxP << " " << FVMCell::NumCells << " " << FVMCell::NumNodes << " " << globalSystem->NumThreads << " " << (t1 - t0)*1000 << " " << Mult * (t3 - t2)*1000 << " " << TreeLevs << " " << TreeSize << " " << AbsErrsStats.Max() << " " << AbsErrsStats.Mean() << " " << AbsErrsStats.StandardDeviation() << " " << RelErrsStats.Max() << " " << RelErrsStats.Mean() << " " << RelErrsStats.StandardDeviation() << " " << sparsity << " " << geometry << endl;
    //    myfile << n << " " << globalSystem->MaxP << " " << Vmean << " " << sqrt(l2 / Posns.size()) << " " << sqrt(l2 / Posns.size())/Vmean << " " << (t1 - t0)*1000 << " " << (t3 - t2)*1000 << endl;
    myfile.close();

#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000 << " seconds" << endl;
    cout << "Output written to output.dat:" << endl;
    cout << "p_max FVMCell::NumCells FVMCell::NumNodes nthreads t_fmm (ms)  t_dir (ms) TreeLevs TreeSize AbsErrsMax AbsErrsMean AbsErrsStd RelErrsMax RelErrsMean RelErrsStd sparse/dense sphere/cube;" << endl;
#endif


}

/**************************************************************/
void WeeAmble() {
    /*  This function is to test the PANEL class, in particular 
     *  in creating influence matrices and calculating 
     *  velocities/perturbation potentials. The geomtry is 
     *  from Penfei Lui's PhD thesis, and this serves roughly 
     *  the same purpose as the appendix in which it's found.
     */
    
    
    Array <PANEL> Pans;

    Pans.push_back(PANEL(Vect3(0, 1, 0), Vect3(0, -0.5, 0.866), Vect3(1, 0, 0), Vect3(0, 1, 0)));
    Pans.push_back(PANEL(Vect3(0, -0.5, 0.866), Vect3(0, -0.5, -0.866), Vect3(1, 0, 0), Vect3(0, -0.5, 0.866)));
    Pans.push_back(PANEL(Vect3(0, -0.5, -0.866), Vect3(0, 1, 0), Vect3(1, 0, 0), Vect3(0, -0.5, -0.866)));
    Pans.push_back(PANEL(Vect3(0, 1, 0), Vect3(-1, 0, 0), Vect3(0, -0.5, 0.866), Vect3(0, 1, 0)));
    Pans.push_back(PANEL(Vect3(0, -0.5, 0.866), Vect3(-1, 0, 0), Vect3(0, -0.5, -0.866), Vect3(0, -0.5, 0.866)));
    Pans.push_back(PANEL(Vect3(0, -0.5, -0.866), Vect3(-1, 0, 0), Vect3(0, 1, 0), Vect3(0, -0.5, -0.866)));
    
    
    


    REAL dlta = 1e-6;
    Vect3 dx(dlta, 0, 0), dy(0, dlta, 0), dz(0, 0, dlta);
    
    
    
    Vect3 Target = Vect3(6.,7.,8.);
    Vect3 VDTarget(0.0,0.0,0.0), VSTarget(0.0,0.0,0.0);
    Vect3 VDTargetD(0.0,0.0,0.0), VSTargetD(0.0,0.0,0.0);
    REAL PhiS = 0.0, PhiD = 0.0;
    

    Array <Array <Array <REAL> > > PhiDs, PhiSs;
    PhiDs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));

    Array <Array <Array <REAL> > > PhiDd, PhiSd;
    PhiDd.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSd.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    
    cout << "------" << endl;
    PANEL::FarField = 1e32;
    Array < Array <REAL> > FiS = UTIL::zeros(6, 6), FiD = UTIL::zeros(6, 6), FiSd = UTIL::zeros(6, 6), FiDd = UTIL::zeros(6, 6);
    for (int i = 0; i < 6; ++i) {
        
        for (int j = 0; j < 6; ++j) {
            PhiS = 0.0, PhiD = 0.0;
            Pans[i].GetNormal();
            Pans[j].GetNormal();
            PANEL *src = &Pans[j], *trg = &Pans[i];
            PANEL::SourceDoubletPotential(src, trg->CollocationPoint, PhiD, PhiS, i, j);
            FiS[i][j] = PhiS;
            FiD[i][j] = PhiD;
        }
        Pans[i].GetNormal();
        Pans[i].GetEdgeInfo();
        PANEL *src = &Pans[i];
        src->Gamma = i + 3.14159;
        src->Sigma = i - 3.14159;
        src->Mu = i + 3.14159;
        PANEL::FarField = 1e32;
        for (int I = -1; I < 2; ++I)
            for (int J = -1; J < 2; ++J)
                for (int K = -1; K < 2; ++K)
                    PANEL::SourceDoubletPotential(src, Target + Vect3(I * dlta, J * dlta, K * dlta), PhiDs[I + 1][J + 1][K + 1], PhiSs[I + 1][J + 1][K + 1], 1, 2);
       
            
        PANEL::FarField = 1;
        for (int I = -1; I < 2; ++I)
            for (int J = -1; J < 2; ++J)
                for (int K = -1; K < 2; ++K)
                    PANEL::SourceDoubletPotential(src, Target + Vect3(I * dlta, J * dlta, K * dlta), PhiDd[I + 1][J + 1][K + 1], PhiSd[I + 1][J + 1][K + 1], 1, 2);


        PANEL::FarField = 1e32;
        src->Mu = src->Gamma = i + 3.14159;
        src->Sigma = 0.0;
        VDTarget = VDTarget + src->VortexPanelVelocity(Target);
        src->Mu = src->Gamma = 0;
        src->Sigma = i - 3.14159;
        VSTarget = VSTarget + src->SourceVel(Target);

        PANEL::FarField = 1;
        src->Mu = src->Gamma = i + 3.14159;
        src->Sigma = 0.0;
        VDTargetD = VDTargetD + src->BodyPanelVelocity(Target);
        src->Mu = src->Gamma = 0;
        src->Sigma = i - 3.14159;
        VSTargetD = VSTargetD + src->BodyPanelVelocity(Target);       

    }
    
    
    cout << "Sigma Panel -------" << endl;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j)
            cout  << " " << FiS[i][j];
        cout << ";" << endl;
    }
    cout << "Sigma by Gaussian Quadrature -------" << endl;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j)
            cout << " " << FiSd[i][j];
        cout << ";" << endl;
    }


    cout << "Mu Panel -------" << endl;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j)
            cout  << " " << FiD[i][j];
        cout << ";" << endl;
    }

    cout << "Mu by Gaussian Quadrature -------" << endl;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j)
            cout  << " " << FiDd[i][j];
        cout << ";" << endl;
    }

    cout << "Vels -------" << endl;


    //        PANEL P(Vect3(-0.91, -1.32, 0.15), Vect3(1.4, -2.1, -0.3), Vect3(1.8, 1.6, 1.25), Vect3(-1.2, 1.8, 0.189));

    


   

    


    Vect3 GraD = Vect3(PhiDs[2][1][1] - PhiDs[0][1][1], PhiDs[1][2][1] - PhiDs[1][0][1], PhiDs[1][1][2] - PhiDs[1][1][0]);
    GraD = GraD / (2. * dlta);
    Vect3 GraS = Vect3(PhiSs[2][1][1] - PhiSs[0][1][1], PhiSs[1][2][1] - PhiSs[1][0][1], PhiSs[1][1][2] - PhiSs[1][1][0]);
    GraS = GraS / (2. * dlta);
    cout << "GradPhiD: " << GraD << "\t <-- From O2 c. diff on phi using phi from SourceDoubletPotential()" << endl;
    cout << "Linear D: " << VDTarget << "\t <-- From VortexPanelVel()" << endl;
    cout << "GradPhiS: " << GraS << "\t <-- From O2 c. diff on phi using phi from SourceDoubletPotential()" << endl;
    cout << "Linear S: " << VSTarget << "\t <-- From SourceVel()" << endl;
    
    
    cout << "--------------" << endl;
//    PANEL P(Vect3(-2.1, -3.2, 4.0), Vect3(2.3, -1.4, 3.0), Vect3(1.5, 2.6, 1.0), Vect3(-1.7, 2.8, 0.0));
//        PANEL P(Vect3(-1.2, -1.3, 0.2), Vect3(1.5, -1, 1), Vect3(1, 1, 0.12), Vect3(-1.23, 0.987, -0.3));
    PANEL P(Vect3(-1.,-1., 0.), Vect3(1., -1., 0.), Vect3(1., 1., 0.), Vect3(-1., 1., 0.));

    P.GetNormal();
    
    

    
    
    PANEL *src = &P;
    P.Sigma = 1.0;
    P.Gamma = 1.0;
    P.Mu = 1.0;
    
    PhiD = 0, PhiS = 0;
    
    Vect3 VT = 0.0;

    REAL PhiDoubletDirect = 0.0;

    PhiDs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));

    PANEL::FarField = 1e32;
    for (int I = -1; I < 2; ++I)
        for (int J = -1; J < 2; ++J)
            for (int K = -1; K < 2; ++K)
                PANEL::SourceDoubletPotential(src, Target + Vect3(I * dlta, J * dlta, K * dlta), PhiDs[I + 1][J + 1][K + 1], PhiSs[I + 1][J + 1][K + 1], 1, 2);
    

    VDTarget = src->VortexPanelVelocity(Target);
    VSTarget = src->SourceVel(Target);
    GraD = Vect3(PhiDs[2][1][1] - PhiDs[0][1][1], PhiDs[1][2][1] - PhiDs[1][0][1], PhiDs[1][1][2] - PhiDs[1][1][0]);
    GraD = GraD / (2. * dlta);
    GraS = Vect3(PhiSs[2][1][1] - PhiSs[0][1][1], PhiSs[1][2][1] - PhiSs[1][0][1], PhiSs[1][1][2] - PhiSs[1][1][0]);
    GraS = GraS / (2. * dlta);
    cout << "GradPhiD: " << GraD << "\t <-- From O2 c. diff on phi using phi from SourceDoubletPotential()" << endl;
    
    cout << "GradPhiS: " << GraS << "\t <-- From O2 c. diff on phi using phi from SourceDoubletPotential()" << endl;
    
    

    Array <Array <Array <REAL> > > PhiDD, PhiSD;
    PhiDD.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSD.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));

    Vect3 VelD(0.0,0.0,0.0), VelS(0.0,0.0,0.0);
    REAL PhiSourceDirect = 0.0;
    int npts = 500;
    {
        REAL Us = 0;
        REAL Vs = 0;
        REAL Ws = 0;
        
        REAL Ud = 0;
        REAL Vd = 0;
        REAL Wd = 0;
        
        
        int n = npts;
        Array < Array < Vect3 > > CP;
        Array < Array < Vect3 > > N;
        Array < Array < REAL > > A;
        src->DivPanel(n, CP, N, A);


        REAL PhiSource = 0.0;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {

                Vect3 DX = (CP[i][j] - Target);
                REAL DXMag = DX.Mag();
                REAL DXMag2 = DXMag * DXMag;
                REAL DXMag3 = DXMag2 * DXMag;
             

                Vect3 MU = P.Mu * A[i][j] * N[i][j];


                REAL DXdotMU = DX.Dot(MU);

                REAL phi_d = -DXdotMU / (two_pi * DXMag3);

                VelD += phi_d * (MU / DXdotMU - 3.0 * DX / DXMag2);
                VelS += A[i][j] * DX / (two_pi * DXMag3);   
                 
                PhiSourceDirect += -A[i][j] / (two_pi * DXMag);
                PhiDoubletDirect += -P.Mu * A[i][j] * N[i][j].Dot(DX) / (two_pi * DXMag3);
                     
                for (int I = -1; I < 2; ++I)
                    for (int J = -1; J < 2; ++J)
                        for (int K = -1; K < 2; ++K) {
                            DX = (CP[i][j] - Target + Vect3(I * dlta, J * dlta, K * dlta));
                            DXMag = DX.Mag();
                            DXdotMU = DX.Dot(MU);
                            DXMag3 = DXMag * DXMag * DXMag;
                            PhiSD[I + 1][J + 1][K + 1] += - A[i][j] / (two_pi * DXMag);
                            PhiDD[I + 1][J + 1][K + 1] += - DXdotMU / (two_pi * DXMag3);
                        }
            }


        GraS = Vect3(PhiSD[2][1][1] - PhiSD[0][1][1], PhiSD[1][2][1] - PhiSD[1][0][1], PhiSD[1][1][2] - PhiSD[1][1][0]);
        GraS = GraS / (2. * dlta);
        cout << "GradPhiS: " << GraS << "\t <-- From O2 c. diff on phi using phi from " << n << " subpanels" << endl;
        cout << "Direct S: " << VelS << "\t <-- From analytical gradiants of phis(x,y,z),  using " << n << " subpanels" << endl;
        cout << "Linear S: " << VSTarget << "\t <-- From SourceVel()" << endl;
        cout << "Direct D: " << VelD << "\t <-- From analytical gradiants of phid(x,y,z),  using " << n << " subpanels" << endl;
        GraD = Vect3(PhiDD[2][1][1] - PhiDD[0][1][1], PhiDD[1][2][1] - PhiDD[1][0][1], PhiDD[1][1][2] - PhiDD[1][1][0]);
        GraD= GraD / (2. * dlta);
        cout << "GradPhiD: " << GraD << "\t <-- From O2 c. diff on phi using phi from " << n << " subpanels" << endl;
        cout << "Linear D: " << VDTarget << "\t <-- From VortexPanelVel()" << endl;

    }
    
   
    PANEL::MaxTheta = 0.01;
    PANEL::NumPans = 1;
    PANEL::MaxRecurse = 2;
    Array <PANEL> daPans;
    PanelRecursiveDivide(P, daPans);

    PhiDs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    VDTarget = VSTarget = Vect3(0.0, 0.0, 0.0);
    PANEL::FarField = 1e32;


    for (int i = 0; i < daPans.size(); ++i) {

        PANEL *src = &daPans[i];
        src->GetNormal();
        src->GetEdgeInfo();
        src->Gamma = src->Mu = src->Sigma = 1.0;

        VDTarget += src->VortexPanelVelocity(Target);
        VSTarget += src->SourceVel(Target);

        for (int I = -1; I < 2; ++I)
            for (int J = -1; J < 2; ++J)
                for (int K = -1; K < 2; ++K)
                    PANEL::SourceDoubletPotential(src, Target + Vect3(I * dlta, J * dlta, K * dlta), PhiDs[I + 1][J + 1][K + 1], PhiSs[I + 1][J + 1][K + 1], 1, 2);

    }


    GraS = Vect3(PhiSs[2][1][1] - PhiSs[0][1][1], PhiSs[1][2][1] - PhiSs[1][0][1], PhiSs[1][1][2] - PhiSs[1][1][0]);
    GraS = GraS / (2. * dlta);
    //cout << "GradPhiS: " << GraS << "\t <-- From O2 c. diff on phi using phi from " << PANEL::NumPans << " recursive subpanels, thetamax = " << PANEL::MaxTheta << endl;
    cout << "Linear D:  " << VDTarget << "\t <-- From VortexPanelVel() on " << PANEL::NumPans << " recursive subpanels, thetamax = " << PANEL::MaxTheta << endl;
    cout << "Linear S:  " << VSTarget << "\t <-- From SourceVel() on " << PANEL::NumPans << " recursive subpanels, thetamax = " << PANEL::MaxTheta << endl;
    GraD = Vect3(PhiDs[2][1][1] - PhiDs[0][1][1], PhiDs[1][2][1] - PhiDs[1][0][1], PhiDs[1][1][2] - PhiDs[1][1][0]);
    GraD = GraD / (2. * dlta);
    //cout << "GradPhiD: " << GraD << "\t <-- From O2 c. diff on phi using phi from " << PANEL::NumPans << " recursive subpanels, thetamax = " << PANEL::MaxTheta << endl;

    PhiDs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    PhiSs.assign(3, Array < Array <REAL> > (3, Array < REAL > (3, 0.0)));
    for (int I = -1; I < 2; ++I)
        for (int J = -1; J < 2; ++J)
            for (int K = -1; K < 2; ++K)
                PhiDs[I + 1][J + 1][K + 1] = src->GetTriTesselatedDoubletPhi(Target + Vect3(I * dlta, J * dlta, K * dlta));

    GraD = Vect3(PhiDs[2][1][1] - PhiDs[0][1][1], PhiDs[1][2][1] - PhiDs[1][0][1], PhiDs[1][1][2] - PhiDs[1][1][0]);
    GraD = GraD / (2. * dlta);
    cout << "GradPhiD: " << GraD << "\t <-- From O2 c. diff on phi using phi from 2 triangular subpanels" << endl;
    cout << "DelPhi D  " << src->DoubletPanelVelocity(Target)   << "\t <-- From analyitcal gradients of quad hyperboloidal doublet panel" << endl; 
    cout << "DelPhi S  " << src->SourcePanelVelocity(Target)   << "\t <-- From analyitcal gradients of quad hyperboloidal source panel" << endl;
    cout << "DelPhiTriS" << src->GetTriTesselatedSourceVel(Target)   << "\t <-- From analyitcal gradients of quad hyperboloidal source panel calculated using tris" << endl;
    REAL PhisSDP = 0.0, TEMP;
    PANEL::SourceDoubletPotential(src, Target, TEMP, PhisSDP, 1, 2);
    cout << setprecision(16) << two_pi*PhiDoubletDirect << " PhiDdirect(x,y,z),  using " << npts << "x" << npts << " points on panel surface" << endl;
    cout << setprecision(16) << two_pi*src->GetTriTesselatedDoubletPhi(Target) <<  " GetTriTesselatedDoubletPhi(x,y,z),  using the Willis method" << endl;
    cout << setprecision(16) << two_pi*src->CurvedDoubletPhi(Target) << " CurvedDoubletPhi(x,y,z),  using the Wang method" << endl;
    cout << setprecision(16) << two_pi*src->HyperboloidDoubletPhi(Target) << " HyperboloidDoubletPhi(x,y,z),  using the hyperboloidal panel (e.g. Vaz) method" << endl;

    cout << setprecision(16) << two_pi*PhiSourceDirect << " PhiSdirect(x,y,z),  using " << npts << "x" << npts << " points on panel surface" << endl;
    cout << setprecision(16) << two_pi*src->CurvedSourcePhi(Target) << " CurvedSourcePhi(x,y,z),  using the Wang method" << endl;
    //cout << PhisSDP << " SourceDoubletPotential(x,y,z),  using the flat panel (e.g. Katz & Plotkin) method" << endl;
    cout << setprecision(16) << two_pi*src->HyperboloidSourcePhi(Target) << " HyperboloidSourcePhi(x,y,z),  using the hyperboloidal panel (e.g. Vaz) method" << endl; 
    
    cout << "Hyper PhiD error  " << abs(src->GetTriTesselatedDoubletPhi(Target) - src->HyperboloidDoubletPhi(Target)) << endl;
    cout << "PhiD error " << abs(src->GetTriTesselatedDoubletPhi(Target) - src->CurvedDoubletPhi(Target)) << endl;
    cout << "PhiS error " << abs(PhiSourceDirect - src->CurvedSourcePhi(Target)) << endl;
    
    return;
    
    PANEL::SourceDoubletPotential(&P, Target, PhiD, PhiS, 1, 2);

    cout << PhiD << " " << PhiS << endl;
    PhiD = 0, PhiS = 0;

    cout << PhiD << " " << PhiS << endl;


    //  Test multipole for panels


    srand((unsigned) time(NULL));


    P.Sigma = P.Gamma = P.Mu = 1.0;

    int n = 1000000;
    Array <Vect3> Targets(n);
    for (int i = 0; i < n; ++i) {
        REAL x = 100.0 * (0.5 - REAL(rand()) / REAL(RAND_MAX));
        REAL y = 100.0 * (0.5 - REAL(rand()) / REAL(RAND_MAX));
        REAL z = 100.0 * (0.5 - REAL(rand()) / REAL(RAND_MAX));
        Targets[i] = Vect3(x, y, z);
    }

    REAL tmp;




    // Test to see if it is quicker to calculate velocities or potentials at what would be cell centroids

    // First calc velocities with no farfield
    PANEL::FarField = 1e32;
    unsigned long int t1 = ticks();
    P.Sigma = 1.0;
    P.Mu = 1.0;


    Array <Vect3> VelsFull(n);
    Array <REAL> PhiDFull(n), PhiSFull(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
        VelsFull[i] = P.BodyPanelVelocity(Targets[i]);
    }

    //  Now calc potential
    unsigned long int t2 = ticks();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
        PhiDFull[i] = 0.0;
        PhiSFull[i] = 0.0;
        PANEL::SourceDoubletPotential(&P, Targets[i], PhiDFull[i], PhiSFull[i], 1, 2);
    }

    unsigned long int t3 = ticks();

    PANEL::FarField = 5.0;
    Array <Vect3> VelsFF(n);
    Array <REAL> PhiDFF(n), PhiSFF(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
        VelsFF[i] = P.BodyPanelVelocity(Targets[i]);
    }

    //  Now calc potential
    unsigned long int t4 = ticks();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
        PhiDFF[i] = 0.0;
        PhiSFF[i] = 0.0;
        PANEL::SourceDoubletPotential(&P, Targets[i], PhiDFF[i], PhiSFF[i], 1, 2);
    }
    unsigned long int t5 = ticks();



    REAL EPhiDMag = 0.0, EPhiSMag = 0.0, EVelMag = 0.0;


    for (int i = 0; i < n; ++i) {
        REAL PhiDerr = abs(PhiDFull[i] - PhiDFF[i]);
        REAL r = (Targets[i] - P.Centroid).Mag() / P.MaxDiagonal;
        if (PhiDerr > EPhiDMag) {
            cout << r << " PhiD " << PhiDerr << " " << PhiDFull[i] << " " << PhiDFF[i] << endl;
            EPhiDMag = PhiDerr;
        }

    }

    for (int i = 0; i < n; ++i) {
        REAL PhiSerr = abs(PhiSFull[i] - PhiSFF[i]);
        REAL r = (Targets[i] - P.Centroid).Mag() / P.MaxDiagonal;


        if (PhiSerr > EPhiSMag) {
            cout << r << " PhiS " << PhiSerr << " " << PhiSFull[i] << " " << PhiSFF[i] << endl;
            EPhiSMag = PhiSerr;
        }


    }

    for (int i = 0; i < n; ++i) {
        REAL VErr = (VelsFull[i] - VelsFF[i]).Mag();
        REAL r = (Targets[i] - P.Centroid).Mag() / P.MaxDiagonal;


        if (VErr > EVelMag) {
            cout << r << " Vels " << VErr << " " << VelsFull[i] << " " << VelsFF[i] << endl;
            EVelMag = VErr;
        }
    }

    cout << "t for full velocity calc: " << t2 - t1 << "\t ms" << endl;
    cout << "t for fast velocity calc: " << t4 - t3 << "\t ms" << endl;
    cout << "t for full potental calc: " << t3 - t2 << "\t ms" << endl;
    cout << "t for fast potental calc: " << t5 - t4 << "\t ms" << endl;
    //        
    //    for (int i = 0; i < 1000; ++i) {
    //
    //        
    //        REAL x = 30.0 *(0.5 - REAL(rand()) / REAL(RAND_MAX));
    //        REAL y = 30.0 *(0.5 - REAL(rand()) / REAL(RAND_MAX));
    //        REAL z = 30.0 *(0.5 - REAL(rand()) / REAL(RAND_MAX));
    //        
    //        Target = Vect3(x,y,z);
    //        Vect3 R = Target - P.CollocationPoint;
    //
    //        Vect3 VS = P.SourceVel(Target);
    //        Vect3 VD = P.BodyPanelVelocity(Target);
    //
    //        REAL mult = P.Area / (2.0 * pi * R.Mag() * R.Mag() * R.Mag());
    //        Vect3 VSp = -mult*R;
    //
    //
    //        
    //        
    //      
    //
    //
    //
    //        REAL dlta = 1e-6;
    //        Vect3 dx(dlta, 0, 0), dy(0, dlta, 0), dz(0, 0, dlta);
    //
    //        REAL PhiDE = 0, PhiDW = 0, PhiDEd = 0, PhiDWd = 0;
    //        REAL PhiSE = 0, PhiSW = 0, PhiSEd = 0, PhiSWd = 0;
    //
    //        P.SubPan(100, Target + dx, 1.0, 1.0, PhiDE, PhiSE, VT);
    //        P.SubPan(100, Target - dx, 1.0, 1.0, PhiDW, PhiSW, VT);
    //        PANEL::SourceDoubletPotential(src, Target + dx, PhiDEd, PhiSEd, 1, 2);
    //        PANEL::SourceDoubletPotential(src, Target - dx, PhiDWd, PhiSWd, 1, 2);
    //
    //        REAL PhiDN = 0, PhiDS = 0, PhiDNd = 0, PhiDSd = 0;
    //        REAL PhiSN = 0, PhiSS = 0, PhiSNd = 0, PhiSSd = 0;
    //
    //        P.SubPan(100, Target + dy, 1.0, 1.0, PhiDN, PhiSN, VT);
    //        P.SubPan(100, Target - dy, 1.0, 1.0, PhiDS, PhiSS, VT);
    //        PANEL::SourceDoubletPotential(src, Target + dy, PhiDNd, PhiSNd, 1, 2);
    //        PANEL::SourceDoubletPotential(src, Target - dy, PhiDSd, PhiSSd, 1, 2);
    //
    //
    //        REAL PhiDT = 0, PhiDB = 0, PhiDTd = 0, PhiDBd = 0;
    //        REAL PhiST = 0, PhiSB = 0, PhiSTd = 0, PhiSBd = 0;
    //
    //        P.SubPan(100, Target + dz, 1.0, 1.0, PhiDT, PhiST, VT);
    //        P.SubPan(100, Target - dz, 1.0, 1.0, PhiDB, PhiSB, VT);
    //        PANEL::SourceDoubletPotential(src, Target + dz, PhiDTd, PhiSTd, 1, 2);
    //        PANEL::SourceDoubletPotential(src, Target - dz, PhiDBd, PhiSBd, 1, 2);
    //
    //
    //        Vect3 GraD100 = Vect3((PhiDE - PhiDW), (PhiDN - PhiDS), (PhiDT - PhiDB));
    //        GraD100 = GraD100 / (2 * dlta);
    //        Vect3 GraD = Vect3((PhiDEd - PhiDWd), (PhiDNd - PhiDSd), (PhiDTd - PhiDBd));
    //        GraD = GraD / (2 * dlta);
    //        Vect3 GraS = Vect3((PhiSEd - PhiSWd), (PhiSNd - PhiSSd), (PhiSTd - PhiSBd));
    //        GraS = GraS / (2 * dlta);
    //        REAL Rmag = R.Mag();
    //        REAL R5 = Rmag*Rmag*Rmag*Rmag*Rmag;
    //        REAL denom = R5;
    //
    //        REAL U = R.x*R.z/denom;
    //        REAL V = R.y*R.z/denom;
    //        REAL W = (R.x*R.x + R.y*R.y - 2*R.z*R.z)/denom;
    //        Vect3 VDp(3 * U * P.Area / two_pi, 3 * V * P.Area / two_pi, -W * P.Area / two_pi);
    //        P.Gamma = 1.0;
    //        cout << "LinD: \t" << P.VortexPanelVelocity(Target) << "\t <-- VortexPanelVelocity()" << endl;
    //        cout << " VDp: \t" << VDp << "\t <-- From PointDbl()" << endl;
    //        P.Sigma = 0.0;
    //        P.Mu = 1.0;
    //        cout << " VDp: \t" << P.BodyPanelVelocity(Target) << "\t <-- From BodyPanelVelocity()" << endl;
    //
    //        cout << "GraD: \t" << GraD << "\t <-- From c.diff on phi using phi from SourceDoubletPotential" << endl;
    //        P.Sigma = 1.0;
    //        P.Mu = 0.0;
    //
    //        cout << "LinO: \t" << P.SourceVel(Target) << "\t <-- From SourceVel()" << endl;
    //        cout << " VSp: \t" << P.BodyPanelVelocity(Target) << "\t <-- From BodyPanelVelocity()" << endl;
    //        cout << " VSp: \t" << VSp << "\t <-- From PointSrc()" << endl;
    //        cout << "GraS: \t" << GraS << "\t <-- From c.diff on phi using phi from SourceDoubletPotential" << endl;
    //        cout << R.Mag()/P.MaxDiagonal << endl;
    //    }


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
/**************************************************************/
REAL PanelMaxTheta(PANEL &P)
{
    Vect3 N1 = ((P.C2-P.C1).Cross(P.C1-P.C4)).Normalise();
    Vect3 N2 = ((P.C3-P.C2).Cross(P.C2-P.C1)).Normalise();
    Vect3 N3 = ((P.C4-P.C3).Cross(P.C3-P.C2)).Normalise();
    Vect3 N4 = ((P.C1-P.C4).Cross(P.C4-P.C3)).Normalise();
    
    
    REAL MaxTheta = max(acosd(N1.Dot(N2)),acosd(N1.Dot(N3)));
    MaxTheta = max(MaxTheta,acosd(N1.Dot(N4)));
    MaxTheta = max(MaxTheta,acosd(N2.Dot(N3)));
    MaxTheta = max(MaxTheta,acosd(N2.Dot(N4)));
    MaxTheta = max(MaxTheta,acosd(N3.Dot(N4)));
    
    return MaxTheta;
    
}
/**************************************************************/
void PanelTriangleDivide(PANEL &P, Array <PANEL> &Output)
{
    Array <PANEL>  tmp(4);
    Vect3 C1 = P.C1, C2 = P.C2, C3 = P.C3, C4 = P.C4;
    Vect3 CP = 0.25*(C1 + C2 + C3 + C4);
   
    //  Subdivide panel
    

    PANEL::NumPans += 4;
    tmp[0] = PANEL(P.C1, P.C2, CP, P.C1);
    tmp[0].GetNormal();
    
    tmp[1] = PANEL(P.C2, P.C3, CP, P.C2);
    tmp[1].GetNormal();
    
    tmp[2] = PANEL(P.C3, P.C4, CP, P.C3);
    tmp[2].GetNormal();
    
    tmp[3] = PANEL(P.C4, P.C1, CP, P.C4);
    tmp[3].GetNormal();
    
    Output.push_back(tmp[0]);
    Output.push_back(tmp[1]);
    Output.push_back(tmp[2]);
    Output.push_back(tmp[3]);

}
/**************************************************************/
void PanelRecursiveDivide(PANEL &P, Array <PANEL> &Output)
{
    
    PANEL::RecurseLev += 1;
    Array <PANEL>  tmp(4);
    //  Get curvature of original panel
    
    Vect3 C1 = P.C1, C2 = P.C2, C3 = P.C3, C4 = P.C4;
    Vect3 CP = 0.25*(C1 + C2 + C3 + C4);
   
    //  Subdivide panel
    
    Vect3 C12 = 0.5*(P.C1 + P.C2);
    Vect3 C23 = 0.5*(P.C2 + P.C3);
    Vect3 C34 = 0.5*(P.C3 + P.C4);
    Vect3 C41 = 0.5*(P.C4 + P.C1);
    PANEL::NumPans += 3;

    tmp[0] = PANEL(P.C1, C12, CP, C41);
    tmp[0].GetNormal();
    REAL th1 = PanelMaxTheta(tmp[0]);
    
    tmp[1] = PANEL(C12, P.C2, C23, CP);
    tmp[1].GetNormal();
    REAL th2 = PanelMaxTheta(tmp[1]);
    
    tmp[2] = PANEL(CP, C23, P.C3, C34);
    tmp[2].GetNormal();
    REAL th3 = PanelMaxTheta(tmp[2]);
    
    tmp[3] = PANEL(C41, CP, C34, P.C4);
    tmp[3].GetNormal();
    REAL th4 = PanelMaxTheta(tmp[3]);
    
    
    


    //cout << th1 << " " << th2 << " " << th3 << " " << th4 << " " << PANEL::NumPans << endl;
    if (th1 > PANEL::MaxTheta)
        if (PANEL::RecurseLev > PANEL::MaxRecurse)
            PanelTriangleDivide(tmp[0], Output);
        else
            PanelRecursiveDivide(tmp[0], Output);
    else
        Output.push_back(tmp[0]);

    if (th2 > PANEL::MaxTheta)
        if (PANEL::RecurseLev > PANEL::MaxRecurse)
            PanelTriangleDivide(tmp[1], Output);
        else
            PanelRecursiveDivide(tmp[1], Output);
    else
        Output.push_back(tmp[1]);

    if (th3 > PANEL::MaxTheta)
        if (PANEL::RecurseLev > PANEL::MaxRecurse)
            PanelTriangleDivide(tmp[2], Output);
        else
            PanelRecursiveDivide(tmp[2], Output);
    else
        Output.push_back(tmp[2]);

    if (th4 > PANEL::MaxTheta)
        if (PANEL::RecurseLev > PANEL::MaxRecurse)
            PanelTriangleDivide(tmp[3], Output);
        else
            PanelRecursiveDivide(tmp[3], Output);
    else
        Output.push_back(tmp[3]);
    
    PANEL::RecurseLev -= 1;
    
}