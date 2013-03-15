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
#include "unit_tests.hpp"
/**************************************************************/
void TEST::TestBEM() {
    /* This function performs a BEM calculation over single or multiple steps.
     * Use to test influence matrix calculalation, subtime-stepping etc.
     */
    
     SYSTEM System(0);

    //  Some default values
    globalSystem->GambitScale = 1.0;
    globalSystem->MaxP = 5;
    globalSystem->Del2 = 0.001;// * globalSystem->GambitScale*globalSystem->GambitScale;
    globalSystem->DS = .3;
    globalSystem->dtInit = 1e6;
    globalSystem->h = 2;
    globalSystem->unscaledVinf = Vect3(0.0);
    globalSystem->NumSubSteps = 1;


    UTIL::cpu_t = ticks();

    UTIL::PreAmble();

    BODY::BodySubStep(1.0e3, 1);

    UTIL::PostAmble(string("Output.mat"));

#ifndef use_NCURSES
    if (WRITE_TO_SCREEN) cout << "CPU time: " << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000 << " seconds" << endl;
#endif
    
    
}
/**************************************************************/
void TEST::TestPanel() {
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
    PANEL::PanelRecursiveDivide(P, daPans);

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
void TEST::TestBiotSavart()
{
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
//        
//        
//        
//    }
//    
}
void TEST::SolveMatfileVels(string fname, int pmax, REAL del2) {

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
void TEST::TestFMM(int argc, char *argv[]) {

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
        cout << "---" << endl << "Vel error percent: " << 100 * (DirectVels[i] - FMMVels[Indices[i]]).Mag() / FMMVels[Indices[i]].Mag() << "%. " << DirectVels[i] << " " << FMMVels[Indices[i]] << endl;
        cout << "xGrads: err " << 100 * (DirectVelGrads[i][0] - globalOctree->AllCells[Indices[i]]->VelGrads[0]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[0].Mag() << "%. " << globalOctree->AllCells[Indices[i]]->VelGrads[0] << " \t" << DirectVelGrads[i][0] << "\t" << NumVelGradients[i][0] << endl;
        cout << "yGrads: err " << 100 * (DirectVelGrads[i][1] - globalOctree->AllCells[Indices[i]]->VelGrads[1]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[1].Mag() << "%. " << globalOctree->AllCells[Indices[i]]->VelGrads[1] << " \t" << DirectVelGrads[i][1] << "\t" << NumVelGradients[i][1] << endl;
        cout << "zGrads: err " << 100 * (DirectVelGrads[i][2] - globalOctree->AllCells[Indices[i]]->VelGrads[2]).Mag() / globalOctree->AllCells[Indices[i]]->VelGrads[2].Mag() << "%. " << globalOctree->AllCells[Indices[i]]->VelGrads[2] << " \t" << DirectVelGrads[i][2] << "\t" << NumVelGradients[i][2] << endl;
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