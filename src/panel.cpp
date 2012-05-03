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


/**************************************************************/
#include "includes.hpp"
#include "types.hpp"
#include "panel.hpp"


#include "panel.hpp"
#include "body.hpp"
#include "pgesv.hpp"

REAL PANEL::FarField = 5.0;
/**************************************************************/
void PANEL::LinearSubPan(PANEL *trg, int n, REAL Mu1, REAL Mu2, REAL &PhiD, Vect3 &V)
{

    V = Vect3(0.0);
    Array <Vect3> C1s = UTIL::globalLinspace(C1,C2,n+1);
    Array <Vect3> C4s = UTIL::globalLinspace(C4,C3,n+1);
    Array <REAL> Mus = UTIL::globalLinspace(Mu1,Mu2,n);
    for (int i = 0; i < n; ++i)
    {
        PANEL T(C1s[i],C1s[i+1],C4s[i+1],C4s[i]);
        T.GetNormal();
        REAL PhiT = 0.0 , tmp;
        PANEL::SourceDoubletPotential(&T, trg->Centroid, PhiT, tmp ,0,1);
        PhiD += (Mus[i]*PhiT);
    }



//    Array < Array < Vect3 > >  CP, N;
//    Array < Array <REAL> > A;
//
//    DivPanel(n, CP, N, A);
//
//    Array <REAL> MuEps1 = UTIL::globalLinspace(MuC1, MuC2, n);
//    Array <REAL> MuEps2 = UTIL::globalLinspace(MuC4, MuC3, n);
//
//    for (int i = 0; i < n; ++i){
//        Array < REAL > MuEta = UTIL::globalLinspace(MuEps1[i],MuEps2[i],n);
//
//        for (int j = 0; j < n; ++j) {
//            Vect3 Mu = MuEta[j] * A[i][j] * N[i][j];
//            Vect3 VT;
//            PANEL::PointDoublet(CP[i][j], P, V, Mu, PhiD);
//        }
//    }
}

/**************************************************************/
void PANEL::DivPanel(int n, Array < Array < Vect3 > > &CP, Array < Array < Vect3 > > &N, Array < Array < REAL > > &A)
{
    Array <Vect3> eps1 =  UTIL::globalLinspace(C1,C2,n+1);
    Array <Vect3> eps2 =  UTIL::globalLinspace(C4,C3,n+1);


    Array < Array <Vect3> > x = UTIL::zerosv(n+1,n+1);
    CP = UTIL::zerosv(n,n), N = UTIL::zerosv(n,n);
    A = UTIL::zeros(n,n);

    for (int i = 0; i < n+1; ++i)
        x[i] = UTIL::globalLinspace(eps1[i],eps2[i],n+1);


    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            CP[i][j] = 0.25*(x[i][j] + x[i+1][j] + x[i+1][j+1] + x[i][j+1]);
            Vect3 dx1 = x[i+1][j+1] - x[i][j];
            Vect3 dx2 = x[i][j+1] - x[i+1][j];

            Vect3 crd = dx1.Cross(dx2);

            A[i][j] = crd.Mag();
            N[i][j] = crd/crd.Mag();

        }

}
/**************************************************************/
REAL PANEL::GetCpD() {
    Vect3 Vinf(globalSystem->unscaledVinf);
    
    Vect3 Vref = Vinf - Vkin; 
    REAL Vref2 = Vref.Dot(Vref); 

    //  Calculate dphi by dt...
    REAL gradT = (Mu - GammaPrev) / (BODY::Time - BODY::TimePrev[0]);
    if (BODY::TimePrev[1] > 0) {
        REAL t0 = BODY::Time;
        REAL t1 = BODY::TimePrev[0];
        REAL t2 = BODY::TimePrev[1];

        Array < Array <REAL> > A = UTIL::zeros(3, 3);

        A[0][0] = 1;
        A[0][1] = t0;
        A[0][2] = t0*t0;
        A[1][0] = 1;
        A[1][1] = t1;
        A[1][2] = t1*t1;
        A[2][0] = 1;
        A[2][1] = t2;
        A[2][2] = t2*t2;

        Array <REAL> RHS(3, 0.0);
        RHS[0] = Mu;
        RHS[1] = PhiPrev[0];
        RHS[2] = PhiPrev[1];

        Array <REAL> x = pgesv(A, RHS);

        gradT = x[1] + 2 * x[2] * t0;

    }
    REAL Vmag = (Vd + VectMultMatrix(TRANS,VCentroid)).Mag();
    REAL Cpress = 1 - (Vmag*Vmag + 0*2 * gradT) / Vref2;
    
    
    return Cpress;
}
/**************************************************************/
REAL PANEL::GetCp() {
    Vect3 Vinf(globalSystem->unscaledVinf);
    //	All panels have four neighbours, just some are on different faces
    PanelNeighbSet <REAL> Ds, DMu;
    PanelNeighbSet <Vect3> Dx;

    if (Neighb[3])
    {    
        Dx[3] = CollocationPoint - Neighb[3]->Centroid;
        DMu[3] = Mu - Neighb[3]->Mu;
    }
    else
    {    
        Dx[3] = Vect3(0.0);
        DMu[3] = 0.0;
    }

    if (Neighb[1])
    {
        Dx[1] = Neighb[1]->CollocationPoint - CollocationPoint;
        DMu[1] = Neighb[1]->Mu - Mu;
    }
    else
    {    
        Dx[1] = Vect3(0.0);
        DMu[1] = 0.0;
    }


    if (Neighb[2])
    {
        Dx[2] = Neighb[2]->CollocationPoint - CollocationPoint;
        DMu[2] = Neighb[2]->Mu - Mu;
    }
    else
    {    
        Dx[2] = Vect3(0.0);
        DMu[2] = 0.0;
    }


    if (Neighb[0])
    {
        Dx[0] = CollocationPoint - Neighb[0]->CollocationPoint;
        DMu[0] = Mu - Neighb[0]->Mu;
    }
    else
    {    
        Dx[0] = Vect3(0.0);
        DMu[0] = 0.0;
    }

        

    Ds[3] = Dx[3].Mag();
    Ds[1] = Dx[1].Mag();
    Ds[2] = Dx[2].Mag();
    Ds[0] = Dx[0].Mag();


    REAL angl = pi / 2;

    if (Theta[3] > angl) {
        Dx[3] = DMu[3] = Ds[3] = 0.;
        DMu[1] = 2 * DMu[1];
    }

    if (Theta[1] > angl) {
        Dx[1] = DMu[1] = Ds[1] = 0.;
        DMu[3] = 2 * DMu[3];
    }

    if (Theta[2] > angl) {
        Dx[2] = DMu[2] = Ds[2] = 0.;
        DMu[0] = 2 * DMu[0];
    }

    if (Theta[0] > angl) {
        Dx[0] = DMu[0] = Ds[0] = 0.;
        DMu[2] = 2 * DMu[2];
    }


    Ds[3] += 1e-16;
    Ds[1] += 1e-16;
    Ds[2] += 1e-16;
    Ds[0] += 1e-16;

    REAL DMuDeta = (DMu[3] / Ds[3] + DMu[1] / Ds[1]) / 2; //(Ds[3] + Ds[1]);
    REAL DMuDxi = (DMu[2] / Ds[2] + DMu[0] / Ds[0]) / 2; //(Ds[2] + Ds[0]);

    //	Now get face vels parallel with zeta and eta
    Vect3 eta = Dx[3] + Dx[1];
    eta = eta / eta.Mag();
    Vect3 xi = Dx[2] + Dx[0];
    xi = xi / xi.Mag();

    Vect3 V;

    //  Iterate over all protowake panels   -   is this necissary or correct here?
//
//        for (int j = 0; j < BODY::AllProtoWakes.size(); ++j){
//            Vect3 U = BODY::AllProtoWakes[j]->WakePanelVelocity(Centroid);
//            V += U;
//        }
//

    REAL Veta = V.Dot(eta) + VCentroid.Dot(eta) + DMuDeta;
    REAL Vxi = V.Dot(xi) + VCentroid.Dot(xi) + DMuDxi;

    
    //Vect3 Vel = Vkin + globalSystem->Vinf;
    Vect3 Vref = Vinf - Vkin; //VCentroid
    REAL Vref2 = Vref.Dot(Vref); //Vkin.Dot(Vkin) + Vinf.Dot(Vinf);

    //  Calculate dphi by dt...
    REAL gradT = (Mu - GammaPrev) / (BODY::Time - BODY::TimePrev[0]);
    if (BODY::TimePrev[1] > 0) {
        REAL t0 = BODY::Time;
        REAL t1 = BODY::TimePrev[0];
        REAL t2 = BODY::TimePrev[1];

        Array < Array <REAL> > A = UTIL::zeros(3, 3);

        A[0][0] = 1;
        A[0][1] = t0;
        A[0][2] = t0*t0;
        A[1][0] = 1;
        A[1][1] = t1;
        A[1][2] = t1*t1;
        A[2][0] = 1;
        A[2][1] = t2;
        A[2][2] = t2*t2;

        Array <REAL> RHS(3, 0.0);
        RHS[0] = Mu;
        RHS[1] = PhiPrev[0];
        RHS[2] = PhiPrev[1];

        Array <REAL> x = pgesv(A, RHS);

        gradT = x[1] + 2 * x[2] * t0;

    }
    REAL Cpress = 1 - (Veta * Veta + Vxi * Vxi + 0*2 * gradT) / Vref2;
    
    //	Check this...
    dF = -Cpress * (.5 * globalSystem->Rho * Vref2) * Area * TRANS[2];
    return Cpress;
}

/**************************************************************/
Vect3 PANEL::VortexPanelVelocity(Vect3 pTarget) {
    Vect3 V;
    V += PANEL::LineVelocity(C1, C2, pTarget, Gamma);
    V += PANEL::LineVelocity(C2, C3, pTarget, Gamma);
    V += PANEL::LineVelocity(C3, C4, pTarget, Gamma);
    V += PANEL::LineVelocity(C4, C1, pTarget, Gamma);
    return V;
}

/**************************************************************/
Vect3 PANEL::BodyPanelVelocity(Vect3 pTarget) {
    Vect3 V;
    //      PV from panel centre to POI in local frame
    Vect3 P = VectMultMatrix(TRANS, pTarget - Centroid);
    
    REAL MagP = P.Mag(), Mult = Sigma / two_pi;

    if (MagP < PANEL::FarField * MaxDiagonal) {

        Vect3 dX1 = P - Xcb[0], dX2 = P - Xcb[1], dX3 = P - Xcb[2], dX4 = P - Xcb[3];
        REAL Pz2 = P.z * P.z;

        Vect4 e(dX1.x * dX1.x + Pz2, dX2.x * dX2.x + Pz2, dX3.x * dX3.x + Pz2, dX4.x * dX4.x + Pz2);
        Vect4 w(dX1.y * dX1.y, dX2.y * dX2.y, dX3.y * dX3.y, dX4.y * dX4.y);

        Vect4 h(dX1.x * dX1.y, dX2.x * dX2.y, dX3.x * dX3.y, dX4.x * dX4.y);

        Vect4 r = sqrt(e + w);

        Vect4 alpha, beta, Pzr;
        REAL aPz = sqrt(Pz2);

        Pzr = P.z * r;
        alpha = (M * e - h) / Pzr;
        beta = (M * permute_1(e) - permute_1(h)) / permute_1(Pzr);
        Vect4 Num = alpha - beta;
        Vect4 Denom = alpha * beta + 1;

        REAL T = 0;

        if (aPz > 0) {
            T = 0;
            if (Num.a != 0.) T += atan2(Num.a, Denom.a);
            if (Num.b != 0.) T += atan2(Num.b, Denom.b);
            if (Num.c != 0.) T += atan2(Num.c, Denom.c);
            if (Num.d != 0.) T += atan2(Num.d, Denom.d);
        }

        Vect4 B((D.a == 0.0) ? 0.0 : log((r.a + r.b - D.a) / (r.a + r.b + D.a)) / D.a,
                (D.b == 0.0) ? 0.0 : log((r.b + r.c - D.b) / (r.b + r.c + D.b)) / D.b,
                (D.c == 0.0) ? 0.0 : log((r.c + r.d - D.c) / (r.c + r.d + D.c)) / D.c,
                (D.d == 0.0) ? 0.0 : log((r.d + r.a - D.d) / (r.d + r.a + D.d)) / D.d);

        //  This return value has been changed......
        Vect3 Vout = VectMultMatrixTranspose(TRANS, Vect3(sum(B * DY) * Mult, sum(B * DX) * Mult, T * Mult));
        V = Vect3(Vout.x, -Vout.y, Vout.z);

        V += PANEL::LineVelocity(C1, C2, pTarget, Mu);
        V += PANEL::LineVelocity(C2, C3, pTarget, Mu);
        V += PANEL::LineVelocity(C3, C4, pTarget, Mu);
        V += PANEL::LineVelocity(C4, C1, pTarget, Mu);
    
    
    } else {
        REAL P3 = (MagP * MagP * MagP);
        Vect3 VS = -Mult * Area * P / P3;


        REAL P5 = MagP * MagP * P3;
        Mult = Area * Mu / (two_pi * P5);

        REAL u = 3 * Mult * P.x * P.z;
        REAL v = 3 * Mult * P.y * P.z;
        REAL w = - Mult * (P.x * P.x + P.y * P.y - 2 * P.z * P.z);
        Vect3 VD(u,v,w);
        V = VectMultMatrixTranspose(TRANS, VS + VD);
    }
        

    return V;
}

/**************************************************************/
Vect3 PANEL::SourceVel(Vect3 pTarget) {
    //      PV from panel centre to POI in local frame
    Vect3 P = VectMultMatrix(TRANS, pTarget - Centroid);
    REAL MagP = P.Mag(), Mult = Sigma / two_pi;

    if (MagP < PANEL::FarField * MaxDiagonal) {

        Vect3 dX1 = P - Xcb[0], dX2 = P - Xcb[1], dX3 = P - Xcb[2], dX4 = P - Xcb[3];
        REAL Pz2 = P.z * P.z;

        Vect4 e(dX1.x * dX1.x + Pz2, dX2.x * dX2.x + Pz2, dX3.x * dX3.x + Pz2, dX4.x * dX4.x + Pz2);
        Vect4 w(dX1.y * dX1.y, dX2.y * dX2.y, dX3.y * dX3.y, dX4.y * dX4.y);

        Vect4 h(dX1.x * dX1.y, dX2.x * dX2.y, dX3.x * dX3.y, dX4.x * dX4.y);

        Vect4 r = sqrt(e + w);

        Vect4 alpha, beta, Pzr;
        REAL aPz = sqrt(Pz2);

        Pzr = P.z * r;
        alpha = (M * e - h) / Pzr;
        beta = (M * permute_1(e) - permute_1(h)) / permute_1(Pzr);
        Vect4 Num = alpha - beta;
        Vect4 Denom = alpha * beta + 1;

        REAL T = 0;

        if (aPz > 0) {
            T = 0;
            if (Num.a != 0.) T += atan2(Num.a, Denom.a);
            if (Num.b != 0.) T += atan2(Num.b, Denom.b);
            if (Num.c != 0.) T += atan2(Num.c, Denom.c);
            if (Num.d != 0.) T += atan2(Num.d, Denom.d);
        }

        Vect4 B((D.a == 0.0) ? 0.0 : log((r.a + r.b - D.a) / (r.a + r.b + D.a)) / D.a,
                (D.b == 0.0) ? 0.0 : log((r.b + r.c - D.b) / (r.b + r.c + D.b)) / D.b,
                (D.c == 0.0) ? 0.0 : log((r.c + r.d - D.c) / (r.c + r.d + D.c)) / D.c,
                (D.d == 0.0) ? 0.0 : log((r.d + r.a - D.d) / (r.d + r.a + D.d)) / D.d);

        //  This return value has been changed......
        Vect3 Vout = VectMultMatrixTranspose(TRANS, Vect3(sum(B * DY) * Mult, sum(B * DX) * Mult, T * Mult));
        return Vect3(Vout.x,-Vout.y,Vout.z);
    } else
        return VectMultMatrixTranspose(TRANS, -Mult * Area * P / (MagP * MagP * MagP));

}

/**************************************************************/
Vect3 PANEL::LineVelocity(Vect3 &lineStart, Vect3 &lineEnd, Vect3 &pTarget, REAL gamma_in) {
    LineVelCnt++;

    Vect3 R1 = pTarget - lineStart;
    Vect3 R2 = pTarget - lineEnd;
    Vect3 C = R1.Cross(R2);
    REAL C2 = C.x * C.x + C.y * C.y + C.z * C.z, MagR1 = R1.Mag(), MagR2 = R2.Mag();
    Vect3 Vout = 0.;



    if ((MagR1 > 0) && (MagR2 > 0) && (C2 > 0)) {
        Vect3 R0 = lineEnd - lineStart;
        REAL Mult = gamma_in / (C2 * two_pi);
        REAL K = Mult * ((R0.Dot(R1) / MagR1) - (R0.Dot(R2) / MagR2));

        Vout += K*C;
    }

//        Vout =  - UTIL::globalDirectVel(pTarget - lineStart, 0.5*gamma_in*(lineEnd - lineStart));
//        Vout =  Vout - UTIL::globalDirectVel(pTarget - lineEnd, 0.5*gamma_in*(lineEnd - lineStart));

    return Vout;
}

/**************************************************************/
void PANEL::GetEdgeInfo() {
    Xcb = Array <Vect3 > (4, Vect3(0.0));
    Vect3 vCorner_g;

    vCorner_g = C1 - Centroid;
    Xcb[0] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C2 - Centroid;
    Xcb[1] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C3 - Centroid;
    Xcb[2] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C4 - Centroid;
    Xcb[3] = VectMultMatrix(TRANS, vCorner_g);

    DX.assign(Xcb[1].x - Xcb[0].x,
            Xcb[2].x - Xcb[1].x,
            Xcb[3].x - Xcb[2].x,
            Xcb[0].x - Xcb[3].x);

    DY.assign(Xcb[1].y - Xcb[0].y,
            Xcb[2].y - Xcb[1].y,
            Xcb[3].y - Xcb[2].y,
            Xcb[0].y - Xcb[3].y);

    M = DY / (DX + _EPS);

    D = sqrt(DX * DX + DY * DY);
}
/**************************************************************/
REAL PANEL::WakePanelPotential(Vect3 target) {

    REAL PhiDoublet = 0.0;
    Vect3 XPg = target - Centroid;
    Vect3 XP = VectMultMatrix(TRANS, XPg);

    if (abs(XP.z) < 1e-12)
        PhiDoublet = 0.5;
    else {
        //  Get source panel LCS

        Vect3 X1 = Xcb[0];
        Vect3 X2 = Xcb[1];
        Vect3 X3 = Xcb[2];
        Vect3 X4 = Xcb[3];

        Vect3 D1 = X2 - X1, D2 = X3 - X2, D3 = X4 - X3, D4 = X1 - X4;
        REAL dx1 = D1.x, dy1 = D1.y,
                dx2 = D2.x, dy2 = D2.y,
                dx3 = D3.x, dy3 = D3.y,
                dx4 = D4.x, dy4 = D4.y;



        
        if (XP.Dot(XP) < (5.0 * MaxDiagonal)*(5.0 * MaxDiagonal)) {

            REAL z = XP.z;

            Vect3 R1 = XP - X1, R2 = XP - X2, R3 = XP - X3, R4 = XP - X4;



            REAL sd1 = D1.Dot(D1), sd2 = D2.Dot(D2), sd3 = D3.Dot(D3), sd4 = D4.Dot(D4);
            REAL r1 = R1.Mag(), r2 = R2.Mag(), r3 = R3.Mag(), r4 = R4.Mag();

            REAL xz1 = (XP.x - X1.x)*(XP.x - X1.x) + XP.z * XP.z;
            REAL xz2 = (XP.x - X2.x)*(XP.x - X2.x) + XP.z * XP.z;
            REAL xz3 = (XP.x - X3.x)*(XP.x - X3.x) + XP.z * XP.z;
            REAL xz4 = (XP.x - X4.x)*(XP.x - X4.x) + XP.z * XP.z;
            REAL xy1 = (XP.x - X1.x)*(XP.y - X1.y);
            REAL xy2 = (XP.x - X2.x)*(XP.y - X2.y);
            REAL xy3 = (XP.x - X3.x)*(XP.y - X3.y);
            REAL xy4 = (XP.x - X4.x)*(XP.y - X4.y);

            REAL Dt1, Dt2, Dt3, Dt4, St1, St2, St3, St4;

            if (sd1 < 1e-12) {
                Dt1 = 0.;
                St1 = 0.;
            } else {
                REAL s11 = dy1 * xz1 - dx1*xy1;
                REAL c11 = r1 * z*dx1;
                REAL s12 = dy1 * xz2 - dx1*xy2;
                REAL c12 = r2 * z*dx1;
                REAL ss1 = s11 * c12 - s12*c11;
                REAL cc1 = c11 * c12 + s11*s12;
                Dt1 = atan2(ss1, cc1);
            }

            if (sd2 < 1e-12) {
                Dt2 = 0.;
                St2 = 0.;
            } else {
                REAL s21 = dy2 * xz2 - dx2*xy2;
                REAL c21 = r2 * z*dx2;
                REAL s22 = dy2 * xz3 - dx2*xy3;
                REAL c22 = r3 * z*dx2;
                REAL ss2 = s21 * c22 - s22*c21;
                REAL cc2 = c21 * c22 + s21*s22;
                Dt2 = atan2(ss2, cc2);
            }

            if (sd3 < 1e-12) {
                Dt3 = 0.;
                St3 = 0.;

            } else {
                REAL s31 = dy3 * xz3 - dx3*xy3;
                REAL c31 = r3 * z*dx3;
                REAL s32 = dy3 * xz4 - dx3*xy4;
                REAL c32 = r4 * z*dx3;
                REAL ss3 = s31 * c32 - s32*c31;
                REAL cc3 = c31 * c32 + s31*s32;
                Dt3 = atan2(ss3, cc3);
            }

            if (sd4 < 1e-12) {
                Dt4 = 0.;
                St4 = 0.;
            } else {
                REAL s41 = dy4 * xz4 - dx4*xy4;
                REAL c41 = r4 * z*dx4;
                REAL s42 = dy4 * xz1 - dx4*xy1;
                REAL c42 = r1 * z*dx4;
                REAL ss4 = s41 * c42 - s42*c41;
                REAL cc4 = c41 * c42 + s41*s42;
                Dt4 = atan2(ss4, cc4);
            }

            PhiDoublet = (Dt1 + Dt2 + Dt3 + Dt4) / two_pi;
        }
        else {
            REAL MagP = XPg.Mag();
            REAL Mult = Area / (4 * two_pi * MagP * MagP * MagP);
            PhiDoublet = XPg.Dot(NC1) + XPg.Dot(NC2) + XPg.Dot(NC3) + XPg.Dot(NC4);
            PhiDoublet *= Mult;
        }
    }
    return PhiDoublet * Gamma;
}
/**************************************************************/
REAL PANEL::BodyPanelPotential(Vect3 target) {
    //  Get source panel LCS
    REAL PhiSource = 0.0, PhiDoublet = 0.0;
    Vect3 XPg = target - Centroid;
    
//    if ((XPg.Dot(XPg) > Mu*ValidRange) && (XPg.Dot(XPg) > Sigma*ValidRange))
//        return 0.0;
    
    if (XPg.Dot(XPg) < ((PANEL::FarField * MaxDiagonal)*(PANEL::FarField * MaxDiagonal)) ){

        Vect3 X1g = C1 - CollocationPoint;
        Vect3 X1 = VectMultMatrix(TRANS, X1g);
        Vect3 X2g = C2 - CollocationPoint;
        Vect3 X2 = VectMultMatrix(TRANS, X2g);
        Vect3 X3g = C3 - CollocationPoint;
        Vect3 X3 = VectMultMatrix(TRANS, X3g);
        Vect3 X4g = C4 - CollocationPoint;
        Vect3 X4 = VectMultMatrix(TRANS, X4g);


        Vect3 D1 = X2 - X1, D2 = X3 - X2, D3 = X4 - X3, D4 = X1 - X4;
        REAL dx1 = D1.x, dy1 = D1.y,
                dx2 = D2.x, dy2 = D2.y,
                dx3 = D3.x, dy3 = D3.y,
                dx4 = D4.x, dy4 = D4.y;



        Vect3 XP = VectMultMatrix(TRANS, XPg);

        REAL x1 = X1.x, x2 = X2.x, x3 = X3.x, x4 = X4.x;
        REAL y1 = X1.y, y2 = X2.y, y3 = X3.y, y4 = X4.y;
        REAL x = XP.x, y = XP.y, z = XP.z;

        Vect3 R1 = XP - X1, R2 = XP - X2, R3 = XP - X3, R4 = XP - X4;



        REAL sd1 = D1.Mag(), sd2 = D2.Mag(), sd3 = D3.Mag(), sd4 = D4.Mag();
        REAL r1 = R1.Mag(), r2 = R2.Mag(), r3 = R3.Mag(), r4 = R4.Mag();

        REAL xz1 = (XP.x - X1.x)*(XP.x - X1.x) + XP.z * XP.z;
        REAL xz2 = (XP.x - X2.x)*(XP.x - X2.x) + XP.z * XP.z;
        REAL xz3 = (XP.x - X3.x)*(XP.x - X3.x) + XP.z * XP.z;
        REAL xz4 = (XP.x - X4.x)*(XP.x - X4.x) + XP.z * XP.z;
        REAL xy1 = (XP.x - X1.x)*(XP.y - X1.y);
        REAL xy2 = (XP.x - X2.x)*(XP.y - X2.y);
        REAL xy3 = (XP.x - X3.x)*(XP.y - X3.y);
        REAL xy4 = (XP.x - X4.x)*(XP.y - X4.y);

        REAL Dt1, Dt2, Dt3, Dt4, St1, St2, St3, St4;

        if (sd1 < 1e-12) {
            Dt1 = 0.;
            St1 = 0.;
        } else {
            REAL s11 = dy1 * xz1 - dx1*xy1;
            REAL c11 = r1 * z*dx1;
            REAL s12 = dy1 * xz2 - dx1*xy2;
            REAL c12 = r2 * z*dx1;
            REAL ss1 = s11 * c12 - s12*c11;
            REAL cc1 = c11 * c12 + s11*s12;
            Dt1 = atan2(ss1, cc1);
            St1 = ((x - x1) * dy1 - (y - y1) * dx1) / sd1 * log((r1 + r2 + sd1) / (r1 + r2 - sd1));
        }

        if (sd2 < 1e-12) {
            Dt2 = 0.;
            St2 = 0.;
        } else {
            REAL s21 = dy2 * xz2 - dx2*xy2;
            REAL c21 = r2 * z*dx2;
            REAL s22 = dy2 * xz3 - dx2*xy3;
            REAL c22 = r3 * z*dx2;
            REAL ss2 = s21 * c22 - s22*c21;
            REAL cc2 = c21 * c22 + s21*s22;
            Dt2 = atan2(ss2, cc2);
            St2 = ((x - x2) * dy2 - (y - y2) * dx2) / sd2 * log((r2 + r3 + sd2) / (r2 + r3 - sd2));
        }

        if (sd3 < 1e-12) {
            Dt3 = 0.;
            St3 = 0.;

        } else {
            REAL s31 = dy3 * xz3 - dx3*xy3;
            REAL c31 = r3 * z*dx3;
            REAL s32 = dy3 * xz4 - dx3*xy4;
            REAL c32 = r4 * z*dx3;
            REAL ss3 = s31 * c32 - s32*c31;
            REAL cc3 = c31 * c32 + s31*s32;
            Dt3 = atan2(ss3, cc3);
            St3 = ((x - x3) * dy3 - (y - y3) * dx3) / sd3 * log((r3 + r4 + sd3) / (r3 + r4 - sd3));
        }

        if (sd4 < 1e-12) {
            Dt4 = 0.;
            St4 = 0.;
        } else {
            REAL s41 = dy4 * xz4 - dx4*xy4;
            REAL c41 = r4 * z*dx4;
            REAL s42 = dy4 * xz1 - dx4*xy1;
            REAL c42 = r1 * z*dx4;
            REAL ss4 = s41 * c42 - s42*c41;
            REAL cc4 = c41 * c42 + s41*s42;
            Dt4 = atan2(ss4, cc4);
            St4 = ((x - x4) * dy4 - (y - y4) * dx4) / sd4 * log((r4 + r1 + sd4) / (r4 + r1 - sd4));
        }
          if (abs(XP.z) < 1e-12)
            PhiDoublet = 0.0;
        else
            PhiDoublet = (Dt1 + Dt2 + Dt3 + Dt4) / four_pi;


        PhiSource = -2 * ((St1 + St2 + St3 + St4) / four_pi - XP.z * PhiDoublet);
        PhiDoublet *= 2;
        
    } else {
        REAL MagP = XPg.Mag();
        REAL Mult = Area / (4 * two_pi * MagP * MagP * MagP);
        PhiDoublet = XPg.Dot(NC1) + XPg.Dot(NC2) + XPg.Dot(NC3) + XPg.Dot(NC4);
        PhiDoublet *= Mult;

        PhiSource = Area / (two_pi * MagP);
    }
    
//    if ((abs(PhiSource) < 1e-6) && (abs(PhiDoublet) < 1e-6))
//        ValidRange = XPg.Mag();
    
    
    return PhiSource*Sigma + PhiDoublet*Mu;
}
/**************************************************************/
void PANEL::SourceDoubletPotential(PANEL *source, Vect3 target, REAL &PhiDoublet, REAL &PhiSource, int i, int j) {
    //  Get source panel LCS
    PhiSource = PhiDoublet = 0.0;
    Vect3 XPg = target - source->Centroid;
    
    if (XPg.Dot(XPg) < ((PANEL::FarField * source->MaxDiagonal)*(PANEL::FarField * source->MaxDiagonal)) ){

        Vect3 X1g = source->C1 - source->CollocationPoint;
        Vect3 X1 = VectMultMatrix(source->TRANS, X1g);
        Vect3 X2g = source->C2 - source->CollocationPoint;
        Vect3 X2 = VectMultMatrix(source->TRANS, X2g);
        Vect3 X3g = source->C3 - source->CollocationPoint;
        Vect3 X3 = VectMultMatrix(source->TRANS, X3g);
        Vect3 X4g = source->C4 - source->CollocationPoint;
        Vect3 X4 = VectMultMatrix(source->TRANS, X4g);


        Vect3 D1 = X2 - X1, D2 = X3 - X2, D3 = X4 - X3, D4 = X1 - X4;
        REAL dx1 = D1.x, dy1 = D1.y,
                dx2 = D2.x, dy2 = D2.y,
                dx3 = D3.x, dy3 = D3.y,
                dx4 = D4.x, dy4 = D4.y;



        Vect3 XP = VectMultMatrix(source->TRANS, XPg);

        REAL x1 = X1.x, x2 = X2.x, x3 = X3.x, x4 = X4.x;
        REAL y1 = X1.y, y2 = X2.y, y3 = X3.y, y4 = X4.y;
        REAL x = XP.x, y = XP.y, z = XP.z;

        Vect3 R1 = XP - X1, R2 = XP - X2, R3 = XP - X3, R4 = XP - X4;



        REAL sd1 = D1.Mag(), sd2 = D2.Mag(), sd3 = D3.Mag(), sd4 = D4.Mag();
        REAL r1 = R1.Mag(), r2 = R2.Mag(), r3 = R3.Mag(), r4 = R4.Mag();

        REAL xz1 = (XP.x - X1.x)*(XP.x - X1.x) + XP.z * XP.z;
        REAL xz2 = (XP.x - X2.x)*(XP.x - X2.x) + XP.z * XP.z;
        REAL xz3 = (XP.x - X3.x)*(XP.x - X3.x) + XP.z * XP.z;
        REAL xz4 = (XP.x - X4.x)*(XP.x - X4.x) + XP.z * XP.z;
        REAL xy1 = (XP.x - X1.x)*(XP.y - X1.y);
        REAL xy2 = (XP.x - X2.x)*(XP.y - X2.y);
        REAL xy3 = (XP.x - X3.x)*(XP.y - X3.y);
        REAL xy4 = (XP.x - X4.x)*(XP.y - X4.y);

        REAL Dt1, Dt2, Dt3, Dt4, St1, St2, St3, St4;

        if (sd1 < 1e-12) {
            Dt1 = 0.;
            St1 = 0.;
        } else {
            REAL s11 = dy1 * xz1 - dx1*xy1;
            REAL c11 = r1 * z*dx1;
            REAL s12 = dy1 * xz2 - dx1*xy2;
            REAL c12 = r2 * z*dx1;
            REAL ss1 = s11 * c12 - s12*c11;
            REAL cc1 = c11 * c12 + s11*s12;
            Dt1 = atan2(ss1, cc1);
            St1 = ((x - x1) * dy1 - (y - y1) * dx1) / sd1 * log((r1 + r2 + sd1) / (r1 + r2 - sd1));
        }

        if (sd2 < 1e-12) {
            Dt2 = 0.;
            St2 = 0.;
        } else {
            REAL s21 = dy2 * xz2 - dx2*xy2;
            REAL c21 = r2 * z*dx2;
            REAL s22 = dy2 * xz3 - dx2*xy3;
            REAL c22 = r3 * z*dx2;
            REAL ss2 = s21 * c22 - s22*c21;
            REAL cc2 = c21 * c22 + s21*s22;
            Dt2 = atan2(ss2, cc2);
            St2 = ((x - x2) * dy2 - (y - y2) * dx2) / sd2 * log((r2 + r3 + sd2) / (r2 + r3 - sd2));
        }

        if (sd3 < 1e-12) {
            Dt3 = 0.;
            St3 = 0.;

        } else {
            REAL s31 = dy3 * xz3 - dx3*xy3;
            REAL c31 = r3 * z*dx3;
            REAL s32 = dy3 * xz4 - dx3*xy4;
            REAL c32 = r4 * z*dx3;
            REAL ss3 = s31 * c32 - s32*c31;
            REAL cc3 = c31 * c32 + s31*s32;
            Dt3 = atan2(ss3, cc3);
            St3 = ((x - x3) * dy3 - (y - y3) * dx3) / sd3 * log((r3 + r4 + sd3) / (r3 + r4 - sd3));
        }

        if (sd4 < 1e-12) {
            Dt4 = 0.;
            St4 = 0.;
        } else {
            REAL s41 = dy4 * xz4 - dx4*xy4;
            REAL c41 = r4 * z*dx4;
            REAL s42 = dy4 * xz1 - dx4*xy1;
            REAL c42 = r1 * z*dx4;
            REAL ss4 = s41 * c42 - s42*c41;
            REAL cc4 = c41 * c42 + s41*s42;
            Dt4 = atan2(ss4, cc4);
            St4 = ((x - x4) * dy4 - (y - y4) * dx4) / sd4 * log((r4 + r1 + sd4) / (r4 + r1 - sd4));
        }



        if (abs(XP.z) < 1e-12)
            PhiDoublet = 0.0;
        else
            PhiDoublet = (Dt1 + Dt2 + Dt3 + Dt4) / four_pi;

        if (i == j)
            PhiDoublet = 0.5;


        PhiSource = -2 * ((St1 + St2 + St3 + St4) / four_pi - XP.z * PhiDoublet);
        PhiDoublet *= 2;
        
    } else {
        REAL MagP = XPg.Mag();
        REAL Mult = source->Area / (4 * two_pi * MagP * MagP * MagP);
        PhiDoublet = XPg.Dot(source->NC1) + XPg.Dot(source->NC2) + XPg.Dot(source->NC3) + XPg.Dot(source->NC4);
        PhiDoublet *= Mult;

        PhiSource = source->Area / (two_pi * MagP);
    }
}
/**************************************************************/

/**************************************************************/
void PANEL::GetNewGlobalPosition() {

    C1 = Owner->CG + VectMultMatrix(Owner->TRANS, C1o - Owner->CGo);
    C2 = Owner->CG + VectMultMatrix(Owner->TRANS, C2o - Owner->CGo);
    C3 = Owner->CG + VectMultMatrix(Owner->TRANS, C3o - Owner->CGo);
    C4 = Owner->CG + VectMultMatrix(Owner->TRANS, C4o - Owner->CGo);




    if (isBound) {
        if (BoundBC == 0) {
            edgeX1 = C1;
            edgeX2 = C2;
        }
        if (BoundBC == 1) {
            edgeX1 = C2;
            edgeX2 = C3;
        }
        if (BoundBC == 2) {
            edgeX1 = C3;
            edgeX2 = C4;
        }
        if (BoundBC == 3) {
            edgeX1 = C4;
            edgeX2 = C1;
        }
    }
    GetNormal();
}

/**************************************************************/
void PANEL::GetNormal() {
    //      Centre point of panel (aka PV of panel centre from global origin)

    Centroid = 0.25 * (C1 + C2 + C3 + C4);

    Vect3 R1 = C2 - C1, R2 = C3 - C2, R3 = C4 - C3, R4 = C1 - C4;
    REAL R1Mag = R1.Dot(R1), R2Mag = R2.Dot(R2), R3Mag = R3.Dot(R3), R4Mag = R4.Dot(R4);


    //   First diagonal
    Vect3 D1 = C3 - C1;
    //   Second diagonal
    Vect3 D2 = C4 - C2;

    Vect3 CR = D1.Cross(D2);

    REAL D1Mag = D1.Mag();
    D1 /= D1Mag;

    REAL D2Mag = D2.Mag();
    D2 /= D2Mag;
    
    
    MaxDiagonal = max(D1Mag, D2Mag);
    REAL MinDiagonal = min(D1Mag, D2Mag);
    //   TRANS[2]

    REAL CRMag = CR.Mag();
    if (MinDiagonal == 0)
        TRANS[2] = Vect3(0.0);
    else
        TRANS[2] = CR/CRMag;


    if (((R1Mag == 0) && (R3Mag == 0)) || ((R2Mag == 0) && (R4Mag == 0))) {
        Centroid = 0.25 * (C1 + C2 + C3 + C4);

        Area = 0.5*(D1.Mag() + D2.Mag());
        
    } else {

            Area = 0.5 * CRMag;

        if (R1Mag < 1e-12)
        Centroid = (C2 + C3 + C4) / 3;

        if (R2Mag < 1e-12)
        Centroid = (C1 + C3 + C4) / 3;

        if (R3Mag < 1e-12)
        Centroid = (C1 + C2 + C4) / 3;

        if (R4Mag < 1e-12)
        Centroid = (C1 + C2 + C3) / 3;
    }   

    R1 /= R1.Mag();
    R2 /= R2.Mag();
    R3 /= R3.Mag();
    R4 /= R4.Mag();
    NC1 = R1.Cross(R4);
    NC2 = R2.Cross(R1);
    NC3 = R3.Cross(R2);
    NC4 = R4.Cross(R3);



    Vect3 VJ = 0.5 * (C3 + C4) - Centroid;


    //  Unit vectors of the same
    
    REAL VJMag = VJ.Mag();
    if (VJMag > 0)
        TRANS[1] = VJ / (VJMag);
    else
        TRANS[1] = Vect3(0.0);
    //    Get local Y
    TRANS[0] = TRANS[1].Cross(TRANS[2]);



    GetEdgeInfo();

    CollocationPoint = Centroid;
    Normal = TRANS[2];

}
/**************************************************************/
void PANEL::CheckNeighb(PANEL *Face) {
    if (Face != this) {
        //  Check C1
        if ((C1 == Face->C1) || (C1 == Face->C2) || (C1 == Face->C3) || (C1 == Face->C4))
            if ((C2 == Face->C1) || (C2 == Face->C2) || (C2 == Face->C3) || (C2 == Face->C4)) {
                Neighb[0] = Face;
                //                cout << "Got Neighb[0] " << endl;
            }

        if ((C2 == Face->C1) || (C2 == Face->C2) || (C2 == Face->C3) || (C2 == Face->C4))
            if ((C3 == Face->C1) || (C3 == Face->C2) || (C3 == Face->C3) || (C3 == Face->C4)) {
                Neighb[1] = Face;
                //                cout << "Got Neighb[1] " << endl;
            }


        if ((C3 == Face->C1) || (C3 == Face->C2) || (C3 == Face->C3) || (C3 == Face->C4))
            if ((C4 == Face->C1) || (C4 == Face->C2) || (C4 == Face->C3) || (C4 == Face->C4)) {
                Neighb[2] = Face;
                //                cout << "Got Neighb[2] " << endl;
            }

        if ((C4 == Face->C1) || (C4 == Face->C2) || (C4 == Face->C3) || (C4 == Face->C4))
            if ((C1 == Face->C1) || (C1 == Face->C2) || (C1 == Face->C3) || (C1 == Face->C4)) {
                Neighb[3] = Face;
                //                cout << "Got Neighb[3] " << endl;
            }

    }


    //        //      int RECIP[4] = { 2, 3, 0, 1 };
    //        Vect3 * LHS[5] = {this->C1, this->C2, this->C3, this->C4, this->C1};
    //        Vect3 * RHS[5] = {Face->C1, Face->C2, Face->C3, Face->C4, Face->C1};
    //
    //        for (int i = 0; i < 4; ++i)
    //            for (int j = 0; j < 4; ++j) {
    //                if (LHS[i] == RHS[j])
    //                    for (int k = 0; k < 5; ++k)
    //                        if (LHS[i + 1] == RHS[k]) {
    //                            cout << "Got neighbour " << i << endl;
    //                            Neighb[i] = Face;
    //                            //NeighbNeighb[i] = Face->Neighb[RECIP[j]];
    //                            return;
    //                        }
    //            }
    //    }


}