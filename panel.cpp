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


/**************************************************************/
#include "includes.hpp"
#include "types.hpp"
#include "panel.hpp"
/***************************************************************
This code accompanies panel.hpp and is used to set up first order
panels on a boundary surface described by the BODY object
 ***************************************************************/

PANEL::~PANEL() {

    for (int i = 0; i < 4; ++i)
        if (Neighb[i])
            for (int j = 0; j < 4; ++j)
                if (Neighb[i]->Neighb[j] == this) Neighb[i]->Neighb[j] = NULL;

    Neighb = NULL;

    delete CollocationPoint;
}

/***************************************************************/
PANEL::PANEL(POINT *P1, POINT *P2, POINT *P3, POINT *P4) : C1(P1), C2(P2), C3(P3), C4(P4) {
    edgeX1 = NULL;
    edgeX2 = NULL;
    Neighb = NULL;
    NeighbNeighb = NULL;
    OtherBoundarySurface = NULL;
    CollocationPoint = new POINT();
    GetCollocationPoint();
    isBound = false;
    Owner = NULL;
    mu_prev = gamma_prev = 0.0;
    Vfmm = Vect3(0.,0.,0.);
    Theta = 0.;
    Ds = Vn = 0.;
}

/**************************************************************/
PANEL::PANEL(const PANEL &src) {
    ID = src.ID;
    Neighb = src.Neighb;
    NeighbNeighb = src.NeighbNeighb;
    if (WRITE_TO_SCREEN) cout << " WARNING!!!!!!!!!!!!!!!!!!!! Someone is trying to copy a panel!" << endl;
}

/**************************************************************/
void PANEL::DecomposeIntoBlobs() {


}

/**************************************************************/
void PANEL::PrintCollocationPointNormal() {
    if (WRITE_TO_SCREEN) cout.setf(ios::fixed, ios::floatfield);
    CollocationPoint->PrintP();
    if (WRITE_TO_SCREEN) cout << TRANS[2] << ";" << endl;
}

/**************************************************************/
void PANEL::PrintCollocationPoint() {
    if (WRITE_TO_SCREEN) cout.setf(ios::fixed, ios::floatfield);
    CollocationPoint->PrintP();
}

/**************************************************************/
void PANEL::PrintCorners() {
    C1->PrintPos();
    C2->PrintPos();
    C3->PrintPos();
    C4->PrintPos();
}

/**************************************************************/
void PANEL::GetEdgeInfo() {

    Vect3 vCorner_g;

    vCorner_g = C1->vP - CollocationPoint->vP;
    Xcb[0] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C2->vP - CollocationPoint->vP;
    Xcb[1] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C3->vP - CollocationPoint->vP;
    Xcb[2] = VectMultMatrix(TRANS, vCorner_g);
    vCorner_g = C4->vP - CollocationPoint->vP;
    Xcb[3] = VectMultMatrix(TRANS, vCorner_g);

    DX.assign(Xcb[1].x - Xcb[0].x,
            Xcb[2].x - Xcb[1].x,
            Xcb[3].x - Xcb[2].x,
            Xcb[0].x - Xcb[3].x);

    DY.assign(Xcb[1].y - Xcb[0].y,
            Xcb[2].y - Xcb[1].y,
            Xcb[3].y - Xcb[2].y,
            Xcb[0].y - Xcb[3].y);

    M = DY / DX;

    D = sqrt(DX * DX + DY * DY);
}

/**************************************************************/
void PANEL::GetNewGlobalPosition() {
    C1->vP = Owner->CG.vP + VectMultMatrix(Owner->TRANS, C1->vO);
    C2->vP = Owner->CG.vP + VectMultMatrix(Owner->TRANS, C2->vO);
    C3->vP = Owner->CG.vP + VectMultMatrix(Owner->TRANS, C3->vO);
    C4->vP = Owner->CG.vP + VectMultMatrix(Owner->TRANS, C4->vO);
}

void PANEL::GetCollocationPoint() {
    //      Centre point of panel (aka PV of panel centre from global origin)

    Centroid = 0.25 * (C1->vP + C2->vP + C3->vP + C4->vP);

    //   First diagonal
    Vect3 D1 = C2->vP - C4->vP;
    //   Second diagonal
    Vect3 D2 = C3->vP - C1->vP;
    //   TRANS[2]
    Vect3 CR = D1.Cross(D2);
    //   local Z (aka unit normal)
    TRANS[2] = CR / CR.Mag();


    Area = 0.5 * CR.Mag();

    Vect3 BR = .5*(C1->vP + C4->vP) - .5*(C2->vP + C3->vP);
    //  Unit vectors of the same
    TRANS[1] = BR / BR.Mag();
    //    Get local X
    TRANS[0] = TRANS[1].Cross(TRANS[2]);

    MaxDiagonal = max(D1.Mag(), D2.Mag());

    CollocationPoint->vP = Centroid - _EPS2 * CR;


    GetEdgeInfo();
}

/**************************************************************/
void PANEL::SourceDoubletPotential(POINT *IP, REAL Mu, REAL Sigma, REAL Phi[])
{
    //      PV from panel centre to POI in local frame
    Vect3 P = VectMultMatrix(TRANS, IP->vP - Centroid);
    REAL MagP = P.Mag();

//    if (MagP < FarField * MaxDiagonal) {
        Vect3 dX1 = P - Xcb[0], dX2 = P - Xcb[1], dX3 = P - Xcb[2], dX4 = P - Xcb[3];
        REAL Pz2 = P.z * P.z;

        Vect4 e(dX1.x * dX1.x + Pz2, dX2.x * dX2.x + Pz2, dX3.x * dX3.x + Pz2, dX4.x * dX4.x + Pz2);
        Vect4 w(dX1.y * dX1.y, dX2.y * dX2.y, dX3.y * dX3.y, dX4.y * dX4.y);
        Vect4 h(dX1.x * dX1.y, dX2.x * dX2.y, dX3.x * dX3.y, dX4.x * dX4.y);


        Vect4 r = sqrt(e + w);

        REAL aPz = sqrt(Pz2);

        Vect4 alpha, beta, Pzr;

        Pzr = P.z * r;

        alpha = (M * e - h) / Pzr;
        beta = (M * permute_1(e) - permute_1(h)) / permute_1(Pzr);
        Vect4 Num = alpha - beta;
        Vect4 Denom = alpha * beta + 1;

        REAL T = 0;

//        if (MagP == 0)
//            Phi[0] = -two_pi;
//        else
//            Phi[0] = Mu / 2;


//        if (aPz > 0) {
//            //  This is basically all the arctan terms rolled into one atan2 call using double angle formulas
//            //  atan2 does not seem to cascade.....
//            Vect4 Pzr = P.z * r;
//            Vect4 alpha = (M * e - h) / Pzr;
//            Vect4 beta = (M * permute_1(e) - permute_1(h)) / permute_1(Pzr);
//            Vect4 delta = alpha - beta;
//            Vect4 A = (alpha - beta) / (1 + alpha * beta);
//            REAL X = (A.a + A.b) / (1 - A.a * A.b);
//            REAL Y = (A.c + A.d) / (1 - A.c * A.d);
//            T = atan((X + Y)/(1 - X * Y));
//            Phi[0] = Mu * T / (four_pi);
//        }


//        REAL temp = T;



        if (aPz > 0) {
            T = 0;
            if (Num.a != 0.) T += atan2(Num.a, Denom.a);
            if (Num.b != 0.) T += atan2(Num.b, Denom.b);
            if (Num.c != 0.) T += atan2(Num.c, Denom.c);
            if (Num.d != 0.) T += atan2(Num.d, Denom.d);
            Phi[0] = Mu * T / (four_pi);
        }

//        if (fabs(temp - T) > 1e-3)  cout << T << " " << temp << endl;
        REAL B = 0;
        B += (D.a == 0. && (r.a + r.b - D.a) > 0.) ? 0.0 : ((dX1.x * DY.a - dX1.y * DX.a) * log((r.a + r.b + D.a) / (r.a + r.b - D.a))) / D.a;
        B += (D.b == 0. && (r.b + r.c - D.b) > 0.) ? 0.0 : ((dX2.x * DY.b - dX2.y * DX.b) * log((r.b + r.c + D.b) / (r.b + r.c - D.b))) / D.b;
        B += (D.c == 0. && (r.c + r.d - D.c) > 0.) ? 0.0 : ((dX3.x * DY.c - dX3.y * DX.c) * log((r.c + r.d + D.c) / (r.c + r.d - D.c))) / D.c;
        B += (D.d == 0. && (r.d + r.a - D.d) > 0.) ? 0.0 : ((dX4.x * DY.d - dX4.y * DX.d) * log((r.d + r.a + D.d) / (r.d + r.a - D.d))) / D.d;

        Phi[1] = Sigma * (B - P.z * T) / (four_pi);
//    } else {
//        Phi[0] = - Mu * Area * fabs(P.z) / (MagP * MagP * MagP * four_pi);
//        Phi[1] = - Sigma * Area / (MagP * four_pi);
//    }
}

/**************************************************************/
Vect3 PANEL::SourceVel(Vect3 pTarget) {
    //      PV from panel centre to POI in local frame
    Vect3 P = VectMultMatrix(TRANS, pTarget - Centroid);
    REAL MagP = P.Mag(), Mult = *sigma / four_pi;

    if (MagP < FarField * MaxDiagonal) {

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
        return VectMultMatrixTranspose(TRANS, Vect3(sum(B * DY) * Mult, sum(B * DX) * Mult, T * Mult));
    } else
        return VectMultMatrixTranspose(TRANS, Mult * Area * P / (MagP * MagP * MagP));

}
/**************************************************************/
void PANEL::WakeNeighbSet() {
    doSide = true;
    gammaSide = gamma;
    
//        if (Neighb.B) {
//            doSide.B = false;
//            Neighb.B->gammaSide.T = Neighb.B->gamma - gamma;
//        }
//        if (Neighb.R) {
//            doSide.R = false;
//            Neighb.R->gammaSide.L = Neighb.R->gamma - gamma;
//        }
//        if (Neighb.L) gammaSide.L -= Neighb.L->gamma;
//        if (Neighb.T) gammaSide.T -= Neighb.T->gamma;
}

/**************************************************************/
Vect3 PANEL::WakePanelVelocity(Vect3 pTarget) {
    Vect3 V;
    V += LineVelocity(C1->vP, C2->vP, pTarget, gammaSide.L);
    V += LineVelocity(C2->vP, C3->vP, pTarget, gammaSide.T);
    if (doSide.R)
        V += LineVelocity(C3->vP, C4->vP, pTarget, gammaSide.R);
    if (doSide.B)
        V += LineVelocity(C4->vP, C1->vP, pTarget, gammaSide.B);
    return V;
}

/**************************************************************/
Vect3 PANEL::BodyPanelVelocity(Vect3 pTarget, REAL mu) {
    return LineVelocity(C1->vP, C2->vP, pTarget, mu) +
            LineVelocity(C2->vP, C3->vP, pTarget, mu) +
            LineVelocity(C3->vP, C4->vP, pTarget, mu) +
            LineVelocity(C4->vP, C1->vP, pTarget, mu);
}
/**************************************************************/
Vect3 PANEL::BodyPanelVelocity(Vect3 pTarget) {
    return BodyPanelVelocity(pTarget, *(this->mu));
}
/**************************************************************/
Vect3 PANEL::LineVelocity(Vect3 lineStart, Vect3 lineEnd, Vect3 pTarget, const REAL gamma_in) {
    Vect3 R1 = pTarget - lineStart, R2 = pTarget - lineEnd, C = R1.Cross(R2);
    REAL MagC = C.Mag(), MagR1 = R1.Mag(), MagR2 = R2.Mag();
    Vect3 Vout;
    if ((MagR1 > 0.) && (MagR2 > 0.) && (MagC > 0.)) {
        Vect3 R0 = lineEnd - lineStart;
        REAL Mult = gamma_in / (MagC * MagC * four_pi);
        REAL K = Mult * ((R0.Dot(R1) / MagR1) - (R0.Dot(R2) / MagR2));

        Vout -= K*C;
    }
    return Vout;
}
/**************************************************************/
void PANEL::CheckNeighb(PANEL *Face) {
    if (Face != this) {
//		int RECIP[4] = { 2, 3, 0, 1 };
		POINT *LHS[5] = { this->C1, this->C2, this->C3, this->C4, this->C1 };
		POINT *RHS[5] = { Face->C1, Face->C2, Face->C3, Face->C4, Face->C1 };

		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j) {
				if (LHS[i]->vP == RHS[j]->vP)
					for (int k = 0; k < 5; ++k)
						if (LHS[i + 1]->vP == RHS[k]->vP) {
							Neighb[i] = Face;
							//NeighbNeighb[i] = Face->Neighb[RECIP[j]];
							return;
						}
			}
	}
}
/**************************************************************/
void PANEL::GetCp(REAL dt) {
	//	All panels have four neighbours, just some are on different faces
	PanelNeighbSet <REAL> Ds, Dmu;
	PanelNeighbSet <Vect3> Dx;

	Dx.L = CollocationPoint->vP - Neighb.L->CollocationPoint->vP;
	Dx.R = Neighb.R->CollocationPoint->vP - CollocationPoint->vP;
	Dx.T = Neighb.T->CollocationPoint->vP - CollocationPoint->vP;
	Dx.B = CollocationPoint->vP - Neighb.B->CollocationPoint->vP;

	Dmu.L = *(mu) - *(Neighb.L->mu);
	Dmu.R = *(Neighb.R->mu) - *(mu);
	Dmu.T = *(Neighb.T->mu) - *(mu);
	Dmu.B = *(mu) - *(Neighb.B->mu);

	Ds.L = Dx.L.Mag();
	Ds.R = Dx.R.Mag();
	Ds.T = Dx.T.Mag();
	Ds.B = Dx.B.Mag();


	if (Theta.L > pi/3)
		Dx.L = Dmu.L = Ds.L = 0.;

	if (Theta.R > pi/3)
		Dx.R = Dmu.R = Ds.R = 0.;

	if (Theta.T > pi/3)
		Dx.T = Dmu.T = Ds.T = 0.;

	if (Theta.B > pi/3)
		Dx.B = Dmu.B = Ds.B = 0.;

	REAL DmuDeta = (Dmu.L + Dmu.R)/(Ds.L + Ds.R);
	REAL DmuDxi = (Dmu.T + Dmu.B)/(Ds.T + Ds.B);

	//	Now get face vels parallel with zeta and eta
	Vect3 eta = Dx.L + Dx.R;
	eta = eta/eta.Mag();
	Vect3 xi = Dx.T + Dx.B;
	xi = xi/xi.Mag();

	REAL Veta = CollocationPoint->vV.Dot(eta) + DmuDeta;
	REAL Vxi = CollocationPoint->vV.Dot(xi) + DmuDxi;

	REAL DphiDt = (*mu - gamma_prev)/dt;

	REAL Vref2 = ((CollocationPoint->vV).Dot(CollocationPoint->vV));

	Cpress = 1 - (Veta*Veta + Vxi*Vxi - 2*DphiDt)/Vref2;
}
