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


#include "body.hpp"
/**************************************************************/
BODY::BODY(Vect3 Origin_, Vect3 Attitude, Vect3 Velocity, Vect3 Rates, string name, SYSTEM *sys) {
    globalSystem = sys;
    Name = name;
    CG.vP = Origin_ * globalSystem->GambitScale;
    EulerAngles = Attitude;
    SetEulerTrans();
    BodyRates = Rates;

    CG.vV = Velocity * globalSystem->GambitScale;
    GetEulerRates();

}

/**************************************************************/
BODY::~BODY() {
    for (int i = 0; i < ProtoWake.size(); ++i) delete ProtoWake[i];
    for (int i = 0; i < ProtoWakePoints.size(); ++i) delete ProtoWakePoints[i];
    for (int i = 0; i < WakePoints.size(); ++i)
        for (int j = 0; j < WakePoints[i].size(); ++j) delete WakePoints[i][j];
    for (int i = 0; i < Faces.size(); ++i) delete Faces[i];


}

/**************************************************************/
void BODY::MoveBody(REAL dt) {

    //    First off calculate new Euler angles
    EulerAngles += dt*EulerRates;
    //    Now get new cg position in global reference frame
    CG.vP += CG.vV*dt;
    //    Now set appropriate body rates, and transforms etc.
    SetEulerTrans();
//    GetBodyRates();   //  Body rates are constant in this case....                        //  Check BodyRates
    //      Now update position of body faces

    for (int i = 0; i < (int) ProtoWake.size(); ++i) {
        ProtoWake[i]->GetNewGlobalPosition();
        ProtoWake[i]->GetCollocationPoint();
    }

    for (int i = 0; i < NumFaces; ++i) {
        Faces[i]->GetNewGlobalPosition();
        Faces[i]->GetCollocationPoint();
    }
}

/**************************************************************/
void BODY::SortWake(REAL dt) {
    if (WRITE_TO_SCREEN) cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    PANEL *src, *trg;
    Vect3 Vinf(globalSystem->uinf, globalSystem->vinf, globalSystem->winf);
    Array <PANEL*> ProtoCopy;


    Array <POINT*> new_points, new_points_t;
    for (int i = 0; i < nSpanPanels; ++i) {
        src = ProtoWake[i];

        src->gamma_prev = src->gamma;

        POINT *X1, *X2, *X3, *X4;
        Vect3 P1 = src->C4->vP, P2 = src->C3->vP;
        if (globalSystem->LiftingLineMode) {
            P1 = src->Shedder->C3->vP;
            P2 = src->Shedder->C2->vP;
        }
        //	Copy protowake to a new panel
        bool test1 = false, test2 = false;
        for (int j = 0; j < (int) new_points.size(); ++j) {
            if (new_points[j]->vP == P1) {
                test1 = true;
                X1 = new_points[j];
            }

            if (new_points[j]->vP == P2) {
                test2 = true;
                X2 = new_points[j];
            }
        }

        if (!test1) {
            X1 = new POINT(P1);
            new_points.push_back(X1);
        }

        if (!test2) {
            X2 = new POINT(P2);
            new_points.push_back(X2);
        }

        if (WakeGlobal.size() == 0) {
            //	Allow the trailing points to convect, including OMEGAxR
            // 	Get point kinematic velocity - rotational part first
            Vect3 Vkin3 = CG.vV + BodyRates.Cross(P2 - CG.vP);
            Vect3 Vkin4 = CG.vV + BodyRates.Cross(P1 - CG.vP);

            Vect3 P3 = P2 + dt * (Vinf + Vkin3);
            Vect3 P4 = P1 + dt * (Vinf + Vkin4);

            bool test3 = false, test4 = false;
            for (int j = 0; j < (int) new_points_t.size(); ++j) {
                if (new_points_t[j]->vP == P3) {
                    test3 = true;
                    X3 = new_points_t[j];
                }

                if (new_points_t[j]->vP == P4) {
                    test4 = true;
                    X4 = new_points_t[j];
                }
            }
            if (!test3) {
                X3 = new POINT(P3);
                new_points_t.push_back(X3);
            }


            if (!test4) {
                X4 = new POINT(P4);
                new_points_t.push_back(X4);
            }
        } else {
            X3 = WakeGlobal[0][i]->C2;
            X4 = WakeGlobal[0][i]->C1;
        }


        trg = new PANEL(X1, X2, X3, X4);
        ProtoCopy.push_back(trg);
        if (!globalSystem->LiftingLineMode)
            trg->gamma = src->gamma;
        else
            trg->gamma = src->Shedder->gamma;

    }


    for (int i = 0; i < (int) ProtoCopy.size(); ++i) {
        for (int j = 0; j < (int) ProtoCopy.size(); ++j)
            ProtoCopy[i]->CheckNeighb(ProtoCopy[j]);

        if (WakeGlobal.size() != 0)
            for (int j = 0; j < (int) WakeGlobal[0].size(); ++j) {
                ProtoCopy[i]->CheckNeighb(WakeGlobal[0][j]);
                WakeGlobal[0][j]->CheckNeighb(ProtoCopy[i]);
            }
    }



    if (WakeGlobal.size() == 0) WakePoints.push_back(new_points_t);
    WakePoints.push_back(new_points);


    for (int i = 0; i < ProtoCopy.size(); ++i)
        ProtoCopy[i]->WakeNeighbSet();


    WakeGlobal.push_front(ProtoCopy);


}

/**************************************************************/
void BODY::InitNascentWake(REAL dt) {
    if (globalSystem->LiftingLineMode) dt = 0.0;
    //     int edge;
    for (int i = 0; i < (int) BoundaryFaces.size(); ++i) {
        if (BoundaryFaces[i]->isTop) {
            Vect3 vVinf(globalSystem->uinf, globalSystem->vinf, globalSystem->winf);
            PANEL *src, *trg;
            src = BoundaryFaces[i];
            POINT *X1, *X2, *X3, *X4;
            X1 = new POINT;
            X2 = new POINT;
            X3 = new POINT;
            X4 = new POINT;
			X1->vV = X2->vV = X3->vV = X4->vV = vVinf;
			//	Rotate to previous position
			MoveBody(-dt * DS);
			Vect3 vx1, vx2, vx3, vx4;

			vx4 = src->edgeX1->vP;
			vx3 = src->edgeX2->vP;

			//	Rotate back
			MoveBody(dt * DS);

			vx1 = src->edgeX1->vP;
			vx2 = src->edgeX2->vP;

			//	Convect x1 and x2
			vx3 += (DS * dt) * vVinf;
			vx4 += (DS * dt) * vVinf;
            //	Put positions into POINTs

            X1->vP = vx1;
            X2->vP = vx2;
            X3->vP = vx3;
            X4->vP = vx4;
            X1->vO = X1->vP - CG.vP;
            X2->vO = X2->vP - CG.vP;
            X3->vO = X3->vP - CG.vP;
            X4->vO = X4->vP - CG.vP;

            //	Make Panel
            trg = new PANEL(X1, X2, X3, X4);
            trg->Shedder = src;
            trg->Owner = this;

            if (src->OtherBoundarySurface)
                src->OtherBoundarySurface->Wake = trg;

            src->Wake = trg;
            ProtoWake.push_back(trg);
            ProtoWakePoints.push_back(X1);
            ProtoWakePoints.push_back(X2);
            ProtoWakePoints.push_back(X3);
            ProtoWakePoints.push_back(X4);
        }
    }

    nSpanPanels = (int) ProtoWake.size();

    for (int i = 0; i < nSpanPanels; ++i)
        for (int j = 0; j < nSpanPanels; ++j)
            ProtoWake[i]->CheckNeighb(ProtoWake[j]);
}

/**************************************************************/
void BODY::GetPanels(Array <int> &GROUP, Array <PANEL*> &PANELS) {
    //  For all panels in the group
    for (int i = 0; i < (int) GROUP.size(); ++i) {
        //  get number of the panel
        int id = GROUP[i];
        Faces.push_back(PANELS[id]);
        PANELS[id]->Owner = this;
        PANELS[id]->C1->Owner = this;
        PANELS[id]->C2->Owner = this;
        PANELS[id]->C3->Owner = this;
        PANELS[id]->C4->Owner = this;
        if (PANELS[id]->isBound)
            BoundaryFaces.push_back(PANELS[id]);

    }
    NumFaces = (int) Faces.size();



    for (int i = 0; i < (int) BoundaryFaces.size(); ++i) {
        Vect3 xs1 = BoundaryFaces[i]->edgeX1->vP;
        Vect3 xs2 = BoundaryFaces[i]->edgeX2->vP;
        for (int j = 0; j < (int) BoundaryFaces.size(); ++j) {
            if (i != j) {
                Vect3 xt1 = BoundaryFaces[j]->edgeX1->vP;
                Vect3 xt2 = BoundaryFaces[j]->edgeX2->vP;
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
#ifndef use_NCURSES
        if (WRITE_TO_SCREEN) cout << "┌\t\t\t\t\t\t\t\t┐" << endl;
        if (WRITE_TO_SCREEN) cout << "|\tNumber of wake shedding elements = number of elements\t|" << endl;
        if (WRITE_TO_SCREEN) cout << "|\tUsing Lifting Line Mode\t\t\t\t\t|" << endl;
        if (WRITE_TO_SCREEN) cout << "└\t\t\t\t\t\t\t\t┘" << endl;
#endif
        globalSystem->LiftingLineMode = true;
        for (int i = 0; i < (int) BoundaryFaces.size(); ++i) {
            BoundaryFaces[i]->isTop = true;
            BoundaryFaces[i]->OtherBoundarySurface = NULL;
        }
    }

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
    REAL c2 = cosphi * sinthe * sinpsi-sinphi*cospsi;
    REAL c3 = cosphi*costhe;
    TRANS[0] = Vect3(a1,b1,c1);
    TRANS[1] = Vect3(a2,b2,c2);
    TRANS[2] = Vect3(a3,b3,c3);
}
/**************************************************************/
#ifndef use_NCURSES

void BODY::PrintSurface() {
    if (WRITE_TO_SCREEN) cout.setf(ios::fixed, ios::floatfield);
    if (WRITE_TO_SCREEN) cout << "hold all" << endl;
    for (int j = 0; j < NumFaces; ++j)
        SURF(Faces[j]->C1->vP, Faces[j]->C2->vP, Faces[j]->C3->vP, Faces[j]->C4->vP, *(Faces[j]->mu));
    if (WRITE_TO_SCREEN) cout << "axis equal" << endl;
}
#endif

/**************************************************************/
void BODY::WriteSurface(ostream& out_stream) {
    if (WRITE_TO_FILE) {
        out_stream.setf(ios::fixed, ios::floatfield);
        out_stream << "hold all" << endl;
        for (int j = 0; j < NumFaces; ++j){
            SURFW((1.0)*Faces[j]->C1->vP, (1.0)*Faces[j]->C2->vP, (1.0)*Faces[j]->C3->vP, (1.0)*Faces[j]->C4->vP, Faces[j]->gamma, out_stream);
//            out_stream << "scatter3(" << Faces[j]->C1->vP.x << "," << Faces[j]->C1->vP.y << "," << Faces[j]->C1->vP.z << ",'s');" << endl;
//            out_stream << "scatter3(" << Faces[j]->C2->vP.x << "," << Faces[j]->C2->vP.y << "," << Faces[j]->C2->vP.z << ",'f');" << endl;
//            out_stream << "scatter3(" << Faces[j]->C3->vP.x << "," << Faces[j]->C3->vP.y << "," << Faces[j]->C3->vP.z << ",'.');" << endl;
//            out_stream << "scatter3(" << Faces[j]->C4->vP.x << "," << Faces[j]->C4->vP.y << "," << Faces[j]->C4->vP.z << ");" << endl;
        }
        out_stream << "axis equal" << endl;
    }
}

/**************************************************************/
void BODY::WriteWake(ostream& out_stream) {
    if (WRITE_TO_FILE) {
        out_stream.setf(ios::fixed, ios::floatfield);
        out_stream << "hold all" << endl;
    }
    if (WRITE_TO_FILE)
        for (int j = 0; j < nSpanPanels; ++j)
            SURFW(ProtoWake[j]->C1->vP, ProtoWake[j]->C2->vP, ProtoWake[j]->C3->vP, ProtoWake[j]->C4->vP, 0.0*ProtoWake[j]->gamma, out_stream);


    out_stream << "X0 = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vP.x << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;


    out_stream << "Y0 = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vP.y << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;

    out_stream << "Z0 = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vP.z << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;

    out_stream << "Ox = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vO.x << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;

    out_stream << "Oy = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vO.y << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;

    out_stream << "Oz = [ ";
    for (int i = 0; i < WakePoints.size(); ++i) {
        for (int j = 0; j < WakePoints[i].size(); ++j)
            out_stream << WakePoints[i][j]->vO.z << " ";
        out_stream << ";" << endl;
    }
    out_stream << "];" << endl;

     out_stream << "scatter3(X0(:),Y0(:),Z0(:),50,sqrt(Ox(:).*Ox(:) + Oy(:).*Oy(:) + Oz(:).*Oz(:)),'filled');" << endl;

    out_stream << "view(3); axis tight; axis equal;drawnow;set(gcf,'Renderer','OpenGL');" << endl;
}
/**************************************************************/

void BODY::PrintBoundary() {
}

/**************************************************************/
void BODY::PrintWake() {
    if (WRITE_TO_SCREEN) cout.setf(ios::fixed, ios::floatfield);
    if (WRITE_TO_SCREEN) cout << "hold all" << endl;


    for (int j = 0; j < nSpanPanels; ++j)
        SURF(ProtoWake[j]->C1->vP, ProtoWake[j]->C2->vP, ProtoWake[j]->C3->vP, ProtoWake[j]->C4->vP, ProtoWake[j]->gamma);

    for (int i = 0; i < (int) WakeGlobal.size(); ++i)
        for (int j = 0; j < (int) WakeGlobal[i].size(); ++j)
            SURF(WakeGlobal[i][j]->C1->vP, WakeGlobal[i][j]->C2->vP, ProtoWake[j]->C3->vP, WakeGlobal[i][j]->C4->vP, WakeGlobal[i][j]->gamma);

    if (WRITE_TO_SCREEN) cout << "axis equal" << endl;
}

/**************************************************************/
Vect3 BODY::GetVel(Vect3 Target) {

    // Why can this section not be multi-threaded????
    Vect3 U, V, W;
    if (!globalSystem->LiftingLineMode)
    for (int i = 0; i < (int) ProtoWake.size(); ++i)
        U -= ProtoWake[i]->WakePanelVelocity(Target);

    for (int i = 0; i < WakePoints.size(); ++i)
        for (int j = 0; j < WakePoints[i].size(); ++j)
            globalDirectVel(Target - WakePoints[i][j]->vP, WakePoints[i][j]->vO, V);


    for (int l = 0; l < (int) Faces.size(); ++l) {
        if (!globalSystem->LiftingLineMode) {
            W -= Faces[l]->SourceVel(Target);
        }
        W -= Faces[l]->BodyPanelVelocity(Target);//, globalSystem->Del2);
    }

    return V;//U + V + W;
}

/**************************************************************/
Vect3 BODY::GetWakeVel(Vect3 Target) {
    Vect3 U;
//    for (i = 0; i < (int) WakeGlobal.size(); ++i)
//        for (j = 0; j < (int) WakeGlobal[i].size(); ++j)
//            U += WakeGlobal[i][j]->WakePanelVelocity(Target);

    for (int i = 0; i < WakePoints.size(); ++i)
        for (int j = 0; j < WakePoints[i].size(); ++j)
            globalDirectVel(WakePoints[i][j]->vP - Target, WakePoints[i][j]->vO, U);

    return U;
}

/**************************************************************/
void BODY::DissolveWake(REAL dt) {
    Array <REAL> R, RR, G, GG, GP, GPGP;
    Array <Vect3> DR, DRDR, V, VV, X, XX;
    Vect3 Start, End;

    Vect3 Vinf(globalSystem->uinf, globalSystem->vinf, globalSystem->winf);
    if (ProtoWake.size() > 0){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int) ProtoWake.size(); ++i) {

        ProtoWake[i]->GetCollocationPoint();
        Vect3 Pos = ProtoWake[i]->CollocationPoint->vP - CG.vP;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = ProtoWake[i]->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = ProtoWake[i]->Owner->CG.vV + Vrot;
        // 	Include freestream
        ProtoWake[i]->CollocationPoint->vV = Vinf - Vkin; // + AllBodyPanels[i]->CollocationPoint->vVfmm;
        //  Iterate over all wake panels
        for (int J = 0; J < (int) globalSystem->Bodies.size(); ++J)
            ProtoWake[i]->CollocationPoint->vV += globalSystem->Bodies[J]->GetVel(Pos + CG.vP);

    }

    for (int i = 0; i < ProtoWake.size(); ++i) {
        PANEL *temp_pan = ProtoWake[i];
        Vect3 S1,S2;
        if (!ProtoWake[i]->Neighb.B) {
            Start = S1 = .5 * (temp_pan->C4->vP + temp_pan->C1->vP);
            S2 = .5 * (temp_pan->C2->vP + temp_pan->C3->vP);
            XX.push_back(S1);
            XX.push_back(S2);
            DRDR.push_back((S2-S1));
            R.push_back((temp_pan->CollocationPoint->vP - Start).Mag());
            G.push_back(temp_pan->gamma);
            GP.push_back(temp_pan->gamma_prev);
            DR.push_back((.5*(temp_pan->C4->vP + temp_pan->C1->vP) - .5*(temp_pan->C2->vP + temp_pan->C3->vP)));
            V.push_back(temp_pan->CollocationPoint->vV);
            X.push_back(temp_pan->CollocationPoint->vP);
            while ((temp_pan->Neighb.T) && (temp_pan->Neighb.T!=temp_pan)) {
                temp_pan = temp_pan->Neighb.T;
                S1 = .5 * (temp_pan->C4->vP + temp_pan->C1->vP);
                S2 = .5 * (temp_pan->C2->vP + temp_pan->C3->vP);
                XX.push_back(S2);
                DRDR.push_back((S2-S1));
                R.push_back(R.back() + (temp_pan->CollocationPoint->vP - temp_pan->Neighb.B->CollocationPoint->vP).Mag());
                G.push_back(temp_pan->gamma);
                GP.push_back(temp_pan->gamma_prev);
                DR.push_back((.5*(temp_pan->C4->vP + temp_pan->C1->vP) - .5*(temp_pan->C2->vP + temp_pan->C3->vP)));
                V.push_back(temp_pan->CollocationPoint->vV);
                X.push_back(temp_pan->CollocationPoint->vP);
            }
            End = .5 * (temp_pan->C2->vP + temp_pan->C3->vP);
        }
    }

    GG.push_back(0);
    RR.push_back(0);
    for (int i = 0; i < V.size()-1; ++i){
    	GG.push_back(.5*(G[i]+G[i+1]));
    	RR.push_back(RR.back() + DRDR[i].Mag());
    }

    GG.push_back(0);
    RR.push_back(RR.back() + DRDR.back().Mag());


        int n = R.size();
        Array <Vect3> dOdt(n), UdelGamma(n), Om(n);
    {
        REAL r[n+1], g[n+1];
        for (int i = 0; i < R.size(); ++i)
            dOdt[i] = ((G[i] - GP[i])/dt)*DR[i];	//	the DR is to give the direction to the vorticity

        for (int i = 0; i < n+1; ++i)
        {
            r[i] = RR[i];
            g[i] = GG[i];
        }


        gsl_interp_accel *acc
                = gsl_interp_accel_alloc();
        gsl_spline *spline
                = gsl_spline_alloc(gsl_interp_akima,n+1);

        gsl_spline_init(spline, r, g, n+1);
        for (int i = 0; i < n; ++i) {
            UdelGamma[i] = gsl_spline_eval_deriv(spline, R[i], acc) * V[i] * globalSystem->GambitScale;
            Om[i] = dt*(UdelGamma[i] - dOdt[i]);
        }

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

    }
        Array <POINT*> new_points;
    //  Release a particle into the freestream
        for (int i = 0; i < X.size(); ++i)
            new_points.push_back(new POINT(X[i] + V[i]*dt,Om[i]));

        WakePoints.push_front(new_points);

    }
}
