/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 36               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-16 17:3#$:  Date of last commit

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


#ifndef PANEL_HPP
#define PANEL_HPP
#include "includes.hpp"
#include "types.hpp"
#include "body.hpp"

#include "utils.hpp"


class PANEL {
public:
    static REAL FarField;

    Vect3 C1, C2, C3, C4, Centroid, CollocationPoint, Normal, Eta, Epsilon, TRANS[3], edgeX1, edgeX2, Vfmm, Vfmm0, Vfmm1, Vkin, VWake, VCentroid, dF;
    Array <Vect3> Xcb;
    Vect4 DX, DY, M, D;
    Vect3 C1o, C2o, C3o, C4o;
    REAL Phi, Sigma, Mu, Gamma, MuPrev, GammaPrev, Cloc, Rloc;
    REAL Area, MaxDiagonal, Vn;
    BODY *Owner;
    bool isBound, isWake, isTop;
    int BoundBC, Group, ID;
    PANEL *OtherBoundarySurface, *AttachedProtoWake, *SheddingSurface;
    PanelNeighbSet <PANEL*> Neighb;
    PanelNeighbSet <REAL> Theta;
    Array <REAL> PhiPrev;

    void CN() {

//        if (!Neighb.L)
//            cout << "%\tMissing Neighb.L" << endl;
//
//        if (!Neighb.R)
//            cout << "%\tMissing Neighb.R" << endl;
//
//        if (!Neighb.T)
//            cout << "%\tMissing Neighb.T" << endl;
//
//        if (!Neighb.B)
//            cout << "%\tMissing Neighb.B" << endl;

        SURFW(C1,C2,C3,C4,Phi,cout);



    }

    PANEL()
    {
        C1 = C2 = C3 = C4 = C1o = C2o = C3o = C4o = edgeX1 = edgeX2 = dF = Vect3(0.0);
        Sigma = Mu = 0.0;
        Gamma = Phi = 0.0;
        Area = MaxDiagonal = Sigma = Mu = Gamma = Vn = MuPrev = GammaPrev = 0.0;
        Centroid = CollocationPoint = Vect3(.25 * (C1 + C2 + C3 + C4));
        Owner = NULL;
        OtherBoundarySurface = AttachedProtoWake = SheddingSurface = NULL;
        Neighb = NULL;
        Theta = 0.0;
        isBound = isTop = isWake = false;
        BoundBC = ID = -1;
        TRANS[0] = TRANS[1] = TRANS[2] = Vect3(0.);
        PhiPrev = Array <REAL> (4,0.0);
        Cloc = Rloc = 0.0;
        Vfmm = Vect3(0.);
    }

    PANEL(Vect3 P1, Vect3 P2, Vect3 P3, Vect3 P4) {
        C1  = P1; C2  = P2; C3  = P3; C4  = P4;
        C1o = C1; C2o = C2; C3o = C3; C4o = C4;
        edgeX1 = edgeX2 = dF = Vect3(0.0);
        Sigma = Mu = 0.0;
        Gamma = Phi = Sigma = Mu = 0.0;
        Area = MaxDiagonal = Gamma = Vn = MuPrev = GammaPrev = 0.0;
        Centroid = Vect3(.25 * (C1 + C2 + C3 + C4));
        Owner = NULL;
        OtherBoundarySurface = AttachedProtoWake = SheddingSurface = NULL;
        Neighb = NULL;
        Theta = 0.0;
        isBound = isTop = isWake = false;
        BoundBC = ID = -1;
        GetNormal();
        PhiPrev = Array <REAL> (4,0.0);
        Cloc = Rloc = 0.0;
        Vfmm = Vect3(0.);
    }

    PANEL(const PANEL &C) {
        C1 =  C.C1; C2 =  C.C2; C3 =  C.C3; C4 =  C.C4;
        C1o = C.C1; C2o = C.C2; C3o = C.C3; C4o = C.C4;
        Sigma = C.Sigma;
        Mu = C.Mu;
        Gamma = C.Gamma;
        Phi = C.Phi;
        Area = C.Area;
        MaxDiagonal = C.MaxDiagonal;
        Gamma = C.Gamma;
        GammaPrev = C.GammaPrev;
        Vn = C.Vn;
        MuPrev = C.MuPrev;
        Centroid = C.Centroid;
        TRANS[0] = C.TRANS[0];
        TRANS[1] = C.TRANS[1];
        TRANS[2] = C.TRANS[2];
        Owner = C.Owner;
        isBound = C.isBound;
        BoundBC = C.BoundBC;
        isTop  = C.isTop;
        Neighb = C.Neighb;
        Theta = C.Theta;
        OtherBoundarySurface = C.OtherBoundarySurface;
        AttachedProtoWake = C.AttachedProtoWake;
        edgeX1 = C.edgeX1;
        edgeX2 = C.edgeX2;
        SheddingSurface = C.SheddingSurface;
        Xcb = C.Xcb;
        Group = C.Group;
        isWake = C.isWake;
        GetNormal();
        PhiPrev = C.PhiPrev;
        VWake = C.VWake;
        Cloc = C.Cloc;
        Rloc = C.Rloc;
        dF = C.dF;
        ID = C.ID;
        Vfmm = C.Vfmm;
    }


    void GetNewGlobalPosition();
    void GetNormal();
    void CheckNeighb(PANEL *Face);
    void SourceDoubletPotential(Vect3 IP, REAL Mu, REAL Sigma, REAL Phi[]);
    void GetEdgeInfo();
    Vect3 BodyPanelVelocity(Vect3);
    Vect3 WakePanelVelocity(Vect3);
    Vect3 SourceVel(Vect3);
    REAL GetCp();
    static Vect3 LineVelocity(Vect3 &lineStart, Vect3 &lineEnd, Vect3 &pTarget, REAL gamma_in);
    static void SourceDoubletPotential(PANEL*, PANEL*, REAL &, REAL &, int, int);
    static void DoubletPotential(PANEL *source, Vect3, REAL &PhiDoublet, int i, int j);
    static void PointDoublet(Vect3 X0, Vect3 XP, Vect3 &V, Vect3 &Mu, REAL &Phi);
    static void PointSource(Vect3 X0, Vect3 XP, Vect3 &V, REAL &Sigma, REAL &Phi);
    static void SourceDoubletPotential(PANEL *source, Vect3 target, REAL &PhiDoublet, REAL &PhiSource, int i, int j);
    void DivPanel(int n, Array < Array < Vect3 > > &CP, Array < Array < Vect3 > > &N, Array < Array < REAL > > &A);

    void SubPan(int n, Vect3 P, REAL  MuNormal, REAL Sigma, REAL &PhiD, REAL &PhiS, Vect3 &V);
    void LinearSubPan(PANEL*, int n,  REAL Mu1, REAL Mu2,  REAL &PhiD, Vect3 &V);
};

#endif
