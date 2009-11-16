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


#ifndef PANEL_HPP
#define PANEL_HPP
#include "includes.hpp"
#include "types.hpp"
#include "body.hpp"




class PANEL {
public:
    long unsigned int ID;
    Vect3 Centroid, TRANS[3], Xcb[4], Vfmm, Vkin; //  Transformation matrix from BODY to PANEL axis; Position of corners
    REAL *sigma, *mu, gamma, Area, MaxDiagonal, mu_prev, gamma_prev, alpha, chord, Uinf, Vn, Cpress;
    POINT *C1, *C2, *C3, *C4, *CollocationPoint, *edgeX1, *edgeX2;
    Vect4 DX, DY, M, D;
    bool isBound, isTop;
    BODY *Owner;
    Vect3 tang1, tang2, CPVel;
    Array <Vect3> VelInterp;

    PanelNeighbSet <PANEL*> Neighb, NeighbNeighb;

    PanelNeighbSet <int> RecipNeighb;

    PanelNeighbSet <bool> doSide;

    PanelNeighbSet <REAL> gammaSide, Ds, Theta;
        
    PANEL *OtherBoundarySurface, *Wake, *Shedder;

    ~PANEL();
    PANEL() {}
    PANEL(POINT *x1, POINT *x2, POINT *x3, POINT *x4);
    PANEL(const PANEL &);
    void CheckNeighb(PANEL *);
    void GetCollocationPoint();
    void GetEdgeInfo();
    void DecomposeIntoBlobs();
    void UpdateCollocationPoint() {
        GetCollocationPoint();
        GetEdgeInfo();
    }
    void WakeNeighbSet();
    void SourceDoubletPotential(POINT *, REAL, REAL, REAL []);
    void PrintCorners();
    void PrintCollocationPointNormal();
    void UpdatePosition(REAL);
    void CalcNormalFlow();
    void GetLatticeVels();
    void PrintCollocationPoint();
    void GetNewGlobalPosition();
    void GetCp(REAL);
    Vect3 LineVelocity(Vect3 lineStart, Vect3 lineEnd, Vect3 pTarget, const REAL gamma_in);
    Vect3 SourceVel(Vect3 pTarget);
    Vect3 WakePanelVelocity(Vect3);
    Vect3 BodyPanelVelocity(Vect3);
    Vect3 BodyPanelVelocity(Vect3, REAL);

    class NoIntersect {
    };
};
#endif
