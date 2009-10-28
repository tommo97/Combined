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


#ifndef BODY_HPP
#define BODY_HPP

#include "includes.hpp"
#include "types.hpp"
#include "system.hpp"
#include "panel.hpp"



class BODY
{
  public:
    ~BODY();
    SYSTEM *globalSystem;
    BODY(SYSTEM *sys) {globalSystem = sys;}
    BODY(Vect3, Vect3, Vect3, Vect3, string, SYSTEM *sys);
    int BodyID;
    POINT CG;   //  ie position of local axis origin in global (inertial) axis x0,y0,z0
    POINT Origin;	// 	PV to the global origin in local axis.....
    Vect3 EulerAngles;		//  ie phi, theta, psi
    Vect3 TRANS[3];  
    Vect3 EulerRates;		//  ie phi_dot, theta_dot, psi_dot
    Vect3 BodyRates; 		//  ie P,Q,R
    Vect3 Rmax;                 //  Point on body furthest from origin
    Array <Vect3> ForceRecords;
    Array <Vect3> TorqueRecords;
    deque < Array <POINT*> > WakePoints;
    Array <POINT*> ProtoWakePoints;
    Vect3 vFORCE;
    Vect3 vTORQUE;
    void BodyVelocity(REAL  P[], REAL V[]);     //  Velocity V at point P due to all elements on this body
    void GetEulerRates();
    void GetBodyRates();
    void SetEulerTrans();
    void Earth2Body(REAL [], REAL []);
    void Body2Earth(REAL [], REAL []);
    void MoveBody(REAL);
    void InitNascentWake(REAL);
    void SortWake(REAL);
    void GetPanels(Array <int> &GROUP, Array <PANEL*> &PANELS);
    void PrintSurface();
    void PrintBoundary();
    void PrintWake();
    void WriteSurface(ostream& out_stream);
    void WriteWake(ostream& out_stream);
    void GetNewWakePanels(BODY &);
    void DissolveWake(REAL);
    Vect3 GetWakeVel(Vect3 Target);
    void PrintVels();   
    Vect3 GetVel(Vect3 Target);
    string Name;
    int NumFaces;   //  Number of faces (ie panels) on body
    Array <PANEL*> Faces, BoundaryFaces, ProtoWake;
    int nSpanPanels;
    deque <Array <PANEL*> > WakeGlobal;
    gsl_matrix * A;     //  A (doublet) influence coefficient matrix for this body
    gsl_matrix * B;     //  B (source) influence coefficient matrix for this body
    gsl_permutation * p;
    gsl_vector * Mu;
    gsl_vector * Sigma;
    gsl_vector * RHS;

};

#endif
