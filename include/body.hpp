/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 2                $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2009-10-28 20:1#$:  Date of last commit

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


#include "utils.hpp"

class BODY {
public:
    static Array <Vect3> CGS;
    Vect3 CG, CGo;
    static Array <Vect3> ATTITUDE;
    Vect3 Attitude;
    static Array <Vect3> VELOCITY;
    Vect3 Velocity;
    static Array <Vect3> RATES;
    static Array <string> NAMES;
    static Array <REAL> TimePrev;
    static Array <REAL> DPhi;
    static Array <REAL> DCp;
    static Array <REAL> DphiDCp;
    static Array <REAL> dDeltaCp_dDeltaPhi;
    static Array <REAL> LiftHist;
    static Array <PANEL*> AllProtoWakes;
    static Array <Vect3> VortonPositions;
    static Array <Vect3> VortonStrengths;
    static REAL Radius, TSR, RHO;
    string Name;
    static Array <REAL> Times;
    Vect3 Rmax;
#ifdef USEGSL
    gsl_matrix * gslA; //  A (doublet) influence coefficient matrix for this body
    gsl_matrix * gslB; //  B (source) influence coefficient matrix for this body
    gsl_permutation * gslPermut;
    gsl_vector * gslMu;
    gsl_vector * gslSigma;
    gsl_vector * gslRHS;
#endif
    static REAL Time;
    static Array <Vect3> TorqueHist;
    static Array <Vect3> ForceHist;
    static Array < Array < Array < int > > > Surfaces;
    static Array < Array < Array < PANEL*> > > ptSurfaces;
    static REAL GambitScale;
    static Array < Array <Array <int> > > PtIDS;
    static Array < Array <Array <int> > > TipInboardUSIDS;
    static Array < Array <Array <int> > > TipInboardLSIDS;
    static Array < Array <Array <int> > > TipOutboardUSIDS;
    static Array < Array <Array <int> > > TipOutboardLSIDS;
    static Array < Array < Array < int > > > InnerTipUSPanelIDS, OuterTipUSPanelIDS;
    static Array < Array < Array < int > > > InnerTipLSPanelIDS, OuterTipLSPanelIDS;
    static Array <Vect3> AllBodyPoints;
    static Array < Array <REAL> > A;
    static Array < Array <REAL> > B;
    static Array < Array <REAL> > C;
    static Array < Array <REAL> > D;
    static Array < Array <REAL> > LocalChordRadius;
    static Array <REAL> RHS;
    static Array <REAL> Mu;
    static Array <REAL> Sigma;
    static int NumBodies;
    static int NumFaces;
    static Array <BODY*> Bodies;
    static bool LiftingLineMode;
    static Array <PANEL*> AllBodyFaces;
    static int BodyPanelIDCounter, BodyPointIDCounter;
    static void BodySubStep(REAL delta_t, int n_steps);
    
    void MakeTmpWake();

    Vect3 BodyRates, BodyAngles, EulerRates, EulerAngles;

	Array <Vect3> AngleHist;

    Array < Array <REAL> > localA;
    Array < Array <REAL> > localB;
    Array <REAL> localSigma, localRHS, localMu;

    Vect3 TRANS[3];
    Array <PANEL> Faces;
    Array <Array <PANEL> > ProtoWakes, ProtoWakes0;
    Array < Array < Array < Vect3 > > > VortonX, VortonOM, VortonVel;
    Array < Array < REAL > > VortonAge;
    Array < int  > VXSizes, VOMSizes, VVelSizes;
    Array <PANEL*> BoundaryFaces, FirstProtoWakes;

    Array < Array < Array < Vect3 > > > WakePoints;
    Array < Array < Array < REAL > > > WakeGamma;
    Array < Array < Array < PANEL > > > WkTmp;

    BODY(Vect3 Or, Vect3 At, Vect3 Ve, Vect3 Ra, string Na) : CG(Or), Attitude(At), Velocity(Ve), Name(Na), BodyRates(Ra){
        EulerAngles = Attitude; //  psi theta phi
        SetEulerTrans();
        BodyAngles = Vect3(0.0); //  omx omy omz
        CGo = CG;
        SetEulerTrans();
        GetEulerRates();
        BODY::LiftingLineMode = false;
//        WkTmp = Array < Array < Array <PANEL> > > ();
    };

    
    static void ReadNeuGetBodies(string neu_file, string name, Vect3 dpos, Vect3 cg, Vect3 vel, Vect3 att, Vect3 rates, bool flip, int plane);
    void GetPanels(Array <PANEL> &PANELS);
    void GetEulerRates();
    void GetBodyRates();
    void SetEulerTrans();
    void MoveBody(REAL);
    void WakeToVortons();
    static bool IPKC;
    static void SetUpProtoWakes(REAL);
    static void PollFaces();
    static void UpdateGlobalInfluenceMatrices();
    static void SplitUpLinearAlgebra();
    static void SetUpInfluenceMatrices();
    Vect3 GetVel(Vect3);
    void MakeWake();
    static void GetLinearRHS();
    static void GetNonLinearRHS();
    static void LinAlg();
    Vect3 GetWakeVel(Vect3);
    REAL GetWakePhi(Vect3);

};
#endif
