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


#ifndef FVMCell_INCL
#define FVMCell_INCL

#include "includes.hpp"
#include "types.hpp"
#include "array.hpp"
#include "vect34.hpp"
#include "node.hpp"
#include "branch.hpp"
#include "system.hpp"

/**************************************************************/
class FVMCell : public Node {
public:

    FVMCell();

    FVMCell(Node *parent, int i, int j, int k);
    
    static Array <FVMCell*> AllCells;
    
    static Array < Array <Vect3> > CellDerivs0, CellDerivs1, CellDerivs;
    
    static Array < Array <int> > MomentInds;
    
    static int MomentIndsSize;
    
    static Array < Array <REAL> > OffsetPows;
    
    static int NumCells;

    NeighbSet <Vect3> FaceVels;

    Vect3  VelHold, OmegaHold;
    
    void AdvanceDt(REAL);
    
    void GetISAVels();
    
    void GetISBVels();

    void GetISAGrads();

    void GetISBGrads();
    
    static REAL (*limiter) (REAL);

    inline static REAL O1UW(REAL r) {
        return 0.0; //O(1)UW
    }

    inline static REAL minmod(REAL r) {
        return max(REAL(0.0), min(REAL(1.0), r));
    };

    inline static REAL superbee(REAL r) {
        return max(max(min(REAL(2.0) * r, REAL(1.0)), min(r, REAL(1.0))), REAL(0.0));
    };

    inline static REAL vanLeer(REAL r) {
        return (r + fabs(r)) / (1.0 + fabs(r));
    };

    inline static REAL Koren(REAL r) {
        return max(0.0, min(min(2.0 * r, (2.0 + r) / 3.0), 2.0));
    };

    inline static REAL UMIST(REAL r) {
        REAL M = 2.0, kappa = 0.5;
        return max(0.0, min(M, min(min(0.5 * (1. + kappa) * r + 0.5 * (1. - kappa), 0.5 * (1. - kappa) * r + 0.5 * (1. + kappa)), 2.0 * r))); //UMIST
        //return max(0.0,min(2.0*r,min(min(0.25+0.75*r,0.75+0.25*r),2.)));
    }

    inline static REAL  MUSCL(REAL r) {
        REAL M = 2.0, kappa = 0.0;
        return max(0.0, min(M, min(min(0.5 * (1. + kappa) * r + 0.5 * (1. - kappa), 0.5 * (1. - kappa) * r + 0.5 * (1. + kappa)), 2.0 * r))); //MUSCL
    }

    inline static REAL  CADA(REAL r) {
        REAL ar = fabs(r);
        return max(-3.0 / ar, min((2.0 + r) / 3.0, min((2.0 * ar) / (0.6 + ar), 5.0 / ar))); //CADA
    }

    inline static REAL  MC(REAL r) {
        return max(0.0, min(min(2.0 * r, (1.0 + r) / 2.0), 2.0)); //
    }



    inline static Vect3 flim(Vect3 r) {
        return Vect3(limiter(r.x), limiter(r.y), limiter(r.z));
    };


    
//    REAL Phi;                         // unused?

    Array <Vect3> VelGrads, VelGradsHold;

    static void InitMomsInds();

    void CheckActive();
    
    void NormaliseObliterate();
    
    void Reset(){};

    void ReportToIO();

    void CountCells();

    void ImageToIO();

    void vRemoveFromNeighbs() {
    };

    void vEvalCapsule(OctreeCapsule &c);

    void vCheckNeighbs();
    
    void vMakeNodeAtTrans(Array < int > &){};
    
    void vReList();

    void vApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    void vApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    void TolDown();

    void SetVelsZero();
    
    void PutInOctreeCellList();

    void Integrate();
    
    void Stretch();
    
    Vect3 Stretch(int);
    
    void Diffuse();
    
    Vect3 Diffuse(int);
    
    void DiffuseX();
    
    void DiffuseY();
    
    void DiffuseZ();
    
    Vect3 ReturnSpectralRadius();

    void ReportSpectralRadius();
    
    void UpdateVelsFromDeltaOmega();
    
    void UpdateGradsFromDeltaOmega();
    
    void vCollapseVField();

    void vPassMmnts2Prnt();

    void CollapseToIP(OctreeCapsule &);

    void SetVelsEqual();

    void Report();

    void ReList();

    void O1UW();
    Vect3 O1UW(int);
    void O2UW();
    Vect3 O2UW(int);
    void O2UWx();
    void O2UWy();
    void O2UWz();
#ifdef USE_MUSCL
    void GetBEV();
    void MUSCL();
#endif
protected:
    virtual ~FVMCell();
};

#endif /* if !defined(FVMCell_INCL) */
