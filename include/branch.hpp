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



#ifndef _BRANCH_HPP
#define	_BRANCH_HPP

#include "node.hpp"
#include "cell.hpp"

/**************************************************************/
class Branch : public Node {
public:

    virtual ~Branch();

    Branch();

    Branch(Node *parent, int i, int j, int k);
    
    static int NumBranches;

    static Array < Array <int> > srcMomentInds, trgMomentInds, srcInheritInds, trgInheritInds;
    
    static int MomentIndsSize, InheritIndsSize;
    
    static Array < Array < Array <REAL> > > BinomMults, InheritMults;
    
    static void InitMomsInds();
    
    void PutInOctreeDownList();
    
    void PutInOctreeUpList();
    
    void vRemoveFromNeighbs() {};

    void vCheckNeighbs() {};

    void vCollapseVField() {};
    
    void vMakeNodeAtTrans(Array < int > &trn) ;
    
    void vReList();

    void vEvalCapsule(OctreeCapsule&);

    void vApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    void vApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    JaggedArray <Vect3> Moments;
    JaggedArray <Vect3> VelField;
//    Array <Vect3> InlineMoments;

    void BranchCount();
    void InheritVField();
    void vPassMmnts2Prnt();
    void GetVelField();
    void SetFieldsZero();
    void ReList();
    void SetVelsZero()
    {
    	Velocity = 0.;
        PanelVel = 0.;
    }
};



#endif	/* _BRANCH_HPP */

