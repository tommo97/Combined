/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 2                $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-10-28 20:1#$:  Date of last commit

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
    
    static unsigned long int NumCells;

    NeighbSet <Vect3> FaceVels;

    Array <Vect3> Deriv, IPs, IPVels;

    int age;
    
    REAL Phi;

    Vect3 VelTensor[3];

    Array <Vect3> Laplacian;

    void CheckActive();
    
    void Reset(){};

    void ReportToIO();

    void CountCells();

    void ImageToIO();

    void vRemoveFromNeighbs() {
    };

    void vEvalCapsule(OctreeCapsule &c);

    void vCheckNeighbs();
    
    void vReList();

    void vApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    void vApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up);

    void TolDown();

    void SetVelsZero();

    void Integrate();

    Vect3 ReturnSpectralRadius();

    void GetLaplacian();

    void GetVelTensor();

    void ReportSpectralRadius();
    
    void vCollapseVField();

    void vPassMmnts2Prnt();

    void CollapseToIP(OctreeCapsule &);

    void SetVelsEqual();

    void Report();

    void ReList();

    void O1UW();

    void O2UW();

    void GetBEV();

    void MUSCL();

protected:
    virtual ~FVMCell();
};

#endif /* if !defined(FVMCell_INCL) */
