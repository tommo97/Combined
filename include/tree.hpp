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



#ifndef _TREE_HPP
#define	_TREE_HPP
#include "includes.hpp"
#include "types.hpp"
#include "node.hpp"
#include "branch.hpp"
#include "cell.hpp"
#include "system.hpp"

class OCTREE {
public:
    ~OCTREE();
    OCTREE();
    Branch *Root;
    int CellCount;              //  These are used for bookkeeping
    Array <int> BranchCount;    //  ditto

    Array <FVMCell*> AllCells;
    static Array <FVMCell*> CellsInOrder;
    Array <Array <Branch*> > AllBranches;
    static Array <Branch*> BranchesInDownOrder, BranchesInUpOrder;
    void ClearNodes();
    void AddToTree();
    void Reset();
    void FVM();
    void Integrate();
    void GetVels();
    string GetDirectVels();
    string GetRecursiveFMMVels();
    string GetPseudoRecursiveFMMVels();
    string GetNonRecursiveFMMVels();
    void GetSRad();
    void Prune();
    void ResetAllVelsAndFields();
    void BulkLoader(Array <Vect3> &, Array <Vect3> &, Array <int> &);
    
    
    void StretchAndAdvance(REAL);
    void DiffuseAndAdvance(REAL);
    void DiffuseXAndAdvance(REAL);
    void DiffuseYAndAdvance(REAL);
    void DiffuseZAndAdvance(REAL);
    void O2UWxAndAdvance(REAL);
    void O2UWyAndAdvance(REAL);
    void O2UWzAndAdvance(REAL);
    void O2UWAndAdvance(REAL);
    Vect3 TreeVel(Vect3);
};


#endif	/* _TREE_HPP */

