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
    Array <Array <Branch*> > AllBranches;
    void ClearNodes();
    void AddToTree();
    void Reset();
    void FVM();
    void Integrate();
    void GetVels();
    void GetSRad();
    void Prune();
    void InitVelsGetLaplacian();
    Vect3 TreeVel(Vect3);
};


#endif	/* _TREE_HPP */
