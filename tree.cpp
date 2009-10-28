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



#include "tree.hpp"

OCTREE::OCTREE() {
    Root = new Branch();
    Node::Root = Root;
    Root->UpdateMomentMults();
}

OCTREE::~OCTREE()
{
    delete Root;
}

void OCTREE::InitVelsGetLaplacian() {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i){
        AllCells[i]->SetVelsZero();
        AllCells[i]->GetLaplacian();
    }
}

void OCTREE::Prune() {
    Root->ApplyRecursively(&Node::DoNothing, &Node::DoNothing, &Node::Prune);
}

void OCTREE::GetSRad() {
    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::ReportSpectralRadius, &Node::DoNothing);
}

void OCTREE::ClearNodes() {
    Root->ClearChildren();
}

void OCTREE::Reset() {
    Prune();
    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::CheckNeighbs, &Node::DoNothing);


    AllCells.clear();
    AllBranches.allocate(12);
    BranchCount.assign(11,0);
    Root->ApplyRecursively(&Branch::BranchCount, &Node::DoNothing, &Node::DoNothing);

    for (int i = 0; i < BranchCount.size(); ++i){
        if (BranchCount[i]>0)
            AllBranches[i].assign(BranchCount[i],NULL);
        BranchCount[i] = 0; //  Recycle for use as a pseudo iterator
    }
    
    CellCount = 0;  //  This is used as a pseudo iterator for assigning into AllCells
    AllCells.assign(globalNum_FVMCELLS, NULL);
    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::ReList, &Branch::ReList);
}

void OCTREE::FVM() {
#ifdef RECURSE
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::O2UW, &Node::DoNothing);
#else

    #ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i) {
        
        AllCells[i]->GetBEV();
    }

//    if (globalTimeStepper->t > 2)
//    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->MUSCL();
//    }
//    else
//    {
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for (int i = 0; i < AllCells.size(); ++i)
//        AllCells[i]->O1UW();
//    }
#endif
}

void OCTREE::Integrate() {
#ifdef RECURSE
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::GetBEV, &Node::DoNothing);
    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::Integrate, &Node::DoNothing); // <- this doesn't work in parallel
    Root->ApplyRecursivelyP(&Node::DoNothing, &Node::DoNothing, &Branch::SetFieldsZero);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i) {
        AllCells[i]->GetVelTensor();
        AllCells[i]->Integrate();
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->CheckActive();
    
    for (int mlev = 0; mlev < AllBranches.size(); ++mlev)
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i) AllBranches[mlev][i]->SetFieldsZero();

#endif
}

void OCTREE::GetVels() {
#ifdef RECURSE
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::SetVelsZero, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::PassMmnts2Prnt, &Node::PassMmnts2Prnt);
    Root->ApplyRecursivelyP(&Branch::GetVelField, &Node::DoNothing, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::CollapseVField, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::SetVelsEqual, &Node::DoNothing);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->PassMmnts2Prnt();

    //  Sweep moments up OCTREE (from cells at L12 -> root at L0)
    for (int mlev = AllBranches.size() - 1; mlev >= 0; --mlev) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i){
            AllBranches[mlev][i]->PassMmnts2Prnt();
        }
    }
    //  Sweep velocity fields down OCTREE
    for (int mlev = 0; mlev < AllBranches.size(); ++mlev) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->InheritVField();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->GetVelField(); //  This seems not to like being paralellised
    }

    //  Collapse velocity fields onto children
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->CollapseVField();
#endif

//    for (int i = 0; i < AllCells.size(); ++i)
//        AllCells[i]->Report();
}

Vect3 OCTREE::TreeVel(Vect3 P) {
    OctreeCapsule C(P, Vect3(0,0,0), false);
    C.IP = true;
    Root->EvalCapsule(C);
    return Vect3(C.Velocity);
}