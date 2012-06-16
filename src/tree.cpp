/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 28               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-09 15:0#$:  Date of last commit

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

/**************************************************************/
OCTREE::OCTREE() {
    Root = new Branch();
    Node::Root = Root;
    Root->UpdateMomentMults();
    FVMCell::InitMomsInds(globalSystem->MaxP);
    Branch::InitMomsInds(globalSystem->MaxP);
}

/**************************************************************/
OCTREE::~OCTREE() {
    delete Root;
}

/**************************************************************/
void OCTREE::InitVelsGetLaplacian() {
#ifdef TIME_STEPS
    long unsigned int t4 = ticks();
#endif
    Root->ApplyRecursively(&Branch::SetVelsZero, &FVMCell::SetVelsZero, &Node::DoNothing);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i) {
        AllCells[i]->GetLaplacian();
    }
#ifdef TIME_STEPS
    long unsigned int t5 = ticks();
    stringstream tmp;
    tmp << "InitVelsGetLaplacian     : " << double(t5 - t4) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void OCTREE::Prune() {
    Root->ApplyRecursively(&Node::Prune, &Node::DoNothing, &Node::DoNothing);
}

/**************************************************************/
void OCTREE::GetSRad() {
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->ReportSpectralRadius();
    //Root->ApplyRecursively(&Node::DoNothing, &FVMCell::ReportSpectralRadius, &Node::DoNothing);
}

/**************************************************************/
void OCTREE::ClearNodes() {
    Root->ClearChildren();
}

/**************************************************************/
void OCTREE::Reset() {
#ifdef TIME_STEPS
    long unsigned int t3 = ticks();
#endif
    //  Everything which changes the shape of the tree must be done recursively
    Root->ApplyRecursively(&Node::MarkWithoutLoad, &Node::MarkWithoutLoad, &Node::DoNothing);
    Root->ApplyRecursively(&Node::DoNothing, &Node::CheckLoad, &Node::DoNothing);

    if (globalTimeStepper->PruneNow) {
        Prune();
        globalTimeStepper->PruneNow = false;
    }

    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::CheckNeighbs, &Node::DoNothing);


    AllCells.clear();
    AllBranches.allocate(OCTREE_LEVS);
    BranchCount.assign(OCTREE_LEVS - 1, 0);
    Root->ApplyRecursively(&Branch::BranchCount, &Node::DoNothing, &Node::DoNothing);

    for (int i = 0; i < BranchCount.size(); ++i) {
        if (BranchCount[i] > 0)
            AllBranches[i].assign(BranchCount[i], NULL);
        BranchCount[i] = 0; //  Recycle for use as a pseudo iterator
    }

    CellCount = 0; //  This is used as a pseudo iterator for assigning into AllCells
    AllCells.assign(FVMCell::NumCells, NULL);
    Node::AllNodes.allocate(Node::NumNodes);
    Node::NodeCount = 0;
    Node::UpList.clear();
    Node::DownList.clear();
    Root->ApplyRecursively(&Node::DoNothing, &Node::ReList, &Node::ReList);
#ifdef TIME_STEPS
    long unsigned int t4 = ticks();
    stringstream tmp;
    tmp << "Calculate FMM: reset()   : " << double(t4 - t3) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void OCTREE::FVM() {
#ifdef TIME_STEPS
    long unsigned int t10 = ticks();
#endif
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
        AllCells[i]->O2UW();
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
#ifdef TIME_STEPS
    long unsigned int t11 = ticks();
    stringstream tmp;
    tmp << "FVM                      : " << double(t11 - t10) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void OCTREE::Integrate() {
#ifdef TIME_STEPS
    long unsigned int t12 = ticks();
#endif

#ifdef RECURSE
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::GetBEV, &Node::DoNothing);
    Root->ApplyRecursively(&Node::DoNothing, &FVMCell::Integrate, &Node::DoNothing); // <- this doesn't work in parallel
    Root->ApplyRecursivelyP(&Node::DoNothing, &Node::DoNothing, &Branch::SetFieldsZero);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i) {
//        AllCells[i]->GetVelTensor();
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
#ifdef TIME_STEPS
    long unsigned int t13 = ticks();
    stringstream tmp;
    tmp << "Integrate                : " << double(t13 - t12) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void OCTREE::GetVels() {
#ifdef TIME_STEPS
    long unsigned int t5 = ticks();
#endif
#ifdef RECURSE
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::SetVelsZero, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::PassMmnts2Prnt, &Node::PassMmnts2Prnt);
    Root->ApplyRecursivelyP(&Branch::GetVelField, &Node::DoNothing, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::CollapseVField, &Node::DoNothing);
    Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::SetVelsEqual, &Node::DoNothing);
#else

    //    cout << "Passing Moments to Parents" << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->PassMmnts2Prnt();

    //  Sweep moments up OCTREE (from cells at L12 -> root at L0)
    //    cout << "Sweeping Moments Up Tree" << endl;
    for (int mlev = AllBranches.size() - 1; mlev >= 0; --mlev) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i) {
            AllBranches[mlev][i]->PassMmnts2Prnt();
        }
    }
    //  Sweep velocity fields down OCTREE

    for (int mlev = 0; mlev < AllBranches.size(); ++mlev) {
        //  Inherit vel fields from parent
        //        cout << "Sweeping Velocity Fields Level " << mlev <<  endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->InheritVField();

        //  Add influence from neighbours
        //        cout << "Adding Influence of Neighbours Level " << mlev <<  endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->GetVelField(); //  This seems not to like being paralellised
    }

    //    cout << "Collapsing Velocity Fields Onto Cells" <<  endl;
    //  Collapse velocity fields onto children
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->CollapseVField();
#endif

    //    for (int i = 0; i < AllCells.size(); ++i)
    //        AllCells[i]->Report();
#ifdef TIME_STEPS
    long unsigned int t6 = ticks();
    stringstream tmp;
    tmp << "GetVels()                : " << double(t6 - t5) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
Vect3 OCTREE::TreeVel(Vect3 P) {
    OctreeCapsule C(P, Vect3(0, 0, 0), false);
    C.IP = true;
    Root->EvalCapsule(C);
    return Vect3(C.Velocity);
}
