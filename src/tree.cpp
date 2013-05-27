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



#include "tree.hpp"
Array <FVMCell*> OCTREE::CellsInOrder;
Array <Branch*> OCTREE::BranchesInDownOrder;
Array <Branch*> OCTREE::BranchesInUpOrder;
/**************************************************************/
OCTREE::OCTREE() {
    Root = new Branch();
    Node::Root = Root;
    Root->UpdateMomentMults();
    FVMCell::InitMomsInds(globalSystem->MaxP);
    Branch::InitMomsInds(globalSystem->MaxP);
    Root->SetUpISBIndices();
}

/**************************************************************/
OCTREE::~OCTREE() {
    delete Root;
}

/**************************************************************/
void OCTREE::ResetAllVelsAndFields() {
#ifdef TIME_STEPS
    long unsigned int t4 = ticks();
#endif
    Root->ApplyRecursively(&Branch::SetVelsZero, &FVMCell::SetVelsZero, &Node::DoNothing);

  
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        globalOctree->AllCells[i]->CheckActive();

    for (int mlev = 0; mlev < globalOctree->AllBranches.size(); ++mlev)
#ifndef NOFMM
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < globalOctree->AllBranches[mlev].size(); ++i) globalOctree->AllBranches[mlev][i]->SetFieldsZero();


    ResetAllVelsAndFields();
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->O2UW();
   
#endif
#ifdef TIME_STEPS
    long unsigned int t11 = ticks();
    stringstream tmp;
    tmp << "FVM                      : " << double(t11 - t10) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}
/**************************************************************/
void OCTREE::DiffuseAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->Diffuse();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::DiffuseXAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->DiffuseX();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::DiffuseYAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->DiffuseY();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::DiffuseZAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->DiffuseZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::StretchAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->GetISAGrads();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->Stretch();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::O2UWxAndAdvance(REAL dt) {
    #ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->GetISAVels();
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->O2UWx();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}

/**************************************************************/
void OCTREE::O2UWyAndAdvance(REAL dt) {
    #ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->GetISAVels();
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->O2UWy();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}

/**************************************************************/
void OCTREE::O2UWzAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->GetISAVels();
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->O2UWz();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
}
/**************************************************************/
void OCTREE::O2UWAndAdvance(REAL dt) {
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->GetISAVels();
#ifdef _OPENMP
#pragma omp parallel for
#endif   
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->SetVelsEqual();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->O2UW();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->AdvanceDt(dt);
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
string OCTREE::GetDirectVels() {
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->Velocity = Vect3(0.0);

#pragma omp parallel for
    for (int i = 0; i < AllCells.size(); ++i)
        for (int j = 0; j < AllCells.size(); ++j)
            AllCells[i]->Velocity += UTIL::globalDirectVel(AllCells[i]->Position - AllCells[j]->Position, AllCells[j]->Omega);
    
    string tmp("");
    return tmp;
}

/**************************************************************/
string OCTREE::GetRecursiveFMMVels() {
    stringstream tmp;
#ifdef TIME_STEPS
    long unsigned int tt1 = ticks();
#endif
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::SetVelsZero, &Node::DoNothing);
#ifdef TIME_STEPS
    long unsigned int tt2 = ticks();
#endif
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::PassMmnts2Prnt, &Node::PassMmnts2Prnt);
#ifdef TIME_STEPS
    long unsigned int tt3 = ticks();
#endif
    Root->ApplyRecursivelyP(&Branch::GetVelField, &Node::DoNothing, &Node::DoNothing);
#ifdef TIME_STEPS
    long unsigned int tt4 = ticks();
#endif
    Root->ApplyRecursivelyP(&Branch::InheritVField, &Node::CollapseVField, &Node::DoNothing);
#ifdef TIME_STEPS
    long unsigned int tt5 = ticks();
#endif
    Root->ApplyRecursivelyP(&Node::DoNothing, &FVMCell::SetVelsEqual, &Node::DoNothing);
#ifdef TIME_STEPS
    long unsigned int tt6 = ticks();
    tmp << "Reset Vels:       " << double(tt2 - tt1) / 1000.0 << endl;
    tmp << "ME:               " << double(tt3 - tt2) / 1000.0 << endl;
    tmp << "Get Vel field     " << double(tt4 - tt3) / 1000.0 << endl;
    tmp << "L2L & Kernel Exp.:" << double(tt5 - tt4) / 1000.0 << endl;
    tmp << "Set Vels Equal:   " << double(tt6 - tt5) / 1000.0 << endl;
#endif
    string output = tmp.str();
    return output;
}
/**************************************************************/
string OCTREE::GetPseudoRecursiveFMMVels() {
    stringstream tmp;
#ifdef TIME_STEPS
    long unsigned int tt1 = ticks();
#endif

    OCTREE::BranchesInDownOrder.clear();
    OCTREE::BranchesInUpOrder.clear();
    OCTREE::CellsInOrder.clear();
    Root->ApplyRecursively(&Branch::PutInOctreeDownList, &FVMCell::PutInOctreeCellList, &Branch::PutInOctreeUpList);
    cout << "here! " << OCTREE::CellsInOrder.size() << " " << AllCells.size() << endl;
#ifdef TIME_STEPS
    long unsigned int tt2 = ticks();
#endif

    for (int i = 0; i < CellsInOrder.size(); ++i) {
        CellsInOrder[i]->SetVelsZero();
        CellsInOrder[i]->PassMmnts2Prnt();
    }

    for (int i = 0; i < BranchesInUpOrder.size(); ++i)
        BranchesInUpOrder[i]->PassMmnts2Prnt();

    for (int i = 0; i < BranchesInDownOrder.size(); ++i)
        BranchesInDownOrder[i]->GetVelField();

    for (int i = 0; i < BranchesInDownOrder.size(); ++i)
        BranchesInDownOrder[i]->InheritVField();

    for (int i = 0; i < CellsInOrder.size(); ++i)
        CellsInOrder[i]->CollapseVField();

    for (int i = 0; i < CellsInOrder.size(); ++i)
        CellsInOrder[i]->SetVelsEqual();

#ifdef TIME_STEPS
    tmp << "Construct lists   " << double(tt2 - tt1) / 1000.0 << endl;
#endif
    string output = tmp.str();
    return output;
}
/**************************************************************/
string OCTREE::GetNonRecursiveFMMVels() {
//      Pass moments to parents
#ifdef TIME_STEPS
    long unsigned int tt1 = ticks();
#endif
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i){
        AllCells[i]->PassMmnts2Prnt();
    }

    //  Sweep moments up OCTREE (from cells at L12 -> root at L0)
#ifdef TIME_STEPS
    long unsigned int tt2 = ticks();
#endif
    for (int mlev = AllBranches.size() - 1; mlev >= 0; --mlev) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i) {
            AllBranches[mlev][i]->PassMmnts2Prnt();
        }
    }
    //  Sweep velocity fields down OCTREE
#ifdef TIME_STEPS
    long unsigned int tt3 = ticks();
    long unsigned int tL2L = 0, tM2L = 0;
#endif
    for (int mlev = 0; mlev < AllBranches.size(); ++mlev) {
        //  Inherit vel fields from parent
#ifdef TIME_STEPS
        long unsigned int ttL2L = ticks();
#endif
#ifdef _OPENMP
//#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->GetVelField(); //  This seems not to like being paralellised

#ifdef TIME_STEPS
        tL2L += ticks() - ttL2L;
        long unsigned int ttM2L = ticks();
#endif

        //  Add influence from colleagues
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < AllBranches[mlev].size(); ++i)
            AllBranches[mlev][i]->InheritVField();

#ifdef TIME_STEPS
        tM2L += ticks() - ttM2L;
#endif
    }
    //  Collapse velocity fields onto children
#ifdef TIME_STEPS
    long unsigned int tt4 = ticks();
#endif
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < AllCells.size(); ++i)
        AllCells[i]->CollapseVField();

#ifdef TIME_STEPS
    long unsigned int tt5 = ticks();
    stringstream tmp;
    tmp << "ME:               " << double(tt2 - tt1) / 1000.0 << endl;
    tmp << "M2M:              " << double(tt3 - tt2) / 1000.0 << endl;
    tmp << "M2L:              " << double(tM2L) / 1000.0 << endl;
    tmp << "L2L:              " << double(tL2L) / 1000.0 << endl;
    tmp << "Kernel Expansion: " << double(tt5 - tt4) / 1000.0 << endl;
#endif
    string output = tmp.str();
    return output;
}
/**************************************************************/
void OCTREE::GetVels() {

#ifdef TIME_STEPS
    long unsigned int t5 = ticks();
#endif
    //#define USE_DIRECT
#ifdef USE_DIRECT
    string time_info = GetDirectVels();
#else
#ifdef RECURSE
    string time_info = GetRecursiveFMMVels();
#else
    string time_info = GetPseudoRecursiveFMMVels();
#endif
#endif
#ifdef TIME_STEPS
    long unsigned int t6 = ticks();
    stringstream tmp;
    tmp << "GetVels()                : " << double(t6 - t5) / 1000.0 << endl;
    globalIO->step_data += tmp.str()  + time_info;
#endif
}

/**************************************************************/
Vect3 OCTREE::TreeVel(Vect3 P) {
    OctreeCapsule C(P, Vect3(0, 0, 0), false);
    C.IP = true;
    Root->EvalCapsule(C);
    return Vect3(C.Velocity);
}
/**************************************************************/
void BulkLoader(Array <Vect3> &X, Array <Vect3> &OM, Array <int> &BodyID)
{
    
    int S = OCTREE_SIZE / 2;
    
    Array <long unsigned int> IDs (X.size(),0);
    Array <Vect3> XCopy = X;
    while (S > 0) {
        S /= 2;
        for (int i = 0; i < X.size(); ++i) {
            int I = (int) (XCopy[i].x > 0), J = (int) (XCopy[i].y > 0), K = (int) (XCopy[i].z > 0);
            IDs[i] = ((((((IDs[i] << 1) + I) << 1) + J) << 1) + K);
            XCopy[i].x -= S * (2 * I - 1);
            XCopy[i].y -= S * (2 * J - 1);
            XCopy[i].z -= S * (2 * K - 1);
        }
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
