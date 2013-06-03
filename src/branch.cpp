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


#include "branch.hpp"
unsigned long int Branch::NumBranches = 0;
Array < Array <int> > Branch::srcMomentInds, Branch::trgMomentInds;
Array < Array < Array <REAL> > > Branch::BinomMults(OCTREE_LEVS);
Array < Array <int> > Branch::srcInheritInds, Branch::trgInheritInds;
Array < Array < Array <REAL> > > Branch::InheritMults(OCTREE_LEVS);
int Branch::MomentIndsSize = 0, Branch::InheritIndsSize = 0;

/**************************************************************/
Branch::~Branch() {
    Branch::NumBranches--;
}

/**************************************************************/
#ifdef NOFMM
Branch::Branch() : Node() {
#else
Branch::Branch() : Node(), Moments(SYSTEM::MaxP), VelField(SYSTEM::MaxP) {
#endif
    Branch::NumBranches++;
}

/**************************************************************/
#ifndef NOFMM
Branch::Branch(Node *parent, int i, int j, int k) : Node(parent, i, j, k), Moments(SYSTEM::MaxP), VelField(SYSTEM::MaxP) {
#else
Branch::Branch(Node *parent, int i, int j, int k) : Node(parent, i, j, k) {
#endif
    Branch::NumBranches++;
}

/**************************************************************/
void Branch::InitMomsInds() {

    Branch::BinomMults.allocate(OCTREE_LEVS);
    Branch::srcMomentInds.clear();
    Branch::trgMomentInds.clear();

    Branch::InheritMults.allocate(OCTREE_LEVS);
    Branch::srcInheritInds.clear();
    Branch::trgInheritInds.clear();

    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                for (int n1 = k1; n1 < SYSTEM::MaxP; n1++)
                    for (int n2 = k2; n1 + n2 < SYSTEM::MaxP; n2++)
                        for (int n3 = k3; n1 + n2 + n3 < SYSTEM::MaxP; n3++) {
                            Array <int> inds(3);
                            inds[0] = k1;
                            inds[1] = k2;
                            inds[2] = k3;
                            Branch::trgInheritInds.push_back(inds);
                            inds[0] = n1;
                            inds[1] = n2;
                            inds[2] = n3;
                            Branch::srcInheritInds.push_back(inds);
                        }

    Branch::InheritIndsSize = Branch::srcInheritInds.size();

    for (int mlev = 1; mlev < OCTREE_LEVS; ++mlev) {
        Branch::InheritMults[mlev].allocate(8);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
                        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                                for (int n1 = k1; n1 < SYSTEM::MaxP; n1++)
                                    for (int n2 = k2; n1 + n2 < SYSTEM::MaxP; n2++)
                                        for (int n3 = k3; n1 + n2 + n3 < SYSTEM::MaxP; n3++) {
                                            Branch::InheritMults[mlev][Node::Indxs[i][j][k]].push_back(Node::InhrtMlt[mlev][i][j][k][k1][k2][k3][n1][n2][n3]);
                                        }
    }

    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                for (int n1 = 0; n1 <= k1; n1++)
                    for (int n2 = 0; n2 <= k2; n2++)
                        for (int n3 = 0; n3 <= k3; n3++) {
                            Array <int> inds(3);
                            inds[0] = k1;
                            inds[1] = k2;
                            inds[2] = k3;
                            Branch::trgMomentInds.push_back(inds);
                            inds[0] = k1 - n1;
                            inds[1] = k2 - n2;
                            inds[2] = k3 - n3;
                            Branch::srcMomentInds.push_back(inds);

                        }

    Branch::MomentIndsSize = Branch::srcMomentInds.size();

    for (int mlev = 1; mlev < OCTREE_LEVS; ++mlev) {
        Branch::BinomMults[mlev].allocate(8);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
                        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                                for (int n1 = 0; n1 <= k1; n1++)
                                    for (int n2 = 0; n2 <= k2; n2++)
                                        for (int n3 = 0; n3 <= k3; n3++) {
                                            Branch::BinomMults[mlev][Node::Indxs[i][j][k]].push_back(Node::BinomMlt[mlev][i][j][k][k1][k2][k3][n1][n2][n3]);
                                        }
    }

};

/**************************************************************/
void Branch::BranchCount() {
    globalOctree->BranchCount[m]++;
}
/**************************************************************/
void Branch::PutInOctreeUpList()
{
    OCTREE::BranchesInUpOrder.push_back(this);
}
/**************************************************************/
void Branch::PutInOctreeDownList()
{
    OCTREE::BranchesInDownOrder.push_back(this);
}
/**************************************************************/
void Branch::vReList() {
    globalOctree->AllBranches[m][globalOctree->BranchCount[m]] = this;
    globalOctree->BranchCount[m]++;
}

/**************************************************************/
void Branch::vMakeNodeAtTrans(Array < int > &trn) {

    int i = Node::deREF[trn.front()][0];
    int j = Node::deREF[trn.front()][1];
    int k = Node::deREF[trn.front()][2];


    if (!Children[i][j][k]) {
        if (trn.size() == 1) {
            Children[i][j][k] = new FVMCell(this, i, j, k);
            return;
        } else {
            Children[i][j][k] = new Branch(this, i, j, k);
        }
    }
    trn.pop_front();
    Children[i][j][k]->MakeNodeAtTrans(trn);
}

/**************************************************************/
void Branch::SetFieldsZero() {
    VelField = Vect3(0.);
    Moments = Vect3(0.);
}

/**************************************************************/
void Branch::vPassMmnts2Prnt() {
    if ((m > 0) && (HasLoad)) {
#ifdef USE_ROLLED_LOOPS

        for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
            for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                    for (int n1 = 0; n1 <= k1; n1++)
                        for (int n2 = 0; n2 <= k2; n2++)
                            for (int n3 = 0; n3 <= k3; n3++)
                                (static_cast<Branch*> (Parent))->Moments[k1][k2][k3] += BinomMlt[m][x][y][z][k1][k2][k3][n1][n2][n3] * Moments[k1 - n1][k2 - n2][k3 - n3];


#else
        for (int i = 0; i < Branch::MomentIndsSize; ++i)
            (static_cast<Branch*> (Parent))->Moments[Branch::trgMomentInds[i][0]][Branch::trgMomentInds[i][1]][Branch::trgMomentInds[i][2]]
                += Branch::BinomMults[m][indx][i] * Moments[Branch::srcMomentInds[i][0]][Branch::srcMomentInds[i][1]][Branch::srcMomentInds[i][2]];


#endif
        Parent->HasLoad = true;
    }

}

/**************************************************************/
void Branch::vEvalCapsule(OctreeCapsule &C) {

    C.S /= 2;
    int i = (int) (C.Position.x > 0), j = (int) (C.Position.y > 0), k = (int) (C.Position.z > 0);

    C.Position.x -= C.S * (2 * i - 1);
    C.Position.y -= C.S * (2 * j - 1);
    C.Position.z -= C.S * (2 * k - 1);

    if ((!Children[i][j][k]) && (!C.checkExist)) {
        if (C.S > 0) {
            Branch *child = new Branch(this, i, j, k);
            Children[i][j][k] = child;
            if (C.IP) {
                child->InheritVField();
                child->GetVelField();
            }
        } else {
            Children[i][j][k] = new FVMCell(this, i, j, k);
        }
        Children[i][j][k]->GetISA();
        Children[i][j][k]->GetISB();
        if (C.toMonitor)
            Children[i][j][k]->toMonitor = true;
        // the above line means that if the branch/cell is being made by a monitoring OctreeCapsule, then it can safely be deleted without affecting the FVM/FMM calcs....
    }
    
    if (Children[i][j][k])
        Children[i][j][k]->EvalCapsule(C);

}

/**************************************************************/
void Branch::GetVelField() {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
        for (int i = 0; i < 216; ++i)
            if (LinISB[i]) {
                for (int n1 = 0; n1 < SYSTEM::MaxP; ++n1)
                    for (int n2 = 0; n1 + n2 < SYSTEM::MaxP; ++n2)
                        for (int n3 = 0; n1 + n2 + n3 < SYSTEM::MaxP; ++n3)
                            for (int k1 = n1; k1 < SYSTEM::MaxP; ++k1)
                                for (int k2 = n2; k1 + k2 < SYSTEM::MaxP; ++k2)
                                    for (int k3 = n3; k1 + k2 + k3 < SYSTEM::MaxP; ++k3) {
                                        VelField[n1][n2][n3] += VlFldMlt[n1][n2][n3][k1][k2][k3] * LinTlrCffts[m][indx][i][k1][k2][k3].Cross((static_cast<Branch*> (LinISB[i]))->Moments[k1 - n1][k2 - n2][k3 - n3]);
                                    }
            }
}
/**************************************************************/
void Branch::InheritVField() {
    if (m > 0) {
//        JaggedArray <Vect3> *Ptr2VField = &(static_cast<Branch*> (Parent))->VelField;
//        JaggedArray <Vect3> &ParentField = *Ptr2VField;

#ifdef USE_ROLLED_LOOPS

        for (int n1 = 0; n1 < SYSTEM::MaxP; ++n1)
            for (int n2 = 0; n2 + n1 < SYSTEM::MaxP; ++n2)
                for (int n3 = 0; n3 + n2 + n1 < SYSTEM::MaxP; ++n3)
                    for (int k1 = n1; k1 < SYSTEM::MaxP; k1++)
                        for (int k2 = n2; k1 + k2 < SYSTEM::MaxP; k2++) 
                            for (int k3 = n3; k1 + k2 + k3 < SYSTEM::MaxP; k3++) 
                                VelField[n1][n2][n3] += Node::InhrtMlt[m][x][y][z][n1][n2][n3][k1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];


#else

        if (m > 0) {
            for (int i = 0; i < Branch::InheritIndsSize; ++i)
                VelField[Branch::trgInheritInds[i][0]][Branch::trgInheritInds[i][1]][Branch::trgInheritInds[i][2]] +=
                    Branch::InheritMults[m][indx][i] * (static_cast<Branch*> (Parent))->VelField[Branch::srcInheritInds[i][0]][Branch::srcInheritInds[i][1]][Branch::srcInheritInds[i][2]];
        }
#endif
    }
}

/**************************************************************/
void Branch::vApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up) {
    (this->*down)();
    if (Children[0][0][0])
        Children[0][0][0]->vApplyRecursively(down, bottom, up);
    if (Children[1][0][0])
        Children[1][0][0]->vApplyRecursively(down, bottom, up);
    if (Children[0][1][0])
        Children[0][1][0]->vApplyRecursively(down, bottom, up);
    if (Children[0][0][1])
        Children[0][0][1]->vApplyRecursively(down, bottom, up);
    if (Children[0][1][1])
        Children[0][1][1]->vApplyRecursively(down, bottom, up);
    if (Children[1][0][1])
        Children[1][0][1]->vApplyRecursively(down, bottom, up);
    if (Children[1][1][0])
        Children[1][1][0]->vApplyRecursively(down, bottom, up);
    if (Children[1][1][1])
        Children[1][1][1]->vApplyRecursively(down, bottom, up);
    (this->*up)();
}

/**************************************************************/
void Branch::vApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up) {
    (this->*down)();
#define MODE_1
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
#pragma omp single
#endif
        {
#ifdef MODE_1
#ifdef _OPENMP
#pragma omp task
#endif
            {
                if (Children[0][0][0])
                    Children[0][0][0]->vApplyRecursivelyP(down, bottom, up);
                if (Children[1][0][0])
                    Children[1][0][0]->vApplyRecursivelyP(down, bottom, up);
                if (Children[0][1][0])
                    Children[0][1][0]->vApplyRecursivelyP(down, bottom, up);
                if (Children[0][0][1])
                    Children[0][0][1]->vApplyRecursivelyP(down, bottom, up);
            }
#ifdef _OPENMP
#pragma omp task
#endif
            {
                if (Children[0][1][1])
                    Children[0][1][1]->vApplyRecursivelyP(down, bottom, up);
                if (Children[1][0][1])
                    Children[1][0][1]->vApplyRecursivelyP(down, bottom, up);
                if (Children[1][1][0])
                    Children[1][1][0]->vApplyRecursivelyP(down, bottom, up);
                if (Children[1][1][1])
                    Children[1][1][1]->vApplyRecursivelyP(down, bottom, up);
            }
#else
            if (Children[0][0][0])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[0][0][0]->vApplyRecursivelyP(down, bottom, up);
            if (Children[1][0][0])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[1][0][0]->vApplyRecursivelyP(down, bottom, up);
            if (Children[0][1][0])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[0][1][0]->vApplyRecursivelyP(down, bottom, up);
            if (Children[0][0][1])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[0][0][1]->vApplyRecursivelyP(down, bottom, up);
            if (Children[0][1][1])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[0][1][1]->vApplyRecursivelyP(down, bottom, up);
            if (Children[1][0][1])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[1][0][1]->vApplyRecursivelyP(down, bottom, up);
            if (Children[1][1][0])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[1][1][0]->vApplyRecursivelyP(down, bottom, up);
            if (Children[1][1][1])
#ifdef _OPENMP
#pragma omp task
#endif
                Children[1][1][1]->vApplyRecursivelyP(down, bottom, up);
#endif
#ifdef _OPENMP
#pragma omp taskwait
#endif
        }
    }
    (this->*up)();
}
