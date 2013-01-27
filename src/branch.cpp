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
Branch::Branch() : Node(), Moments(globalSystem->MaxP), VelField(globalSystem->MaxP) {
    Branch::NumBranches++;
}

/**************************************************************/
Branch::Branch(Node *parent, int i, int j, int k) : Node(parent, i, j, k), Moments(globalSystem->MaxP), VelField(globalSystem->MaxP) {
    Branch::NumBranches++;
}

/**************************************************************/
void Branch::InitMomsInds(int MaxP) {

    Branch::BinomMults.allocate(OCTREE_LEVS);
    Branch::srcMomentInds.clear();
    Branch::trgMomentInds.clear();

    Branch::InheritMults.allocate(OCTREE_LEVS);
    Branch::srcInheritInds.clear();
    Branch::trgInheritInds.clear();

    for (int k1 = 0; k1 < globalSystem->MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < globalSystem->MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < globalSystem->MaxP; ++k3)
                for (int n1 = k1; n1 < globalSystem->MaxP; n1++)
                    for (int n2 = k2; n1 + n2 < globalSystem->MaxP; n2++)
                        for (int n3 = k3; n1 + n2 + n3 < globalSystem->MaxP; n3++) {
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
                    for (int k1 = 0; k1 < MaxP; ++k1)
                        for (int k2 = 0; k2 + k1 < MaxP; ++k2)
                            for (int k3 = 0; k3 + k2 + k1 < MaxP; ++k3)
                                for (int n1 = k1; n1 < globalSystem->MaxP; n1++)
                                    for (int n2 = k2; n1 + n2 < globalSystem->MaxP; n2++)
                                        for (int n3 = k3; n1 + n2 + n3 < globalSystem->MaxP; n3++) {
                                            Branch::InheritMults[mlev][Node::Indxs[i][j][k]].push_back(Node::InhrtMlt[mlev][i][j][k][k1][k2][k3][n1][n2][n3]);
                                        }
    }

    for (int k1 = 0; k1 < MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < MaxP; ++k3)
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
                    for (int k1 = 0; k1 < MaxP; ++k1)
                        for (int k2 = 0; k2 + k1 < MaxP; ++k2)
                            for (int k3 = 0; k3 + k2 + k1 < MaxP; ++k3)
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

#ifdef USE_ROLLED_LOOPS
    if ((m >= 1) && (HasLoad)) {
        for (int k1 = 0; k1 < globalSystem->MaxP; ++k1)
            for (int k2 = 0; k2 + k1 < globalSystem->MaxP; ++k2)
                for (int k3 = 0; k3 + k2 + k1 < globalSystem->MaxP; ++k3)
                    for (int n1 = 0; n1 <= k1; n1++)
                        for (int n2 = 0; n2 <= k2; n2++)
                            for (int n3 = 0; n3 <= k3; n3++)
                                (static_cast<Branch*> (Parent))->Moments[k1][k2][k3] += BinomMlt[m][x][y][z][k1][k2][k3][n1][n2][n3] * Moments[k1 - n1][k2 - n2][k3 - n3];


        Parent->HasLoad = true;
    }
#else
    if ((m >= 1) && (HasLoad)) {
        for (int i = 0; i < Branch::MomentIndsSize; ++i)
            (static_cast<Branch*> (Parent))->Moments[Branch::trgMomentInds[i][0]][Branch::trgMomentInds[i][1]][Branch::trgMomentInds[i][2]]
                += Branch::BinomMults[m][indx][i] * Moments[Branch::srcMomentInds[i][0]][Branch::srcMomentInds[i][1]][Branch::srcMomentInds[i][2]];

        Parent->HasLoad = true;
    }
#endif


}

/**************************************************************/
void Branch::vEvalCapsule(OctreeCapsule &C) {

    C.S /= 2;
    int i = (int) (C.Position.x > 0), j = (int) (C.Position.y > 0), k = (int) (C.Position.z > 0);

    if (C.has_load)
        HasLoad = true;

    Omega += C.Omega;

    C.Position.x -= C.S * (2 * i - 1);
    C.Position.y -= C.S * (2 * j - 1);
    C.Position.z -= C.S * (2 * k - 1);

    if (!Children[i][j][k]) {
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
    }

    Children[i][j][k]->EvalCapsule(C);
}

/**************************************************************/
void Branch::GetVelField() {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if (m > 0) {
        for (int ix = 0; ix < 3; ++ix)
            for (int iy = 0; iy < 3; ++iy)
                for (int iz = 0; iz < 3; ++iz)
                    if (ISA[ix][iy][iz]) ISA[ix][iy][iz]->skip_here[tid] = true; //  This is it! I think this means it is not thread safe!

        for (int ix = 0; ix < 3; ++ix)
            for (int iy = 0; iy < 3; ++iy)
                for (int iz = 0; iz < 3; ++iz)
                    if (Parent->ISA[ix][iy][iz])
                        for (int ay = 0; ay < 2; ++ay)
                            for (int be = 0; be < 2; ++be)
                                for (int ce = 0; ce < 2; ++ce)
                                    if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]) {
                                        if (!Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->skip_here[tid] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad) {
                                            int MaxP = globalSystem->MaxP;
                                            Array < Array < Array < Vect3 > > > *Ptr2Coeffts = &TlrCffts[m][x][y][z][ix][iy][iz][ay][be][ce];
                                            JaggedArray <Vect3> *Ptr2Moms = &(static_cast<Branch*> (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]))->Moments;
                                            Array < Array < Array <Vect3> > > &a = *Ptr2Coeffts;
                                            JaggedArray <Vect3> &m = *Ptr2Moms;
                                            for (int n1 = 0; n1 < MaxP; ++n1)
                                                for (int n2 = 0; n1 + n2 < MaxP; ++n2)
                                                    for (int n3 = 0; n1 + n2 + n3 < MaxP; ++n3)
                                                        for (int k1 = n1; k1 < MaxP; ++k1)
                                                            for (int k2 = n2; k1 + k2 < MaxP; ++k2)
                                                                for (int k3 = n3; k1 + k2 + k3 < MaxP; ++k3) {
                                                                    VelField[n1][n2][n3] += VlFldMlt[n1][n2][n3][k1][k2][k3] * a[k1][k2][k3].Cross(m[k1-n1][k2-n2][k3-n3]);
                                                                    //REAL mult = (pow(-1, n1) * pow(-1, n2) * pow(-1, n3)) * (REAL(globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3])) / (REAL(globalFactorial[k1 - n1] * globalFactorial[k2 - n2] * globalFactorial[k3 - n3]));
                                                                    //VelField[n1][n2][n3].x = VelField[n1][n2][n3].x + mult * (a[k1][k2][k3].y * m[k1-n1][k2-n2][k3-n3].z - a[k1][k2][k3].z * m[k1-n1][k2-n2][k3-n3].y);
                                                                    //VelField[n1][n2][n3].y = VelField[n1][n2][n3].y + mult * (a[k1][k2][k3].z * m[k1-n1][k2-n2][k3-n3].x - a[k1][k2][k3].x * m[k1-n1][k2-n2][k3-n3].z);
                                                                    //VelField[n1][n2][n3].z = VelField[n1][n2][n3].z + mult * (a[k1][k2][k3].x * m[k1-n1][k2-n2][k3-n3].y - a[k1][k2][k3].y * m[k1-n1][k2-n2][k3-n3].x);
                                                                }
                                            
                                            
                                            
                                            
                                            
//                                            
//                                            
//                                            for (int n1 = 0; n1 < globalSystem->MaxP; ++n1)
//                                                for (int n2 = 0; n2 + n1 < globalSystem->MaxP; ++n2)
//                                                    for (int n3 = 0; n3 + n2 + n1 < globalSystem->MaxP; ++n3)
//                                                        for (int k1 = n1; k1 < globalSystem->MaxP; k1++)
//                                                            for (int k2 = n2; k1 + k2 < globalSystem->MaxP; k2++)
//                                                                for (int k3 = n3; k1 + k2 + k3 < globalSystem->MaxP; k3++) {
////                                                                    REAL mult = (pow(-1,n1)*pow(-1,n2)*pow(-1,n3)) * (REAL(globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3] ))/(REAL (globalFactorial[k1 - n1] * globalFactorial[k2 - n2] * globalFactorial[k3 - n3]));
//                                                                    VelField[n1][n2][n3] += (VlFldMlt[n1][n2][n3][k1][k2][k3] *
//                                                                            TlrCffts[m][x][y][z][ix][iy][iz][ay][be][ce][k1][k2][k3].Cross((static_cast<Branch*> (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]))->Moments[k1 - n1][k2 - n2][k3 - n3]));
//                                                                }
//                                        
//                                        
//                                        
//                                        
//                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        } else {
                                            Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->skip_here[tid] = false;
                                        }
                                    }
    }
}

/**************************************************************/
void Branch::InheritVField() {
    if (m > 0) {
        JaggedArray <Vect3> *Ptr2VField = &(static_cast<Branch*> (Parent))->VelField;
        JaggedArray <Vect3> &ParentField = *Ptr2VField;

#ifdef USE_ROLLED_LOOPS

        for (int n1 = 0; n1 < globalSystem->MaxP; ++n1)
            for (int n2 = 0; n2 + n1 < globalSystem->MaxP; ++n2)
                for (int n3 = 0; n3 + n2 + n1 < globalSystem->MaxP; ++n3)
                    for (int k1 = n1; k1 < globalSystem->MaxP; k1++)
                        for (int k2 = n2; k1 + k2 < globalSystem->MaxP; k2++) 
                            for (int k3 = n3; k1 + k2 + k3 < globalSystem->MaxP; k3++) 
                                VelField[n1][n2][n3] += Node::InhrtMlt[m][x][y][z][n1][n2][n3][k1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];


#else

        if (m > 0) {
            for (int i = 0; i < Branch::InheritIndsSize; ++i)
                VelField[Branch::trgInheritInds[i][0]][Branch::trgInheritInds[i][1]][Branch::trgInheritInds[i][2]] +=
                    Branch::InheritMults[m][indx][i] * ParentField[Branch::srcInheritInds[i][0]][Branch::srcInheritInds[i][1]][Branch::srcInheritInds[i][2]];
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
