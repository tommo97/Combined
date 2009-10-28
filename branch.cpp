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


#include "branch.hpp"

/**************************************************************/
Branch::~Branch() {
    globalNum_BRANCHES--;
}

/**************************************************************/
Branch::Branch() : Node(), Moments(globalSystem->MaxP), VelField(globalSystem->MaxP) {
    globalNum_BRANCHES++;
}

/**************************************************************/
Branch::Branch(Node *parent, int i, int j, int k) : Node(parent, i, j, k), Moments(globalSystem->MaxP), VelField(globalSystem->MaxP) {
    globalNum_BRANCHES++;
}
/**************************************************************/
void Branch::BranchCount() {
        globalOctree->BranchCount[m]++;
}
/**************************************************************/
void Branch::ReList() {
        globalOctree->AllBranches[m][globalOctree->BranchCount[m]] = this;
        globalOctree->BranchCount[m]++;
}
/**************************************************************/
void Branch::SetFieldsZero() {
    VelField = Vect3(0.);
    Moments = Vect3(0.);
}

/**************************************************************/
void Branch::vPassMmnts2Prnt() {
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
        if (C.S > 0)
        {
            Branch *child = new Branch(this, i, j, k);
            Children[i][j][k] = child;
            if (C.IP)
            {
                child->InheritVField();
                child->GetVelField();
            }
        }

        else
        {
            Children[i][j][k] = new FVMCell(this, i, j, k);
        }
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
                    if (ISA[ix][iy][iz]) ISA[ix][iy][iz]->skip_here[tid] = true;     //  This is it! I think this means it is not thread safe!

        for (int ix = 0; ix < 3; ++ix)
            for (int iy = 0; iy < 3; ++iy)
                for (int iz = 0; iz < 3; ++iz)
                    if (Parent->ISA[ix][iy][iz])
                        for (int ay = 0; ay < 2; ++ay)
                            for (int be = 0; be < 2; ++be)
                                for (int ce = 0; ce < 2; ++ce)
                                    if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]) {
                                        if (!Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->skip_here[tid] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad) {
                                            for (int k1 = 0; k1 < globalSystem->MaxP; ++k1)
                                                for (int k2 = 0; k2 + k1 < globalSystem->MaxP; ++k2)
                                                    for (int k3 = 0; k3 + k2 + k1 < globalSystem->MaxP; ++k3)
                                                        for (int n1 = k1; n1 < globalSystem->MaxP; n1++)
                                                            for (int n2 = k2; n1 + n2 < globalSystem->MaxP; n2++)
                                                                for (int n3 = k3; n1 + n2 + n3 < globalSystem->MaxP; n3++){
                                                                    VelField[k1][k2][k3] += VlFldMlt[k1][k2][k3][n1][n2][n3] *
                                                                        TlrCffts[m][x][y][z][ix][iy][iz][ay][be][ce][n1][n2][n3].Cross((static_cast<Branch*> (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]))->Moments[n1 - k1][n2 - k2][n3 - k3]);
                                                                }
                                        } else {
                                            Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->skip_here[tid] = false;
                                        }
                                    }
    }
}

/**************************************************************/
void Branch::InheritVField() {
    if (m > 0) {
        //        Vect3 Dx = Parent->Position - Position;
        for (int k1 = 0; k1 < globalSystem->MaxP; ++k1)
            for (int k2 = 0; k2 + k1 < globalSystem->MaxP; ++k2)
                for (int k3 = 0; k3 + k2 + k1 < globalSystem->MaxP; ++k3)
                    for (int n1 = k1; n1 < globalSystem->MaxP; n1++)
                        for (int n2 = k2; n1 + n2 < globalSystem->MaxP; n2++)
                            for (int n3 = k3; n1 + n2 + n3 < globalSystem->MaxP; n3++)
                                VelField[k1][k2][k3] += Node::InhrtMlt[m][x][y][z][k1][k2][k3][n1][n2][n3] * (static_cast<Branch*> (Parent))->VelField[n1][n2][n3];
                            
    }
}
//#define MODE_1

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
