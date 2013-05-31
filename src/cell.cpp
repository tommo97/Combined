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


#include "vect34.hpp"
#include "node.hpp"
#include "types.hpp"
#include "cell.hpp"
unsigned long int FVMCell::NumCells = 0;
Array < Array <int> > FVMCell::MomentInds;
int FVMCell::MomentIndsSize;
Array < Array <REAL> > FVMCell::OffsetPows;
Array <FVMCell*> FVMCell::AllCells;
Array < Array < Vect3> > FVMCell::CellDerivs0, FVMCell::CellDerivs1, FVMCell::CellDerivs;
/**************************************************************/
FVMCell::FVMCell() : Node(), FaceVels(0.) {
    if (WRITE_TO_SCREEN) cout << "Warning - uninitialised cell" << endl;
    FVMCell::NumCells++;
}
/**************************************************************/
FVMCell::FVMCell(Node *parent, int i, int j, int k) : Node(parent, i, j, k), FaceVels(0.) {
    FVMCell::NumCells++;
    VelGrads.assign(3, Vect3(0.0, 0.0, 0.0));
#ifdef USE_MUSCL
    BEV.assign(SYSTEM::NumTransVars, Vect3());
#endif
}
/**************************************************************/
void FVMCell::InitMomsInds() {
    FVMCell::MomentInds.clear();
    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                Array <int> inds;
                inds.push_back(k1);
                inds.push_back(k2);
                inds.push_back(k3);
                FVMCell::MomentInds.push_back(inds);
            }
    FVMCell::MomentIndsSize = FVMCell::MomentInds.size();
    FVMCell::OffsetPows.allocate(8);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                Vect3 Dx = Node::Offset[i][j][k];
                for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
                    for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                        for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                            REAL mult = 0.0;
                            for (int ii = 0; ii < UTIL::QuadPts.size(); ++ii)
                                for (int jj = 0; jj < UTIL::QuadPts.size(); ++jj)
                                    for (int kk = 0; kk < UTIL::QuadPts.size(); ++kk) {
                                        REAL wt = UTIL::QuadWts[ii] * UTIL::QuadWts[jj] * UTIL::QuadWts[kk];
                                        Vect3 DXQ = Dx + Vect3(UTIL::QuadPts[ii], UTIL::QuadPts[jj], UTIL::QuadPts[kk]);
                                        mult += wt * pow(DXQ, k1, k2, k3);
                                    }
                            FVMCell::OffsetPows[Node::Indxs[i][j][k]].push_back(mult);
                        }
            }
    
    
//    This is for when there is only a single point
//    FVMCell::OffsetPows.allocate(8);
//    for (int i = 0; i < 2; ++i)
//        for (int j = 0; j < 2; ++j)
//            for (int k = 0; k < 2; ++k) {
//                Vect3 Dx = Node::Offset[i][j][k];
//                for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
//                    for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
//                        for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
//                            FVMCell::OffsetPows[Node::Indxs[i][j][k]].push_back(pow(Dx, k1, k2, k3));
//            }
//    
//    
    
};

/**************************************************************/
void FVMCell::Integrate() {
    globalTimeStepper->Integrate(this);   
    NormaliseObliterate();
}
/**************************************************************/
void FVMCell::NormaliseObliterate() {
    HasLoad = false;
    Omega = Vect3(0., 0., 0.);


    //  Normalise/obliterate the vorticities in the cell
    for (int q1 = 0; q1 < SYSTEM::NumTransVars; ++q1)
        for (int q2 = 0; q2 < SYSTEM::NumTransVars; ++q2)
            if (q1 != q2) {
                if (SIGN(TransVars[q1].x) != SIGN(TransVars[q2].x)) {
                    if (fabs(TransVars[q1].x) > fabs(TransVars[q2].x)) {
                        TransVars[q1].x += TransVars[q2].x;
                        TransVars[q2].x = 0;
                    }
                }
                if (SIGN(TransVars[q1].y) != SIGN(TransVars[q2].y)) {
                    if (fabs(TransVars[q1].y) > fabs(TransVars[q2].y)) {
                        TransVars[q1].y += TransVars[q2].y;
                        TransVars[q2].y = 0;
                    }
                }
                if (SIGN(TransVars[q1].z) != SIGN(TransVars[q2].z)) {
                    if (fabs(TransVars[q1].z) > fabs(TransVars[q2].z)) {
                        TransVars[q1].z += TransVars[q2].z;
                        TransVars[q2].z = 0;
                    }
                }

            }

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        if ((TransVars[q].x * TransVars[q].x) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
            TransVars[q].x = 0.0;
        if ((TransVars[q].y * TransVars[q].y) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
            TransVars[q].y = 0.0;
        if ((TransVars[q].z * TransVars[q].z) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
            TransVars[q].z = 0.0;
        Omega += TransVars[q];
    }
    if (Omega.Dot(Omega) > (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
        HasLoad = true;
}
/**************************************************************/
void FVMCell::CheckActive() {
    if (!HasLoad) {
        for (int i = 0; i < 6; ++i) {
            if (Neighb[i] && Neighb[i]->Omega.Dot(Neighb[i]->Omega) > (VORTICITY_CUTOFF * VORTICITY_CUTOFF)) {
                HasLoad = true;
                return;
            }
        }
    }
}

/**************************************************************/
void FVMCell::ReportToIO() {
    if (WRITE_TO_FILE) {
        for (int q = 0; q < SYSTEM::NumTransVars; ++q)
            if (TransVars[q].Mag())
                globalIO->currentOfstream << Position + .5 << " " << TransVars[q] << endl;
    }
}

/**************************************************************/
void FVMCell::ImageToIO() {
#ifdef _PNGWRITER
    if (Position.z == -0.5) {

        REAL MAXM = 2., MINM = -2, u = Omega.z;
        REAL b = .0, g = .0, r = .0, gradient = 1. / 20.;

        //  set saturation
        if (u > MAXM) u = MAXM;
        if (u < MINM) u = MINM;

        REAL transf = (u + .5 * (MAXM - MINM))*(63 / (MAXM - MINM)) + 1;
        b = g = r = 0.0;
        if (transf <= 64) {
            b = 1;
            g = (-gradient) * transf + 3.5;
            r = 0;

        }
        if (transf < 50) {
            b = (gradient) * transf - 1.5;
            g = 1;
            r = (-gradient) * transf + 2.5;
        }
        if (transf < 30) {
            b = 0;
            g = (gradient) * transf - 0.5;
            r = 1;
        }
        if (transf < 10) {
            b = 0.;
            g = 0.;
            r = (gradient) * transf + 0.5;
        }
        globalIO->CurrentPNG.plot(int (Position.x + .5), int (Position.y + .5), r, g, b);
    }
#endif
}

/**************************************************************/
void FVMCell::vEvalCapsule(OctreeCapsule &C) {
    if (C.IP)
        CollapseToIP(C);
    
    if (C.toMonitor){
        C.Ptr2CellVelocity = &Velocity;
        C.Ptr2CellPosition = &Position;
    }

    if (C.has_load) TransVars[C.AssociatedBody] += C.Omega;

}

/**************************************************************/
void FVMCell::CountCells() {
    globalOctree->CellCount++;
}

/**************************************************************/
void FVMCell::vReList() {
    FVMCell::AllCells[globalOctree->CellCount] = this;
    globalOctree->CellCount++;
}

/**************************************************************/
void FVMCell::Report() {
    Vect3 Vel, VelN, VelS, VelE, VelW, VelT, VelB, Vinf(globalSystem->scaledVinf.x, globalSystem->scaledVinf.y, globalSystem->scaledVinf.z), Vfmm = Vinf;


    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                Vfmm += CllpsMlt[x][y][z][k1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];

    for (int ix = 0; ix < 3; ++ix)
        for (int iy = 0; iy < 3; ++iy)
            for (int iz = 0; iz < 3; ++iz)
                if (Parent->ISA[ix][iy][iz])
                    for (int ay = 0; ay < 2; ++ay)
                        for (int be = 0; be < 2; ++be)
                            for (int ce = 0; ce < 2; ++ce)
                                if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad) {
                                    Vect3 PV = Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Position - Position;
                                    Vfmm += UTIL::globalDirectVel(PV, Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Omega);
                                }
    for (int i = 0; i < FVMCell::AllCells.size(); ++i) {
        Vect3 PV = FVMCell::AllCells[i]->Position - Position;
        Vel += UTIL::globalDirectVel(PV, FVMCell::AllCells[i]->Omega);
        VelN += UTIL::globalDirectVel(PV - Vect3(.0, .5, .0), FVMCell::AllCells[i]->Omega);
        VelS += UTIL::globalDirectVel(PV + Vect3(.0, .5, .0), FVMCell::AllCells[i]->Omega);
        VelE += UTIL::globalDirectVel(PV - Vect3(.5, .0, .0), FVMCell::AllCells[i]->Omega);
        VelW += UTIL::globalDirectVel(PV + Vect3(.5, .0, .0), FVMCell::AllCells[i]->Omega);
        VelT += UTIL::globalDirectVel(PV - Vect3(.0, .0, .5), FVMCell::AllCells[i]->Omega);
        VelB += UTIL::globalDirectVel(PV + Vect3(.0, .0, .5), FVMCell::AllCells[i]->Omega);
    }
    if (WRITE_TO_SCREEN) {
        cout << " x-> " << Position << " FMM: " << Vfmm << " Dir: " << Vel + Vinf << endl;
        cout << "N fmm: " << FaceVels.N << " dir: " << VelN + Vinf;
        if (Neighb.N)
            cout << " Neighb val " << static_cast<FVMCell*> (Neighb.N)->FaceVels.S;
        cout << endl;
        cout << "S fmm: " << FaceVels.S << " dir: " << VelS + Vinf;
        if (Neighb.S)
            cout << " Neighb val " << static_cast<FVMCell*> (Neighb.S)->FaceVels.N;
        cout << endl;
        cout << "E fmm: " << FaceVels.E << " dir: " << VelE + Vinf;
        if (Neighb.E)
            cout << " Neighb val " << static_cast<FVMCell*> (Neighb.E)->FaceVels.W;
        cout << endl;
        cout << "W fmm: " << FaceVels.W << " dir: " << VelW + Vinf;
        if (Neighb.W) cout << " Neighb val " << static_cast<FVMCell*> (Neighb.W)->FaceVels.E;
        cout << endl;
        cout << "T fmm: " << FaceVels.T << " dir: " << VelT + Vinf;
        if (Neighb.T) cout << " Neighb val " << static_cast<FVMCell*> (Neighb.T)->FaceVels.B;
        cout << endl;
        cout << "B fmm: " << FaceVels.B << " dir: " << VelB + Vinf;
        if (Neighb.B) cout << " Neighb val " << static_cast<FVMCell*> (Neighb.B)->FaceVels.T;
        cout << endl;
        cout << "--------- " << ID << " ProcessID " << Parent->ID << " " << Omega << endl;
    }
}

/**************************************************************/
void FVMCell::vApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up) {
    (this->*bottom)();
    down = up = NULL;
}

/**************************************************************/
void FVMCell::vApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom, BranchFuncPtr up) {
    (this->*bottom)();
    down = up = NULL;
}

/**************************************************************/
void FVMCell::SetVelsZero() {
    FaceVels = Velocity = Vect3(globalSystem->scaledVinf.x, globalSystem->scaledVinf.y, globalSystem->scaledVinf.z);
    VelGrads[0] = VelGrads[1] = VelGrads[2] = Vect3(0.0);
}
/**************************************************************/
void FVMCell::PutInOctreeCellList()
{
    OCTREE::CellsInOrder.push_back(this);
}
/**************************************************************/
FVMCell::~FVMCell() {
    FVMCell::NumCells--;
}

/**************************************************************/
void FVMCell::vCheckNeighbs() {
    if (HasLoad) {
        //  Should make some neighbours here for the FVM
        //  Quickest and easiest way is to drop some more capsules in to octree
        if (!Neighb.N) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position + Vect3(.0, 1.0, .0);
            Root->EvalCapsule(C);
        }

        if (!Neighb.S) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position - Vect3(.0, 1.0, .0);
            Root->EvalCapsule(C);
        }

        if (!Neighb.E) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position + Vect3(1.0, .0, .0);
            Root->EvalCapsule(C);
        }

        if (!Neighb.W) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position - Vect3(1.0, .0, .0);
            Root->EvalCapsule(C);
        }
#ifdef MODE_3D
        if (!Neighb.T) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position + Vect3(.0, .0, 1.0);
            Root->EvalCapsule(C);
        }

        if (!Neighb.B) {
            OctreeCapsule C;
            C.has_load = false;
            C.Position = C.tPosition = Position - Vect3(.0, .0, 1.0);
            Root->EvalCapsule(C);
        }
#endif
    }
}

/**************************************************************/

void FVMCell::SetVelsEqual() {


#ifndef COLLAPSE_TO_FACES
    FaceVels = Velocity;
    if (Neighb.S) FaceVels.S = .5 * (Velocity + static_cast<FVMCell*> (Neighb.S)->Velocity);
    if (Neighb.W) FaceVels.W = .5 * (Velocity + static_cast<FVMCell*> (Neighb.W)->Velocity);
    if (Neighb.B) FaceVels.B = .5 * (Velocity + static_cast<FVMCell*> (Neighb.B)->Velocity);
    if (Neighb.N) FaceVels.N = .5 * (Velocity + static_cast<FVMCell*> (Neighb.N)->Velocity);
    if (Neighb.E) FaceVels.E = .5 * (Velocity + static_cast<FVMCell*> (Neighb.E)->Velocity);
    if (Neighb.T) FaceVels.T = .5 * (Velocity + static_cast<FVMCell*> (Neighb.T)->Velocity);
#endif

}

/**************************************************************/
Vect3 FVMCell::ReturnSpectralRadius() {
    return Vect3(max(fabs(FaceVels.E.x), fabs(FaceVels.W.x)),
            max(fabs(FaceVels.N.y), fabs(FaceVels.S.y)),
            max(fabs(FaceVels.T.z), fabs(FaceVels.B.z)));
}

/**************************************************************/
void FVMCell::ReportSpectralRadius() {
    Vect3 srad = ReturnSpectralRadius();
    globalTimeStepper->srad.x = max(globalTimeStepper->srad.x, srad.x);
    globalTimeStepper->srad.y = max(globalTimeStepper->srad.y, srad.y);
    globalTimeStepper->srad.z = max(globalTimeStepper->srad.z, srad.z);
}

/**************************************************************/
void FVMCell::UpdateGradsFromDeltaOmega() {
    //  Updates the velocity gradients due to the changes in the nearby vorticity field

    VelGrads = VelGradsHold;
    for (int ix = 0; ix < 3; ++ix)
        for (int iy = 0; iy < 3; ++iy)
            for (int iz = 0; iz < 3; ++iz)
                if (Parent->ISA[ix][iy][iz])
                    for (int ay = 0; ay < 2; ++ay)
                        for (int be = 0; be < 2; ++be)
                            for (int ce = 0; ce < 2; ++ce)
                                if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad) {
                                    Vect3 PV = Position - Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Position;

                                    UTIL::globalCubicDirectVelGrads(PV, (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Omega - OmegaHold), VelGrads);
                                }
}

/**************************************************************/
void FVMCell::UpdateVelsFromDeltaOmega() {
    //  Updates the velocity due to the changes in the nearby vorticity field
    Velocity = VelHold;
    for (int ix = 0; ix < 3; ++ix)
        for (int iy = 0; iy < 3; ++iy)
            for (int iz = 0; iz < 3; ++iz)
                if (Parent->ISA[ix][iy][iz])
                    for (int ay = 0; ay < 2; ++ay)
                        for (int be = 0; be < 2; ++be)
                            for (int ce = 0; ce < 2; ++ce)
                                if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad) {
                                    Vect3 PV = Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Position - Position;

                                    Velocity += Node::DirVelMultsX[x][y][z][ix][iy][iz][ay][be][ce] * (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Omega - OmegaHold).x;
                                    Velocity += Node::DirVelMultsY[x][y][z][ix][iy][iz][ay][be][ce] * (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Omega - OmegaHold).y;
                                    Velocity += Node::DirVelMultsZ[x][y][z][ix][iy][iz][ay][be][ce] * (Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Omega - OmegaHold).z;
                                }
    FaceVels = Velocity;
}

/**************************************************************/

void FVMCell::vCollapseVField() {
#ifdef COLLAPSE_TO_FACES
    for (int i = 0; i < 6; ++i)
        if ((i == 0) || (i == 2) || (i == 4) || !Neighb[i])
#endif
        {
            {
                int k1 = 0;
                for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                    for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                        Velocity += Node::CllpsMlt[x][y][z][k1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];

                REAL dxkn = 1;
                for (k1 = 1; k1 < SYSTEM::MaxP; ++k1) {
                    for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2) {
                        for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                            Velocity += Node::CllpsMlt[x][y][z][k1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];
                            // du/dx, dv/dx, dw/dx NB this is += rather than -= since am reusing the CllpsMlt, which is -ve what the mult should be
                            VelGrads[0] += Node::CllpsMlt[x][y][z][k1 - 1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];
                            // du/dy, dv/dy, dw/dy
                            VelGrads[1] += Node::CllpsMlt[x][y][z][k1 - 1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1 - 1][k2 + 1][k3];
                            // du/dz, dv/dz, dw/dz
                            VelGrads[2] += Node::CllpsMlt[x][y][z][k1 - 1][k2][k3] * (static_cast<Branch*> (Parent))->VelField[k1 - 1][k2][k3 + 1];
                        }
                    }
                }
            }


            for (int i = 0; i < 216; ++i)
                if (LinISB[i]) {

                    VelGrads[0] += Node::ISBGradMultsX[indx][i][0] * LinISB[i]->Omega.x;
                    VelGrads[0] += Node::ISBGradMultsY[indx][i][0] * LinISB[i]->Omega.y;
                    VelGrads[0] += Node::ISBGradMultsZ[indx][i][0] * LinISB[i]->Omega.z;

                    VelGrads[1] += Node::ISBGradMultsX[indx][i][1] * LinISB[i]->Omega.x;
                    VelGrads[1] += Node::ISBGradMultsY[indx][i][1] * LinISB[i]->Omega.y;
                    VelGrads[1] += Node::ISBGradMultsZ[indx][i][1] * LinISB[i]->Omega.z;

                    VelGrads[2] += Node::ISBGradMultsX[indx][i][2] * LinISB[i]->Omega.x;
                    VelGrads[2] += Node::ISBGradMultsY[indx][i][2] * LinISB[i]->Omega.y;
                    VelGrads[2] += Node::ISBGradMultsZ[indx][i][2] * LinISB[i]->Omega.z;

                    Velocity += Node::ISBDirMultsX[indx][i] * LinISB[i]->Omega.x;
                    Velocity += Node::ISBDirMultsY[indx][i] * LinISB[i]->Omega.y;
                    Velocity += Node::ISBDirMultsZ[indx][i] * LinISB[i]->Omega.z;

                }

            for (int i = 0; i < 27; ++i)
                if (LinISA[i]) {
                    VelGrads[0] += Node::ISAGradMultsX[i][0] * LinISA[i]->Omega.x;
                    VelGrads[0] += Node::ISAGradMultsY[i][0] * LinISA[i]->Omega.y;
                    VelGrads[0] += Node::ISAGradMultsZ[i][0] * LinISA[i]->Omega.z;

                    VelGrads[1] += Node::ISAGradMultsX[i][1] * LinISA[i]->Omega.x;
                    VelGrads[1] += Node::ISAGradMultsY[i][1] * LinISA[i]->Omega.y;
                    VelGrads[1] += Node::ISAGradMultsZ[i][1] * LinISA[i]->Omega.z;

                    VelGrads[2] += Node::ISAGradMultsX[i][2] * LinISA[i]->Omega.x;
                    VelGrads[2] += Node::ISAGradMultsY[i][2] * LinISA[i]->Omega.y;
                    VelGrads[2] += Node::ISAGradMultsZ[i][2] * LinISA[i]->Omega.z;

                    Velocity += Node::ISADirMultsX[i] * LinISA[i]->Omega.x;
                    Velocity += Node::ISADirMultsY[i] * LinISA[i]->Omega.y;
                    Velocity += Node::ISADirMultsZ[i] * LinISA[i]->Omega.z;


                }
        }
#ifndef COLLAPSE_TO_FACES
    FaceVels = Velocity;
#endif
    VelHold = Velocity;
    VelGradsHold = VelGrads;
}

/**************************************************************/
void FVMCell::vPassMmnts2Prnt() {


    //    JaggedArray <Vect3> *Ptr2Moms = &(static_cast<Branch*> (Parent))->Moments;
    //    JaggedArray <Vect3> &ParentMoments = *Ptr2Moms;
    if (HasLoad) {
#ifdef USE_ROLLED_LOOPS 

        Vect3 Dx = Node::Offset[x][y][z];

        for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
            for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
                for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                    (static_cast<Branch*> (Parent))->Moments[k1][k2][k3] += Omega * pow(Dx, k1, k2, k3);

#else

        for (int i = 0; i < FVMCell::MomentIndsSize; ++i)
            (static_cast<Branch*> (Parent))->Moments[FVMCell::MomentInds[i][0]][FVMCell::MomentInds[i][1]][FVMCell::MomentInds[i][2]]
                += Omega * FVMCell::OffsetPows[indx][i];

#endif
        Parent->HasLoad = true;
    }

}

/**************************************************************/
void FVMCell::O1UW() {
    REAL FPE = max(FaceVels.E.x, (REAL) .0);
    REAL FPW = max(-FaceVels.W.x, (REAL) .0);
    REAL FPN = max(FaceVels.N.y, (REAL) .0);
    REAL FPS = max(-FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FPT = max(FaceVels.T.z, (REAL) .0);
    REAL FPB = max(-FaceVels.B.z, (REAL) .0);
#endif
    REAL FME = max(-FaceVels.E.x, (REAL) .0);
    REAL FMW = max(FaceVels.W.x, (REAL) .0);
    REAL FMN = max(-FaceVels.N.y, (REAL) .0);
    REAL FMS = max(FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FMT = max(-FaceVels.T.z, (REAL) .0);
    REAL FMB = max(FaceVels.B.z, (REAL) .0);
#endif
    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 fe = FPE * TransVars[q] - FME * (*Neighb_Val[q].E);
        Vect3 fw = FPW * TransVars[q] - FMW * (*Neighb_Val[q].W);
        Vect3 fn = FPN * TransVars[q] - FMN * (*Neighb_Val[q].N);
        Vect3 fs = FPS * TransVars[q] - FMS * (*Neighb_Val[q].S);
#ifdef MODE_3D
        Vect3 ft = FPT * TransVars[q] - FMT * (*Neighb_Val[q].T);
        Vect3 fb = FPB * TransVars[q] - FMB * (*Neighb_Val[q].B);
#else
        Vect3 ft, fb;
#endif

        Vect3 convection = fe + fw + fn + fs + ft + fb;
        TransDerivs[q] = -convection;
    }
}

/**************************************************************/
Vect3 FVMCell::O1UW(int q) {
    REAL FPE = max(FaceVels.E.x, (REAL) .0);
    REAL FPW = max(-FaceVels.W.x, (REAL) .0);
    REAL FPN = max(FaceVels.N.y, (REAL) .0);
    REAL FPS = max(-FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FPT = max(FaceVels.T.z, (REAL) .0);
    REAL FPB = max(-FaceVels.B.z, (REAL) .0);
#endif
    REAL FME = max(-FaceVels.E.x, (REAL) .0);
    REAL FMW = max(FaceVels.W.x, (REAL) .0);
    REAL FMN = max(-FaceVels.N.y, (REAL) .0);
    REAL FMS = max(FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FMT = max(-FaceVels.T.z, (REAL) .0);
    REAL FMB = max(FaceVels.B.z, (REAL) .0);
#endif
    Vect3 fe = FPE * TransVars[q] - FME * (*Neighb_Val[q].E);
    Vect3 fw = FPW * TransVars[q] - FMW * (*Neighb_Val[q].W);
    Vect3 fn = FPN * TransVars[q] - FMN * (*Neighb_Val[q].N);
    Vect3 fs = FPS * TransVars[q] - FMS * (*Neighb_Val[q].S);
#ifdef MODE_3D
    Vect3 ft = FPT * TransVars[q] - FMT * (*Neighb_Val[q].T);
    Vect3 fb = FPB * TransVars[q] - FMB * (*Neighb_Val[q].B);
#else
    Vect3 ft, fb;
#endif

    return -1.0 * (fe + fw + fn + fs + ft + fb);
}
/**************************************************************/
void FVMCell::GetISBVels() {
    if (globalSystem->useFMM) {
        Velocity = VelHold;
        for (int i = 0; i < 216; ++i)
            if (LinISB[i]) {
                Vect3 DeltaOmega = LinISB[i]->Omega - (static_cast<FVMCell*> (LinISB[i]))->OmegaHold;
                Velocity += Node::ISBDirMultsX[indx][i] * DeltaOmega.x;
                Velocity += Node::ISBDirMultsY[indx][i] * DeltaOmega.y;
                Velocity += Node::ISBDirMultsZ[indx][i] * DeltaOmega.z;
            }
    }
}

/**************************************************************/
void FVMCell::GetISAVels() {
    if (globalSystem->useFMM) {
        Velocity = VelHold;
        for (int i = 0; i < 27; ++i)
            if (LinISA[i]) {
                Vect3 DeltaOmega = LinISA[i]->Omega - (static_cast<FVMCell*> (LinISA[i]))->OmegaHold;
                Velocity += Node::ISADirMultsX[i] * DeltaOmega.x;
                Velocity += Node::ISADirMultsY[i] * DeltaOmega.y;
                Velocity += Node::ISADirMultsZ[i] * DeltaOmega.z;
            }
    }
}

/**************************************************************/
void FVMCell::GetISBGrads() {
    if (globalSystem->useFMM) {
        VelGrads = VelGradsHold;
        for (int i = 0; i < 216; ++i)
            if (LinISB[i]) {
                Vect3 DeltaOmega = LinISA[i]->Omega - (static_cast<FVMCell*> (LinISB[i]))->OmegaHold;
                VelGrads[0] += Node::ISBGradMultsX[indx][i][0] * DeltaOmega.x;
                VelGrads[0] += Node::ISBGradMultsY[indx][i][0] * DeltaOmega.y;
                VelGrads[0] += Node::ISBGradMultsZ[indx][i][0] * DeltaOmega.z;

                VelGrads[1] += Node::ISBGradMultsX[indx][i][1] * DeltaOmega.x;
                VelGrads[1] += Node::ISBGradMultsY[indx][i][1] * DeltaOmega.y;
                VelGrads[1] += Node::ISBGradMultsZ[indx][i][1] * DeltaOmega.z;

                VelGrads[2] += Node::ISBGradMultsX[indx][i][2] * DeltaOmega.x;
                VelGrads[2] += Node::ISBGradMultsY[indx][i][2] * DeltaOmega.y;
                VelGrads[2] += Node::ISBGradMultsZ[indx][i][2] * DeltaOmega.z;
            }
    }
}

/**************************************************************/
void FVMCell::GetISAGrads() {
    if (globalSystem->useFMM) {
        VelGrads = VelGradsHold;
        for (int i = 0; i < 27; ++i)
            if (LinISA[i]) {
                Vect3 DeltaOmega = LinISA[i]->Omega - (static_cast<FVMCell*> (LinISA[i]))->OmegaHold;
                VelGrads[0] += Node::ISAGradMultsX[i][0] * DeltaOmega.x;
                VelGrads[0] += Node::ISAGradMultsY[i][0] * DeltaOmega.y;
                VelGrads[0] += Node::ISAGradMultsZ[i][0] * DeltaOmega.z;

                VelGrads[1] += Node::ISAGradMultsX[i][1] * DeltaOmega.x;
                VelGrads[1] += Node::ISAGradMultsY[i][1] * DeltaOmega.y;
                VelGrads[1] += Node::ISAGradMultsZ[i][1] * DeltaOmega.z;

                VelGrads[2] += Node::ISAGradMultsX[i][2] * DeltaOmega.x;
                VelGrads[2] += Node::ISAGradMultsY[i][2] * DeltaOmega.y;
                VelGrads[2] += Node::ISAGradMultsZ[i][2] * DeltaOmega.z;
            }
    }
}
/**************************************************************/
void FVMCell::Stretch() {
    Vect3 GradU(VelGrads[0].x, VelGrads[1].x, VelGrads[2].x);
    Vect3 GradV(VelGrads[0].y, VelGrads[1].y, VelGrads[2].y);
    Vect3 GradW(VelGrads[0].z, VelGrads[1].z, VelGrads[2].z);
    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 StretchDeriv = Vect3(TransVars[q].Dot(GradU), TransVars[q].Dot(GradV), TransVars[q].Dot(GradW));
        TransDerivs[q] = StretchDeriv;
    }
}
/**************************************************************/
Vect3 FVMCell::Stretch(int q) {
    Vect3 GradU(VelGrads[0].x, VelGrads[1].x, VelGrads[2].x);
    Vect3 GradV(VelGrads[0].y, VelGrads[1].y, VelGrads[2].y);
    Vect3 GradW(VelGrads[0].z, VelGrads[1].z, VelGrads[2].z);
    return Vect3(TransVars[q].Dot(GradU), TransVars[q].Dot(GradV), TransVars[q].Dot(GradW));
}
/**************************************************************/
void FVMCell::Diffuse() {

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 dw = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].W)) - ((*Neighb_Val[q].E) + TransVars[q] - (*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / 12;
        Vect3 de = (4. / 3.)*((*Neighb_Val[q].E) - TransVars[q]) - ((*Neighb_Neighb_Val[q].E) + (*Neighb_Val[q].E) - TransVars[q] - (*Neighb_Val[q].W)) / 12;

        Vect3 ds = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].S)) - ((*Neighb_Val[q].N) + TransVars[q] - (*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / 12;
        Vect3 dn = (4. / 3.)*((*Neighb_Val[q].N) - TransVars[q]) - ((*Neighb_Neighb_Val[q].N) + (*Neighb_Val[q].N) - TransVars[q] - (*Neighb_Val[q].S)) / 12;

        Vect3 db = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].B)) - ((*Neighb_Val[q].T) + TransVars[q] - (*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / 12;
        Vect3 dt = (4. / 3.)*((*Neighb_Val[q].T) - TransVars[q]) - ((*Neighb_Neighb_Val[q].T) + (*Neighb_Val[q].T) - TransVars[q] - (*Neighb_Val[q].B)) / 12;
        Vect3 ViscDeriv = globalSystem->Nu * ((de - dw) + (dn - ds) + (dt - db));
        TransDerivs[q] = ViscDeriv;
    }
}
/**************************************************************/
Vect3 FVMCell::Diffuse(int q) {
        Vect3 dw = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].W)) - ((*Neighb_Val[q].E) + TransVars[q] - (*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / 12;
        Vect3 de = (4. / 3.)*((*Neighb_Val[q].E) - TransVars[q]) - ((*Neighb_Neighb_Val[q].E) + (*Neighb_Val[q].E) - TransVars[q] - (*Neighb_Val[q].W)) / 12;

        Vect3 ds = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].S)) - ((*Neighb_Val[q].N) + TransVars[q] - (*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / 12;
        Vect3 dn = (4. / 3.)*((*Neighb_Val[q].N) - TransVars[q]) - ((*Neighb_Neighb_Val[q].N) + (*Neighb_Val[q].N) - TransVars[q] - (*Neighb_Val[q].S)) / 12;

        Vect3 db = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].B)) - ((*Neighb_Val[q].T) + TransVars[q] - (*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / 12;
        Vect3 dt = (4. / 3.)*((*Neighb_Val[q].T) - TransVars[q]) - ((*Neighb_Neighb_Val[q].T) + (*Neighb_Val[q].T) - TransVars[q] - (*Neighb_Val[q].B)) / 12;
        return globalSystem->Nu * ((de - dw) + (dn - ds) + (dt - db));
}
/**************************************************************/
void FVMCell::DiffuseX() {

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 dw = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].W)) - ((*Neighb_Val[q].E) + TransVars[q] - (*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / 12;
        Vect3 de = (4. / 3.)*((*Neighb_Val[q].E) - TransVars[q]) - ((*Neighb_Neighb_Val[q].E) + (*Neighb_Val[q].E) - TransVars[q] - (*Neighb_Val[q].W)) / 12;
        Vect3 ViscDeriv = globalSystem->Nu * ((de - dw));
        TransDerivs[q] = ViscDeriv;
    }
}

/**************************************************************/
void FVMCell::DiffuseY() {

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 ds = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].S)) - ((*Neighb_Val[q].N) + TransVars[q] - (*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / 12;
        Vect3 dn = (4. / 3.)*((*Neighb_Val[q].N) - TransVars[q]) - ((*Neighb_Neighb_Val[q].N) + (*Neighb_Val[q].N) - TransVars[q] - (*Neighb_Val[q].S)) / 12;
        Vect3 ViscDeriv = globalSystem->Nu * ((dn - ds));
        TransDerivs[q] = ViscDeriv;
    }
}

/**************************************************************/
void FVMCell::DiffuseZ() {

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 db = (4. / 3.)*(TransVars[q] - (*Neighb_Val[q].B)) - ((*Neighb_Val[q].T) + TransVars[q] - (*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / 12;
        Vect3 dt = (4. / 3.)*((*Neighb_Val[q].T) - TransVars[q]) - ((*Neighb_Neighb_Val[q].T) + (*Neighb_Val[q].T) - TransVars[q] - (*Neighb_Val[q].B)) / 12;
        Vect3 ViscDeriv = globalSystem->Nu * ((dt - db));
        TransDerivs[q] = ViscDeriv;
    }
}

/**************************************************************/
void FVMCell::O2UWx() {

    REAL FPE = max(FaceVels.E.x, (REAL) .0);
    REAL FPW = max(-FaceVels.W.x, (REAL) .0);

    REAL FME = max(-FaceVels.E.x, (REAL) .0);
    REAL FMW = max(FaceVels.W.x, (REAL) .0);


    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 ep = ((*Neighb_Val[q].E) - TransVars[q]);
        Vect3 pw = (TransVars[q] - (*Neighb_Val[q].W));

        Vect3 rep = (pw / (ep + 1e-16));
        Vect3 rwp = (ep / (pw + 1e-16));

        Vect3 rem = (((*Neighb_Neighb_Val[q].E) - (*Neighb_Val[q].E)) / (ep + 1e-16));
        Vect3 rwm = (((*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / (pw + 1e-16));

        Vect3 spsip_e = flim(rep);
        Vect3 spsip_w = flim(rwp);


        Vect3 spsim_e = flim(rem);
        Vect3 spsim_w = flim(rwm);

        Vect3 fe = FPE * (TransVars[q] + .5 * spsip_e * ep) - FME * ((*Neighb_Val[q].E) - .5 * spsim_e * ep);
        Vect3 fw = FPW * (TransVars[q] - .5 * spsip_w * pw) - FMW * ((*Neighb_Val[q].W) + .5 * spsim_w * pw);


        Vect3 convection = fe + fw;
        TransDerivs[q] = -convection;
    }
}

/**************************************************************/
void FVMCell::O2UWy() {


    REAL FPN = max(FaceVels.N.y, (REAL) .0);
    REAL FPS = max(-FaceVels.S.y, (REAL) .0);

    REAL FMN = max(-FaceVels.N.y, (REAL) .0);
    REAL FMS = max(FaceVels.S.y, (REAL) .0);


    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {

        Vect3 np = ((*Neighb_Val[q].N) - TransVars[q]);
        Vect3 ps = (TransVars[q] - (*Neighb_Val[q].S));




        Vect3 rnp = (ps / (np + 1e-16));
        Vect3 rsp = (np / (ps + 1e-16));

        Vect3 rnm = (((*Neighb_Neighb_Val[q].N) - (*Neighb_Val[q].N)) / (np + 1e-16));
        Vect3 rsm = (((*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / (ps + 1e-16));
        //

        Vect3 spsip_n = flim(rnp);
        Vect3 spsip_s = flim(rsp);

        Vect3 spsim_n = flim(rnm);
        Vect3 spsim_s = flim(rsm);

        Vect3 fn = FPN * (TransVars[q] + .5 * spsip_n * np) - FMN * ((*Neighb_Val[q].N) - .5 * spsim_n * np);
        Vect3 fs = FPS * (TransVars[q] - .5 * spsip_s * ps) - FMS * ((*Neighb_Val[q].S) + .5 * spsim_s * ps);

        Vect3 convection = fn + fs;
        TransDerivs[q] = -convection;
    }
}

/**************************************************************/
void FVMCell::O2UWz() {
    REAL FPT = max(FaceVels.T.z, (REAL) .0);
    REAL FPB = max(-FaceVels.B.z, (REAL) .0);

    REAL FMT = max(-FaceVels.T.z, (REAL) .0);
    REAL FMB = max(FaceVels.B.z, (REAL) .0);


    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 tp = ((*Neighb_Val[q].T) - TransVars[q]);
        Vect3 pb = (TransVars[q] - (*Neighb_Val[q].B));



        Vect3 rtp = (pb / (tp + 1e-16));
        Vect3 rbp = (tp / (pb + 1e-16));

        Vect3 rtm = (((*Neighb_Neighb_Val[q].T) - (*Neighb_Val[q].T)) / (tp + 1e-16));
        Vect3 rbm = (((*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / (pb + 1e-16));

        Vect3 spsip_t = flim(rtp);
        Vect3 spsip_b = flim(rbp);

        Vect3 spsim_t = flim(rtm);
        Vect3 spsim_b = flim(rbm);

        Vect3 ft = FPT * (TransVars[q] + .5 * spsip_t * tp) - FMT * ((*Neighb_Val[q].T) - .5 * spsim_t * tp);
        Vect3 fb = FPB * (TransVars[q] - .5 * spsip_b * pb) - FMB * ((*Neighb_Val[q].B) + .5 * spsim_b * pb);

        Vect3 convection = ft + fb;
        TransDerivs[q] = -convection;
    }
}

/**************************************************************/
void FVMCell::O2UW() {
    REAL FPE = max(FaceVels.E.x, (REAL) .0);
    REAL FPW = max(-FaceVels.W.x, (REAL) .0);
    REAL FPN = max(FaceVels.N.y, (REAL) .0);
    REAL FPS = max(-FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FPT = max(FaceVels.T.z, (REAL) .0);
    REAL FPB = max(-FaceVels.B.z, (REAL) .0);
#endif
    REAL FME = max(-FaceVels.E.x, (REAL) .0);
    REAL FMW = max(FaceVels.W.x, (REAL) .0);
    REAL FMN = max(-FaceVels.N.y, (REAL) .0);
    REAL FMS = max(FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FMT = max(-FaceVels.T.z, (REAL) .0);
    REAL FMB = max(FaceVels.B.z, (REAL) .0);
#endif

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 ep = ((*Neighb_Val[q].E) - TransVars[q]);
        Vect3 pw = (TransVars[q] - (*Neighb_Val[q].W));
        Vect3 np = ((*Neighb_Val[q].N) - TransVars[q]);
        Vect3 ps = (TransVars[q] - (*Neighb_Val[q].S));
#ifdef MODE_3D
        Vect3 tp = ((*Neighb_Val[q].T) - TransVars[q]);
        Vect3 pb = (TransVars[q] - (*Neighb_Val[q].B));
#endif



        Vect3 rep = (pw / (ep + 1e-16));
        Vect3 rwp = (ep / (pw + 1e-16));
        Vect3 rnp = (ps / (np + 1e-16));
        Vect3 rsp = (np / (ps + 1e-16));
#ifdef MODE_3D
        Vect3 rtp = (pb / (tp + 1e-16));
        Vect3 rbp = (tp / (pb + 1e-16));
#endif
        Vect3 rem = (((*Neighb_Neighb_Val[q].E) - (*Neighb_Val[q].E)) / (ep + 1e-16));
        Vect3 rwm = (((*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / (pw + 1e-16));
        Vect3 rnm = (((*Neighb_Neighb_Val[q].N) - (*Neighb_Val[q].N)) / (np + 1e-16));
        Vect3 rsm = (((*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / (ps + 1e-16));
#ifdef MODE_3D
        Vect3 rtm = (((*Neighb_Neighb_Val[q].T) - (*Neighb_Val[q].T)) / (tp + 1e-16));
        Vect3 rbm = (((*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / (pb + 1e-16));
#endif
        //

        Vect3 spsip_e = flim(rep);
        Vect3 spsip_w = flim(rwp);
        Vect3 spsip_n = flim(rnp);
        Vect3 spsip_s = flim(rsp);
#ifdef MODE_3D
        Vect3 spsip_t = flim(rtp);
        Vect3 spsip_b = flim(rbp);
#endif
        Vect3 spsim_e = flim(rem);
        Vect3 spsim_w = flim(rwm);
        Vect3 spsim_n = flim(rnm);
        Vect3 spsim_s = flim(rsm);
#ifdef MODE_3D
        Vect3 spsim_t = flim(rtm);
        Vect3 spsim_b = flim(rbm);
#endif
        Vect3 fe = FPE * (TransVars[q] + .5 * spsip_e * ep) - FME * ((*Neighb_Val[q].E) - .5 * spsim_e * ep);
        Vect3 fw = FPW * (TransVars[q] - .5 * spsip_w * pw) - FMW * ((*Neighb_Val[q].W) + .5 * spsim_w * pw);
        Vect3 fn = FPN * (TransVars[q] + .5 * spsip_n * np) - FMN * ((*Neighb_Val[q].N) - .5 * spsim_n * np);
        Vect3 fs = FPS * (TransVars[q] - .5 * spsip_s * ps) - FMS * ((*Neighb_Val[q].S) + .5 * spsim_s * ps);
#ifdef MODE_3D
        Vect3 ft = FPT * (TransVars[q] + .5 * spsip_t * tp) - FMT * ((*Neighb_Val[q].T) - .5 * spsim_t * tp);
        Vect3 fb = FPB * (TransVars[q] - .5 * spsip_b * pb) - FMB * ((*Neighb_Val[q].B) + .5 * spsim_b * pb);
#else
        Vect3 ft, fb;
#endif

        Vect3 convection = fe + fw + fn + fs + ft + fb;
        TransDerivs[q] = -convection;
    }
}

/**************************************************************/
Vect3 FVMCell::O2UW(int q) {
    REAL FPE = max(FaceVels.E.x, (REAL) .0);
    REAL FPW = max(-FaceVels.W.x, (REAL) .0);
    REAL FPN = max(FaceVels.N.y, (REAL) .0);
    REAL FPS = max(-FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FPT = max(FaceVels.T.z, (REAL) .0);
    REAL FPB = max(-FaceVels.B.z, (REAL) .0);
#endif
    REAL FME = max(-FaceVels.E.x, (REAL) .0);
    REAL FMW = max(FaceVels.W.x, (REAL) .0);
    REAL FMN = max(-FaceVels.N.y, (REAL) .0);
    REAL FMS = max(FaceVels.S.y, (REAL) .0);
#ifdef MODE_3D
    REAL FMT = max(-FaceVels.T.z, (REAL) .0);
    REAL FMB = max(FaceVels.B.z, (REAL) .0);
#endif

    Vect3 ep = ((*Neighb_Val[q].E) - TransVars[q]);
    Vect3 pw = (TransVars[q] - (*Neighb_Val[q].W));
    Vect3 np = ((*Neighb_Val[q].N) - TransVars[q]);
    Vect3 ps = (TransVars[q] - (*Neighb_Val[q].S));
#ifdef MODE_3D
    Vect3 tp = ((*Neighb_Val[q].T) - TransVars[q]);
    Vect3 pb = (TransVars[q] - (*Neighb_Val[q].B));
#endif



    Vect3 rep = (pw / (ep + 1e-16));
    Vect3 rwp = (ep / (pw + 1e-16));
    Vect3 rnp = (ps / (np + 1e-16));
    Vect3 rsp = (np / (ps + 1e-16));
#ifdef MODE_3D
    Vect3 rtp = (pb / (tp + 1e-16));
    Vect3 rbp = (tp / (pb + 1e-16));
#endif
    Vect3 rem = (((*Neighb_Neighb_Val[q].E) - (*Neighb_Val[q].E)) / (ep + 1e-16));
    Vect3 rwm = (((*Neighb_Val[q].W) - (*Neighb_Neighb_Val[q].W)) / (pw + 1e-16));
    Vect3 rnm = (((*Neighb_Neighb_Val[q].N) - (*Neighb_Val[q].N)) / (np + 1e-16));
    Vect3 rsm = (((*Neighb_Val[q].S) - (*Neighb_Neighb_Val[q].S)) / (ps + 1e-16));
#ifdef MODE_3D
    Vect3 rtm = (((*Neighb_Neighb_Val[q].T) - (*Neighb_Val[q].T)) / (tp + 1e-16));
    Vect3 rbm = (((*Neighb_Val[q].B) - (*Neighb_Neighb_Val[q].B)) / (pb + 1e-16));
#endif
    //

    Vect3 spsip_e = flim(rep);
    Vect3 spsip_w = flim(rwp);
    Vect3 spsip_n = flim(rnp);
    Vect3 spsip_s = flim(rsp);
#ifdef MODE_3D
    Vect3 spsip_t = flim(rtp);
    Vect3 spsip_b = flim(rbp);
#endif
    Vect3 spsim_e = flim(rem);
    Vect3 spsim_w = flim(rwm);
    Vect3 spsim_n = flim(rnm);
    Vect3 spsim_s = flim(rsm);
#ifdef MODE_3D
    Vect3 spsim_t = flim(rtm);
    Vect3 spsim_b = flim(rbm);
#endif
    Vect3 fe = FPE * (TransVars[q] + .5 * spsip_e * ep) - FME * ((*Neighb_Val[q].E) - .5 * spsim_e * ep);
    Vect3 fw = FPW * (TransVars[q] - .5 * spsip_w * pw) - FMW * ((*Neighb_Val[q].W) + .5 * spsim_w * pw);
    Vect3 fn = FPN * (TransVars[q] + .5 * spsip_n * np) - FMN * ((*Neighb_Val[q].N) - .5 * spsim_n * np);
    Vect3 fs = FPS * (TransVars[q] - .5 * spsip_s * ps) - FMS * ((*Neighb_Val[q].S) + .5 * spsim_s * ps);
#ifdef MODE_3D
    Vect3 ft = FPT * (TransVars[q] + .5 * spsip_t * tp) - FMT * ((*Neighb_Val[q].T) - .5 * spsim_t * tp);
    Vect3 fb = FPB * (TransVars[q] - .5 * spsip_b * pb) - FMB * ((*Neighb_Val[q].B) + .5 * spsim_b * pb);
#else
    Vect3 ft, fb;
#endif

    return -1.0 * (fe + fw + fn + fs + ft + fb);
}
/**************************************************************/
inline static Vect3 minmod(Vect3 x, Vect3 y);

inline static Vect3 minmod(Vect3 x, Vect3 y) {

    return 0.5 * (sign(x) + sign(y)) * min(fabs(x), fabs(y));
}

inline static Vect3 minmod3(Vect3 x, Vect3 y, Vect3 z);

inline static Vect3 minmod3(Vect3 x, Vect3 y, Vect3 z) {

    return minmod(x, minmod(y, z));
}
#ifdef USE_MUSCL
/**************************************************************/
void FVMCell::GetBEV() {

    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        /*Get Boundary Extrapolated Values...*/
        Vect3 SLOPEx = minmod3(BETA * ((*Neighb_Val[q].E) - TransVars[q]), .5 * ((*Neighb_Val[q].E) - (*Neighb_Val[q].W)), BETA * (TransVars[q] - (*Neighb_Val[q].W)));
        Vect3 SLOPEy = minmod3(BETA * ((*Neighb_Val[q].N) - TransVars[q]), .5 * ((*Neighb_Val[q].N) - (*Neighb_Val[q].S)), BETA * (TransVars[q] - (*Neighb_Val[q].S)));
        Vect3 SLOPEz = minmod3(BETA * ((*Neighb_Val[q].T) - TransVars[q]), .5 * ((*Neighb_Val[q].T) - (*Neighb_Val[q].B)), BETA * (TransVars[q] - (*Neighb_Val[q].B)));


        BEV[q].E = TransVars[q] + 0.5 * SLOPEx;
        BEV[q].W = TransVars[q] - 0.5 * SLOPEx;
        BEV[q].N = TransVars[q] + 0.5 * SLOPEy;
        BEV[q].S = TransVars[q] - 0.5 * SLOPEy;
        BEV[q].T = TransVars[q] + 0.5 * SLOPEz;
        BEV[q].B = TransVars[q] - 0.5 * SLOPEz;
    }

}

/**************************************************************/
void FVMCell::MUSCL() {
    /*    Approximation of the spectral radii with face velocities */
    Vect3 srad = ReturnSpectralRadius();
    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        Vect3 Ew, We, Ns, Sn, Tb, Bt; //  BEVs at current cell's interfaces
        Vect3 Pn = BEV[q].N;
        Vect3 Ps = BEV[q].S;
        Vect3 Pe = BEV[q].E;
        Vect3 Pw = BEV[q].W;
        Vect3 Pt = BEV[q].T;
        Vect3 Pb = BEV[q].B;
        //
        //  BEVs on other side of current cell interfaces
        if (Neighb.E)
            Ew = Neighb.E->BEV[q].W;

        if (Neighb.W)
            We = Neighb.W->BEV[q].E;

        if (Neighb.N)
            Ns = Neighb.N->BEV[q].S;

        if (Neighb.S)
            Sn = Neighb.S->BEV[q].N;

        if (Neighb.T)
            Tb = Neighb.T->BEV[q].B;

        if (Neighb.B)
            Bt = Neighb.B->BEV[q].T;



        Vect3 fPe = Pe * FaceVels.E.x;
        Vect3 fPw = Pw * FaceVels.W.x;

        Vect3 gPn = Pn * FaceVels.N.y;
        Vect3 gPs = Ps * FaceVels.S.y;

        Vect3 hPt = Pt * FaceVels.T.z;
        Vect3 hPb = Pb * FaceVels.B.z;


        Vect3 fEw = Ew * FaceVels.E.x;
        Vect3 fWe = We * FaceVels.W.x;

        Vect3 gNs = Ns * FaceVels.N.y;
        Vect3 gSn = Sn * FaceVels.S.y;

        Vect3 hTb = Tb * FaceVels.T.z;
        Vect3 hBt = Bt * FaceVels.B.z;


        Vect3 Hx_e = 0.5 * ((fEw + fPe) - srad.x * (Ew - Pe));
        Vect3 Hx_w = 0.5 * ((fPw + fWe) - srad.x * (Pw - We));

        Vect3 Hy_n = 0.5 * ((gNs + gPn) - srad.y * (Ns - Pn));
        Vect3 Hy_s = 0.5 * ((gPs + gSn) - srad.y * (Ps - Sn));

        Vect3 Hz_t = 0.5 * ((hTb + hPt) - srad.z * (Tb - Pt));
        Vect3 Hz_b = 0.5 * ((hPb + hBt) - srad.z * (Pb - Bt));


        Vect3 convection = (Hx_e - Hx_w) + (Hy_n - Hy_s) + (Hz_t - Hz_b);
        TransDerivs[q] = -convection;
    }
}
#endif
/**************************************************************/
void FVMCell::CollapseToIP(OctreeCapsule &IP) {

    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
                IP.Velocity -= (pow(IP.tPosition - Parent->Position, k1, k2, k3) /
                    (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3])) *
                (static_cast<Branch*> (Parent))->VelField[k1][k2][k3];

    for (int ix = 0; ix < 3; ++ix)
        for (int iy = 0; iy < 3; ++iy)
            for (int iz = 0; iz < 3; ++iz)
                if (Parent->ISA[ix][iy][iz])
                    for (int ay = 0; ay < 2; ++ay)
                        for (int be = 0; be < 2; ++be)
                            for (int ce = 0; ce < 2; ++ce)
                                if (Parent->ISA[ix][iy][iz]->Children[ay][be][ce] && Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->HasLoad)
                                    IP.Velocity += UTIL::globalCubicDirectVel(IP.tPosition - Parent->ISA[ix][iy][iz]->Children[ay][be][ce]->Position,
                                        static_cast<FVMCell*> (Parent->ISA[ix][iy][iz]->Children[ay][be][ce])->Omega);
}


/**************************************************************/
void FVMCell::AdvanceDt(REAL dt) {

    // Advance transport variables...
    for (int q = 0; q < TransVars.size(); ++q) {
        TransVars[q] += dt * TransDerivs[q];
    }




    //    //  Normalise/obliterate the vorticities in the cell
    //    for (int q1 = 0; q1 < SYSTEM::NumTransVars; ++q1)
    //        for (int q2 = 0; q2 < SYSTEM::NumTransVars; ++q2)
    //            if (q1 != q2) {
    //                if (SIGN(TransVars[q1].x) != SIGN(TransVars[q2].x)) {
    //                    if (fabs(TransVars[q1].x) > fabs(TransVars[q2].x)) {
    //                        TransVars[q1].x += TransVars[q2].x;
    //                        TransVars[q2].x = 0;
    //                    }
    //                }
    //                if (SIGN(TransVars[q1].y) != SIGN(TransVars[q2].y)) {
    //                    if (fabs(TransVars[q1].y) > fabs(TransVars[q2].y)) {
    //                        TransVars[q1].y += TransVars[q2].y;
    //                        TransVars[q2].y = 0;
    //                    }
    //                }
    //                if (SIGN(TransVars[q1].z) != SIGN(TransVars[q2].z)) {
    //                    if (fabs(TransVars[q1].z) > fabs(TransVars[q2].z)) {
    //                        TransVars[q1].z += TransVars[q2].z;
    //                        TransVars[q2].z = 0;
    //                    }
    //                }
    //            }

    Omega = Vect3(0., 0., 0.);
    for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
        //        if ((TransVars[q].x * TransVars[q].x) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
        //            TransVars[q].x = 0.0;
        //        if ((TransVars[q].y * TransVars[q].y) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
        //            TransVars[q].y = 0.0;
        //        if ((TransVars[q].z * TransVars[q].z) < (VORTICITY_CUTOFF * VORTICITY_CUTOFF))
        //            TransVars[q].z = 0.0;
        Omega += TransVars[q];
    }


}