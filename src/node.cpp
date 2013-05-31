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



#include "node.hpp"
#include "cell.hpp"
#include "vect34.hpp"
unsigned long int Node::NumNodes = 0;
unsigned long int Node::NodeCount = 0; 
int Node::MomentSize = 0;
Array <Array <int> > Node::CellMomentIndex(3, Array <int> ());
Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > Node::DirVelMultsX;
Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > Node::DirVelMultsY;
Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > Node::DirVelMultsZ;
Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > Node::DirGradMultsX;
Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > Node::DirGradMultsY;
Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > Node::DirGradMultsZ;
Array <Node*> Node::AllNodes, Node::UpList, Node::DownList;
Array < Array <int> > Node::ISBRecipInds(8, Array <int> (216,-1));
Array <int> Node::ISARecipInds(27,-1);
Array < Array < Array <Vect3> > > Node::ISBGradMultsX(8, Array < Array <Vect3 > > (216, Array <Vect3> (3, Vect3(0.0))));
Array < Array < Array <Vect3> > > Node::ISBGradMultsY(8, Array < Array <Vect3 > > (216, Array <Vect3> (3, Vect3(0.0))));
Array < Array < Array <Vect3> > > Node::ISBGradMultsZ(8, Array < Array <Vect3 > > (216, Array <Vect3> (3, Vect3(0.0))));
Array < Array <Vect3> > Node::ISBDirMultsX(8, Array <Vect3 > (216, Vect3(0.)));
Array < Array <Vect3> > Node::ISBDirMultsY(8, Array <Vect3 > (216, Vect3(0.)));
Array < Array <Vect3> > Node::ISBDirMultsZ(8, Array <Vect3 > (216, Vect3(0.)));

Array < Array <Vect3> > Node::ISAGradMultsX(27, Array <Vect3 > (3, Vect3(0.)));
Array < Array <Vect3> > Node::ISAGradMultsY(27, Array <Vect3 > (3, Vect3(0.)));
Array < Array <Vect3> > Node::ISAGradMultsZ(27, Array <Vect3 > (3, Vect3(0.)));
Array <Vect3> Node::ISADirMultsX(Array < Vect3 > (27, Vect3(0.)));
Array <Vect3> Node::ISADirMultsY(Array < Vect3 > (27, Vect3(0.)));
Array <Vect3> Node::ISADirMultsZ(Array < Vect3 > (27, Vect3(0.)));
    
    

/**************************************************************/
Node::~Node() {
    RemoveFromNeighbs();
    ClearChildren();
    if (Parent)
        Parent->Children[x][y][z] = NULL;
    Node::NumNodes--;
    Trans.clear();
}

/**************************************************************/
Node::Node() : Parent(NULL), x(-1), y(-1), z(-1), m(0), indx(-1), ID(00), size(OCTREE_SIZE),
SkipAdd2List(false), Neighb(NULL), HasLoad(false), toMonitor(false) {
    Neighb_Val.assign(SYSTEM::NumTransVars, &ZERO);
    Neighb_Neighb_Val.assign(SYSTEM::NumTransVars, &ZERO);
    TransVars.assign(SYSTEM::NumTransVars, Vect3());
    TransVarsHold.assign(SYSTEM::NumTransVars, Vect3());
    TransVars0.assign(SYSTEM::NumTransVars, Vect3());
    TransDerivs.assign(SYSTEM::NumTransVars, Vect3());
    SetISANull();
    SetKidsNull();
    Node::NumNodes++;
    Trans.push_back(-1);
}

/**************************************************************/
Node::Node(Node *parent, int i, int j, int k) : Parent(parent), x(i), y(j), z(k), m(parent->m + 1),indx(Indxs[i][j][k]),
ID(((((((parent->ID << 1) + i) << 1) + j) << 1) + k)), size(.5 * parent->size), SkipAdd2List(false),
Position(parent->Position + .5 * Node::Offset[i][j][k] * parent->size), Neighb(NULL), HasLoad(false), toMonitor(false) {
    Neighb_Val.assign(SYSTEM::NumTransVars, &ZERO);
    Neighb_Neighb_Val.assign(SYSTEM::NumTransVars, &ZERO);
    TransVars.assign(SYSTEM::NumTransVars, Vect3());
    TransVarsHold.assign(SYSTEM::NumTransVars, Vect3());
    TransVars0.assign(SYSTEM::NumTransVars, Vect3());
    TransDerivs.assign(SYSTEM::NumTransVars, Vect3());
    Trans = parent->Trans;
    Trans.push_back(REF[i][j][k]);
    SetISANull();
    SetKidsNull();
//    GetISA();
    Node::NumNodes++;
    skip_here.assign(globalSystem->NumThreads, false);


}
/**************************************************************/
const Vect3 Node::Offset[][2][2] = {
    {
        {Vect3(-.5, -.5, -.5), Vect3(-.5, -.5, .5)},
        {Vect3(-.5, .5, -.5), Vect3(-.5, .5, .5)}
    },
    {
        {Vect3(.5, -.5, -.5), Vect3(.5, -.5, .5)},
        {Vect3(.5, .5, -.5), Vect3(.5, .5, .5)}
    }
};

const int Node::Indxs[][2][2] = {
    {
        {0, 1},
        {2, 3}
    },
    {
        {4, 5},
        {6, 7}
    }
};

const Vect3 Node::NeighbOffset[6] = {Vect3(0.0, 0.5, 0.0), Vect3(0.0, -.5, 0.0),
    Vect3(0.5, 0.0, 0.0), Vect3(-.5, 0.0, 0.0),
    Vect3(0.0, 0.0, 0.5), Vect3(0.0, 0.0, -.5)};

Vect3 Node::ZERO = Vect3(0.);
Node* Node::Root = NULL;
int Node::RootSize = OCTREE_SIZE;

int Node::REF[][2][2] = {
    {
        {0, 4},
        {1, 5}
    },
    {
        {2, 6},
        {3, 7}
    }
};
/**************************************************************/
ARRAY13(Vect3) Node::TlrCffts(OCTREE_LEVS);
ARRAY6(REAL) Node::VlFldMlt; //  This is modified (once) later
ARRAY10(Vect3) Node::TlrCfftsdx(OCTREE_LEVS);
ARRAY6(NeighbSet <REAL>) Node::FaceCllpsMlt(2);
ARRAY10(REAL) Node::BinomMlt(OCTREE_LEVS);
ARRAY10(REAL) Node::InhrtMlt(OCTREE_LEVS);
ARRAY6(REAL) Node::CllpsMlt(2);
ARRAY6 (Vect3) Node::LinTlrCffts(OCTREE_LEVS);
const int Node::trans_neighb[][5] = {
    {0, 2, 4, 6, 1},
    {1, 3, 5, 7, -1},
    {0, 1, 4, 5, 2},
    {2, 3, 6, 7, -2},
    {0, 1, 2, 3, 4},
    {4, 5, 6, 7, -4}
};
const int Node::deREF[][3] = {
    {0, 0, 0},
    {0, 1, 0},
    {1, 0, 0},
    {1, 1, 0},
    {0, 0, 1},
    {0, 1, 1},
    {1, 0, 1},
    {1, 1, 1}
};
int Node::Op[6] = {1, 0, 3, 2, 5, 4};

/**************************************************************/
void Node::RemoveFromNeighbs() {
    for (int i = 0; i < 6; ++i) {
        if (Neighb[i]) {
            for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
                Neighb[i]->Neighb_Val[q][Op[i]] = &ZERO;
                if (Neighb[i]->Neighb[i]) {
                    Neighb[i]->Neighb[i]->Neighb_Neighb_Val[q][Op[i]] = &ZERO;
                }
            }
            Neighb[i]->Neighb[Op[i]] = NULL;
        }
    }
    for (int i = 0, I = 2; i < 3; ++i, --I)
        for (int j = 0, J = 2; j < 3; ++j, --J)
            for (int k = 0, K = 2; k < 3; ++k, --K)
                if (ISA[i][j][k])
                    ISA[i][j][k]->ISA[I][J][K] = NULL;

    for (int i = 0; i < 216; ++i)
        if (LinISB[i]) {
            LinISB[i]->LinISB[Node::ISBRecipInds[indx][i]] = NULL;
        }

    for (int i = 0; i < 27; ++i)
        if (LinISA[i])
            LinISA[i]->LinISA[ISARecipInds[i]] = NULL;
}

/**************************************************************/
void Node::Recurrance(JaggedArray <REAL> &coeffts, int k1, int k2, int k3, Vect3 diff, REAL R2) {
    //  Checked
    REAL absK = REAL(k1 + k2 + k3), m1 = (2 * absK - 1) / (R2 * absK), m2 = (absK - 1) / (R2 * absK);
    if (k1 > 0) {
        coeffts[k1][k2][k3] += m1 * diff.x * coeffts[k1 - 1][k2][k3];
        if (k1 > 1)
            coeffts[k1][k2][k3] -= m2 * coeffts[k1 - 2][k2][k3];
    }
    if (k2 > 0) {
        coeffts[k1][k2][k3] += m1 * diff.y * coeffts[k1][k2 - 1][k3];
        if (k2 > 1)
            coeffts[k1][k2][k3] -= m2 * coeffts[k1][k2 - 2][k3];
    }
    if (k3 > 0) {
        coeffts[k1][k2][k3] += m1 * diff.z * coeffts[k1][k2][k3 - 1];
        if (k3 > 1)
            coeffts[k1][k2][k3] -= m2 * coeffts[k1][k2][k3 - 2];
    }
}

/**************************************************************/
void Node::SetUpISBIndices() {
    //  This is a bit of a cludge...



    {
        int count = 0;
        for (int i = -1; i < 2; ++i)
            for (int j = -1; j < 2; ++j)
                for (int k = -1; k < 2; ++k) {
                    Vect3 R = Vect3(1.0 * REAL(i), 1.0 * REAL(j), 1.0 * REAL(k)); //   PV from source  centroid to target ISA centroid   

                    {
                        Vect3 Vx = UTIL::globalCubicDirectVel(-1.0*R, Vect3(1., 0., 0.));
                        Vect3 Vy = UTIL::globalCubicDirectVel(-1.0*R, Vect3(0., 1., 0.));
                        Vect3 Vz = UTIL::globalCubicDirectVel(-1.0*R, Vect3(0., 0., 1.));
                        Node::ISADirMultsX[count] = Vx;
                        Node::ISADirMultsY[count] = Vy;
                        Node::ISADirMultsZ[count] = Vz;


                        Array <Vect3 > GradsX(3, Vect3(0.0));
                        UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(1., 0., 0.), GradsX);
                        Array <Vect3 > GradsY(3, Vect3(0.0));
                        UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(0., 1., 0.), GradsY);
                        Array <Vect3 > GradsZ(3, Vect3(0.0));
                        UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(0., 0., 1.), GradsZ);

                        Node::ISAGradMultsX[count] = GradsX;
                        Node::ISAGradMultsY[count] = GradsY;
                        Node::ISAGradMultsZ[count] = GradsZ;
                    }



                    int count2 = 0;
                    for (int I = -1; I < 2; ++I)
                        for (int J = -1; J < 2; ++J)
                            for (int K = -1; K < 2; ++K) {
                                Vect3 R2 = Vect3(1.0 * REAL(I), 1.0 * REAL(J), 1.0 * REAL(K)); //   PV from source parent centroid to target ISA centroid 
                                if ((R + R2).Mag() < 1e-9) // i.e. if one PV is the equal and opposite of the other
                                    Node::ISARecipInds[count] = count2;

                                count2++;
                            }

                    count++;
                }
    }


    //  For all 8 nodes in parent, find all ISBs, then for each ISB, go through all of its ISBs until the original node is found - this gives the reciprocal
    Node::ISBRecipInds = Array < Array <int > > (8, Array <int> (216));
    for (int ix = 0; ix < 2; ++ix)
        for (int iy = 0; iy < 2; ++iy)
            for (int iz = 0; iz < 2; ++iz) { //   This has looped over the 8 nodes inside
                
                Vect3 R1 = -1.0 * Offset[ix][iy][iz]; //     PV to target from parent centroid
                int count = 0;
                
                for (int i = -1; i < 2; ++i)
                    for (int j = -1; j < 2; ++j)
                        for (int k = -1; k < 2; ++k) { //     This has looped over 27 of the parent ISA
                            
                            Vect3 R2 = Vect3(2.0 * REAL(i), 2.0 * REAL(j), 2.0 * REAL(k)); //   PV from source parent centroid to target parent centroid     
                            for (int jx = 0; jx < 2; ++jx)
                                for (int jy = 0; jy < 2; ++jy)
                                    for (int jz = 0; jz < 2; ++jz) { //   This has looped over 8 nodes of each parent ISA
                                        
                                        Vect3 R3 = 1.0 * Offset[jx][jy][jz]; //      PV from source cell to source parent.
                                        { //  Need to exclude ISA of initial node -- these are cells which are |R1 + R2 + R3| < (1, sqrt(2) or sqrt(3)) away
                                            {
                                                Vect3 R = R1 + R2 + R3;
                                                Vect3 Vx = UTIL::globalCubicDirectVel(-1.0 * R, Vect3(1., 0., 0.));
                                                Vect3 Vy = UTIL::globalCubicDirectVel(-1.0 * R, Vect3(0., 1., 0.));
                                                Vect3 Vz = UTIL::globalCubicDirectVel(-1.0 * R, Vect3(0., 0., 1.));
                                                Node::ISBDirMultsX[Indxs[ix][iy][iz]][count] = Vx;
                                                Node::ISBDirMultsY[Indxs[ix][iy][iz]][count] = Vy;
                                                Node::ISBDirMultsZ[Indxs[ix][iy][iz]][count] = Vz;


                                                Array <Vect3 > GradsX(3, Vect3(0.0));
                                                UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(1., 0., 0.), GradsX);
                                                Array <Vect3 > GradsY(3, Vect3(0.0));
                                                UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(0., 1., 0.), GradsY);
                                                Array <Vect3 > GradsZ(3, Vect3(0.0));
                                                UTIL::globalCubicDirectVelGrads(-1.0*R, Vect3(0., 0., 1.), GradsZ);

                                                Node::ISBGradMultsX[Indxs[ix][iy][iz]][count] = GradsX;
                                                Node::ISBGradMultsY[Indxs[ix][iy][iz]][count] = GradsY;
                                                Node::ISBGradMultsZ[Indxs[ix][iy][iz]][count] = GradsZ;
                                                
                                                
                                                
                                                
                                                
                                                
                                                
                                                Vect3 RR1 = -1.0*Offset[jx][jy][jz]; //     PV to target from parent centroid
                                                //      Now have position of ISB node - now work tho' its ISB until reaching original node
                                                int count2 = 0;
                                                for (int I = -1; I < 2; ++I)
                                                    for (int J = -1; J < 2; ++J)
                                                        for (int K = -1; K < 2; ++K) {
                                                            Vect3 RR2 = Vect3(2.0 * REAL(I), 2.0 * REAL(J), 2.0 * REAL(K)); //   PV from source parent centroid to target parent centroid     
                                                            for (int JX = 0; JX < 2; ++JX)
                                                                for (int JY = 0; JY < 2; ++JY)
                                                                    for (int JZ = 0; JZ < 2; ++JZ) {
                                                                        Vect3 RR3 = 1.0 * Offset[JX][JY][JZ]; //      PV from source cell to source parent.
                                                                            if ((R + RR1 + RR2 + RR3).Mag() < 1e-9)
                                                                                Node::ISBRecipInds[Indxs[ix][iy][iz]][count] = count2;
                                                                            count2++;
                                                                    }
                                                        }
                                                count++;

                                            }
                                        }
                                    }
                        }
            }
}
/**************************************************************/
void Node::CompCoeffts(Vect3 diff, JaggedArray <REAL> &coeffts) {
//#define VERSION_2
#ifdef VERSION_2
    REAL R2 = diff.Dot(diff) +  globalSystem->Del2;

    // base case
#ifdef MODE_3D
    coeffts[0][0][0] = sqrt(1 / (R2)) / four_pi; // divide by 4*Pi to get Phi(z) ORIGINAL
#else
    coeffts[0][0][0] = sqrt(1 / (R2)) / two_pi; // divide by 4*Pi to get Phi(z) ORIGINAL
#endif

    /* two of the indices are zero */
    Recurrance(coeffts, 0, 0, 1, diff, R2);
    Recurrance(coeffts, 0, 1, 0, diff, R2);
    Recurrance(coeffts, 1, 0, 0, diff, R2);
    for (int i = 2; i <= SYSTEM::MaxP; ++i) {
        Recurrance(coeffts, 0, 0, i, diff, R2);
        Recurrance(coeffts, 0, i, 0, diff, R2);
        Recurrance(coeffts, i, 0, 0, diff, R2);
    }

    /* one index = 0, one index = 1, other index >= 1 */
    Recurrance(coeffts, 0, 1, 1, diff, R2);
    Recurrance(coeffts, 1, 0, 1, diff, R2);
    Recurrance(coeffts, 1, 1, 0, diff, R2);
    for (int i = 2; i < SYSTEM::MaxP; ++i) {
        Recurrance(coeffts, 0, 1, i, diff, R2);
        Recurrance(coeffts, 1, 0, i, diff, R2);
        Recurrance(coeffts, 0, i, 1, diff, R2);
        Recurrance(coeffts, 1, i, 0, diff, R2);
        Recurrance(coeffts, i, 1, 0, diff, R2);
        Recurrance(coeffts, i, 0, 1, diff, R2);
    }

    /* one index = 0, other indices >= 2 */
    for (int i = 2; i <= SYSTEM::MaxP - 2; ++i)
        for (int j = 2; i + j <= SYSTEM::MaxP; ++j) {
            Recurrance(coeffts, 0, i, j, diff, R2);
            Recurrance(coeffts, i, 0, j, diff, R2);
            Recurrance(coeffts, i, j, 0, diff, R2);
        }

    /* two indices = 1, other index >= 1 */
    Recurrance(coeffts, 1, 1, 1, diff, R2);
    for (int i = 2; i <= SYSTEM::MaxP - 2; ++i) {
        Recurrance(coeffts, 1, 1, i, diff, R2);
        Recurrance(coeffts, 1, i, 1, diff, R2);
        Recurrance(coeffts, i, 1, 1, diff, R2);
    }

    /* one index = 1, other indices >= 2 */
    for (int i = 2; i <= SYSTEM::MaxP - 3; ++i)
        for (int j = 2; i + j < SYSTEM::MaxP; ++j) {
            Recurrance(coeffts, 1, i, j, diff, R2);
            Recurrance(coeffts, i, 1, j, diff, R2);
            Recurrance(coeffts, i, j, 1, diff, R2);
        }


    /* all indices >= 2 */
    for (int i = 2; i <= SYSTEM::MaxP - 4; ++i)
        for (int j = 2; i + j <= SYSTEM::MaxP - 2; ++j)
            for (int k = 2; i + j + k <= SYSTEM::MaxP; ++k)
                Recurrance(coeffts, i, j, k, diff, R2);



#else
    // use this instead
			Vect3 Diff = -1.0*diff;
    REAL del2 = globalSystem->Del2;
    int MAX_P = SYSTEM::MaxP;
    int p = MAX_P;
    Array <REAL> cff1(p + 1), cff2(p + 1);
    REAL x12, x22, x32, mult;
   int i, j, k;

   for (i = 1; i <= p; i++) {
      cff1[i] = 1 - (REAL) 1 / (2 * i);
      cff2[i] = 1 - (REAL) 1 / i;
   }
   x12 = Diff.x * 2;
   x22 = Diff.y * 2;
   x32 = Diff.z * 2;
   mult = -1 / (del2 + Diff.x * Diff.x + Diff.y * Diff.y + Diff.z * Diff.z);

   /* base case */
   coeffts[0][0][0] = sqrt (-mult)/four_pi;

   /* two of the indices are zero */
   coeffts[0][0][1] = Diff.z * (coeffts[0][0][0] * mult);
   coeffts[0][1][0] = Diff.y * (coeffts[0][0][0] * mult);
   coeffts[1][0][0] = Diff.x * (coeffts[0][0][0] * mult);
   for (i = 2; i <= p; i++) {
      coeffts[0][0][i] = (x32 * cff1[i] * coeffts[0][0][i - 1] + cff2[i] * coeffts[0][0][i - 2]) * mult;
      coeffts[0][i][0] = (x22 * cff1[i] * coeffts[0][i - 1][0] + cff2[i] * coeffts[0][i - 2][0]) * mult;
      coeffts[i][0][0] = (x12 * cff1[i] * coeffts[i - 1][0][0] + cff2[i] * coeffts[i - 2][0][0]) * mult;
   }

   /* one index = 0, one index = 1, other index >= 1 */
   coeffts[0][1][1] = 3 * Diff.z * (coeffts[0][1][0] * mult);
   coeffts[1][0][1] = 3 * Diff.z * (coeffts[1][0][0] * mult);
   coeffts[1][1][0] = 3 * Diff.y * (coeffts[1][0][0] * mult);
   for (i = 2; i < p; i++) {
      coeffts[0][1][i] = (Diff.y * coeffts[0][0][i] + x32 * coeffts[0][1][i - 1] + coeffts[0][1][i - 2]) * mult;
      coeffts[1][0][i] = (Diff.x * coeffts[0][0][i] + x32 * coeffts[1][0][i - 1] + coeffts[1][0][i - 2]) * mult;
      coeffts[0][i][1] = (x22 * coeffts[0][i - 1][1] + coeffts[0][i - 2][1] + Diff.z * coeffts[0][i][0]) * mult;
      coeffts[1][i][0] = (Diff.x * coeffts[0][i][0] + x22 * coeffts[1][i - 1][0] + coeffts[1][i - 2][0]) * mult;
      coeffts[i][1][0] = (x12 * coeffts[i - 1][1][0] + coeffts[i - 2][1][0] + Diff.y * coeffts[i][0][0]) * mult;
      coeffts[i][0][1] = (x12 * coeffts[i - 1][0][1] + coeffts[i - 2][0][1] + Diff.z * coeffts[i][0][0]) * mult;
   }

   /* one index = 0, other indices >= 2 */
   for (i = 2; i <= p - 2; i++)
      for (j = 2; i + j <= p; j++) {
	 coeffts[0][i][j] =
	    (x22 * cff1[i] * coeffts[0][i - 1][j] + cff2[i] * coeffts[0][i - 2][j] + x32 * coeffts[0][i][j - 1] + coeffts[0][i][j - 2]) * mult;
	 coeffts[i][0][j] =
	    (x12 * cff1[i] * coeffts[i - 1][0][j] + cff2[i] * coeffts[i - 2][0][j] + x32 * coeffts[i][0][j - 1] + coeffts[i][0][j - 2]) * mult;
	 coeffts[i][j][0] =
	    (x12 * cff1[i] * coeffts[i - 1][j][0] + cff2[i] * coeffts[i - 2][j][0] + x22 * coeffts[i][j - 1][0] + coeffts[i][j - 2][0]) * mult;
      }

   /* two indices = 1, other index >= 1 */
   coeffts[1][1][1] = 5 * Diff.z * mult * coeffts[1][1][0];
   for (i = 2; i <= p - 2; i++) {
      coeffts[1][1][i] = (Diff.x * coeffts[0][1][i] + x22 * coeffts[1][0][i] + x32 * coeffts[1][1][i - 1] + coeffts[1][1][i - 2]) * mult;
      coeffts[1][i][1] = (Diff.x * coeffts[0][i][1] + x22 * coeffts[1][i - 1][1] + coeffts[1][i - 2][1] + x32 * coeffts[1][i][0]) * mult;
      coeffts[i][1][1] = (x12 * coeffts[i - 1][1][1] + coeffts[i - 2][1][1] + Diff.y * coeffts[i][0][1] + x32 * coeffts[i][1][0]) * mult;
   }

   /* one index = 1, other indices >= 2 */
   for (i = 2; i <= p - 3; i++)
      for (j = 2; i + j < p; j++) {
	 coeffts[1][i][j] =
	    (Diff.x * coeffts[0][i][j] + x22 * coeffts[1][i - 1][j] + coeffts[1][i - 2][j] + x32 * coeffts[1][i][j - 1] + coeffts[1][i][j - 2]) * mult;
	 coeffts[i][1][j] =
	    (x12 * coeffts[i - 1][1][j] + coeffts[i - 2][1][j] + Diff.y * coeffts[i][0][j] + x32 * coeffts[i][1][j - 1] + coeffts[i][1][j - 2]) * mult;
	 coeffts[i][j][1] =
	    (x12 * coeffts[i - 1][j][1] + coeffts[i - 2][j][1] + x22 * coeffts[i][j - 1][1] + coeffts[i][j - 2][1] + Diff.z * coeffts[i][j][0]) * mult;
      }

   /* all indices >= 2 */
   for (i = 2; i <= p - 4; i++)
      for (j = 2; i + j <= p - 2; j++)
	 for (k = 2; i + j + k <= p; k++)
	    coeffts[i][j][k] =
	       (x12 * cff1[i] * coeffts[i - 1][j][k] + cff2[i] * coeffts[i - 2][j][k] + x22 * coeffts[i][j - 1][k] + coeffts[i][j - 2][k] +
		x32 * coeffts[i][j][k - 1] + coeffts[i][j][k - 2]) * mult;
 
#endif
}

/**************************************************************/
void Node::GetISA() {
    
    
    if (m > 0)
    {
        int count = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k){
                    if (!LinISA[count])
                        LinISA[count] = ReturnNeighb(i, j, k);
                    if (LinISA[count])
                        LinISA[count]->LinISA[Node::ISARecipInds[count]] = this;
                    count ++;
                }
        
        
        
    }
    
    
    
    
    
    
    
    ISA[1][1][1] = this;
    if (m > 0) {
        for (int i = 0, I = 2; i < 3; ++i, --I)
            for (int j = 0, J = 2; j < 3; ++j, --J)
                for (int k = 0, K = 2; k < 3; ++k, --K) {
                    if (!ISA[i][j][k]) {
                        ISA[i][j][k] = ReturnNeighb(i, j, k);
                        if (ISA[i][j][k])
                            ISA[i][j][k]->ISA[I][J][K] = this;
                    }
                }
        Neighb.E = ISA[2][1][1];
        Neighb.W = ISA[0][1][1];
        Neighb.N = ISA[1][2][1];
        Neighb.S = ISA[1][0][1];
        Neighb.T = ISA[1][1][2];
        Neighb.B = ISA[1][1][0];
      
//    Neighb.E = ReturnNeighb(2, 1, 1);
//    Neighb.W = ReturnNeighb(0, 1, 1);
//    Neighb.N = ReturnNeighb(1, 2, 1);
//    Neighb.S = ReturnNeighb(1, 0, 1);
//    Neighb.T = ReturnNeighb(1, 1, 2);
//    Neighb.B = ReturnNeighb(1, 1, 0);
        }
    for (int i = 0; i < 6; ++i) {
        if (Neighb[i]) {
            Neighb[i]->Neighb[Op[i]] = this;
            for (int q = 0; q < SYSTEM::NumTransVars; ++q) {
                Neighb_Val[q][i] = &(Neighb[i]->TransVars[q]);
                Neighb[i]->Neighb_Val[q][Op[i]] = &TransVars[q];
                if (Neighb[i]->Neighb[i]) {
                    Neighb_Neighb_Val[q][i] = &(Neighb[i]->Neighb[i]->TransVars[q]);
                    Neighb[i]->Neighb[i]->Neighb_Neighb_Val[q][Op[i]] = &TransVars[q];
                }
            }
        }
    }
    
}

/**************************************************************/
void Node::GetISB() {
    if (m > 1) {
        int count = 0;
        for (int i = 0; i < 27; ++i)
            if (LinISA[i])
                LinISA[i]->SkipAdd2List = true;

        for (int i = 0; i < 27; ++i)
            for (int ay = 0; ay < 2; ++ay)
                for (int be = 0; be < 2; ++be)
                    for (int ce = 0; ce < 2; ++ce) {
                        if (Parent->LinISA[i] && Parent->LinISA[i]->Children[ay][be][ce] && (!Parent->LinISA[i]->Children[ay][be][ce]->SkipAdd2List)) {
                            LinISB[count] = Parent->LinISA[i]->Children[ay][be][ce];
                            LinISB[count]->LinISB[Node::ISBRecipInds[indx][count]] = this;
                        }
                        count++;
                    }
        for (int i = 0; i < 27; ++i)
            if (LinISA[i])
                LinISA[i]->SkipAdd2List = false;
    }
}
/**************************************************************/
Node* Node::ReturnNeighb(int i, int j, int k) {
    Array <int> v;
    v = this->Trans;
    if (i == 0) Trans2Neighb(v, 3);
    else if (i == 2) Trans2Neighb(v, 2);
    if (j == 0) Trans2Neighb(v, 1);
    else if (j == 2) Trans2Neighb(v, 0);
    if (k == 0) Trans2Neighb(v, 5);
    else if (k == 2) Trans2Neighb(v, 4);
    if (v[0] != -1) return NULL;
#ifndef USE_ARRAY
    Array <int> temp;
    for (int i = 1; i < v.size(); ++i) temp.push_back(v[i]);
    v = temp;
#else
    v.pop_front();
#endif

    return Root->MoveDownTree(v);
}

/**************************************************************/
void Node::Trans2Neighb(Array <int> &v, int dirn) {
    int M = (int) v.size() - 1;
    bool isdone = false;
    while ((M >= 0) && (!isdone)) {
        if ((v[M] == trans_neighb[dirn][0]) || (v[M] == trans_neighb[dirn][1]) || (v[M] == trans_neighb[dirn][2])
                || (v[M] == trans_neighb[dirn][3])) {
            v[M] += trans_neighb[dirn][4];
            isdone = true;
        } else {
            v[M] -= trans_neighb[dirn][4];
            M--;
        }
    }
}

/**************************************************************/
Node* Node::MoveDownTree(Array <int> &Directions) {
    int ind = Directions[0];
    if (ind == -1) return NULL;
    if (Directions.size() > 1) {
#ifndef USE_ARRAY
        Array <int> temp;
        for (int i = 1; i < Directions.size(); ++i) temp.push_back(Directions[i]);
        Directions = temp;
#else
        Directions.pop_front();
#endif

        if (Children[deREF[ind][0]][deREF[ind][1]][deREF[ind][2]])
            return Children[deREF[ind][0]][deREF[ind][1]][deREF[ind][2]]->MoveDownTree(Directions);
        else
            return NULL;
    } else
        return Children[deREF[ind][0]][deREF[ind][1]][deREF[ind][2]];
}

/**************************************************************/
void Node::EvalCapsule(OctreeCapsule &c) {
    if (c.has_load) {
        HasLoad = true;
        Omega += c.Omega;
    }
    vEvalCapsule(c);
}

/**************************************************************/
void Node::UpdateMomentMults() {
    
    for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1)
        for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2)
            for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3)
            {
                Node::MomentSize++;
                Node::CellMomentIndex[0].push_back(k1);
                Node::CellMomentIndex[1].push_back(k2);
                Node::CellMomentIndex[2].push_back(k3);

            }

    cout << "MomentSize " << Node::MomentSize << endl;
    
    DirVelMultsX = ARRAY9(Vect3) (2);
    DirVelMultsY = ARRAY9(Vect3) (2);
    DirVelMultsZ = ARRAY9(Vect3) (2);
    DirGradMultsX = ARRAY10(Vect3) (2);
    DirGradMultsY = ARRAY10(Vect3) (2);
    DirGradMultsZ = ARRAY10(Vect3) (2);

    for (int ix = 0; ix < 2; ++ix) {
        DirVelMultsX[ix] = ARRAY8(Vect3) (2);
        DirVelMultsY[ix] = ARRAY8(Vect3) (2);
        DirVelMultsZ[ix] = ARRAY8(Vect3) (2);
        DirGradMultsX[ix] = ARRAY9(Vect3) (2);
        DirGradMultsY[ix] = ARRAY9(Vect3) (2);
        DirGradMultsZ[ix] = ARRAY9(Vect3) (2);
        for (int iy = 0; iy < 2; ++iy) {
            DirVelMultsX[ix][iy] = ARRAY7(Vect3) (2);
            DirVelMultsY[ix][iy] = ARRAY7(Vect3) (2);
            DirVelMultsZ[ix][iy] = ARRAY7(Vect3) (2);
            DirGradMultsX[ix][iy] = ARRAY8(Vect3) (2);
            DirGradMultsY[ix][iy] = ARRAY8(Vect3) (2);
            DirGradMultsZ[ix][iy] = ARRAY8(Vect3) (2);            
            for (int iz = 0; iz < 2; ++iz) {
                DirVelMultsX[ix][iy][iz] = ARRAY6(Vect3) (3);
                DirVelMultsY[ix][iy][iz] = ARRAY6(Vect3) (3);
                DirVelMultsZ[ix][iy][iz] = ARRAY6(Vect3) (3);
                DirGradMultsX[ix][iy][iz] = ARRAY7(Vect3) (3);
                DirGradMultsY[ix][iy][iz] = ARRAY7(Vect3) (3);
                DirGradMultsZ[ix][iy][iz] = ARRAY7(Vect3) (3);
                Vect3 R1 = -1.0*Offset[ix][iy][iz]; //     PV to target from parent centroid
                for (int i = -1, I = 0; i < 2; ++i, ++I) {
                    DirVelMultsX[ix][iy][iz][I] = ARRAY5(Vect3) (3);
                    DirVelMultsY[ix][iy][iz][I] = ARRAY5(Vect3) (3);
                    DirVelMultsZ[ix][iy][iz][I] = ARRAY5(Vect3) (3);
                    DirGradMultsX[ix][iy][iz][I] = ARRAY6(Vect3) (3);
                    DirGradMultsY[ix][iy][iz][I] = ARRAY6(Vect3) (3);
                    DirGradMultsZ[ix][iy][iz][I] = ARRAY6(Vect3) (3);
                    for (int j = -1, J = 0; j < 2; ++j, ++J) {
                        DirVelMultsX[ix][iy][iz][I][J] = ARRAY4(Vect3) (3);
                        DirVelMultsY[ix][iy][iz][I][J] = ARRAY4(Vect3) (3);
                        DirVelMultsZ[ix][iy][iz][I][J] = ARRAY4(Vect3) (3);
                        DirGradMultsX[ix][iy][iz][I][J] = ARRAY5(Vect3) (3);
                        DirGradMultsY[ix][iy][iz][I][J] = ARRAY5(Vect3) (3);
                        DirGradMultsZ[ix][iy][iz][I][J] = ARRAY5(Vect3) (3);
                        for (int k = -1, K = 0; k < 2; ++k, ++K) {
                            DirVelMultsX[ix][iy][iz][I][J][K] = ARRAY3(Vect3) (2);
                            DirVelMultsY[ix][iy][iz][I][J][K] = ARRAY3(Vect3) (2);
                            DirVelMultsZ[ix][iy][iz][I][J][K] = ARRAY3(Vect3) (2);
                            DirGradMultsX[ix][iy][iz][I][J][K] = ARRAY4(Vect3) (2);
                            DirGradMultsY[ix][iy][iz][I][J][K] = ARRAY4(Vect3) (2);
                            DirGradMultsZ[ix][iy][iz][I][J][K] = ARRAY4(Vect3) (2);
                            Vect3 R2 = Vect3(2.0 * REAL(i), 2.0 * REAL(j), 2.0 * REAL(k)); //   PV from source parent centroid to target parent centroid     
                            for (int jx = 0; jx < 2; ++jx) {
                                DirVelMultsX[ix][iy][iz][I][J][K][jx] = ARRAY2(Vect3) (2);
                                DirVelMultsY[ix][iy][iz][I][J][K][jx] = ARRAY2(Vect3) (2);
                                DirVelMultsZ[ix][iy][iz][I][J][K][jx] = ARRAY2(Vect3) (2);
                                DirGradMultsX[ix][iy][iz][I][J][K][jx] = ARRAY3(Vect3) (2);
                                DirGradMultsY[ix][iy][iz][I][J][K][jx] = ARRAY3(Vect3) (2);
                                DirGradMultsZ[ix][iy][iz][I][J][K][jx] = ARRAY3(Vect3) (2);

                                for (int jy = 0; jy < 2; ++jy) {
                                    DirVelMultsX[ix][iy][iz][I][J][K][jx][jy] = Array <Vect3> (2);
                                    DirVelMultsY[ix][iy][iz][I][J][K][jx][jy] = Array <Vect3> (2);
                                    DirVelMultsZ[ix][iy][iz][I][J][K][jx][jy] = Array <Vect3> (2);
                                    DirGradMultsX[ix][iy][iz][I][J][K][jx][jy] = Array < Array <Vect3> > (2);
                                    DirGradMultsY[ix][iy][iz][I][J][K][jx][jy] = Array < Array <Vect3> > (2);
                                    DirGradMultsZ[ix][iy][iz][I][J][K][jx][jy] = Array < Array <Vect3> > (2);
                                    for (int jz = 0; jz < 2; ++jz) {
                                        Vect3 R3 = 1.0 * Offset[jx][jy][jz]; //      PV from source cell to source parent.
                                        Vect3 Vx = UTIL::globalCubicDirectVel(-1.0*(R1 + R2 + R3), Vect3(1.,0.,0.));
                                        Vect3 Vy = UTIL::globalCubicDirectVel(-1.0*(R1 + R2 + R3), Vect3(0.,1.,0.));
                                        Vect3 Vz = UTIL::globalCubicDirectVel(-1.0*(R1 + R2 + R3), Vect3(0.,0.,1.));
                                        
                                        DirVelMultsX[ix][iy][iz][I][J][K][jx][jy][jz] = Vx;
                                        DirVelMultsY[ix][iy][iz][I][J][K][jx][jy][jz] = Vy;
                                        DirVelMultsZ[ix][iy][iz][I][J][K][jx][jy][jz] = Vz;
                                                
                                        DirGradMultsX[ix][iy][iz][I][J][K][jx][jy][jz] = Array <Vect3 > (3);
                                        DirGradMultsY[ix][iy][iz][I][J][K][jx][jy][jz] = Array <Vect3 > (3);
                                        DirGradMultsZ[ix][iy][iz][I][J][K][jx][jy][jz] = Array <Vect3 > (3);
                                        Array <Vect3 > GradsX(3, Vect3(0.0));
                                        UTIL::globalCubicDirectVelGrads(-1.0*(R1 + R2 + R3), Vect3(1., 0., 0.), GradsX);
                                        Array <Vect3 > GradsY(3, Vect3(0.0));
                                        UTIL::globalCubicDirectVelGrads(-1.0*(R1 + R2 + R3), Vect3(0., 1., 0.), GradsY);
                                        Array <Vect3 > GradsZ(3, Vect3(0.0));
                                        UTIL::globalCubicDirectVelGrads(-1.0*(R1 + R2 + R3), Vect3(0., 0., 1.), GradsZ);
                                        DirGradMultsX[ix][iy][iz][I][J][K][jx][jy][jz] = GradsX;
                                        DirGradMultsY[ix][iy][iz][I][J][K][jx][jy][jz] = GradsY;
                                        DirGradMultsZ[ix][iy][iz][I][J][K][jx][jy][jz] = GradsZ;
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }
                    
                
            
        
    
    REAL Size = Node::RootSize;
    Node::VlFldMlt = ARRAY6(REAL) (SYSTEM::MaxP);
    //	Go down OCTREE and do a test run on each level

    for (int mlev = 1; mlev < OCTREE_LEVS; ++mlev) {
        TlrCffts[mlev] = ARRAY12(Vect3) (2);
        LinTlrCffts[mlev] = ARRAY5(Vect3) (8);
        TlrCfftsdx[mlev] = ARRAY9(Vect3) (2);
        BinomMlt[mlev] = InhrtMlt[mlev] = ARRAY9(REAL) (2);

        for (int ix = 0; ix < 2; ++ix) {
            TlrCffts[mlev][ix] = ARRAY11(Vect3) (2);
            TlrCfftsdx[mlev][ix] = ARRAY8(Vect3) (2);
            BinomMlt[mlev][ix] = InhrtMlt[mlev][ix] = ARRAY8(REAL) (2);

            for (int iy = 0; iy < 2; ++iy) {
                TlrCffts[mlev][ix][iy] = ARRAY10(Vect3) (2);
                TlrCfftsdx[mlev][ix][iy] = ARRAY7(Vect3) (2);
                BinomMlt[mlev][ix][iy] = InhrtMlt[mlev][ix][iy] = ARRAY7(REAL) (2);

                for (int iz = 0; iz < 2; ++iz) {
                    BinomMlt[mlev][ix][iy][iz] = InhrtMlt[mlev][ix][iy][iz] = ARRAY6(REAL) (SYSTEM::MaxP);
                    Vect3 Dx = .5 * Offset[ix][iy][iz] * Size;

                    for (int n1 = 0; n1 < SYSTEM::MaxP; ++n1) {
                        BinomMlt[mlev][ix][iy][iz][n1] = InhrtMlt[mlev][ix][iy][iz][n1] = ARRAY5(REAL) (SYSTEM::MaxP);

                        for (int n2 = 0; n2 + n1 < SYSTEM::MaxP; ++n2) {
                            BinomMlt[mlev][ix][iy][iz][n1][n2] = InhrtMlt[mlev][ix][iy][iz][n1][n2] = ARRAY4(REAL) (SYSTEM::MaxP);

                            for (int n3 = 0; n3 + n2 + n1 < SYSTEM::MaxP; ++n3) {
                                BinomMlt[mlev][ix][iy][iz][n1][n2][n3] = ARRAY3(REAL) (n1 + 1);

                                for (int k1 = 0; k1 <= n1; ++k1) {
                                    BinomMlt[mlev][ix][iy][iz][n1][n2][n3][k1] = ARRAY2(REAL) (n2 + 1);

                                    for (int k2 = 0; k2 <= n2; ++k2) {
                                        BinomMlt[mlev][ix][iy][iz][n1][n2][n3][k1][k2] = Array <REAL > (n3 + 1, 0.);

                                        for (int k3 = 0; k3 <= n3; ++k3) {
                                            BinomMlt[mlev][ix][iy][iz][n1][n2][n3][k1][k2][k3] = (REAL) (globalFactorial[n1] * globalFactorial[n2] * globalFactorial[n3]) * (REAL) pow(Dx, k1, k2, k3) /
                                                    (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3] * globalFactorial[n1 - k1] * globalFactorial[n2 - k2] * globalFactorial[n3 - k3]);
                                        }
                                    }
                                }
                                InhrtMlt[mlev][ix][iy][iz][n1][n2][n3] = ARRAY3(REAL) (SYSTEM::MaxP);

                                for (int k1 = n1; k1 < SYSTEM::MaxP; ++k1) {
                                    InhrtMlt[mlev][ix][iy][iz][n1][n2][n3][k1] = ARRAY2(REAL) (SYSTEM::MaxP);

                                    for (int k2 = n2; k2 + k1 < SYSTEM::MaxP; ++k2) {
                                        InhrtMlt[mlev][ix][iy][iz][n1][n2][n3][k1][k2] = Array <REAL > (SYSTEM::MaxP, 0.);

                                        for (int k3 = n3; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                                            unsigned long int k_nfactorial = globalFactorial[k1 - n1] * globalFactorial[k2 - n2] * globalFactorial[k3 - n3];
                                            //  Velocity Field Inheritance Multipliers
                                            InhrtMlt[mlev][ix][iy][iz][n1][n2][n3][k1][k2][k3] = pow(Dx, (k1 - n1), (k2 - n2), (k3 - n3)) / REAL(k_nfactorial);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //  Now Prep Taylor coefficients etc.
                    //  First do a dummy run at Parents Level
                    TlrCfftsdx[mlev][ix][iy][iz] = ARRAY6(Vect3) (3);
                    TlrCffts[mlev][ix][iy][iz] = ARRAY9(Vect3) (3);
                    LinTlrCffts[mlev][Indxs[ix][iy][iz]] = ARRAY4(Vect3) (216);
                    int count = 0;
                    for (int i = -1, I = 0; i < 2; ++i, ++I) {
                        TlrCfftsdx[mlev][ix][iy][iz][I] = ARRAY5(Vect3) (3);
                        TlrCffts[mlev][ix][iy][iz][I] = ARRAY8(Vect3) (3);

                        for (int j = -1, J = 0; j < 2; ++j, ++J) {
                            TlrCfftsdx[mlev][ix][iy][iz][I][J] = ARRAY4(Vect3) (3);
                            TlrCffts[mlev][ix][iy][iz][I][J] = ARRAY7(Vect3) (3);

                            for (int k = -1, K = 0; k < 2; ++k, ++K) {
                                TlrCfftsdx[mlev][ix][iy][iz][I][J][K] = ARRAY3(Vect3) (2);
                                TlrCffts[mlev][ix][iy][iz][I][J][K] = ARRAY6(Vect3) (2);
                                

                                for (int jx = 0; jx < 2; ++jx) {
                                    TlrCfftsdx[mlev][ix][iy][iz][I][J][K][jx] = ARRAY2(Vect3) (2);
                                    TlrCffts[mlev][ix][iy][iz][I][J][K][jx] = ARRAY5(Vect3) (2);

                                    for (int jy = 0; jy < 2; ++jy) {
                                        TlrCfftsdx[mlev][ix][iy][iz][I][J][K][jx][jy] = Array <Vect3 > (2);
                                        TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy] = ARRAY4(Vect3) (2);

                                        for (int jz = 0; jz < 2; ++jz) {
                                            Vect3 diff = Dx - Size * Vect3(REAL(i), REAL(j), REAL(k)) - 0.5 * Offset[jx][jy][jz] * Size; // check signs here
                                            diff = -1.0*diff;
                                            JaggedArray <REAL> coeffts(SYSTEM::MaxP + 1);
                                            coeffts = 0.0;
                                            CompCoeffts(diff, coeffts);
                                            TlrCfftsdx[mlev][ix][iy][iz][I][J][K][jx][jy][jz] = diff;
                                            TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz] = ARRAY3(Vect3) (SYSTEM::MaxP);
                                            LinTlrCffts[mlev][Indxs[ix][iy][iz]][count] = ARRAY3(Vect3) (SYSTEM::MaxP);
                                            for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1) {
                                                TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz][k1] = ARRAY2(Vect3) (SYSTEM::MaxP);
                                                LinTlrCffts[mlev][Indxs[ix][iy][iz]][count][k1] = ARRAY2(Vect3) (SYSTEM::MaxP);
                                                for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2) {
                                                    TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz][k1][k2] = Array <Vect3 > (SYSTEM::MaxP);
                                                    LinTlrCffts[mlev][Indxs[ix][iy][iz]][count][k1][k2] = Array <Vect3 > (SYSTEM::MaxP);
                                                    for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                                                        //TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz][k1][k2][k3].x = pow(-1,k1+k2+k3) * (k1 + 1) * coeffts[k1 + 1][k2][k3];
                                                        //TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz][k1][k2][k3].y = pow(-1,k1+k2+k3) * (k2 + 1) * coeffts[k1][k2 + 1][k3];
                                                        //TlrCffts[mlev][ix][iy][iz][I][J][K][jx][jy][jz][k1][k2][k3].z = pow(-1,k1+k2+k3) * (k3 + 1) * coeffts[k1][k2][k3 + 1];
                                                        LinTlrCffts[mlev][Indxs[ix][iy][iz]][count][k1][k2][k3].x  = pow(-1.,k1+k2+k3) * (k1 + 1) * coeffts[k1 + 1][k2][k3];
                                                        LinTlrCffts[mlev][Indxs[ix][iy][iz]][count][k1][k2][k3].y  = pow(-1.,k1+k2+k3) * (k2 + 1) * coeffts[k1][k2 + 1][k3];
                                                        LinTlrCffts[mlev][Indxs[ix][iy][iz]][count][k1][k2][k3].z  = pow(-1.,k1+k2+k3) * (k3 + 1) * coeffts[k1][k2][k3 + 1];
                                                    }
                                                }
                                            }
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Size /= 2;
    }



    for (int ix = 0; ix < 2; ++ix) {
        CllpsMlt[ix] = ARRAY5(REAL) (2);
        FaceCllpsMlt[ix] = ARRAY5(NeighbSet <REAL>) (2);

        for (int iy = 0; iy < 2; ++iy) {
            CllpsMlt[ix][iy] = ARRAY4(REAL) (2);
            FaceCllpsMlt[ix][iy] = ARRAY4(NeighbSet <REAL>) (2);

            for (int iz = 0; iz < 2; ++iz) {
                CllpsMlt[ix][iy][iz] = ARRAY3(REAL) (SYSTEM::MaxP);
                FaceCllpsMlt[ix][iy][iz] = ARRAY3(NeighbSet <REAL>) (SYSTEM::MaxP);

                for (int k1 = 0; k1 < SYSTEM::MaxP; ++k1) {
                    CllpsMlt[ix][iy][iz][k1] = ARRAY2(REAL) (SYSTEM::MaxP);
                    FaceCllpsMlt[ix][iy][iz][k1] = ARRAY2(NeighbSet <REAL>) (SYSTEM::MaxP);

                    for (int k2 = 0; k2 + k1 < SYSTEM::MaxP; ++k2) {
                        CllpsMlt[ix][iy][iz][k1][k2] = Array <REAL > (SYSTEM::MaxP);
                        FaceCllpsMlt[ix][iy][iz][k1][k2] = Array < NeighbSet <REAL> > (SYSTEM::MaxP);

                        for (int k3 = 0; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                            CllpsMlt[ix][iy][iz][k1][k2][k3] = pow(Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].N = pow(Vect3(0., 0.5, .0) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].S = pow(Vect3(0., -.5, .0) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].E = pow(Vect3(0.5, .0, .0) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].W = pow(Vect3(-.5, .0, .0) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].T = pow(Vect3(0., .0, 0.5) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                            FaceCllpsMlt[ix][iy][iz][k1][k2][k3].B = pow(Vect3(0., .0, -.5) + Offset[ix][iy][iz], k1, k2, k3) / (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]);
                        }
                    }
                }
            }
        }
    }


    for (int n1 = 0; n1 < SYSTEM::MaxP; ++n1) {
        VlFldMlt[n1] = ARRAY5(REAL) (SYSTEM::MaxP);
        for (int n2 = 0; n2 + n1 < SYSTEM::MaxP; ++n2) {
            VlFldMlt[n1][n2] = ARRAY4(REAL) (SYSTEM::MaxP);
            for (int n3 = 0; n3 + n2 + n1 < SYSTEM::MaxP; ++n3) {
                REAL mult = pow(-1.0, n1) * pow(-1.0, n2) * pow(-1.0, n3);
                VlFldMlt[n1][n2][n3] = ARRAY3(REAL) (SYSTEM::MaxP);
                for (int k1 = n1; k1 < SYSTEM::MaxP; ++k1) {
                    VlFldMlt[n1][n2][n3][k1] = ARRAY2(REAL) (SYSTEM::MaxP);
                    for (int k2 = n2; k2 + k1 < SYSTEM::MaxP; ++k2) {
                        VlFldMlt[n1][n2][n3][k1][k2] = Array <REAL > (SYSTEM::MaxP);
                        for (int k3 = n3; k3 + k2 + k1 < SYSTEM::MaxP; ++k3) {
                            REAL k_nfactorial = (REAL) (globalFactorial[k1 - n1] * globalFactorial[k2 - n2] * globalFactorial[k3 - n3]);
                            //  Velocity Field Multipliers
                            VlFldMlt[n1][n2][n3][k1][k2][k3] = mult * (REAL) (globalFactorial[k1] * globalFactorial[k2] * globalFactorial[k3]) / k_nfactorial;
                        }
                    }
                }
            }
        }
    }
}

/**************************************************************/
void Node::RecursivePanelVel(PANEL* Pan){
    Vect3 P = Position;
    
    Vect3 R = (Pan->CollocationPoint - P);
    REAL R2 = R.Dot(R);
    bool DoHere = true;

    if (R2 < (0.75*size*size)) {                //      If panel centroid is inside node, do at children
        if (Children[0][0][0]) {
            DoHere = false;
            Children[0][0][0]->RecursivePanelVel(Pan);
        }
        if (Children[0][0][1]) {
            DoHere = false;
            Children[0][0][1]->RecursivePanelVel(Pan);
        }
        if (Children[0][1][0]) {
            DoHere = false;
            Children[0][1][0]->RecursivePanelVel(Pan);
        }
        if (Children[0][1][1]) {
            DoHere = false;
            Children[0][1][1]->RecursivePanelVel(Pan);
        }
        if (Children[1][0][0]) {
            DoHere = false;
            Children[1][0][0]->RecursivePanelVel(Pan);
        }
        if (Children[1][0][1]) {
            DoHere = false;
            Children[1][0][1]->RecursivePanelVel(Pan);
        }
        if (Children[1][1][0]) {
            DoHere = false;
            Children[1][1][0]->RecursivePanelVel(Pan);
        }
        if (Children[1][1][1]) {
            DoHere = false;
            Children[1][1][1]->RecursivePanelVel(Pan);
        }
    }
    if (DoHere) {
        PanelVel += Pan->BodyPanelVelocity(P);
    }
}
