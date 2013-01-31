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

#ifndef NODE_INCL
#define NODE_INCL

#include "includes.hpp"
#include "types.hpp"


/**************************************************************/
#define callMemberFunction(object,ptrToMember)  ((object).*(ptrToMember))

class Node {
public:

    virtual ~Node();

    Node();

    Node(Node *parent, int i, int j, int k);
    static Node *Root;
    static unsigned long int NumNodes, NodeCount;
    static int MomentSize;
    static Array <Array <int> > CellMomentIndex;
    static Array < Array <Array <REAL> > > CellMomentDisplacementPowers;
    static Array <Node*> AllNodes, UpList, DownList;
    static Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > DirVelMultsX;
    static Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > DirVelMultsY;
    static Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > DirVelMultsZ;
    
    static Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > DirGradMultsX;
    static Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > DirGradMultsY;
    static Array < Array <Array <Array < Array < Array < Array < Array < Array <Array < Vect3 > > > > > > > > > > DirGradMultsZ;
    
    static Array < Array <int> > ISBRecipInds;
    static Array <int> ISARecipInds;

    static Array < Array < Array <Vect3> > > ISBGradMultsX;
    static Array < Array < Array <Vect3> > > ISBGradMultsY;
    static Array < Array < Array <Vect3> > > ISBGradMultsZ;
    static Array < Array <Vect3> > ISBDirMultsX;
    static Array < Array <Vect3> > ISBDirMultsY;
    static Array < Array <Vect3> > ISBDirMultsZ;

    static Array < Array <Vect3> > ISAGradMultsX;
    static Array < Array <Vect3> > ISAGradMultsY;
    static Array < Array <Vect3> > ISAGradMultsZ;
    static Array <Vect3> ISADirMultsX;
    static Array <Vect3> ISADirMultsY;
    static Array <Vect3> ISADirMultsZ;
    
    static const Vect3 Offset[2][2][2];
    static const int Indxs[2][2][2];
    static const Vect3 NeighbOffset[6];
    static const int deREF[][3];
    static ARRAY13(Vect3) TlrCffts;
    static ARRAY6 (Vect3) LinTlrCffts;
    static ARRAY6(REAL) VlFldMlt;
    static ARRAY10(Vect3) TlrCfftsdx;
    static ARRAY6(NeighbSet <REAL>) FaceCllpsMlt;
    static ARRAY10(REAL) BinomMlt;
    static ARRAY10(REAL) InhrtMlt;
    static ARRAY6(REAL) CllpsMlt;
    static int REF[][2][2];
    static int Op[6];
    static const int trans_neighb[][5];
    static Vect3 ZERO;
    static int RootSize;

    
    //  Node class definately needs
    Node *Parent;
    int m, indx;
    long unsigned int ID;
    REAL size;
    Array<bool> skip_here;
    bool SkipAdd2List;

    bool HasLoad;
    bool to_report;
    Node *Children[2][2][2];
    Node *ISA[3][3][3]; // this should replace Neighb?
    Array <Node*> LinISB, LinISA;
//    Node *ISB[6][6][6];
    Array<int> Trans;
    
    
    
    //  Not really needed any more
    int x, y, z;
    
    
    //   possibly devolve the following members solely to cells
    Vect3 Position, Omega, Velocity, PanelVel;
    NeighbSet<Node*> Neighb;
    Array<NeighbSet<Vect3> > BEV;
    Array<NeighbSet<Vect3*> > Neighb_Val, Neighb_Neighb_Val;
    Array<Vect3> TransVars, TransVarsHold, TransVars0;
    Array <Array<Vect3> > TransDerivs;
    
    
    void UpdateMomentMults();
    
    void SetUpISBIndices();

    void GetISA();

    void Prune() {
        HasLoad = false;
        Omega = 0.;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                    if (Children[i][j][k]) {
                        if (!Children[i][j][k]->HasLoad) {
                            delete Children[i][j][k];
                            Children[i][j][k] = NULL;
                        } else {
                            HasLoad = true;
                            Omega += Children[i][j][k]->Omega;
                        }
                    }
    }
    
    
    void GetISB();
    
    void MarkWithoutLoad() {
        HasLoad = false;
    }

    void MarkWithLoad() {
        HasLoad = true;
        if (Parent && !Parent->HasLoad)
            Parent->MarkWithLoad();
    }

    void CheckLoad() {

        HasLoad = false;

        //        if (Omega.Mag() > VORTICITY_CUTOFF) {
        //            MarkWithLoad();
        //            return;
        //        }

        for (int q = 0; q < TransVars.size(); ++q)
            if (TransVars[q].Dot(TransVars[q]) > (VORTICITY_CUTOFF*VORTICITY_CUTOFF)) {
                MarkWithLoad();
                return;
            }



    }
    void vPruneChild(int i, int j, int k) {
        delete Children[i][j][k];
        Children[i][j][k] = NULL;
    }

    void SetISANull() {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    ISA[i][j][k] = NULL;
            
            LinISB = Array <Node*> (216,NULL);
            LinISA = Array <Node*> (27,NULL);
    }

    void SetKidsNull() {
        Children[0][0][0] = Children[0][1][0] = Children[0][0][1]
                = Children[0][1][1] = NULL;
        Children[1][0][0] = Children[1][1][0] = Children[1][0][1]
                = Children[1][1][1] = NULL;
    }

    void ClearChildren() {
        delete Children[0][0][0];
        delete Children[0][0][1];
        delete Children[0][1][0];
        delete Children[0][1][1];
        delete Children[1][0][0];
        delete Children[1][0][1];
        delete Children[1][1][0];
        delete Children[1][1][1];
        SetKidsNull();
    }


    Node* MoveDownTree(Array<int> &Directions);

    void MoveUpTree() {
        if (Parent->m > 0)
            Parent->MoveUpTree();
    }

    void CompCoeffts(Vect3 diff, JaggedArray<REAL> &);
    void Recurrance(JaggedArray<REAL> &, int, int, int, Vect3, REAL);

    //    void GetVelField();

    void PrintNeighb() {
        if (WRITE_TO_SCREEN)
            cout << Position << " ";
        for (int i = 0; i < 6; ++i)
            if (Neighb[i])
                if (WRITE_TO_SCREEN)
                    cout << Neighb[i]->Position << " ";
        if (WRITE_TO_SCREEN)
            cout << endl;
    }

    void DoNothing() {
    }

    void ReportSomething() {
    }


    void RemoveFromNeighbs();

    void Trans2Neighb(Array<int> &v, int dirn);

    Node* ReturnNeighb(int, int, int);

    virtual void vEvalCapsule(OctreeCapsule&) = 0;

    void EvalCapsule(OctreeCapsule &c);

    void RecursivePanelVel(PANEL*);

    void RecursivePassPanelVelsDown() {
        if (Children[0][0][0]) {
            Children[0][0][0]->PanelVel += PanelVel;
            Children[0][0][0]->RecursivePassPanelVelsDown();
        }
        if (Children[0][0][1]) {
            Children[0][0][1]->PanelVel += PanelVel;
            Children[0][0][1]->RecursivePassPanelVelsDown();
        }
        if (Children[0][1][0]) {
            Children[0][1][0]->PanelVel += PanelVel;
            Children[0][1][0]->RecursivePassPanelVelsDown();
        }
        if (Children[0][1][1]) {
            Children[0][1][1]->PanelVel += PanelVel;
            Children[0][1][1]->RecursivePassPanelVelsDown();
        }
        if (Children[1][0][0]) {
            Children[1][0][0]->PanelVel += PanelVel;
            Children[1][0][0]->RecursivePassPanelVelsDown();
        }
        if (Children[1][0][1]) {
            Children[1][0][1]->PanelVel += PanelVel;
            Children[1][0][1]->RecursivePassPanelVelsDown();
        }
        if (Children[1][1][0]) {
            Children[1][1][0]->PanelVel += PanelVel;
            Children[1][1][0]->RecursivePassPanelVelsDown();
        }
        if (Children[1][1][1]) {
            Children[1][1][1]->PanelVel += PanelVel;
            Children[1][1][1]->RecursivePassPanelVelsDown();
        }
    }

    virtual void
    vApplyRecursively(BranchFuncPtr, FVMCellFuncPtr, BranchFuncPtr) = 0;

    void ApplyRecursively(BranchFuncPtr down, FVMCellFuncPtr bottom,
            BranchFuncPtr up) {
        vApplyRecursively(down, bottom, up);
    }


    virtual void
    vApplyRecursivelyP(BranchFuncPtr, FVMCellFuncPtr, BranchFuncPtr) = 0;

    void ApplyRecursivelyP(BranchFuncPtr down, FVMCellFuncPtr bottom,
            BranchFuncPtr up) {
        vApplyRecursivelyP(down, bottom, up);
    }


    virtual void vPassMmnts2Prnt() = 0;

    void PassMmnts2Prnt() {
        vPassMmnts2Prnt();
    }


    virtual void vCollapseVField() = 0;

    void CollapseVField() {
        vCollapseVField();
    }


    virtual void vCheckNeighbs() = 0;

    void CheckNeighbs() {
        vCheckNeighbs();
    }


    virtual void vMakeNodeAtTrans(Array < int > &trn) = 0;

    void MakeNodeAtTrans(Array < int > &trn) {
        vMakeNodeAtTrans(trn);
    };
    
    virtual void vReList() = 0;

    void ReList() {
        Node::AllNodes[(int) Node::NodeCount] = this;
        Node::NodeCount++;
        vReList();
    }

    

    void Trans2Neighb(Array <int> &v, int M, int dirn) {
        /* Find the neighbours of an octree cluster in a given direction
            Binary transformations are as follows:
            n 1   2   3   4   5   6   7   8
            x 0   0   1   1   0   0   1   1
            y 0   1   0   1   0   1   0   1
            z 0   0   0   0   1   1   1   1
            WSB WNB ESB ENB WST WNT EST ENT
         */

        int a = 0;

        a += (v[M] == trans_neighb[dirn][0]);
        a += (v[M] == trans_neighb[dirn][1]);
        a += (v[M] == trans_neighb[dirn][2]);
        a += (v[M] == trans_neighb[dirn][3]);
        if (a != 0) {
            v[M] += trans_neighb[dirn][4];
        } else {
            v[M] -= trans_neighb[dirn][4];
            Trans2Neighb(v, M - 1, dirn);
        }
    };

};

#endif /* if !defined(NODE_INCL) */

