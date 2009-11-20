/*
 This file is part of the Combined Wake Modelling Code Version 1.0

 VTM Code Copyright Tom McCombes 2009
 This code solves the 3D unsteady incompressible
 Navier-Stokes equations in velociy vorticity form


 $Rev:: 28               $:  Revision of last commit
 $Author:: tom           $:  Author of last commit
 $Date:: 2009-11-09 15:0#$:  Date of last commit

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
	Node *Parent;
	int x, y, z, m;
	long unsigned int ID;
	REAL size;
	Array<bool> skip_here;
	bool InList;
	Vect3 Position, Omega, Velocity;
	NeighbSet<Node*> Neighb;
	Array<NeighbSet<Vect3> > BEV;
	Array<NeighbSet<Vect3*> > Neighb_Val, Neighb_Neighb_Val;
	bool HasLoad;
	bool to_report;
	Node *Children[2][2][2];
	Node *ISA[3][3][3];
	Array<int> Trans;
	static int RootSize;
	static const Vect3 Offset[2][2][2];
	static const Vect3 NeighbOffset[6];
	static const int deREF[][3];

	Array<Vect3> TransVars, TransDerivs;

	void UpdateMomentMults();

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


	void vPruneChild(int i, int j, int k) {
		delete Children[i][j][k];
		Children[i][j][k] = NULL;
	}


	void SetISANull() {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					ISA[i][j][k] = NULL;
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
		if (Parent->ID > 0)
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

	void RecursivePanelVel(PANEL&);

	void RecursivePassPanelVelsDown()
	{
		if (Children[0][0][0]){
			Children[0][0][0]->Velocity += Velocity;
			Children[0][0][0]->RecursivePassPanelVelsDown();
		}
		if (Children[0][0][1]) {
			Children[0][0][1]->Velocity += Velocity;
			Children[0][0][1]->RecursivePassPanelVelsDown();
		}
		if (Children[0][1][0]) {
			Children[0][1][0]->Velocity += Velocity;
			Children[0][1][0]->RecursivePassPanelVelsDown();
		}
		if (Children[0][1][1]) {
			Children[0][1][1]->Velocity += Velocity;
			Children[0][1][1]->RecursivePassPanelVelsDown();
		}
		if (Children[1][0][0]) {
			Children[1][0][0]->Velocity += Velocity;
			Children[1][0][0]->RecursivePassPanelVelsDown();
		}
		if (Children[1][0][1]) {
			Children[1][0][1]->Velocity += Velocity;
			Children[1][0][1]->RecursivePassPanelVelsDown();
		}
		if (Children[1][1][0]) {
			Children[1][1][0]->Velocity += Velocity;
			Children[1][1][0]->RecursivePassPanelVelsDown();
		}
		if (Children[1][1][1]) {
			Children[1][1][1]->Velocity += Velocity;
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


	static ARRAY13(Vect3) TlrCffts;
	static ARRAY6(REAL) VlFldMlt;
	static ARRAY10(Vect3) TlrCfftsdx;
	static ARRAY6(NeighbSet <REAL>) FaceCllpsMlt;
	static ARRAY10(REAL) BinomMlt;
	static ARRAY10(REAL) InhrtMlt;
	static ARRAY6(REAL) CllpsMlt;
	static int REF[][2][2];
	static int Op[6];
	static const int trans_neighb[][5];
	long unsigned int NID;
	static Vect3 ZERO;

};

#endif /* if !defined(NODE_INCL) */
