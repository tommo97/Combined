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


#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "includes.hpp"
#include "types.hpp"
#include "io.hpp"
#include "tree.hpp"
#include "body.hpp"


class SYSTEM
{
public:
  int MaxP, NumBodyPanels, ProcessID, NumBodies, NumThreads, NumSubSteps, SysDumpInterval, NumTransVars;
  bool LiftingLineMode, ZeroBodyRelativeMotion, useBodies;
  REAL dtInit, GambitScale, Del2, DS, h;
  REAL Mu,Nu,Rho,Temp,g;
  string NeuFile, CaseName, WorkingDir;
  Vect3 scaledVinf, unscaledVinf;
  SYSTEM(int);
  ~SYSTEM();
  bool PanelMode;
  int num_out;
  string InputStr;

  Array <Vect3> ATTITUDE, VELOCITY, ORIGIN, RATES;
  Array <Vect3> VortonsXs, VortonOmegas;

#ifdef USEGSL
  gsl_matrix * globalA;     //  A (doublet) influence coefficient matrix for this body
  gsl_matrix * globalB;     //  B (source) influence coefficient matrix for this body
  gsl_matrix * panelU;
  gsl_matrix * panelV;
  gsl_matrix * panelW;
  gsl_permutation * globalP;
  gsl_vector * globalMu;
  gsl_vector * globalSigma;
  gsl_vector * globalRHS;
#else
  REAL **A;
  REAL **panelU;
  REAL **panelV;
  REAL **panelW;
  REAL *rhs;
  REAL *mu;
  REAL *sigma;

  REAL **B;
  REAL **B_x_axis;
  REAL **B_y_axis;
  REAL **B_z_axis;
  int *ipiv;
#endif
	
  void Initialise();
  void ReadNeuGetBodies();
  void SetupGlobalInfluenceMatrices();
  void UpdateGlobalInfluenceMatrices();
  void BodySubStep(REAL, int);
  void GetGlobalRHS();
  void PrintInfluenceMatrices();
  void LinAlg();
  void PrintBodiesAndWakes();
  void DecomposePanelsIntoTree();
  void InitialiseMembers();
  void TimeStep();
  void PrintBodiesVels();
  void WriteBodiesAndWakes(ostream& out_stream);
  void GetCellsCalcCurlV();
  void GetPressures(REAL);
  void PutWakesInTree();
  void AddVortonsToTree(Array <Vect3> &, Array <Vect3> &, Array <int> &);
  void GetFaceVels();
  void GetPanelFMMVelocities(REAL);
  void LinearWaves();
  void MoveBodies(REAL, bool);
  void WriteDomain();
  void WriteData();
  void WritePanelVels();
  void WriteBodies();
    enum exception {
        GENERAL_ERROR
    };

};

#endif
