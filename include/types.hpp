/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 35               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-16 00:1#$:  Date of last commit

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



#ifndef TYPES_INCL
#define TYPES_INCL



#define STRINGIZE(X) #X
#define StringFromMakefile(X) (STRINGIZE(X))

#define WRITE_TO_SCREEN true
#define WRITE_TO_FILE true

//  Enable GNU Scientific Library linear algebra - uses GSL matrix/Array classes. Using LAPACK/ATLAS is faster...
//#define USEGSL
//#define RECURSE
using namespace std;
#define MODE_3D
//#define _PNGWRITER
//  Precision
#define DOUBLE_PRECISION

#define OCTREE_SIZE 4096
#define OCTREE_LEVS 13
#define LIMITER Koren

//#define COLLAPSE_TO_FACES

#define NX 300
#define NY 300
//  Disables bounds checking in Array class... use with caution. Make sure code runs, and doesn't change...
#define ARRAY_NO_CHECK
//#define USE_ROLLED_LOOPS

//  Some simulation parameters
#define USE_ARRAY
//#define TIME_STEPS


/**************************************************************/

class Branch;
class FVMCell;
class Node;
class OctreeCapsule;



/**************************************************************/
// Some multi-dimensional array macros
#define ARRAY2(type) Array < Array < type > >
#define ARRAY3(type) ARRAY2(Array < type >)
#define ARRAY4(type) ARRAY3(Array < type >)
#define ARRAY5(type) ARRAY4(Array < type >)
#define ARRAY6(type) ARRAY5(Array < type >)
#define ARRAY7(type) ARRAY6(Array < type >)
#define ARRAY8(type) ARRAY7(Array < type >)
#define ARRAY9(type) ARRAY8(Array < type >)
#define ARRAY10(type) ARRAY9(Array < type >)
#define ARRAY11(type) ARRAY10(Array < type >)
#define ARRAY12(type) ARRAY11(Array < type >)
#define ARRAY13(type) ARRAY12(Array < type >)
/**************************************************************/

#define USE_MATRIX_INVERSE

#ifdef DOUBLE_PRECISION
#define REAL double
#define _EPS 1e-3
#define _EPS2 1e-6
#define _TOL 1e-12
#define VORTICITY_CUTOFF 1e-6
#else
#define REAL float
#define _EPS 1e-3
#define _EPS2 1e-6
#define _TOL 1e-6
#define VORTICITY_CUTOFF 1e-3
#endif

#include "includes.hpp"
#include "array.hpp"
#include "vect34.hpp"


//  Tells code to include the tilting term (not required for 2D cases)
// #define INCLUDE_TILTING_TERM            

//  MINMOD3 parameter: 1 results in MIN-MOD limiter, 2 results in SUPERBEE           
#define BETA 2.0
//  Defines whether to use teams mode for openMP
//  #define TEAMS_MODE
//  Use when just using an input file (patch)            
// #define TEST_WITHOUT_BODIES
//  Use to perform an additional velocity calc at the cell centre            
// #define COLLAPSE_TO_CELL_CENTRE         

//  Change to QUIET or DEBUG for more verbosity - can slow things down and result in LOTS of screen output            
#define SILENT        
//  Use Ncurses for display output. Don't use on HPC or when dumping screen to a file (otherwise output will be mince)     
//#define use_NCURSES

//  how many output lines are printed before a header is printed to screen                
#define HEADER_OUTPUT 25

//  display CPU/memory information during run
#define TOP                        
/**************************************************************/
//  Some constants that should probably not be fiddled with...   
#define pi 3.1415926535897932384626433832795
#define two_pi 6.2831853071795864769252867665590
#define four_pi 12.566370614359172953850573533118


#define SUPERBEE(r) max(max(min(REAL(2.0)*r, REAL(1.0)), min(r, REAL(1.0))), REAL(0.0))

#define MINMOD(r) max(REAL(0.0), min(REAL(1.0), r))

#define VANLEER(r) (r+abs(r))/(1+abs(r))

#define Koren(r) max(0.0,min(min(2.0*r,(2.0+r)/3.0),2.0))

#define flim(r) Vect3(LIMITER(r.x), LIMITER(r.y), LIMITER(r.z))

#define MINMOD2(x,y) 0.5*(sign(x) + sign(y))*min(fabs(x),fabs(y))

#define MINMOD3(x,y,z) MINMOD2(x,MINMOD2(y,z))


#define NVARS 3
#define CELL_ISA_NUMBER 26
#define NFACES 6            
// number of transport variables
//#define N 8
// #ifndef MODE_IS_2D

// #else
// #define CELL_ISA_NUMBER 8
// #define NFACES 4
// #endif



#ifdef _OPENMP
#define running std::endl <<  "\t\t\t- \033[1;32mrunning openmp if enabled\033[0m"
#endif
#ifndef _OPENMP
#define running std::endl << "\t\t\t- \033[1;31mnot running openmp\033[0m"
#endif




#ifdef USE_SQRT_APPROX
#define SQRT fast_sqrt_approx
#else
#define SQRT sqrt
#endif                    

/* REAL is the floating point type to be used. If double is used, then fast approximations for atan2, invsqrt and SQRT cannot be used */
#define BIGINT long unsigned int

#define DEG2RAD(a) two_pi*a/360

#define RAD2DEG(a) 360*a/two_pi

#define acosd(a) RAD2DEG(acos(a))

#define asind(a) RAD2DEG(asin(a))

#define atand(a) RAD2DEG(atan(a))

#define sind(a) sin(DEG2RAD(a))

#define cosd(a) cos(DEG2RAD(a))

#define tand(a) tan(DEG2RAD(a))

#define SURF(X1,X2,X4,X3,C) if (WRITE_TO_SCREEN) std::cout << "surf([" << X1[0] << ", " << X2[0] << "; " << X3[0] << ", " << X4[0] << "],[" << X1[1] << ", " << X2[1] << "; " << X3[1] << ", " << X4[1] << "],[" << X1[2] << ", " << X2[2] << ";" << X3[2] << ", " << X4[2] << "]," << C << ");\n";

#define SURFW(X1,X2,X4,X3,C,W) W << "surf(plot_ax,[" << X1[0] << ", " << X2[0] << "; " << X3[0] << ", " << X4[0] << "],[" << X1[1] << ", " << X2[1] << "; " << X3[1] << ", " << X4[1] << "],[" << X1[2] << ", " << X2[2] << ";" << X3[2] << ", " << X4[2] << "]," << C << ");\n";

        

/**************************************************************/
typedef void (Branch::*BranchReportPtr)(ostream &);
typedef void (FVMCell::*FVMCellReportPtr)(ostream &);
typedef void (Branch::*BranchFuncPtr)();
typedef void (FVMCell::*FVMCellFuncPtr)();

/**************************************************************/

#ifdef _LP64
#define VERS "64 bit"
#else
#define VERS "32 bit"
#endif

/**************************************************************/
/* Some Pre-definitions */

class LINE;

class POINT;

class FVMCell;

class BODY;

class PANEL;

class IO;

class TIME_STEPPER;

class OCTREE;

class SYSTEM;

/**************************************************************/
/* global variable declarations */
extern long unsigned int globalFactorial[24];

extern SYSTEM *globalSystem;

extern IO *globalIO;

extern TIME_STEPPER *globalTimeStepper;

extern OCTREE *globalOctree;

extern long int globalTimer;

/*
*  Global Functions
*/
long int ticks();

Array <REAL> globalLinspace(REAL start, REAL end, int n);

Array <Vect3> globalLinspace(Vect3 start, Vect3 end, int n);

string globalGetStdoutFromCommand(string cmd);

/**************************************************************/
static struct _ticks_global {
struct timeval start;

_ticks_global() {
	gettimeofday(&start, NULL);
}
} ticks_global;

/**************************************************************/

class POINT
{
    public:
	BODY *Owner;
	Vect3 vP, vV, vVfmm[2], vO, vVstar, vXstar;
	POINT(Vect3 POS, Vect3 VEL, Vect3 OMEGA) : vP(POS), vV(VEL), vO(OMEGA){};
	POINT(Vect3 POS, Vect3 OMEGA) : vP(POS), vO(OMEGA){};
	POINT(Vect3 POS) : vP(POS){};
	POINT(){};
	void PrintPos() { if (WRITE_TO_SCREEN) std::cout << vP << std::endl;}
	void PrintP() { if (WRITE_TO_SCREEN) std::cout << vP << std::endl;}
	void PrintV() { if (WRITE_TO_SCREEN) std::cout << vV << std::endl;}
	void PrintO() { if (WRITE_TO_SCREEN) std::cout << vO << std::endl;}
};
/**************************************************************/
class OctreeCapsule
{
public:
    int S, AssociatedBody;
    Vect3 Position, Omega, tPosition, Velocity;

    bool has_load, IP, toMonitor;

    OctreeCapsule() : S(OCTREE_SIZE / 2), has_load(false), IP(false), toMonitor(false) {
    };

    OctreeCapsule(Vect3 P, Vect3 O, bool t) : S(OCTREE_SIZE / 2), Position(P), Omega(O), tPosition(P), has_load(t), IP(false), toMonitor(false) {
    };

    OctreeCapsule(const OctreeCapsule &A) : S(A.S), Position(A.Position), Omega(A.Omega), tPosition(A.tPosition), has_load(A.has_load), IP(A.IP), toMonitor(A.toMonitor) {
    }
    inline bool operator>(const OctreeCapsule B) {
        return Position > B.Position;
    }

    inline bool operator<(const OctreeCapsule B) {
        return Position < B.Position;
    }

    inline bool operator >=(const OctreeCapsule B) {
        return Position >= B.Position;
    }

    inline bool operator <=(const OctreeCapsule B) {
        return Position <= B.Position;
    }

    inline bool operator ==(const OctreeCapsule B) {
        return Position == B.Position;
    }

    inline bool operator !=(const OctreeCapsule B) {
        return Position != B.Position;
    }
};


/**************************************************************/
void globalDirectVel(Vect3 diff, Vect3 omega, Vect3 & vel);
Vect3 globalDirectVel(Vect3 diff, Vect3 omega);

#endif /* if !defined(TYPES_INCL) */
