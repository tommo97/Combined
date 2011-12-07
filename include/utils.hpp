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



#ifndef UTILS_INCL
#define UTILS_INCL
#include "types.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <complex>
#include <matio.h>
#include <omp.h>
#include <algorithm>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>

using namespace std;



#define EPS 1e-16
#define USE_MATRIX_INVERSE
#define TEST_MODE
#define USE_ARRAY
#define ARRAY_NO_CHECK
#define USE_FULL_BLADE_LOOPS
#define RETURN_WAKES_AT_END
#define REAL double
#define TOL 1e-6
//#define SURFW(X1,X2,X4,X3,C,W) W << "surf([" << X1.x << ", " << X2.x << "; " << X3.x << ", " << X4.x << "],[" << X1.y << ", " << X2.y << "; " << X3.y << ", " << X4.y << "],[" << X1.z << ", " << X2.z << ";" << X3.z << ", " << X4.z << "]," << C << ");\n";
//#define SURF(X1,Y1,Z1,X2,Y2,Z2,X4,Y4,Z4,X3,Y3,Z3,C,W) W << "surf([" << X1 << ", " << X2 << "; " << X3 << ", " << X4 << "],[" << Y1 << ", " << Y2 << "; " << Y3 << ", " << Y4 << "],[" << Z1 << ", " << Z2 << ";" << Z3 << ", " << Z4 << "]," << C << ");\n";
//#define pi 3.141592653589793
#define piover2 1.570796326794897
//#define two_pi 6.283185307179586
//#define four_pi 12.56637061435917

#include "array.hpp"
#include "vect34.hpp"
#include "includes.hpp"

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
class BODY;
class PANEL;
extern unsigned long int LineVelCnt;
extern int NumThreads;

/**************************************************************/
class UTIL {
public:

    static unsigned long int cpu_t;
    static inline Array < Array <double> > zeros(int a, int b) {
        return Array < Array <double> > (a, Array <double> (b, 0.0));
    }

    static inline Array < Array <Array <double> > > zeros(int a, int b, int c) {
        return Array < Array <Array <double> > > (a, Array <Array <double> > (b, Array <double> (c, 0.0)));
    }

    static inline Array < Array <Vect3> > zerosv(int a, int b) {
        return Array < Array <Vect3> > (a, Array <Vect3 > (b, 0.0));
    }

    template<class T>
    inline static string toString(const T& t) {
        stringstream ss;
        ss << t;
        return ss.str();
    }

    static double FastInverseSqrt(double);

    static int read_neu(string infname,
            Array<Vect3> &X,
            Array<Array<int> > &PNLS,
            Array<Array<int> > &GROUPS,
            Array<Array<int> > &BCS,
            Array<string> &NAMES,
            Array <Array < Array < int > > > &Surfaces,
            Array < Array <REAL> > &CRLocal,
            Array < Array < Array < int > > > &PtIDS,
            Array < Array <Array <int> > > &TipInboardUSIDS,
            Array < Array <Array <int> > > &TipInboardLSIDS,
            Array < Array <Array <int> > > &TipOutboardUSIDS,
            Array < Array <Array <int> > > &TipOutboardLSIDS,
            Array < Array <Array <int> > > &InnerTipUSPanelIDS,
            Array < Array <Array <int> > > &InnerTipLSPanelIDS,
            Array < Array <Array <int> > > &OuterTipUSPanelIDS,
            Array < Array <Array <int> > > &OuterTipLSPanelIDS);


    static void write2D(string varname, string fname, Array<Array<double> > &input, int m,
            int n);
    static void write2D(string varname, string fname, Array<Array<int> > &input, int m,
            int n);
    static void write1D(string varname, string fname, Array<double> &input, int m);
    static void write1D(string varname, string fname, Array<int> &input, int m);
    static void WriteMATLABMatrix2D(string vname, string fname,
            Array<Array<int> > &data);
    static void WriteMATLABMatrix2D(string vname, string fname,
            Array<Array<double> > &data);
    static void WriteMATLABMatrix2DVect3(string vname, string fname,
            Array<Array<Vect3> > &data);
    static void WriteMATLABMatrix1D(string vname, string fname, double data);
    static void WriteMATLABMatrix1D(string vname, string fname, Array<double> &data);
    static void WriteMATLABMatrix1D(string vname, string fname, Array<int> &data);
    static void WriteMATLABMatrix1DVect3(string vname, string fname,
            Array<Vect3> &data);
    static void WriteMATLABMatrix1DVect3(string vname, string fname,
            Vect3 &data);

    static void write1D(string varname, string fname, string &input, int m);
    static void WriteMATLABString(string vname, string fname, string data);
    
    static void ReadBinaryVect3(Array <Vect3> &, string);

    inline static double rad2deg(double theta) {
        return theta / 0.017453292519943;
    }

    inline static double deg2rad(double theta) {
        return theta * 0.017453292519943;
    }


    static double interp2(Array<Array<double> > &X, Array<Array<double> > &Y, Array<Array<
            double> > &Z, double Xi, double Yi);

    static double interp3(Array < Array<Array<double> > > &X,  Array < Array<Array<double> > > &Y, Array <Array<Array<
            double> > > &Z, double Xi, double Yi, double Zi);
    
    static Array <REAL> globalLinspace(REAL start, REAL end, int n);
    static Array <Vect3> globalLinspace(Vect3 start, Vect3 end, int n);
    static Vect3 globalDirectVel(Vect3 diff, Vect3 omega, REAL del2);
    static void PreAmble();
    static void PostAmble(string);
};
/**************************************************************/
 
#endif	/* if !defined(UTILS_INCL) */
