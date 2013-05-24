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


#include "array.hpp"
#include "vect34.hpp"
#include "includes.hpp"

/**************************************************************/
extern unsigned long int LineVelCnt;
extern int NumThreads;

/**************************************************************/
class UTIL {
public:

    static unsigned long int cpu_t;
    static int NumCellGaussPts;
    static inline Array < Array <double> > zeros(int a, int b) {
        return Array < Array <double> > (a, Array <double> (b, 0.0));
    }

    static inline Array < Array <Array <double> > > zeros(int a, int b, int c) {
        return Array < Array <Array <double> > > (a, Array <Array <double> > (b, Array <double> (c, 0.0)));
    }

    static inline Array < Array <Vect3> > zerosv(int a, int b) {
        return Array < Array <Vect3> > (a, Array <Vect3 > (b, 0.0));
    }


    template <typename T> static ARRAY3(T) zeros(int a, int b, int c) {
        return ARRAY3(T) (a, ARRAY2(T) (b, Array <T> (c, 0.0)));
    }

    static inline Array < Array < Array <Vect3> > > zerosv(int a, int b, int c) {
        return Array < Array <Array <Vect3> > > (a, Array <Array <Vect3> > (b, Array <Vect3> (c, Vect3(0.0))));
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
    static void write1D(string varname, string fname, Array<long unsigned > &input, int m);

    static int readmat(string fname, string varname, Array <REAL> &data, Array <int> &dims, bool verbose);
    static void WriteMATLABMatrix2D(string vname, string fname,
            Array<Array<int> > &data);
    static void WriteMATLABMatrix2D(string vname, string fname,
            Array<Array<double> > &data);
    static void WriteMATLABMatrix2DVect3(string vname, string fname,
            Array<Array<Vect3> > &data);
    static void WriteMATLABMatrix1D(string vname, string fname, double data);
    static void WriteMATLABMatrix1D(string vname, string fname, Array<double> &data);
    static void WriteMATLABMatrix1D(string vname, string fname, Array<int> &data);
    static void WriteMATLABMatrix1D(string vname, string fname, Array<long unsigned int> &data);

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
    static Array <PANEL> Pans;
    static Array <Vect3> globalLinspace(Vect3 start, Vect3 end, int n);
    static Vect3 globalDirectVel(Vect3 diff, Vect3 omega);
    static Vect3 globalCubicDirectVel(Vect3 diff, Vect3 omega);
    static void globalDirectVelGrads(Vect3 diff, Vect3 omega, Array <Vect3> &Grads);
    static void globalCubicDirectVelGrads(Vect3 diff, Vect3 omega, Array <Vect3> &Grads);
    static void PreAmble();
    static void PostAmble(string);
    static void GetCellPans();
    static Array <REAL> QuadPts, QuadWts;
    static void lgwt(int N, Array <REAL> &x, Array <REAL> &w);
    template <typename T> static T interp3(ARRAY3(Vect3) &X, ARRAY3(T) &U, Vect3 Xi);
    
    
    
};
/**************************************************************/
 
#endif	/* if !defined(UTILS_INCL) */
