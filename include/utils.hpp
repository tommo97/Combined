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


#include "array.hpp"
#include "vect34.hpp"
#include "includes.hpp"

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
