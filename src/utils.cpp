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
#include "utils.hpp"
#include "panel.hpp"

unsigned long int LineVelCnt = 0;
int NumThreads = 1;

long unsigned int UTIL::cpu_t = 0;
int UTIL::NumCellGaussPts;
Array <REAL> UTIL::QuadPts, UTIL::QuadWts;
Array <PANEL> UTIL::Pans;
/**************************************************************/
/*      This piece of code is for fast reciprocal square root approximation as used in Quake*/
double UTIL::FastInverseSqrt( double number )
{
        long i;
        double x2, y;
        const double threehalfs = 1.5F;
 
        x2 = number * 0.5F;
        y  = number;
        i  = * ( long * ) &y;                       // evil floating point bit level hacking
        i  = 0x5fe6eb50c7aa19f9 - ( i >> 1 );       // what the fuck?
        y  = * ( double * ) &i;
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
 
        return y;
}
/**************************************************************/
void UTIL::GetCellPans() {
    Array <Vect3> Verts(8);
    Verts[0] = 0.5 * Vect3(-1, -1, -1);
    Verts[1] = 0.5 * Vect3(-1, -1, 1);
    Verts[2] = 0.5 * Vect3(1, -1, 1);
    Verts[3] = 0.5 * Vect3(1, -1, -1);
    Verts[4] = 0.5 * Vect3(-1, 1, -1);
    Verts[5] = 0.5 * Vect3(-1, 1, 1);
    Verts[6] = 0.5 * Vect3(1, 1, 1);
    Verts[7] = 0.5 * Vect3(1, 1, -1);
    

    UTIL::Pans.push_back(PANEL(Verts[0], Verts[1], Verts[2], Verts[3]));
    UTIL::Pans.push_back(PANEL(Verts[3], Verts[2], Verts[6], Verts[7]));
    UTIL::Pans.push_back(PANEL(Verts[7], Verts[6], Verts[5], Verts[4]));
    UTIL::Pans.push_back(PANEL(Verts[4], Verts[5], Verts[1], Verts[0]));
    UTIL::Pans.push_back(PANEL(Verts[1], Verts[5], Verts[6], Verts[2]));
    UTIL::Pans.push_back(PANEL(Verts[4], Verts[0], Verts[3], Verts[7]));

    for (int i = 0; i < 6; ++i) {
        UTIL::Pans[i].GetNormal();
        UTIL::Pans[i].Sigma = UTIL::Pans[i].Gamma = UTIL::Pans[i].Mu = 1.0;
    }
}
/**************************************************************/
Vect3 UTIL::globalDirectVel(Vect3 diff, Vect3 omega) {

    REAL mult, nrm;
    nrm = sqrt(SYSTEM::Del2 + diff.Dot(diff));
    mult = -1. / (four_pi * nrm * nrm * nrm);
    return mult * diff.Cross(omega);
}
/**************************************************************/
Vect3 UTIL::globalCubicDirectVel(Vect3 diff, Vect3 omega) {
   if (diff.Dot(diff) > 4) {
        Vect3 VelQ(0., 0., 0.);
        for (int i = 0; i < UTIL::QuadPts.size(); ++i)
            for (int j = 0; j < UTIL::QuadPts.size(); ++j)
                for (int k = 0; k < UTIL::QuadPts.size(); ++k)
                    VelQ += UTIL::QuadWts[i] * UTIL::QuadWts[j] * UTIL::QuadWts[k] *
                        UTIL::globalDirectVel(diff - Vect3(UTIL::QuadPts[i], UTIL::QuadPts[j], UTIL::QuadPts[k]), omega);
        return VelQ; 
    } else {
        Vect3 VelP(0., 0., 0.);
        for (int i = 0; i < 6; ++i) {
            REAL Phi = UTIL::Pans[i].HyperboloidSourcePhi(diff);
            VelP -= UTIL::Pans[i].TRANS[2].Cross(omega) * Phi / 2.0;
        }
        return VelP;
    }
}

/**************************************************************/
void UTIL::globalDirectVelGrads(Vect3 diff, Vect3 omega, Array <Vect3> &Grads) {
    //  Diff here is Xtarget - Xsource (hence -ve (-=) below)
    REAL R, R3, R5;
    
    R = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z + SYSTEM::Del2);
    R3 = R*R*R;
    R5 = R3*R*R;
    
    Grads[0].x  -= -(3.0*diff.x*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5); 
    Grads[0].y  -= (3.0*diff.x*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5) - omega.z/(four_pi*R3);
    Grads[0].z  -= omega.y/(four_pi*R3) - (3.0*diff.x*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    Grads[1].x  -= omega.z/(four_pi*R3) - (3.0*diff.y*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5);
    Grads[1].y  -= (3.0*diff.y*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5);
    Grads[1].z  -= -omega.x/(four_pi*R3) - (3.0*diff.y*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    Grads[2].x  -= -omega.y/(four_pi*R3) - (3.0*diff.z*(omega.z*diff.y - omega.y*diff.z))/(four_pi*R5);
    Grads[2].y  -= omega.x/(four_pi*R3) + (3.0*diff.z*(omega.z*diff.x - omega.x*diff.z))/(four_pi*R5);
    Grads[2].z  -= -(3.0*diff.z*(omega.y*diff.x - omega.x*diff.y))/(four_pi*R5);
    

}
/**************************************************************/
void UTIL::globalCubicDirectVelGrads(Vect3 diff, Vect3 Omega, Array <Vect3> &Grads) {
    //  Diff here is Xtarget - Xsource (hence -ve (-=) below)


    if (diff.Dot(diff) > 4) {
        for (int i = 0; i < UTIL::QuadPts.size(); ++i)
            for (int j = 0; j < UTIL::QuadPts.size(); ++j)
                for (int k = 0; k < UTIL::QuadPts.size(); ++k)
                    UTIL::globalDirectVelGrads(diff - Vect3(UTIL::QuadPts[i], UTIL::QuadPts[j], UTIL::QuadPts[k]), UTIL::QuadWts[i] * UTIL::QuadWts[j] * UTIL::QuadWts[k] * Omega, Grads);
    } else {
        for (int i = 0; i < 6; ++i) {
            Vect3 V = UTIL::Pans[i].SourcePanelVelocity(diff) * two_pi / four_pi;
            Vect3 C = UTIL::Pans[i].TRANS[2].Cross(Omega);
            Grads[0] += V.x*C;
            Grads[1] += V.y*C;  // NB this is -ve because the return value for V.y in SourceVel is -ve (why?)
            Grads[2] += V.z*C;
        }
    }
}
/**************************************************************/
void UTIL::lgwt(int N, Array <REAL> &x, Array <REAL> &w)
{
/*
 * Translation of lgwt.m -- see original header below
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
*/

    REAL a = -1.0, b = 1.0;
    N = N - 1;
    int N1 = N + 1, N2 = N + 2;
    Array <REAL> xu = globalLinspace(-1.0, 1.0, N1);

    // Initial guess
    Array <REAL> y(N1);
    x.allocate(N1);
    w.allocate(N1);
    for (int i = 0; i < N1; ++i)
        y[i] = cos((2 * i + 1) * pi / (2 * N + 2))+(0.27 / N1) * sin(pi * xu[i] * N / N2);


    // Legendre-Gauss Vandermonde Matrix
    Array <Array <REAL> > L = UTIL::zeros(N1, N2);


    // Derivative of LGVM
     Array <REAL> Lp = Array <REAL> (N1,0.0);

    // Compute the zeros of the N+1 Legendre Polynomial
    // using the recursion relation and the Newton-Raphson method


// Iterate until new points are uniformly within epsilon of old points
    REAL Linf = 1.0e6;
    while (Linf > 1e-10) {

        for (int i = 0; i < N1; ++i) {
            L[i][0] = 1.0;
            Lp[i] = 0.0;
            L[i][1] = y[i];
        }

        for (int i = 0; i < N1; ++i)
            for (int k = 0; k < N; ++k)
                L[i][k + 2] = (((2 * (k + 2)) - 1) * y[i] * L[i][k + 1]-(k + 1) * L[i][k]) / (k + 2);

         
        for (int i = 0; i < N1; ++i)
            Lp[i] = N2 * (L[i][N] - y[i] * L[i][N1]) / (1.0 - y[i] * y[i]);


        Array <REAL> y0 = y;
        for (int i = 0; i < N1; ++i)
            y[i] = y0[i] - L[i][N1] / Lp[i];

        Linf = 0.0;
        for (int i = 0; i < N1; ++i)
            Linf = max(Linf, abs(y[i] - y0[i]));

    }

    x.allocate(N1);
    w.allocate(N1);
    for (int i = 0; i < N1; ++i) {
        // Linear map from[-1, 1] to [a, b]
        x[i] = (a * (1. - y[i]) + b * (1. + y[i])) / 2.0;

        // Compute the weights
        w[i] = (b-a)/((1.-y[i]*y[i])*Lp[i]*Lp[i])*(1.0*N2*N2)/(1.0*N1*N1);

    }
    
}
/**************************************************************/

int UTIL::read_neu(string infname,
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
        Array < Array <Array <int> > > &OuterTipLSPanelIDS) {
    ifstream input;
    input.open(infname.c_str());
    if (!input) {
        // #ifndef use_NCURSES
        //         if (WRITE_TO_SCREEN)
        cout << "%\t" << "Unable to open file: " << infname << endl;
        //         throw NO_FILE;
        return -1;
        // #endif
    }
    string line, WRD, fname;

    int nodeNum, elemNum, NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL;
    bool EOS = false;
    if (input.is_open()) {
        getline(input, line); //  Read control info header
        getline(input, line); //  Read filetype specification
        getline(input, line); //  Read filename
        {
            istringstream strm(line);
            strm >> fname;
#ifndef use_NCURSES
            //             if (WRITE_TO_SCREEN)
            //                cout << "%\t"<< fname << endl;
#endif
        }
        getline(input, line); //  Read generating program
        getline(input, line); //  Read date and time of generation
        getline(input, line); //  Read mesh data headers
        getline(input, line); //  Read mesh data
        {
            istringstream strm(line);
            strm >> NUMNP >> NELEM >> NGRPS >> NBSETS >> NDFCD >> NDFVL;
        }

        Array<int> ptemp;
        ptemp.assign(4, 0);
        X.assign(NUMNP, Vect3());
        PNLS.assign(NELEM, ptemp);

        getline(input, line); //  Read ENDOFSECTION
        getline(input, line); //  Read NODAL COORDINATES 2.4.6
        //  Read NODAL coordinates
        while ((!input.eof()) && (!EOS)) {
            getline(input, line);
            istringstream strm(line);
            // READ LINE. If line is ENDOFSECTION then set EOS = true
            strm >> WRD;
            if (WRD == "ENDOFSECTION")
                EOS = true;
            else {
                nodeNum = atoi(WRD.c_str()) - 1;
                strm >> X[nodeNum].x >> X[nodeNum].y >> X[nodeNum].z;
            }
        }
        getline(input, line); //  Read ELEMENTS/CELL S2.4.6
        EOS = false;
        int t1, t2, x1, x2, x3, x4;
        while ((!input.eof()) && (!EOS)) {
            getline(input, line);
            istringstream strm(line);
            // READ LINE. If line is ENDOFSECTION then set EOS = true
            strm >> WRD;
            if (WRD == "ENDOFSECTION")
                EOS = true;
            else {
                elemNum = atoi(WRD.c_str()) - 1;
                strm >> t1 >> t2 >> x1 >> x2 >> x3 >> x4;
                PNLS[elemNum][0] = x1 - 1;
                PNLS[elemNum][1] = x2 - 1;
                PNLS[elemNum][2] = x3 - 1;
                PNLS[elemNum][3] = x4 - 1;
            }
        }
        string NAME;

        for (int i = 0; i < NGRPS; ++i) {
            Array<int> GRP;
            getline(input, line); //  Read ELEMENTS GROUP 2.4.6
            getline(input, line); //  Read Group ID; Number of elements, material and flags
            {
                istringstream strm(line);
                strm >> t1 >> t1 >> t1 >> NELEM;
            }

            EOS = false;
            getline(input, line); //  Read Group Name
            {
                istringstream strm(line);
                strm >> NAME;
                NAMES.push_back(NAME);
            }
            getline(input, line); //  Read 0
            while ((!input.eof()) && (!EOS)) {
                getline(input, line);
                istringstream strm(line);
                // READ LINE. If line is ENDOFSECTION then set EOS = true
                strm >> WRD;
                if (WRD == "ENDOFSECTION")
                    EOS = true;
                else {
                    //LINE.push_back(atoi(WRD.c_str()));
                    GRP.push_back(atoi(WRD.c_str()) - 1);
                    for (int j = 0; j < 9; ++j) {
                        strm >> t1;
                        if (t2 == t1)
                            break;
                        //LINE.push_back(t1 - 1);
                        GRP.push_back(t1 - 1);
                        t2 = t1;
                    }
                    //GRP.push_back(LINE);
                }

            }
            GROUPS.push_back(GRP);
        }

        //cout << "%\tAttempting to read " << NBSETS << " of groups..." << endl;
        for (int i = 0; i < 1; ++i) {

            getline(input, line); //  Read BOUNDARY CONDITIONS 2.4.6
            getline(input, line); //  Read Group Name; something, number of elements and some other stuff
            {
                istringstream strm(line);
                strm >> NAME >> t2 >> NELEM; // Read group name and number of elements
            }
#ifndef use_NCURSES
            //             if (WRITE_TO_SCREEN)
            //                cout << "%\t"<< NAME << " " << t2 << " " << NELEM << endl;
#endif
            {
                EOS = false;

                while ((!input.eof()) && (!EOS)) {
                    getline(input, line);
                    Array<int> GRP;
                    istringstream strm(line);
                    // READ LINE. If line is ENDOFSECTION then set EOS = true
                    strm >> WRD;
                    if (WRD == "ENDOFSECTION")
                        EOS = true;
                    else {
                        //LINE.push_back(atoi(WRD.c_str()));
                        GRP.push_back(atoi(WRD.c_str()) - 1);
                        strm >> t1 >> t2;
                        //LINE.push_back(t2 - 1);
                        GRP.push_back(t2 - 1);
                        BCS.push_back(GRP);
                        //GRP.push_back(LINE);
                    }

                }

            }

        }
    }


    //  Now read surface data   -   not standard NEU format hence commented
    bool cont = false;
    if (!input.eof()) {
        getline(input, line);
        istringstream strm(line);
        strm >> WRD;
        if (WRD == "/SURFACE")
            cont = true;
    }



//    InnerTipUSPanelIDS = Array < Array < Array < int > > > ();
//    OuterTipUSPanelIDS = Array < Array < Array < int > > > ();
//    InnerTipLSPanelIDS = Array < Array < Array < int > > > ();
//    OuterTipLSPanelIDS = Array < Array < Array < int > > > ();

    if (cont) {
        cout << "%\tReading connectivity..." << endl;

        Array <int> ln;
        int mode = 0;
        int sx = 0, sy = 0, bn = 0;
        {
            EOS = false;
            while (!input.eof() && (!EOS)) {
                getline(input, line);
                istringstream strm(line);
                strm >> WRD;
                if (WRD == "/%") {
                    EOS = true;
                } else {
                    if ((mode != 0) && (sx > 0)) {
                        //                        cout << "------ " << sy << " " << mode <<  endl;
                        ln = Array <int> (sy, 0);
                        for (int i = 0; i < sy; ++i)
                            strm >> ln[i];

                        if (mode == 1)
                            Surfaces.back().push_back(ln);
                        if (mode == 2)
                            InnerTipUSPanelIDS.back().push_back(ln);
                        if (mode == 3)
                            InnerTipLSPanelIDS.back().push_back(ln);
                        if (mode == 4)
                            OuterTipUSPanelIDS.back().push_back(ln);
                        if (mode == 5)
                            OuterTipLSPanelIDS.back().push_back(ln);

                        sx--;
                    }


                    if (WRD == "/P") {
                        mode = 1;
                        //                        cout << WRD << endl;
                        strm >> bn >> sx >> sy;
                        Surfaces.push_back(Array < Array <int> > ());
                        InnerTipUSPanelIDS.push_back(Array < Array <int> > ());
                        OuterTipUSPanelIDS.push_back(Array < Array <int> > ());
                        InnerTipLSPanelIDS.push_back(Array < Array <int> > ());
                        OuterTipLSPanelIDS.push_back(Array < Array <int> > ());

                    }
                    if (WRD == "/pTUSI") {
                        mode = 2;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/pTLSI") {
                        mode = 3;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/pTUSO") {
                        mode = 4;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/pTLSO") {
                        mode = 5;
                        strm >> bn >> sx >> sy;
                    }

                }
            }
        }
        if (Surfaces.size() > 1)
            cout << "%\tRead file " << infname << " got " << Surfaces.size() << " surfaces of size " << Surfaces[0].size() << " by " << Surfaces[0][0].size() << endl;
        else
            cout << "%\tRead file " << infname << " got " << Surfaces.size() << " surface of size " << Surfaces[0].size() << " by " << Surfaces[0][0].size() << endl;
    }



    //  Now read surface data   -   not standard NEU format hence commented
    cont = false;


    if (!input.eof()) {
        getline(input, line);
        istringstream strm(line);
        strm >> WRD;
        if (WRD == "/LOCALCHORD")
            cont = true;
    }


    if (cont) {
        cout << "%\tReading local chord/radius ordinates..." << endl;
//        CRLocal = Array < Array < REAL > > ();
        {
            EOS = false;
            while (!input.eof() && (!EOS)) {
                getline(input, line);
                istringstream strm(line);
                strm >> WRD;
                if (WRD == "/%") {
                    EOS = true;
                } else {
                    Array <REAL> ln(2, 0.0);
                    int l;
                    strm >> l;
                    for (int i = 0; i < 2; ++i)
                        strm >> ln[i];
                    CRLocal.push_back(ln);
                }
            }
        }
        cout << "%\tRead file " << infname << " got " << CRLocal.size() << " points" << endl;

    }
    cont = false;
    if (!input.eof()) {
        getline(input, line);
        istringstream strm(line);
        strm >> WRD;
        //        cout << "---> " << WRD << endl;
        if (WRD == "/MESHING")
            cont = true;
    }



//    PtIDS = Array < Array <Array <int> > > ();
//    TipInboardUSIDS = Array < Array <Array <int> > > ();
//    TipInboardLSIDS = Array < Array <Array <int> > > ();
//    TipOutboardUSIDS = Array < Array <Array <int> > > ();
//    TipOutboardLSIDS = Array < Array <Array <int> > > ();
    if (cont) {
        cout << "%\tReading local panel mesh points..." << endl;

        Array <int> ln;
        int mode = 0;
        int sx = 0, sy = 0, bn = 0;
        {
            EOS = false;
            while (!input.eof() && (!EOS)) {
                getline(input, line);
                istringstream strm(line);
                strm >> WRD;
                if (WRD == "/%") {
                    EOS = true;
                } else {
                    if ((mode != 0) && (sx > 0)) {
                        //                        cout << "------ " << sy << " " << mode <<  endl;
                        ln = Array <int> (sy, 0);
                        for (int i = 0; i < sy; ++i)
                            strm >> ln[i];

                        if (mode == 1)
                            PtIDS.back().push_back(ln);
                        if (mode == 2)
                            TipInboardUSIDS.back().push_back(ln);
                        if (mode == 3)
                            TipInboardLSIDS.back().push_back(ln);
                        if (mode == 4)
                            TipOutboardUSIDS.back().push_back(ln);
                        if (mode == 5)
                            TipOutboardLSIDS.back().push_back(ln);

                        sx--;
                    }


                    if (WRD == "/M") {
                        mode = 1;
                        //                        cout << WRD << endl;
                        strm >> bn >> sx >> sy;
                        PtIDS.push_back(Array < Array <int> > ());
                        TipInboardUSIDS.push_back(Array < Array <int> > ());
                        TipInboardLSIDS.push_back(Array < Array <int> > ());
                        TipOutboardUSIDS.push_back(Array < Array <int> > ());
                        TipOutboardLSIDS.push_back(Array < Array <int> > ());

                    }
                    if (WRD == "/TUSI") {
                        mode = 2;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/TLSI") {
                        mode = 3;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/TUSO") {
                        mode = 4;
                        strm >> bn >> sx >> sy;
                    }
                    if (WRD == "/TLSO") {
                        mode = 5;
                        strm >> bn >> sx >> sy;
                    }

                }
            }
        }
        //cout << "%\tRead file " << infname << " got " << TipInboardUSIDS.front().size() << " points" << endl;
    }





    //    please work pretty pretty please
    return NGRPS;
}

/**************************************************************/
void UTIL::write2D(string varname, string fname, Array<Array<double> > &input, int m,
        int n) {
    size_t dims[2] = { size_t(m), size_t(n)};
    //double d[m * n];

    double *d = new double[m * n];

    mat_t *mat;
    matvar_t *matvar;
    int k = 0;

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; i++) {
            d[k] = input[i][j];
            k++;
        }

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2,
                 dims, d, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }

    delete[] d;
}

/**************************************************************/
void UTIL::write2D(string varname, string fname, Array<Array<int> > &input, int m,
        int n) {
    size_t dims[2] = { size_t(m), size_t(n)};
    //double d[m * n];

    int *d = new int[m * n];

    mat_t *mat;
    matvar_t *matvar;
    int k = 0;

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; i++) {
            d[k] = input[i][j];
            k++;
        }

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_INT32, MAT_T_INT32, 2,
                dims, d, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }

    delete[] d;
}


void UTIL::WriteMATLABString(string vname, string fname, string data) {
    int m = (int) strlen(data.c_str());
    if (m == 0)
        return;
    write1D(vname, fname, data, m);
}


void UTIL::write1D(string varname, string fname, string &input, int m) {
    
    m+=1;
    
    size_t  dims[2];
    
    char *str = new char[strlen(input.c_str())+1];

    mat_t *mat;
    matvar_t *matvar;
    strcpy(str,input.c_str());

    dims[0] = 1;
    dims[1] = (int) strlen(str);
    

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_CHAR, MAT_T_INT8, 2, dims, str, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] str;
}

void UTIL::write1D(string varname, string fname, Array<double> &input, int m) {
    size_t dims[2] = {size_t(m), 1};
    //double d[m];
    double *d = new double[m];

    mat_t *mat;
    matvar_t *matvar;

    for (int i = 0; i < m; i++)
        d[i] = input[i];

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2,
                dims, d, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] d;
}


void UTIL::write1D(string varname, string fname, Array<int> &input, int m) {
    size_t dims[2] = {size_t(m), 1};
    //double d[m];
    int *d = new int[m];

    mat_t *mat;
    matvar_t *matvar;

    for (int i = 0; i < m; i++)
        d[i] = input[i];

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_INT32, MAT_T_INT32, 2,
                dims, d, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] d;
}

void UTIL::write1D(string varname, string fname, Array<long unsigned int> &input, int m) {
    size_t dims[2] = {size_t(m), 1};
    //double d[m];
    long unsigned int *d = new long unsigned int [m];

    mat_t *mat;
    matvar_t *matvar;

    for (int i = 0; i < m; i++)
        d[i] = input[i];

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_UINT64, MAT_T_UINT64, 2,
                dims, d, 0);
        Mat_VarWrite(mat, matvar, (matio_compression) 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] d;
}

int UTIL::readmat(string fname, string varname, Array <REAL> &data, Array <int> &dims, bool verbose) {
    data.clear();
    dims.clear();
    int err = 0;
    mat_t *mat;
    matvar_t *matvar;
    char *varname_cstr = new char [varname.size() + 1];
    strcpy(varname_cstr, varname.c_str());
    mat = Mat_Open(fname.c_str(), MAT_ACC_RDONLY);
    if (mat) {
        matvar = Mat_VarRead(mat, varname_cstr);
        if (matvar == NULL) {
            err = 1;
        } else {
            int numel = int(matvar->dims[0]);
            for (int i = 0; i < matvar->rank; ++i) {
                if (i > 0)
                    numel *= int(matvar->dims[i]);
                dims.push_back(int(matvar->dims[i]));
            }
            
            data.allocate(numel);
            for (int i = 0; i < numel; ++i)
                data[i] = ((double*)matvar->data)[i];
            
            if (verbose) {
                Mat_VarPrint(matvar, 1);
                cout << "Name of data: " << matvar->name << endl;
                cout << "Rank of input: " << matvar->rank << endl;
                cout << "Dimensions [";
                for (int i = 0; i < matvar->rank; ++i)
                    cout << matvar->dims[i] << " ";
                cout << "]" << endl;
                cout << "Number of elements: " << numel << endl;
            }
        }
        Mat_VarFree(matvar);
        Mat_Close(mat);
    } else {
        err = 1;
    }
    return err;
}



void UTIL::WriteMATLABMatrix1DVect3(string vname, string fname,
        Array<Vect3> &data) {
    int m = data.size();
    if (m == 0)
        return;

    Array < REAL > tmp(m);

    for (int i = 0; i < m; ++i)
        tmp[i] = data[i].x;

    UTIL::WriteMATLABMatrix1D(vname + "_x", fname, tmp);

    for (int i = 0; i < m; ++i)
        tmp[i] = data[i].y;

    UTIL::WriteMATLABMatrix1D(vname + "_y", fname, tmp);

    for (int i = 0; i < m; ++i)
        tmp[i] = data[i].z;

    UTIL::WriteMATLABMatrix1D(vname + "_z", fname, tmp);
}

void UTIL::WriteMATLABMatrix1DVect3(string vname, string fname,
        Vect3 &data) {

    UTIL::WriteMATLABMatrix1D(vname + "_x", fname, data.x);
    UTIL::WriteMATLABMatrix1D(vname + "_y", fname, data.y);
    UTIL::WriteMATLABMatrix1D(vname + "_z", fname, data.z);
}

void UTIL::WriteMATLABMatrix2DVect3(string vname, string fname,
        Array<Array<Vect3> > &data) {
    int m = data.size();
    if (m == 0)
        return;
    int n = data[0].size();
    if (n == 0)
        return;
    Array < Array < REAL > > tmp(m, n);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            tmp[i][j] = data[i][j].x;

    UTIL::write2D(vname + "_x", fname, tmp, m, n);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            tmp[i][j] = data[i][j].y;

    UTIL::write2D(vname + "_y", fname, tmp, m, n);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            tmp[i][j] = data[i][j].z;

    UTIL::write2D(vname + "_z", fname, tmp, m, n);
}

void UTIL::WriteMATLABMatrix2D(string vname, string fname,
        Array<Array<int> > &data) {
    int m = data.size();
    if (m == 0)
        return;
    int n = data[0].size();
    if (n == 0)
        return;

    UTIL::write2D(vname, fname, data, m, n);
}

void UTIL::WriteMATLABMatrix2D(string vname, string fname,
        Array<Array<double> > &data) {
    int m = data.size();
    if (m == 0)
        return;
    int n = data[0].size();
    if (n == 0)
        return;

    UTIL::write2D(vname, fname, data, m, n);
}

void UTIL::WriteMATLABMatrix1D(string vname, string fname, double data) {
    Array <double> temp(1, data);
    UTIL::write1D(vname, fname, temp, 1);
}

void UTIL::WriteMATLABMatrix1D(string vname, string fname, Array<double> &data) {
    int m = data.size();
    if (m == 0)
        return;
    UTIL::write1D(vname, fname, data, m);
}

void UTIL::WriteMATLABMatrix1D(string vname, string fname, Array<int> &data) {
    int m = data.size();
    if (m == 0)
        return;
    UTIL::write1D(vname, fname, data, m);

}

void UTIL::WriteMATLABMatrix1D(string vname, string fname, Array<long unsigned int> &data) {
    int m = data.size();
    if (m == 0)
        return;
    UTIL::write1D(vname, fname, data, m);

}

void UTIL::ReadBinaryVect3(Array <Vect3> &tmp, string fname) {

    int length;
    double * buffer;

    ifstream infile;

    infile.open(fname.c_str(), ios::binary);

    if (infile.is_open()) {
        // get length of file:
        infile.seekg(0, ios::end);
        length = int(infile.tellg());
        infile.seekg(0, ios::beg);

        // allocate memory:
        buffer = new double [length / (sizeof (double))];

        // read data as a block:
        infile.read((char *) buffer, length);
        infile.close();
        int sz = int(length / (3 * sizeof (double)));
        cout << "Reading " << length << " bytes, equiv to " << length / (sizeof (double)) << " doubles or " << length / (3 * sizeof (double)) << " Vect3 " << endl;

        tmp.allocate(sz);

        for (int i = 0; i < sz; ++i)
            tmp[i] = Vect3(buffer[i], buffer[i + sz], buffer[i + 2 * sz]);


        delete[] buffer;
    } else {
        cout << "Error opening file " << fname << endl;
    }
}
/**************************************************************/
Array <REAL> UTIL::globalLinspace(REAL start, REAL end, int n) {
    Array <REAL> output(n);

    REAL d = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i)
        output[i] = start + i * d;

    return output;

}

/**************************************************************/
Array <Vect3> UTIL::globalLinspace(Vect3 start, Vect3 end, int n) {
    Array <Vect3> output(n);

    Vect3 d = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i)
        output[i] = start + i * d;

    return output;

}

/**************************************************************/
double UTIL::interp2(Array<Array<double> > &X, Array<Array<double> > &Y, Array<Array<
        double> > &Z, double Xi, double Yi) {
    //	Array <double> XatYi(Y.size(), 0.0), ZatYi(Y.size(),0.0);
    //	for (int i = 0; i < Y.size(); ++i)
    //	{
    //		XatYi[i] = X[i][0];
    //		ZatYi[i] = interp1(Y[i],Z[i],Yi);
    //	}
    //
    //	double out = interp1(XatYi,ZatYi,Xi);

    //	Step 1:	Get upper and lower bounds for patch (assuming rectangular)
    int ni = X.size(), nj = X[0].size();
    double xl = X[0][0], xu = X[ni - 1][nj - 1], yl = Y[0][0], yu =
            Y[ni - 1][nj - 1];
    //	Step 2: If Xi or Yi are outwith these bounds return error
    //if (Xi >= xu || Xi <= xl || Yi >= yu || Yi <= yl)
    //	Do something

    //	Step 3: Convert (Xi, Yi) to (i,j) -> get neighbours and remainders
    double ix = ni * (Xi - xl) / (xu - xl) - 1, jy = nj * (Yi - yl) / (yu - yl)
            - 1;
    int ixp = int(ix) + 1, ixm = int(ix), jyp = int(jy) + 1, jym = int(jy);

    //	Step 4:	Check to see if point is on boundary, and if so use interior values instead
    //	Do something


    //cout << ix << " " << jy << " " << ixp << " " << jyp << " " << ixm << " " << jym << endl;
    //	Step 5:	calculate estimate
    double denom = (X[ixp][0] - X[ixm][0]) * (Y[0][jyp] - Y[0][jym]);
    double est = 0.0;
    est += Z[ixm][jym] * (X[ixp][0] - Xi) * (Y[0][jyp] - Yi) / denom;
    est += Z[ixp][jym] * (Xi - X[ixm][0]) * (Y[0][jyp] - Yi) / denom;
    est += Z[ixm][jyp] * (X[ixp][0] - Xi) * (Yi - Y[0][jym]) / denom;
    est += Z[ixp][jyp] * (Xi - X[ixm][0]) * (Yi - Y[0][jym]) / denom;

    return est;
}
/**************************************************************/

Vect3 UTIL::ODE4(Vect3 (*functocall)(Vect3, Vect3),Array <REAL> t, Array <Vect3> &Y, Vect3 ICs, Vect3 Params)
{
  Array <double> h = diff(t);
  Y.allocate(t.size());
  Y[0] = ICs;
  Vect3 F1, F2, F3, F4;
  for (int i = 1; i < t.size(); ++i)
  {
    double ti = t[i-1];
    double hi = h[i-1];
    Vect3 yi = Y[i-1];
    
    F1 = functocall(yi, Params);
    F2 = functocall(yi+0.5*hi*F1, Params);
    F3 = functocall(yi+0.5*hi*F2, Params);
    F4 = functocall(yi+hi*F3, Params);
   
    
    
    Y[i] = (yi + (hi/6.0)*(F1 + 2.0*F2 + 2.0*F3 + F4));
  }
  
  return Y.back();
}
/**************************************************************/

Vect3 UTIL::ODE4Final(Vect3 (*functocall)(Vect3, Vect3),REAL t0, REAL tmax, REAL dt, Vect3 ICs, Vect3 Params)
{
  Vect3 Y = ICs;
  Vect3 F1, F2, F3, F4;
  int n = ceil((tmax - t0)/dt) + 1;
  dt = (tmax - t0)/n;
  for (int i = 1; i < n; ++i)
  {
    F1 = functocall(Y, Params);
    F2 = functocall(Y+0.5*dt*F1, Params);
    F3 = functocall(Y+0.5*dt*F2, Params);
    F4 = functocall(Y+dt*F3, Params);
    Y += (dt/6.0)*(F1 + 2.0*F2 + 2.0*F3 + F4);
  }
  
  return Y;
}
