/*
 This file is part of the Combined Wake Modelling Code Version 1.0

 V3D Code Copyright Tom McCombes 2011
 This code solves the 3D unsteady incompressible
 Navier-Stokes equations in velociy vorticity form


 $Rev:: 29               $:  Revision of last commit
 $Author:: tom           $:  Author of last commit
 $Date:: 2011-11-09 21:5#$:  Date of last commit

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

unsigned long int LineVelCnt = 0;
int NumThreads = 1;

long unsigned int UTIL::cpu_t = 0;


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
        cout << "%\tRead file " << infname << " got " << Surfaces.size() << " surfaces of size " << Surfaces[0].size() << " by " << Surfaces[0][0].size() << endl;
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
    int dims[2] = {m, n};
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
        Mat_VarWrite(mat, matvar, 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }

    delete[] d;
}

/**************************************************************/
void UTIL::write2D(string varname, string fname, Array<Array<int> > &input, int m,
        int n) {
    int dims[2] = {m, n};
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
        Mat_VarWrite(mat, matvar, 0);
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
    int dims[2];
    
    char *str = new char[strlen(input.c_str())+1];

    mat_t *mat;
    matvar_t *matvar;
    strcpy(str,input.c_str());

    dims[0] = 1;
    dims[1] = (int) strlen(str);
    

    mat = Mat_Open(fname.c_str(), MAT_ACC_RDWR);
    if (mat) {
        matvar = Mat_VarCreate(varname.c_str(), MAT_C_CHAR, MAT_T_INT8, 2, dims, str, 0);
        Mat_VarWrite(mat, matvar, 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] str;
}

void UTIL::write1D(string varname, string fname, Array<double> &input, int m) {
    int dims[2] = {m, 1};
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
        Mat_VarWrite(mat, matvar, 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] d;
}

void UTIL::write1D(string varname, string fname, Array<int> &input, int m) {
    int dims[2] = {m, 1};
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
        Mat_VarWrite(mat, matvar, 0);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
    delete[] d;
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

