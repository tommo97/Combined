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

#ifndef UNIT_TESTS_HPP
#define	UNIT_TESTS_HPP

#include "includes.hpp"
#include "types.hpp"
#include "io.hpp"
#include "tree.hpp"
#include "body.hpp"




class TEST {
public:
    static void TestFMM(int argc, char *argv[]);
    static void TestPanel();
    static void SolveMatfileVels(string fname, int pmax, REAL del2);
    static void TestBiotSavart();
    static void TestBEM();
    static void SimpleTestPanel();
    static void TestBEMFVM();
    static void TestBulkLoader(int n);
    

    
};

class RunningStat {
public:

    RunningStat() : m_n(0), mmax(-1e32), mmin(1e32) {
    }

    void Clear() {
        m_n = 0;
    }

    void Push(REAL x) {
        m_n++;

        if (x > mmax) mmax = x;
        if (x < mmin) mmin = x;
        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1) {
            m_oldM = m_newM = x;
            m_oldS = 0.0;
        } else {
            m_newM = m_oldM + (x - m_oldM) / m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

            // set up for next iteration
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
    }

    int NumDataValues() const {
        return m_n;
    }

    REAL Mean() const {
        return (m_n > 0) ? m_newM : 0.0;
    }

    REAL Variance() const {
        return ( (m_n > 1) ? m_newS / (m_n - 1) : 0.0);
    }

    REAL StandardDeviation() const {
        return sqrt(Variance());
    }

    REAL Max() const {
        return mmax;
    }

    REAL Min() const {
        return mmin;
    }
private:
    int m_n;
    REAL mmax, mmin;
    REAL m_oldM, m_newM, m_oldS, m_newS;
};

#endif	/* UNIT_TESTS_HPP */

