/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2013
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form

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

#ifndef WAVES_HPP
#define	WAVES_HPP

#include "includes.hpp"
#include "types.hpp"
#include "io.hpp"
#include "tree.hpp"
#include "body.hpp"


class ModalFreqNotSet {
public:
    ModalFreqNotSet() {
        cout << "Modal Frequency Not Set ... Aborting." << endl;
    }
};

class HSigNotSet {
public:
    HSigNotSet() {
        cout << "Significant Wave Height Not Set ... Aborting." << endl;
    }
};

class PeakFreqNotSet {
public:
    PeakFreqNotSet() {
        cout << "Peak Frequency Not Set ... Aborting." << endl;
    }
};
    
class WaveField {
public:

    
    virtual ~WaveField();

    WaveField();
    
    void getWaveFeld(void (WaveField::*spectrum_func)(), REAL omega_min, REAL omega_max, int n_omega);
    
    void Bretschneider();
    
    void JONSWAP();
    
    void PiersonMoskowitz();
    
    REAL getPhi2D(Vect3);
    
    void SetPeakFreq(REAL in)
    {
        PeakFreq = in;
        PeakFreqSet = true;
    }
    
    void SetHSig(REAL in)
    {
        HSig = in;
        HSigSet = true;
    }
    
        void SetModalFreq(REAL PF)
    {
        ModalFreq = PF;
        ModalFreqSet = true;
    }
        
        Array <REAL> Spectrum()
        {
            return S;
        }

    static REAL Depth;
private:
    bool PeakFreqSet, ModalFreqSet, HSigSet, WindDataSet;
    static REAL Gravity;
    Array <REAL> S, Freqs, Amplitudes, Phases, WaveNumbers, WaveLengths, Periods;
    REAL HSig, ModalFreq, PeakFreq;

    
};
#endif	/* WAVES_HPP */
