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


#include "system.hpp"
#include "types.hpp"
#include "waves.hpp"

REAL WaveField::Depth;
REAL WaveField::Gravity = 9.80665;
/**************************************************************/
WaveField::~WaveField(){
   
}
/**************************************************************/
WaveField::WaveField(){
    PeakFreqSet = ModalFreqSet = HSigSet = WindDataSet = false;
}
/**************************************************************/
void WaveField::getWaveFeld(void (WaveField::*spectrum_func)(), REAL T_min, REAL T_max, int n_omega) {
    Freqs = globalLinspace(1./(2.*pi*T_max), 1./(2.*pi*T_min), n_omega);
    WaveNumbers = Array <REAL> (n_omega,0.0);
    WaveLengths = Array <REAL> (n_omega,0.0);
    Amplitudes = Array <REAL> (n_omega,0.0);
    S = Array <REAL> (n_omega,0.0);
    ((this)->*(spectrum_func))();
    for (int i = 0; i < Freqs.size(); ++i) {
        WaveNumbers[i] = WaveField::Depth;
        REAL eps = 1;
        int count = 0;
        
        //      Newton-Raphson for dispersion relation...
        while (eps > 1e-6) {
            REAL val = WaveNumbers[i] * tanh(Depth * WaveNumbers[i]) - Freqs[i] * Freqs[i] / Gravity;
            REAL deriv = Depth * WaveNumbers[i] * sech(Depth * WaveNumbers[i]) * sech(Depth * WaveNumbers[i]) + tanh(Depth * WaveNumbers[i]);
            REAL new_val = WaveNumbers[i] - (val / deriv);
            eps = fabs(WaveNumbers[i] - new_val);
            WaveNumbers[i] = new_val;
            count++;

            if (count > 1e6)
                eps = 0;
        }
        WaveLengths[i] = 2.0*pi/WaveNumbers[i];
    }
    
}    /**************************************************************/
void WaveField::Bretschneider() {

    if (!ModalFreqSet)
        throw ModalFreqNotSet();

    if (!HSigSet)
        throw HSigNotSet();

    REAL dFreqs = Freqs[1] - Freqs[0];
    Amplitudes = Array <REAL > (Freqs.size(), 0.0);
    for (int i = 0; i < Freqs.size(); ++i) {
        S[i] = (5. / 16.) * (pow(ModalFreq, 4) / pow(Freqs[i], 5)) * HSig * HSig * exp(-5.0 * pow(ModalFreq, 4) / 4.0 / pow(Freqs[i], 4));
        Amplitudes[i] = sqrt(2. * S[i] * dFreqs);
    }
}

/**************************************************************/
void WaveField::JONSWAP() {

    if (!PeakFreqSet)
        throw PeakFreqNotSet();

    if (!HSigSet)
        throw HSigNotSet();
    
    REAL sigma;
    REAL alpha = 0.0081;
    REAL gamma = 3.3;
    REAL beta = 1.25;
    
    REAL dFreqs = Freqs[1] - Freqs[0];
    Amplitudes = Array <REAL > (Freqs.size(), 0.0);
    for (int i = 0; i < Freqs.size(); ++i) {
        if (Freqs[i]>PeakFreq)
            sigma = 0.09;
        else
            sigma = 0.07;
        
        REAL a = exp(-(Freqs[i] - PeakFreq)*(Freqs[i] - PeakFreq)/(2.0*PeakFreq*PeakFreq*sigma*sigma));
        S[i] = alpha*Gravity/pow(Freqs[i],5) *  exp(-beta * pow(PeakFreq, 4) / pow(Freqs[i], 4))*pow(gamma,a);
        Amplitudes[i] = sqrt(2. * S[i] * dFreqs);
    }
    
    
}
/**************************************************************/
REAL WaveField::getPhi2D(Vect3 XP)
{
    
    
    
    return 0.0;
}