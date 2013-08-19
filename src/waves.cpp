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
WaveField WaveField::Cnoidal;
WaveField WaveField::Linear;
REAL WaveField::Period, WaveField::Height, WaveField::WaveLength, WaveField::WaveNumber, WaveField::Celerity;
/**************************************************************/
WaveField::~WaveField(){
   
}
/**************************************************************/
WaveField::WaveField(){
    PeakFreqSet = ModalFreqSet = HSigSet = WindDataSet = false;
}
/**************************************************************/
void WaveField::getWaveFeld(void (WaveField::*spectrum_func)(), REAL T_min, REAL T_max, int n_omega) {
    Freqs = globalLinspace(1./(2*pi*T_max), 1./(2*pi*T_min), n_omega);
    WaveNumbers = Array <REAL> (n_omega,0.0);
    WaveLengths = Array <REAL> (n_omega,0.0);
    Amplitudes = Array <REAL> (n_omega,0.0);
    S = Array <REAL> (n_omega,0.0);
    ((this)->*(spectrum_func))();
    for (int i = 0; i < Freqs.size(); ++i) {
        WaveNumbers[i] = WaveField::Depth;
        REAL eps = 1;
        int count = 0;
        
        //      Newton-Raphson for dispersion relation
        while (eps > 1e-6) {
            REAL val = WaveNumbers[i] * tanh(Depth * WaveNumbers[i]) - Freqs[i] * Freqs[i] / SYSTEM::g;
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
        S[i] = alpha*SYSTEM::g/pow(Freqs[i],5) *  exp(-beta * pow(PeakFreq, 4) / pow(Freqs[i], 4))*pow(gamma,a);
        Amplitudes[i] = sqrt(2. * S[i] * dFreqs);
    }
    
    
}

void WaveField::WaveFieldCnoidal() {


    REAL T = WaveField::Period;
    REAL H = WaveField::Height;
    REAL d = WaveField::Depth;
    REAL omega = two_pi/T;
    
    //  Initial guess for wavenumber based on linear wave theory
    REAL k = d;
    {
        REAL om2_g = omega * omega / SYSTEM::g;
        int count = 0;
        REAL eps = 1.0;
        while (eps > 1e-12) {
            REAL val = k * tanh(d * k) - om2_g;
            REAL deriv = d * k * sech(d * k) * sech(d * k) + tanh(d * k);

            REAL new_val = k - (val / deriv);
            eps = fabs(k - new_val);
            k = abs(new_val);
            count += 1;
        }
        lambda = 2 * pi / k;
        cout << "Initial guess" << endl;
    cout << "k: "  << k << " lambda: " << lambda << " # iterations: " << count << endl;
    }

    

    REAL knew = k;


    REAL eps_k = 10;
    int numits = 0;

    Array <REAL> c(17), s(17);
    
    while (eps_k > 1e-12) {
        k = 0.5 * (k + knew);
        c[0] = 1.0;
        s[0] = 1.0;
        for (int i = 1; i < 17; ++i) {
            c[i] = c[i - 1] * cosh(k * d);
            s[i] = s[i - 1] * sinh(k * d);
            
        }

        B33 = 3 * (8 * c[6] + 1) / (64 * s[6]);
        B35 = (88128 * c[14] - 208224 * c[12] + 70848 * c[10] + 54000 * c[8] - 21816 * c[6] + 6264 * c[4] - 54 * c[2]) / ((12288 * s[12])*(6 * c[2] - 1));
        B55 = (192000 * c[16] - 262720 * c[14] + 83680 * c[12] + 20160 * c[10] - 7280 * c[8] + 7260 * c[6] - 1800 * c[4] - 1050 * c[2] + 225) / 
                ((12288 * s[10])*(6 * c[2] - 1)*(8 * c[4] - 11 * c[2] + 3));

        C1 = (8 * c[4] - 8 * c[2] + 9) / (8 * s[4]);
        C2 = (3840 * c[12] - 4096 * c[10] + 2592 * c[8] - 1008 * c[6] + 5944 * c[4] - 1830 * c[2] + 147) / ((512 * s[10])*(6 * c[2] - 1));

        
        // Now iterate internally for lambda
        REAL eps_lambda = 10;

        while (eps_lambda > 1e-12) {


            REAL f = 0.5 * k * H - (lambda + lambda * lambda * lambda * B33 + lambda * lambda * lambda * lambda * lambda * (B35 + B55));
            REAL f_dash = -1 - 3 * B33 * lambda * lambda - 5 * (B35 + B55) * lambda * lambda * lambda*lambda;

            REAL lambdanew = lambda - (f / f_dash);
            REAL lambdaold = lambda;
            lambda = 0.5 * (lambda + lambdanew);

            eps_lambda = abs((lambda - lambdaold) / lambda);

        }
        
        lambda2 = lambda * lambda;
        lambda3 = lambda2 * lambda;
        lambda4 = lambda3 * lambda;
        lambda5 = lambda4 * lambda;
        REAL F = omega * omega - SYSTEM::g * k * tanh(d * k)*(1 + lambda * lambda * C1 + lambda * lambda * lambda * lambda * C2);
        REAL F_dash = -SYSTEM::g * d * k * (1 + C1 * lambda * lambda + C2 * lambda * lambda * lambda * lambda) * sech(d * k) * sech(d * k) - SYSTEM::g * (1 + C1 * lambda * lambda + C2 * lambda * lambda * lambda * lambda) * tanh(d * k);
        knew = k - (F / F_dash);
        REAL kold = k;
        eps_k = abs((knew - kold) / kold);
        numits = numits + 1;
    }

    cout << "Final converged solution" << endl;
    cout << "k: "  << k  << " lambda: " << lambda << " # iterations: " << numits << endl;





    
    
    
    WaveField::WaveNumber = k;
    WaveField::WaveLength = 2 * pi / WaveField::WaveNumber;
    REAL L = WaveField::WaveLength = 2 * pi / k;
    WaveField::Celerity = 2*pi/k/T;
    
    
    beta = 2 * pi / L;
    

    A11 = 1 / s[1];
    A13 = -c[2]*(5 * c[2] + 1) / (8 * s[5]);
    A15 = -(1184 * c[10] - 1440 * c[8] - 1992 * c[6] + 2641 * c[4] - 249 * c[2] + 18) / (1536 * s[11]);

    A22 = 3 / (8 * s[4]);
    A24 = (192 * c[8] - 424 * c[6] - 312 * c[4] + 480 * c[2] - 17) / (768 * s[10]);

    A33 = (13 - 4 * c[2]) / (64 * s[7]);
    A35 = (512 * c[12] + 4224 * c[10] - 6800 * c[8] - 12808 * c[6] + 16704 * c[4] - 3154 * c[2] + 107) / (4096 * s[13]*(6 * c[2] - 1));
    A44 = (80 * c[6] - 816 * c[4] + 1338 * c[2] - 197) / (1536 * s[10]*(6 * c[2] - 1));
    A55 = -(2880 * c[10] - 72480 * c[8] + 324000 * c[6] - 432000 * c[4] + 163470 * c[2] - 16245) / (61440 * s[11]*(6 * c[2] - 1)*(8 * c[4] - 11 * c[2] + 3));




    B22 = (2 * c[2] + 1) * c[1] / (4 * s[3]);
    B24 = c[1]*(272 * c[8] - 504 * c[6] - 192 * c[4] + 322 * c[2] + 21) / (384 * s[9]);
    B33 = 3 * (8 * c[6] + 1) / (64 * s[6]);
    B35 = (88128 * c[14] - 208224 * c[12] + 70848 * c[10] + 54000 * c[8] - 21816 * c[6] + 6264 * c[4] - 54 * c[2]) / ((12288 * s[12])*(6 * c[2] - 1));
    B44 = c[1]*(768 * c[10] - 448 * c[8] - 48 * c[6] + 48 * c[4] + 106 * c[2] - 21) / (384 * s[9]*(6 * c[2] - 1));
    B55 = (192000 * c[16] - 262720 * c[14] + 83680 * c[12] + 20160 * c[10] - 7280 * c[8] + 7260 * c[6] - 1800 * c[4] - 1050 * c[2] + 225)/
            ((12288 * s[10])*(6 * c[2] - 1)*(8 * c[4] - 11 * c[2] + 3));


    C02 = SYSTEM::g * (tanh(beta * d));
    C1 = (8 * c[4] - 8 * c[2] + 9) / (8 * s[4]);
    C2 = (3840 * c[12] - 4096 * c[10] + 2592 * c[8] - 1008 * c[6] + 5944 * c[4] - 1830 * c[2] + 147) / ((512 * s[10])*(6 * c[2] - 1));
    C3 = -1 / (4 * s[1] * c[1]);
    C4 = (12 * c[8] + 36 * c[6] - 162 * c[4] + 141 * c[2] - 27) / (192 * c[1] * s[9]);

    cout << "------------------------ Constants ------------------------" << endl;
    cout << "B22 = " << B22 << "    B33 = " << B33 << "    B44 = " << B44 << endl;
    cout << "B24 = " << B24 << "    B35 = " << B35 << "    B55 = " << B55 << endl;
    cout << "A11 = " << A11 << "    A22 = " << A22 << "    A35 = " << A35 << endl;
    cout << "A13 = " << A13 << "    A24 = " << A24 << "    A44 = " << A44 << endl;
    cout << "A15 = " << A15 << "    A33 = " << A33 << "    A55 = " << A55 << endl;
    cout << "C1 = " << C1 << "    C3 = " << C3 << "    sinh(beta d) = " << s[1] << endl;
    cout << "C2 = " << C2 << "    C4 = " << C4 << "    cosh(beta d) = " << c[1] << endl;
    cout << "------------------------ Results ------------------------" << endl;
    cout << "Wave length:   L       = " << WaveField::WaveLength << " m" << endl;
    cout << "Wave celerity: C       = " << WaveField::Celerity << " m/s" << endl;

}



void WaveField::WaveFieldLinear() {


    REAL T = WaveField::Period;
    REAL d = WaveField::Depth;
    REAL omega = two_pi/T;
    
    //  Initial guess for wavenumber based on linear wave theory
    REAL k = d;
    {
        REAL om2_g = omega * omega / SYSTEM::g;
        int count = 0;
        REAL eps = 1.0;
        while (eps > 1e-12) {
            REAL val = k * tanh(d * k) - om2_g;
            REAL deriv = d * k * sech(d * k) * sech(d * k) + tanh(d * k);

            REAL new_val = k - (val / deriv);
            eps = fabs(k - new_val);
            k = abs(new_val);
            count += 1;
        }
        lambda = 2 * pi / k;
        cout << "Converged solution" << endl;
    cout << "k: "  << k << " lambda: " << lambda << " # iterations: " << count << endl;
    }

    
    
    WaveField::WaveNumber = k;
    WaveField::WaveLength = 2 * pi / WaveField::WaveNumber;
    WaveField::Celerity = 2*pi/k/T;
    
    
    
    cout << "------------------------ Results ------------------------" << endl;
    cout << "Wave length:   L       = " << WaveField::WaveLength << " m" << endl;
    cout << "Wave celerity: C       = " << WaveField::Celerity << " m/s" << endl;

}

Vect3 WaveField::LinearVelocity(Vect3 X, REAL t) {
    REAL a = WaveField::Height, om = two_pi / WaveField::Period, d = WaveField::Depth, g = 9.80665, k = WaveField::WaveNumber;
    REAL UWave = -(a * g * k * sin(X.x * k - om * t) * cosh(k * ((X.z - 0.8) + d))) / (om * cosh(d * k));
    REAL WWave = (a * g * k * cos(X.x * k - om * t) * sinh(k * ((X.z - 0.8) + d))) / (om * cosh(d * k));
    return Vect3(UWave, 0.0, WWave);
}

REAL WaveField::LinearPerturbationPotential(Vect3 X, REAL t) {
    REAL a = WaveField::Height, om = two_pi / WaveField::Period, d = WaveField::Depth, g = 9.80665, k = WaveField::WaveNumber;
    return (a * g / om) * cosh(k * (d + (X.z - 0.8))) * cos(om * t - k * X.x) / cosh(k * d);
}

REAL WaveField::CnoidalPerturbationPotential(Vect3 X, REAL time) {
    REAL px = X.x;
    REAL pz = X.z;

    REAL C = WaveField::Celerity;
    REAL d = WaveField::Depth;
    REAL S = pz + d;
    REAL THETAs = beta * (px - C * time);
    REAL L1 = (lambda * A11 + lambda3 * A13 + lambda5 * A15) * cosh(beta * S) * sin(THETAs);
    REAL L2 = (lambda2 * A22 + lambda4 * A24) * cosh(2 * beta * S) * sin(2 * THETAs);
    REAL L3 = (lambda3 * A33 + lambda5 * A35) * cosh(3 * beta * S) * sin(3 * THETAs);
    REAL L4 = (lambda4 * A44) * cosh(4 * beta * S) * sin(4 * THETAs);
    REAL L5 = (lambda5 * A55) * cosh(5 * beta * S) * sin(5 * THETAs);

    REAL betaphi_C = L1 + L2 + L3 + L4 + L5;
    return betaphi_C * C / beta;
}

Vect3 WaveField::CnoidalVelocity(Vect3 X, REAL time)
{
    REAL px = X.x;
    REAL pz = X.z;
    
    REAL C = WaveField::Celerity;
    REAL d = WaveField::Depth;
    
    REAL pu = (C*(beta*(A11*lambda + A13*lambda3 + A15*lambda5)*
        cos(beta*(-(C*time) + px))*cosh(beta*(d + pz)) + 
        2*beta*(A22*lambda2 + A24*lambda4)*cos(2*beta*(-(C*time) + px))*
        cosh(2*beta*(d + pz)) + 
        3*beta*(A33*lambda3 + A35*lambda5)*cos(3*beta*(-(C*time) + px))*
        cosh(3*beta*(d + pz)) + 
        4*A44*beta*lambda4*cos(4*beta*(-(C*time) + px))*
        cosh(4*beta*(d + pz)) + 
        5*A55*beta*lambda5*cos(5*beta*(-(C*time) + px))*
        cosh(5*beta*(d + pz))))/beta;
    
    REAL pw = (C*(beta*(A11*lambda + A13*lambda3 + A15*lambda5)*
        sin(beta*(-(C*time) + px))*sinh(beta*(d + pz)) + 
        2*beta*(A22*lambda2 + A24*lambda4)*sin(2*beta*(-(C*time) + px))*
        sinh(2*beta*(d + pz)) + 
        3*beta*(A33*lambda3 + A35*lambda5)*sin(3*beta*(-(C*time) + px))*
        sinh(3*beta*(d + pz)) + 
        4*A44*beta*lambda4*sin(4*beta*(-(C*time) + px))*
        sinh(4*beta*(d + pz)) + 
        5*A55*beta*lambda5*sin(5*beta*(-(C*time) + px))*
        sinh(5*beta*(d + pz))))/beta;
    
    return Vect3(pu, 0.0, pw);
}
///**************************************************************/
//REAL WaveField::getPhi2D(Vect3 XP)
//{
//    
//    
//    
//    return 0.0;
//}