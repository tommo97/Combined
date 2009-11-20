/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 35               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2009-11-16 00:1#$:  Date of last commit

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




/**************************************************************/
#include "time_integrator.hpp"
#include <gsl/gsl_sf_bessel.h>
/**************************************************************/
TIME_STEPPER::TIME_STEPPER() {
    dt_prev = 1e16;
    dx = dy = dz = 1.0;
    ChangeOver = false;
    n = -1;
    RKStep = 0;
    t = substep_time = 0.0;
    cfl_lim = 0.45;
    dt_out = .1;
    t_out = dt_out;
    max_t = 20;
    lambda = mu = nu = 0.0;
    cpu_t = ticks();
    dump_next = false;
    first_step = true;
}
/**************************************************************/
void trig_cone(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL magnitude, REAL radius);

void trig_cone(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL magnitude, REAL radius) {
    for (int x = 1; x <= NX; ++x) {
        for (int y = 1; y <= NY; ++y) {
            REAL r = sqrt((centre.x - x)*(centre.x - x)+ (centre.y - y)*(centre.y - y));
            if (r / radius < 1) {
                Omega.push_back(Vect3(0., 0., magnitude * sin((pi * 0.5 * (1 + (r / radius))))));
                X.push_back(Vect3(x, y, 0.));
            }

        }
    }
}
/**************************************************************/
void cylinder_dipole(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL offset, REAL amplitude, REAL radius);

void cylinder_dipole(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL offset, REAL amplitude, REAL radius) {
    REAL c_cone[3] = {centre.x + offset, centre.x - offset, centre.y};
    REAL r1 = 0, r2 = 0;

    for (int x = 1; x <= 300; ++x) {
        for (int y = 1; y <= 300; ++y) {
            r1 = sqrt((c_cone[0] - x)*(c_cone[0] - x) + (c_cone[2] - y)*(c_cone[2] - y));
            r2 = sqrt((c_cone[1] - x)*(c_cone[1] - x) + (c_cone[2] - y)*(c_cone[2] - y));
            if (r1 <= radius) {
                Omega.push_back(Vect3(0., 0., (0. * amplitude) - amplitude));
                X.push_back(Vect3(x, y, 0.));
            }
            if (r2 <= radius) {
                Omega.push_back(Vect3(0., 0., (0. * amplitude) + amplitude));
                X.push_back(Vect3(x, y, 0.));
            }
        }
    }
}
/**************************************************************/
void lamb_dipole(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL theta_dipole, REAL radius, REAL a, REAL k);

void lamb_dipole(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL theta_dipole, REAL radius, REAL a, REAL k) {
    //  Produce a Lamb Dipole, eg
    //        lamb_dipole(Vect3(100, 50, 0), X, OMEGA, 30, pi / 2, 30, 30, 3.81 / 30);


    REAL U;
    for (int x = 1; x <= NX; ++x) {
        for (int y = 1; y <= NY; ++y) {
            U = 0.;
            if (x < centre.x) U = -amplitude;
            if (x > centre.x) U = amplitude;

            REAL r = sqrt((centre[0] - x)*(centre[0] - x) + (centre[1] - y)*(centre[1] - y));
            if (r < radius) {
                REAL theta = atan2((centre[1] - y)*(centre[1] - y) + 1e-16, (centre[0] - x)*(centre[0] - x) + 1e-16);
                REAL numerator = gsl_sf_bessel_J1(k * r);
                REAL denominator = gsl_sf_bessel_J0(k * a);
                Omega.push_back(Vect3(0.0, 0.0, -U * 2 * k * sin(theta - theta_dipole) * numerator / denominator));
                X.push_back(Vect3((REAL) x, (REAL) y, 0.0));
            }
        }
    }
}
/**************************************************************/
void sheet(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL radius, REAL k);

void sheet(Vect3 centre, Array <Vect3> &X, Array <Vect3> &Omega, REAL amplitude, REAL radius, REAL k) {
    int nr = 1000, nt = 360;
    amplitude += k;
    REAL r = 0;
    for (int i = 0; i < nr; ++i) {
        REAL gamma = 2 * radius * sqrt(1 - ((r / radius) * (r / radius)));
        REAL theta = 0;
        for (int j = 0; j < nt; ++j) {
            X.push_back(centre + Vect3(r * cos(theta), r * sin(theta), 0.0));
            Omega.push_back(Vect3(-gamma * sin(theta), gamma * cos(theta), 0.0));
            theta += two_pi / nt;
        }
        r += radius / nr;
    }

    Array <Vect3> Xn, Omn;
    Xn.push_back(X[0]);
    Omn.push_back(Omega[0]);
    for (int i = 1; i < X.size(); ++i) {
        X[i] = floor(X[i]) + .5;
        bool put_in = true;
        for (int j = 0; j < Xn.size(); ++j)
            if (X[i].x == Xn[j].x && X[i].y == Xn[j].y && X[i].z == Xn[j].z) put_in = false;

        if (put_in) {
            Xn.push_back(floor(X[i]) + .5);
            Omn.push_back(Omega[i]);
        }

    }

    X = Xn;
    Omega = Omn;
}

/**************************************************************/

void TIME_STEPPER::time_loop() {

    if (first_step) {
    	globalSystem->NumSubSteps = 25;
    	globalSystem->dtInit = 5;
		for (int i = 0; i < globalSystem->NumBodies; ++i)
			globalSystem->Bodies[i]->InitNascentWake(globalSystem->dtInit / globalSystem->NumSubSteps);
		globalSystem->SetupGlobalInfluenceMatrices();
		dt = globalSystem->dtInit;
        globalSystem->BodySubStep(dt, globalSystem->NumSubSteps);
        globalSystem->GetPressures(dt);
        globalIO->write_m();
        globalSystem->WriteBodies();
//        globalSystem->PutWakesInTree();
//        globalOctree->Reset();
//        globalOctree->InitVelsGetLaplacian();
//        globalOctree->GetVels();
        first_step = false;
    }

//    time_step();
//    globalIO->stat_step();
//
//    globalOctree->FVM(); //  t = t0
//    globalOctree->Integrate(); //  t = t0 -> t1
//
//    globalSystem->BodySubStep(dt, globalSystem->NumSubSteps);
//
//    globalSystem->PutWakesInTree();
//    globalOctree->Reset();
//    globalOctree->InitVelsGetLaplacian();
//    globalOctree->GetVels();
//    globalSystem->GetPanelFMMVelocities();  //  t = t1
//    globalSystem->GetFaceVels(); //  What do we do if this pushes it over the CFL limit?
//
//    if (globalTimeStepper->dump_next){
////        globalSystem->WriteDomain();
//        globalSystem->WriteVorticity();
////        globalOctree->Reset();
//    }

}
/**************************************************************/
void TIME_STEPPER::time_step() {

    //  Need a sensible way to figure out how long to make the globla Eulerian time step
    //  and then the number and length of the Lagrangian time steps

    //  First off calculate the Eulerian time-step length
    srad = 0.0;

    globalOctree->GetSRad();

    REAL dt_euler = cfl_lim / srad.Mag(); //(srad.x + srad.y + srad.z);
    if (srad.Mag() == 0.0) dt_euler = globalSystem->dtInit;

    dt = dt_euler;

    dump_next = false;

    if ((t + dt >= t_out) || (t + dt >= max_t)) {
        dump_next = true;
        dt = dt_euler = min(t_out - t, max_t - t);
        t_out += dt_out;
    }


    //  Now calculate the Lagrangian time-step length
    REAL OmRMax = 0;
    for (int i = 0; i < globalSystem->NumBodies; ++i)
        OmRMax = max(OmRMax, max(fabs(globalSystem->Bodies[i]->CG.vV + globalSystem->Bodies[i]->BodyRates.Cross(globalSystem->Bodies[i]->Rmax))));

    //  If Lagrangian time-step is infinite (ie body is not moving) use a sensible number of sub-steps
    REAL dt_lagrange = min(dt_euler / 10, cfl_lim / OmRMax);

    int nss = ceil(dt_euler / dt_lagrange);
    
    CFL = srad * dt;

    globalSystem->NumSubSteps = nss;

//    if (n == 0) dump_next = true;

    if (n > 0)
    {
        dt_prev = dt;
        t = t + dt;
    }


//    if ((ChangeOver == false) && (t >= 1))
//    {
//
//        for (int i = 0; i < globalSystem->NumBodies; ++i)
//            globalSystem->Bodies[i]->CG.vV = - globalSystem->Vinf;
//
//        globalSystem->uinf = 0;
//        globalSystem->vinf = 0;
//        globalSystem->winf = 0;
//        globalSystem->Vinf = 0;
//        ChangeOver = true;
//    }


    n++;
}

/**************************************************************/
void TIME_STEPPER::Integrate(FVMCell * cell) {

    Euler(cell);

}

/**************************************************************/
void TIME_STEPPER::Euler(FVMCell * cell) {

    for (int q = 0; q < globalSystem->NumTransVars; ++q)
    {
        cell->TransVars[q] += dt * cell->TransDerivs[q];
        cell->TransDerivs[q] = 0.;
    }
}

/**************************************************************/
void TIME_STEPPER::RK2(FVMCell * cell) {
    //  Heun's Method for RK2
    if (RKStep == 0) { //  Predictor - an Euler step
        cell->Omega += dt * cell->Deriv[0];
    } else { //  Corrector - trapezoidal step
        cell->Omega += .5 * dt * (cell->Deriv[1] - cell->Deriv[0]);
        cell->Deriv.clear();
    }
}

/**************************************************************/
void TIME_STEPPER::RK4(FVMCell * cell) {
    switch (RKStep) {
        case 0:
        {
            cell->Omega += dt * .5 * cell->Deriv[0];
            break;
        }

        case 1:
        {
            cell->Omega += .5 * dt * (cell->Deriv[1] - cell->Deriv[0]);
            break;
        }
        case 2:
        {
            cell->Omega += dt * (cell->Deriv[2] - .5 * cell->Deriv[1]);
            break;
        }
        case 3:
        {
            cell->Omega += dt * ((cell->Deriv[0] + 2 * (cell->Deriv[1] + cell->Deriv[2]) + cell->Deriv[3]) / 6 - cell->Deriv[2]);
            cell->Deriv.clear();
            break;
        }
    }
}

/**************************************************************/
void TIME_STEPPER::ABM4(FVMCell * cell) {
    if (cell->age < 3) {
        if (RKStep == 0)
            cell->Omega += dt * cell->Deriv[cell->age];
        else {
            cell->Omega += .5 * dt * (cell->Deriv[cell->age + 1] - cell->Deriv[cell->age]);
            cell->Deriv[cell->age] = .5 * (cell->Deriv[cell->age] + cell->Deriv[cell->age + 1]);
            cell->Deriv.pop_back();
            cell->age++;
        }
    } else {
        if (RKStep == 0) //  Predictor
        {
            cell->Omega += (dt / 24) * (55 * cell->Deriv[3] - 59 * cell->Deriv[2] + 37 * cell->Deriv[1] - 9 * cell->Deriv[0]);
        } else //  Corrector
        {
            cell->Omega += (dt / 24) * (9 * cell->Deriv[4] + 19 * cell->Deriv[3] - 5 * cell->Deriv[2] + cell->Deriv[1]);
            cell->Omega -= (dt / 24) * (55 * cell->Deriv[3] - 59 * cell->Deriv[2] + 37 * cell->Deriv[1] - 9 * cell->Deriv[0]);
            // Trim of oldest derivative and also derivative of the predictor


#ifndef USE_ARRAY
            Array <Vect3> temp;
            for (int i = 1; i < cell->Deriv.size() - 1; ++i) temp.push_back(cell->Deriv[i]);
            cell->Deriv = temp;
#else
            cell->Deriv.pop_front();
            cell->Deriv.pop_back();
#endif
            cell->age++;
        }
    }
}