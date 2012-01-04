/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 35               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-16 00:1#$:  Date of last commit

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
REAL TIME_STEPPER::MaxTime = 0;
REAL TIME_STEPPER::SimTime = 0;
REAL TIME_STEPPER::SubStepTime = 0;
/**************************************************************/
TIME_STEPPER::TIME_STEPPER() {
    dt_prev = 1e16;
    dx = dy = dz = 1.0;
    ChangeOver = last_step = PruneNow = false;
    n = -1;
    RKStep = 0;
    t = substep_time = sim_time = 0.0;
    cfl_lim = 0.45;
    dt_out = 1;
    t_out = dt_out;
    lambda = mu = nu = 0.0;
    cpu_t = ticks();
    dump_next = show_rundata = false;
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
        dt = globalSystem->dtInit;
        cout << "### " << dt << " " << globalSystem->NumSubSteps << endl;
        long unsigned int t1 = ticks();
        BODY::BodySubStep(globalSystem->dtInit, globalSystem->NumSubSteps);
        long unsigned int t2 = ticks();
        //globalSystem->GetPressures(dt);
        //globalIO->write_m();
        //globalSystem->WriteBodies();
        
        Array < Vect3 > Posns, Vdir, VdirTree;
        Posns.push_back(Vect3(10,20,30));
        Posns.push_back(Vect3(20,-10,20));
        Posns.push_back(Vect3(-30,30,10));
        Posns.push_back(Vect3(40,40,0));
        Posns.push_back(Vect3(50,-10,-10));
        Posns.push_back(Vect3(60,20,-20));

        for (int i = 0; i < Posns.size(); ++i) {
            Vect3 Vel(0.0, 0.0, 0.0);
            for (int j = 0; j < BODY::VortexPositions.size(); ++j)
                Vel += globalDirectVel(Posns[i] - BODY::VortexPositions[j], BODY::VortexOmegas[j]);
            Vdir.push_back(Vel);

        }
        long unsigned int t3 = ticks();
        globalSystem->PutWakesInTree();
        long unsigned int t4 = ticks();
        globalOctree->Reset();
        globalOctree->InitVelsGetLaplacian();
        globalOctree->GetVels();
        

        for (int i = 0; i < Posns.size(); ++i) {
            Vect3 Vel(0.0, 0.0, 0.0);
            
            for (int j = 0; j < globalOctree->AllCells.size(); ++j)
                Vel += globalDirectVel(Posns[i] - globalOctree->AllCells[j]->Position, globalOctree->AllCells[j]->Omega);
            VdirTree.push_back(Vel);

        }
        

        
        


        
        for (int i = 0; i < Posns.size(); ++i)
            cout << "VD: " << Vdir[i] << ", VDTree: " << VdirTree[i] << " VFMM: " << globalOctree->TreeVel(Posns[i]) << endl;
        
        
        
        
        
        
        //        globalSystem->GetPanelFMMVelocities(0.0); //  t = t1
        first_step = false;
    } else {
        
        
        //      Bin panel wake into tree
        long unsigned int t1 = ticks();
        globalSystem->PutWakesInTree();
        long unsigned int t2 = ticks();
        cpu_sort_t = t2-t1;
        
        
        //      Calculate FMM
        long unsigned int t3 = ticks();
        if (!fmod((REAL) n, 100.0))
            PruneNow = true;
        globalOctree->Reset();
        long unsigned int t4 = ticks();
        globalOctree->InitVelsGetLaplacian();
        long unsigned int t5 = ticks();
        globalOctree->GetVels();


        //      Produce Output
        if (globalTimeStepper->dump_next) {
            
            if (TIME_STEPPER::SimTime < 1)
            for (int i = 0; i < globalOctree->AllCells.size(); ++i)
                if (((globalOctree->AllCells[i]->Position - BODY::Bodies[0]->CG).Mag() > (globalSystem->GambitScale * 2.0 * 5.03))) {
                    globalOctree->AllCells[i]->Omega = Vect3(0.0, 0.0, 0.0);
                    globalOctree->AllCells[i]->TransVars = Vect3(0.0, 0.0, 0.0);
                    globalOctree->AllCells[i]->HasLoad = false;
                }
            globalSystem->WriteData();
        }



        //      Calculate Panel contribution to FVM face fluxes
        long unsigned int t6 = ticks();
        globalSystem->GetFaceVels(); //  What do we do if this pushes it over the CFL limit?
        
        
        //      Get Timestep length
        long unsigned int t7 = ticks();
        time_step();
        
        //      Display Status
        long unsigned int t8 = ticks();
        globalIO->stat_step();
        
        
        //      Get FMM Vels on Panels @ t and t+dt
        long unsigned int t9 = ticks();
        globalSystem->GetPanelFMMVelocities(dt); //  t = t1
        
        //      Advance FVM to t + dt
        long unsigned int t10 = ticks();
        globalOctree->FVM(); //  t = t0
        long unsigned int t11 = ticks();
        BODY::BodySubStep(dt, globalSystem->NumSubSteps);
        long unsigned int t12 = ticks();
        cpu_ss_t = t12-t11;
        
        
        //      Evolve dw/dt
        globalOctree->Integrate(); //  t = t0 -> t1
        long unsigned int t13 = ticks();

        


        cout << "# Vortons                : " << BODY::VortexPositions.size() << endl;
        cout << "# faces                  : " << BODY::AllBodyFaces.size() << endl;

        
        cout << "Bin panel wake into tree : " << double(t2-t1)/1000.0 << endl;
        cout << "Calculate FMM: reset()   : " << double(t4-t3)/1000.0 << endl;
        cout << "InitVelsGetLaplacian     : " << double(t5-t4)/1000.0 << endl;
        cout << "GetVels()                : " << double(t6-t5)/1000.0 << endl;
        cout << "GetFaceVels()            : " << double(t7-t6)/1000.0 << endl;
        cout << "time_step()              : " << double(t8-t7)/1000.0 << endl;
        cout << "stat_step                : " << double(t9-t8)/1000.0 << endl;
        cout << "GetPanelFMMVelocities(.) : " << double(t10-t9)/1000.0 << endl;
        cout << "FVM                      : " << double(t11-t10)/1000.0 << endl;
        cout << "BodySubStep              : " << double(t12-t11)/1000.0 << endl;
        cout << "Integrate                : " << double(t13-t12)/1000.0 << endl;
        cout << "Total....................: " << double(t13-t1)/1000.0 << endl;
        
    }
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

    //  Calculate timestep length such that no body travels further than a single cell
    REAL OmRMax = 0.0;
    
    for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
        Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;
        // 	Get point kinematic velocity - rotational part first
        Vect3 Vrot = BODY::AllBodyFaces[i]->Owner->BodyRates.Cross(Pos);
        // 	Add to translational velocity....
        Vect3 Vkin = BODY::AllBodyFaces[i]->Owner->Velocity + Vrot;
    
        OmRMax = max(Vkin.Mag(),OmRMax);

    }
    
   
    
//    dt = min(dt_euler,cfl_lim/OmRMax);
    
   
    //  Check to see if this takes us over a time when we should be writing some output
    dump_next = false;

    if ((TIME_STEPPER::SimTime + dt >= t_out) || (TIME_STEPPER::SimTime + dt >= TIME_STEPPER::MaxTime)) {
        if (TIME_STEPPER::SimTime + dt >= TIME_STEPPER::MaxTime)
            last_step = true;
        dump_next = true;
        dt = dt_euler = min(t_out - TIME_STEPPER::SimTime, TIME_STEPPER::MaxTime - TIME_STEPPER::SimTime);
        t_out += dt_out;
    }


    //  If Lagrangian time-step is infinite (ie body is not moving) use a sensible number of sub-steps
    REAL dt_lagrange = dt_euler / 10 ;//min(dt_euler / 25, cfl_lim / (OmRMax));

    int nss = ceil(dt_euler / dt_lagrange);

    CFL = srad * dt;

    globalSystem->NumSubSteps = nss;
    

    //    if (n == 0) dump_next = true;

    if (n > 0) {
        dt_prev = dt;
        t += dt;
    }


    //    if ((ChangeOver == false) && (t >= 1))
    //    {
    //
    //        for (int i = 0; i < globalSystem->NumBodies; ++i)
    //            globalSystem->Bodies[i]->CG.vV = - globalSystem->Vinf;
    //
    //        globalSystem->Vinf.x = 0;
    //        globalSystem->Vinf.y = 0;
    //        globalSystem->Vinf.z = 0;
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

    for (int q = 0; q < globalSystem->NumTransVars; ++q) {
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
