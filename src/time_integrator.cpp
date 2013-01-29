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
int TIME_STEPPER::RKStep = 0;
bool TIME_STEPPER::RK2Mode = false;

/**************************************************************/
TIME_STEPPER::TIME_STEPPER() {
    dt_prev = 1e16;
    dx = dy = dz = 1.0;
    ChangeOver = last_step = PruneNow = false;
    n = -1;
    RKStep = 0;
    t = substep_time = sim_time = 0.0;
    cfl_lim = 0.45;
    dt_out = 1.0;
    t_out = dt_out;
    lambda = mu = nu = 0.0;
    cpu_t = ticks();
    dump_next = show_rundata = false;
    first_step = dump_next = true;
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

void TIME_STEPPER::DoFMM() {
    globalOctree->ResetAllVelsAndFields();
    globalOctree->GetVels();
}

/**************************************************************/
void TIME_STEPPER::TimeAdvance() {
#ifdef USE_SWSS
    //  This is a SWSS
    
    //  First 
    
    //  t0: Calc FMM and get DT
    DoFMM();
    //  t0: calculate face velocities due to body
    globalSystem->GetFaceVels();
    //  t0: calculate face velocities due to body
    time_step();
    //  t0: get panel FMM Vels
    if (globalSystem->useBodies) {
        globalSystem->GetPanelFMMVelocities(dt);

        //  t0: get time derivatives at t0

        //  t0: advance innter (body) timestep

        if (TIME_STEPPER::RK2Mode)
            BODY::BodySubStep(dt / 2.0, globalSystem->NumSubSteps);
        else
            BODY::BodySubStep(dt, globalSystem->NumSubSteps);

    }
    
    globalSystem->GetFaceVels();
    
    /* Order of integrations is as follows:
     * 0) Panel  
     * 1) stretch
     * 2) diffuse
     * 3) O2x
     * 4) O2y
     * 5) O2z
     * 
     * 
     * Other sweep is reverse
     */
    
    
    
    
    
    
    
    
    
    
#endif
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    TIME_STEPPER::RKStep = 0;
    //  t0: Calc FMM and get DT
    DoFMM();
    //  t0: calculate face velocities due to body
    globalSystem->GetFaceVels();
    time_step();
    //  t0: get panel FMM Vels
    if (globalSystem->useBodies) {
        globalSystem->GetPanelFMMVelocities(dt);

        //  t0: get time derivatives at t0

        //  t0: advance innter (body) timestep

        if (TIME_STEPPER::RK2Mode)
            BODY::BodySubStep(dt / 2.0, globalSystem->NumSubSteps);
        else
            BODY::BodySubStep(dt, globalSystem->NumSubSteps);

    }
    globalSystem->GetFaceVels();
    globalOctree->FVM(); //  dom_dt(t0)
    //  t0: advance outer wake to t*
    globalOctree->Integrate(); //  t = t0 -> t*

    if (TIME_STEPPER::RK2Mode) {
        TIME_STEPPER::RKStep = 1;
        //  t*: Do the FMM again
        DoFMM();
        if (globalSystem->useBodies) {
            globalSystem->GetPanelFMMVelocities(dt / 2);
            BODY::BodySubStep(dt / 2.0, globalSystem->NumSubSteps); //  t = t0 -> t*
            //  t*: calculate face velocities due to body
        }
        globalSystem->GetFaceVels();
        //  t*: calculate gradients at t*
        globalOctree->FVM(); //  dom_dt(t*)
        //  t*: advance to t1
        globalOctree->Integrate(); //  t = t0 -> t1
    }


    //  Put the wake in the tree from the bodytimestep
    if (globalSystem->useBodies) {
        globalSystem->PutWakesInTree();
    } else
        TIME_STEPPER::SimTime += dt;
    //  Clean up
    //if (!fmod((REAL) n, 100.0))
    //    PruneNow = true;
    globalOctree->Reset();




}

/**************************************************************/
void TIME_STEPPER::time_loop() {

    if (first_step) {
        dt = globalSystem->dtInit;
        if (globalSystem->useBodies) {
            BODY::BodySubStep(globalSystem->dtInit, globalSystem->NumSubSteps);
            globalSystem->PutWakesInTree();
        }
        globalOctree->Reset();
        
//        globalOctree->GetVels();
          for (int i = 0; i < globalOctree->AllCells.size(); ++i)
        globalOctree->AllCells[i]->Velocity = Vect3(0.0);
#pragma omp parallel for
        for (int i = 0; i < globalOctree->AllCells.size(); ++i)
            for (int j = 0; j < globalOctree->AllCells.size(); ++j)
                globalOctree->AllCells[i]->Velocity += globalDirectVel(globalOctree->AllCells[j]->Position - globalOctree->AllCells[i]->Position, globalOctree->AllCells[j]->Omega);
        first_step = false;
    } else {
#ifdef TIME_STEPS
        long unsigned int t1 = ticks();
#endif
        //        
        //        //      Bin panel wake into tree
        //        globalSystem->PutWakesInTree();
        //        //      Calculate FMM
        //        
        //        if (!fmod((REAL) n, 100.0))
        //            PruneNow = true;
        //        globalOctree->Reset();
        //        globalOctree->InitVelsGetLaplacian();
        //        globalOctree->GetVels();

        TimeAdvance();
        //      Produce Output
        if (globalTimeStepper->dump_next) {
            globalSystem->WriteData();
        }
        //        //      Calculate Panel contribution to FVM face fluxes
        //        globalSystem->GetFaceVels(); //  What do we do if this pushes it over the CFL limit?
        //        //      Get Timestep length
        //        time_step();
        //
        //        //      Get FMM Vels on Panels @ t and t+dt
        //        globalSystem->GetPanelFMMVelocities(dt); //  t = t1
        //        //      Advance FVM to t + dt
        //        globalOctree->FVM(); //  t = t0
        //        BODY::BodySubStep(dt, globalSystem->NumSubSteps);
        //        
        //        
        //        //      Evolve dw/dt
        //        globalOctree->Integrate(); //  t = t0 -> t1





#ifdef TIME_STEPS
        long unsigned int t13 = ticks();
        stringstream tmp;
        tmp << "Total....................: " << double(t13 - t1) / 1000.0 << endl;
        globalIO->step_data += tmp.str();
#endif
        //      Display Status
        globalIO->stat_step();
    }
}

/**************************************************************/
void TIME_STEPPER::time_step() {
#ifdef TIME_STEPS
    long unsigned int t7 = ticks();
#endif
    //  Need a sensible way to figure out how long to make the globla Eulerian time step
    //  and then the number and length of the Lagrangian time steps

    //  First off calculate the Eulerian time-step length
    srad = 0.0;

    globalOctree->GetSRad();

    REAL dt_euler = cfl_lim / srad.Mag(); //(srad.x + srad.y + srad.z);
    if (srad.Mag() == 0.0) dt_euler = globalSystem->dtInit;

    dt = dt_euler;

    //  Calculate timestep length such that no body travels further than a single cell
    REAL OmRMax = 0.0, MaxRadius = 0.0;
    //    REAL MaxX, MaxY, MaxZ; MaxX = MaxY = MaxZ = -1e32;
    //    REAL MinX, MinY, MinZ; MinX = MinY = MinZ = 1e32;
    if (globalSystem->useBodies) {
        for (int i = 0; i < BODY::AllBodyFaces.size(); ++i) {
            Vect3 Pos = BODY::AllBodyFaces[i]->CollocationPoint - BODY::AllBodyFaces[i]->Owner->CG;
            // 	Get point kinematic velocity - rotational part first
            Vect3 Vrot = BODY::AllBodyFaces[i]->Owner->BodyRates.Cross(Pos);
            // 	Add to translational velocity....
            Vect3 Vkin = BODY::AllBodyFaces[i]->Owner->Velocity + Vrot;

            OmRMax = max(Vkin.Mag(), OmRMax);


            MaxRadius = max(MaxRadius, Pos.Mag());

            //
            //        MaxX = max(MaxX, BODY::AllBodyFaces[i]->CollocationPoint.x);
            //        MinX = min(MinX, BODY::AllBodyFaces[i]->CollocationPoint.x);
            //        MaxY = max(MaxY, BODY::AllBodyFaces[i]->CollocationPoint.y);
            //        MinY = min(MinY, BODY::AllBodyFaces[i]->CollocationPoint.y);
            //        MaxZ = max(MaxZ, BODY::AllBodyFaces[i]->CollocationPoint.z);
            //        MinZ = min(MinZ, BODY::AllBodyFaces[i]->CollocationPoint.z);



        }
    }


    //    int DX = ceil(MaxX) - floor(MinX);
    //    int DY = ceil(MaxY) - floor(MinY);
    //    int DZ = ceil(MaxZ) - floor(MinZ);
    //    cout << "Maximum Radius " << MaxRadius << " Enclosing Box size: " << DX;
    //    cout << " " << DY << " " << DZ << " " << DX*DY*DZ << endl;
    //    cout << "--------------- Calcing Inerp Vels --------------" << endl;
    //    cout << globalOctree->AllCells.size()*BODY::AllBodyFaces.size() << " Interactions " << endl;
    //    Array < Array < Array < Vect3 > > > Posns, Vels;
    //            
    //    Posns = Array < Array < Array <Vect3> > > (DX, Array < Array < Vect3 > > (DY, Array < Vect3 > (DZ, Vect3(0.0))));
    //    Vels = Posns;
    //    
    //    for (int i = 0; i < DX; ++i)
    //        for (int j = 0; j < DY; ++j)
    //            for (int k = 0; k < DZ; ++k)
    //            {
    //                Vect3 XP(MinX + i, MinY + j, MinZ + k);
    //                Posns[i][j][k] = XP;
    //                for (int l = 0; l < globalOctree->AllCells.size(); ++l)
    //                    Vels[i][j][k] += UTIL::globalDirectVel(globalOctree->AllCells[l]->Position - XP, globalOctree->AllCells[l]->Omega, globalSystem->Del2);
    //                    
    //            }


    //    dt = min(dt_euler,cfl_lim/OmRMax);


    //  Check to see if this takes us over a time when we should be writing some output
    dump_next = false;

    if ((TIME_STEPPER::SimTime + dt >= t_out) || (TIME_STEPPER::SimTime + dt >= TIME_STEPPER::MaxTime)) {
        if (TIME_STEPPER::SimTime + dt >= TIME_STEPPER::MaxTime)
            last_step = true;
        dump_next = true;
        //dt = dt_euler = min(t_out - TIME_STEPPER::SimTime, TIME_STEPPER::MaxTime - TIME_STEPPER::SimTime);
        t_out += dt_out;
    }

    if (globalSystem->useBodies) {
        //  If Lagrangian time-step is infinite (ie body is not moving) use a sensible number of sub-steps
        REAL dt_lagrange = min(dt_euler / 10, cfl_lim / (OmRMax)); //

        int nss = ceil(dt_euler / dt_lagrange);


        globalSystem->NumSubSteps = nss;
    }
    CFL = srad * dt;

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
#ifdef TIME_STEPS
    long unsigned int t8 = ticks();
    stringstream tmp;
    tmp << "time_step()              : " << double(t8 - t7) / 1000.0 << endl;
    globalIO->step_data += tmp.str();
#endif
}

/**************************************************************/
void TIME_STEPPER::Integrate(FVMCell * cell) {

    Euler(cell);
    //    RK2(cell);

}

/**************************************************************/
void TIME_STEPPER::Euler(FVMCell * cell) {

    for (int q = 0; q < globalSystem->NumTransVars; ++q) {
        cell->TransVars[q] += dt * cell->TransDerivs[TIME_STEPPER::RKStep][q];
        cell->TransDerivs[TIME_STEPPER::RKStep][q] = Vect3(0.);
    }
}

/**************************************************************/
void TIME_STEPPER::RK2(FVMCell * cell) {
    //  Heun's Method for RK2
    if (RKStep == 0) { //  Predictor - an Euler step
        for (int q = 0; q < globalSystem->NumTransVars; ++q) {
            cell->TransVars[q] += dt * cell->TransDerivs[TIME_STEPPER::RKStep][q];
        }
    } else { //  Corrector - trapezoidal step
        for (int q = 0; q < globalSystem->NumTransVars; ++q) {
            cell->TransVars[q] += 0.5 * dt * (cell->TransDerivs[1][q] - cell->TransDerivs[0][q]);
            cell->TransDerivs[0][q] = Vect3(0.);
            cell->TransDerivs[1][q] = Vect3(0.);
        }
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
