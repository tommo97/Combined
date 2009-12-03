/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 2                $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2009-10-28 20:1#$:  Date of last commit

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


#ifndef TIMESTEP_INCL
#define TIMESTEP_INCL

#include "includes.hpp"
#include "types.hpp"
#include "array.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "io.hpp"



class TIME_STEPPER
{
    public:
        TIME_STEPPER();
        int n;
        REAL substep_time, t, t_out, cfl_lim, dt_out, dt_prev, sim_time;
        REAL max_t, dt, dx, dy, dz;
        REAL *A;
        REAL lambda, mu, nu;
        bool dump_next, first_step, ChangeOver;
        Vect3 centre;
        int num_images, RKStep;
        long int cpu_t;
        Vect3 srad, CFL;

#ifdef TIMER
        long int cpu_t_fmm_e;
        long int cpu_t_fvm_e;
        long int cpu_t_octree_e;
        long int cpu_t_IO_e;
        long int cpu_t_fmm_s;
        long int cpu_t_fvm_s;
        long int cpu_t_octree_s;
        long int cpu_t_IO_s;
        long int cpu_tot_s;
#endif
        void time_loop();
        void evolution(int);
        void reconstruction(int);
        void time_step();
        void PanelOctreeInterface();
        void boundary();
        void boundary_init(int);
        void Euler(FVMCell *);
        void RK2(FVMCell *);
        void RK4(FVMCell *);
        void ABM4(FVMCell *);
        void Integrate(FVMCell *);

	//~TIME_STEPPER();
};

long int ticks();
#endif /* if !defined(TIME_STEPPER) */
