/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2011
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 27               $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2011-11-07 22:3#$:  Date of last commit

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


#ifndef IO_HPP
#define IO_HPP
#include "includes.hpp"
#include "types.hpp"
#include "array.hpp"
#include "system.hpp"

#include "time_integrator.hpp"



class IO {
public:
    IO();

    int num_file;
    int num_dat;
    int num_forces;
    int num_pressures;
    int num_m;
    int num_msh;
    int num_images;
    int ID;
//    vector <REAL> x, y, z;
//    Array < Array <REAL> > domain;
//    Gnuplot *Plot;
    string pad, filename, directory, image_type, file_type, out_name, suffix, OS, Home;
    string top_header, top_data, top_command;
    string ProcessID, header, step_data, step_header;
    ofstream currentOfstream;
    Array <int> NumFiles;
    Array <string> FileNames;
    #ifdef _PNGWRITER
    pngwriter CurrentPNG;
    #endif
    list <string> lines;

    void read_neu(string infname, Array <Vect3> &X, Array <Array <int> > &PNLS, Array < Array <int> > &GROUPS, Array < Array < int > > &BCS, Array <string> &);

    void print_header();

    void print_line();

    void PrepOutputDir();

    void print_header(ostream&);

    void write_ISA();

    void WriteBinary();

    string StringToLower(string);

    string StringToUpper(string);

    void read_dat(char* infname, Array <Vect3> &x, Array <Vect3> &omega);

    void write_file(stringstream &, string outfname, string ext, bool);
    
    void write_2D_mat(Array < Array <REAL> > &outdata, string vname, string OutName, bool disp);

    void write_dat();

    void write_vels_and_pressures();

    void read_input(string infname);

    void write_m();

    void write_forces();

    void write_GMSH();

    template<class T>
    inline string to_string(const T& t) {
    	stringstream ss;
    	ss << t;
    	return ss.str();
    }



#ifdef _PNGWRITER
    void write_image();

    void create_image(string);
#endif
    void print_stat_step(ostream &);

    void display_stat_step();

    void stat_step();
#ifdef use_NCURSES
    void start_curses();

    void end_curses();
#endif      

    ~IO();

    enum exception {
        NO_FILE
    };

};


#endif
