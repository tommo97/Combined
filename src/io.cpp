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

#include "io.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
//#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communication with Gnuplot

IO::IO() {
    num_file = 0;
    num_dat = 0;
    num_m = 0;
    num_msh = 0;
    num_forces = 0;
    num_images = 0;
    //    Plot = NULL;
    //    domain.assign(NX,Array<REAL> (NY,0.));

#ifdef use_NCURSES
    start_curses();
#endif
    stringstream temp, out_stream;
    print_header(out_stream);

    temp << getpid();
    ProcessID = temp.str();
    OS = getenv("OSTYPE");
    Home = getenv("HOME");
    temp >> globalSystem->ProcessID;
#ifdef use_NCURSES
    top_header = "\t\t%%CPU %%MEM COMMAND             TIME   PID";
#else
    top_header = "\t\t%CPU %MEM COMMAND             TIME   PID";
#endif
    top_command = "ps -p " + ProcessID + " -o pcpu,pmem,comm,time,pid | grep "
            + ProcessID;
    top_data = " ";
    step_header
            = "\rstep   s/steps\t  dt /sim t  \tWC t /CPU t\tCFL\t# cells";

    string out = out_stream.str();
    suffix = "";
#ifdef use_NCURSES
    mvprintw(0, 0, "%s", out.c_str());
    refresh();
#else
    if (WRITE_TO_SCREEN)
        cout << out.c_str();
    if (WRITE_TO_SCREEN)
        cout << "Version " << StringFromMakefile(REV) << endl << "Compiled on " << StringFromMakefile(DATE_TIME) << endl;
        cout << "Platform identified itself as: " << OS;
    if (WRITE_TO_SCREEN)
        cout << " " << VERS << "; Number of threads: "
            << globalSystem->NumThreads << endl;
#endif
    latest_file.clear();
}

/**************************************************************/
IO::~IO() {
#ifdef use_NCURSES
    end_curses();
#endif
}
/**************************************************************/
#ifdef use_NCURSES

void IO::start_curses() {
    initscr(); /* Start curses mode 		*/
}

/**************************************************************/
void IO::end_curses() {
    endwin(); /* End curses mode		  */
}
/**************************************************************/
#endif

void IO::print_header() {
    if (!setlocale(LC_CTYPE, "")) {
        fprintf(stderr, "Can't set the specified locale! "
                "Check LANG, LC_CTYPE, LC_ALL.\n");
    }

    if (WRITE_TO_SCREEN)
        cout << setfill('*') << setw(80) << "*" << endl;
    if (WRITE_TO_SCREEN)
        cout << "*" << VERS << " ω-V Code " << StringFromMakefile(REV)  << " Copyright© Tom McCombes " << StringFromMakefile(DATE_YEAR) << "\t\t       *" << endl;

    if (WRITE_TO_SCREEN)
        cout << setfill('*') << setw(80) << "*" << endl;

}

/**************************************************************/
void IO::print_line() {
    if (WRITE_TO_SCREEN)
        cout << setfill('*') << setw(80) << "*" << endl;
}

/**************************************************************/
void IO::print_header(ostream& out_stream) {
    if (!setlocale(LC_CTYPE, "")) {
        fprintf(stderr, "Can't set the specified locale! "
                "Check LANG, LC_CTYPE, LC_ALL.\n");
    }
    if (WRITE_TO_FILE)
        out_stream << setfill('*') << setw(81) << "*" << endl << "*";
#ifndef use_NCURSES      
    if (WRITE_TO_FILE)
        out_stream << "\t\t" <<  VERS << "  ω-V Code \t Copyright© Tom McCombes " << StringFromMakefile(DATE_YEAR);
#else 
    if (WRITE_TO_FILE) out_stream << VERS << "  Wake Code v." << REV << " Copyright© Tom McCombes 2013";
#endif
    if (WRITE_TO_FILE)
        out_stream << "\t\t*" << endl << setfill('*') << setw(81) << "*" << endl;

}

/**************************************************************/
string StringToLower(string strToConvert) {//change each element of the string to lower case
    //    locale loc;
    //
    //    for (unsigned int i = 0; i < strToConvert.length(); i++) {
    //        strToConvert[i] = tolower(strToConvert[i], loc);
    //    }
    return strToConvert; //return the converted string
}

/**************************************************************/
string StringToUpper(string strToConvert) {//change each element of the string to lower case
    //    locale loc;
    //
    //    for (unsigned int i = 0; i < strToConvert.length(); i++) {
    //        strToConvert[i] = toupper(strToConvert[i], loc);
    //    }
    return strToConvert; //return the converted string
}

/**************************************************************/
void IO::read_neu(string infname, Array<Vect3> &X, Array<Array<int> > &PNLS,
        Array<Array<int> > &GROUPS, Array<Array<int> > &BCS,
        Array<string> &NAMES) {
    ifstream input;
    input.open(infname.c_str());
    if (!input) {
#ifndef use_NCURSES  
        if (WRITE_TO_SCREEN)
            cout << "Unable to open file: " << infname << endl;
        throw NO_FILE;
        return;
#endif
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
            if (WRITE_TO_SCREEN)
                cout << fname << endl;
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

        for (int i = 0; i < NBSETS; ++i) {

            getline(input, line); //  Read BOUNDARY CONDITIONS 2.4.6
            getline(input, line); //  Read Group Name; something, number of elements and some other stuff
            {
                istringstream strm(line);
                strm >> NAME >> t2 >> NELEM; // Read group name and number of elements
            }
#ifndef use_NCURSES     
            if (WRITE_TO_SCREEN)
                cout << NAME << " " << t2 << " " << NELEM << endl;
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
}
/**************************************************************/
static Array<Vect3> ReadVectorsFromLine(string);
Array<Vect3> ReadVectorsFromLine(string line) {
    Array<Vect3> output;
    size_t open = line.find_first_of("[");
    while (open != string::npos) {
        //  Read to first closing bracket
        string vect;
        size_t close = line.find_first_of("]", open + 1);
        for (int i = int(open) + 1; i < int(close); ++i)
            vect.push_back(line[i]);

        line.erase(0, close + 1);
        istringstream strm(vect);
        Vect3 temp;
        strm >> temp.x >> temp.y >> temp.z;
        output.push_back(temp);
        open = line.find_first_of("[");
    }
    return output;
}

static string ChopLine(string);

string ChopLine(string line) {
    line.erase(0, line.find_first_of(":") + 1);
    line.erase(line.find_first_of(";"), line.length());
    return line;
}
/**************************************************************/
void IO::PrepOutputDir() {
    stringstream out_stream;

    directory = globalSystem->WorkingDir + "output/";

    string command = "mkdir -p " + directory + "binary_output/";
    string output = globalGetStdoutFromCommand(command);
    cout << command << " " << output << endl;
    command = "mkdir -p " + directory;

    out_name = "_output_pid" + to_string(globalSystem->ProcessID);
    file_type = ".dat";
    image_type = ".png";

    string check = "touch -c " + directory + "test";
    int j = system(check.c_str());
    if (j == 0) {
#ifndef use_NCURSES
        if (WRITE_TO_FILE)
            out_stream << "Success!" << endl;
#endif
    } else {
        if (WRITE_TO_FILE)
            out_stream << endl << "Unable to access output directory, aborting"
                << endl;
        throw NO_FILE;
    }
}
/**************************************************************/
//void IO::read_input(string infname) {
//    ifstream input;
//    input.open(infname.c_str());
//    if (!input) {
//        if (WRITE_TO_SCREEN)
//            cout << "Unable to open input file: " << infname << endl;
//        return;
//    }
//    if (input.is_open()) {
//        while (!input.eof()) {
//            string line;
//            getline(input, line);
//            istringstream strm;
//            if (!line.empty()) {
//                if ((line[0] != '%') && (line[0] != '/')) {
//                    if ((line[0] == '#') && WRITE_TO_SCREEN)
//                        cout << line << endl;
//                    else
//                        strm.str(ChopLine(line));
//
//                    size_t found;
//                    found=line.find_first_of(":");
//
//                    string ln = line.substr(0,found);
//                    if (line.compare(0,found,"Input Neu File") == 0)
//                        strm >> globalSystem->NeuFile;
//                    if (line.compare(0,found,"max P") == 0) {
//                        strm >> SYSTEM::MaxP;
//                    }
//                    if (line.compare(0,found,"Scale") == 0)
//                        strm >> SYSTEM::GambitScale;
//                    if (line.compare(0,found,"dt Init") == 0)
//                        strm >> globalSystem->dtInit;
//                    if (line.compare(0,found,"ds") == 0){
//                        strm >> globalSystem->DS;
//                    }
//                    if (line.compare(0,found,"Num Sub-Steps") == 0){
//                        strm >> globalSystem->NumSubSteps;
//                    }
//                    if (line.compare(0,found,"Case Name") == 0)
//                        strm >> globalSystem->CaseName;
//                    if (line.compare(0,found,"Num Bodies") == 0) {
//                        strm >> globalSystem->NumBodies;
//                    }
//                    if (line.compare(0,found,"Body CGs") == 0)
//                        globalSystem->ORIGIN = ReadVectorsFromLine(ChopLine(
//                            line));
//                    if (line.compare(0,found,"Body Rates") == 0)
//                        globalSystem->RATES = ReadVectorsFromLine(
//                            ChopLine(line));
//                    if (line.compare(0,found,"Body Kin-Vels") == 0)
//                        globalSystem->VELOCITY = ReadVectorsFromLine(ChopLine(
//                            line));
//                    if (line.compare(0,found,"Body Attitudes") == 0)
//                        globalSystem->ATTITUDE = ReadVectorsFromLine(ChopLine(line));
//                    if (line.compare(0,found,"Vinf") == 0) {
//                        Array <Vect3> temp = ReadVectorsFromLine(ChopLine(line));
//
//                        globalSystem->Vinf = temp[0];
//                    }
//                }
//            }
//        }
//    }
//}

/**************************************************************/
void IO::read_dat(char* infname, Array<Vect3> &x, Array<Vect3> &omega) {
    ifstream input;
    input.open(infname);
    if (!input) {
        if (WRITE_TO_SCREEN)
            cout << "Unable to open input file: " << infname << endl;
        return;
    }

    if (input.is_open()) {
        while (!input.eof()) {
            string line;
            getline(input, line);
            istringstream strm(line);
            if (!line.empty()) {
                Vect3 X, OM;
                strm >> X.x >> X.y >> X.z >> OM.x >> OM.y >> OM.z;
                //                X = floor(X) + .5;
                x.push_back(X);
                omega.push_back(OM);
            }
        }

    }
}

/**************************************************************/
#ifdef _PNGWRITER

void IO::create_image(string outname) {
    // Create the PNGwriter instance.

    CurrentPNG = pngwriter(NX, NY, 1.0, outname.c_str());
    //    CurrentPNG.filledsquare_blend(5, 5, NX-5, NY-5, 1., 0.875, 1., 0.125);
    //    x.clear(); y.clear(); z.clear();
    //    domain.assign(NX,Array<REAL> (NY,0.));
    //    globalOctree->Root->ApplyRecursively(&Node::DoNothing, &FVMCell::ImageToIO, &Node::DoNothing);
    //
    //    if (Plot)
    //    {
    //        Plot->reset_plot();
    //        Plot->unset_grid();
    //    }
    //    else
    //        Plot = new Gnuplot("lines");
    //
    //    Plot->plot_xyz(x,y,z,"something");
    //
    ////    write_file();
    //
    //    // Write the file to disk.
    //    CurrentPNG.plot(floor(globalTimeStepper->centre.x), floor(globalTimeStepper->centre.y)+1, 1, 1, 1);
    //    CurrentPNG.plot(floor(globalTimeStepper->centre.x)+1, floor(globalTimeStepper->centre.y), 1, 1, 1);
    //    CurrentPNG.plot(floor(globalTimeStepper->centre.x), floor(globalTimeStepper->centre.y)-1, 1, 1, 1);
    //    CurrentPNG.plot(floor(globalTimeStepper->centre.x)-1, floor(globalTimeStepper->centre.y), 1, 1, 1);
    //    CurrentPNG.plot(floor(globalTimeStepper->centre.x), floor(globalTimeStepper->centre.y), 1, 1, 1);
    CurrentPNG.close();
    suffix = "Written image: " + outname;
}
#endif

/**************************************************************/
void IO::stat_step() {
    #ifdef TIME_STEPS
    long unsigned int t8 = ticks();
#endif
    stringstream out_stream;

#ifndef use_NCURSES
    if ((globalTimeStepper->n % HEADER_OUTPUT == 0) || (globalTimeStepper->show_rundata)) {
#ifdef TOP
        top_data = "\t\t" + globalGetStdoutFromCommand(top_command);
        out_stream << top_header << endl << top_data << endl;
#endif
        out_stream << step_header << endl;
        globalTimeStepper->show_rundata = false;
    }
#endif
    out_stream.setf(ios::fixed, ios::floatfield);
    out_stream.precision(4);
    out_stream << globalTimeStepper->n << "\t" << globalSystem->NumSubSteps
            << "\t" << globalTimeStepper->dt << "/" << TIME_STEPPER::SimTime;
    out_stream.precision(3);
    out_stream << "\t" << (REAL) (ticks() - globalTimer) / 1000;
    out_stream << "/" << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    globalTimeStepper->cpu_t = ticks();

    out_stream << "\t" << globalTimeStepper->CFL.Mag();

    out_stream << "\t" << FVMCell::NumCells;
    if (globalTimeStepper->dump_next) {
        out_stream << "\t   *";
    } else {
        out_stream << "\t   ";
    }
    step_data += out_stream.str();

    if (WRITE_TO_SCREEN)
        display_stat_step();

    double MEM_PERCENT, temp;
    stringstream psdata;
    psdata << top_data;
    psdata >> temp >> MEM_PERCENT;

    if (MEM_PERCENT > 95) {
        cout << setfill('!') << setw(80) << "!" << endl;

        cout << "Out of memory. Quitting to avoid swapping to disk." << endl;
        cout << "Memory used: " << MEM_PERCENT << "%" << endl;

        cout << setfill('!') << setw(80) << "!" << endl;
        throw OutOfMemory();


    }


    //            if (WRITE_TO_FILE) write_stat_step();
    #ifdef TIME_STEPS
    long unsigned int t9 = ticks();

    stringstream tmp;
    tmp << "stat_step                : " << double(t9 - t8) / 1000.0 << endl;
    step_data += tmp.str();
#endif
}

/**************************************************************/
void IO::display_stat_step() {
    if (WRITE_TO_SCREEN) {

#ifdef use_NCURSES
        if (lines.size() > HEADER_OUTPUT) lines.pop_back();
        lines.push_front(step_data);
        mvprintw(5, 0, "%s", step_header->c_str());
        int count = 0;
        for (list<string>::iterator it = lines.begin(); it != lines.end(); ++it) {
            mvprintw(6 + count, 0, "%s", it->c_str());
            count++;
        }
        mvprintw(6 + HEADER_OUTPUT + 3, 0, "%s", suffix.c_str());

#ifdef TOP
        mvprintw(19, 0, top_header.c_str());
        mvprintw(20, 0, top_data.c_str());
#endif

        //for (int i = 0; i < lines.size(); ++i) mvprintw(4+i,0,"%s",out.c_str());
        //mvprintw(4,0,"%s",out.c_str());
        refresh();
#else
        if (WRITE_TO_SCREEN) {
            cout << step_data.c_str() << " " << latest_file << endl;
            latest_file.clear();
            step_data.clear();
        }
#endif
    }

}

/**************************************************************/
void IO::print_stat_step(ostream & out) {
    out << "*" << endl;
    out.setf(ios::fixed, ios::floatfield);
    out << "*\tstep:\t" << globalTimeStepper->n << " \tsim time:\t"
            << globalTimeStepper->sim_time << "\tdt\t" << globalTimeStepper->dt
            << endl;
    out << "*\tcpu_t:\t" << (REAL) (ticks() - globalTimeStepper->cpu_t) / 1000;
    out << "\tt:\t" << (REAL) (ticks() - globalTimer) / 1000;
    out << "\t# cells:\t" << FVMCell::NumCells << endl;
    out << "*" << endl;
}

/**************************************************************/
void IO::write_file(stringstream &outstream, string OutName, string ext,
        bool disp) {

    int ID = -1;
    for (int i = 0; i < NumFiles.size(); ++i)
        if (FileNames[i] == OutName)
            ID = i;

    if (ID == -1) {
        FileNames.push_back(OutName);
        NumFiles.push_back(0);
        ID = NumFiles.size() - 1;
    }

    string fname = directory + FileNames[ID];
    stringstream num;
    num << NumFiles[ID];

    int pad = 6 - int(num.str().length());
    fname += "_" + globalSystem->CaseName + string(pad, '0') + num.str() + "."
            + ext;
    num_file++;
    if (WRITE_TO_SCREEN && disp)
        cout << fname << endl;
    fstream filestr;

    filestr.open(fname.c_str(), fstream::in | fstream::out | fstream::app);
    if (filestr.is_open())
        filestr << outstream.str() << endl;
    else {
        if (WRITE_TO_SCREEN)
            cout << "Unable to open output file: " << fname << endl;
        throw NO_FILE;
    }

    NumFiles[ID]++;
    filestr.close();
}
/**************************************************************/
void IO::writeMATLABOutputStruct(MATLABOutputStruct &outdata, string OutName) {

    int ID = -1;
    for (int i = 0; i < NumFiles.size(); ++i)
        if (FileNames[i] == OutName)
            ID = i;

    if (ID == -1) {
        FileNames.push_back(OutName);
        NumFiles.push_back(0);
        ID = NumFiles.size() - 1;
    }

    string fname = directory + FileNames[ID];
    stringstream num;
    num << NumFiles[ID];

    int pad = 6 - int(num.str().length());
    fname += "_" + globalSystem->CaseName + string(pad, '0') + num.str() + ".mat";
    num_file++;
    latest_file = "<-W " + fname;

    fstream filestr;

    for (int i = 0; i < outdata.Int1DArrays.size(); ++i)
        UTIL::WriteMATLABMatrix1D(outdata.Int1DArrayStrings[i], fname, outdata.Int1DArrays[i]);

    for (int i = 0; i < outdata.Double1DArrays.size(); ++i)
        UTIL::WriteMATLABMatrix1D(outdata.Double1DArrayStrings[i], fname, outdata.Double1DArrays[i]);

    for (int i = 0; i < outdata.Vect1DArrays.size(); ++i)
    {
        Array < Array < REAL> > data = UTIL::zeros(outdata.Vect1DArrays[i].size(), 3);
        for (int j = 0; j < outdata.Vect1DArrays[i].size(); ++j)
        {
            data[j][0] = outdata.Vect1DArrays[i][j].x;
            data[j][1] = outdata.Vect1DArrays[i][j].y;
            data[j][2] = outdata.Vect1DArrays[i][j].z;
        }
        UTIL::WriteMATLABMatrix2D(outdata.Vect1DArrayStrings[i], fname, data);
    }

    for (int i = 0; i < outdata.Int2DArrays.size(); ++i)
        UTIL::WriteMATLABMatrix2D(outdata.Int2DArrayStrings[i], fname, outdata.Int2DArrays[i]);

    for (int i = 0; i < outdata.Double2DArrays.size(); ++i)
        UTIL::WriteMATLABMatrix2D(outdata.Double2DArrayStrings[i], fname, outdata.Double2DArrays[i]);

    for (int i = 0; i < outdata.Vect2DArrays.size(); ++i)
        UTIL::WriteMATLABMatrix2DVect3(outdata.Vect2DArrayStrings[i], fname, outdata.Vect2DArrays[i]);

    if (globalSystem->useBodies) {
        UTIL::PostAmble(fname);
    }
    UTIL::WriteMATLABMatrix1D("Del2", fname, globalSystem->Del2);
    UTIL::WriteMATLABMatrix1D("h", fname, SYSTEM::M4Radius_in_cells);
    UTIL::WriteMATLABMatrix1D("GambitScale", fname, SYSTEM::GambitScale);
    UTIL::WriteMATLABMatrix1D("MaxP", fname, SYSTEM::MaxP);
    UTIL::WriteMATLABMatrix1D("KinematicViscosity", fname, globalSystem->Nu);
    UTIL::WriteMATLABMatrix1D("LastFVMDeltaT", fname, globalTimeStepper->dt);
    NumFiles[ID]++;
    filestr.close();
}
/**************************************************************/
#ifdef _PNGWRITER

void IO::write_image() {
    stringstream temp;
    temp << num_images;
    string filename_padding, prefix = temp.str(), image;
    if (num_images < 100000) filename_padding = "0";
    if (num_images < 10000) filename_padding = "00";
    if (num_images < 1000) filename_padding = "000";
    if (num_images < 100) filename_padding = "0000";
    if (num_images < 10) filename_padding = "00000";
    num_images++;

    image = directory + filename_padding + prefix /*+ out_name*/ + image_type;
    create_image(image);
}
#endif

/**************************************************************/
void IO::write_m() {
//    ofstream out_stream;
//    stringstream temp;
//    temp << num_m;
//    string filename_padding, prefix = temp.str();
//    if (num_m < 100000)
//        filename_padding = "dump_0";
//    if (num_m < 10000)
//        filename_padding = "dump_00";
//    if (num_m < 1000)
//        filename_padding = "dump_000";
//    if (num_m < 100)
//        filename_padding = "dump_0000";
//    if (num_m < 10)
//        filename_padding = "dump_00000";
//    num_m++;
//    string f_type = ".m";
//    filename = directory + filename_padding + prefix + f_type;
//    out_stream.open(filename.c_str());
//    if (!out_stream) {
//        if (WRITE_TO_SCREEN)
//            cout << "Unable to open output file: " << filename << endl;
//        throw NO_FILE;
//    }
//    if (!setlocale(LC_CTYPE, "")) {
//        fprintf(stderr, "Can't set the specified locale! "
//                "Check LANG, LC_CTYPE, LC_ALL.\n");
//    }
//
//    globalSystem->WriteBodiesAndWakes(out_stream);
//
//    out_stream.close();
}

/**************************************************************/
void IO::write_forces() {
}

/**************************************************************/
void IO::write_vels_and_pressures() {
}

/**************************************************************/
void IO::write_GMSH() {
}

/**************************************************************/
void IO::WriteBinary() {
    /*
            struct Body {
                    char Name[80];
                    int NPanels;
                    Vect3 Position;
                    Vect3 Velocity;
                    Vect3 Attitude;
                    Vect3 Rates;
            };

            Array<Body> Bodies(123);

            for (int i = 0; i < Bodies.size(); ++i)
                    for (int j = 0; j < 3; ++j) {
                            string Name = "Body " + to_string(i);
                            Name.insert(Name.length(), 80 - Name.length(), ' ');
                            strcpy(Bodies[i].Name, Name.c_str());
                            Bodies[i].Velocity[j] = j + .5;
                            Bodies[i].Attitude[j] = 10 * (i + 1);
                            Bodies[i].Position[j] = j + pi;
                            Bodies[i].Rates[j] = 10 * (j + 1);
                            Bodies[i].NPanels = 123 + i;
                    }

            ofstream fout("bodies.bin", ios::out | ios::binary | ios::ate);
            if (!fout) {
                    cout << "Cannot open file.\n";
            }
            int NumBodies = Bodies.size();

            fout.write((char *) &NumBodies, sizeof(int));

            for (int i = 0; i < NumBodies; ++i) {
                    fout.write((char *) &Bodies[i].Name, sizeof(char[80]));
                    fout.write((char *) &Bodies[i].NPanels, sizeof(int));
                    fout.write((char *) &Bodies[i].Position, sizeof(Vect3));
                    fout.write((char *) &Bodies[i].Velocity, sizeof(Vect3));
                    fout.write((char *) &Bodies[i].Attitude, sizeof(Vect3));
                    fout.write((char *) &Bodies[i].Rates, sizeof(Vect3));
            }

            int x = 50, y = 75, z = 100;
            fout.write((char *) &x, sizeof(int));
            fout.write((char *) &y, sizeof(int));
            fout.write((char *) &z, sizeof(int));
            Array<Array<Array<REAL> > > domain(x);
            Array<Array<Array<Vect3> > > geom(x);
            int count = 0;
            for (int i = 0; i < x; ++i) {
                    domain[i].allocate(y);
                    geom[i].allocate(y);
                    for (int j = 0; j < y; ++j) {
                            domain[i][j].allocate(z);
                            geom[i][j].allocate(z);
                            for (int k = 0; k < z; ++k) {
                                    domain[i][j][k] = count + pi;
                                    geom[i][j][k] = Vect3(i + .5, j + .5, k + .5);
                                    count++;
                            }
                    }
            }

            for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                            for (int k = 0; k < z; ++k)
                                    fout.write((char *) &domain[i][j][k], sizeof(REAL));

            for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                            for (int k = 0; k < z; ++k)
                                    fout.write((char *) &geom[i][j][k].x, sizeof(REAL));

            for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                            for (int k = 0; k < z; ++k)
                                    fout.write((char *) &geom[i][j][k].y, sizeof(REAL));

            for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                            for (int k = 0; k < z; ++k)
                                    fout.write((char *) &geom[i][j][k].z, sizeof(REAL));

            fout.close();
     */
}
