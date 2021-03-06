clear all; close all
cd ~/Desktop/Workspace/Combined
dir1 = './neu_files/';
dir2 = './case_files/';
dir3 = './run_files/';
commit = false;
%%  Geometry Parameters
cutout = 0;
numbodies = 1;
nblades = [1 1 1 1 1 1];
Reflect = [-1 1 1 -1 -1 -1];
%%  Simulation Parameters
num_procs = 4;
MaxP = 3;
Scale = 5;          %   Scaling is done in the simulation




Vels{1} = [-1.70 0 0];
Origin{1} = [0 0 0];
Attitudes{1} = [0 deg2rad(8.5) 0];
Rates{1} = [0 0 0];


%Vels{2} = [-10 0 0];
%Vels{3} = [-10 0 0];
%
% Origin{2} = [4 0 0];
% Origin{3} = [4 0 0];
% 
% Attitudes{2} = [pi/2 0 deg2rad(120)];
% Attitudes{3} = [pi/2 0 deg2rad(240)];
% 
% Rates{2} = [0 0 7.5];
% Rates{3} = [0 0 7.5];
% 
% Vels{4} = [-10 0 0];
% Vels{5} = [-10 0 0];
% Vels{6} = [-10 0 0];
% Origin{4} = [4 0 0];
% Origin{5} = [4 0 0];
% Origin{6} = [4 0 0];
% Attitudes{4} = [pi/2 0 deg2rad(0)];
% Attitudes{5} = [pi/2 0 deg2rad(120)];
% Attitudes{6} = [pi/2 0 deg2rad(240)];
% Rates{4} = [0 0 -7.5];
% Rates{5} = [0 0 -7.5];
% Rates{6} = [0 0 -7.5];


name = '0015_';
fname = [name '.neu'];
fid = fopen(fname, 'wt');
fid_cas = fopen([dir2 name '.cas'],'wt');
fid_run = fopen([dir3 name],'wt');

v = ver;
for k=1:length(v)
    if strfind(v(k).Name, 'MATLAB')
        version = v(k).Version;
        break
    end
end

%%  Define blade geometry
%  NREL Data
RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
    2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];

CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
    0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
    3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];


%   Span lies along y axis, chord parallel to x axis
n = 128;

RADIUS = linspace(0, 4,n);
CHORD = ones(size(RADIUS));
THETA = zeros(size(RADIUS));


radius = max(RADIUS);


th0 = 0;%-(3-THETA(end));
pitch_axis = 0.3; % chord


M1 = [];
M2 = [];
cnt = 0;
count = 0;
%%
for q = 1:numbodies
    dpsi = 360/nblades(q);
    origin = Origin{q};
    g{q}.start = cnt + 1;
    EulerAngles = Attitudes{q};
    cosphi = cos(EulerAngles(1)); costhe = cos(EulerAngles(2)); cospsi = cos(EulerAngles(3));
    sinphi = sin(EulerAngles(1)); sinthe = sin(EulerAngles(2)); sinpsi = sin(EulerAngles(3));
    a1 = costhe*cospsi;
    a2 = costhe*sinpsi;
    a3 = -sinthe;
    b1 = sinphi*sinthe*cospsi - cosphi*sinpsi;
    b2 = sinphi*sinthe*sinpsi + cosphi*cospsi;
    b3 = sinphi*costhe;
    c1 = cosphi*sinthe*cospsi + sinphi*sinpsi;
    c2 = cosphi*sinthe*sinpsi - sinphi*cospsi;
    c3 = cosphi*costhe;
    TRANS(1,:) = [a1 b1 c1];
    TRANS(2,:) = [a2 b2 c2];
    TRANS(3,:) = [a3 b3 c3];
    for i = 1:nblades(q)
        count = count + 1;
        
        rads = linspace(cutout,radius,n)';
        chord = interp1(RADIUS,CHORD,rads,'cubic');        
        
        LEx = pitch_axis*chord;
        TEx = -(1-pitch_axis)*chord;
        th = th0 + interp1(RADIUS,THETA,rads,'cubic');
        
        xdata0 = [LEx TEx];
        ydata0 = [rads rads];
        zdata0 = zeros(size(ydata0));
        
        th = [th th];
        xdatab = Reflect(q)*(xdata0.*cosd(th) + zdata0.*sind(th)) + origin(1);
        ydatab = ydata0 + origin(2);
        zdatab = -xdata0.*sind(th) + zdata0.*cosd(th) + origin(3);
        
        
        xdata = xdatab*TRANS(1,1) + ydatab*TRANS(1,2) + zdatab*TRANS(1,3);
        ydata = xdatab*TRANS(2,1) + ydatab*TRANS(2,2) + zdatab*TRANS(2,3);
        zdata = xdatab*TRANS(3,1) + ydatab*TRANS(3,2) + zdatab*TRANS(3,3);
        
%         if q <4
%             ydata = ydata - 6;
%         else
%             ydata = ydata + 6;
%         end
        
        hold all
        surf(xdata,ydata,zdata)
        axis equal
        
        
        C = xdata;
        
        C(:) = (1:numel(xdata)) + (count-1) * numel(xdata);
        
        
        X = [xdata(:) ydata(:) zdata(:)];
        
        
        
        
        c1 = C(1:end-1,1);
        c2 = C(1:end-1,2);
        c3 = C(2:end,2);
        c4 = C(2:end,1);
        
        M1 = [M1; X];   %   This holds the positions of the points
        M2 = [M2; [c1 c2 c3 c4]];   %   This hols the panel corners
        
        cnt = cnt + numel(c1);
    end
    
    g{q}.end = cnt;
end

scatter3(0,0,0); axis tight; view(3);
nPts = size(M1,1);
nPnls = size(M2,1);
nGrps = numbodies;
nBCs = 1;
Empty = 0;

% Attitudes{1} = [0 0 0];
% Attitudes{2} = [0 0 0];
% Attitudes{3} = [0 0 0];
% Origin{1} = [0 -6 0];
% Origin{2} = [0 -6 0];
% Origin{3} = [0 -6 0];
% 
% Attitudes{4} = [0 0 0];
% Attitudes{5} = [0 0 0];
% Attitudes{6} = [0 0 0];
% Origin{4} = [0 6 0];
% Origin{5} = [0 6 0];
% Origin{6} = [0 6 0];

fprintf(fid,'        CONTROL INFO 2.4.6\n** GAMBIT NEUTRAL FILE\n%s\n',name);
fprintf(fid,'PROGRAM:                MATLAB     VERSION:  %s\n', version);
fprintf(fid,'%s\n', datestr(now));
fprintf(fid,'\tNUMNP\tNELEM\tNGRPS\tNBSETS\tNDFCD\tNDFVL\n');
fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\n',nPts,nPnls,nGrps,nBCs,Empty,Empty);
fprintf(fid,'ENDOFSECTION\n   NODAL COORDINATES 2.4.6\n');


str = [repmat(' ',size(M1,1),10 - size(num2str(size(M1,1)),2)) num2str([1:size(M1,1)]','%g    ') repmat('   ',size(M1,1),1) num2str(M1,'%-1.12e    ')];

for i = 1:size(M1,1)
    fprintf(fid,'%s\n',str(i,:));
    %    fprintf(fid,'\t%g\t%-1.12e\t%-1.12e\t%-1.12e\n',i,M1(i,1),M1(i,2),M1(i,3));
end
fprintf(fid,'ENDOFSECTION\n      ELEMENTS/CELLS 2.4.6\n');
for i = 1:size(M2,1) - 1
    fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',i,2,4,M2(i,1),M2(i,2),M2(i,3),M2(i,4));
end
fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g',size(M2,1),2,4,M2(end,1),M2(end,2),M2(end,3),M2(end,4));

for n = 1:nGrps
    fprintf(fid,'\nENDOFSECTION\n       ELEMENT GROUP 2.4.6\n');
    fprintf(fid,'GROUP:\t\t%g ELEMENTS:\t\t%g MATERIAL:\t\t4 NFLAGS:\t\t1\n\t\t\t%s%g\n\t0\n',n,1 + g{n}.end - g{n}.start,name,n);
    pans = g{n}.start:g{n}.end;
    M3 = zeros(10,floor(numel(pans)/10));
    M3(:) = pans(1:numel(M3));
    M3 = M3';
    M4 = pans(1+numel(M3):end);
    
    for i = 1:size(M3,1)
        fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',M3(i,:));
    end
    
    for i = 1:size(M4,2)
        fprintf(fid,'\t%g',M4(i));
    end
end


fprintf(fid,'\nENDOFSECTION\n BOUNDARY CONDITIONS 2.4.6\n\t\t\t\t%s\t1\t%g\t0\t6\n','WAKE',nPnls);
fprintf(fid,'\t%g\t2\t2\n',1:nPnls);
fprintf(fid,'ENDOFSECTION\n');
fclose(fid);

%%  Write case file
CaseData = struct([]);
CaseData = AddVal2Case(CaseData,'INPUT',fname,'string');
CaseData = AddVal2Case(CaseData,'PMAX',MaxP,'int');
CaseData = AddVal2Case(CaseData,'SCALE',Scale,'double');
CaseData = AddVal2Case(CaseData,'BODYNUM',numbodies,'int');
CaseData = AddVal2Case(CaseData,'ATTITUDES', Attitudes,'double');
CaseData = AddVal2Case(CaseData,'CGBODIES',Origin, 'double');
CaseData = AddVal2Case(CaseData,'RATEBODIES',Rates,'double');
CaseData = AddVal2Case(CaseData,'VELBODIES',Vels,'double');
CaseData = AddVal2Case(CaseData,'NAME',name,'string');


fprintf(fid_cas,'#          Input file generated by MATLAB Version:  %s\n', version);
fprintf(fid_cas,'#          %s\n', datestr(now));


for i = 1:CaseData.NumVals
    fprintf(fid_cas,'%s:\t%s;\n',CaseData.Data{i}.Name, CaseData.Data{i}.Val);
end

fclose(fid_cas);

system(['cat ' dir2 name '.cas']);

%%  Make a runfile
fprintf(fid_run,'%s\n','#!/bin/sh');
fprintf(fid_run,'%s%g\n','#PBS -l nodes=1:ppn=',num_procs);
fprintf(fid_run,'%s\n','#PBS -m bea');
fprintf(fid_run,'%s\n\n','#PBS -M tom.mccombes@strath.ac.uk');
fprintf(fid_run,'%s\n','cd $PBS_O_WORKDIR');
fprintf(fid_run,'%s\n','nprocs=`wc -l $PBS_NODEFILE | awk ''{ print $1 }''`');
fprintf(fid_run,'%s%s%s%s\n','/home/lap05140/./main ',name,'.cas > dump_',name);


%%  Commit them
if commit
    system(['svn add ' dir1 name '.neu ' dir2 name '.cas ' dir3 name]);
    system(['svn commit ' dir1 name '.neu ' dir2 name '.cas ' dir3 name ' -m "MATLAB'...
        ' automatic commit of case/data/runfile for ' name '"']);
end
