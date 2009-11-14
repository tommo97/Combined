clear all; close all; clc;
cd ~/Desktop/Workspace/Combined

dir2 = './case_files/';
dir3 = './run_files/';
dir4 = './mat_files/';
commit = false;


split = true;


%%  Simulation Parameters
num_procs = 4;
MaxP = 3;
Scale = 1;          %   Scaling is done in the simulation




name = 'Straight_';
fname = [name '.neu'];


%%  Geometry Parameters - Positions/Vels/Attitudes etc

th0 = 180;

Vels{1} = [10 0 -2];
Origin{1} = [0 0 0];
Attitudes{1} = [0 0 pi];
Rates{1} = [0 0 0];

Vels{2} = [-10 0 0];
Origin{2} = [0 0 0];
Attitudes{2} = [deg2rad(180)+pi/2 0 pi];
Rates{2} = [7.5 0 0];

Vels{3} = [-10 0 0];
Origin{3} = [0 0 0];
Attitudes{3} = [deg2rad(240) 0 pi];
Rates{3} = [7.5 0 0];


%%  Geometry Parameters - make a template/skeleton blade
NRELBlade.Attitude = [0 pi/2 0];
NRELBlade.Origin = [0 0 0];
NRELBlade.th0 =  th0;
NRELBlade.PitchAxis = 0.3;
%   Aerofoil
NRELBlade.NChord = 25;
NRELBlade.NSpan = 26;
%%  Aerofoil

x = BellShape(0,1,NRELBlade.NChord,5);
%x = linspace(0,1,NRELBlade.NChord);
[Aerofoil z] = NRELFoil(x);

NRELBlade.FOIL = z.S809;

%%  NRELBlade -- NREL data
RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
    2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
    0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
    3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];


RADIUS = linspace(-5,5);
THETA = 10*zeros(size(RADIUS));


CHORD = 1*ones(size(RADIUS));

th = linspace(0,pi);
%CHORD = 2*sqrt(sin(linspace(0,pi))) + .1;

NRELBlade.RADIUS = RADIUS;
NRELBlade.CHORD = CHORD;
NRELBlade.THETA = THETA;




%%  Make rotor systems
hold on
NumPanels = 0;
NumPoints = 0;
for i = 1:1
    Blade = NRELBlade;
    Blade.Attitude = Attitudes{i};
    Blade.Velocity = Vels{i};
    Blade.Rates = Rates{i};
    Blade.Origin = Origin{i};
    Bodies{i} = MakeBlade(Blade);
    S = zeros(size(Bodies{i}.X(Bodies{i}.N.Local)));
    for j = 1:size(S,2)
        S(:,j) = j;
    end
    surf(Bodies{i}.X(Bodies{i}.N.Local),Bodies{i}.Y(Bodies{i}.N.Local),Bodies{i}.Z(Bodies{i}.N.Local),S);
    hold on
     surf(Bodies{i}.X(Bodies{i}.Tip.Inboard.US.N.Local),Bodies{i}.Y(Bodies{i}.Tip.Inboard.US.N.Local),Bodies{i}.Z(Bodies{i}.Tip.Inboard.US.N.Local));
     surf(Bodies{i}.X(Bodies{i}.Tip.Inboard.LS.N.Local),Bodies{i}.Y(Bodies{i}.Tip.Inboard.LS.N.Local),Bodies{i}.Z(Bodies{i}.Tip.Inboard.LS.N.Local));
     surf(Bodies{i}.X(Bodies{i}.Tip.Outboard.US.N.Local),Bodies{i}.Y(Bodies{i}.Tip.Outboard.US.N.Local),Bodies{i}.Z(Bodies{i}.Tip.Outboard.US.N.Local));
     surf(Bodies{i}.X(Bodies{i}.Tip.Outboard.LS.N.Local),Bodies{i}.Y(Bodies{i}.Tip.Outboard.LS.N.Local),Bodies{i}.Z(Bodies{i}.Tip.Outboard.LS.N.Local));
    %
    Bodies{i}.N.Global = Bodies{i}.N.Local + NumPoints;
    Bodies{i}.Panels.c1.Global = Bodies{i}.Panels.c1.Local + NumPoints;
    Bodies{i}.Panels.c2.Global = Bodies{i}.Panels.c2.Local + NumPoints;
    Bodies{i}.Panels.c3.Global = Bodies{i}.Panels.c3.Local + NumPoints;
    Bodies{i}.Panels.c4.Global = Bodies{i}.Panels.c4.Local + NumPoints;
    Bodies{i}.Panels.WakeShedders.US.Global = Bodies{i}.Panels.WakeShedders.US.Local + NumPanels;
    Bodies{i}.Panels.WakeShedders.LS.Global = Bodies{i}.Panels.WakeShedders.LS.Local + NumPanels;
    NumPanels = NumPanels + numel(Bodies{i}.Panels.c1.Local);
    NumPoints = NumPoints + numel(Bodies{i}.X);
    Attitudes{i} = [0 0 0];
    %   The following are some TINY errors, they serve to make the
    %   influence coefficient matrices non-singular. Why and how are
    %   unknown to me...
    %     Bodies{i}.X = Bodies{i}.X + rand*1e-16;
    %     Bodies{i}.Y = Bodies{i}.Y + rand*1e-16;
    %     Bodies{i}.Z = Bodies{i}.Z + rand*1e-16;
    
    
    
    
    
    
    
    
    
    
    Bodies{i}.Faces.C1.Body = [Bodies{i}.X(Bodies{i}.Panels.c1.Global) Bodies{i}.Y(Bodies{i}.Panels.c1.Global) Bodies{i}.Z(Bodies{i}.Panels.c1.Global)];
    Bodies{i}.Faces.C2.Body = [Bodies{i}.X(Bodies{i}.Panels.c2.Global) Bodies{i}.Y(Bodies{i}.Panels.c2.Global) Bodies{i}.Z(Bodies{i}.Panels.c2.Global)];
    Bodies{i}.Faces.C3.Body = [Bodies{i}.X(Bodies{i}.Panels.c3.Global) Bodies{i}.Y(Bodies{i}.Panels.c3.Global) Bodies{i}.Z(Bodies{i}.Panels.c3.Global)];
    Bodies{i}.Faces.C4.Body = [Bodies{i}.X(Bodies{i}.Panels.c4.Global) Bodies{i}.Y(Bodies{i}.Panels.c4.Global) Bodies{i}.Z(Bodies{i}.Panels.c4.Global)];
    
    Bodies{i}.Faces.CP.Body = .25*(Bodies{i}.Faces.C1.Body + Bodies{i}.Faces.C2.Body +...
        Bodies{i}.Faces.C3.Body + Bodies{i}.Faces.C4.Body);
    Bodies{i}.Faces.D1.Body = Bodies{i}.Faces.C3.Body - Bodies{i}.Faces.C1.Body;
    Bodies{i}.Faces.D2.Body = Bodies{i}.Faces.C4.Body - Bodies{i}.Faces.C2.Body;
    
    
    ly =  .5*(Bodies{i}.Faces.C1.Body+Bodies{i}.Faces.C4.Body) - .5*(Bodies{i}.Faces.C2.Body+Bodies{i}.Faces.C3.Body);
    lymag = sqrt(dot(ly,ly,2));
    lz = cross(Bodies{i}.Faces.C2.Body-Bodies{i}.Faces.C4.Body, Bodies{i}.Faces.C3.Body-Bodies{i}.Faces.C1.Body);
    lzmag = sqrt(dot(lz,lz,2));
    Bodies{i}.Faces.Area = .5*lzmag;
    Bodies{i}.Faces.LocalAxis.Y.Body = [ly(:,1)./lymag, ly(:,2)./lymag, ly(:,3)./lymag];
    Bodies{i}.Faces.LocalAxis.Z.Body = [lz(:,1)./lzmag, lz(:,2)./lzmag, lz(:,3)./lzmag];
    Bodies{i}.Faces.LocalAxis.X.Body = cross(Bodies{i}.Faces.LocalAxis.Y.Body, Bodies{i}.Faces.LocalAxis.Z.Body);
%    scatter3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3));
%     quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
%         Bodies{i}.Faces.LocalAxis.X.Body(:,1),Bodies{i}.Faces.LocalAxis.X.Body(:,2),Bodies{i}.Faces.LocalAxis.X.Body(:,3),'green');
%     quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
%         Bodies{i}.Faces.LocalAxis.Y.Body(:,1),Bodies{i}.Faces.LocalAxis.Y.Body(:,2),Bodies{i}.Faces.LocalAxis.Y.Body(:,3),'red');
%     quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
%         Bodies{i}.Faces.LocalAxis.Z.Body(:,1),Bodies{i}.Faces.LocalAxis.Z.Body(:,2),Bodies{i}.Faces.LocalAxis.Z.Body(:,3),'blue');
%     
    
    
    
    
    
    
    
end

%%  Prepare output to neutral file
MakeNEU(Bodies,name,split);

if split
    numbodies = length(Bodies);
else
    numbodies = 1;
end
%%  Write case file
fid_cas = fopen([dir2 name '.cas'],'wt');

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
fid_run = fopen([dir3 name],'wt');
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

save([dir4 name '.mat'],'Bodies')

axis equal tight
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');