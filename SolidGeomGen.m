clear all; close all
cd ~/Desktop/Workspace/Combined
dir1 = './neu_files/';
dir2 = './case_files/';
dir3 = './run_files/';
commit = false;
%%  Geometry Parameters - Aerofoil
NChord = 25;

x = cumtrapz(1 - cos(linspace(0,pi,NChord)).^2)/max(cumtrapz(1 - cos(linspace(0,pi,NChord)).^2));


[Aerofoil z] = NRELFoil(x);

hold on
scatter(x,z.N0012.US);
scatter(x,z.N0012.LS);
axis equal
