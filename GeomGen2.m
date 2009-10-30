clear all; close all
cd ~/Desktop/Combined
dir1 = '~/Desktop/Combined/neu_files/'



%%  Geometry Parameters
cutout = 1.257;
numbodies = 2;
nblades = [3 4];
Reflect = [1 -1];
%%  Simulation Parameters
MaxP = 3;
Scale = 3;          %   Scaling is done in the simulation
Vels{1} = [-10 0 0];
Vels{2} = [-10 0 0];
Origin{1} = [0 0 0];
Origin{2} = [2.5 0 0];
Attitudes{1} = [0 0 0];
Attitudes{2} = [0 0 0];
Rates{1} = [0 0 0];
Rates{2} = [0 0 0];




name = 'coax';
fname = [name '.neu']
fid = fopen(fname, 'wt');
fid_cas = fopen([name '.cas'],'wt');


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
% RADIUS = 0:5;%
% CHORD = ones(size(RADIUS));
% THETA = zeros(size(RADIUS));

%   Span lies along y axis, chord parallel to x axis
n = 30;
radius = max(RADIUS);


th0 = (3-THETA(end));
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
    for i = 1:nblades(q)
        count = count + 1;
        psi = (i-1)*dpsi;
        
        rads = linspace(cutout,radius,n)';
        C = interp1(RADIUS,CHORD,rads,'cubic');
        LEx = -pitch_axis*C;
        TEx = (1-pitch_axis)*C;
        th = th0 + interp1(RADIUS,THETA,rads,'cubic');
        
        xdata0 = [LEx TEx];
        ydata0 = [rads rads];
        zdata0 = zeros(size(ydata0));
        
        th = [th th];
        xdata1 = Reflect(q)*xdata0.*cosd(th);
        ydata1 = ydata0;
        zdata1 = xdata0.*sind(th);
        
        xdata2 = xdata1*cosd(psi) - ydata1*sind(psi);
        ydata2 = xdata1*sind(psi) + ydata1*cosd(psi);
        zdata2 = zdata1;
        
        
        zdata = xdata2 + origin(3);
        ydata = ydata2 + origin(2);
        xdata = zdata2 + origin(1);
        
        
        
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
CaseData = AddVal2Case(CaseData,'NUMBODIES',numbodies,'int');
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


