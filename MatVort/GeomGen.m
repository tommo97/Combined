clear all
%close all
hold all
cd ~/Desktop/Combined/neu_files
name = 'blade';
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

cutout = 1.257;
th0 = (3-THETA(end));
pitch_axis = 0.3; % chord
nblades = 2;
dpsi = 360/nblades;
M1 = [];
M2 = [];

%%
for i = 1:nblades
psi = (i-1)*dpsi;   

rads = linspace(cutout,radius,n)';
C = interp1(RADIUS,CHORD,rads,'cubic');
LEx = -pitch_axis*C;
TEx = (1-pitch_axis)*C;
th = th0 + interp1(RADIUS,THETA,rads,'cubic');

xdata = [LEx TEx];
ydata = [rads rads];
zdata = zeros(size(ydata));

th = [th th];
xdata = xdata.*cosd(th) - zdata.*sind(th);
zdata = xdata.*sind(th) + zdata.*cosd(th);


Zdata = xdata*cosd(psi) - ydata*sind(psi);
Ydata = xdata*sind(psi) + ydata*cosd(psi);
Xdata = zdata;
hold all
surf(Xdata,Ydata,Zdata)
axis equal



C = Xdata;
C(:) = (1:numel(Xdata)) + (i-1) * numel(Xdata);


X = [Xdata(:) Ydata(:) Zdata(:)];
c1 = C(1:end-1,1);
c2 = C(1:end-1,2);
c3 = C(2:end,2);
c4 = C(2:end,1);

M1 = [M1; X];
M2 = [M2; [c1 c2 c3 c4]];
end
scatter3(0,0,0); axis tight; view(3);
nPts = nblades*numel(Xdata);
nPnls = nblades*numel(c1);
nGrps = 1;
nBCs = 1;
Empty = 0;



fid = fopen('line.neu', 'wt');
fprintf(fid,'        CONTROL INFO 2.4.6\n** GAMBIT NEUTRAL FILE\n%s\n',name);
fprintf(fid,'PROGRAM:                MATLAB     VERSION:  %s\n', version);
fprintf(fid,'%s\n', datestr(now));
fprintf(fid,'\tNUMNP\tNELEM\tNGRPS\tNBSETS\tNDFCD\tNDFVL\n');
fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\n',nPts,nPnls,nGrps,nBCs,Empty,Empty);
fprintf(fid,'ENDOFSECTION\n   NODAL COORDINATES 2.4.6\n');


str = [repmat(' ',size(M1,1),10 - size(num2str(size(M1,1)),2)) num2str([1:size(M1,1)]','%g    ') repmat('   ',size(M1,1),1) num2str(M1,'%-1.12e    ')]

for i = 1:size(M1,1)
    fprintf(fid,'%s\n',str(i,:));
%    fprintf(fid,'\t%g\t%-1.12e\t%-1.12e\t%-1.12e\n',i,M1(i,1),M1(i,2),M1(i,3));
end
fprintf(fid,'ENDOFSECTION\n      ELEMENTS/CELLS 2.4.6\n');
for i = 1:size(M2,1)
    fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',i,2,4,M2(i,1),M2(i,2),M2(i,3),M2(i,4));
end  
fprintf(fid,'ENDOFSECTION\n       ELEMENT GROUP 2.4.6\n');
fprintf(fid,'GROUP:\t\t1 ELEMENTS:\t\t%g MATERIAL:\t\t4 NFLAGS:\t\t1\n\t\t\t%s\n\t0\n',nPnls,name);

M3 = zeros(10,floor(nPnls/10));
M3(:) = 1:numel(M3);
M3 = M3';
M4 = numel(M3)+1:nPnls;

for i = 1:size(M3,1)
    
    
    
    fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',M3(i,:));
end  

for i = 1:size(M4,2)
    fprintf(fid,'\t%g',M4(i));
end
fprintf(fid,'\nENDOFSECTION\n BOUNDARY CONDITIONS 2.4.6\n\t\t\t\t%s\t1\t%g\t0\t6\n','WAKE',nPnls);
fprintf(fid,'\t%g\t2\t2\n',1:nPnls);     
fprintf(fid,'ENDOFSECTION\n');
fclose(fid);

cd ..
