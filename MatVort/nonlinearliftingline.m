%%  Non linear lifting line
%clf;
clear all;
clc;

%%  Turbine stats
RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
    2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];

CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
    0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
    3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];

croot= 1;%input('Please input root chord   ');
ctip= 1;%input('Please input tip chord   ');
span= 10;%input('Please input span   ');
thetaroot= 12;%20.04;%input('Please input Root twist angle in degrees   ');
thetatip= 12;%1.815;%input('Please input tip twist angle in degrees   ');

%%  Define some points

%   Span lies along y axis, chord parallel to x axis
n = 101;
cutout = 0;%1.257;




RADIUS = linspace(0,span,n);

CHORD = linspace(croot,ctip,n);

THETA = linspace(thetaroot,thetatip,n);

radius = max(RADIUS);
th0 = 0;% (3-THETA(end));
pitch_axis = 0.3; % chord

rads = linspace(cutout,radius,n)';
chord = interp1(RADIUS,CHORD,rads,'cubic');
LEx = pitch_axis*chord;
TEx = -(1-pitch_axis)*chord;
th = th0 + interp1(RADIUS,THETA,rads,'cubic');

xdata0 = [LEx TEx];
ydata0 = [rads rads];
zdata0 = zeros(size(ydata0));

th = [th th];
xdata = xdata0.*cosd(th) + zdata0.*sind(th);
ydata = ydata0;
zdata = -xdata0.*sind(th) + zdata0.*cosd(th);
C = xdata;
C(:) = (1:numel(xdata));


X = [xdata(:) ydata(:) zdata(:)];




c1 = C(1:end-1,1);
c2 = C(1:end-1,2);
c3 = C(2:end,2);
c4 = C(2:end,1);

x1 = X(c1,:);
x2 = X(c2,:);
x3 = X(c3,:);
x4 = X(c4,:);

CP = .5*(.5*(x1+x4) + .25*(x1 + x2 + x3 + x4));

R = sqrt(dot(CP,CP,2));
C = interp1(rads,chord,R,'cubic','extrap');
alpha_g = interp1(rads,th(:,1),R,'cubic','extrap');

dr = R(2) - R(1);


localx =  .5*(x1+x4) - .5*(x2+x3);
localx = [localx(:,1)./sqrt(dot(localx,localx,2)) localx(:,2)./sqrt(dot(localx,localx,2)) localx(:,3)./sqrt(dot(localx,localx,2))];
localz = cross(x4-x2,x3-x1);
localz = [localz(:,1)./(sqrt(dot(localz,localz,2))) localz(:,2)./(sqrt(dot(localz,localz,2))) localz(:,3)./(sqrt(dot(localz,localz,2)))];
localy = cross(localz,localx);
% hold all
% scatter3(CP(:,1),CP(:,2),CP(:,3));
% quiver3(CP(:,1),CP(:,2),CP(:,3),localz(:,1),localz(:,2),localz(:,3));
% quiver3(CP(:,1),CP(:,2),CP(:,3),localx(:,1),localx(:,2),localx(:,3));
% quiver3(CP(:,1),CP(:,2),CP(:,3),localy(:,1),localy(:,2),localy(:,3));
% surf(xdata,ydata,zdata);
% axis equal
%%  Assume a lift distribution

uinf = [-1 0 0];

Uinf = ones(size(localx));
Uinf(:,1) = uinf(1);
Uinf(:,2) = uinf(2);
Uinf(:,3) = uinf(3);
Vxlocal = dot(Uinf,localx,2);
Vzlocal = dot(Uinf,localz,2);

alpha_i = -atand(Vzlocal./Vxlocal);
alpha_e = alpha_i;

[Cl Cd] = SectionCoefftsBahajPostStall(alpha_e);


isdone = false;
count = 1;
maxits = 1000;
res = ones(5,1);

uinfmag = sqrt(dot(uinf,uinf));
mult = 1/(4*pi*uinfmag);
gamma_input = 0.5*uinfmag.*C.*Cl;
gamma_input =  sin(linspace(0,pi,n-1))' + .1;
while ~isdone
    %%  Calculate a new induced alpha
    %   Simpson quadrature
   
    abcissa = R;
   
    
    %f = gradient(gamma_input,dr);
    
    temp_g = interp1(R,gamma_input,RADIUS,'cubic',0);
    f = diff(temp_g)./diff(RADIUS);
    f = f';
    
    alpha_i = zeros(size(alpha_g));
    for i = 1:length(R);
        
        DR = R(i) - R;
        fmod = f./DR;
        %%  remove singularity
        if i==1
            fmod(i) = fmod(i+1);
        elseif i==length(R)
            fmod(i) = fmod(i-1);
        else
            fmod(i) = 0.5*(fmod(i+1) + fmod(i-1));
        end
        
        alpha_i(i) =  mult*Simpson(fmod,dr);
    end
    %%  Recalculate gamma
    alpha_e =  alpha_g - rad2deg(alpha_i);
    [Cl Cd] = SectionCoefftsBahajPostStall(alpha_e);
    gamma_new = 0.5*uinfmag.*C.*Cl;
    
    
    if count <= 5
        res(count) = max(abs(gamma_new - gamma_input));
    else
        res(1:end-1) = res(2:end);
        res(end) = max(abs(gamma_new - gamma_input));
    end

    
    
    isdone = (max(res) < 1e-3) || (count > maxits);
    
    count = count + 1;
    
    gamma_input = gamma_input + 0.01*(gamma_new - gamma_input);
    disp([count res(end)]);
end

hold all
plot(2*(.5-R./max(R)),Cl);
axis([0 1 0 2])



