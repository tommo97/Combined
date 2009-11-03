% Lifting Line code in MATLAB
% Coded by L. sankar, December 1998
%
croot= 0.737;%input('Please input root chord   ');
ctip= 0.355;%input('Please input tip chord   ');
span= 5.25;%input('Please input span   ');
thetaroot= 20.04;%input('Please input Root twist angle in degrees   ');
thetatip= -1.815;%input('Please input tip twist angle in degrees   ');
a0root= 2*pi;%input('Please input root lift curve slope in units/ radian   ');
a0tip= 2*pi;%input('Please input lift curve slope at the tip, in units/radian   ');
alpha= 0;%input('Please input angle of attack, in degrees   ');
alpha0root= 0;%input('Please input zero-lift angle at the root   ');
alpha0tip= 0;%input('Please input zero lift angle at the tip   ');

croot= 1;%input('Please input root chord   ');
ctip= 1;%input('Please input tip chord   ');
span= 10;%input('Please input span   ');
thetaroot= 12;%20.04;%input('Please input Root twist angle in degrees   ');
thetatip= 12;%1.815;%input('Please input tip twist angle in degrees   ');

thetaroot=deg2rad(thetaroot);
thetatip=deg2rad(thetatip);
alpha = deg2rad(alpha);
alpha0root= deg2rad(alpha0root);
alpha0tip= deg2rad(alpha0tip);
n = 100;
theta=zeros(1,n);
y=zeros(1,n);
c=zeros(1,n);
cl=zeros(1,n);
alp=zeros(1,n);
rhs=zeros(n,1);
b=zeros(n,n);
a=zeros(n,1);
%
% Define properties at n span stations
%
for i=1:n
    theta(i) = (i) * pi/(2. * n);
    y(i) = span * 0.5 * cos(theta(i));
    c(i) = croot+(ctip-croot)*y(i)*2./span;
    alp(i) = alpha+thetaroot-(alpha0root+(alpha0tip-alpha0root+thetaroot-thetatip)*y(i)*2./span);
    a(i) = a0root+(a0tip-a0root)*y(i)*2./span;
end

% Set up 2n x 2n system of equations for A1, A3 , ... A2n-1
for j=1:n
    mu = c(j)* a(j) / (4. * span);
    rhs(j,1)=alp(j)*sin(theta(j))*c(j)*a(j)/(4*span);
    for i=1:n
        l = 2 * i-1;
        b(j,i)=sin(l*theta(j))*(mu*l+sin(theta(j)));
    end
end
%
% Solve for the Fourier Coefficients
%
a=b\rhs;
% Compute wing area and aspect ratio
S=(croot+ctip)/2.*span;
AR=span*span/S;
% Compute CL
CL=a(1)*pi*AR;
% Compute CD
CD0=1.;
for i=2:n
    CD0=CD0+(2.*i-1)*a(i)*a(i)/(a(1)*a(1));
end
CD = pi * AR * a(1) * a(1) * CD0;
% Compute spanwise load distribution, normalized by Freestream velocity times Span
gamma=zeros(1,n);
for i=1:n
    gamma(i)=0.0;
    y(i)=y(i)*2./span;
    for j=1:n
        gamma(i)=gamma(i)+2.*a(j)*sin((2*j-1)*theta(i));
    end
    gamma(i) = gamma(i) *  span;
    cl(i) = gamma(i)/(0.5*c(i));
end
cl
%
% Plot the circulation distribution from root to tip with labels
%
plot(y,gamma);
xlabel('y/Semi-Span');
ylabel('Cl');
title('Sectional Lift Coefficient distribution');
