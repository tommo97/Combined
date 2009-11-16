%function test_ends()
%%  Aerofoil
close all; clear all;
n = 25;
x = cumtrapz(1 - cos(linspace(0,pi,n)).^2)/max(cumtrapz(1 - cos(linspace(0,pi,n)).^2));
x = linspace(0,1,n);
[Aerofoil z] = NRELFoil(x);
d = cos(linspace(0,pi,n-1));

zu = z.S809.US;
zl = z.S809.LS;

%   first attempt
xxx = [x;[x(1) x(2:end)+linspace(1,0,n-1).*d.*diff(x)];x];
zzz = [zu;.5*(zu+zl);zl];

surf(xxx,zzz,zeros(size(xxx)))
view(2)
axis equal
