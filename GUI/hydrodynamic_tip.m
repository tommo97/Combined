function hydrodynamic_tip()
cla
load ../GeomInputMatFiles/DMTP4119Approximation.mat
h = 1e-3;
Radius(end) = 0.995;
Chord(end) = 0.1;
Radius = [Radius 1];
Chord = [Chord 0];


N = 100;

r = linspace(min(Radius),max(Radius),N);

c = interp1(Radius,Chord,r,'cubic');


%   Get gradient d(c/2)/dr at r=0.97R
r1 = min(r);
r2 = 0.97;
c1 = -0.0625; 
ds_dr_1 = 0.;

spans =  linspace(0,1,10);


plot(r,c/2);
hold on

plot(r,-c/2);

r_ = linspace(r1,r2,N);

[c_, c2, dc_dr, coeffts] = getLineDriver(Chord, Radius,r_, r1, r2, c1, ds_dr_1, h);

s_ = [0 cumsum(sqrt(diff(c_).^2 + diff(r_).^2))]; 
s_ = s_/s_(end);

rsp_ = interp1(s_,r_,spans,'cubic');
csp_ = interp1(s_,c_,spans,'cubic');
scatter(rsp_,csp_)

plot(r_,c_);


%   Now need to cast a ray along a radius perpendicular to the line r_,c_
%   until it meets the upper surface








s_full = [0 cumsum(sqrt(diff(c).^2 + diff(r).^2))]; s_full_max = s_full(end);
s_full = s_full/s_full_max;



su = interp1(r,s_full,linspace(r1,r2, N),'cubic');


rsu = interp1(s_full,r,su,'cubic');
csu = interp1(s_full,c,su,'cubic');

rsup = interp1(s_full,r,max(su)*spans,'cubic');
csup = interp1(s_full,c,max(su)*spans,'cubic');





r0 = rsup(3);
c0 = csup(3)/2;


gradnt = 3*coeffts(1)*r.^2 + 2*coeffts(2)*r + coeffts(3);

dydx = (c0 - c)./(r0 - r);

figure

plot(r,dydx);
hold all
plot(r,-1./gradnt);
ax = axis;
axis([ax(1) ax(2) -10 10])
return




















global AA BB CC DD X1 Y1
AA = coeffts(1); BB = coeffts(2); CC = coeffts(3); DD = coeffts(4);


x0 = X1;           % Make a starting guess at the solution
options=optimset('Display','iter');   % Option to display output
[x,fval] = fsolve(@myfun,r2,options)  % Call solver


[x2 res] = lsqnonlin(@myfun,r2,r1,r2);










%[x1 count] = getRay(rsup(3),csup(3)/2,coeffts);

plot([rsup(3) x2],[csup(3)/2 AA.*x2.*x2.*x2 + BB.*x2.*x2 + CC.*x2 + DD],'r')


y2 = AA.*x2.*x2.*x2 + BB.*x2.*x2 + CC.*x2 + DD;
m1 = (y2 - Y1)./(x2 - X1 + 1e-16);
m2 = 3*AA.*x2.*x2 + 2*BB.*x2 + CC;
disp(num2str([m1 m2]))




plot(rsu,csu/2,'r')
scatter(rsup,csup/2,'s')



cl = -interp1(r,c,linspace(r1,r2,N))/2;

su = [0 sqrt(diff(c).^2 + diff(r).^2)]; su = su/max(su);





grads = 3*coeffts(1)*spans.^2 + 2*coeffts(2)*spans + coeffts(3);



scatter(r2,c2);

plot(r2 + 5*[-h h], c2 + 5*[-h h]*dc_dr)
plot(r2 + 5*[-h h], c2 + -5*[-h h]/dc_dr)
plot(r1 + 5*[-h h], c1 + 5*[-h h]*ds_dr_1)




axis equal



function F = myfun(x2)
global AA BB CC DD X1 Y1
y2 = AA.*x2.*x2.*x2 + BB.*x2.*x2 + CC.*x2 + DD;
m1 = (y2 - Y1)./(x2 - X1 + 1e-16);
m2 = 3*AA.*x2.*x2 + 2*BB.*x2 + CC;
F = abs(m2 - m1);

  
  
function [x1 count] = getRay(x0,y0,coeffts)
a = coeffts(1); b = coeffts(2); c = coeffts(3); d = coeffts(4);
x1 = x0;
eps = 1;
count = 0;
while eps > 1e-6
    val = d - y0 - (x0 - x1)/(3*a*x1^2 + 2*b*x1 + c) + c*x1 + a*x1^3 + b*x1^2;
    deriv = c + 2*b*x1 + 3*a*x1^2 + 1/(3*a*x1^2 + 2*b*x1 + c) + ((2*b + 6*a*x1)*(x0 - x1))/(3*a*x1^2 + 2*b*x1 + c)^2;
    
    new_val = val - (val/deriv);
    eps = abs(val - new_val);
    x1 = abs(new_val);
    
    disp(num2str([x1 count]));
    count = count + 1;
end




function [c, c3, dc_dr, coeffts] = getLineDriver(Chord, Radius,r , x1, x2, y1, dy_dx_1, h)
y2 = interp1(Radius,Chord,x2,'cubic')/2;
[dc_dr, c3] = getGradient(Radius,Chord, x2,h);



%   Gradient of new sweep line @ r/R = rPt is thus -1/dc_dr
ds_dr_2 = -1/dc_dr;
dy_dx_2 = ds_dr_2;
%   Gradient of new sweep line @ r/R = Rmin/R is 0


%   Two points and two gradients - can thus specify a quadratic line
[c, coeffts] = getLine(x1, y1, x2, y2, dy_dx_1, dy_dx_2, r);


function [c, coeffts] = getLine(x1, y1, x2, y2, dy_dx_1, dy_dx_2, r)
A = [x1^3 x1^2 x1 1; x2^3 x2^2 x2 1; 3*x1^2 2*x1 1 0; 3*x2^2 2*x2 1 0];
b = [y1;y2;dy_dx_1;dy_dx_2];
coeffts = A\b;
c = coeffts(1)*r.^3 + coeffts(2)*r.^2 + coeffts(3)*r + coeffts(4);

function [dc_dr, c3] = getGradient(Radius,Chord, rPt, h)

r1 = rPt - 2*h;
r2 = rPt - h;
r3 = rPt;
r4 = rPt + h;
r5 = rPt + 2*h;

c1 = interp1(Radius,Chord,r1,'cubic')/2;
c2 = interp1(Radius,Chord,r2,'cubic')/2;
c3 = interp1(Radius,Chord,r3,'cubic')/2;
c4 = interp1(Radius,Chord,r4,'cubic')/2;
c5 = interp1(Radius,Chord,r5,'cubic')/2;

r21 = r2-r1; r32 = r3-r2; r43 = r4-r3; r54 = r5-r4;
r31 = r3-r1; r42 = r4-r2; r53 = r5-r3;
r41 = r4-r1; r52 = r5-r2; r51 = r5-r1;
w21 = -r32*r43*r53/r31/r41/r51;
w32 = (r31*r41*(r31+r52)+r32*r52*(r31+r42))*r43*r53/r31/r42/r41/r52/r51;
w43 = (r53*r52*(r53+r41)+r43*r41*(r53+r42))*r32*r31/r53/r42/r41/r52/r51;
w54 = -r32*r43*r31/r53/r52/r51;
dc_dr = w21*(c2-c1)/r21+w32*(c3-c2)/r32+w43*(c4-c3)/r43+w54*(c5-c4)/r54;