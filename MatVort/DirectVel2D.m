function [V1 V2] = DirectVel2D(dx, omega1, omega2)
scale = 10;
DELa = 0.;
DELb = scale*DELa;
dx1 = dx;
dx2 = scale*dx;
nrm1 = sqrt(DELa + dx1);
nrm2 = sqrt(DELb + dx2); 
mult1 = -1 / (4 * pi * nrm1 * nrm1 * nrm1);
malt = -1/(4*pi*dx*dx);
mult2 = -1 / (4 * pi * nrm2 * nrm2 * nrm2);



V1 = mult1 * dx1*omega1;
V2 = mult2 * dx2*omega2;
V3 = malt * sign(dx)*omega1
