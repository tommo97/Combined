clear all
clc
load RunData_000001.mat
hold all
c = 'b';
s = 'o';
scale = GambitScale;

OMEGA = [BodyRates0_x BodyRates0_y BodyRates0_z];
VEL = [1 0 0] * scale;

R = 0.4 * scale;

denom = 0.5*1027*VEL(1)*VEL(1)*2*pi*R*R;

CT = 2*Force_x/denom;

CP = -2*BodyRates0_x*Torque_x/(VEL(1)*denom);



plot(Times,(smooth(CP,25)),[c '-']);

scatter(0:0.1:max(Times),interp1(Times,smooth(CP,25),0:0.1:max(Times)),[c s],'filled')

