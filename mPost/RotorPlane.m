function RotorPlane()
files = dir('R*.mat');
close all
s = size(files,1);

scale = 50;
val = 5
fname = files(s).name;
load(fname)


r = .4*ones(1,1000);
theta = linspace(0,2*pi,1000);
[x y] = pol2cart(theta,r);

r = linspace(.02,0.4);
th1 = 0 + 15*Time;
th2 = deg2rad(120) + 15*Time;
th3 = deg2rad(240) + 15*Time;


Dx = 1; 
Dy = 1;
Dz = 1;

figure

%   Rotor Diameter Upstream

BodyPhi(:,:,1) = DomainSliceXNormalBodyPhiMinus1;
BodyPhi(:,:,2) = DomainSliceXNormalBodyPhi1;
BodyPhi(:,:,3) = DomainSliceXNormalBodyPhiPlus1;

ProtoWakePhi(:,:,1) = DomainSliceXNormalProtoWakePhiMinus1;
ProtoWakePhi(:,:,2) = DomainSliceXNormalProtoWakePhi1;
ProtoWakePhi(:,:,3) = DomainSliceXNormalProtoWakePhiPlus1;

[dphidy dphidz dphidx] = gradient(BodyPhi + ProtoWakePhi,Dy,Dz,Dx);



subplot(3,3,1)
Vel = DomainSliceXNormalWakeVel1_x/scale;
contour(DomainSliceXNormalPosn1_y/scale,DomainSliceXNormalPosn1_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);


subplot(3,3,2)
Vel = (dphidx(:,:,2)/2)/scale;
contour(DomainSliceXNormalPosn1_y/scale,DomainSliceXNormalPosn1_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);

subplot(3,3,3)
Vel = 1+(dphidx(:,:,2)/2 + DomainSliceXNormalWakeVel1_x)/scale;
contour(DomainSliceXNormalPosn1_y/scale,DomainSliceXNormalPosn1_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);

%   Rotor Plane

BodyPhi(:,:,1) = DomainSliceXNormalBodyPhiMinus2;
BodyPhi(:,:,2) = DomainSliceXNormalBodyPhi2;
BodyPhi(:,:,3) = DomainSliceXNormalBodyPhiPlus2;

ProtoWakePhi(:,:,1) = DomainSliceXNormalProtoWakePhiMinus2;
ProtoWakePhi(:,:,2) = DomainSliceXNormalProtoWakePhi2;
ProtoWakePhi(:,:,3) = DomainSliceXNormalProtoWakePhiPlus2;

[dphidy dphidz dphidx] = gradient(BodyPhi + ProtoWakePhi,Dy,Dz,Dx);



subplot(3,3,4)
Vel = DomainSliceXNormalWakeVel2_x/scale;
contour(DomainSliceXNormalPosn2_y/scale,DomainSliceXNormalPosn2_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);


subplot(3,3,5)
Vel = (dphidx(:,:,2)/2)/scale;
contour(DomainSliceXNormalPosn2_y/scale,DomainSliceXNormalPosn2_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);

subplot(3,3,6)
Vel = 1+(dphidx(:,:,2)/2 + DomainSliceXNormalWakeVel2_x)/scale;
contour(DomainSliceXNormalPosn2_y/scale,DomainSliceXNormalPosn2_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);



%   1D DS

BodyPhi(:,:,1) = DomainSliceXNormalBodyPhiMinus3;
BodyPhi(:,:,2) = DomainSliceXNormalBodyPhi3;
BodyPhi(:,:,3) = DomainSliceXNormalBodyPhiPlus3;

ProtoWakePhi(:,:,1) = DomainSliceXNormalProtoWakePhiMinus3;
ProtoWakePhi(:,:,2) = DomainSliceXNormalProtoWakePhi3;
ProtoWakePhi(:,:,3) = DomainSliceXNormalProtoWakePhiPlus3;

[dphidy dphidz dphidx] = gradient(BodyPhi + ProtoWakePhi,Dy,Dz,Dx);



subplot(3,3,7)
Vel = DomainSliceXNormalWakeVel3_x/scale;
contour(DomainSliceXNormalPosn3_y/scale,DomainSliceXNormalPosn3_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);


subplot(3,3,8)
Vel = (dphidx(:,:,2)/2)/scale;
contour(DomainSliceXNormalPosn3_y/scale,DomainSliceXNormalPosn3_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);

subplot(3,3,9)
Vel = 1+(dphidx(:,:,2)/2 + DomainSliceXNormalWakeVel3_x)/scale;
contour(DomainSliceXNormalPosn3_y/scale,DomainSliceXNormalPosn3_z/scale,Vel,23);
axis([-.5 .5 -.5 .5]); axis square; axis square
box on
grid on
hold on
plot(x,y,'k-','LineWidth',2);
plot(r*cos(th1),r*sin(th1),'k-','LineWidth',2);
plot(r*cos(th2),r*sin(th2),'k-','LineWidth',2);
plot(r*cos(th3),r*sin(th3),'k-','LineWidth',2);


