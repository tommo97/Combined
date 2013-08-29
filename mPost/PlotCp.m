%close all
clear all
clc
figure
rhoH2O = 1027;
load('Output.mat','Cp','BodyPointsX','BodyPointsY','BodyPointsZ','Rloc','Cloc','BodySurface0','CollocPts_x','CollocPts_y','CollocPts_z')
load('Output.mat','Cp','VCollocPts_x','VCollocPts_y','VCollocPts_z','Mu')
C1 = [BodyPointsX(:,1) BodyPointsY(:,1) BodyPointsZ(:,1)];
C2 = [BodyPointsX(:,2) BodyPointsY(:,2) BodyPointsZ(:,2)];
C3 = [BodyPointsX(:,3) BodyPointsY(:,3) BodyPointsZ(:,3)];
C4 = [BodyPointsX(:,4) BodyPointsY(:,4) BodyPointsZ(:,4)];
Norms = cross(C3-C1,C4-C2);
Norms = [Norms(:,1)./sqrt(dot(Norms,Norms,2)) Norms(:,2)./sqrt(dot(Norms,Norms,2)) Norms(:,3)./sqrt(dot(Norms,Norms,2))];
CollocPts = 0.25*(C1 + C2 + C3 + C4);


Rads = Rloc;
Chrd = Cloc;
R0 = Rads(BodySurface0);
C0 = Chrd(BodySurface0);


Xcp0 = CollocPts_x(BodySurface0);
Ycp0 = CollocPts_y(BodySurface0);
Zcp0 = CollocPts_z(BodySurface0);
PressChord = zeros(size(C0(1,:)));

inds = BodySurface0;
Cp = Cp(inds);
q = 0.5.*rhoH2O.*sqrt(VCollocPts_x(inds).^2 + VCollocPts_y(inds).^2 + VCollocPts_z(inds).^2);

Mucp0 = Mu(BodySurface0);
CpCp0 = Cp(BodySurface0);

r = .95;



CPress = zeros(size(Cp(1,:)));

for j = 1:size(R0,2)
    PressChord(j) = interp1(R0(:,j),C0(:,j),r,'cubic');
    CPress(j) = interp1(R0(:,j),-CpCp0(:,j),r,'cubic');
end
clf
plot(PressChord,CPress);



axis([0 1 -1 1.6])
axis square
grid on
box on

title(['C_p(x) at r/R = ' num2str(r)])
xlabel('x');
ylabel('C_p','rotation',0)

drawnow;

