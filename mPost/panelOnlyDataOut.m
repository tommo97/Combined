clear all; clear mex; clc;
close all
figure
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');

fname = 'Output.mat';
load(fname)


hold all

scale = 1;


C1 = [BodyPointsX(:,1) BodyPointsY(:,1) BodyPointsZ(:,1)]/scale;
C2 = [BodyPointsX(:,2) BodyPointsY(:,2) BodyPointsZ(:,2)]/scale;
C3 = [BodyPointsX(:,3) BodyPointsY(:,3) BodyPointsZ(:,3)]/scale;
C4 = [BodyPointsX(:,4) BodyPointsY(:,4) BodyPointsZ(:,4)]/scale;
BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);


p1 = patch('Vertices',BodyPanPts,...
    'Faces',PtIDS,'FaceVertexCData',-Cp(:),...
    'FaceColor','flat','EdgeColor','k');

BodyPointsX = BodyPointsX;
C1 = [WakePanC1_x WakePanC1_y WakePanC1_z]/scale;
C2 = [WakePanC2_x WakePanC2_y WakePanC2_z]/scale;
C3 = [WakePanC3_x WakePanC3_y WakePanC3_z]/scale;
C4 = [WakePanC4_x WakePanC4_y WakePanC4_z]/scale;
WakePanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);


p2 = patch('Vertices',WakePanPts,...
    'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
    'FaceColor','flat','EdgeColor','k');

box on
grid on
view(3)

axis equal tight;% lighting phong; camlight right;

set(gca,'Projection','perspective')


