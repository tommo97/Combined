clear all
close all

clc
load Output.mat

C1 = [BodyPointsX(:,1) BodyPointsY(:,1) BodyPointsZ(:,1)];
C2 = [BodyPointsX(:,2) BodyPointsY(:,2) BodyPointsZ(:,2)];
C3 = [BodyPointsX(:,3) BodyPointsY(:,3) BodyPointsZ(:,3)];
C4 = [BodyPointsX(:,4) BodyPointsY(:,4) BodyPointsZ(:,4)];
figure
BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);
p = patch('Vertices',BodyPanPts,...
'Faces',PtIDS,'FaceVertexCData',Cp,...
'FaceColor','flat','EdgeColor','none');
set(gcf,'Renderer','OpenGL')
hold all

view(3)
scatter3(BodyCG0_x,BodyCG0_y,BodyCG0_z)



WakePanPts = [WakePanC1_x WakePanC1_y WakePanC1_z;
WakePanC2_x WakePanC2_y WakePanC2_z;
WakePanC3_x WakePanC3_y WakePanC3_z;
WakePanC4_x WakePanC4_y WakePanC4_z];
PtIDS = 1:4*length(WakePanC1_x);
PtIDS = reshape(PtIDS,length(WakePanC1_x),4);
p2 = patch('Vertices',WakePanPts,...
'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
'FaceColor','flat','EdgeColor','k');


axis equal tight