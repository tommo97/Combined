clear all
close all
clc
load Output.mat


% 
% for i = 1:size(CollocPtHistory0_x,2)
%     
% %x1 = CollocPtHistory1_x(i,:);
% %z1 = CollocPtHistory1_z(i,:);
% %y1 = CollocPtHistory1_y(i,:);
% 
% x0 = CollocPtHistory0_x(:,i);
% z0 = CollocPtHistory0_z(:,i);
% y0 = CollocPtHistory0_y(:,i);
% 
% 
% %scatter3(x1,y1,z1);
% hold all
% scatter3(x0,y0,z0);
% 
% axis equal
% 
% view(3)
% 
% pause(0.1)
% 
% end








%return
C1 = [BodyPointsX(:,1) BodyPointsY(:,1) BodyPointsZ(:,1)];
C2 = [BodyPointsX(:,2) BodyPointsY(:,2) BodyPointsZ(:,2)];
C3 = [BodyPointsX(:,3) BodyPointsY(:,3) BodyPointsZ(:,3)];
C4 = [BodyPointsX(:,4) BodyPointsY(:,4) BodyPointsZ(:,4)];
figure
BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);
p = patch('Vertices',BodyPanPts,...
'Faces',PtIDS,'FaceVertexCData',1,...
'FaceColor','flat','EdgeColor','none');
set(gcf,'Renderer','OpenGL')
hold all
WakePanPts = [WakePanC1_x WakePanC1_y WakePanC1_z;
WakePanC2_x WakePanC2_y WakePanC2_z;
WakePanC3_x WakePanC3_y WakePanC3_z;
WakePanC4_x WakePanC4_y WakePanC4_z];
PtIDS = 1:4*length(WakePanC1_x);
PtIDS = reshape(PtIDS,length(WakePanC1_x),4);
p2 = patch('Vertices',WakePanPts,...
'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
'FaceColor','flat','EdgeColor','k');
% ProtoWakePanPts = [ProtoWakePointsX(:,1) ProtoWakePointsY(:,1) ProtoWakePointsZ(:,1);
%     ProtoWakePointsX(:,2) ProtoWakePointsY(:,2) ProtoWakePointsZ(:,2);
%     ProtoWakePointsX(:,3) ProtoWakePointsY(:,3) ProtoWakePointsZ(:,3);
%     ProtoWakePointsX(:,4) ProtoWakePointsY(:,4) ProtoWakePointsZ(:,4)];
%
% PtIDS = 1:4*length(ProtoWakePointsX(:,1));
% PtIDS = reshape(PtIDS,length(ProtoWakePointsX(:,1)),4);
%
% p3 = patch('Vertices',ProtoWakePanPts,...
%    'Faces',PtIDS,'FaceVertexCData',ProtoWakeGamma(:),...
%    'FaceColor','flat','EdgeColor','k');
axis equal tight
scatter3(BodyCG0_x,BodyCG0_y,BodyCG0_z)
view(3)