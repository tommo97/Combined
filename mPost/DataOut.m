%function DataOut()


clear all
clc
%close all
files = dir('R*.mat');

count = 0;
for i = 1:length(files)
    load(files(i).name);
    
    
    
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
    
    
    for j = 1:size(CpHistoryAll,1)
        count = count + 1;
        Cp = CpHistoryAll(j,:)';
        F = -[Area.*Cp.*Norms(:,1) Area.*Cp.*Norms(:,2) Area.*Cp.*Norms(:,3)];
        
        M = cross(F,[CollocPts_x CollocPts_y CollocPts_z]);
        Mx(j) = sum(M(:,1));
        
        CLift = sum(F(:,3))*cosd(8.5) - sum(F(:,1))*sind(8.5);
        CL(j) = sum(CLift)./12;
        CDrag = Norms(:,1).*Area.*(Cp);
        CD(j) = sum(CDrag)./12;
        
        
        
        
        Mucp0 = Mu(BodySurface0);
        CpCp0 = Cp(BodySurface0);
        
        r = 0;
        
        
        
        CPress = zeros(size(Cp(1,:)));
        
        for k = 1:size(R0,2)
            PressChord(k) = interp1(R0(:,k),C0(:,k),r,'cubic');
            CPress(k) = interp1(R0(:,k),-CpCp0(:,k),r,'cubic');
        end
       
        
        Cl(count) = trapz(PressChord,-CPress);
        %clf
        %plot(PressChord,CPress);
        disp([i j Cl(count)])
        %drawnow
    end
end
% clf
% 
% 
% BodyPanPts = [C1;C2;C3;C4];
% PtIDS = 1:4*length(C1);
% PtIDS = reshape(PtIDS,length(C1),4);
% 
% p = patch('Vertices',BodyPanPts,...
%     'Faces',PtIDS,'FaceVertexCData',Cp(:),...
%     'FaceColor','flat','EdgeColor','none');
% set(gcf,'Renderer','OpenGL')
% hold all
% 
% WakePanPts = [WakePanC1_x WakePanC1_y WakePanC1_z;
%     WakePanC2_x WakePanC2_y WakePanC2_z;
%     WakePanC3_x WakePanC3_y WakePanC3_z;
%     WakePanC4_x WakePanC4_y WakePanC4_z];
% 
% PtIDS = 1:4*length(WakePanC1_x);
% PtIDS = reshape(PtIDS,length(WakePanC1_x),4);
% 
% p2 = patch('Vertices',WakePanPts,...
%     'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
%     'FaceColor','flat','EdgeColor','k');
% 
% 
% % ProtoWakePanPts = [ProtoWakePointsX(:,1) ProtoWakePointsY(:,1) ProtoWakePointsZ(:,1);
% %     ProtoWakePointsX(:,2) ProtoWakePointsY(:,2) ProtoWakePointsZ(:,2);
% %     ProtoWakePointsX(:,3) ProtoWakePointsY(:,3) ProtoWakePointsZ(:,3);
% %     ProtoWakePointsX(:,4) ProtoWakePointsY(:,4) ProtoWakePointsZ(:,4)];
% % 
% % 
% % 
% % 
% % PtIDS = 1:4*length(ProtoWakePointsX(:,1));
% % PtIDS = reshape(PtIDS,length(ProtoWakePointsX(:,1)),4);
% % 
% % p3 = patch('Vertices',ProtoWakePanPts,...
% %    'Faces',PtIDS,'FaceVertexCData',ProtoWakeGamma(:),...
% %    'FaceColor','flat','EdgeColor','k');
% 
% axis equal tight
% view(3)
% 


%scatter3(VortonX_x(:),VortonX_y(:),VortonX_z(:),'+')
% 
%figure
plot(SubTimes(2:end),Cl(2:end))
%plot(Times,CL,'-b')
hold all
%plot(Times,Cl,'-k')
%plot(Times,2*pi*AlphaHistory,'-r')