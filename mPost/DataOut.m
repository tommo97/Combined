clear all
close all
clc

set(0,'defaulttextinterpreter','none','defaultaxesposition',[0.10    0.10    .89    .8]);
files = dir('Output.mat');


%   Find out the interpolation weights
load(files(end).name);
r = 0.9;
Rads = Rloc;
Chrd = Cloc;
R0 = Rads(BodySurface0);
C0 = Chrd(BodySurface0);


Xcp0 = CollocPts_x(BodySurface0);
Ycp0 = CollocPts_y(BodySurface0);
Zcp0 = CollocPts_z(BodySurface0);
PressChord = zeros(size(C0(1,:)));
Cp = CpHistory(1,:)';
CPress = zeros(size(Cp(1,:)));

load FDC0015.mat

chord_mlt = zeros(size(R0));

for i = 1:size(R0,2)
    Ri = R0(:,i);
    
    ind = find(Ri==r);
    
    if isempty(ind)
        ind_below = find(Ri<r);
        Rib = Ri(max(ind_below));
        
        ind_above = find(Ri>r);
        Riu = Ri(min(ind_above));
        
        ratio = (r - Rib)/(Riu - Rib);
        
        C1 = (1-ratio);
        C2 = ratio;
        chord_mlt(max(ind_below),i) = C1;
        chord_mlt(min(ind_above),i) = C2;
    else
        
        chord_mlt(ind) = 1;
    end
    
    
    
end


PressChord = sum(chord_mlt.*C0,1);


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
    
    
    for j = 1:size(CpHistoryAll,1)
        count = count + 1;
        Cp = CpHistoryAll(j,:)';
        CpD = CpHistoryAllD(j,:)';
        F = -[Area.*Cp.*Norms(:,1) Area.*Cp.*Norms(:,2) Area.*Cp.*Norms(:,3)];
        
        M = cross(F,[CollocPts_x CollocPts_y CollocPts_z]);
        Mx(j) = sum(M(:,1));
        
        CLift = sum(F(:,3))*cosd(8.5) - sum(F(:,1))*sind(8.5);
        CL(j) = sum(CLift)./4;
        CDrag = Norms(:,1).*Area.*(Cp);
        CD(j) = sum(CDrag)./4;
        
        
        
        
        Mucp0 = Mu(BodySurface0);
        CpCp0 = Cp(BodySurface0);
        CPress = sum(-chord_mlt.*CpCp0,1);
        Cl(count) = trapz(PressChord,-CPress);
        
        
        CpCp0D = CpD(BodySurface0);
        CPressD = sum(-chord_mlt.*CpCp0D,1);
        
        clf
        plot(PressChord,CPress);
        hold all
        plot(PressChord,CPressD);
        scatter(BTZSOTONpt9(:,1),BTZSOTONpt9(:,2));
        disp([i j Cl(count)])
        drawnow
    end
end
figure
% plot(PressChord,CPress);
% figure
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
a = 1;
b =  ones(100,1)/100;

plot(SubTimes(2:end),Cl(2:end))
hold all
Cl = smooth(Cl,250,'rloess');
plot(SubTimes(2:end),Cl(2:end))
figure
semilogy(0.5*(SubTimes(2:end-1) + SubTimes(3:end)),abs(diff(Cl(2:end))))
%plot(Times,CL,'-b')

%plot(Times,Cl,'-k')
%plot(Times,2*pi*AlphaHistory,'-r')