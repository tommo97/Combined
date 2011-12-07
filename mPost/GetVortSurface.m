clear all; clear mex; clc;
files = dir('R*.mat');
close all
s = size(files,1);

scale = 50;
val = 5
fname = files(s).name;
load(fname)

XCG = scale * Time * 1.0;
THETA = Time * BodyRates0_x;

[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);
YI = YI + XCG/scale;
p1 = patch(isosurface(XI,YI,ZI,VI,val));
isonormals(XI,YI,ZI,VI,p1);
isocolors(XI,YI,ZI,VIx,p1);
set(p1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')


isocolors(XI,YI,ZI,VIx,p1);


set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
hold all

return
drawnow


n = 50;


y = min(Domain(:,2)):1:max(Domain(:,2));
z = min(Domain(:,3)):1:max(Domain(:,3));
[rotorplaney rotorplanez] = meshgrid(y,z);
rotorplaneu = zeros(size(rotorplaney));
rotorplanev = zeros(size(rotorplaney));
rotorplanew = zeros(size(rotorplaney));
DEL2 = 0.25;

for i = 1:numel(rotorplaney)
    
    y = rotorplaney(i);
    x = Time * -1 * scale;
    z = rotorplanez(i);
    
    dx = x - Domain(:,1);
    dy = y - Domain(:,2);
    dz = z - Domain(:,3);
    
    
    
    
    nrm = sqrt(DEL2 + dx.*dx + dy.*dy + dz.*dz);
    
    mult = 1 ./ (4 * pi .* nrm .* nrm .* nrm);
    
    omx = Domain(:,4);
    omy = Domain(:,5);
    omz = Domain(:,6);
    
    Vx = mult.*(omy.*dz - omz.*dy);
    Vy = mult.*(omz.*dx - omx.*dz);
    Vz = mult.*(omx.*dy - omy.*dx);
    
    rotorplaneu(i) = rotorplaneu(i) + sum(Vx);
    rotorplanev(i) = rotorplanev(i) + sum(Vy);
    rotorplanew(i) = rotorplanew(i) + sum(Vz);
    
    disp(i)
end


figure
[C,h] = contour(rotorplaney/scale,rotorplanez/scale,(rotorplanew/scale));
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
axis equal tight
grid on
box on
xlabel('$y$[m]')
yl = ylabel('$z$[m]')
set(yl,'Rotation',0)
mlf2pdf(gcf,'ZVel64x64jet')



ViewParams = GetViewParams;
return

load(fname)

BodyPointsX = BodyPointsX + XCG;
C1 = [BodyPointsY(:,1) BodyPointsX(:,1) BodyPointsZ(:,1)]/scale;
C2 = [BodyPointsY(:,2) BodyPointsX(:,2) BodyPointsZ(:,2)]/scale;
C3 = [BodyPointsY(:,3) BodyPointsX(:,3) BodyPointsZ(:,3)]/scale;
C4 = [BodyPointsY(:,4) BodyPointsX(:,4) BodyPointsZ(:,4)]/scale;
BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);


p2 = patch('Vertices',BodyPanPts,...
    'Faces',PtIDS,'FaceVertexCData',-Cp(:),...
    'FaceColor','flat','EdgeColor','none');

box on
grid on
view(3)

axis equal tight; lighting phong; camlight right;

set(gca,'Projection','perspective')






return







WakePanPts = [  WakePanC1_y WakePanC1_x WakePanC1_z;
    WakePanC2_y WakePanC2_x WakePanC2_z;
    WakePanC3_y WakePanC3_x WakePanC3_z;
    WakePanC4_y WakePanC4_x WakePanC4_z]/scale;

PtIDS = 1:4*length(WakePanC1_x);
PtIDS = reshape(PtIDS,length(WakePanC1_x),4);

p3 = patch('Vertices',WakePanPts,...
    'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
    'FaceColor','flat','EdgeColor','k');

cd3 = get(p3,'CData');
cax3 = [min(WakePanGamma(:)) max(WakePanGamma(:))];




i1 = 1;
cmap = [cmapX];
i2 = length(cmap);
i3 = i2 + 1;
cmap = [cmap;cmapY];
i4 = length(cmap);
i5 = i4 + 1;
cmap = [cmap;cmapZ];
i6 = length(cmap);

colormap(cmap)
set(p1,'CData',i1:1:i2);
%set(p2,'CData',i3:1:i4);
%set(p3,'CData',i5:1:i6);
%scatter3(VortonX_y/scale,VortonX_x/scale,VortonX_z/scale,'.');
%
axis equal tight
view(3)


VIs = VI;
figure
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
VIs(VIs>5*std(VIs(:))) = 5*std(VIs(:));
slice(XI,YI,ZI,VIs,[0],[ ],[])
axis equal
shading flat
box on
view(3)

axis equal tight
set(gca,'Projection','perspective')