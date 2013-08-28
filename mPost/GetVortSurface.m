clear all; clear mex; clc;
files = dir('RunData*8.mat');
close all
figure
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale','Times','BodyRates0_x','BodyPointsX','BodyPointsY','BodyPointsZ')
load(fname,'Time','XCG')

scale = GambitScale;
val = 10;

XCG = 0.0;%scale * Time * 10.0;
THETA = Time * BodyRates0_x;

[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);

%[curlx,curly,curlz,cav] = curl(y,x,z,V,U,W) ;
%OM = sqrt(curlx.^2 + curly.^2 + curlz.^2);
YIMoved = YI + XCG/scale;
p1 = patch(isosurface(XI,YIMoved ,ZI,VI,val));
isonormals(XI,YIMoved ,ZI,VI,p1);
%isocolors(XI,YI,ZI,VIx,p1);
set(p1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')


%isocolors(XI,YI,ZI,VIx,p1);


hold all




ViewParams = GetViewParams;

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
    'FaceColor','flat','EdgeColor','k');

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

[sx sy sz] = meshgrid(0,-0.4:0.1:.4,-0.4:.1:0.4);
x = x - min(x(:));
h = streamline(y/scale,x/scale,z/scale,V/scale,1+U/scale,W/scale,sx,sy,sz);
%h = streamline(x,y,z,u,v,w,sx,sy,sz);
set(h,'Color','red')
VIs = VI;
figure
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
VIs(VIs>5*std(VIs(:))) = 5*std(VIs(:));
slice(XI,YI,ZI,VIs,[0],-80+[10 20 30 40 50 60],[])
axis equal
shading flat
box on
view(3)

axis equal tight
set(gca,'Projection','perspective')