clear all; clear mex; clc;
files = dir('RunData*2.mat');
close all
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale')


scale = GambitScale;
val = 30;



[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);

% hold all
% set(gcf,'renderer','opengl')
% p1 = patch(isosurface(XI,YI,ZI,VI,val/3));
% isonormals(XI,YI,ZI,VI,p1);
% set(p1,'FaceColor',[1 0.8 1],'EdgeColor','none','faceAlpha',.2)
% % 
% p2 = patch(isosurface(XI,YI,ZI,VI,2*val/3));
% isonormals(XI,YI,ZI,VI,p2);
% set(p2,'FaceColor',[1 0.8 1],'EdgeColor','none','faceAlpha',.2)

p3 = patch(isosurface(XI,YI,ZI,VI,val));
isonormals(XI,YI,ZI,VI,p3);
set(p3,'FaceColor',[1 0.8 1],'EdgeColor','none');%,'faceAlpha',.2)

set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
hold all




ViewParams = GetViewParams;


box on
grid on
view(3)

axis equal tight; lighting phong; camlight right;
axis([-1.5 1.5 -1.5 1.5 -2.5 0.5])
set(gca,'Projection','perspective')
set(gca,'fontsize',14)
return
figure

slice(XI,YI,ZI,VI,[0],[0],[0])
axis equal
shading flat
box off
axis off
grid off
view(3)

axis equal tight
set(gca,'Projection','perspective')