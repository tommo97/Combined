clear all; clear mex; clc;
files = dir('RunData*00.mat');
close all
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale')


scale = GambitScale;
val = 100;

col{1} = 'r';
col{2} = 'b';

[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);

% hold all
set(gcf,'renderer','opengl')
% p1 = patch(isosurface(XI,YI,ZI,VI,val/3));
% isonormals(XI,YI,ZI,VI,p1);
% set(p1,'FaceColor',[1 0.8 1],'EdgeColor','none','faceAlpha',.2)
% % 
% p2 = patch(isosurface(XI,YI,ZI,VI,2*val/3));
% isonormals(XI,YI,ZI,VI,p2);
% set(p2,'FaceColor',[1 0.8 1],'EdgeColor','none','faceAlpha',.2)
hold all
for i = 1:size(VIx,2)
    vi = sqrt(VIx{i}.^2 + VIy{i}.^2 + VIz{i}.^2);
    p{i} = patch(isosurface(XI,YI,ZI,vi,val));
    isonormals(XI,YI,ZI,vi,p{i});
    set(p{i},'FaceColor',col{i},'EdgeColor','none');%,'faceAlpha',.2)
end
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');





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