clear all; clear mex; close all; clc;
files = dir('R*.mat');

s = size(files,1);
% while s < 100
%     files = dir('f*.dat');
%     s = size(files,1);
%     disp(s)
%     pause(60);
% end

scale = 1;
val = 20;

fname = files(s).name;

[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);

p1 = patch(isosurface(XI,YI,ZI,VI,val));
isonormals(XI,YI,ZI,VI,p1);
isocolors(XI,YI,ZI,VIx,p1);
set(p1,'FaceColor','interp','EdgeColor','none')
cax1 = caxis;

VIs = VIx;

VIs = VIs - cax1(1);
VIs = 1 + 31*VIs/(cax1(2) - cax1(1));
VIs = floor(VIs);

isocolors(XI,YI,ZI,VIs,p1);


set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
hold all
%scatter3([0 0 0 0],[0 0 410 410],[-20 20 -20 20],[1 1 1 1],'white');
%view([0 90]);
%view(3);
axis equal tight;

%set(gca,'CameraUpVector',[-1 0 0])
axis equal tight; lighting phong; camlight right;

%zoom(1.8);
drawnow
ViewParams = GetViewParams;
%camlight right
%view(38.5,16);

makemovie = false;
if makemovie
    aviobj = avifile(['temp.avi'],'fps',10);
    
    for i = 1:s
        clf;
        fname = files(i).name;
        %     extract(fname,val);
        %     campos = ViewParams.campos;
        %     camtarget = ViewParams.camtarget;
        %     camva = ViewParams.camva;
        %     camproj = ViewParams.camproj;
        extract(fname,val);
        set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
        hold all
        %scatter3([0 0 0 0],[0 0 410 410],[-20 20 -20 20],[1 1 1 1],'white');
        view([90 0]);
        view(38.5,16);
        axis equal tight off; lighting phong; camlight right;
        %zoom(1.5);
        %view([0 90]);
        
        axis([-30 50 0 250 -75 75]);
        axis equal on;
        
        %set(gca,'CameraUpVector',[-1 0 0])
        drawnow
        caxis([val 10]);
        view(38.5,16);
        drawnow;
        frame = getframe(fig);
        aviobj = addframe(aviobj,frame);
    end
    aviobj = close(aviobj); close all
    command = ['mencoder temp.avi -oac mp3lame -lameopts cbr:br=128:vol=0 -srate 48000 -aid 1 -ovc lavc -sws 0 -lavcopts  threads=2:vcodec=mpeg4:vbitrate=940:keyint=240:vqmin=2:vqmax=15 -ofps 15.00 -vf scale=1024:768 -noodml -ffourcc DIVX  -o output.avi'];
    system(command);
    
end


load(fname)


C1 = [BodyPointsY(:,1) BodyPointsX(:,1) BodyPointsZ(:,1)]/scale;
C2 = [BodyPointsY(:,2) BodyPointsX(:,2) BodyPointsZ(:,2)]/scale;
C3 = [BodyPointsY(:,3) BodyPointsX(:,3) BodyPointsZ(:,3)]/scale;
C4 = [BodyPointsY(:,4) BodyPointsX(:,4) BodyPointsZ(:,4)]/scale;
BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);




cax2 = [min(Cp(:)) max(Cp(:))];

ctemp = linspace(cax2(1),cax2(2),31) - cax2(1);
ctemp = 1+31*ctemp/max(ctemp);
Cps = Cp;
Cps = Cps - cax2(1);
Cps = 1 + 31*Cps/(cax2(2) - cax2(1));
Cps = floor(Cps);
p2 = patch('Vertices',BodyPanPts,...
    'Faces',PtIDS,'FaceVertexCData',Cp(:),...
    'FaceColor','flat','EdgeColor','none');

box on
grid on
view(3)

axis equal tight
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
slice(XI,YI,ZI,VIs,[0],[ -5 -10 -15 -20 -25 -30],[])
axis equal
shading flat
box on
view(3)

axis equal tight
set(gca,'Projection','perspective')