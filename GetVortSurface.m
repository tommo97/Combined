close all; clear all; clear mex;
files = dir('f_be*.dat');
val = 7.5;
s = size(files,1);

% while s < 100
%     files = dir('f*.dat');
%     s = size(files,1);
%     disp(s)
%     pause(60);
% end



makemovie = false;
fname = files(s).name;
fullscreen = get(0,'ScreenSize');
fig = figure('Position',[0 0 1200 1200*16/10]);
set(fig,'DoubleBuffer','on');
extract(fname,val);
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
hold all
scatter3([0 0 0 0],[0 0 410 410],[-20 20 -20 20],[1 1 1 1],'white');
%view([0 90]);
%view(3);
axis([-30 50 0 250 -75 75]);
axis equal off;

%set(gca,'CameraUpVector',[-1 0 0])
axis equal tight; lighting phong; camlight right;
%zoom(1.8);
caxis([val 10])
drawnow
ViewParams = GetViewParams;
%camlight right
view(38.5,16);


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
        scatter3([0 0 0 0],[0 0 410 410],[-20 20 -20 20],[1 1 1 1],'white');
        view([90 0]);
        view(38.5,16);
        axis equal tight off; lighting phong; camlight right;
        %zoom(1.5);
        %view([0 90]);
        
axis([-30 50 0 250 -75 75]);
axis equal off;

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
