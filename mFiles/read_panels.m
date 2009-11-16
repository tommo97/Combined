clear all
close all
s = 0;
while s < 300
files = dir('*.m');
s = size(files,1);
disp(s);
%pause(10);
end

fullscreen = get(0,'ScreenSize');
fig = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
set(fig,'DoubleBuffer','on');
aviobj = avifile(['temp.avi']);

set(0,'DefaultSurfaceLine','none')

for i = 1:s
    clf;
    a = files(i).name(1:end-2);
    eval(a);
     
    %camlight left;
   
    axis([-10 40 -15 15 0 4]);
    zoom(1);
    zoom(2);
    view(-44, 18);
    set(gcf,'Renderer','zbuffer');
    set(gcf,'Renderer','OpenGL')
    lighting phong; axis off;
    drawnow;
    frame = getframe(fig);
    aviobj = addframe(aviobj,frame);
    disp(i);
end

aviobj = close(aviobj); close all
command = ['mencoder temp.avi -oac mp3lame -lameopts cbr:br=128:vol=0 -srate 48000 -aid 1 -ovc lavc -sws 0 -lavcopts  threads=2:vcodec=mpeg4:vbitrate=940:keyint=240:vqmin=2:vqmax=15 -ofps 15.00 -vf scale=1024:768 -noodml -ffourcc DIVX  -o video.avi'];
i = system(command);
set(0,'DefaultSurfaceLine','-')
if i==0 
    delete temp.avi
end