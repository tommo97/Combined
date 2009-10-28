close all; clear all; clear mex;
files = dir('f*.dat');
val = 1;
s = size(files,1);
makemovie = false;
fname = files(s).name;
fullscreen = get(0,'ScreenSize');
fig = figure('Position',[0 -50 fullscreen(3)/4 fullscreen(4)/2]);
set(fig,'DoubleBuffer','on');
extract(fname,val);
set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
view(-38.5,16); axis equal off; lighting phong; camlight headlight; drawnow;
ViewParams = GetViewParams;



if makemovie
    aviobj = avifile(['temp.avi']);

for i = 1:s
    clf;
    fname = files(i).name;
    extract(fname,val);
    campos = ViewParams.campos;
    camtarget = ViewParams.camtarget;
    camva = ViewParams.camva;
    camproj = ViewParams.camproj;
    axis(ViewParams.axis);
    axis equal off;
    view(ViewParams.view)
    lighting phong; camlight headlight;
    drawnow;
    frame = getframe(fig);
    aviobj = addframe(aviobj,frame);
end
aviobj = close(aviobj); close all
command = ['mencoder temp.avi -oac mp3lame -lameopts cbr:br=128:vol=0 -srate 48000 -aid 1 -ovc lavc -sws 0 -lavcopts  threads=2:vcodec=mpeg4:vbitrate=940:keyint=240:vqmin=2:vqmax=15 -ofps 15.00 -vf scale=800:600 -noodml -ffourcc DIVX  -o output.avi'];
system(command);
delete temp.avi;
end