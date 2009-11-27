function gotonewest(handles)
cdir = pwd;
cd ~/Desktop/Workspace/Combined
a = dir('~/output');
s = size(a,1);
dates = zeros(1,s-2);
for i = 3:s
    dates(-2+i) = a(i).datenum * a(i).isdir;
end

newest = find(dates==max(dates)) + 2;

disp(a(newest).date);
expr = ['~/output/' a(newest).name];
cd(expr);
plot_ax = handles.goto_axes;
cla(plot_ax,'reset');
a = dir;
dump_000000;
cd(cdir);