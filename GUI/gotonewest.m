function gotonewest(handles)
cdir = pwd;
cd(handles.bin_dir);
cd ..;

plot_ax = handles.goto_axes;
cla(plot_ax,'reset');

dump_000000;
cd(cdir);