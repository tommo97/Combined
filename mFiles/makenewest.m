clear all
close all
cd ~/Desktop/Combined
a = dir('output');
s = size(a,1);
dates = zeros(1,s-2);
close all
for i = 3:s
    dates(-2+i) = a(i).datenum * a(i).isdir;
end

newest = find(dates==max(dates)) + 2;

disp(a(newest).date);
expr = ['output/' a(newest).name];
cd(expr);
a = dir('*.m');
expr = a(size(a,1)).name;
eval(expr(1:end-2));
%surf(X0,Y0,Z0,sqrt(Ox.^2+Oy.^2+Oz.^2));
disp(expr);
set(gcf,'Renderer','OpenGL');
%extract;
cd ..;
cd ..;