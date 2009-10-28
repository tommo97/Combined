cd ~/Desktop/Combined
a = dir('output');
s = size(a,1);
dates = zeros(1,s-2);
for i = 3:s
    dates(-2+i) = a(i).datenum * a(i).isdir;
end

newest = find(dates==max(dates)) + 2;

disp(a(newest).date);
expr = ['output/' a(newest).name];
cd(expr);
clear all;