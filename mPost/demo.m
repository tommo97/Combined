clc
clear

in=round(10*randn(40,30,20));
[y x z]=meshgrid(1:.3:40,1:.2:30,1:.1:20);

start_time1 = clock;
out=interp3(in,x,y,z,'*linear',0);
end_time_1 = clock;
start_time2 = clock;
out2=trilinear(in,x,y,z);
end_time_2 = clock;
err=abs((out) - (out2));
sum(err(:))

time1 = etime(end_time_1,start_time1)
time2 = etime(end_time_2,start_time2)

time_difference = etime(end_time_1,start_time1) - etime(end_time_2,start_time2)

time2/time1 