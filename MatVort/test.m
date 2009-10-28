dif = rand;
om = rand;



for i = 1:1000;
    [V1 V2] = DirectVel2D(dif,om,i*dif/10);  
    sc(i) = i/10;
    rat(i) = V2./V1;
end
close all
plot(sc,rat)
disp(mean(rat./sc));