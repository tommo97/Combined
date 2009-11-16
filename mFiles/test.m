clear all
P0x = 1:10;
P0y = ones(size(P0x));
P0z = zeros(size(P0x));

close all

Vp = [0 1 0];
Vi = [0 0 0];

dt = 1;
scatter3(P0x,P0y,P0z);
TEx = P0x;
TEy = P0y;
TEz = P0z;
Wx = TEx;
Wy = TEy;
Wz = TEz;

hold on
drawnow
for i = 1:10
    hold off

    Wx(1,:) = Wx(1,:) + .666*dt*Vp(1);
    Wy(1,:) = Wy(1,:) + .666*dt*Vp(2);
    Wz(1,:) = Wz(1,:) + .666*dt*Vp(3);
    Wx(2:end,:) = Wx(2:end,:) + dt*Vi(1);
    Wy(2:end,:) = Wy(2:end,:) + dt*Vi(2);
    Wz(2:end,:) = Wz(2:end,:) + dt*Vi(3);
    P0x = P0x + Vp(1)*dt;
    P0y = P0y + Vp(2)*dt;
    P0z = P0z + Vp(3)*dt;
    TEx = P0x;
    TEy = P0y;
    TEz = P0z;


        Wx = [TEx;Wx];
    Wy = [TEy;Wy];
    Wz = [TEz;Wz];


    scatter3(P0x,P0y,P0z);
    hold on
    surf(Wx,Wy, zeros(size(Wx)));
    drawnow
end
