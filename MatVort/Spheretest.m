clear all; close all; clc;
[x,y,z] = sphere;
Sphere.r = 10;
Sphere.c = [0 0 0];
Sphere.x = Sphere.r*x + Sphere.c(1);
Sphere.y = Sphere.r*y + Sphere.c(2);
Sphere.z = Sphere.r*z + Sphere.c(3);

for i = 1:size(Sphere.x,1) - 1
    for j = 1:size(Sphere.x,2) - 1
        c1 = [Sphere.x(i,j) Sphere.y(i,j) Sphere.z(i,j)];
        c2 = [Sphere.x(i+1,j) Sphere.y(i+1,j) Sphere.z(i+1,j)];
        c3 = [Sphere.x(i+1,j+1) Sphere.y(i+1,j+1) Sphere.z(i+1,j+1)];
        c4 = [Sphere.x(i,j+1) Sphere.y(i,j+1) Sphere.z(i,j+1)];
        cp = .25*(c1 + c2 + c3 + c4);
        Sphere.CP.x(i,j) = cp(1);
        Sphere.CP.y(i,j) = cp(2);
        Sphere.CP.z(i,j) = cp(3);
        
        Sphere.c1.x(i,j) = c1(1);
        Sphere.c1.y(i,j) = c1(2);
        Sphere.c1.z(i,j) = c1(3);
        
        Sphere.c2.x(i,j) = c2(1);
        Sphere.c2.y(i,j) = c2(2);
        Sphere.c2.z(i,j) = c2(3);
        
        Sphere.c3.x(i,j) = c3(1);
        Sphere.c3.y(i,j) = c3(2);
        Sphere.c3.z(i,j) = c3(3);
        
        Sphere.c4.x(i,j) = c4(1);
        Sphere.c4.y(i,j) = c4(2);
        Sphere.c4.z(i,j) = c4(3);
        
        t1 =  -.5*(.5*(c1+c2) -.5*(c3+c4));
        t1 = t1/sqrt(dot(t1,t1));
        n = cross(c4-c2,c3-c1);
        n = n/(sqrt(dot(n,n)));
        t2 = cross(t1,n);
        Sphere.n.x(i,j) = n(1);
        Sphere.n.y(i,j) = n(2);
        Sphere.n.z(i,j) = n(3);
        Sphere.t1.x(i,j) = t1(1);
        Sphere.t1.y(i,j) = t1(2);
        Sphere.t1.z(i,j) = t1(3);
        Sphere.t2.x(i,j) = t2(1);
        Sphere.t2.y(i,j) = t2(2);
        Sphere.t2.z(i,j) = t2(3);
        
    end
end

A = zeros(numel(Sphere.CP.x),numel(Sphere.CP.x));
B = A;
hold all
scatter3(Sphere.CP.x(:),Sphere.CP.y(:),Sphere.CP.z(:));
quiver3(Sphere.CP.x(:),Sphere.CP.y(:),Sphere.CP.z(:),Sphere.n.x(:),Sphere.n.y(:),Sphere.n.z(:));
quiver3(Sphere.CP.x(:),Sphere.CP.y(:),Sphere.CP.z(:),Sphere.t1.x(:),Sphere.t1.y(:),Sphere.t1.z(:));
quiver3(Sphere.CP.x(:),Sphere.CP.y(:),Sphere.CP.z(:),Sphere.t2.x(:),Sphere.t2.y(:),Sphere.t2.z(:));
axis equal

for i = 1:numel(Sphere.CP.x)
    trg = [Sphere.CP.x(i) Sphere.CP.y(i) Sphere.CP.z(i)];
    tnm = [Sphere.n.x(i) Sphere.n.y(i) Sphere.n.z(i)];
    for j = 1:numel(Sphere.CP.x)
        src = [Sphere.CP.x(j) Sphere.CP.y(j) Sphere.CP.z(j)];
        snm = [Sphere.n.x(j) Sphere.n.y(j) Sphere.n.z(j)];
        src = src + .25*snm;
        ome = snm;
        
        c1 = [Sphere.c1.x(j) Sphere.c1.y(j) Sphere.c1.z(j)];
        c2 = [Sphere.c2.x(j) Sphere.c2.y(j) Sphere.c2.z(j)];
        c3 = [Sphere.c3.x(j) Sphere.c3.y(j) Sphere.c3.z(j)];
        c4 = [Sphere.c4.x(j) Sphere.c4.y(j) Sphere.c4.z(j)];
        V = DirectVel(src - trg, ome);
        V1 = LinVel(c1,c2,trg,1);
        V2 = LinVel(c2,c3,trg,1);
        V3 = LinVel(c3,c4,trg,1);
        V4 = LinVel(c4,c1,trg,1);
        V = V1 + V2 + V3 + V4;
        %B(i,j) = sqrt(V(1)*V(1) + V(2)*V(2) + V(3)*V(3));
        A(i,j) = n(1)*V(1) + n(2)*V(2) + n(3)*V(3);
    end
    disp([i numel(Sphere.CP.x)]);
end

Sphere.RHS = zeros(size(Sphere.x,1) - 1,size(Sphere.y,1) - 1);
Vinf = [10 0 0];
for i = 1:size(Sphere.x,1) - 1
    for j = 1:size(Sphere.x,2) - 1
        Sphere.RHS(i,j) = ...
            sqrt((dot(Vinf,[Sphere.t1.x(i,j) Sphere.t1.y(i,j) Sphere.t1.z(i,j)])^2) + ...
            (dot(Vinf,[Sphere.t2.x(i,j) Sphere.t2.y(i,j) Sphere.t2.z(i,j)]))^2);       
    end
end
    
E = A;
E(end,:) = 1;
    
rhs = Sphere.RHS(:);

RHS = rhs;
RHS(end) = 0;
gam = E\RHS;

Sphere.Gamma = zeros(size(Sphere.RHS));
Sphere.Gamma(:) = gam;
surf(Sphere.x, Sphere.y, Sphere.z, Sphere.Gamma);
