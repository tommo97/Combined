clear all
clc
load RunData_000000.mat;

add_padding = true;
pad_cells = 2*[5 5 5];


Posns = CellPos;
Omegas = CellOms;

Xoffset = - min(Posns(:,1)) + 1;
Yoffset = - min(Posns(:,2)) + 1;
Zoffset = - min(Posns(:,3)) + 1;

Posns(:,1) = Posns(:,1) + Xoffset;
Posns(:,2) = Posns(:,2) + Yoffset;
Posns(:,3) = Posns(:,3) + Zoffset;



maxs = max(Posns);

if add_padding
Posns(:,1) = Posns(:,1) + pad_cells(1);
Posns(:,2) = Posns(:,2) + pad_cells(2);
Posns(:,3) = Posns(:,3) + pad_cells(3);

Xoffset = Xoffset + pad_cells(1);
Yoffset = Yoffset + pad_cells(2);
Zoffset = Zoffset + pad_cells(3);

maxs = max(Posns);
maxs = maxs + pad_cells;
end

Omx = zeros(maxs);
Omy = zeros(maxs);
Omz = zeros(maxs);


inds = sub2ind(maxs,Posns(:,1),Posns(:,2),Posns(:,3));


Omx(inds) = CellOms(:,1);
Omy(inds) = CellOms(:,2);
Omz(inds) = CellOms(:,3);

[YS XS ZS] = meshgrid(1:maxs(2),1:maxs(1),1:maxs(3));


xs = XS - Xoffset;
ys = YS - Yoffset;
zs = ZS - Zoffset;

Omegas = [Omx(:) Omy(:) Omz(:)];



[interp3(ys,xs,zs,Omx,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))...
interp3(ys,xs,zs,Omy,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))...
interp3(ys,xs,zs,Omz,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))]
CellOms(1807,:)


Posns = [xs(:) ys(:) zs(:)] + 0.5;

%save('InArrayA.mat','Posns');
%save('InArrayB.mat','Omegas');
%return
%return
load('OutArrayA.mat')
X = [Posns_x Posns_y Posns_z];
mins = min(X);


X(:,1) = round(X(:,1) - mins(1) + 1);
X(:,2) = round(X(:,2) - mins(2) + 1);
X(:,3) = round(X(:,3) - mins(3) + 1);

maxs = [max(X(:,1)) max(X(:,2)) max(X(:,3))];


u = zeros(maxs);
v = zeros(maxs);
w = zeros(maxs);
x = zeros(maxs);
y = zeros(maxs);
z = zeros(maxs);

omx = x;
omy = y;
omz = z;
inds = sub2ind(maxs,X(:,1),X(:,2),X(:,3));

x(inds) = X(:,1);
y(inds) = X(:,2);% + XCG/scale;
z(inds) = X(:,3);

u(inds) = Vels_x;%/GambitScale;
v(inds) = Vels_y;%/GambitScale;
w(inds) = Vels_z;%/GambitScale;

omx(inds) = Omegas_x;
omy(inds) = Omegas_y;
omz(inds) = Omegas_z;


[interp3(ys,xs,zs,omx,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))...
interp3(ys,xs,zs,omy,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))...
interp3(ys,xs,zs,omz,  CellPos(1807,2),  CellPos(1807,1),  CellPos(1807,3))]

xs = xs/(GambitScale*3.3); % 3.3 = 6.6/2. 6.6 is the aspect ratio of the McAlister blade, chordlength 1. This makes semispan +-1.
ys = ys/(GambitScale*3.3);
zs = zs/(GambitScale*3.3);
close all
isosurface(ys,xs,zs,sqrt(omx.*omx + omy.*omy + omz.*omz),5);
axis equal
view(3)



clf
linex = linspace(min(ys(:)),-0.25);%max(ys(:)))
XCG = GambitScale * Time * -10.0 / (3.3*GambitScale);
TE = XCG + cosd(12)*0.75/3.3; % foil TE is 0.75*c*cos(12) from 0.0
liney = ones(size(linex))*(TE + 4*1.0/3.3); % 0.2*c from foil TE...
hold all
for z = -0.05%0:-0.01:-0.3
linez = z*ones(size(linex));


W = interp3(ys,xs,zs,w,linex,liney,linez,'*cubic');

plot(linex+1,W/50);
end


