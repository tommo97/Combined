clear all
clc
load RunData_000002.mat;

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

save('InArrayA.mat','Posns');
save('InArrayB.mat','Omegas');

