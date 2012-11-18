clear all
load RunData_000005.mat;

Posns = CellPos;
Omegas = CellOms;

Posns(:,1) = Posns(:,1) - min(Posns(:,1)) + 1;
Posns(:,2) = Posns(:,2) - min(Posns(:,2)) + 1;
Posns(:,3) = Posns(:,3) - min(Posns(:,3)) + 1;

maxs = max(Posns);

Omx = zeros(maxs);
Omy = zeros(maxs);
Omz = zeros(maxs);

inds = sub2ind(maxs,Posns(:,1),Posns(:,2),Posns(:,3));


Omx(inds) = CellOms(:,1);
Omy(inds) = CellOms(:,2);
Omz(inds) = CellOms(:,3);
[YS XS ZS] = meshgrid(1:maxs(2),1:maxs(1),1:maxs(3));

Posns = [XS(:) YS(:) ZS(:)];
Omegas = [Omx(:) Omy(:) Omz(:)];

save('InArrayA.mat','Posns');
save('InArrayB.mat','Omegas');

