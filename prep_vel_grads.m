
load RunData_000005.mat

scale = GambitScale;

minx = min(CellPos(:,1));
miny = min(CellPos(:,2));
minz = min(CellPos(:,3));

maxx = max(CellPos(:,1));
maxy = max(CellPos(:,2));
maxz = max(CellPos(:,3));

DomOmx = zeros(1 + maxx-minx,1 + maxy-miny,1 + maxz-minz);
DomOmy = zeros(1 + maxx-minx,1 + maxy-miny,1 + maxz-minz);
DomOmz = zeros(1 + maxx-minx,1 + maxy-miny,1 + maxz-minz);

subs = 1 + [CellPos(:,1) - minx CellPos(:,2) - miny CellPos(:,3) - minz];

inds = sub2ind(size(Domain),subs(:,1),subs(:,2),subs(:,3));

DomOmx(inds) = CellOms(:,1);
DomOmy(inds) = CellOms(:,2);
DomOmz(inds) = CellOms(:,3);

[XI YI ZI] = meshgrid([min(CellPos(:,2)):1:max(CellPos(:,2))]/scale,...
    [min(CellPos(:,1)):1:max(CellPos(:,1))]/scale,...
    [min(CellPos(:,3)):1:max(CellPos(:,3))]/scale);

Posns = [XI(:) YI(:) ZI(:)];
Omegas = [DomOmx(:) DomOmy(:) DomOmz(:)]/(scale*scale);