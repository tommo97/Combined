clear all; clear mex; clc;
files = dir('RunData*5.mat');
close all
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale')


Data = load(fname,'CellOms','GambitScale','CellPos','TransVars_x','TransVars_y','TransVars_z');

clear Data2;

Data.subs = Data.CellPos(:,1:3);
Data.subs(:,1) = Data.subs(:,1) - min(Data.CellPos(:,1));
Data.subs(:,2) = Data.subs(:,2) - min(Data.CellPos(:,2));
Data.subs(:,3) = Data.subs(:,3) - min(Data.CellPos(:,3));
Data.subs = Data.subs+1;
Data.inds = sub2ind(max(Data.subs),Data.subs(:,1),Data.subs(:,2),Data.subs(:,3));
Data.do_multiple_vorts = false;

Data.VI.mag = zeros(max(Data.subs));
Data.VI.mag(Data.inds) = sqrt(Data.CellOms(:,1).^2 + Data.CellOms(:,2).^2 + Data.CellOms(:,3).^2);
Data.VI.x = zeros(max(Data.subs));
Data.VI.y = zeros(max(Data.subs));
Data.VI.z = zeros(max(Data.subs));

Data.P.x = zeros(max(Data.subs));
Data.P.y = zeros(max(Data.subs));
Data.P.z = zeros(max(Data.subs));

Data.P.x(Data.inds) = Data.CellPos(:,1);
Data.P.y(Data.inds) = Data.CellPos(:,2);
Data.P.z(Data.inds) = Data.CellPos(:,3);
Data.VI.x(Data.inds) = Data.CellOms(:,1);
Data.VI.y(Data.inds) = Data.CellOms(:,2);
Data.VI.z(Data.inds) = Data.CellOms(:,3);

Xlocus = -5;


[xi yi] = meshgrid(min(Data.P.z(:)):max(Data.P.z(:)),min(Data.P.y(:)):max(Data.P.y(:)));


slicedata = slice(Data.VI.x,[],10,[]);

slicevals = get(slicedata, 'CData');
close all

u = zeros(size(xi));
v = zeros(size(xi));


del2 = 0.25;

for i = 1:numel(xi)
    
    
    for j = 1:numel(xi)
        dy =  (yi(j)-yi(i));
        dx =  (xi(j)-xi(i));
        denom = slicevals(j)/(2*pi*(dx*dx + dy*dy + del2));
        u(i) = u(i) - denom*dy;
        v(i) = v(i) + denom*dx;
        
    end
    
end


load('OutArrayA.mat')
X = [Posns_x Posns_y Posns_z];
mins = min(X);


X(:,1) = round(X(:,1) - mins(1) + 1);
X(:,2) = round(X(:,2) - mins(2) + 1);
X(:,3) = round(X(:,3) - mins(3) + 1);

maxs = [max(X(:,2)) max(X(:,1)) max(X(:,3))];



u = zeros(maxs);
v = zeros(maxs);
w = zeros(maxs);
x = zeros(maxs);
y = zeros(maxs);
z = zeros(maxs);
omx = x;
omy = y;
omz = z;
inds = sub2ind(maxs,X(:,2),X(:,1),X(:,3));

x(inds) = X(:,1);
y(inds) = X(:,2);% + XCG/scale;
z(inds) = X(:,3);

u(inds) = Vels_x;%/GambitScale;
v(inds) = Vels_y;%/GambitScale;
w(inds) = Vels_z;%/GambitScale;

omx(inds) = Omegas_x;
omy(inds) = Omegas_y;
omz(inds) = Omegas_z;

slicedata = slice(omz,[],[1:20:400],[]);
sliceplane = get(slicedata,'CData');
shading flat
axis equal