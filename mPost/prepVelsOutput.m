clear all; clear mex; clc;
files = dir('RunData*5.mat');
close all
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale','Times','BodyRates0_x','BodyPointsX','BodyPointsY','BodyPointsZ')
load(fname,'Time','XCG')

scale = GambitScale;
val = 10;

XCG = scale * Time * 10;
THETA = Time * BodyRates0_x;

val = 50;
col(1) = 'r';
col(2) = 'b';


Data = load(fname,'CellOms','GambitScale','CellPos','TransVars_x','TransVars_y','TransVars_z');


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
Data.VI.x(Data.inds) = Data.CellOms(:,1);
Data.VI.y(Data.inds) = Data.CellOms(:,2);
Data.VI.z(Data.inds) = Data.CellOms(:,3);

Data.VI.max = max(Data.VI.mag(:));
Data.VI.mean = mean(Data.VI.mag(:));
Data.VI.total = sum(Data.VI.mag(:));
if exist('TransVars_x');
    Data.NumTransVars = size(Data.TransVars_x,2);
    
    for j = 1:Data.NumTransVars
        Data.V(j).mag = zeros(max(Data.subs));
        Data.V(j).mag(Data.inds) = sqrt(Data.TransVars_x(:,j) + Data.TransVars_y(:,j).^2 + Data.TransVars_z(:,j).^2);
        
    end
end
disp(['Domain size: ' num2str(size(Data.VI.mag)) '; i.e. ' num2str(numel(Data.VI.mag)) ' cells' ]);
disp(['Number of Vorticity Cells: ' num2str(length(Data.subs(:,1)))]);
disp(['Occupancy Ratio: ' num2str(length(Data.subs(:,1))/numel(Data.VI.mag))]);


[Data.XI Data.YI Data.ZI] = meshgrid([min(Data.CellPos(:,2)):1:max(Data.CellPos(:,2))]/Data.GambitScale,...
    [min(Data.CellPos(:,1)):1:max(Data.CellPos(:,1))]/Data.GambitScale,...
    [min(Data.CellPos(:,3)):1:max(Data.CellPos(:,3))]/Data.GambitScale);






YIMoved = Data.YI + XCG/scale;

slice(Data.XI,Data.YI+XCG/scale,Data.ZI,Data.VI.mag,[0],[10],[]);

return
%%  Try and calculate the velocity at this slice plane...



slicedata = slice(Data.XI,Data.YI+XCG/scale,Data.ZI,Data.XI,[],[10],[]);
sliceplaney = get(slicedata,'CData');
slicedata = slice(Data.XI,Data.YI+XCG/scale,Data.ZI,Data.YI,[],[10],[]);
sliceplanex = get(slicedata,'CData');
slicedata = slice(Data.XI,Data.YI+XCG/scale,Data.ZI,Data.ZI,[],[10],[]);
sliceplanez = get(slicedata,'CData');


sliceplanex = sliceplanex(2:2:end,2:2:end);
sliceplaney = sliceplaney(2:2:end,2:2:end);
sliceplanez = sliceplanez(2:2:end,2:2:end);



XS = Data.XI;
YS = Data.YI+XCG/scale;
ZS = Data.ZI;

OMX = Data.VI.x;
OMY = Data.VI.y;
OMZ = Data.VI.z;

Active = find(Data.VI.mag);


Vx = zeros(size(sliceplanez));
Vy = Vx;
Vz = Vx;

XSa = XS(Active);
YSa = YS(Active);
ZSa = ZS(Active);
Oms = [OMX(Active) OMY(Active) OMZ(Active)];

for i = 1:numel(sliceplanex)
    
    Dx = XSa - sliceplanex(i);
    Dy = YSa - sliceplaney(i);
    Dz = ZSa - sliceplanez(i);
    
    
    R = [Dx Dy Dz];
    Rmag2 = Dx.*Dx + Dy.*Dy +Dz.*Dz;
    den = sqrt(Rmag2 + 0.25);
    den3 = den.*den.*den;
    denom = [den3 den3 den3];
    
    K = -R./denom;
    
    Vs = sum(cross(K,Oms))/(4*pi);
    
    Vx(i) = Vx(i) + (Vs(1));
    Vy(i) = Vy(i) + (Vs(2));
    Vz(i) = Vz(i) + (Vs(3));
    
    disp(i)
end













return


axis equal
shading flat
box on
view(3)

Posns = [Data.XI(:) Data.YI(:) Data.ZI(:)];

Omegas = [Data.VI.x(:) Data.VI.y(:) Data.VI.z(:)]/(GambitScale.^2);

save('InArrayA.mat','Posns');
save('InArrayB.mat','Omegas');
XCG = scale * Time * 1.0;
pause

load('OutArrayA.mat')
%   Get Domain Bounds
X = DomainData(:,1:3);
mins = min(X);


X(:,1) = X(:,1) - mins(1) + 1;
X(:,2) = X(:,2) - mins(2) + 1;
X(:,3) = X(:,3) - mins(3) + 1;

maxs = [max(X(:,2)) max(X(:,1)) max(X(:,3))];



u = zeros(maxs);
v = zeros(maxs);
w = zeros(maxs);
x = zeros(maxs);
y = zeros(maxs);
z = zeros(maxs);

inds = sub2ind(maxs,X(:,2),X(:,1),X(:,3));

x(inds) = X(:,1);
y(inds) = X(:,2) + XCG/scale;
z(inds) = X(:,3);

u(inds) = DomainData(:,7)/GambitScale;
v(inds) = DomainData(:,8)/GambitScale;
w(inds) = DomainData(:,9)/GambitScale;

ZS = min(Data.YI(:) + XCG/scale) + 46.0;%5.03*[1 2 3 4 5];

for I = 1:length(ZS);
    sl_data = slice(Data.XI,Data.YI + XCG/scale,Data.ZI,w.*w + v.*v,[],[ZS(I)],[]);
    sliceplane = get(sl_data(1),'CData');
    sl_x = slice(Data.XI,Data.YI + XCG/scale,Data.ZI,Data.XI,[],[ZS(I)],[]);
    slicex = get(sl_x(1),'CData');
    sl_y = slice(Data.XI,Data.YI + XCG/scale,Data.ZI,Data.ZI,[],[ZS(I)],[]);
    slicey = get(sl_y(1),'CData');
    close all
    figure
    contourf(slicex,slicey,sliceplane)
    axis equal
    
    [THETA,RHO] = cart2pol(slicex,slicey);
    
    
    thetas = linspace(0,2*pi,101);
    rhos = linspace(0,6.5,99);
    [THETAS,RHOS] = meshgrid(thetas,rhos);
    [xs, ys] = pol2cart(THETAS,RHOS);
    
    vis = interp2(slicey,slicex,sliceplane,ys,xs);
    
    MeanVel{I} = mean(vis,2)';
    MeanRad{I} = mean(RHOS,2)';
    
    
end

close all
hold all
for i = 1:length(ZS)
    plot([fliplr(-MeanRad{i}) MeanRad{i}], [fliplr(MeanVel{i}) MeanVel{i}])
    leg{i} = ['$z=' num2str(ZS(i)) '$m']
end

axis equal
legend(leg)









figure
h = polar(xs,ys);
hold on;
contourf(xs,ys,vis);
% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image

