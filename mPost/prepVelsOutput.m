clear all; clear mex; clc;
files = dir('RunData*9.mat');
close all
s = size(files,1);
fname = files(s).name;
load(fname,'GambitScale','Times','BodyRates0_x','BodyPointsX','BodyPointsY','BodyPointsZ')
load(fname,'Time','XCG')

scale = GambitScale;
val = 10;

XCG = scale * Time * 10;
THETA = Time * BodyRates0_x;

[XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val,scale);

YIMoved = YI + XCG/scale;

slice(XI,YI+XCG/scale,ZI,VI,[0],[],[]);
axis equal
shading flat
box on
view(3)

Posns = [XI(:) YI(:) ZI(:)];

Omegas = [VIx(:) VIy(:) VIz(:)]/(GambitScale.^2);

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

ZS = min(YI(:) + XCG/scale) + 5.03*[1 2 3 4 5];

for I = 1:length(ZS);
    sl_data = slice(XI,YI + XCG/scale,ZI,v+10,[],[ZS(I)],[]);
    sliceplane = get(sl_data(1),'CData');
    sl_x = slice(XI,YI + XCG/scale,ZI,XI,[],[ZS(I)],[]);
    slicex = get(sl_x(1),'CData');
    sl_y = slice(XI,YI + XCG/scale,ZI,ZI,[],[ZS(I)],[]);
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

