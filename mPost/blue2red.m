function M = blue2red(datamin, datamax,n)
%   We have a maximum and a minimum value, 
%   want zero to be white, and want n colours in map.

inds = 1:n;
range = linspace(datamin,datamax,n);

%   Get extended range
s = max(abs(datamax),abs(datamin));
range_ex = linspace(-s,s,n);

new_inds = interp1(range_ex,inds,range);

nx = round(2*n*s/(datamax - datamin));

nx = 2+(2*floor(nx/2));
mapoid =  redblue(nx);


m(:,1) = interp1(linspace(1,n,nx),mapoid(:,1),new_inds,'nearest');
m(:,2) = interp1(linspace(1,n,nx),mapoid(:,2),new_inds,'nearest');
m(:,3) = interp1(linspace(1,n,nx),mapoid(:,3),new_inds,'nearest');

m(m(:)>1) = 1.;

M = m;
return
%   We have a maximum and a minimum value, 
%   want zero to be white, and want n colours in map.


range = linspace(datamin,datamax,n);

%   Get extended range
s = max(abs(datamax),abs(datamin));
nx = 2*s*n/(datamax - datamin);
range_ex = linspace(-s,s,nx);
inds = 1:nx;
new_inds = interp1(range_ex,inds,range,'nearest')+1;

mapoid = redblue(nx+1);
m(:,1) = mapoid(new_inds,1);
m(:,2) = mapoid(new_inds,2);
m(:,3) = mapoid(new_inds,3);
m(m(:)>1) = 1.;

M = m;
return