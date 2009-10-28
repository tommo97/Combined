mult = 5;
x = mult*X0(:);
y = mult*Y0(:);
z = mult*Z0(:);

x = floor(x);
y = floor(y);
z = floor(z);
x = x - min(x) + 1;
y = y - min(y) + 1;
z = z - min(z) + 1;
lims = [min(x) max(x) min(y) max(y) min(z) max(z)];
om = sqrt(Ox(:).*Ox(:) + Oy(:).*Oy(:) + Oz(:).*Oz(:));
[X Y Z] = ndgrid([lims(1):lims(2)],[lims(3):lims(4)],[lims(5):lims(6)]);
domain = zeros(size(X));

for i = 1:length(x)
    domain(x(i),y(i),z(i)) = domain(x(i),y(i),z(i)) + om(i);
end
domain = smooth3(domain,'gaussian',[mult mult mult],.95);
val = 10;%mean(om);
[f v] = isosurface(X,Y,Z,domain,val/mult);
n = isonormals(Y,X,Z,domain,v);

[fig p] = createfigure(n, Y, X, Z, v, f,true);