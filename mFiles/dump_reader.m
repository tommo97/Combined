clear
close all;
dump_000150;
set(gcf,'Renderer','OpenGL');
close all;
x = X0(:);
y = Y0(:);
z = Z0(:);
x = floor(x) + .5;
y = floor(y) + .5;
z = floor(z) + .5;
x = x-min(x) + 1;
y = y-min(y) + 1;
z = z-min(z) + 1;
lims = [min(x) max(x) min(y) max(y) min(z) max(z)];
[X Y Z] = ndgrid(lims(1):lims(2),lims(3):lims(4),lims(5):lims(6));
domain = zeros(size(X));
Oms = sqrt(Ox(:).*Ox(:) + Oy(:).*Oy(:) + Oz(:).*Oz(:));
for i = 1:length(x)
domain(x(i),y(i),z(i)) = Oms(i);
end
val = mean(Oms(:));
domain = smooth3(domain);
[f v] = isosurface(X,Y,Z,domain,val);
n = isonormals(Y,X,Z,domain,v);
[fig p] = createfigure(n, X, Y, Z, v, f,true);
% imagesc(domain(:,:,2));
% axis square
%slice(Y,X,Z,domain,[30],[30],[3])

