figure
files = dir('data*.dat');
val = 5;
s = size(files,1);
data = dlmread(files(s).name);
subs = data(:,1:3);
max(subs)
domain = zeros(max(subs));
U = domain;
V = domain;
W = domain;
X = domain;
Y = domain;
Z = domain;
for i = 1:size(subs,1)
    X(subs(i,1),subs(i,2),subs(i,3)) = data(i,4);
    Y(subs(i,1),subs(i,2),subs(i,3)) = data(i,5);
    Z(subs(i,1),subs(i,2),subs(i,3)) = data(i,6);
    U(subs(i,1),subs(i,2),subs(i,3)) = data(i,7);
    V(subs(i,1),subs(i,2),subs(i,3)) = data(i,8);
    W(subs(i,1),subs(i,2),subs(i,3)) = data(i,9);
end
    
U = U + 10;

X = X - min(X(:));

%-----Define viewing and lighting
hold all

sx = [-20:1:20];
sy = min(X(:))*ones(size(sx));
sz = -0.5*ones(size(sx));

verts = stream3(Y,X,Z,V,U,W,sx,sy,sz);
[curlx,curly,curlz,cav] = curl(Y,X,Z,V,U,W);
spd = sqrt(U.^2 + V.^2 + W.^2);
streamribbon(verts,Y,X,Z,cav,spd,2);        
domain = sqrt(curlx.*curlx + curly.*curly + curlz.*curlz);         
axis tight
shading interp


cdata = 10*curlz;

p = patch(isosurface(Y,X,Z,domain,val));
isonormals(Y,X,Z,domain,p);
isocolors(Y,X,Z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none')


set(gcf,'Color',[1,1,1],'Renderer','OpenGL')
view(-38.5,16); axis equal off; lighting phong; camlight headlight;
