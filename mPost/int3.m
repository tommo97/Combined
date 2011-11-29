function out = int3(Xv,Yv,Zv,Vol,xi,yi,zi)

[l m n] = size(Vol);

minx = Xv(1);
maxx = Xv(end);
dx = (maxx - minx)/l;
deltax = (xi - minx)/dx;
indx = floor(deltax)+1;
xd = deltax - indx;


miny = Yv(1);
maxy = Yv(end);
dy = (maxy - miny)/m;
deltay = (yi - miny)/dy;
indy = floor(deltay)+1;
yd = deltay - indy;


minz = Zv(1);
maxz = Zv(end);
dz = (maxz - minz)/n;
deltaz = (zi - minz)/dz;
indz = floor(deltaz)+1;
zd = deltaz - indz;

i1 = Vol(indx,indy,indz)*(1-zd) + Vol(indx,indy,indz+1)*zd;
i2 = Vol(indx,indy+1,indz)*(1-zd) + Vol(indx,indy+1,indz+1)*zd;
j1 = Vol(indx+1,indy,indz)*(1-zd) + Vol(indx+1,indy,indz+1)*zd;
j2 = Vol(indx+1,indy+1,indz)*(1-zd) + Vol(indx+1,indy+1,indz+1)*zd;

w1 = i1*(1-yd) + i2*yd;
w2 = j1*(1-yd) + j2*yd;

out = w1*(1-xd) + w2*xd;



