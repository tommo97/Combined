syms m12 m23 m34 m41 e1 e2 e3 e4 r1 r2 r3 r4 h1 h2 h3 h4 d12 d23 d34 d41
syms x1 x2 x3 x4 y1 y2 y3 y4 x y z 
syms num denom

d12 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
d23 = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2));
d34 = sqrt((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3));
d41 = sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4));

m12 = (y2-y1) / (x2-x1);
m23 = (y3-y2) / (x3-x2);
m34 = (y4-y3) / (x4-x3);
m41 = (y1-y4) / (x1-x4);

r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
r3 = sqrt((x-x3)^2 + (y-y3)^2 + z^2);
r4 = sqrt((x-x4)^2 + (y-y4)^2 + z^2);

e1 = (x-x1)^2 + z^2;
e2 = (x-x2)^2 + z^2;
e3 = (x-x3)^2 + z^2;
e4 = (x-x4)^2 + z^2;

h1 = (x-x1)*(y-y1);
h2 = (x-x2)*(y-y2);
h3 = (x-x3)*(y-y3);
h4 = (x-x4)*(y-y4);

alpha = (m12*e1 - h1)/(z*r1);
beta = (m12*e2 - h2)/(z*r2);

num = alpha - beta;
denom = 1 + alpha*beta;






syms Numa Numb Numc Numd Dena Denb Denc Dend Aa Ab Ac Ad X Y Z
Aa = (Numa+Dena)/(1-Numa*Dena);
Ab = (Numb+Denb)/(1-Numb*Denb);
Ac = (Numc+Denc)/(1-Numc*Denc);
Ad = (Numd+Dend)/(1-Numd*Dend);

X = (Aa + Ab)/(1 - Aa*Ab);
Y = (Ac + Ad)/(1 - Ac*Ad);

Z = (X + Y)/(1 - X*Y)






