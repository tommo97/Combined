function V = DirVel(X,Y,Z,P,Omega)

R(:,:,1) = P(1) - X;
R(:,:,2) = P(2) - Y;
R(:,:,3) = P(3) - Z;

DEL2 = 1e-6;

R2 = (R(:,:,1).*R(:,:,1) + R(:,:,2).*R(:,:,2) + R(:,:,3).*R(:,:,3));
nrm = sqrt(DEL2 + R2);
mult = -1./(4*pi*nrm.*nrm.*nrm);
C(:,:,1) = R(:,:,2).*Omega(:,:,3) - R(:,:,3).*Omega(:,:,2);
C(:,:,2) = R(:,:,3).*Omega(:,:,1) - R(:,:,1).*Omega(:,:,3);
C(:,:,3) = R(:,:,1).*Omega(:,:,2) - R(:,:,2).*Omega(:,:,1);
Vel(:,:,1) = mult.*C(:,:,1);
Vel(:,:,2) = mult.*C(:,:,2);
Vel(:,:,3) = mult.*C(:,:,3);

Vx = Vel(:,:,1);
Vy = Vel(:,:,2);
Vz = Vel(:,:,3);

V = [sum(Vx(:)) sum(Vy(:)) sum(Vz(:))];