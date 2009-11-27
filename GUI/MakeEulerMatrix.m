function TRANS = MakeEulerMatrix(Attitude)
EulerAngles = Attitude;
cosphi = cos(EulerAngles(3)); costhe = cos(EulerAngles(2)); cospsi = cos(EulerAngles(1));
sinphi = sin(EulerAngles(3)); sinthe = sin(EulerAngles(2)); sinpsi = sin(EulerAngles(1));
a1 = costhe*cospsi;
a2 = costhe*sinpsi;
a3 = -sinthe;
b1 = sinphi*sinthe*cospsi - cosphi*sinpsi;
b2 = sinphi*sinthe*sinpsi + cosphi*cospsi;
b3 = sinphi*costhe;
c1 = cosphi*sinthe*cospsi + sinphi*sinpsi;
c2 = cosphi*sinthe*sinpsi - sinphi*cospsi;
c3 = cosphi*costhe;
TRANS(1,:) = [a1 a2 a3];
TRANS(2,:) = [b1 b2 b3];
TRANS(3,:) = [c1 c2 c3];