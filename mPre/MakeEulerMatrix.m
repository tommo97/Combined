function TRANS = MakeEulerMatrix(Attitude)
EulerAngles = Attitude;
cosphi = cos(EulerAngles(1)); costhe = cos(EulerAngles(2)); cospsi = cos(EulerAngles(3));
sinphi = sin(EulerAngles(1)); sinthe = sin(EulerAngles(2)); sinpsi = sin(EulerAngles(3));
a1 = costhe*cospsi;
a2 = costhe*sinpsi;
a3 = -sinthe;
b1 = sinphi*sinthe*cospsi - cosphi*sinpsi;
b2 = sinphi*sinthe*sinpsi + cosphi*cospsi;
b3 = sinphi*costhe;
c1 = cosphi*sinthe*cospsi + sinphi*sinpsi;
c2 = cosphi*sinthe*sinpsi - sinphi*cospsi;
c3 = cosphi*costhe;
TRANS(1,:) = [a1 b1 c1];
TRANS(2,:) = [a2 b2 c2];
TRANS(3,:) = [a3 b3 c3];