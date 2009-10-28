function V = LinVel(Start,End,Target,Gamma)

R1 = Target - Start;
R2 = Target - End;

C(:,1) = R1(:,2).*R2(:,3) - R1(:,3).*R2(:,2);
C(:,2) = R1(:,3).*R2(:,1) - R1(:,1).*R2(:,3);
C(:,3) = R1(:,1).*R2(:,2) - R1(:,2).*R2(:,1);
MagC = sqrt(C(:,1).*C(:,1) + C(:,2).*C(:,2) + C(:,3).*C(:,3));
MagR1 = sqrt(R1(:,1).*R1(:,1) + R1(:,2).*R1(:,2) + R1(:,3).*R1(:,3));
MagR2 = sqrt(R2(:,1).*R2(:,1) + R2(:,2).*R2(:,2) + R2(:,3).*R2(:,3));

R0 = End-Start;
% disp('-');
% disp([R0; R1; R2; C])
% disp('-');
R0R1 = R0(:,1).*R1(:,1) + R0(:,2).*R1(:,2) + R0(:,3).*R1(:,3);
R0R2 = R0(:,1).*R2(:,1) + R0(:,2).*R2(:,2) + R0(:,3).*R2(:,3);
Mult = Gamma./(MagC.*MagC*4*pi);
K = Mult.*((R0R1./MagR1) - (R0R2./MagR2));
K(MagC < 1e-16) = 0;
K(MagR1 < 1e-16) = 0;
K(MagR2 < 1e-16) = 0;
V = -[K K K].*C;