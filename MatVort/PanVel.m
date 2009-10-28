function V = PanVel(C1,C2,C3,C4,P,Gamma)

V1 = LinVel(C1,C2,P,Gamma);
V2 = LinVel(C2,C3,P,Gamma);
V3 = LinVel(C3,C4,P,Gamma);
V4 = LinVel(C4,C1,P,Gamma);

V = V1 + V2 + V3 + V4;
% disp('!!!')
% disp([C1;C2;C3;C4;P;V]);
% disp(Gamma);
% disp('---')