clear all
mode1 = false;
uinf = [10 0 -1];

infty = 1e6;

if mode1
    [Xdata PNLS] = ReadNeu('0012.neu');
else
    [Xdata PNLS] = ReadNeu('0012_flat.neu');
end


C1 = Xdata(PNLS(:,1),:);
C2 = Xdata(PNLS(:,2),:);
C3 = Xdata(PNLS(:,3),:);
C4 = Xdata(PNLS(:,4),:);


CP = 0.25*(C1+C2+C3+C4);

Vinf = CP;
Vinf(:,1) = uinf(1);
Vinf(:,2) = uinf(2);
Vinf(:,3) = uinf(3);

if ~mode1
    P1 = 0.5*(0.5*(C1+C2) + C1);
    P4 = 0.5*(0.5*(C3+C4) + C4);
    
    CP = 0.5*(CP + 0.5*(C2+C3));
    Disp = infty*Vinf;
    P2 = P1 + Disp;
    P3 = P4 + Disp;
end

D1 = C3 - C1;
D2 = C4 - C2;
CR = cross(D1,D2);
CRmag = sqrt(dot(CR,CR,2));
n = CR./[CRmag CRmag CRmag];
Area = 0.5*CRmag;
A = zeros(length(CP));
TRG = CP;
NRM = CP;



for i = 1:length(CP)
    TRG(:,1) = CP(i,1);
    TRG(:,2) = CP(i,2);
    TRG(:,3) = CP(i,3);
    NRM(:,1) = n(i,1);
    NRM(:,2) = n(i,2);
    NRM(:,3) = n(i,3);
    
    
    if mode1
        V = PanVel(C1,C2,C3,C4,TRG,ones(length(CP),1));
    else
        %V = LinVel(P1,P2,TRG,1);
        V = PanVel(P1,P2,P3,P4,TRG,ones(length(CP),1));
    end
    Vn = dot(NRM,V,2);
    
    A(i,:) = Vn';
    disp(i)
end

RHS = dot(Vinf,n,2);
if mode1
    A(end,:) = 1;
    
    RHS(end) = 0;
end
GAMMA = A\RHS;

clf
hold on
for i = 1:length(GAMMA)
    SurfPan(C1(i,:), C4(i,:), C2(i,:), C3(i,:), GAMMA(i));
end
