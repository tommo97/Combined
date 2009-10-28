%function USLiftingLineBlob
clear all
%close all

tic
NumSteps = 50;
dt = .02;

%%  Define blade geometry
%   Span lies along y axis, chord parallel to x axis
n = 50;
b = 5;
C = 2;
[Xdata Ydata Zdata] = meshgrid(-.5*C:.5*C, 10*linspace(-b,b,n), 0:0);

C1 = [Xdata(1:end-1,1) Ydata(1:end-1,1) Zdata(1:end-1,1)];
C2 = [Xdata(1:end-1,2) Ydata(1:end-1,2) Zdata(1:end-1,2)];
C3 = [Xdata(2:end,2) Ydata(2:end,2) Zdata(2:end,2)];
C4 = [Xdata(2:end,1) Ydata(2:end,1) Zdata(2:end,1)];

Centroid = 0.25*(C1+C2+C3+C4);

CP = Centroid;%0.5*(Centroid + 0.5*(C2+C3));

dr = diff(CP(1:2,:));

D1 = C3 - C1;
D2 = C4 - C2;
CR = cross(D1,D2);
CRmag = sqrt(dot(CR,CR,2));
Norm = CR./[CRmag CRmag CRmag];
Area = 0.5*CRmag;
A = zeros(length(CP));
TRG = CP;
NRM = CP;

%%  Some inflow type parameters
Uinf = 10*[10 0 2];
UinfMag = sqrt(dot(Uinf,Uinf,2));
Vinf = [Uinf(1) * ones(n-1,1) Uinf(2) * ones(n-1,1) Uinf(3) * ones(n-1,1)];

%%  Place a vortex ring of about chord length starting at c/4
P1 = 0.5*(0.5*(C1+C2) + C1);
P4 = 0.5*(0.5*(C3+C4) + C4);
Disp = C*Vinf/UinfMag;
P2 = P1 + Disp;
P3 = P4 + Disp;

%%  Get Influence Matrix

for i = 1:length(CP)
    TRG(:,1) = CP(i,1);
    TRG(:,2) = CP(i,2);
    TRG(:,3) = CP(i,3);
    NRM(:,1) = Norm(i,1);
    NRM(:,2) = Norm(i,2);
    NRM(:,3) = Norm(i,3);
    %V = LinVel(P1,P4,TRG,1);
    V = PanVel(P1,P2,P3,P4,TRG,ones(length(CP),1));
    
    Vn = dot(NRM,V,2);
    
    A(i,:) = Vn';
    disp(i)
end

%%  Get Initial RHS and solve for GAMMA

RHS = dot(-Vinf,Norm,2);

GAMMA = A\RHS;

GAMMA_prev = zeros(size(GAMMA));
dgammadt = (GAMMA - GAMMA_prev)/dt;
GAMMA_prev = GAMMA;


V = Vinf;


%   Interpolate up

R = CP(:,2);
r = linspace(-b,b);
gamma = interp1(R,GAMMA,r,'cubic','extrap');
dgammadr = diff(gamma)./diff(r);
delGamma = interp1(r(1:end-1) + .5*diff(r),dgammadr, R, 'cubic','extrap');

dOmdt(:,1) = dgammadt*dr(1);
dOmdt(:,2) = dgammadt*dr(2);
dOmdt(:,3) = dgammadt*dr(3);

uDelOm(:,1) = V(:,1).*delGamma;
uDelOm(:,2) = V(:,2).*delGamma;
uDelOm(:,3) = V(:,3).*delGamma;
Om = dt*(uDelOm + dOmdt);

%%  Shed a vortex blob

WakePtsX = [];
WakePtsY = [];
WakePtsZ = [];
WakeOMEGAx = [];
WakeOMEGAy = [];
WakeOMEGAz = [];
%%  Shed new wake blob
WakePtsX = [CP(:,1)+.5*C   WakePtsX ];
WakePtsY = [CP(:,2)     WakePtsY ];
WakePtsZ = [CP(:,3)     WakePtsZ ];
WakeOMEGAx = [Om(:,1) WakeOMEGAx];
WakeOMEGAy = [Om(:,2) WakeOMEGAy];
WakeOMEGAz = [Om(:,3) WakeOMEGAz];


counter = 0;
%%  Start of time marching loop
for i = 0:NumSteps
    %%  Convect the wake
    %   Get velocities at points due to lifting line and other wake panels
    WakeVelX = zeros(size(WakePtsX));
    WakeVelY = zeros(size(WakePtsY));
    WakeVelZ = zeros(size(WakePtsZ));
    OM = [];
    OM(:,:,1) = WakeOMEGAx;
    OM(:,:,2) = WakeOMEGAy;
    OM(:,:,3) = WakeOMEGAz;
    
    
    for k = 1:numel(WakePtsX)
        trg(1) = WakePtsX(k);
        trg(2) = WakePtsY(k);
        trg(3) = WakePtsZ(k);
        TRG(:,1) = trg(1);
        TRG(:,2) = trg(2);
        TRG(:,3) = trg(3);
        
        V = PanVel(P1,P2,P3,P4,TRG,GAMMA);
        
        WakeVel = DirVel(WakePtsX,WakePtsY,WakePtsZ,trg,OM);
        
        WakeVelX(k) = Uinf(1) + sum(V(:,1)) - WakeVel(1);
        WakeVelY(k) = Uinf(2) + sum(V(:,2)) - WakeVel(2);
        WakeVelZ(k) = Uinf(3) + sum(V(:,3)) - WakeVel(3);
        
    end
    
    
    
    WakePtsX = WakePtsX + dt*WakeVelX;
    WakePtsY = WakePtsY + dt*WakeVelY;
    WakePtsZ = WakePtsZ + dt*WakeVelZ;
    
    counter = counter + numel(WakePtsX);
    %%  Calculate downwash on blade and solve linear algebra
    V = Vinf;
    for j = 1:size(CP,1)
        trg(1) = CP(j,1);
        trg(2) = CP(j,2);
        trg(3) = CP(j,3);
        WakeVel = DirVel(WakePtsX,WakePtsY,WakePtsZ,trg,OM);
        
        V(j,:) = V(j,:) - WakeVel;
        
        
    end
    
    Vn = dot(Norm,-V,2);
    
    GAMMA = A\Vn;
    
    dgammadt = (GAMMA - GAMMA_prev)/dt;
    GAMMA_prev = GAMMA;
    
    %%  Shed new wake element
    
    gamma = interp1(R,GAMMA,r,'cubic','extrap');
    dgammadr = diff(gamma)./diff(r);
    delGamma = interp1(r(1:end-1) + .5*diff(r),dgammadr, R, 'cubic','extrap');
    dOmdt(:,1) = dgammadt*dr(1);
    dOmdt(:,2) = dgammadt*dr(2);
    dOmdt(:,3) = dgammadt*dr(3);
    
    uDelOm(:,1) = V(:,1).*delGamma;
    uDelOm(:,2) = V(:,2).*delGamma;
    uDelOm(:,3) = V(:,3).*delGamma;
    
    
    Om = dt*(uDelOm + dOmdt);
    
    WakePtsX = [CP(:,1) + .5*C WakePtsX ];
    WakePtsY = [CP(:,2)     WakePtsY ];
    WakePtsZ = [CP(:,3)     WakePtsZ ];
    WakeOMEGAx = [Om(:,1) WakeOMEGAx];
    WakeOMEGAy = [Om(:,2) WakeOMEGAy];
    WakeOMEGAz = [Om(:,3) WakeOMEGAz];
    c = sqrt(WakeOMEGAx.^2 + WakeOMEGAy.^2 + WakeOMEGAz.^2);
    %clf
    scatter3(WakePtsX(:),WakePtsY(:),WakePtsZ(:),25, c(:), 'filled');
    %hold on
    %surf(WakePtsX,WakePtsY,WakePtsZ,c);
    view(3)
    %plot(GAMMA);
    axis equal
    drawnow
    disp(i)
end
toc

BlobWake.x = WakePtsX;
BlobWake.y = WakePtsY;
BlobWake.z = WakePtsZ;
BlobWake.Omega.x = WakeOMEGAx;
BlobWake.Omega.y = WakeOMEGAy;
BlobWake.Omega.z = WakeOMEGAz;
save('BlobWake','BlobWake');