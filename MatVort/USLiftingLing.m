%function USLiftingLing
clear all
close all

tic
NumSteps = 25;
dt = .02;

%%  Define blade geometry
%   Span lies along y axis, chord parallel to x axis
n = 25;
b = 5;
C = 1;
[Xdata Ydata Zdata] = meshgrid(-.5*C:.5*C, linspace(-b,b,n), 0:0);

C1 = [Xdata(1:end-1,1) Ydata(1:end-1,1) Zdata(1:end-1,1)];
C2 = [Xdata(1:end-1,2) Ydata(1:end-1,2) Zdata(1:end-1,2)];
C3 = [Xdata(2:end,2) Ydata(2:end,2) Zdata(2:end,2)];
C4 = [Xdata(2:end,1) Ydata(2:end,1) Zdata(2:end,1)];

Centroid = 0.25*(C1+C2+C3+C4);

CP = Centroid;

D1 = C3 - C1;
D2 = C4 - C2;
CR = cross(D1,D2);
CRmag = sqrt(dot(CR,CR,2));
Norm = CR./[CRmag CRmag CRmag];
Area = 0.5*CRmag;
A = zeros(n-1);
TRG = CP;
NRM = CP;

%%  Some inflow type parameters
Uinf = [10 0 1];
UinfMag = sqrt(dot(Uinf,Uinf,2));
Vinf = [Uinf(1) * ones(n-1,1) Uinf(2) * ones(n-1,1) Uinf(3) * ones(n-1,1)];

%%  Place a vortex ring of about chord length starting at c/4
P1 = 0.5*(0.5*(C1+C2) + C1);
P4 = 0.5*(0.5*(C3+C4) + C4);
Disp = C*Vinf/UinfMag;
Disp(:,1) = 1;
Disp(:,2:3) = 0;
P2 = P1 + Disp;
P3 = P4 + Disp;

%%  Get Influence Matrix

for i = 1:n-1
    TRG(:,1) = CP(i,1);
    TRG(:,2) = CP(i,2);
    TRG(:,3) = CP(i,3);
    NRM(:,1) = Norm(i,1);
    NRM(:,2) = Norm(i,2);
    NRM(:,3) = Norm(i,3);
    %V = LinVel(P1,P4,TRG,1);
    %disp(TRG)
    V = PanVel(P1,P2,P3,P4,TRG,ones(n-1,1));
    
    Vn = dot(NRM,V,2);
    
    A(i,:) = Vn';
    disp(i)
end

%%  Get Initial RHS and solve for GAMMA

RHS = dot(-Vinf,Norm,2);

GAMMA = A\RHS;
disp([RHS GAMMA])
%disp(GAMMA')
%%  Shed a vortex panel

WakePtsX = [];
WakePtsY = [];
WakePtsZ = [];
WakeGAMMA = [];
%%  Shed new wake element
WakePtsX = [[P2(1,1) ; P3(:,1)] WakePtsX ];
WakePtsY = [[P2(1,2) ; P3(:,2)] WakePtsY ];
WakePtsZ = [[P2(1,3) ; P3(:,3)] WakePtsZ ];

counter = 0;
%%  Start of time marching loop
for i = 0:26
    %%  Convect the wake
    %   Get velocities at points due to lifting line and other wake panels
    WakeVelX = Uinf(1) * ones(size(WakePtsX));
    WakeVelY = Uinf(2) * ones(size(WakePtsY));
    WakeVelZ = Uinf(3) * ones(size(WakePtsZ));
    A1x = WakePtsX(1:end-1, 1:end-1,:);
    A2x = WakePtsX(1:end-1, 2:end,:);
    A3x = WakePtsX(2:end,   2:end,:);
    A4x = WakePtsX(2:end,   1:end-1,:);
    A1y = WakePtsY(1:end-1, 1:end-1,:);
    A2y = WakePtsY(1:end-1, 2:end,:);
    A3y = WakePtsY(2:end,   2:end,:);
    A4y = WakePtsY(2:end,   1:end-1,:);
    A1z = WakePtsZ(1:end-1, 1:end-1,:);
    A2z = WakePtsZ(1:end-1, 2:end,:);
    A3z = WakePtsZ(2:end,   2:end,:);
    A4z = WakePtsZ(2:end,   1:end-1,:);
    for j = 1:numel(WakeVelX)
        P = zeros(size(CP));
        P(:,1) = WakePtsX(j);
        P(:,2) = WakePtsY(j);
        P(:,3) = WakePtsZ(j);
        
        V = PanVel(P1,P2,P3,P4,P,GAMMA);
        WakeVelX(j) = Uinf(1) + sum(V(:,1));
        WakeVelY(j) = Uinf(2) + sum(V(:,2));
        WakeVelZ(j) = Uinf(3) + sum(V(:,3));
        P = zeros(numel(A1x),3);
        P(:,1) = WakePtsX(j);
        P(:,2) = WakePtsY(j);
        P(:,3) = WakePtsZ(j);
        %         if i > 0
        PanelVel = PanVel([A1x(:) A1y(:) A1z(:)],[A2x(:) A2y(:) A2z(:)],[A3x(:) A3y(:) A3z(:)],[A4x(:) A4y(:) A4z(:)],P,WakeGAMMA(:));
        V =  sum(PanelVel,1);
        %             V = PatchVel(WakePtsX,WakePtsY,WakePtsZ,[WakePtsX(j) WakePtsY(j) WakePtsZ(j)],WakeGAMMA);
        WakeVelX(j) = WakeVelX(j) + V(1);
        WakeVelY(j) = WakeVelY(j) + V(2);
        WakeVelZ(j) = WakeVelZ(j) + V(3);
        %         end
    end
    WakePtsX = WakePtsX + dt*WakeVelX;
    WakePtsY = WakePtsY + dt*WakeVelY;
    WakePtsZ = WakePtsZ + dt*WakeVelZ;
    
    %     if i == 0
    %         WakePtsX = [[P2(1,1) ; P3(:,1)] WakePtsX ];
    %         WakePtsY = [[P2(1,2) ; P3(:,2)] WakePtsY ];
    %         WakePtsZ = [[P2(1,3) ; P3(:,3)] WakePtsZ ];
    %         WakeGAMMA = GAMMA;
    %     end
    counter = counter + numel(WakePtsX);
    %%  Calculate downwash on blade and solve linear algebra
    V = Vinf;
    for j = 1:size(CP,1)
        P = zeros(numel(A1x),3);
        P(:,1) = CP(j,1);
        P(:,2) = CP(j,2);
        P(:,3) = CP(j,3);
        PanelVel = PanVel([A1x(:) A1y(:) A1z(:)],[A2x(:) A2y(:) A2z(:)],[A3x(:) A3y(:) A3z(:)],[A4x(:) A4y(:) A4z(:)],P,WakeGAMMA(:));
        V(j,:) = V(j,:) + sum(PanelVel,1);
    end
    
    Vn = dot(Norm,-V,2);
    
    GAMMA = A\Vn;
    disp([Vn GAMMA]);
    
    


    %%  Shed new wake element
    WakePtsX = [[P2(1,1) ; P3(:,1)] WakePtsX ];
    WakePtsY = [[P2(1,2) ; P3(:,2)] WakePtsY ];
    WakePtsZ = [[P2(1,3) ; P3(:,3)] WakePtsZ ];
    
    %     for l = 1:n-1
    %         for m = 1:i-1
    %         p1 = [WakePtsX(l,m)  WakePtsY(l,m) WakePtsZ(l,m)];
    %         p2 = [WakePtsX(l+1,m)  WakePtsY(l+1,m) WakePtsZ(l+1,m)];
    %         p3 = [WakePtsX(l+1,m+1)  WakePtsY(l+1,m+1) WakePtsZ(l+1,m+1)];
    %         p4 = [WakePtsX(l,m+1)  WakePtsY(l,m+1) WakePtsZ(l,m+1)];
    %
    %         disp(p1);
    %         disp(p2);
    %         disp(p3);
    %         disp(p4);
    %         V = PanVel(p1,p2,p3,p4,[5.75 5 0.5],WakeGAMMA(l,m));
    %         disp(V)
    %         disp(WakeGAMMA(l,m));
    %         disp(WakeGAMMA(l,m));
    %         disp('---');
    %         end
    %     end
    WakeGAMMA = [GAMMA WakeGAMMA];
    
    disp(i)
    if i > 0
        surf(WakePtsX,WakePtsY,WakePtsZ,WakeGAMMA);
        view(3)
        axis equal
        %         plot(GAMMA);
        drawnow
    end
end
toc
PanelWake.x = WakePtsX;
PanelWake.y = WakePtsY;
PanelWake.z = WakePtsZ;
PanelWake.gamma = WakeGAMMA;
save('PanelWake','PanelWake');
