function Cp = VelsAndForces(C,Mu,Uinf)

x = C.x;
y = C.y;
z = C.z;


C1 = cat(3,x(1:end-1,1:end-1),y(1:end-1,1:end-1),z(1:end-1,1:end-1));
C2 = cat(3,x(1:end-1,2:end),y(1:end-1,2:end),z(1:end-1,2:end));
C3 = cat(3,x(2:end,2:end),y(2:end,2:end),z(2:end,2:end));
C4 = cat(3,x(2:end,1:end-1),y(2:end,1:end-1),z(2:end,1:end-1));
CP = .25*(C1 + C2 + C3 + C4);

dXi = CP(2:end,:) - CP(1:end-1,:);
dXj = CP(:,2:end) - CP(:,1:end-1);
ri = sqrt(dot(dXi,dXi,3));
rj = sqrt(dot(dXj,dXj,3));

dMui = Mu(2:end,:) - Mu(1:end-1,:);
dMuj = Mu(:,2:end) - Mu(:,1:end-1);




%%  Need to get the freestream/wake/kinematic velocities at the Collocation
%   point into local frames parallel with the span and chordwise directions








[M,N] = size(C.x);
M = M-1;
N = N-1;





K = 0;
Area = zeros(M,N);

A = C3 - C1;
B = C4 - C2;

C = cross(A,B,3);
CMag = sqrt(dot(C,C,3));
Area = CMag/2;
PanelNormal(:,:,1) = -C(:,:,1)./CMag;
PanelNormal(:,:,2) = -C(:,:,2)./CMag;
PanelNormal(:,:,3) = -C(:,:,3)./CMag;



X = CP(:,:,1);
Y = CP(:,:,2);
Z = CP(:,:,3);

% calculation of beta(longitudinal), zeta(transversal) and eta(perpendicular)
% vectors

beta = .5*(C3 + C4 - C1 - C2);
betamag = sqrt(dot(beta,beta,3));
beta(:,:,1) = beta(:,:,1)./betamag
beta(:,:,2) = beta(:,:,2)./betamag
beta(:,:,3) = beta(:,:,3)./betamag

zeta = .5*(C2 + C3 - C1 - C4);
zetamag = sqrt(dot(zeta,zeta,3));
zeta(:,:,1) = zeta(:,:,1)./zetamag;
zeta(:,:,2) = zeta(:,:,2)./zetamag;
zeta(:,:,3) = zeta(:,:,3)./zetamag;

eta = cross(beta,PanelNormal,3);







vx = zeros(M,N);
vy = zeros(M,N);
vz = zeros(M,N);
v = zeros(M,N);
Cp = zeros(M,N);
q = 0.5*1.226*dot(Uinf,Uinf);



for i = 1:M-1
    for j = 1:N-1
        % calculation of induced speeds
        % boundary panels - first order interpolation
        if i == 1
            dx = X(i+1,j)-X(i,j);
            dy = Y(i+1,j)-Y(i,j);
            dz = Z(i+1,j)-Z(i,j);
            dr = sqrt(dx^2+dy^2+dz^2);
            qu = (Mu(i+1,j)-Mu(i,j))/dr;
        end
        if i == M
            dx = X(i,j)-X(i-1,j);
            dy = Y(i,j)-Y(i-1,j);
            dz = Z(i,j)-Z(i-1,j);
            dr = sqrt(dx^2+dy^2+dz^2);
            qu = (Mu(i,j)-Mu(i-1,j))/dr;
        end
        if j == 1
            dx = X(i,j+1)-X(i,j);
            dy = Y(i,j+1)-Y(i,j);
            dz = Z(i,j+1)-Z(i,j);
            dr = sqrt(dx^2+dy^2+dz^2);
            qp = (Mu(i,j+1)-Mu(i,j))/dr;
        end
        if j == N
            dx = X(i,j)-X(i,j-1);
            dy = Y(i,j)-Y(i,j-1);
            dz = Z(i,j)-Z(i,j-1);
            dr = sqrt(dx^2+dy^2+dz^2);
            qp = (Mu(i,j)-Mu(i,j-1))/dr;
        end
        if i ~= 1 && i ~= M
            % inner panels - second order interpolation
            % longitudinal
            dx = X(i+1,j)-X(i,j);
            dy = Y(i+1,j)-Y(i,j);
            dz = Z(i+1,j)-Z(i,j);
            drupwind = sqrt(dx^2+dy^2+dz^2);
            dx = X(i,j)-X(i-1,j);
            dy = Y(i,j)-Y(i-1,j);
            dz = Z(i,j)-Z(i-1,j);
            drdownwind = sqrt(dx^2+dy^2+dz^2);
            qu = (Mu(i+1,j) - Mu(i-1,j))./(drupwind + drdownwind);
        end
        if j ~= 1 && j ~= N
            % transversal
            dx = X(i,j+1)-X(i,j);
            dy = Y(i,j+1)-Y(i,j);
            dz = Z(i,j+1)-Z(i,j);
            drupwind = sqrt(dx^2+dy^2+dz^2);
            dx = X(i,j)-X(i,j-1);
            dy = Y(i,j)-Y(i,j-1);
            dz = Z(i,j)-Z(i,j-1);
            drdownwind = sqrt(dx^2+dy^2+dz^2);
            qp = (Mu(i,j+1) - Mu(i,j-1))./(drupwind + drdownwind);
        end
        
        qo = (zeta(i,j,1)*eta(i,j,1) + zeta(i,j,2)*eta(i,j,2) + zeta(i,j,3)*eta(i,j,3))*qp;
        
        gu = beta(i,j,1)*Uinf(1) + beta(i,j,2)*Uinf(2) + beta(i,j,3)*Uinf(3);
        go = eta(i,j,1)*Uinf(1) + eta(i,j,2)*Uinf(2) + eta(i,j,3)*Uinf(3);
        
        vx(i,j) = (-qu + gu)*beta(i,j,1) + (-qo + go)*eta(i,j,1);
        vy(i,j) = (-qu + gu)*beta(i,j,2) + (-qo + go)*eta(i,j,2);
        vz(i,j) = (-qu + gu)*beta(i,j,3) + (-qo + go)*eta(i,j,3);
        v(i,j) = realsqrt(vx(i,j)^2+vy(i,j)^2+vz(i,j)^2);
        
        Cp(i,j) = 1-v(i,j)^2/dot(Uinf,Uinf);

    end
end