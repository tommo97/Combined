function V = WakeVel(X,Y,Z,G)

%%  This does not seem to be any faster.

R(:,:,1) = X;
R(:,:,2) = Y;
R(:,:,3) = Z;

[m n o] = size(R);

Xt(:,:,1) = repmat(X(:),[1 (m-1)*(n-1)])';
Xt(:,:,2) = repmat(Y(:),[1 (m-1)*(n-1)])';
Xt(:,:,3) = repmat(Z(:),[1 (m-1)*(n-1)])';



[m n o] = size(R);


Gma = repmat(G(:),[1 numel(Z)]);
A = {};

A{1} = R(1:end-1, 1:end-1,:);
A{2} = R(1:end-1, 2:end,:);
A{3} = R(2:end,   2:end,:);
A{4} = R(2:end,   1:end-1,:);
A{5} = R(1:end-1, 1:end-1,:);
for i = 1:4
x = A{i}(:,:,1);
y = A{i}(:,:,2);
z = A{i}(:,:,3);
Xs1{i}(:,:,1) = repmat(x(:),[1 numel(X)]);
Xs1{i}(:,:,2) = repmat(y(:),[1 numel(Y)]);
Xs1{i}(:,:,3) = repmat(z(:),[1 numel(Z)]);
x = A{i}(:,:,1);
y = A{i}(:,:,2);
z = A{i}(:,:,3);
Xs2{i}(:,:,1) = repmat(x(:),[1 numel(X)]);
Xs2{i}(:,:,2) = repmat(y(:),[1 numel(Y)]);
Xs2{i}(:,:,3) = repmat(z(:),[1 numel(Z)]);



Vel{i} = PatchLineVel(Xs1{i},Xs2{i},Xt,Gma);

end


V = Vel{1} + Vel{2} + Vel{3} + Vel{4};



% Vel = Vel + PatchLineVel(A2,A3,Target,G);
% Vel = Vel + PatchLineVel(A3,A4,Target,G);
% Vel = Vel + PatchLineVel(A4,A1,Target,G);


% Vx = Vel(:,:,1);
% Vy = Vel(:,:,2);
% Vz = Vel(:,:,3);
% 
% V = [sum(Vx(:)) sum(Vy(:)) sum(Vz(:))];


end


function V = PatchLineVel(Start,End,Target,Gamma)
V = zeros(size(Target));
R1 = Target - Start;
R2 = Target - End;

C(:,:,1) = R1(:,:,2).*R2(:,:,3) - R1(:,:,3).*R2(:,:,2);
C(:,:,2) = R1(:,:,3).*R2(:,:,1) - R1(:,:,1).*R2(:,:,3);
C(:,:,3) = R1(:,:,1).*R2(:,:,2) - R1(:,:,2).*R2(:,:,1);
MagC = sqrt(C(:,:,1).*C(:,:,1) + C(:,:,2).*C(:,:,2) + C(:,:,3).*C(:,:,3));
MagR1 = sqrt(R1(:,:,1).*R1(:,:,1) + R1(:,:,2).*R1(:,:,2) + R1(:,:,3).*R1(:,:,3));
MagR2 = sqrt(R2(:,:,1).*R2(:,:,1) + R2(:,:,2).*R2(:,:,2) + R2(:,:,3).*R2(:,:,3));

R0 = End-Start;
R0R1 = R0(:,:,1).*R1(:,:,1) + R0(:,:,2).*R1(:,:,2) + R0(:,:,3).*R1(:,:,3);
R0R2 = R0(:,:,1).*R2(:,:,1) + R0(:,:,2).*R2(:,:,2) + R0(:,:,3).*R2(:,:,3);
Mult = Gamma./(MagC.*MagC*4*pi);
K = Mult.*((R0R1./MagR1) - (R0R2./MagR2));
K(MagC < 1e-16) = 0;
K(MagR1 < 1e-16) = 0;
K(MagR2 < 1e-16) = 0;
V(:,:,1) = -K.*C(:,:,1);
V(:,:,2) = -K.*C(:,:,2);
V(:,:,3) = -K.*C(:,:,3);
end

