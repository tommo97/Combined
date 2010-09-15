function Blade = MakeBlade(Blade)
%%  Aerofoil
x = Blade.FOIL.X;
USz = Blade.FOIL.US;
LSz = Blade.FOIL.LS;
Blade.NChord = length(x);
%%   Blade

%   Attitude
EulerAngles = Blade.Attitude;
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
Blade.TRANS(1,:) = [a1 b1 c1];
Blade.TRANS(2,:) = [a2 b2 c2];
Blade.TRANS(3,:) = [a3 b3 c3];


Rmin = min(Blade.RADIUS);
Rmax = max(Blade.RADIUS);


y = BellShape(Rmin,Rmax,Blade.NSpan,5)';
%y = linspace(Rmin,Rmax,Blade.NSpan)'



%%  Make blade surfaces

UpperS.x = repmat(x,[Blade.NSpan 1]);% - Blade.PitchAxis;
UpperS.y = repmat(y,[1 Blade.NChord]);
UpperS.z = repmat(USz,[Blade.NSpan 1]);

LowerS.x = repmat(x,[Blade.NSpan 1]) - Blade.PitchAxis;
LowerS.y = repmat(y,[1 Blade.NChord]);
LowerS.z = repmat(LSz,[Blade.NSpan 1]);

%%   Scale and twist
thetas = Blade.th0 + repmat(interp1(Blade.RADIUS,Blade.THETA,y,'cubic'),[1 Blade.NChord]);
chords = repmat(interp1(Blade.RADIUS,Blade.CHORD,y,'cubic'),[1 Blade.NChord]);

Blade.Upper.x = chords.*(UpperS.x.*cosd(thetas) - UpperS.z.*sind(thetas));
Blade.Upper.y = UpperS.y;
Blade.Upper.z = chords.*(UpperS.x.*sind(thetas) + UpperS.z.*cosd(thetas));

Blade.Lower.x = chords.*(LowerS.x.*cosd(thetas) - LowerS.z.*sind(thetas));
Blade.Lower.y = LowerS.y;
Blade.Lower.z = chords.*(LowerS.x.*sind(thetas) + LowerS.z.*cosd(thetas));

Blade.Upper.x = Blade.Upper.x - min(min(Blade.Upper.x(:),min(Blade.Lower.x(:))));
Blade.Lower.x = Blade.Lower.x - min(min(Blade.Upper.x(:),min(Blade.Lower.x(:))));
%%  Close ends - make some caps
d = cos(linspace(0,pi,Blade.NChord-1));

xi = .5*(Blade.Upper.x(1,:) + Blade.Lower.x(1,:));
yi = .5*(Blade.Upper.y(1,:) + Blade.Lower.y(1,:));% - 0.1*sin(linspace(0,pi,Blade.NChord));
zi = .5*(Blade.Upper.z(1,:) + Blade.Lower.z(1,:));

xo = .5*(Blade.Upper.x(end,:) + Blade.Lower.x(end,:));
yo = .5*(Blade.Upper.y(end,:) + Blade.Lower.y(end,:));% + 0.1*sin(linspace(0,pi,Blade.NChord));
zo = .5*(Blade.Upper.z(end,:) + Blade.Lower.z(end,:));

xi = [xi(1) xi(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xi)];
xo = [xo(1) xo(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xo)];

% 
% xi = [xi(1) xi(1)+1.25*(xi(2)-xi(1)) xi(3:end-2) xi(end)+1.25*(xi(end-1)-xi(end)) xi(end)];
% xo = [xo(1) xo(1)+1.25*(xo(2)-xo(1)) xo(3:end-2) xo(end)+1.25*(xo(end-1)-xo(end)) xo(end)];
% 

Blade.US.Local.x = [Blade.Upper.x;xo];
Blade.US.Local.y = [Blade.Upper.y;yo];
Blade.US.Local.z = [Blade.Upper.z;zo];

Blade.LS.Local.x = [Blade.Lower.x];
Blade.LS.Local.y = [Blade.Lower.y];
Blade.LS.Local.z = [Blade.Lower.z];


%%  Put into attitude specified by Euler angles
Blade.US.Global.x = Blade.Origin(1) + Blade.US.Local.x*Blade.TRANS(1,1) + Blade.US.Local.y*Blade.TRANS(1,2) + Blade.US.Local.z*Blade.TRANS(1,3);
Blade.US.Global.y = Blade.Origin(2) + Blade.US.Local.x*Blade.TRANS(2,1) + Blade.US.Local.y*Blade.TRANS(2,2) + Blade.US.Local.z*Blade.TRANS(2,3);
Blade.US.Global.z = Blade.Origin(3) + Blade.US.Local.x*Blade.TRANS(3,1) + Blade.US.Local.y*Blade.TRANS(3,2) + Blade.US.Local.z*Blade.TRANS(3,3);

Blade.LS.Global.x = Blade.Origin(1) + Blade.LS.Local.x*Blade.TRANS(1,1) + Blade.LS.Local.y*Blade.TRANS(1,2) + Blade.LS.Local.z*Blade.TRANS(1,3);
Blade.LS.Global.y = Blade.Origin(2) + Blade.LS.Local.x*Blade.TRANS(2,1) + Blade.LS.Local.y*Blade.TRANS(2,2) + Blade.LS.Local.z*Blade.TRANS(2,3);
Blade.LS.Global.z = Blade.Origin(3) + Blade.LS.Local.x*Blade.TRANS(3,1) + Blade.LS.Local.y*Blade.TRANS(3,2) + Blade.LS.Local.z*Blade.TRANS(3,3);

%%  Weld seams together -- points first
%   These are the indices
US.N = zeros(size(Blade.US.Global.x));
US.N(:) = 1:numel(US.N);

LS.N = zeros(size(Blade.LS.Global.x));
LS.N(:) = numel(LS.N)+[1:numel(LS.N)];
%   Now get unique points
X = [[Blade.US.Global.x(:);Blade.LS.Global.x(:)] [Blade.US.Global.y(:);Blade.LS.Global.y(:)] [Blade.US.Global.z(:);Blade.LS.Global.z(:)]];
N = [US.N(:);LS.N(:)];
[b, m, n] = unique(X,'rows');
disp('--------')
disp(size(X));
disp(size(unique(X,'rows')));

ind1 = [1:length(m)]';
ind2 = ind1(n);

US.UN = zeros(size(Blade.US.Global.x));
US.UN(:) = ind2(US.N(:));

LS.UN = zeros(size(Blade.LS.Global.x));
LS.UN(:) = ind2(LS.N(:));

US.UN(1,1) = LS.UN(2,2);
LS.UN(1,1) = US.UN(2,2);

US.UN(end,1) = LS.UN(end-1,2);
LS.UN(end,1) = US.UN(end-1,2);


US.UN(1,end) = LS.UN(2,end-1);
LS.UN(1,end) = US.UN(2,end-1);

US.UN(end,end) = LS.UN(end-1,end-1);
LS.UN(end,end) = US.UN(end-1,end-1);

Blade.X = X(m,1);
Blade.Y = X(m,2);
Blade.Z = X(m,3);


%%  Prepare for export
Blade.N.Local = [fliplr(LS.UN(2:end-1,2:end)) US.UN(2:end-1,:)];
%   Tips

Blade.Tip.Inboard.US.N.Local = US.UN(1:2,:);
Blade.Tip.Outboard.US.N.Local = US.UN(end-1:end,:);
Blade.Tip.Inboard.LS.N.Local = LS.UN(1:2,2:end-1);
Blade.Tip.Outboard.LS.N.Local = LS.UN(end-1:end,2:end-1);

%   Now mini closure panel bit

MainPans = zeros(size(Blade.N.Local) - 1);

MainPans(:) = 1:numel(MainPans);

[tc1 tc2 tc3 tc4] = fcorner(Blade.N.Local);
[tiu1 tiu2 tiu3 tiu4] = fcorner(Blade.Tip.Inboard.US.N.Local);
[til4 til3 til2 til1] = fcorner(Blade.Tip.Inboard.LS.N.Local);
[tou1 tou2 tou3 tou4] = fcorner(Blade.Tip.Outboard.US.N.Local);
[tol4 tol3 tol2 tol1] = fcorner(Blade.Tip.Outboard.LS.N.Local);



Blade.Panels.MainPans = MainPans;
Blade.Panels.c1.Local = [tc1(:);tiu1(:);til1(:);tou1(:);tol1(:)];
Blade.Panels.c2.Local = [tc2(:);tiu2(:);til2(:);tou2(:);tol2(:)];
Blade.Panels.c3.Local = [tc3(:);tiu3(:);til3(:);tou3(:);tol3(:)];
Blade.Panels.c4.Local = [tc4(:);tiu4(:);til4(:);tou4(:);tol4(:)];

disp(size(Blade.Panels.c1.Local));
Blade.nPnls = numel(Blade.Panels.c1.Local);
Blade.nPts = numel(Blade.X);
Blade.Panels.WakeShedders.LS.Local = MainPans(:,1);
Blade.Panels.WakeShedders.US.Local = MainPans(:,end);


function [Corner1 Corner2 Corner3 Corner4] = fcorner(Pans)
Corner1 = Pans(1:end-1,1:end-1);
Corner2 = Pans(1:end-1,2:end);
Corner3 = Pans(2:end,2:end);
Corner4 = Pans(2:end,1:end-1);
