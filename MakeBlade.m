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

y = linspace(Rmin,Rmax,Blade.NSpan)';
%%  Make blade surfaces
UpperS.x = repmat(x,[Blade.NSpan 1]) - Blade.PitchAxis;
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

%%  Close ends
% Blade.Upper.x = [.5*(Upper.x(1,:) + Lower.x(1,:));...
%     Upper.x; .5*(Upper.x(end,:) + Lower.x(end,:))];
% Blade.Upper.y = [.5*(Upper.y(1,:) + Lower.y(1,:));...
%     Upper.y; .5*(Upper.y(end,:) + Lower.y(end,:))];
% Blade.Upper.z = [.5*(Upper.z(1,:) + Lower.z(1,:));...
%     Upper.z; .5*(Upper.z(end,:) + Lower.z(end,:))];
% 
% Blade.Lower.x = [.5*(Upper.x(1,:) + Lower.x(1,:));...
%     Lower.x; .5*(Upper.x(end,:) + Lower.x(end,:))];
% Blade.Lower.y = [.5*(Upper.y(1,:) + Lower.y(1,:));...
%     Lower.y; .5*(Upper.y(end,:) + Lower.y(end,:))];
% Blade.Lower.z = [.5*(Upper.z(1,:) + Lower.z(1,:));...
%     Lower.z; .5*(Upper.z(end,:) + Lower.z(end,:))];
%%  Wrap around
Blade.Surface.Local.z = [fliplr(Blade.Lower.z) Blade.Upper.z(:,2:end-1)];
Blade.Surface.Local.x = [fliplr(Blade.Lower.x) Blade.Upper.x(:,2:end-1)];
Blade.Surface.Local.y = [fliplr(Blade.Lower.y) Blade.Upper.y(:,2:end-1)];

%%  Put into attitude specified by Euler angles
Blade.Surface.Global.x = Blade.Origin(1) + Blade.Surface.Local.x*Blade.TRANS(1,1) + Blade.Surface.Local.y*Blade.TRANS(1,2) + Blade.Surface.Local.z*Blade.TRANS(1,3);
Blade.Surface.Global.y = Blade.Origin(2) + Blade.Surface.Local.x*Blade.TRANS(2,1) + Blade.Surface.Local.y*Blade.TRANS(2,2) + Blade.Surface.Local.z*Blade.TRANS(2,3);
Blade.Surface.Global.z = Blade.Origin(3) + Blade.Surface.Local.x*Blade.TRANS(3,1) + Blade.Surface.Local.y*Blade.TRANS(3,2) + Blade.Surface.Local.z*Blade.TRANS(3,3);

%%  Prepare for export
Blade.X = Blade.Surface.Global.x(:);
Blade.Y = Blade.Surface.Global.y(:);
Blade.Z = Blade.Surface.Global.z(:);
Blade.N.Local = zeros(size(Blade.Surface.Global.x));
Blade.N.Local(:) = 1:numel(Blade.N.Local);
%   Include "missing" panel points
Blade.N.Local = [Blade.N.Local Blade.N.Local(:,1)];
%%   Get panel corner points for a single blade
Blade.Panels.ID.Local = zeros(size(Blade.N.Local) - 1);
Blade.Panels.ID.Local(:) = 1:numel(Blade.Panels.ID.Local);
Blade.Panels.c1.Local = Blade.N.Local(1:end-1,1:end-1);
Blade.Panels.c2.Local = Blade.N.Local(1:end-1,2:end);
Blade.Panels.c3.Local = Blade.N.Local(2:end,2:end);
Blade.Panels.c4.Local = Blade.N.Local(2:end,1:end-1);
Blade.nPnls = numel(Blade.Panels.ID.Local);
Blade.nPts = numel(Blade.X);
Blade.Panels.WakeShedders.LS.Local = Blade.Panels.ID.Local(1:end,1);
Blade.Panels.WakeShedders.US.Local = Blade.Panels.ID.Local(1:end,end);

