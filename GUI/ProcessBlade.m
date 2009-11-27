function Blade = ProcessBlade(Blade)
%%  Make blade surfaces 
%
%   .n is the chordwise position between 0(leading edge) and 1(trailing
%   edge);

TUpperS.x = repmat(Blade.Section.Tip.X,[Blade.NSpan 1]) - Blade.PitchAxis;
TUpperS.y = repmat(Blade.Radius,[1 Blade.NChord]);
TUpperS.z = repmat(Blade.Section.Tip.US,[Blade.NSpan 1]);
TUpperS.n = repmat(Blade.Section.Tip.X,[Blade.NSpan 1]);

RUpperS.x = repmat(Blade.Section.Root.X,[Blade.NSpan 1]) - Blade.PitchAxis;
RUpperS.y = repmat(Blade.Radius,[1 Blade.NChord]);
RUpperS.z = repmat(Blade.Section.Root.US,[Blade.NSpan 1]);
RUpperS.n = repmat(Blade.Section.Root.X,[Blade.NSpan 1]);

TLowerS.x = repmat(Blade.Section.Tip.X,[Blade.NSpan 1]) - Blade.PitchAxis;
TLowerS.y = repmat(Blade.Radius,[1 Blade.NChord]);
TLowerS.z = repmat(Blade.Section.Tip.LS,[Blade.NSpan 1]);
TLowerS.n = repmat(Blade.Section.Tip.X,[Blade.NSpan 1]);

RLowerS.x = repmat(Blade.Section.Root.X,[Blade.NSpan 1]) - Blade.PitchAxis;
RLowerS.y = repmat(Blade.Radius,[1 Blade.NChord]);
RLowerS.z = repmat(Blade.Section.Root.LS,[Blade.NSpan 1]);
RLowerS.n = repmat(Blade.Section.Root.X,[Blade.NSpan 1]);

RootBlendCoefft = repmat(linspace(1,0,Blade.NSpan)', [1 Blade.NChord]);
TipBlendCoefft = 1-RootBlendCoefft;


UpperS.x = RootBlendCoefft.*RUpperS.x + TipBlendCoefft.*TUpperS.x;
UpperS.y = RootBlendCoefft.*RUpperS.y + TipBlendCoefft.*TUpperS.y;
UpperS.z = RootBlendCoefft.*RUpperS.z + TipBlendCoefft.*TUpperS.z;
UpperS.n = RootBlendCoefft.*RUpperS.n + TipBlendCoefft.*TUpperS.n;

LowerS.x = RootBlendCoefft.*RLowerS.x + TipBlendCoefft.*TLowerS.x;
LowerS.y = RootBlendCoefft.*RLowerS.y + TipBlendCoefft.*TLowerS.y;
LowerS.z = RootBlendCoefft.*RLowerS.z + TipBlendCoefft.*TLowerS.z;
LowerS.n = RootBlendCoefft.*RLowerS.n + TipBlendCoefft.*TLowerS.n;
%%   Scale and twist
thetas = repmat(Blade.Theta,[1 Blade.NChord]);
chords = repmat(Blade.Chord,[1 Blade.NChord]);


Blade.Upper.x = chords.*(UpperS.x.*cosd(thetas) - UpperS.z.*sind(thetas));
Blade.Upper.y = UpperS.y;
Blade.Upper.z = chords.*(UpperS.x.*sind(thetas) + UpperS.z.*cosd(thetas));
Blade.Upper.n = UpperS.n;

Blade.Lower.x = chords.*(LowerS.x.*cosd(thetas) - LowerS.z.*sind(thetas));
Blade.Lower.y = LowerS.y;
Blade.Lower.z = chords.*(LowerS.x.*sind(thetas) + LowerS.z.*cosd(thetas));
Blade.Lower.n = LowerS.n;


%%  Close ends - make some caps
d = cos(linspace(0,pi,Blade.NChord-1));

xi = .5*(Blade.Upper.x(1,:) + Blade.Lower.x(1,:));
yi = .5*(Blade.Upper.y(1,:) + Blade.Lower.y(1,:));% - 0.1*sin(linspace(0,pi,Blade.NChord));
zi = .5*(Blade.Upper.z(1,:) + Blade.Lower.z(1,:));
xin = .5*(Blade.Upper.n(1,:) + Blade.Lower.n(1,:));

xo = .5*(Blade.Upper.x(end,:) + Blade.Lower.x(end,:));
yo = .5*(Blade.Upper.y(end,:) + Blade.Lower.y(end,:));% + 0.1*sin(linspace(0,pi,Blade.NChord));
zo = .5*(Blade.Upper.z(end,:) + Blade.Lower.z(end,:));
xon = .5*(Blade.Upper.n(end,:) + Blade.Lower.n(end,:));

xi = [xi(1) xi(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xi)];
xin = [xin(1) xin(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xin)];
xo = [xo(1) xo(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xo)];
xon = [xon(1) xon(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xon)];
% 
% xi = [xi(1) xi(1)+1.25*(xi(2)-xi(1)) xi(3:end-2) xi(end)+1.25*(xi(end-1)-xi(end)) xi(end)];
% xo = [xo(1) xo(1)+1.25*(xo(2)-xo(1)) xo(3:end-2) xo(end)+1.25*(xo(end-1)-xo(end)) xo(end)];
% 

Blade.US.Local.x = [xi;Blade.Upper.x;xo];
Blade.US.Local.y = [yi;Blade.Upper.y;yo];
Blade.US.Local.z = [zi;Blade.Upper.z;zo];
Blade.US.Local.n = [xin;Blade.Upper.n;xon];

Blade.LS.Local.x = [xi;Blade.Lower.x;xo];
Blade.LS.Local.y = [yi;Blade.Lower.y;yo];
Blade.LS.Local.z = [zi;Blade.Lower.z;zo];
Blade.LS.Local.n = [xin;Blade.Lower.n;xon];

%%  Put into attitude specified by Euler angles
Blade.US.Global.x = Blade.US.Local.x;
Blade.US.Global.y = Blade.US.Local.y;
Blade.US.Global.z = Blade.US.Local.z;
Blade.US.Global.n = Blade.US.Local.n;

Blade.LS.Global.x = Blade.LS.Local.x;
Blade.LS.Global.y = Blade.LS.Local.y;
Blade.LS.Global.z = Blade.LS.Local.z;
Blade.LS.Global.n = Blade.LS.Local.n;

%%  Weld seams together -- points first
%   These are the indices
US.N = zeros(size(Blade.US.Global.x));
US.N(:) = 1:numel(US.N);

LS.N = zeros(size(Blade.LS.Global.x));
LS.N(:) = numel(LS.N)+[1:numel(LS.N)];
%   Now get unique points
X = [[Blade.US.Global.x(:);Blade.LS.Global.x(:)] [Blade.US.Global.y(:);Blade.LS.Global.y(:)] [Blade.US.Global.z(:);Blade.LS.Global.z(:)]];
P = [Blade.US.Global.n(:);Blade.LS.Global.n(:)];
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
Blade.n = P(m);

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



% [n1 n2 n3 n4] = fcorner(Blade.N.Local);
% Blade.Panels.c1.Local.n = n1;
% Blade.Panels.c2.Local.n = n2;
% Blade.Panels.c3.Local.n = n3;
% Blade.Panels.c4.Local.n = n4;



Blade.nPnls = numel(Blade.Panels.c1.Local);
Blade.nPts = numel(Blade.X);
Blade.Panels.WakeShedders.LS.Local = MainPans(:,1);
Blade.Panels.WakeShedders.US.Local = MainPans(:,end);


function [Corner1 Corner2 Corner3 Corner4] = fcorner(Pans)
Corner1 = Pans(1:end-1,1:end-1);
Corner2 = Pans(1:end-1,2:end);
Corner3 = Pans(2:end,2:end);
Corner4 = Pans(2:end,1:end-1);