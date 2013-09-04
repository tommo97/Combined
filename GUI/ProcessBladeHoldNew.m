function Blade = ProcessBlade(Blade)
%%  Make blade surfaces
%
%   .n is the chordwise position between 0(leading edge) and 1(trailing
%   edge);
clc

closeTips = true;
disp(Blade);


SpanRepSize = [Blade.NSpan 1];
ChordRepSize = [1 Blade.NChord];

TUpperS.x = repmat(Blade.Section.Tip.X,SpanRepSize) - Blade.PitchAxis;
TUpperS.y = repmat(Blade.Radius,ChordRepSize);
TUpperS.z = repmat(Blade.Section.Tip.US,SpanRepSize);
TUpperS.n = repmat(Blade.Section.Tip.X,SpanRepSize);
TUpperS.m = repmat(Blade.Radius,ChordRepSize)./max(Blade.Radius(:));

RUpperS.x = repmat(Blade.Section.Root.X,SpanRepSize) - Blade.PitchAxis;
RUpperS.y = repmat(Blade.Radius,ChordRepSize);
RUpperS.z = repmat(Blade.Section.Root.US,SpanRepSize);
RUpperS.n = repmat(Blade.Section.Root.X,SpanRepSize);
RUpperS.m = repmat(Blade.Radius,ChordRepSize)./max(Blade.Radius(:));




TLowerS.x = repmat(Blade.Section.Tip.X,SpanRepSize) - Blade.PitchAxis;
TLowerS.y = repmat(Blade.Radius,ChordRepSize);
TLowerS.z = repmat(Blade.Section.Tip.LS,SpanRepSize);
TLowerS.n = repmat(Blade.Section.Tip.X,SpanRepSize);
TLowerS.m = repmat(Blade.Radius,ChordRepSize)./max(Blade.Radius(:));


RLowerS.x = repmat(Blade.Section.Root.X,SpanRepSize) - Blade.PitchAxis;
RLowerS.y = repmat(Blade.Radius,ChordRepSize);
RLowerS.z = repmat(Blade.Section.Root.LS,SpanRepSize);
RLowerS.n = repmat(Blade.Section.Root.X,SpanRepSize);
RLowerS.m = repmat(Blade.Radius,ChordRepSize)./max(Blade.Radius(:));


% if (Blade.Reverse)
%     TUpperS.x = -TUpperS.x;
%     RUpperS.x = -RUpperS.x;
%     TLowerS.x = -TLowerS.x;
%     RLowerS.x = -RLowerS.x;
%
% end




RootBlendCoefft = repmat(linspace(1,0,Blade.NSpan)', ChordRepSize);

if ~isempty(Blade.Thickness)
    th = Blade.Thickness - min(Blade.Thickness);
    RootBlendCoefft = repmat(th./max(th)', ChordRepSize);
end


if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
    %RootBlendCoefft = repmat(Blade.TransitionPiece', ChordRepSize);
    if isempty(RootBlendCoefft)
        RootBlendCoefft = repmat(linspace(1,0,Blade.NSpan)', ChordRepSize);
    end
    
  
end



% if Blade.isSOTON
%     RootBlendCoefft = repmat((Blade.Thickness - Blade.Thickness(end)), ChordRepSize)/(Blade.Thickness(1) - Blade.Thickness(end));
% end

if (~isempty(Blade.Thickness))
    RootBlendCoefft = repmat((Blade.Thickness - Blade.Thickness(end)), ChordRepSize)/(Blade.Thickness(1) - Blade.Thickness(end));
end


TipBlendCoefft = 1-RootBlendCoefft;

UpperS.x = RootBlendCoefft.*RUpperS.x + TipBlendCoefft.*TUpperS.x;
UpperS.y = RootBlendCoefft.*RUpperS.y + TipBlendCoefft.*TUpperS.y;
UpperS.z = RootBlendCoefft.*RUpperS.z + TipBlendCoefft.*TUpperS.z;
UpperS.n = RootBlendCoefft.*RUpperS.n + TipBlendCoefft.*TUpperS.n;
UpperS.m = RootBlendCoefft.*RUpperS.m + TipBlendCoefft.*TUpperS.m;

LowerS.x = RootBlendCoefft.*RLowerS.x + TipBlendCoefft.*TLowerS.x;
LowerS.y = RootBlendCoefft.*RLowerS.y + TipBlendCoefft.*TLowerS.y;
LowerS.z = RootBlendCoefft.*RLowerS.z + TipBlendCoefft.*TLowerS.z;
LowerS.n = RootBlendCoefft.*RLowerS.n + TipBlendCoefft.*TLowerS.n;
LowerS.m = RootBlendCoefft.*RLowerS.m + TipBlendCoefft.*TLowerS.m;



if ~isempty(Blade.Thickness)
    chords = 1;
    LocalThicknesses = repmat(max(UpperS.z' - LowerS.z')',ChordRepSize);
    DesiredThicknesses = repmat(Blade.Thickness,ChordRepSize);
    LowerS.z = LowerS.z.*DesiredThicknesses./LocalThicknesses;
    UpperS.z = UpperS.z.*DesiredThicknesses./LocalThicknesses;
end


if ~isempty(Blade.Camber)
    MeanLineTip = 0.5*(UpperS.z + LowerS.z) + Blade.Camber  * Blade.Section.Tip.Camber;
    MeanLineRoot = 0.5*(UpperS.z + LowerS.z) + Blade.Camber  * Blade.Section.Root.Camber;
    ThicknessDist = (UpperS.z - LowerS.z);
    UpperS.z = (RootBlendCoefft.*MeanLineRoot + TipBlendCoefft.*MeanLineTip) + 0.5*ThicknessDist;
    LowerS.z = (RootBlendCoefft.*MeanLineRoot + TipBlendCoefft.*MeanLineTip) - 0.5*ThicknessDist;
end
% 
% %%  Close ends - make some caps
% if ~Blade.RoundTips
%     d = cos(linspace(0,pi,Blade.NChord-1));
%     
%     xi = .5*(UpperS.x(1,:) + LowerS.x(1,:));
%     yi = .5*(UpperS.y(1,:) + LowerS.y(1,:));% - 0.1*sin(linspace(0,pi,Blade.NChord));
%     zi = .5*(UpperS.z(1,:) + LowerS.z(1,:));
%     xin = .5*(UpperS.n(1,:) + LowerS.n(1,:));
%     xim = .5*(UpperS.m(1,:) + LowerS.m(1,:));
%     
%     xo = .5*(UpperS.x(end,:) + LowerS.x(end,:));
%     yo = .5*(UpperS.y(end,:) + LowerS.y(end,:));% + 0.1*sin(linspace(0,pi,Blade.NChord));
%     zo = .5*(UpperS.z(end,:) + LowerS.z(end,:));
%     xon = .5*(UpperS.n(end,:) + LowerS.n(end,:));
%     xom = .5*(UpperS.m(end,:) + LowerS.m(end,:));
%     
%     xi = [xi(1) xi(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xi)];
%     xin = [xin(1) xin(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xin)];
%     xo = [xo(1) xo(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xo)];
%     xon = [xon(1) xon(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xon)];
%     xom = [xom(1) xom(2:end)+linspace(1,0,Blade.NChord-1).*d.*diff(xom)];
%     %
%     % xi = [xi(1) xi(1)+1.25*(xi(2)-xi(1)) xi(3:end-2) xi(end)+1.25*(xi(end-1)-xi(end)) xi(end)];
%     % xo = [xo(1) xo(1)+1.25*(xo(2)-xo(1)) xo(3:end-2) xo(end)+1.25*(xo(end-1)-xo(end)) xo(end)];
%     %
%     
%     if (isempty(Blade.Thickness)) || (Blade.Thickness(end) > 0) %%  If outer thickness is > 0 then use caps
%         UpperS.x = [xi;UpperS.x;xo];
%         UpperS.y = [yi;UpperS.y;yo];
%         UpperS.z = [zi;UpperS.z;zo];
%         UpperS.n = [xin;UpperS.n;xon];
%         UpperS.m = [xim;UpperS.m;xom];
%         
%         LowerS.x = [xi;LowerS.x;xo];
%         LowerS.y = [yi;LowerS.y;yo];
%         LowerS.z = [zi;LowerS.z;zo];
%         LowerS.n = [xin;LowerS.n;xon];
%         LowerS.m = [xim;LowerS.m;xom];
%     else     %    just seal edge
%         closeTips = false;
%         UpperS.x = [xi;UpperS.x];
%         UpperS.y = [yi;UpperS.y];
%         UpperS.z = [zi;UpperS.z];
%         UpperS.n = [xin;UpperS.n];
%         UpperS.m = [xim;UpperS.m];
%         
%         LowerS.x = [xi;LowerS.x];
%         LowerS.y = [yi;LowerS.y];
%         LowerS.z = [zi;LowerS.z];
%         LowerS.n = [xin;LowerS.n];
%         LowerS.m = [xim;LowerS.m];
%     end
% 
%     
%     
%     
%     
%     
%     
% else
%     
%     
%     %   first we need to make a blending coefficient matrix for upper and lower
%     %   surfaces
%     %Blade.num_tip_pans = 3;
%     num_pts_j = 2*(Blade.num_tip_pans);
%     num_pts_i = size(UpperS.x(1,:),2);
%     
%     
%     blendi = repmat(linspace(0,1,num_pts_j)',[1,num_pts_i]);
%     
%     blendi = flipud(blendi(2:end-1,:));
%     
%     %   Minus below since lower surface is negative and we want the mean...
%     profilei = blendi.*repmat(UpperS.z(1,:),[num_pts_j-2,1]) - (flipud(blendi)).*repmat(LowerS.z(1,:),[num_pts_j-2,1]);
%     
%     profileo = blendi.*repmat(UpperS.z(end,:),[num_pts_j-2,1]) - (flipud(blendi)).*repmat(LowerS.z(end,:),[num_pts_j-2,1]);
%     
%     
%     
%     %   Now need to get the angle of each new surface
%     
%     ZMI = 0.5*(UpperS.z(1,:) + LowerS.z(1,:));
%     ZMO = 0.5*(UpperS.z(end,:) + LowerS.z(end,:));
%     
%     XMI = UpperS.x(1,:);
%     XMO = UpperS.x(end,:);
%     
%     YMI = 0.5*(UpperS.y(1,:) + LowerS.y(1,:)) - 0.5*(UpperS.z(1,:) - LowerS.z(1,:));
%     YMO = 0.5*(UpperS.y(end,:) + LowerS.y(end,:)) + 0.5*(UpperS.z(end,:) - LowerS.z(end,:));
%     
%     
%     
%     thetas = linspace(0,-90,1+num_pts_j/2);
%     thetas = thetas(2:end-1);
%     
%     sz = size(thetas,2);
%     
%     ri = profilei(1:sz,:);
%     ti = repmat(thetas',[1,num_pts_i]);
%     xi = repmat(UpperS.x(1,:),[sz,1]);
%     [YI,ZI,XI] = pol2cart(deg2rad(ti),ri,xi);
%     YI = repmat(UpperS.y(1,:),[sz,1]) - YI;
%     
%     
%     thetas = linspace(-90,0,1+num_pts_j/2);
%     thetas = thetas(2:end-1);
%     ro = profileo(1:sz,:);
%     to = repmat(thetas',[1,num_pts_i]);
%     xo = repmat(UpperS.x(end,:),[sz,1]);
%     [YO,ZO,XO] = pol2cart(deg2rad(to),ro,xo);
%     YO = repmat(UpperS.y(end,:),[sz,1]) + YO;
%     
%     
%     
%     UpperS.x = [XMI;XI;UpperS.x;XO;XMO];
%     UpperS.y = [YMI;YI;UpperS.y;YO;YMO];
%     UpperS.z = [ZMI;repmat(ZMI,[sz,1])-ZI;UpperS.z;repmat(ZMO,[sz,1])-ZO;ZMO];
%     
%     
%     
%     thetas = linspace(0,90,1+num_pts_j/2);
%     thetas = thetas(2:end-1);
%     
%     sz = size(thetas,2);
%     
%     ri = flipud(profilei((sz+1):end,:));
%     ti = repmat(thetas',[1,num_pts_i]);
%     xi = repmat(LowerS.x(1,:),[sz,1]);
%     [YI,ZI,XI] = pol2cart(deg2rad(ti),ri,xi);
%     YI = repmat(LowerS.y(1,:),[sz,1]) - YI;
%     
%     
%     thetas = linspace(90,0,1+num_pts_j/2);
%     thetas = thetas(2:end-1);
%     ro = profileo((sz+1):end,:);
%     to = repmat(thetas',[1,num_pts_i]);
%     xo = repmat(LowerS.x(end,:),[sz,1]);
%     [YO,ZO,XO] = pol2cart(deg2rad(to),ro,xo);
%     YO = repmat(LowerS.y(end,:),[sz,1]) + YO;
%     
%     
%     
%     LowerS.x = [XMI;XI;LowerS.x;XO;XMO];
%     LowerS.y = [YMI;YI;LowerS.y;YO;YMO];
%     LowerS.z = [ZMI;repmat(ZMI,[sz,1])-ZI;LowerS.z;repmat(ZMO,[sz,1])-ZO;ZMO];
%     
%     
%     
%     UpperS.m = [repmat(UpperS.m(1,:),[sz+1,1]); UpperS.m; repmat(UpperS.m(end,:),[sz+1,1])];
%     UpperS.n = [repmat(UpperS.n(1,:),[sz+1,1]); UpperS.n; repmat(UpperS.n(end,:),[sz+1,1])];
%     LowerS.m = [repmat(LowerS.m(1,:),[sz+1,1]); LowerS.m; repmat(LowerS.m(end,:),[sz+1,1])];
%     LowerS.n = [repmat(LowerS.n(1,:),[sz+1,1]); LowerS.n; repmat(LowerS.n(end,:),[sz+1,1])];
%     
% %     figure; surf(UpperS.x, UpperS.y, UpperS.z); axis equal; view(3)
% %     hold all
% %     surf(LowerS.x, LowerS.y, LowerS.z); axis equal; view(3)
%     
% 
% end


rakes = zeros(size(UpperS.x));
skews = rakes;
thetas = rakes;
chords = rakes;

skews(:) = interp1(Blade.Radius,Blade.Skew,UpperS.y(:),'cubic','extrap');
rakes(:) = interp1(Blade.Radius,Blade.Rake,UpperS.y(:),'cubic','extrap');
thetas(:) = interp1(Blade.Radius,Blade.Theta,UpperS.y(:),'cubic','extrap');
chords(:) = interp1(Blade.Radius,Blade.Chord,UpperS.y(:),'cubic','extrap');

UpperS.x = chords.*UpperS.x;
UpperS.y = UpperS.y;
UpperS.z = chords.*UpperS.z;

LowerS.x = chords.*LowerS.x;
LowerS.y = LowerS.y;
LowerS.z = chords.*LowerS.z;


%   A circle of same chord as root, centered at the root
if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
LXs = LowerS.x(1,:);
UXs = UpperS.x(1,:);
CX = 0.5*(min(LXs) + 0.5*(max(LXs) - min(LXs)) +  min(UXs) + 0.5*(max(UXs) - min(UXs)));
RADIUS = 0.25*((max(LXs) - min(LXs)) + (max(UXs) - min(UXs)));

xu = UXs - CX;
xl = LXs - CX;
zu = sqrt(max(RADIUS*RADIUS - xu.*xu,0));
zl = -sqrt(max(RADIUS*RADIUS - xl.*xl,0));


UZs = repmat(zu,[Blade.n2+1,1]);
LZs = repmat(zl,[Blade.n2+1,1]);

[a b ] = size(LowerS.y(1:(Blade.n2+1),:))
scales = repmat(linspace(1,0,a),[b 1]);
scales = repmat([ones(1,Blade.n1) linspace(1,0,a-Blade.n1)],[b 1]);

URADII =  UpperS.x(1:(Blade.n2+1),:) - CX;
LRADII =  LowerS.x(1:(Blade.n2+1),:) - CX;
uinds = find(abs(URADII) <= RADIUS);
linds = find(abs(LRADII) <= RADIUS);

tu = UpperS.z(1:(Blade.n2+1),:);
tl = LowerS.z(1:(Blade.n2+1),:);

tu(uinds) = max(tu(uinds),UZs(uinds));
tl(linds) = min(tl(linds)+LZs(linds));


UpperS.z(1:(Blade.n2+1),:) = scales.^0.5'.*UZs + (1-scales.^0.5').*UpperS.z(1:(Blade.n2+1),:);

LowerS.z(1:(Blade.n2+1),:) = scales.^0.5'.*LZs + (1-scales.^0.5').*LowerS.z(1:(Blade.n2+1),:);

Zs = 0.5*(LowerS.z(1,:) + UpperS.z(1,:));

LowerS.z(1,:) = Zs;
UpperS.z(1,:) = Zs;
end
% 
% CircLower.z = -repmat(sqrt(0.25 - (Blade.Section.Root.X - 0.5).^2),SpanRepSize);
% CircUpper.z = repmat(sqrt(0.25 - (Blade.Section.Root.X - 0.5).^2),SpanRepSize) ;
% if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
%     Blade.CircSection = zeros(size(Blade.TransitionPiece));
%     Blade.CircSection(Blade.n1:Blade.n2) = linspace(1,0,Blade.n2 - Blade.n1 +1);
%     Blade.CircSection(1:Blade.n1) = 1;
%     CircSectionBlendCoefft = repmat(Blade.CircSection', ChordRepSize);
%     
%     NRELscaleR = [0 0.508 0.660 0.883 1.008  1.067 1.133 1.257 1.343 5.532];
%     NRELscaleS = [0.218/0.218 0.218/0.218 0.218/0.218 0.183/0.183 0.163/0.349 0.154/0.442 0.154/0.544 0.154/0.738 0.2095 0.2095];%/0.2095;
%     
%   %if Blade.isNREL
%   %    Supper = interp1(NRELscaleR,NRELscaleS,UpperS.y);
%   %    Slower = interp1(NRELscaleR,NRELscaleS,LowerS.y);
%   %else
%       Supper = ones(size(UpperS.y));
%       Slower = ones(size(UpperS.y));
%   %end
%     
%     UpperS.z(2:end-1,:) = CircSectionBlendCoefft.*CircUpper.z + (1-CircSectionBlendCoefft).*UpperS.z(2:end-1,:);
%     UpperS.z = UpperS.z.*Supper;
%     LowerS.z(2:end-1,:) = CircSectionBlendCoefft.*CircLower.z + (1-CircSectionBlendCoefft).*LowerS.z(2:end-1,:);
%     LowerS.z = LowerS.z.*Slower;
%     
% end
% 


%%   Scale and twist



%UpperS.z(1,:)./sqrt(UpperS.x(1,:).^2 + UpperS.z(1,:).^2);
%LowerS.z(1,:)./sqrt(LowerS.x(1,:).^2 + LowerS.z(1,:).^2);
   
            
utx = (UpperS.x.*cosd(thetas) - UpperS.z.*sind(thetas));
utz = (UpperS.x.*sind(thetas) + UpperS.z.*cosd(thetas)) - 0.0;%414;
uty = UpperS.y;

ltx = (LowerS.x.*cosd(thetas) - LowerS.z.*sind(thetas));
ltz = (LowerS.x.*sind(thetas) + LowerS.z.*cosd(thetas)) - 0.0;%414;
lty = LowerS.y;

%%  Now, if a prop then need to bend around the hub - this might wreck the 
%   skew and rake angles, so check carefully the correct order of this. For
%   the 4119 prop there is no skew nor rake so is ok

if (Blade.isProp)
    if (closeTips)
        R = repmat(Blade.Radius([ 1:end ]),[1 size(utx,2)]);
    else
        R = repmat(Blade.Radius([ 1:end]),[1 size(utx,2)]);
    end
    thta = atan(utx./(uty + 1e-16));
    utx = R.*sin(thta);
    uty = R.*cos(thta);
    
    if (closeTips)
        R = repmat(Blade.Radius([ 1:end ]),[1 size(ltx,2)]);
    else
        R = repmat(Blade.Radius([ 1:end]),[1 size(ltx,2)]);
    end
    thta = atan(ltx./(lty + 1e-16));
    ltx = R.*sin(thta);
    lty = R.*cos(thta);
end


Nhub = 5;

ty = uty(1,:);
tx = utx(1,:);
thta1 = atan(tx./ty);
rs1 = sqrt(tx.*tx + ty.*ty);

ty = lty(1,:);
tx = ltx(1,:);
thta2 = atan(tx./ty);
rs2 = sqrt(tx.*tx + ty.*ty);
thta_mean = 0.5*(thta1 + thta2);
rs = repmat(0.5*(rs1 + rs2), [Nhub 1]);

mult = repmat(linspace(0,1,Nhub + 1)',[1 size(thta1,2) ]);
mult = mult(1:Nhub,:);

thta = repmat(thta1,[Nhub 1]).*mult + (1-mult).*(repmat(thta_mean,[Nhub 1])-pi/3);


utx = [rs .* sin(thta);utx];
uty = [rs .* cos(thta);uty];

thta = repmat(thta2,[Nhub 1]).*mult + (1-mult).*(repmat(thta_mean,[Nhub 1]) + pi/3);
ltx = [rs .* sin(thta);ltx];
lty = [rs .* cos(thta);lty];

tz = 0.5*(utz(1,:) + ltz(1,:));

tzu = repmat(utz(1,:),[Nhub 1]).*mult + (1-mult).*(repmat(tz,[Nhub 1]));
utz = [tzu;utz];
tzl = repmat(ltz(1,:),[Nhub 1]).*mult + (1-mult).*(repmat(tz,[Nhub 1]));
ltz = [tzl;ltz];


thetas = zeros(size(utx));
skews = thetas;
rakes = thetas;
thetas(:) = interp1(Blade.Radius,Blade.Theta,uty,'cubic','extrap');
ty = [repmat(UpperS.y(1,:),[Nhub 1]);UpperS.y];
skews(:) = interp1(Blade.Radius,Blade.Skew,ty,'cubic','extrap');
rakes(:) = interp1(Blade.Radius,Blade.Rake,ty,'cubic','extrap');
xin = .5*(UpperS.n(1,:) + LowerS.n(1,:));
xim = .5*(UpperS.m(1,:) + LowerS.m(1,:));
xn = xin;
        xm = xim;
        for i = 1:Nhub
            xin = [xn;xin];
            xim = [xm;xim];
        end
        UpperS.n = [xin;UpperS.n];
        UpperS.m = [xim;UpperS.m];
        
     
        LowerS.n = [xin;LowerS.n];
        LowerS.m = [xim;LowerS.m];

%%  Put points into attitude specified by Skew angles


cosphi = cos(0); costhe = cos(0); cospsi = cos(skews);
sinphi = sin(0); sinthe = sin(0); sinpsi = sin(skews);
a1 = costhe*cospsi;
a2 = costhe*sinpsi;
a3 = -sinthe;
b1 = sinphi*sinthe*cospsi - cosphi*sinpsi;
b2 = sinphi*sinthe*sinpsi + cosphi*cospsi;
b3 = sinphi*costhe;
c1 = cosphi*sinthe*cospsi + sinphi*sinpsi;
c2 = cosphi*sinthe*sinpsi - sinphi*cospsi;
c3 = cosphi*costhe;


UpperS.x = utx.*a1 + uty.*a2 + utz.*a3;
UpperS.y = utx.*b1 + uty.*b2 + utz.*b3;
UpperS.z = utx.*c1 + uty.*c2 + utz.*c3;

LowerS.x = ltx.*a1 + lty.*a2 + ltz.*a3;
LowerS.y = ltx.*b1 + lty.*b2 + ltz.*b3;
LowerS.z = ltx.*c1 + lty.*c2 + ltz.*c3;

Blade.Upper = UpperS;
Blade.US.Local = UpperS;
Blade.Lower = LowerS;
Blade.LS.Local = LowerS;

clear Blade.Upper Blade.Lower Blade.US.Local Blade.LS.Local


Blade.US.Global.x = Blade.US.Local.x;
Blade.US.Global.y = Blade.US.Local.y;
Blade.US.Global.z = rakes + Blade.US.Local.z;
Blade.US.Global.n = Blade.US.Local.n;
Blade.US.Global.m = Blade.US.Local.m;

Blade.LS.Global.x = Blade.LS.Local.x;
Blade.LS.Global.y = Blade.LS.Local.y;
Blade.LS.Global.z = rakes + Blade.LS.Local.z;
Blade.LS.Global.n = Blade.LS.Local.n;
Blade.LS.Global.m = Blade.LS.Local.m;


%%  Weld seams together -- points first
%   These are the indices
US.N = zeros(size(Blade.US.Global.x));
US.N(:) = 1:numel(US.N);

LS.N = zeros(size(Blade.LS.Global.x));
LS.N(:) = numel(LS.N)+[1:numel(LS.N)];
%   Now get unique points
X = [[Blade.US.Global.x(:);Blade.LS.Global.x(:)] [Blade.US.Global.y(:);Blade.LS.Global.y(:)] [Blade.US.Global.z(:);Blade.LS.Global.z(:)]];
P = [Blade.US.Global.n(:);Blade.LS.Global.n(:)];
Q = [Blade.US.Global.m(:);Blade.LS.Global.m(:)];
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


if ~Blade.RoundTips
    US.UN(1,1) = LS.UN(2,2);
    LS.UN(1,1) = US.UN(2,2);
    
    US.UN(end,1) = LS.UN(end-1,2);
    LS.UN(end,1) = US.UN(end-1,2);
    
    
    US.UN(1,end) = LS.UN(2,end-1);
    LS.UN(1,end) = US.UN(2,end-1);
    
    US.UN(end,end) = LS.UN(end-1,end-1);
    LS.UN(end,end) = US.UN(end-1,end-1);
    Blade.N.Local = [fliplr(LS.UN(2:end-1,2:end)) US.UN(2:end-1,:)];
else
    Blade.N.Local = [fliplr(LS.UN(:,2:end)) US.UN(:,:)];
end
Blade.X = X(m,1);
Blade.Y = X(m,2);
Blade.Z = X(m,3);
Blade.n = P(m);
Blade.m = Q(m);

%%  Prepare for export

%   Tips
if ~Blade.RoundTips
    Blade.Tip.Inboard.US.N.Local = US.UN(1:2,:);
    Blade.Tip.Outboard.US.N.Local = US.UN(end-1:end,:);
    Blade.Tip.Inboard.LS.N.Local = LS.UN(1:2,2:end-1);
    Blade.Tip.Outboard.LS.N.Local = LS.UN(end-1:end,2:end-1);
else
    Blade.Tip.Inboard.US.N.Local = [];
    Blade.Tip.Outboard.US.N.Local = [];
    Blade.Tip.Inboard.LS.N.Local = [];
    Blade.Tip.Outboard.LS.N.Local = [];
end
%   Now mini closure panel bit

MainPans = zeros(size(Blade.N.Local) - 1);

MainPans(:) = 1:numel(MainPans);
[tc1 tc2 tc3 tc4] = fcorner(Blade.N.Local);
if ~Blade.RoundTips
    
    [tiu1 tiu2 tiu3 tiu4] = fcorner(Blade.Tip.Inboard.US.N.Local);
    [til4 til3 til2 til1] = fcorner(Blade.Tip.Inboard.LS.N.Local);
    [tou1 tou2 tou3 tou4] = fcorner(Blade.Tip.Outboard.US.N.Local);
    [tol4 tol3 tol2 tol1] = fcorner(Blade.Tip.Outboard.LS.N.Local);
end


Blade.Panels.MainPans = MainPans;
if ~Blade.RoundTips
    Blade.Panels.c1.Local = [tc1(:);tiu1(:);til1(:);tou1(:);tol1(:)];
    Blade.Panels.c2.Local = [tc2(:);tiu2(:);til2(:);tou2(:);tol2(:)];
    Blade.Panels.c3.Local = [tc3(:);tiu3(:);til3(:);tou3(:);tol3(:)];
    Blade.Panels.c4.Local = [tc4(:);tiu4(:);til4(:);tou4(:);tol4(:)];
else
    
    Blade.Panels.c1.Local = tc1(:);
    Blade.Panels.c2.Local = tc2(:);
    Blade.Panels.c3.Local = tc3(:);
    Blade.Panels.c4.Local = tc4(:);
    
end


% [n1 n2 n3 n4] = fcorner(Blade.N.Local);
% Blade.Panels.c1.Local.n = n1;
% Blade.Panels.c2.Local.n = n2;
% Blade.Panels.c3.Local.n = n3;
% Blade.Panels.c4.Local.n = n4;



Blade.nPnls = numel(Blade.Panels.c1.Local);
Blade.nPts = numel(Blade.X);
Blade.Panels.WakeShedders.LS.Local = MainPans(:,1);
Blade.Panels.WakeShedders.US.Local = MainPans(:,end);

if Blade.RoundTips
    n = Blade.num_tip_pans+1;
    Blade.Panels.WakeShedders.LS.Local = MainPans(n:(end-n),1);
    Blade.Panels.WakeShedders.US.Local = MainPans(n:(end-n),end);
    if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
            Blade.Panels.WakeShedders.LS.Local = MainPans((n+Blade.n2):(end-n),1);
        Blade.Panels.WakeShedders.US.Local = MainPans((n+Blade.n2):(end-n),end);
    end
end
    
if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
    Blade.Panels.WakeShedders.LS.Local = MainPans(Blade.n2:end,1);
    Blade.Panels.WakeShedders.US.Local = MainPans(Blade.n2:end,end);
end
Mp = zeros(size(Blade.N.Local) - 1);
Mp(:) = 1:numel(Mp);
if ~Blade.RoundTips
    t1 = zeros(size(Blade.Tip.Inboard.US.N.Local) - 1); t1(:) = 1:numel(t1); t1 = t1 + max(Mp(:));
    t2 = zeros(size(Blade.Tip.Inboard.LS.N.Local) - 1); t2(:) = 1:numel(t2); t2 = t2 + max(t1(:));
    t3 = zeros(size(Blade.Tip.Outboard.US.N.Local) - 1); t3(:) = 1:numel(t3); t3 = t3 + max(t2(:));
    t4 = zeros(size(Blade.Tip.Outboard.LS.N.Local) - 1); t4(:) = 1:numel(t4); t4 = t4 + max(t3(:));
end
Blade.Panels.MainSurf = Mp;
if ~Blade.RoundTips
    Blade.Panels.TipInnerUS = t1;
    Blade.Panels.TipInnerLS = t2;
    Blade.Panels.TipOuterUS = t3;
    Blade.Panels.TipOuterLS = t4;
else
    Blade.Panels.TipInnerUS = [];
    Blade.Panels.TipInnerLS = [];
    Blade.Panels.TipOuterUS = [];
    Blade.Panels.TipOuterLS = [];
end

function [Corner1 Corner2 Corner3 Corner4] = fcorner(Pans)
Corner1 = Pans(1:end-1,1:end-1);
Corner2 = Pans(1:end-1,2:end);
Corner3 = Pans(2:end,2:end);
Corner4 = Pans(2:end,1:end-1);