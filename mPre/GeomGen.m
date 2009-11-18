function Blade = GeomGen(Blade)


Blade = MakeBlade(Blade);

S = zeros(size(Blade.X(Blade.N.Local)));
for j = 1:size(S,2)
    S(:,j) = j;
end

surf(Blade.axes,Blade.X(Blade.N.Local),Blade.Y(Blade.N.Local),Blade.Z(Blade.N.Local),S);
hold(Blade.axes,'on');
surf(Blade.axes,Blade.X(Blade.Tip.Inboard.US.N.Local),Blade.Y(Blade.Tip.Inboard.US.N.Local),Blade.Z(Blade.Tip.Inboard.US.N.Local));
surf(Blade.axes,Blade.X(Blade.Tip.Inboard.LS.N.Local),Blade.Y(Blade.Tip.Inboard.LS.N.Local),Blade.Z(Blade.Tip.Inboard.LS.N.Local));
surf(Blade.axes,Blade.X(Blade.Tip.Outboard.US.N.Local),Blade.Y(Blade.Tip.Outboard.US.N.Local),Blade.Z(Blade.Tip.Outboard.US.N.Local));
surf(Blade.axes,Blade.X(Blade.Tip.Outboard.LS.N.Local),Blade.Y(Blade.Tip.Outboard.LS.N.Local),Blade.Z(Blade.Tip.Outboard.LS.N.Local));
%
hold(Blade.axes,'off');
% Blade.N.Global = Blade.N.Local + NumPoints;
Blade.Panels.c1.Global = Blade.Panels.c1.Local + NumPoints;
Blade.Panels.c2.Global = Blade.Panels.c2.Local + NumPoints;
Blade.Panels.c3.Global = Blade.Panels.c3.Local + NumPoints;
Blade.Panels.c4.Global = Blade.Panels.c4.Local + NumPoints;
Blade.Panels.WakeShedders.US.Global = Blade.Panels.WakeShedders.US.Local + NumPanels;
Blade.Panels.WakeShedders.LS.Global = Blade.Panels.WakeShedders.LS.Local + NumPanels;
% NumPanels = NumPanels + numel(Blade.Panels.c1.Local);
% NumPoints = NumPoints + numel(Blade.X);


Blade.C.x = Blade.X(Blade.N.Local);
Blade.C.y = Blade.Y(Blade.N.Local);
Blade.C.z = Blade.Z(Blade.N.Local);


Blade.Faces.C1.Body = [Blade.X(Blade.Panels.c1.Global) Blade.Y(Blade.Panels.c1.Global) Blade.Z(Blade.Panels.c1.Global)];
Blade.Faces.C2.Body = [Blade.X(Blade.Panels.c2.Global) Blade.Y(Blade.Panels.c2.Global) Blade.Z(Blade.Panels.c2.Global)];
Blade.Faces.C3.Body = [Blade.X(Blade.Panels.c3.Global) Blade.Y(Blade.Panels.c3.Global) Blade.Z(Blade.Panels.c3.Global)];
Blade.Faces.C4.Body = [Blade.X(Blade.Panels.c4.Global) Blade.Y(Blade.Panels.c4.Global) Blade.Z(Blade.Panels.c4.Global)];

Blade.Faces.CP.Body = .25*(Blade.Faces.C1.Body + Blade.Faces.C2.Body +...
    Blade.Faces.C3.Body + Blade.Faces.C4.Body);
Blade.Faces.D1.Body = Blade.Faces.C3.Body - Blade.Faces.C1.Body;
Blade.Faces.D2.Body = Blade.Faces.C4.Body - Blade.Faces.C2.Body;


ly =  .5*(Blade.Faces.C1.Body+Blade.Faces.C4.Body) - .5*(Blade.Faces.C2.Body+Blade.Faces.C3.Body);
lymag = sqrt(dot(ly,ly,2));
lz = cross(Blade.Faces.C2.Body-Blade.Faces.C4.Body, Blade.Faces.C3.Body-Blade.Faces.C1.Body);
lzmag = sqrt(dot(lz,lz,2));
Blade.Faces.Area = .5*lzmag;
Blade.Faces.LocalAxis.Y.Body = [ly(:,1)./lymag, ly(:,2)./lymag, ly(:,3)./lymag];
Blade.Faces.LocalAxis.Z.Body = [lz(:,1)./lzmag, lz(:,2)./lzmag, lz(:,3)./lzmag];
Blade.Faces.LocalAxis.X.Body = cross(Blade.Faces.LocalAxis.Y.Body, Blade.Faces.LocalAxis.Z.Body);


