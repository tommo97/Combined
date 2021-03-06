function Blade = GeomProcess(handles)

%%  Get Blade as expected by GeomGen

Blade.Radius = handles.Span.Radius';
Blade.Chord = handles.Span.Chord';
Blade.Theta = handles.Span.Theta';
Blade.Skew = handles.Span.Skew';
Blade.RoundTips = handles.Span.RoundTips;
Blade.Rake = handles.Span.Rake';
Blade.isProp = handles.Span.isProp;
Blade.Thickness = handles.Span.Thickness';
Blade.Camber = handles.Span.Camber';

if Blade.RoundTips
    Blade.num_tip_pans = handles.Span.num_tip_pans;
    %disp(Blade.num_tip_pans);
end

Blade.Reverse = handles.Span.REVERSE;

Blade.isNREL = handles.Span.isNREL;
Blade.isSOTON = handles.Span.isSOTON;
Blade.isBarltrop = handles.Span.isBarltrop;

if Blade.isNREL || Blade.isSOTON || Blade.isBarltrop
    Blade.n1 = handles.Span.n1;
    Blade.n2 = handles.Span.n2;
end

Blade.TransitionPiece = handles.Span.TransitionPiece;



Blade.Attitude = [0 0 0];
Blade.Origin = [0 0 0];
Blade.th0 = 0.;
Blade.PitchAxis = str2num(get(handles.pitch_axis,'String'));
Blade.NSpan = length(Blade.Radius);


Blade.Section.Tip.US = handles.Chord.Tip.SectionShape.US;
Blade.Section.Tip.LS = handles.Chord.Tip.SectionShape.LS;
Blade.Section.Tip.X = handles.Chord.Tip.SectionShape.X;
Blade.Section.Tip.Camber = handles.Chord.Tip.SectionShape.Camber;


Blade.Section.Root.US = handles.Chord.Root.SectionShape.US;
Blade.Section.Root.LS = handles.Chord.Root.SectionShape.LS;
Blade.Section.Root.X = handles.Chord.Root.SectionShape.X;
Blade.Section.Root.Camber = handles.Chord.Root.SectionShape.Camber;

Blade.NChord = length(Blade.Section.Root.X);
Blade.Axes = handles.blade_surf_axes;
%%  Pass Blade to GeomGen

Blade = ProcessBlade(Blade);
S = zeros(size(Blade.X(Blade.N.Local)));
for j = 1:size(S,2)
    S(:,j) = j;
end
ax = Blade.Axes;
surf(ax,Blade.X(Blade.N.Local),Blade.Y(Blade.N.Local),Blade.Z(Blade.N.Local),S,'facecolor',[0.75 0.75 1],'edgecolor','w');
hold(ax,'on');
surf(ax,Blade.X(Blade.Tip.Inboard.US.N.Local),Blade.Y(Blade.Tip.Inboard.US.N.Local),Blade.Z(Blade.Tip.Inboard.US.N.Local),'facecolor',[0.75 0.75 1],'edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Inboard.LS.N.Local),Blade.Y(Blade.Tip.Inboard.LS.N.Local),Blade.Z(Blade.Tip.Inboard.LS.N.Local),'facecolor',[0.75 0.75 1],'edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Outboard.US.N.Local),Blade.Y(Blade.Tip.Outboard.US.N.Local),Blade.Z(Blade.Tip.Outboard.US.N.Local),'facecolor',[0.75 0.75 1],'edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Outboard.LS.N.Local),Blade.Y(Blade.Tip.Outboard.LS.N.Local),Blade.Z(Blade.Tip.Outboard.LS.N.Local),'facecolor',[0.75 0.75 1],'edgecolor','b');

surf(ax,Blade.X(Blade.N.Local),Blade.Y(Blade.N.Local),Blade.Z(Blade.N.Local),S,'facecolor','none','edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Inboard.US.N.Local),Blade.Y(Blade.Tip.Inboard.US.N.Local),Blade.Z(Blade.Tip.Inboard.US.N.Local),'facecolor','none','edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Inboard.LS.N.Local),Blade.Y(Blade.Tip.Inboard.LS.N.Local),Blade.Z(Blade.Tip.Inboard.LS.N.Local),'facecolor','none','edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Outboard.US.N.Local),Blade.Y(Blade.Tip.Outboard.US.N.Local),Blade.Z(Blade.Tip.Outboard.US.N.Local),'facecolor','none','edgecolor','b');
surf(ax,Blade.X(Blade.Tip.Outboard.LS.N.Local),Blade.Y(Blade.Tip.Outboard.LS.N.Local),Blade.Z(Blade.Tip.Outboard.LS.N.Local),'facecolor','none','edgecolor','b');
%shading interp
%surf(ax,Blade.X(Blade.N.Local),Blade.Y(Blade.N.Local),Blade.Z(Blade.N.Local),S,'facecolor','none');
%surf(ax,Blade.X(Blade.Tip.Inboard.US.N.Local),Blade.Y(Blade.Tip.Inboard.US.N.Local),Blade.Z(Blade.Tip.Inboard.US.N.Local),'facecolor','none');
%surf(ax,Blade.X(Blade.Tip.Inboard.LS.N.Local),Blade.Y(Blade.Tip.Inboard.LS.N.Local),Blade.Z(Blade.Tip.Inboard.LS.N.Local),'facecolor','none');
%surf(ax,Blade.X(Blade.Tip.Outboard.US.N.Local),Blade.Y(Blade.Tip.Outboard.US.N.Local),Blade.Z(Blade.Tip.Outboard.US.N.Local),'facecolor','none');
%surf(ax,Blade.X(Blade.Tip.Outboard.LS.N.Local),Blade.Y(Blade.Tip.Outboard.LS.N.Local),Blade.Z(Blade.Tip.Outboard.LS.N.Local),'facecolor','none');
view(ax,3);
axis(ax,'equal','tight');
hold(ax,'off');
lighting phong
camlight(-90,0)
camlight(90,0)
