function Blade = GeomProcess(handles)

%%  Get Blade as expected by GeomGen

Blade.Radius = handles.Span.Radius';
Blade.Chord = handles.Span.Chord';
Blade.Theta = handles.Span.Theta';

Blade.MakeCaps = handles.Span.MakeCaps;
Blade.Thickness = handles.Span.Thickness';

Blade.Attitude = [0 0 0];
Blade.Origin = [0 0 0];
Blade.th0 = 0.;
Blade.PitchAxis = str2num(get(handles.pitch_axis,'String'));
Blade.NSpan = length(Blade.Radius);


Blade.Section.Tip.US = handles.Chord.Tip.SectionShape.US;
Blade.Section.Tip.LS = handles.Chord.Tip.SectionShape.LS;
Blade.Section.Tip.X = handles.Chord.Tip.SectionShape.X;

Blade.Section.Root.US = handles.Chord.Root.SectionShape.US;
Blade.Section.Root.LS = handles.Chord.Root.SectionShape.LS;
Blade.Section.Root.X = handles.Chord.Root.SectionShape.X;

Blade.NChord = length(Blade.Section.Root.X);
Blade.Axes = handles.blade_surf_axes;
%%  Pass Blade to GeomGen

Blade = ProcessBlade(Blade);
S = zeros(size(Blade.X(Blade.N.Local)));
for j = 1:size(S,2)
    S(:,j) = j;
end
ax = Blade.Axes;
surf(ax,Blade.X(Blade.N.Local),Blade.Y(Blade.N.Local),Blade.Z(Blade.N.Local),S);

USx = Blade.US.Local.x';
USy = Blade.US.Local.y';
USz = Blade.US.Local.z';
LSx = Blade.LS.Local.x';
LSy = Blade.LS.Local.y';
LSz = Blade.LS.Local.z';


A = [USx(:) USy(:) USz(:)];
B = [LSx(:) LSy(:) LSz(:)];
save US.dat A -ASCII;
save LS.dat B -ASCII;


hold(ax,'on');
surf(ax,Blade.X(Blade.Tip.Inboard.US.N.Local),Blade.Y(Blade.Tip.Inboard.US.N.Local),Blade.Z(Blade.Tip.Inboard.US.N.Local));
surf(ax,Blade.X(Blade.Tip.Inboard.LS.N.Local),Blade.Y(Blade.Tip.Inboard.LS.N.Local),Blade.Z(Blade.Tip.Inboard.LS.N.Local));
surf(ax,Blade.X(Blade.Tip.Outboard.US.N.Local),Blade.Y(Blade.Tip.Outboard.US.N.Local),Blade.Z(Blade.Tip.Outboard.US.N.Local));
surf(ax,Blade.X(Blade.Tip.Outboard.LS.N.Local),Blade.Y(Blade.Tip.Outboard.LS.N.Local),Blade.Z(Blade.Tip.Outboard.LS.N.Local));
view(ax,3);
axis(ax,'equal','tight');
hold(ax,'off');



