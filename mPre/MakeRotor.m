function Rotor = MakeRotor(handles)
Rotor.NumBlades = str2num(get(handles.nblades,'String'));
Rotor.Axes = handles.full_turbine_axes;
Rotor.Attitude = [  str2num(get(handles.psi,'String'))...
                    str2num(get(handles.theta,'String'))...
                    str2num(get(handles.phi,'String'))];
Rotor.Origin = [str2num(get(handles.Ox,'String'))...
                str2num(get(handles.Oy,'String'))...
                str2num(get(handles.Oz,'String'))];
Rotor.th0 = deg2rad(str2num(get(handles.theta_hub,'String')));
ax = Rotor.Axes;
hold(ax,'off');
cla(ax);
azimuth = ([1:Rotor.NumBlades] - 1)*2*pi/Rotor.NumBlades;
for i = 1:Rotor.NumBlades
    Rotor.Blade{i} = handles.Blade;
    Rotor.Blade{i}.Azimuth = azimuth(i);
    x0 = Rotor.Blade{i}.X;
    y0 = Rotor.Blade{i}.Y;
    z0 = Rotor.Blade{i}.Z;  
    % Now put this rotor blade into appropriate azimuthal position -
    % equivilent to rotation about body roll axis
    TRANS2 = MakeEulerMatrix(Rotor.Attitude);
    TRANS = MakeEulerMatrix([azimuth(i) Rotor.th0 0]);
    x1 = x0*TRANS(1,1) + y0*TRANS(1,2) + z0*TRANS(1,3);
    y1 = x0*TRANS(2,1) + y0*TRANS(2,2) + z0*TRANS(2,3);
    z1 = x0*TRANS(3,1) + y0*TRANS(3,2) + z0*TRANS(3,3);
    
    Rotor.Blade{i}.X = Rotor.Origin(1) + x1*TRANS2(1,1) + y1*TRANS2(1,2) + z1*TRANS2(1,3);
    Rotor.Blade{i}.Y = Rotor.Origin(2) + x1*TRANS2(2,1) + y1*TRANS2(2,2) + z1*TRANS2(2,3);
    Rotor.Blade{i}.Z = Rotor.Origin(3) + x1*TRANS2(3,1) + y1*TRANS2(3,2) + z1*TRANS2(3,3);
    hold(ax,'on');
    surf(ax,Rotor.Blade{i}.X(Rotor.Blade{i}.N.Local),Rotor.Blade{i}.Y(Rotor.Blade{i}.N.Local),Rotor.Blade{i}.Z(Rotor.Blade{i}.N.Local));

    surf(ax,Rotor.Blade{i}.X(Rotor.Blade{i}.Tip.Inboard.US.N.Local),Rotor.Blade{i}.Y(Rotor.Blade{i}.Tip.Inboard.US.N.Local),Rotor.Blade{i}.Z(Rotor.Blade{i}.Tip.Inboard.US.N.Local));
    surf(ax,Rotor.Blade{i}.X(Rotor.Blade{i}.Tip.Inboard.LS.N.Local),Rotor.Blade{i}.Y(Rotor.Blade{i}.Tip.Inboard.LS.N.Local),Rotor.Blade{i}.Z(Rotor.Blade{i}.Tip.Inboard.LS.N.Local));
    surf(ax,Rotor.Blade{i}.X(Rotor.Blade{i}.Tip.Outboard.US.N.Local),Rotor.Blade{i}.Y(Rotor.Blade{i}.Tip.Outboard.US.N.Local),Rotor.Blade{i}.Z(Rotor.Blade{i}.Tip.Outboard.US.N.Local));
    surf(ax,Rotor.Blade{i}.X(Rotor.Blade{i}.Tip.Outboard.LS.N.Local),Rotor.Blade{i}.Y(Rotor.Blade{i}.Tip.Outboard.LS.N.Local),Rotor.Blade{i}.Z(Rotor.Blade{i}.Tip.Outboard.LS.N.Local));
    view(ax,3);
end
axis(ax,'equal','tight');
Rotor.Split = get(handles.blades_as_bodies,'Value') ;

hold(ax,'off');