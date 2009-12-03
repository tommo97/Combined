function handles = Kinematics(handles)

for i = 1:handles.BodyData.numbodies
    handles.BodyData.attitudes{i}   = [str2num(get(handles.psi,'String')) str2num(get(handles.theta,'String')) str2num(get(handles.phi,'String'))];
    handles.BodyData.origin{i}      = [str2num(get(handles.Ox,'String')) str2num(get(handles.Oy,'String')) str2num(get(handles.Oz,'String'))];
    handles.BodyData.rates{i}       = [str2num(get(handles.p,'String')) str2num(get(handles.q,'String')) str2num(get(handles.r,'String'))];
    handles.BodyData.vels{i}        = [str2num(get(handles.u,'String')) str2num(get(handles.v,'String')) str2num(get(handles.w,'String'))];
    axis = [str2num(get(handles.rotor_axis_x,'String')) str2num(get(handles.rotor_axis_y,'String')) str2num(get(handles.rotor_axis_z,'String'))];
    axis = axis./(sqrt(dot(axis,axis)));
    handles.BodyData.RotorAxis{i} = axis;
    rotorspeed = 2*pi*str2num(get(handles.rpm,'String'))/60;
    handles.BodyData.rates{i} = handles.BodyData.rates{i} + rotorspeed*axis;
    handles.BodyData.names{i} = [get(handles.name,'String') '_' num2str(i)];
end