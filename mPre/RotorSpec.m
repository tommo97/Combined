function handles = RotorSpec(handles)

handles.BodyData.numbodies = 1;
handles.BodyData.numblades = str2num(get(handles.nblades,'String'));

if get(handles.blades_as_bodies,'Value')
    numbodies = int16(str2num(get(handles.nblades,'String')));
end

