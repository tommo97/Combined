function handle = Text2Slider(handle)
sliderValue = get(handle.thickness_edit_text,'String');
if (isempty(str2num(sliderValue)) || str2num(sliderValue) < 0)
    set(handle.thickness_slider,'Value',0);
    set(handle.thickness_edit_text,'String','0');
    handle.thickness = 0;
elseif (str2num(sliderValue) > 1)
    set(handle.thickness_slider,'Value',1);
    set(handle.thickness_edit_text,'String','1');
    handle.thickness = 1;
else
    set(handle.thickness_slider,'Value',str2num(sliderValue));
    handle.thickness = str2num(sliderValue);
end


        