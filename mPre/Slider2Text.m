function handle = Slider2Text(handle)
handle.thickness = get(handle.thickness_slider,'Value');
set(handle.thickness_edit_text,'String', num2str(handle.thickness));