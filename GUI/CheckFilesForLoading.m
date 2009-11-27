function handles = CheckFilesForLoading(handles)

handles = load_listbox('../bin_files', handles);


function handles = load_listbox(dir_path, handles)
temp = pwd;
cd(dir_path);
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
set(handles.inputFiles_listbox,'String',handles.file_names,...
	'Value',get(handles.inputFiles_listbox,'Value'));
cd(temp);
%set(handles.inputFilesDir,'String',pwd)