function handles = CheckFilesForLoading(handles)

handles = load_listbox('../bin_files', handles);


function handles = load_listbox(dir_path, handles)
temp = pwd;
cd(dir_path);
dir_struct = dir([dir_path '/*.bin']);
if length({dir_struct.name}) > 0
names = {dir_struct.name};
for i = 1:length({dir_struct.name})
    num(i) = str2num(names{i}(9:end-4));
end 

[sorted_nums,sorted_index] = sort(num);
handles.file_names = [];
for i = 1:length(num)
    handles.file_names{i} = names{sorted_index(i)};
end

handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
set(handles.inputFiles_listbox,'String',handles.file_names,...
	'Value',get(handles.inputFiles_listbox,'Value'));
else
    set(handles.inputFiles_listbox,'Value',1);
    set(handles.inputFiles_listbox,'String','Empty',...
	'Value',get(handles.inputFiles_listbox,'Value'));
end
cd(temp);

%set(handles.inputFilesDir,'String',pwd)