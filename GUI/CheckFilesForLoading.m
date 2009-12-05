function handles = CheckFilesForLoading(handles)
a = dir('../scratch');
s = size(a,1);
dates = zeros(1,s-2);
for i = 3:s
    dates(-2+i) = a(i).datenum * a(i).isdir;
end

newest = find(dates==max(dates)) + 2;
handles.bin_dir = ['../scratch/' a(newest).name '/output/binary_output/'];

handles = load_listbox(handles.bin_dir, handles);


function handles = load_listbox(dir_path, handles)
%temp = pwd;
%cd(dir_path);
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
%cd(temp);

%set(handles.inputFilesDir,'String',pwd)