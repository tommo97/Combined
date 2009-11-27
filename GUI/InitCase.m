function handles = InitCase(handles)
%%  Set up Case panel defaults
set(handles.output_case,'Value',true);
set(handles.output_neu,'Value',true);
set(handles.output_runfile,'Value',true);
set(handles.commit_to_head,'enable','off');
set(handles.output_svn,'enable','off');
set(handles.timestamp,'Value',true);