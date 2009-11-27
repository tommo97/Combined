function handles = GetCase(handles)
%%  Set up Case panel defaults
handles.Case.output_case = get(handles.output_case,'Value');
handles.Case.output_neu = get(handles.output_neu,'Value');
handles.Case.output_runfile = get(handles.output_runfile,'Value');
handles.Case.commit_to_head = get(handles.commit_to_head,'Value');
handles.Case.output_svn = get(handles.output_svn,'Value');
handles.Case.timestamp = get(handles.timestamp,'Value');