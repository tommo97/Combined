function handles = BodyDataOut(handles)
tic
%[pathstr, name] = fileparts(mfilename('fullpath'));
%cd(pathstr)
%%  Read in binary data
files = get(handles.inputFiles_listbox,'String');

binfile = strcat('../bin_files/', files(get(handles.inputFiles_listbox,'Value')));

fid = fopen(binfile{1}, 'r');

handles.PostComp.CaseName = fread(fid,fread(fid,1,'int'),'char=>char')';

set(handles.casename,'String',handles.PostComp.CaseName);

load(['../mat_files/' handles.PostComp.CaseName '.mat']);

handles.PostComp.Bodies = Bodies;

NumBodies = fread(fid,1,'int');

for i = 1:NumBodies
    %   The following line reads the integer with the length of the name,
    %   then reads that number of character for the name
    handles.PostComp.Bodies{i}.Name = fread(fid,fread(fid,1,'int'),'char=>char')';
    handles.PostComp.Bodies{i}.NumPanels = fread(fid,1,'int');
    handles.PostComp.Bodies{i}.Position = fread(fid,3,'double');
    handles.PostComp.Bodies{i}.Velocity = fread(fid,3,'double');
    handles.PostComp.Bodies{i}.Attitude = fread(fid,3,'double');
    handles.PostComp.Bodies{i}.Rates = fread(fid,3,'double');
    handles.PostComp.Bodies{i}.Faces.C1.Body = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C2.Body = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C3.Body = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C4.Body = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C1.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C2.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C3.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.C4.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.LocalAxis.X.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.LocalAxis.Y.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.LocalAxis.Z.Global = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.Area1 = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double');
    handles.PostComp.Bodies{i}.Faces.Gamma = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double');
    handles.PostComp.Bodies{i}.Faces.Gamma_Prev = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double');
    handles.PostComp.Bodies{i}.Faces.Mu = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double');
    handles.PostComp.Bodies{i}.Faces.Sigma = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double');
    handles.PostComp.Bodies{i}.Faces.Cpress = fread(fid,handles.PostComp.Bodies{i}.NumPanels,'double')';
    handles.PostComp.Bodies{i}.Faces.dF = fread(fid,[3 handles.PostComp.Bodies{i}.NumPanels],'double')';
    handles.PostComp.Bodies{i}.Faces.CP.Body = .25*(handles.PostComp.Bodies{i}.Faces.C1.Body + handles.PostComp.Bodies{i}.Faces.C2.Body + handles.PostComp.Bodies{i}.Faces.C3.Body + handles.PostComp.Bodies{i}.Faces.C4.Body);
end


set(handles.span_pan,'String',num2str(size(handles.PostComp.Bodies{i}.Panels.MainPans,1)));
toc