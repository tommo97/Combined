function handles = AllDataOut(handles)
hold(handles.goto_axes,'off');
handles = CheckFilesForLoading(handles);
tic
files = get(handles.inputFiles_listbox,'String');
t = zeros(length(files),1);
Fx = t;
Fy = t;
Fz = t;
Mx = t;
My = t;
Mz = t;
for i = 1:length(files)
    binfile = strcat('../bin_files/', files(i));
    fid = fopen(binfile{1}, 'r');
    fread(fid,fread(fid,1,'int'),'char=>char')';
    t(i) = fread(fid,1,'double');
    NumBodies = fread(fid,1,'int');

for j = 1:NumBodies
    %   The following line reads the integer with the length of the name,
    %   then reads that number of character for the name
    Name = fread(fid,fread(fid,1,'int'),'char=>char')';
    NumPanels = fread(fid,1,'int');
    Position = fread(fid,3,'double');
    Velocity = fread(fid,3,'double');
    Attitude = fread(fid,3,'double');
    Rates = fread(fid,3,'double');
    Faces.C1.Body = fread(fid,[3 NumPanels],'double')';
    Faces.C2.Body = fread(fid,[3 NumPanels],'double')';
    Faces.C3.Body = fread(fid,[3 NumPanels],'double')';
    Faces.C4.Body = fread(fid,[3 NumPanels],'double')';
    Faces.C1.Global = fread(fid,[3 NumPanels],'double')';
    Faces.C2.Global = fread(fid,[3 NumPanels],'double')';
    Faces.C3.Global = fread(fid,[3 NumPanels],'double')';
    Faces.C4.Global = fread(fid,[3 NumPanels],'double')';
    Faces.LocalAxis.X.Global = fread(fid,[3 NumPanels],'double')';
    Faces.LocalAxis.Y.Global = fread(fid,[3 NumPanels],'double')';
    Faces.LocalAxis.Z.Global = fread(fid,[3 NumPanels],'double')';
    Faces.Area1 = fread(fid,NumPanels,'double');
    Faces.Gamma = fread(fid,NumPanels,'double');
    Faces.Gamma_Prev = fread(fid,NumPanels,'double');
    Faces.Mu = fread(fid,NumPanels,'double');
    Faces.Sigma = fread(fid,NumPanels,'double');
    Faces.Cpress = fread(fid,NumPanels,'double')';
    Faces.dF = fread(fid,[3 NumPanels],'double')';
    Faces.CP.Body = .25*(Faces.C1.Body + Faces.C2.Body + Faces.C3.Body + Faces.C4.Body);
    Fx(i,j) = sum(Faces.dF(:,1));
    Fy(i,j) = sum(Faces.dF(:,2));
    Fz(i,j) = sum(Faces.dF(:,3));
    M = cross(Faces.CP.Body,Faces.dF,2);
    Mx(i,j) = sum(M(:,1));
    My(i,j) = sum(M(:,2));
    Mz(i,j) = sum(M(:,3));
end
end
cla(handles.goto_axes,'reset');
hold(handles.goto_axes,'all');
for i = 1:size(Fx,2);
plot(handles.goto_axes,t,Fx(:,i),'r');
plot(handles.goto_axes,t,Fy(:,i),'g');
plot(handles.goto_axes,t,Fz(:,i),'b');
plot(handles.goto_axes,t,Mx(:,i),'k');
plot(handles.goto_axes,t,My(:,i),'y');
plot(handles.goto_axes,t,Mz(:,i),'p');
legend(handles.goto_axes,'Fx','Fy','Fz','Mx','My','Mz');
end

%plot(handles.goto_axes,t,Fy);
%plot(handles.goto_axes,t,Fz/1000);
toc