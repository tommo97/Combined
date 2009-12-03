function handles = AllDataOut(handles,keep)
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
    Vinf = fread(fid,[3 1],'double');
    rho = fread(fid,1,'double');
    NumBodies = fread(fid,1,'int');

for j = 1:NumBodies
    MaxRad = max(handles.PostComp.Bodies{j}.Radius);
    
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
    Faces.CP.Global = .25*(Faces.C1.Global + Faces.C2.Global + Faces.C3.Global + Faces.C4.Global);
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
    M = cross(Faces.CP.Global,Faces.dF,2);
    Mx(i,j) = sum(M(:,1));
    My(i,j) = sum(M(:,2));
    Mz(i,j) = sum(M(:,3));
    ax = handles.PostComp.Bodies{j}.RotorAxis;
    Pwr(i,j) = sum(handles.PostComp.Bodies{j}.Rates'.*sum(M,1).*ax);
    V = ax*(Vinf - handles.PostComp.Bodies{j}.Velocity);
    Cp(i,j) = Pwr(i,j)/(.5*V*V*V*rho*pi*MaxRad*MaxRad);
    Ct(i,j) = Fz(i,j)/(0.5*rho*V*V*pi*MaxRad*MaxRad);
end
fclose(fid);
end
cla(handles.coefft_axes,'reset');
cla(handles.goto_axes,'reset');
hold(handles.goto_axes,'all');
hold(handles.coefft_axes,'all');

%   Some basic processing

CtMu3 = mean(sum(Ct,2));
CpMu3 = mean(sum(Cp,2));
CtSigma3 = std(sum(Ct,2));
CpSigma3 = std(sum(Cp,2));

outliersCt = abs(sum(Ct,2) - CtMu3) > 5*CtSigma3;
outliersCp = abs(sum(Cp,2) - CpMu3) > 5*CpSigma3;
Ctm = Ct;
Cpm = Cp;
Ctm(outliersCt,:) = NaN;
Cpm(outliersCp,:) = NaN;

plot(handles.coefft_axes,t,sum(Cpm,2),'.b');
plot(handles.coefft_axes,t,sum(Ctm,2),'.r');
for i = 1:size(Fx,2);

% plot(handles.goto_axes,t,Fx(:,i),'r');
% plot(handles.goto_axes,t,Fy(:,i),'g');

if i==1
plot(handles.goto_axes,t,Fz(:,i),'-.b');
plot(handles.goto_axes,t,Mz(:,i),'-.r');
plot(handles.goto_axes,t,Pwr(:,i),'-.g');

else
    plot(handles.goto_axes,t,Fz(:,i),'-ob');
    plot(handles.goto_axes,t,Mz(:,i),'-or');
    plot(handles.goto_axes,t,Pwr(:,i),'-og');
end
% plot(handles.goto_axes,t,Mx(:,i),'k');
% plot(handles.goto_axes,t,My(:,i),'y');
%plot(handles.goto_axes,t,Mz(:,i),'p');
legend(handles.goto_axes,'Fz','Mz','Location','SouthWest');


end

%plot(handles.goto_axes,t,Fy);
%plot(handles.goto_axes,t,Fz/1000);
drawnow


toc

