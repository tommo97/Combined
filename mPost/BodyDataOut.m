fclose all;clear all;clc; close all
tic
[pathstr, name] = fileparts(mfilename('fullpath'));
cd(pathstr)
%%  Read in binary data
fid = fopen('../bin_files/BodyData50.bin', 'r');
hold all
load('../mat_files/straight.mat')
NumBodies = fread(fid,1,'int');

for i = 1:NumBodies
    %   The following line reads the integer with the length of the name,
    %   then reads that number of character for the name
    Bodies{i}.Name = fread(fid,fread(fid,1,'int'),'char=>char')';
    Bodies{i}.NumPanels = fread(fid,1,'int');
    Bodies{i}.Position = fread(fid,3,'double');
    Bodies{i}.Velocity = fread(fid,3,'double');
    Bodies{i}.Attitude = fread(fid,3,'double');
    Bodies{i}.Rates = fread(fid,3,'double');
    Bodies{i}.Faces.C1.Body = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C2.Body = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C3.Body = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C4.Body = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C1.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C2.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C3.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.C4.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.LocalAxis.X.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.LocalAxis.Y.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.LocalAxis.Z.Global = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.Area1 = fread(fid,Bodies{i}.NumPanels,'double');
    Bodies{i}.Faces.Gamma = fread(fid,Bodies{i}.NumPanels,'double');
    Bodies{i}.Faces.Gamma_Prev = fread(fid,Bodies{i}.NumPanels,'double');
    Bodies{i}.Faces.Mu = fread(fid,Bodies{i}.NumPanels,'double');
    Bodies{i}.Faces.Sigma = fread(fid,Bodies{i}.NumPanels,'double');
    Bodies{i}.Faces.Cpress = fread(fid,Bodies{i}.NumPanels,'double')';
    Bodies{i}.Faces.dF = fread(fid,[3 Bodies{i}.NumPanels],'double')';
    Bodies{i}.Faces.CP.Body = .25*(Bodies{i}.Faces.C1.Body + Bodies{i}.Faces.C2.Body + Bodies{i}.Faces.C3.Body + Bodies{i}.Faces.C4.Body);
end




set(gcf,'Renderer','OpenGL')

%%  Begin Postprocessing
for i = 1:NumBodies
    inds = [1:max(Bodies{i}.Panels.MainPans(:))]';


    Bodies{i}.CPdist = zeros(size(Bodies{1}.Panels.MainPans));
    
    Bodies{i}.Fxdist = zeros(size(Bodies{1}.Panels.MainPans));
    
    Bodies{i}.Fydist = zeros(size(Bodies{1}.Panels.MainPans));
        
    Bodies{i}.Fzdist = zeros(size(Bodies{1}.Panels.MainPans));

    Bodies{i}.GammaDist = zeros(size(Bodies{1}.Panels.MainPans));
    

    Bodies{i}.CPdist(:) = Bodies{i}.Faces.Cpress(inds);
    Bodies{i}.Fxdist(:) = Bodies{i}.Faces.dF(inds,1);
    Bodies{i}.Fydist(:) = Bodies{i}.Faces.dF(inds,2);
    Bodies{i}.Fzdist(:) = Bodies{i}.Faces.dF(inds,3);
    Bodies{i}.GammaDist(:) = Bodies{i}.Faces.Gamma(inds);
    Bodies{i}.Panels.CP.Body.x = zeros(size(Bodies{i}.Panels.MainPans));
    Bodies{i}.Panels.CP.Body.x(:) = Bodies{i}.Faces.CP.Body(inds,1);
    
    Bodies{i}.Panels.CP.Body.y = zeros(size(Bodies{i}.Panels.MainPans));
    Bodies{i}.Panels.CP.Body.y(:) = Bodies{i}.Faces.CP.Body(inds,2);
    
    Bodies{i}.Panels.CP.Body.z = zeros(size(Bodies{i}.Panels.MainPans));
    Bodies{i}.Panels.CP.Body.z(:) = Bodies{i}.Faces.CP.Body(inds,3);
    
    %   Now we can get gradients in a "chordwise" and "spanwise" direction
    
    slice = 12;
    CPmid = Bodies{i}.CPdist(slice,:);
    x = Bodies{i}.Panels.CP.Body.x;
    z = Bodies{i}.Panels.CP.Body.z;
    xmid = x(slice,:);
    cmid = (xmid - min(xmid));
    cmid = cmid/max(abs(cmid));

    fx = Bodies{i}.Fxdist(slice,:);
    fy = Bodies{i}.Fydist(slice,:);
    fz = Bodies{i}.Fzdist(slice,:);
    plot(cmid,-CPmid,'--s')
    %hold all
    %quiver(x(slice,:),z(slice,:),fx,fz)
    %clf
    
    Bodies{i}.GammaDist(slice,:) = 0;
    %surf(Bodies{i}.X(Bodies{i}.N.Local),Bodies{i}.Y(Bodies{i}.N.Local),Bodies{i}.Z(Bodies{i}.N.Local),Bodies{i}.GammaDist);
    %hold all
    %scatter3(Bodies{1}.Panels.CP.Body.x(:),Bodies{1}.Panels.CP.Body.y(:),Bodies{1}.Panels.CP.Body.z(:));
end
axis([-.1 1 -1 4])
axis square


toc