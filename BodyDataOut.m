fclose all;clear all;clc;


%%  Read in binary data
fid = fopen('BodyData.bin', 'r');
hold all

NumBodies = fread(fid,1,'int');
Bodies = cell(NumBodies,1);
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
end




set(gcf,'Renderer','OpenGL')

%%  Begin Postprocessing
for i = 1:NumBodies
Bodies{i}.Faces.CP.Body = (Bodies{i}.Faces.C1.Body + Bodies{i}.Faces.C2.Body +...
    Bodies{i}.Faces.C3.Body + Bodies{i}.Faces.C4.Body);
Bodies{i}.Faces.D1.Body = Bodies{i}.Faces.C3.Body - Bodies{i}.Faces.C1.Body;
Bodies{i}.Faces.D2.Body = Bodies{i}.Faces.C4.Body - Bodies{i}.Faces.C2.Body;


lx =  .5*(Bodies{i}.Faces.C1.Body+Bodies{i}.Faces.C4.Body) - .5*(Bodies{i}.Faces.C2.Body+Bodies{i}.Faces.C3.Body);
lxmag = sqrt(dot(lx,lx,2));
lz = cross(Bodies{i}.Faces.C4.Body-Bodies{i}.Faces.C2.Body, Bodies{i}.Faces.C3.Body-Bodies{i}.Faces.C1.Body);
lzmag = sqrt(dot(lz,lz,2));
Bodies{i}.Faces.Area = .5*lzmag;
Bodies{i}.Faces.LocalAxis.X.Body = [lx(:,1)./lxmag, lx(:,2)./lxmag, lx(:,3)./lxmag];
Bodies{i}.Faces.LocalAxis.Z.Body = [lz(:,1)./lzmag, lz(:,2)./lzmag, lz(:,3)./lzmag];
Bodies{i}.Faces.LocalAxis.Y.Body = cross(Bodies{i}.Faces.LocalAxis.X.Body, Bodies{i}.Faces.LocalAxis.Z.Body);
scatter3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3));
quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
    Bodies{i}.Faces.LocalAxis.X.Body(:,1),Bodies{i}.Faces.LocalAxis.X.Body(:,2),Bodies{i}.Faces.LocalAxis.X.Body(:,3));
quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
    Bodies{i}.Faces.LocalAxis.Y.Body(:,1),Bodies{i}.Faces.LocalAxis.Y.Body(:,2),Bodies{i}.Faces.LocalAxis.Y.Body(:,3));
quiver3(Bodies{i}.Faces.CP.Body(:,1),Bodies{i}.Faces.CP.Body(:,2),Bodies{i}.Faces.CP.Body(:,3),...
    Bodies{i}.Faces.LocalAxis.Z.Body(:,1),Bodies{i}.Faces.LocalAxis.Z.Body(:,2),Bodies{i}.Faces.LocalAxis.Z.Body(:,3));
end


axis equal


