fclose all;clear all;clc;

fid = fopen('bodies.bin', 'r');

NumBodies = fread(fid,1,'int');
Bodies = cell(NumBodies,1);
for i = 1:NumBodies
    Bodies{i}.Name = fread(fid,80,'char=>char')';
    Bodies{i}.NumPanels = fread(fid,1,'int');
    Bodies{i}.Position = fread(fid,3,'double');
    Bodies{i}.Velocity = fread(fid,3,'double');
    Bodies{i}.Attitude = fread(fid,3,'double');
    Bodies{i}.Rates = fread(fid,3,'double');
end
dims = fread(fid,3,'int');
domain = zeros(dims');
x = domain;
y = domain;
z = domain;
for i = 1:dims(3)
    domain(:,:,i) = fread(fid,[dims(1) dims(2)],'double');
end

for i = 1:dims(3)
    x(:,:,i) = fread(fid,[dims(1) dims(2)],'double');
end

for i = 1:dims(3)
    y(:,:,i) = fread(fid,[dims(1) dims(2)],'double');
end

for i = 1:dims(3)
    z(:,:,i) = fread(fid,[dims(1) dims(2)],'double');
end
fclose all;
whos

