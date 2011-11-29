clear all
close all
clc

fid = fopen('numbers');

%%  First find length of file

isdone = false;
position_old = -1e32;
while ~isdone
    fseek(fid,1,'cof');
    
    position = ftell(fid);
    if (position==position_old)
        isdone = true;
    end
    position_old = position;
end

%%  Rewind to start and read data
frewind(fid);
data_in = fread(fid, position, 'double');


%%  Close file
fclose all;

data_out = rand(12345,3);

fid = fopen('numbers2','w');
fwrite(fid, data_out,'double');
fclose(fid);