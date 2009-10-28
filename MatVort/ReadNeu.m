function [Xdata PNLS] = ReadNeu(neufile)

fid = fopen(neufile);


fgetl(fid); %  Read control info header
fgetl(fid); %  Read filetype specification
fname = fgetl(fid); %  Read filename
disp(fname);
fgetl(fid); %  Read generating program
fgetl(fid); %  Read date and time of generation
fgetl(fid); %  Read mesh data headers

data = textscan(fgetl(fid),'%d');   %  Read mesh data

NUMNP = data{1}(1);
NELEM = data{1}(2);
NGRPS = data{1}(3);
NBSETS = data{1}(4);
NDFCD = data{1}(5);
NDFVL = data{1}(6);

Xdata = zeros(NUMNP,3);
PNLS = zeros(NELEM,4);

tline = fgetl(fid); %  Read ENDOFSECTION
tline = fgetl(fid); %  Read NODAL COORDINATES 2.4.6

fod = fopen('temp.tmp','wt');
EOS = false;
while ischar(tline) && ~EOS
    tline = fgetl(fid);
    EOS = sum(tline(1:12) =='ENDOFSECTION');
    if ~EOS
        tmp = textscan(tline(1:12),'%d');
        fprintf(fod, '%s\n', tline);
        i = tmp{1}(1);
    end
end
fclose(fod);

data = dlmread('temp.tmp');

if size(data,2) == 3
    data = [data zeros(size(data,1),1)];
end
Xdata(data(:,1),:) = data(:,2:4);

fgetl(fid); %  Read ELEMENTS/CELL S2.4.6

fod = fopen('temp.tmp','wt');
EOS = false;
while ischar(tline) && ~EOS
    tline = fgetl(fid);
    EOS = sum(tline(1:12) =='ENDOFSECTION');
    if ~EOS
        tmp = textscan(tline(1:12),'%d');
        fprintf(fod, '%s\n', tline);
        i = tmp{1}(1);
    end
end

fclose(fod);
data = dlmread('temp.tmp');

PNLS(data(:,1),:) = data(:,4:7);

fclose(fid);