%function A = ReadBin(fname)
%%  Read in binary data
clear all mex;
fname = 'binfile.bin';
fclose all
fid = fopen(fname, 'r');
fseek(fid, -64, 'eof');
A = fread(fid,64,'*char')';
n = strfind(A,'CHAMELEON');


keywordPos = n - 64 - 1;

fseek(fid, keywordPos, 'eof');
B = fread(fid,18,'*char')'


fseek(fid, keywordPos - 5, 'eof');
footersize = fread(fid,1,'int')
