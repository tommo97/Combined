function D = struct_read(fid,format)
%% A function to read a C structure from a binary file.
% struct_read simplifies the reading of a C structure from a binary file by
% providing a functionality similar to the Python 'struct' module. The
% elements of the structure are specified as a series of characters (and
% optional multipliers) in a format string. struct_read then reads these binary
% values returning a cell array of the results. Results are returned either
% as character arrays or doubles since these are the default MATLAB types.
%
% When reading a file, the function checks to see if the "eof" flag has
% been set and if so sets the last element of D to 'eof' and skips
% subsequent reads. Therefore one may check the final value in D to see if
% the end of the file was reached.
%
% ARGUMENTS:
%
%   fid         File handle of an open binary file.
%   format      Format string specifying the C structure. (see below)
%
% FORMAT SPECIFIERS:
% 
% Specifier     Bytes:        Notes
% c             2             char returned as a string character
% C             2             char returned as a number (double)
% b             2             signed char
% B             2             unsigned char
% h             2             short
% H             4             unsigned short
% i             4             int 
% I             4             unsigned int
% l             4             long
% L             4             unsigned long
% q             8             long long
% f             4             float
% d             8             double
% s             char[]        char[] returns character string
%
% MULTIPLIERS:
% Any of the specifiers above may be prepended with an integer indicating
% the number of typed values to be returned as a single entity. For
% example, an 8 character string may be specified by '8c'. 
%
% EXAMPLES:
%
% Given the following C structure definition for a depth measurement from an
% echo-sounder:
%
% typedef struct 
% { 
%     float     depth1; //Depth1 
%     char      units;  // m for meters, f for feet
%     double    time_stamp; //Time stamp, number of seconds elapsed since midnight (00:00:00), January 1, 1970 
%     char[256] navstring; // GPS Position String
% } ECHO; 
%
% fid = fopen('mybinaryfile','r');
% D = struct_read(fid,'fcd256c')
%
% D:
% { [23.456],'m',[1242238192.870432012], ...
%  '$GPGGA,182503.00,7121.11801,N,15651.53926,W,1,12,0.8,20.49,M,-048,M,,*57' }
%
% Val Schmidt
% University of New Hampshire
% Center for Coastal and Ocean Mapping
% 2009

%% Extract the multipliers
tmp = regexp(format,'\D','split');
N = length(tmp);

% While MATLAB is notorious for poorly performing loops, this one executes
% much faster than following commented out single line.
z=1;
mult_idx = [];
for i=1:N
    if ~isempty(tmp{i})
        mult_idx(z) = i;
        z=z+1;
    end
end
%mult_idx = find(cellfun(@(x) ~isempty(x), tmp));

multipliers = ones(1,N);
multipliers(mult_idx) = cellfun(@str2num,tmp(mult_idx));

%% Extract the types
tmp = regexp(format,'\d','split');
types = strcat(tmp{:});

D = cell(1,length(types));
for i = 1:length(types)
    
    % Check for the eof flag. If it's set, return 'eof' as the last element
    % in D and break out of the loop.
    if feof(fid)
        D{end} = 'eof';
        break;
    else
        switch types(i)
            
            case 'c'
                % char as char
                D{i} = fread(fid,multipliers(i),'char=>char');
                
            case 'C'
                % char as number
                D{i} = fread(fid,multipliers(i),'char');
                
            case 'b'
                % signed char
                D{i} = fread(fid,multipliers(i),'schar');
                
            case 'B'
                % unsigned char
                D{i} = fread(fid,multipliers(i),'uchar');
            case 'h'
                % short
                D{i} = fread(fid,multipliers(i),'short');
                
            case 'H'
                % unsigned short
                D{i} = fread(fid,multipliers(i),'ushort');
                
            case 'i'
                % int
                D{i} = fread(fid,multipliers(i),'int');
                
            case 'I'
                % unsigned int
                D{i} = fread(fid,multipliers(i),'uint');
                
            case 'l'
                % long
                D{i} = fread(fid,multipliers(i),'long');
                
            case 'L'
                % unsigned long
                D{i} = fread(fid,multipliers(i),'ulong');
                
            case 'q'
                % long long
                
            case 'Q'
                % unsigned long long
                
            case 'f'
                % float
                D{i} = fread(fid,multipliers(i),'float');
            case 'd'
                % double
                D{i} = fread(fid,multipliers(i),'double');
                
            case 's'
                % char[]  ??? does this require adding 1 to mult for \null?
                D{i} = fread(fid,multipliers(i),'char=>char');
                
        end
    end
end

    

end
