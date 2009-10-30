function CaseOut = AddVal2Case(CaseData,Name,Val,Type)
if isfield(CaseData,'NumVals')
    CaseOut = CaseData;
    n = CaseData.NumVals + 1;
else
    n = 1;
    CaseOut.NumVals = 0;
    CaseOut.Data = struct([]);
end

[x y] = size(Val);

CaseOut.Data{n}.Name = Name;
str = [];
Format = 's';
if (strcmp(Type, 'int') || strcmp(Type, 'integer'))
    Format = '%g';
end

if (strcmp(Type, 'string'))
    Format = '%s';
end
if (strcmp(Type, 'real') || strcmp(Type, 'double') || strcmp(Type, 'float'))
    Format = '%-1.12e';
end
str = [];
if iscell(Val)
    for i = 1:size(Val,2)
        str = [str '[' sprintf(Format,Val{i}) ']'];
    end
else
    str = [str sprintf(Format,Val)];
end

disp(str);


CaseOut.Data{n}.Val = str;
CaseOut.Data{n}.Type = Format;
CaseOut.NumVals = CaseOut.NumVals + 1;
