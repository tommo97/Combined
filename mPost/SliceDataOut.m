clear all; clear mex; clc;
files = dir('R*.mat');

close all
s = size(files,1);

scale = 16;

fname = files(s).name;
Data = load(fname);


names{1} = '0pt1';
names{2} = '0pt2';
names{3} = '0pt5';
names{4} = '1pt0';
names{5} = '2pt0';
names{6} = '4pt0';
names{7} = '6pt0';

for i = 1:length(names)
    eval(['Vel{i}.Wake.x = Data.SliceWakeVel' names{i} '_x;']);
    eval(['Vel{i}.Wake.y = Data.SliceWakeVel' names{i} '_y;']);
    eval(['Vel{i}.Wake.z = Data.SliceWakeVel' names{i} '_z;']);
    
    eval(['Vel{i}.Body.x = Data.SliceBodyVel' names{i} '_x;']);
    eval(['Vel{i}.Body.y = Data.SliceBodyVel' names{i} '_y;']);
    eval(['Vel{i}.Body.z = Data.SliceBodyVel' names{i} '_z;']);
    
    eval(['Vel{i}.Vorton.x = Data.SliceVortonWakeVel' names{i} '_x;']);
    eval(['Vel{i}.Vorton.y = Data.SliceVortonWakeVel' names{i} '_y;']);
    eval(['Vel{i}.Vorton.z = Data.SliceVortonWakeVel' names{i} '_z;']);
    
    eval(['Vel{i}.ProtoWake.x = Data.SliceProtoWakeVel' names{i} '_x;']);
    eval(['Vel{i}.ProtoWake.y = Data.SliceProtoWakeVel' names{i} '_y;']);
    eval(['Vel{i}.ProtoWake.z = Data.SliceProtoWakeVel' names{i} '_z;']);
    
    eval(['Phi{i}.ProtoWake = Data.SliceProtoWakePhi' names{i} ';']);
    eval(['Phi{i}.Body = Data.SliceBodyPhi' names{i} ';']);
    
    %[Vel{i}.Body.y Vel{i}.Body.z] = gradient(Phi{i}.Body);
    Vel{i}.x = -Vel{i}.ProtoWake.x + Vel{i}.Wake.x - Vel{i}.Vorton.x - Vel{i}.Body.x;
    Vel{i}.y = -Vel{i}.ProtoWake.y + Vel{i}.Wake.y - Vel{i}.Vorton.y - Vel{i}.Body.y;
    Vel{i}.z = -Vel{i}.ProtoWake.z + Vel{i}.Wake.z - Vel{i}.Vorton.z - Vel{i}.Body.z;
    
    eval(['Pos{i}.x = Data.SlicePosn' names{i} '_x;']);
    eval(['Pos{i}.y = Data.SlicePosn' names{i} '_y;']);
    eval(['Pos{i}.z = Data.SlicePosn' names{i} '_z;']);
    [Vel{i}.CurlX,Vel{i}.CavX]= curl(Pos{i}.y',Pos{i}.z',Vel{i}.y', Vel{i}.z');
    
    hold all
    plot(Pos{1}.y(1:80,80)/scale + 3.3,Vel{i}.z(1:80,80)/4)
    axis([-2 4 -.9 .9])
    grid on
    box on
    leg{i}=num2str(Data.Time + Pos{i}.x(1)/scale);

end
legend(leg)

