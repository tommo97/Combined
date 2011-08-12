function handles = MakeNEU(handles)
Bodies = handles.Rotor.Blade;
name = handles.fullname;
split = handles.Rotor.Split;
dir1 = '../neu_files/';
fname = [dir1 name '.neu'];
fid = fopen(fname, 'wt');
Empty = 0;
numBodies = length(Bodies);

nPts = 0;
nPnls = 0;

if split
    nGrps = length(Bodies);
    nBCs = nGrps;
else
    nGrps = 1;
    nBCs = 1;
end
PtsX = [];
nC = [];
mC = [];
C1 = [];
C2 = [];
C3 = [];
C4 = [];
W = [];
F = [];
f{1} = 2;
f{2} = 4;
NumPanels = 0;
NumPoints = 0;
for i = 1:numBodies
    Bodies{i}.Panels.Start = nPnls + 1;
    nPts = nPts + Bodies{i}.nPts;
    nPnls = nPnls + Bodies{i}.nPnls;
    Bodies{i}.Panels.Finish = nPnls;
    disp([Bodies{i}.Panels.Start Bodies{i}.Panels.Finish]);
    Bodies{i}.N.Global = Bodies{i}.N.Local + NumPoints;
    
    
    Bodies{i}.Tip.Inboard.US.N.Global = Bodies{i}.Tip.Inboard.US.N.Local + NumPoints;
    Bodies{i}.Tip.Inboard.LS.N.Global = Bodies{i}.Tip.Inboard.LS.N.Local + NumPoints;
    Bodies{i}.Tip.Outboard.US.N.Global = Bodies{i}.Tip.Outboard.US.N.Local + NumPoints;
    Bodies{i}.Tip.Outboard.LS.N.Global = Bodies{i}.Tip.Outboard.LS.N.Local + NumPoints;
    
    
    
    Bodies{i}.Panels.c1.Global = Bodies{i}.Panels.c1.Local + NumPoints;
    Bodies{i}.Panels.c2.Global = Bodies{i}.Panels.c2.Local + NumPoints;
    Bodies{i}.Panels.c3.Global = Bodies{i}.Panels.c3.Local + NumPoints;
    Bodies{i}.Panels.c4.Global = Bodies{i}.Panels.c4.Local + NumPoints;
    Bodies{i}.Panels.WakeShedders.US.Global = Bodies{i}.Panels.WakeShedders.US.Local + NumPanels;
    Bodies{i}.Panels.WakeShedders.LS.Global = Bodies{i}.Panels.WakeShedders.LS.Local + NumPanels;
    NumPanels = NumPanels + numel(Bodies{i}.Panels.c1.Local);
    NumPoints = NumPoints + numel(Bodies{i}.X);
    PtsX = [PtsX;[Bodies{i}.X Bodies{i}.Y Bodies{i}.Z]];
    nC = [nC;Bodies{i}.n];
    mC = [mC;Bodies{i}.m];
    C1 = [C1;Bodies{i}.Panels.c1.Global(:)];
    C2 = [C2;Bodies{i}.Panels.c2.Global(:)];
    C3 = [C3;Bodies{i}.Panels.c3.Global(:)];
    C4 = [C4;Bodies{i}.Panels.c4.Global(:)];
    f0 = 2*ones(size(Bodies{i}.Panels.WakeShedders.US.Global(:)));
    f1 = f{1}*ones(size(Bodies{i}.Panels.WakeShedders.US.Global(:)));
    f2 = f{2}*ones(size(Bodies{i}.Panels.WakeShedders.LS.Global(:)));
    W = [W;[Bodies{i}.Panels.WakeShedders.US.Global(:)  f0 f1]];
    W = [W;[Bodies{i}.Panels.WakeShedders.LS.Global(:)  f0 f2]];
end



handles.neumsg = ['Writing ' num2str(nGrps) ' groups with ' num2str(nPts) ' points and ' num2str(nPnls) ' panels.'];


fprintf(fid,'        CONTROL INFO 2.4.6\n** GAMBIT NEUTRAL FILE\n%s\n',name);
fprintf(fid,'PROGRAM:                MATLAB     VERSION:  %s\n', version);
fprintf(fid,'%s\n', datestr(now));
fprintf(fid,'\tNUMNP\tNELEM\tNGRPS\tNBSETS\tNDFCD\tNDFVL\n');
fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\n',nPts,nPnls,nGrps,nBCs,Empty,Empty);
fprintf(fid,'ENDOFSECTION\n   NODAL COORDINATES 2.4.6\n');

str = [repmat(' ',nPts,10 - size(num2str(nPts),2)) num2str([1:nPts]','%g    ') repmat('   ',nPts,1) num2str(PtsX,'%-1.12e    ')];


for i = 1:nPts
    fprintf(fid,'%s\n',str(i,:));
    %    fprintf(fid,'\t%g\t%-1.12e\t%-1.12e\t%-1.12e\n',i,M1(i,1),M1(i,2),M1(i,3));
end
fprintf(fid,'ENDOFSECTION\n      ELEMENTS/CELLS 2.4.6\n');
for i = 1:size(C1,1) - 1
    fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',i,2,4,C1(i),C2(i),C3(i),C4(i));
end
fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g',size(C1,1),2,4,C1(end),C2(end),C3(end),C4(end));

if split
    for n = 1:nGrps
        fprintf(fid,'\nENDOFSECTION\n       ELEMENT GROUP 2.4.6\n');
        fprintf(fid,'GROUP:\t\t%g ELEMENTS:\t\t%g MATERIAL:\t\t4 NFLAGS:\t\t1\n\t\t\t%s%g\n\t0\n',n,1 + Bodies{n}.Panels.Finish - Bodies{n}.Panels.Start,name,n);
        pans = Bodies{n}.Panels.Start:Bodies{n}.Panels.Finish;
        M3 = zeros(10,floor(numel(pans)/10));
        M3(:) = pans(1:numel(M3));
        M3 = M3';
        M4 = pans(1+numel(M3):end);
        
        for i = 1:size(M3,1)
            fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',M3(i,:));
        end
        
        for i = 1:size(M4,2)
            fprintf(fid,'\t%g',M4(i));
        end
    end
else

        fprintf(fid,'\nENDOFSECTION\n       ELEMENT GROUP 2.4.6\n');
        fprintf(fid,'GROUP:\t\t%g ELEMENTS:\t\t%g MATERIAL:\t\t4 NFLAGS:\t\t1\n\t\t\t%s\n\t0\n',1,nPnls,name);
        pans = 1:nPnls;
        M3 = zeros(10,floor(numel(pans)/10));
        M3(:) = pans(1:numel(M3));
        M3 = M3';
        M4 = pans(1+numel(M3):end);
        
        for i = 1:size(M3,1)
            fprintf(fid,'\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',M3(i,:));
        end
        
        for i = 1:size(M4,2)
            fprintf(fid,'\t%g',M4(i));
        end
    
end

fprintf(fid,'\nENDOFSECTION\n BOUNDARY CONDITIONS 2.4.6\n\t\t\t\t%s\t1\t%g\t0\t6\n','WAKE',numel(W));
fprintf(fid,'\t%g\t%g\t%g\n',W');
fprintf(fid,'ENDOFSECTION\n');



%   Write connectivity of US/LS here


fprintf(fid,'/SURFACE \n');


for i = 1:numBodies

    pans = Bodies{i}.Panels.Start:Bodies{i}.Panels.Finish;
    
    [sx sy] = size(Bodies{i}.Panels.MainPans);
    
    fprintf(fid,'/P %g\t%g\t%g\n',i,sx,sy);
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],pans(Bodies{i}.Panels.MainPans(j,:)));
    end
    
    
    [sx sy] = size(Bodies{i}.Panels.TipInnerUS);
    
    fprintf(fid,'/pTUSI %g\t%g\t%g\n',i,sx,sy);
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],pans(Bodies{i}.Panels.TipInnerUS(j,:)));
    end
    
      
    [sx sy] = size(Bodies{i}.Panels.TipInnerLS);
    
    fprintf(fid,'/pTLSI %g\t%g\t%g\n',i,sx,sy);
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],pans(Bodies{i}.Panels.TipInnerLS(j,:)));
    end
    
    [sx sy] = size(Bodies{i}.Panels.TipInnerUS);
    
    fprintf(fid,'/pTUSO %g\t%g\t%g\n',i,sx,sy);
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],pans(Bodies{i}.Panels.TipOuterUS(j,:)));
    end
    
      
    [sx sy] = size(Bodies{i}.Panels.TipInnerLS);
    
    fprintf(fid,'/pTLSO %g\t%g\t%g\n',i,sx,sy);
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],pans(Bodies{i}.Panels.TipOuterLS(j,:)));
    end

    
    
end
fprintf(fid,'/%%% \n');
fprintf(fid,'\n/LOCALCHORD\t%g \n',nPts);

str = [repmat(' ',nPts,10 - size(num2str(nPts),2)) num2str([1:nPts]','%g    ') repmat('   ',nPts,1) num2str([nC mC],'%-1.12e    ')];


for i = 1:nPts
    fprintf(fid,'/\t%s\n',str(i,:));
    %    fprintf(fid,'\t%g\t%-1.12e\t%-1.12e\t%-1.12e\n',i,M1(i,1),M1(i,2),M1(i,3));
end







fprintf(fid,'/%%% \n');

fprintf(fid,'\n/MESHING \n');


for i = 1:numBodies

    MainSurf = Bodies{i}.N.Global;
    
    [sx sy] = size(MainSurf);
    
    
      fprintf(fid,'/M %g\t%g\t%g\n',i,sx,sy);
    
    
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],MainSurf(j,:));
    end

    

    
    Tip = Bodies{i}.Tip.Inboard.US.N.Global;
    [sx sy] = size(Tip);
    fprintf(fid,'/TUSI %g\t%g\t%g\n',i,sx,sy);
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],Tip(j,:));
    end
    
    Tip = Bodies{i}.Tip.Inboard.LS.N.Global;
    [sx sy] = size(Tip);
    fprintf(fid,'/TLSI %g\t%g\t%g\n',i,sx,sy);
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],Tip(j,:));
    end
    
    Tip = Bodies{i}.Tip.Outboard.US.N.Global;
    [sx sy] = size(Tip);
    fprintf(fid,'/TUSO %g\t%g\t%g\n',i,sx,sy);
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],Tip(j,:));
    end
    
    Tip = Bodies{i}.Tip.Outboard.LS.N.Global;
    [sx sy] = size(Tip);
    fprintf(fid,'/TLSO %g\t%g\t%g\n',i,sx,sy);
    for j = 1:sx
        strs = repmat('\t%g',[1 sy]);
        fprintf(fid,['/' strs '\n'],Tip(j,:));
    end
end




fprintf(fid,'/%%% \n');



fclose(fid);

handles.Bodies = Bodies;