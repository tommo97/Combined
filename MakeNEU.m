function MakeNEU(Bodies,name,split)
dir1 = './neu_files/';
fname = [dir1 name '.neu'];
fid = fopen(fname, 'wt');
Empty = 0;
nGrps = length(Bodies);
nPts = 0;
nPnls = 0;

if split
    nBCs = nGrps;
else
    nBCs = 1;
end
PtsX = [];
C1 = [];
C2 = [];
C3 = [];
C4 = [];
W = [];
F = [];
f{1} = 2;
f{2} = 4;
for i = 1:nGrps
    Bodies{i}.Panels.Start = nPnls + 1;
    nPts = nPts + Bodies{i}.nPts;
    nPnls = nPnls + Bodies{i}.nPnls;
    Bodies{i}.Panels.Finish = nPnls;
    disp([Bodies{i}.Panels.Start Bodies{i}.Panels.Finish]);
    
    PtsX = [PtsX;[Bodies{i}.X Bodies{i}.Y Bodies{i}.Z]];
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



disp(['Writing ' num2str(nGrps) ' groups with ' num2str(nPts) ' points and ' num2str(nPnls) ' panels.'])


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
disp(W);
fprintf(fid,'ENDOFSECTION\n');
fclose(fid);