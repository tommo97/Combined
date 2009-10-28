clear all;
C = [1 1 1];
FACES{1} = C + [1 0 0]; %E
FACES{2} = C - [1 0 0]; %W
FACES{3} = C + [0 1 0]; %N
FACES{4} = C - [0 1 0]; %S
FACES{5} = C + [0 0 1]; %T
FACES{6} = C - [0 0 1]; %B

EDGES{1}  = C + [1   1  0];  %NEC
EDGES{2}  = C + [-1  1  0];  %NWC
EDGES{3}  = C + [1  -1  0];  %SEC
EDGES{4}  = C + [-1 -1  0];  %SWC
EDGES{5}  = C + [0   1  1];  %NT
EDGES{6}  = C + [0   1 -1];  %NB
EDGES{7}  = C + [0  -1  1];  %ST
EDGES{8}  = C + [0  -1 -1];  %SB
EDGES{9}  = C + [1   0  1];  %ET
EDGES{10} = C + [1   0 -1];  %EB
EDGES{11} = C + [-1  0  1];  %WT
EDGES{12} = C + [-1  0 -1];  %WB

CORNERS{1} = C + [ 1  1  1];    %NET
CORNERS{2} = C + [-1  1  1];    %NWT
CORNERS{3} = C + [1  -1  1];    %SET
CORNERS{4} = C + [-1 -1  1];    %SWT
CORNERS{5} = C + [ 1  1 -1];    %NEB
CORNERS{6} = C + [-1  1 -1];    %NWB
CORNERS{7} = C + [1  -1 -1];    %SEB
CORNERS{8} = C + [-1 -1 -1];    %SWB


close all
hold all
for i = 1:6
    scatter3(FACES{i}(1),FACES{i}(2),FACES{i}(3),100,'filled');
    
    disp(['if (ISA[' num2str(FACES{i}(1)) '][' num2str(FACES{i}(2)) '][' num2str(FACES{i}(3)) ']) FACES[' num2str(i-1) '] = ISA[' num2str(FACES{i}(1)) '][' num2str(FACES{i}(2)) '][' num2str(FACES{i}(3)) ']->Omega;'])  
    
end
for i = 1:8
    scatter3(CORNERS{i}(1),CORNERS{i}(2),CORNERS{i}(3),100);
        disp(['if (ISA[' num2str(CORNERS{i}(1)) '][' num2str(CORNERS{i}(2)) '][' num2str(CORNERS{i}(3)) ']) CORNERS[' num2str(i-1) '] = ISA[' num2str(CORNERS{i}(1)) '][' num2str(CORNERS{i}(2)) '][' num2str(CORNERS{i}(3)) ']->Omega;'])  
end
for i = 1:12
    scatter3(EDGES{i}(1),EDGES{i}(2),EDGES{i}(3),100,'s');
        disp(['if (ISA[' num2str(EDGES{i}(1)) '][' num2str(EDGES{i}(2)) '][' num2str(EDGES{i}(3)) ']) EDGES[' num2str(i-1) '] = ISA[' num2str(EDGES{i}(1)) '][' num2str(EDGES{i}(2)) '][' num2str(EDGES{i}(3)) ']->Omega;'])  
end