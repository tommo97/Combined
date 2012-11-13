clear all
figure

load('RunData_000020.mat')


c = 10;
try
    mins = min(CellPos);
catch
    CellPos = Domain(:,1:3);
    CellOms = Domain(:,4:6);
    mins = min(CellPos);
end
CellInds(:,1) = CellPos(:,1) - mins(1) + 1;
CellInds(:,2) = CellPos(:,2) - mins(2) + 1;
CellInds(:,3) = CellPos(:,3) - mins(3) + 1;
maxs = max(CellInds);


Domain = zeros(maxs);


linearInd = sub2ind(maxs, CellInds(:,1), CellInds(:,2), CellInds(:,3));

Oms = sqrt(sum(CellOms.^2,2));

Domain(linearInd) = Oms;

isosurface(Domain,c);

axis equal

[X,Y,Z] = meshgrid(1:maxs(2),1:maxs(1),1:maxs(3));
[F,V,col] = MarchingCubes(X,Y,Z,Domain,c);

tic

figure
set(gcf,'Renderer','OpenGL')
hold on
%   Find all voxels corner values

[l,m,n] = size(Domain);
count = 0;
cnt = 1;
for i = 2:l-1
    for j = 2:m-1
        for k = 2:n-1
            
            NET = 0.125*(Domain(i,j,k)+Domain(i+1,j,k)+Domain(i+1,j+1,k)+Domain(i,j+1,k)+...
                Domain(i,j,k+1)+Domain(i+1,j,k+1)+Domain(i+1,j+1,k+1)+Domain(i,j+1,k+1));
            NEB = 0.125*(Domain(i,j,k)+Domain(i+1,j,k)+Domain(i+1,j+1,k)+Domain(i,j+1,k)+...
                Domain(i,j,k-1)+Domain(i+1,j,k-1)+Domain(i+1,j+1,k-1)+Domain(i,j+1,k-1));
            
            
            SET = 0.125*(Domain(i,j,k)+Domain(i+1,j,k)+Domain(i+1,j-1,k)+Domain(i,j-1,k)+...
                Domain(i,j,k+1)+Domain(i+1,j,k+1)+Domain(i+1,j-1,k+1)+Domain(i,j-1,k+1));
            SEB = 0.125*(Domain(i,j,k)+Domain(i+1,j,k)+Domain(i+1,j+1,k)+Domain(i,j+1,k)+...
                Domain(i,j,k-1)+Domain(i+1,j,k-1)+Domain(i+1,j+1,k-1)+Domain(i,j+1,k-1));
            
            NWT = 0.125*(Domain(i,j,k)+Domain(i-1,j,k)+Domain(i-1,j+1,k)+Domain(i,j+1,k)+...
                Domain(i,j,k+1)+Domain(i-1,j,k+1)+Domain(i-1,j+1,k+1)+Domain(i,j+1,k+1));
            NWB = 0.125*(Domain(i,j,k)+Domain(i-1,j,k)+Domain(i-1,j+1,k)+Domain(i,j+1,k)+...
                Domain(i,j,k-1)+Domain(i-1,j,k-1)+Domain(i-1,j+1,k-1)+Domain(i,j+1,k-1));
            
            SWT = 0.125*(Domain(i,j,k)+Domain(i-1,j,k)+Domain(i-1,j-1,k)+Domain(i,j-1,k)+...
                Domain(i,j,k+1)+Domain(i-1,j,k+1)+Domain(i-1,j-1,k+1)+Domain(i,j-1,k+1));
            SWB = 0.125*(Domain(i,j,k)+Domain(i-1,j,k)+Domain(i-1,j-1,k)+Domain(i,j-1,k)+...
                Domain(i,j,k-1)+Domain(i-1,j,k-1)+Domain(i-1,j-1,k-1)+Domain(i,j-1,k-1));
            
            if ~(NET==NEB==SET==SEB==NWT==NWB==SWT==SWB)
                
                sNE = (NET-c)*(c-NEB);
                sNW = (NWT-c)*(c-NWB);
                sSE = (SET-c)*(c-SEB);
                sSW = (SWT-c)*(c-SWB);
                
                sNT = (NET-c)*(c-NWT);
                sNB = (NEB-c)*(c-NWB);
                sST = (SET-c)*(c-SWT);
                sSB = (SEB-c)*(c-SWB);
                
                sET = (NET-c)*(c-SET);
                sEB = (NEB-c)*(c-SEB);
                sWT = (NWT-c)*(c-SWT);
                sWB = (NWB-c)*(c-SWB);
                
                S(1,1) = sNE; S(1,2) = sNW; S(1,3) = sSE; S(1,4) = sSW;
                S(2,1) = sNT; S(2,2) = sNB; S(2,3) = sST; S(2,4) = sSB;
                S(3,1) = sET; S(3,2) = sEB; S(3,3) = sWT; S(3,4) = sWB;
                
                S(S==0) = 1;
                S(S<0) = 0;
                
                if (sum(S(:))>0)
                    cnt = cnt+1;
                    h = cube([i,j,k]-.5, 1);
                    set(h,'FaceColor','r','LineStyle','none');%,'FaceAlpha',0.5*iinds(i)/64);
                end
            end
            count = count + 1;
            
        end
    end
    disp(num2str([count numel(Domain)]));
end
toc
            