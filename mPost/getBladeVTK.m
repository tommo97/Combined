clear all
close all
clc
fname = 'RunData_000002.mat';
load(fname,'CellPos','CellOms','GambitScale');
Origin = min(CellPos)/GambitScale;

minx = min(CellPos(:,1));
miny = min(CellPos(:,2));
minz = min(CellPos(:,3));

maxx = max(CellPos(:,1));
maxy = max(CellPos(:,2));
maxz = max(CellPos(:,3));

Domain = zeros(1 + maxx-minx,1 + maxy-miny,1 + maxz-minz);

subs = 1 + [CellPos(:,1) - minx CellPos(:,2) - miny CellPos(:,3) - minz];

inds = sub2ind(size(Domain),subs(:,1),subs(:,2),subs(:,3));

Domain(inds) = sqrt(CellOms(:,1).^2 + CellOms(:,2).^2 + CellOms(:,3).^2);

writeVTK(Domain,'Wake100_3sec.full.vtk',1/GambitScale,Origin);
load(fname,'Cp','NParts','AllBodyPoints_x','AllBodyPoints_y','AllBodyPoints_z');
 i = 1;
eval(['load(fname,''TRANS1_' num2str(i-1) '_x'');']);
eval(['load(fname,''TRANS2_' num2str(i-1) '_x'');']);
eval(['load(fname,''TRANS3_' num2str(i-1) '_x'');']);

eval(['load(fname,''TRANS1_' num2str(i-1) '_y'');']);
eval(['load(fname,''TRANS2_' num2str(i-1) '_y'');']);
eval(['load(fname,''TRANS3_' num2str(i-1) '_y'');']);

eval(['load(fname,''TRANS1_' num2str(i-1) '_z'');']);
eval(['load(fname,''TRANS2_' num2str(i-1) '_z'');']);
eval(['load(fname,''TRANS3_' num2str(i-1) '_z'');']);

eval(['TRANS = [TRANS1_' num2str(i-1) '_x TRANS2_' num2str(i-1) '_x TRANS3_' num2str(i-1) '_x; TRANS1_' num2str(i-1) '_y TRANS2_' num2str(i-1) '_y TRANS3_0_y; TRANS1_' num2str(i-1) '_z TRANS2_' num2str(i-1) '_z TRANS3_' num2str(i-1) '_z]']);


for i = 1:NParts
    tic
    eval(['load(fname,''BodyMainPointIDS' num2str(i-1) ''')']);
    eval(['load(fname,''BodySurface' num2str(i-1) ''')']);
    Cps = eval(['Cp(BodySurface' num2str(i-1) ')']);
    
    [ni nj] = size(Cps);
    
    CPs = [Cps(end,end) Cps(end,:) Cps(end,1);Cps(:,end) Cps Cps(:,1);Cps(1,end) Cps(1,:) Cps(1,1)];
    
    
    [m n] = meshgrid(0:nj+1,0:ni+1);
    
    mi = 0.25*(m(1:end-1,1:end-1) + m(1:end-1,2:end) + m(2:end,1:end-1) + m(2:end,2:end));
    ni = 0.25*(n(1:end-1,1:end-1) + n(1:end-1,2:end) + n(2:end,1:end-1) + n(2:end,2:end));
    
    Body{i}.CPi = interp2(m,n,CPs,mi,ni);
    
    
    
    Body{i}.X = eval(['AllBodyPoints_x(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    Body{i}.Y = eval(['AllBodyPoints_y(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    Body{i}.Z = eval(['AllBodyPoints_z(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    Body{i}.Cp = eval(['Cp(BodySurface' num2str(i-1) ')']);
    
    eval(['vtkname = ''BladeData' num2str(i-1) '.vtk''']);
    
   
    
    xi = Body{i}.X * TRANS(1,1) + Body{i}.Y * TRANS(2,1) + Body{i}.Z * TRANS(3,1) + Origin(1);
    yi = Body{i}.X * TRANS(1,2) + Body{i}.Y * TRANS(2,2) + Body{i}.Z * TRANS(3,2);% + Origin(2);
    zi = Body{i}.X * TRANS(1,3) + Body{i}.Y * TRANS(2,3) + Body{i}.Z * TRANS(3,3);% + Origin(3);
    val = Body{i}.CPi;
    
    
    %---get cells
    fprintf('Arranging cells\n')
    ncellX = size(zi,1)-1;
    ncellY = size(zi,2)-1;
    
    %---we save the data (val) as point data, but need to know the cell data to
    %correctly create the cells, so find average val for each cell.
    
        rval = meanRes(val,ncellY,ncellX,1);
        
        %---it may occur that a cell may have a NaN value in just 1-3 vertices and
        %not all. In this case the cell is used so those vertices need an
        %elevation.
        %zi= fillElevation(zi,ncellX,ncellY,0,rval);
        index  = zeros(ncellX+1,ncellY+1);
        k      = 0;
        
        %---index matvalx.
        for j = 1:ncellY+1
            for i = 1:ncellX+1
                if isfinite(zi(i,j))
                    index(i,j) = k;
                    k=k+1;
                end
            end
        end
        
        
        %---work through cells. If there is no value for a cell, discard it.
        cells = cell(ncellX,ncellY);
        for i = 1:ncellX
            for j = 1:ncellY
                if(isfinite(rval(i,j)))
                    cP1=index(i,j);
                    cP2=index(i+1,j);
                    cP3=index(i+1,j+1);
                    cP4=index(i,j+1);
                    cells(i,j) = {[4 cP1 cP2 cP3 cP4]};
                end
            end
        end
        
        %---shuffle the data to vtk format.
        cells=cells(:);
        cells=cell2mat(cells);
        ncells = size(cells,1);
        xi(isnan(zi))=[];
        yi(isnan(zi))=[];
        zi(isnan(zi))=[];
        points = [xi(:) yi(:) zi(:)];
        npoints = length(zi(:));
        val(isnan(val))=[];
        
        %---wvalte vtk file. remeber to change this to binary for faster reading
        fprintf('Wvalting VTK file\n')
        vtk_file=fopen(vtkname,'w');
        
        fprintf(vtk_file,'# vtk DataFile Version 3.0\n3D export of map data\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d float\n',npoints);
        fprintf(vtk_file,'%i %i %i\n',points.');
        fprintf(vtk_file,'CELLS %i %i\n',ncells,ncells*5);
        fprintf(vtk_file,[repmat('%i ', 1, size(cells, 2)), '\n'],cells.');
        fprintf(vtk_file,'CELL_TYPES %i\n',ncells);
        fprintf(vtk_file,'%i\n',(ones(ncells,1)*9));
        fprintf(vtk_file,'POINT_DATA %u \nSCALARS data float 1 \nLOOKUP_TABLE default \n ',npoints);
        fprintf(vtk_file,[repmat('%.2f ', 1, size(val, 2)), '\n'],val);
        fclose all;
        fprintf('VTK file written\n')
        toc
        
        
    end