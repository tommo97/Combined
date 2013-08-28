function ConvertVTK(infname)
fname = infname;
vtkname = [infname(1:(length(infname) - 3)) 'vtk'];
dirname = [infname(1:(length(infname) - 4)) '_VTKdata'];
system(['rm -rf ' dirname]);
system(['mkdir ' dirname]);
load(fname,'CellPos','CellOms','GambitScale','Time','NBodies','NumTransVars', 'TransVars_x','TransVars_y','TransVars_z');
Origin = min(CellPos)/GambitScale;

minx = min(CellPos(:,1));
miny = min(CellPos(:,2));
minz = min(CellPos(:,3));

maxx = max(CellPos(:,1));
maxy = max(CellPos(:,2));
maxz = max(CellPos(:,3));
  

for i = 1:NumTransVars

    Domain = zeros(1 + maxx-minx,1 + maxy-miny,1 + maxz-minz);
    
    subs = 1 + [CellPos(:,1) - minx CellPos(:,2) - miny CellPos(:,3) - minz];
    
    inds = sub2ind(size(Domain),subs(:,1),subs(:,2),subs(:,3));
        
    Domain(inds) = sqrt(TransVars_x(:,i).^2 + TransVars_y(:,i).^2 + TransVars_z(:,i).^2);
    
    wakeVTKname = [dirname '/Wake_Scale=' num2str(GambitScale) '_Time=' num2str(Time) '_' num2str(i) '_' vtkname];
    
    
    writeDomainVTK(Domain,wakeVTKname,1/GambitScale,Origin);
end
load(fname,'Cp','NParts','AllBodyPoints_x','AllBodyPoints_y','AllBodyPoints_z');



for i = 1:NParts
    tic
   
    eval(['load(fname,''BodyMainPointIDS' num2str(i-1) ''')']);
    eval(['load(fname,''BodySurface' num2str(i-1) ''')']);

    Cps = eval(['Cp(BodySurface' num2str(i-1) ')']);
    
    [ni, nj] = size(Cps);
    
    CPs = [Cps(end,end) Cps(end,:) Cps(end,1);Cps(:,end) Cps Cps(:,1);Cps(1,end) Cps(1,:) Cps(1,1)];
    
    
    [m, n] = meshgrid(0:nj+1,0:ni+1);
    
    mi = 0.25*(m(1:end-1,1:end-1) + m(1:end-1,2:end) + m(2:end,1:end-1) + m(2:end,2:end));
    ni = 0.25*(n(1:end-1,1:end-1) + n(1:end-1,2:end) + n(2:end,1:end-1) + n(2:end,2:end));
    
    val = interp2(m,n,CPs,mi,ni);
    
    
    
    xi = eval(['AllBodyPoints_x(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    yi = eval(['AllBodyPoints_y(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    zi = eval(['AllBodyPoints_z(BodyMainPointIDS' num2str(i-1) ')/GambitScale' ]);
    str = ['vtkname = [''' dirname '/BodyData_Scale=' num2str(GambitScale) '_Time=' num2str(Time) '_Body=' num2str(i-1) '.vtk'']'];
    eval(str);
  
    
    
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
    

function [] = writeDomainVTK(vol,vtkfile,spacing,origin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions
volinfo = whos('vol');

sz = volinfo.size;
X = sz(1); Y = sz(2); Z = 1;
if( length(sz) == 3 )
  Z = sz(3);
end

% open file (OBS! big endian format)
fid = fopen(vtkfile,'w','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.0');
fprintf(fid, '%s\n', 'created by writeVTK (Matlab implementation by Erik Vidholm)');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', X, ' ', Y, ' ', Z);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', origin(1), ' ', origin(2), ' ', origin(3)); 
fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', spacing, ' ', spacing, ' ', spacing); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', X*Y*Z);

tp = volinfo.class;
if( strcmp(tp, 'uint8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_char');
elseif( strcmp(tp, 'int8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data char');
elseif( strcmp(tp, 'uint16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_short');
elseif( strcmp(tp, 'int16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data short');
elseif( strcmp(tp, 'uint32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_int');
elseif( strcmp(tp, 'int32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data int');
elseif( strcmp(tp, 'single') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data float');
elseif( strcmp(tp, 'double') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data double');
end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,vol,tp);

% close file
fclose(fid);