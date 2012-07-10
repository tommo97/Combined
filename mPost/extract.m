function [XI,YI,ZI,VI, VIx, VIy, VIz, x,y,z, U, V, W] = extract(fname,val, scale)
load(fname,'CellPos');

subs = CellPos(:,1:3);
subs(:,1) = subs(:,1) - min(CellPos(:,1));
subs(:,2) = subs(:,2) - min(CellPos(:,2));
subs(:,3) = subs(:,3) - min(CellPos(:,3));
subs = subs+1;

VI = zeros(max(subs));
U = [];%VI;
V = [];%VI;
W = [];%VI;
x = [];%VI;
y = [];%VI;
z = [];%VI;
Vx = [];%zeros(max(subs));
Vy = [];%zeros(max(subs));
Vz = [];%zeros(max(subs));
VIx = [];
VIy = [];
VIz = [];

disp(['Domain size: ' num2str(size(VI)) '; i.e. ' num2str(numel(VI)) ' cells' ]);
disp(['Number of Vorticity Cells: ' num2str(length(subs(:,1)))]);
disp(['Occupancy Ratio: ' num2str(length(subs(:,1))/numel(VI))]);


[XI YI ZI] = meshgrid([min(CellPos(:,2)):1:max(CellPos(:,2))]/scale,...
    [min(CellPos(:,1)):1:max(CellPos(:,1))]/scale,...
    [min(CellPos(:,3)):1:max(CellPos(:,3))]/scale);


sz = size(CellPos,2);
artcols = (sz-3):sz;
viscols = artcols-6;
load(fname,'CellViscDeriv');
load(fname,'CellOms');
load(fname,'CellVel');
Data = CellOms;


data = sqrt(Data(:,1).*Data(:,1) + Data(:,2).*Data(:,2) + Data(:,3).*Data(:,3));

%data = CellVel(:,1);

for i = 1:size(subs,1)
    VI(subs(i,1),subs(i,2), subs(i,3)) =  data(i);
end

return
num_vort = (size(Domain,2) - 6)/3;
colour{1} = [1 0 0];
colour{2} = [0 1 0];
colour{3} = [0 0 1];








for q = 1:num_vort
    %val = val - q;
    cols = [18:20];%[7 + (q-1)*3:7+(q*3)-1];
    vort = Domain(:,cols);
    
    vals = sqrt(vort(:,1).*vort(:,1) + vort(:,2).*vort(:,2) + vort(:,3).*vort(:,3));
    
    max(vals(:))
    
    for i = 1:size(subs,1)
        VI(subs(i,1),subs(i,2), subs(i,3)) =  vals(i);
        %Vx(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,1);
        %Vy(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,2);
        %Vz(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,3);
    end
    %Y = Y - min(Y(:));
    %VI = VI.^3;
    %VI = smooth3(VI);%,'gaussian',[1 1 1],0.65);
    
    
    %VIx = Vx;
    %VIy = Vy;
    %VIz = Vz;
    
    %XI = X;
    %YI = Y;
    %ZI = Z;
    
    
    
    
    %     [x,y,z,VI] = subvolume(XI,YI,ZI,VI,[nan,nan,nan,nan,nan,0]);
    %     p1 = patch(isosurface(x,y,z,VI, val),...
    %         'FaceColor',colour{q},'EdgeColor','none');
    %     isonormals(x,y,z,VI,p1);
    %     %VI = VI.^(1/3);
    %     %VI = smooth3(VI,'gaussian',[5 5 5],0.65);
    %     p2 = patch(isocaps(x,y,z,VI, val),...
    %         'FaceColor','interp','EdgeColor','none','DiffuseStrength',1.0);
    
    
    %     p = patch(isosurface(XI,YI,ZI,VI,val));%,'FaceColor',colour{q},...
    %         %'EdgeColor','none','FaceLighting','phong');
    %     isonormals(XI,YI,ZI,VI,p);
    %     isocolors(XI,YI,ZI,VIx,p);
    %     set(p,'FaceColor','interp','EdgeColor','none')
end;
return
try
    Domain = [CellPositions_x CellPositions_y CellPositions_z CellVelocities_x CellVelocities_y CellVelocities_z];
    subs = Domain(:,1:3);
    subs(:,1) = subs(:,1) - min(Domain(:,1));
    subs(:,2) = subs(:,2) - min(Domain(:,2));
    subs(:,3) = subs(:,3) - min(Domain(:,3));
    subs = subs+1;
    
    for i = 1:size(subs,1)
        x(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,1);
        y(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,2);
        z(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,3);        
        U(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,4);
        V(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,5);
        W(subs(i,1),subs(i,2),subs(i,3)) =  Domain(i,6);
    end
catch
    
end