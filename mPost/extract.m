function [XI,YI,ZI,VI, VIx, VIy, VIz] = extract(fname,val, scale)
load(fname);
data = Domain;
clear Domain;
subs = data(:,1:3);
subs(:,1) = subs(:,1) - min(data(:,1));
subs(:,2) = subs(:,2) - min(data(:,2));
subs(:,3) = subs(:,3) - min(data(:,3));
subs = subs+1;

V = zeros(max(subs));
Vx = zeros(max(subs));
Vy = zeros(max(subs));
Vz = zeros(max(subs));

disp(['Domain size: ' num2str(size(V)) '; i.e. ' num2str(numel(V)) ' cells' ]);
disp(['Number of Vorticity Cells: ' num2str(length(subs(:,1)))]);
disp(['Occupancy Ratio: ' num2str(length(subs(:,1))/numel(V))]);


[X Y Z] = meshgrid([min(data(:,2)):1:max(data(:,2))]/scale,...
    [min(data(:,1)):1:max(data(:,1))]/scale,...
    [min(data(:,3)):1:max(data(:,3))]/scale);





num_vort = (size(data,2) - 6)/3;
colour{1} = [1 0 0];
colour{2} = [0 1 0];
colour{3} = [0 0 1];

for q = 1:num_vort
    %val = val - q;
    cols = [7 + (q-1)*3:7+(q*3)-1];
    vort = data(:,cols);
    
    vals = sqrt(vort(:,1).*vort(:,1) + vort(:,2).*vort(:,2) + vort(:,3).*vort(:,3));
    
    max(vals(:))
    
    for i = 1:size(subs,1)
        V(subs(i,1),subs(i,2), subs(i,3)) =  vals(i);
        Vx(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,1);
        Vy(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,2);
        Vz(subs(i,1),subs(i,2),subs(i,3)) =  vort(i,3);
    end
    %Y = Y - min(Y(:));
    %V = V.^3;
    %V = smooth3(V);%,'gaussian',[1 1 1],0.65);   

        VI = V;
        VIx = Vx;
        VIy = Vy;
        VIz = Vz;
        XI = X;
        YI = Y;
        ZI = Z;
        
 
    
    
%     [x,y,z,v] = subvolume(XI,YI,ZI,VI,[nan,nan,nan,nan,nan,0]);
%     p1 = patch(isosurface(x,y,z,v, val),...
%         'FaceColor',colour{q},'EdgeColor','none');
%     isonormals(x,y,z,v,p1);
%     %v = v.^(1/3);
%     %v = smooth3(v,'gaussian',[5 5 5],0.65);   
%     p2 = patch(isocaps(x,y,z,v, val),...
%         'FaceColor','interp','EdgeColor','none','DiffuseStrength',1.0);
    
 
%     p = patch(isosurface(XI,YI,ZI,VI,val));%,'FaceColor',colour{q},...
%         %'EdgeColor','none','FaceLighting','phong');
%     isonormals(XI,YI,ZI,VI,p);
%     isocolors(XI,YI,ZI,VIx,p);
%     set(p,'FaceColor','interp','EdgeColor','none')
end;