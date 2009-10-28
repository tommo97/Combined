function extract(fname,val)
data = dlmread(fname);
subs = data(:,1:3);
subs(:,1) = subs(:,1) - min(data(:,1));
subs(:,2) = subs(:,2) - min(data(:,2));
subs(:,3) = subs(:,3) - min(data(:,3));
subs = subs+1;
max(subs)
V = zeros(max(subs));
[X Y Z] = meshgrid([min(data(:,2)):1:max(data(:,2))],...
    [min(data(:,1)):1:max(data(:,1))],...
    [min(data(:,3)):1:max(data(:,3))]);


n = 1;

if n ~= 1
[XI YI ZI] = meshgrid([min(data(:,2)):1/n:max(data(:,2))],...
    [min(data(:,1)):1/n:max(data(:,1))],...
    [min(data(:,3)):1/n:max(data(:,3))]);
VI = zeros(size(XI));
end

num_vort = (size(data,2) - 6)/3;
colour{1} = [1 0 0];
colour{2} = [1 1 0];

for q = 1:num_vort
    %val = val - q;
    cols = [7 + (q-1)*3:7+(q*3)-1];
    vort = data(:,cols);
    
    vals = sqrt(vort(:,1).*vort(:,1) + vort(:,2).*vort(:,2) + vort(:,3).*vort(:,3));
    max(vals(:))
    
    for i = 1:size(subs,1)
        V(subs(i,1),subs(i,2), subs(i,3)) =  vals(i);
    end
    Y = Y - min(Y(:));

    V = smooth3(V,'gaussian',[5 5 5],0.65);   
    if n ~= 1
        YI = YI - min(YI(:));
        VI =  interp3(X,Y,Z,V,XI,YI,ZI,'cubic');
    else
        VI = V;
        XI = X;
        YI = Y;
        ZI = Z;
    end
 
    p = patch(isosurface(XI,YI,ZI,VI,val),'FaceColor',colour{q},...
        'EdgeColor','none','FaceLighting','phong');
    isonormals(XI,YI,ZI,VI,p);
end;