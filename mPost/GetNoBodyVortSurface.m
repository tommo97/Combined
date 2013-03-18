clear all; clear mex; clc; close all

val = .1;
col(1) = 'r';
col(2) = 'b';


IDs = [5];
RANGE = 1:length(IDs);
%matlabpool close force local
%matlabpool(7)

for i = RANGE
    
    
    ID = IDs(i);
    files = dir(['RunData*0' num2str(ID) '.mat']);
    
    s = size(files,1);
    fname = files(s).name;
    Data(i) = load(fname,'CellOms','GambitScale','CellPos','TransVars_x','TransVars_y','TransVars_z');
end

for i = RANGE
    
    
    Data(i).subs = Data(i).CellPos(:,1:3);
    Data(i).subs(:,1) = Data(i).subs(:,1) - min(Data(i).CellPos(:,1));
    Data(i).subs(:,2) = Data(i).subs(:,2) - min(Data(i).CellPos(:,2));
    Data(i).subs(:,3) = Data(i).subs(:,3) - min(Data(i).CellPos(:,3));
    Data(i).subs = Data(i).subs+1;
    Data(i).inds = sub2ind(max(Data(i).subs),Data(i).subs(:,1),Data(i).subs(:,2),Data(i).subs(:,3));
    Data(i).do_multiple_vorts = false;
    
    Data(i).VI.mag = zeros(max(Data(i).subs));
    Data(i).VI.mag(Data(i).inds) = sqrt(Data(i).CellOms(:,1).^2 + Data(i).CellOms(:,2).^2 + Data(i).CellOms(:,3).^2);
    
    Data(i).VI.max = max(Data(i).VI.mag(:));
    Data(i).VI.mean = mean(Data(i).VI.mag(:));
    Data(i).VI.total = sum(Data(i).VI.mag(:));
    Data(i).NumTransVars = size(Data(i).TransVars_x,2);
    
    for j = 1:Data(i).NumTransVars
        Data(i).V(j).mag = zeros(max(Data(i).subs));
        Data(i).V(j).mag(Data(i).inds) = sqrt(Data(i).TransVars_x(:,j).^2 + Data(i).TransVars_y(:,j).^2 + Data(i).TransVars_z(:,j).^2);
        
    end
    disp(['Domain size: ' num2str(size(Data(i).VI.mag)) '; i.e. ' num2str(numel(Data(i).VI.mag)) ' cells' ]);
    disp(['Number of Vorticity Cells: ' num2str(length(Data(i).subs(:,1)))]);
    disp(['Occupancy Ratio: ' num2str(length(Data(i).subs(:,1))/numel(Data(i).VI.mag))]);
    
    
    [Data(i).XI Data(i).YI Data(i).ZI] = meshgrid([min(Data(i).CellPos(:,2)):1:max(Data(i).CellPos(:,2))]/Data(i).GambitScale,...
        [min(Data(i).CellPos(:,1)):1:max(Data(i).CellPos(:,1))]/Data(i).GambitScale,...
        [min(Data(i).CellPos(:,3)):1:max(Data(i).CellPos(:,3))]/Data(i).GambitScale);
    
    
    %Data(i).fig = [];
    %Data(i).p = []; Data(i).q = []; Data(i).r = [];
end

for i = RANGE
    
    figure
    for j = 1:Data(i).NumTransVars
        Data(i).p(j) = patch(isosurface(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,val));
        isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,Data(i).p(j));
        set(Data(i).p(j),'FaceColor',col(j),'EdgeColor','none');%,'faceAlpha',1)
        
%         Data(i).q(j) = patch(isosurface(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,0.125*Data(i).VI.mean));
%         isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,Data(i).q(j));
%         set(Data(i).q(j),'FaceColor',col(j),'EdgeColor','none','faceAlpha',.1)
%         
%         Data(i).r(j) = patch(isosurface(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,0.0625*Data(i).VI.mean));
%         isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,Data(i).r(j));
%         set(Data(i).r(j),'FaceColor',col(j),'EdgeColor','none','faceAlpha',.05)
     end
    
    
    
    
    axis equal tight;
    lighting phong
    box on
    grid on
    view(3)
    axis([-4 4 -2 2 -5 -0]);
    axis([-4 4 -2 2 -4 4]);
    camlight('infinite');
    
    set(gcf,'Color',[1,1,1],'Renderer','OpenGL');
    set(gca,'Projection','perspective')
    axis off
    zoom(1.75)
    drawnow
    
    print(gcf,['normal-ring-' num2str(i) '.jpg'],'-r300','-djpeg')
end

axis tight



%matlabpool close