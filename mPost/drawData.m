clear all; close all; clc;
col(1) = 'r';
col(2) = 'b';
load Data.mat
for i = size(Data,2)
    
    figure
    set(gcf,'renderer','OpenGL');
    for j = 1:Data(i).NumTransVars
        p(j) = patch(Data(i).p(j)); 
        isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,p(j));
        set(p(j),'FaceColor',col(j),'EdgeColor','none','faceAlpha',1)
        q(j) = patch(Data(i).q(j)); 
        isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,q(j));
        set(q(j),'FaceColor',col(j),'EdgeColor','none','faceAlpha',.1)
        r(j) = patch(Data(i).r(j)); 
        isonormals(Data(i).XI,Data(i).YI,Data(i).ZI,Data(i).V(j).mag,r(j));
        set(r(j),'FaceColor',col(j),'EdgeColor','none','faceAlpha',.05)
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
        zoom(2)
        drawnow
%     
%     print(gcf,['oblique-ring-' num2str(i) '.jpg'],'-r300','-djpeg')
end