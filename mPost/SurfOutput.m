



TRANS = [TRANS1_0_x TRANS1_0_y TRANS1_0_z; TRANS2_0_x TRANS2_0_y TRANS2_0_z; TRANS3_0_x TRANS3_0_y TRANS3_0_z];
TRANS_INV = TRANS';
infile = '64x64.mat'
load(infile)
set(0,'defaultaxesposition',[0 0 1 1],'defaulttextinterpreter','none')
hold all
nc=24;
for i = 1:NParts
    Cps = eval(['Cp(BodySurface' num2str(i-1) ')']);
    
    [ni nj] = size(Cps);
    
    CPs = [Cps(end,end) Cps(end,:) Cps(end,1);Cps(:,end) Cps Cps(:,1);Cps(1,end) Cps(1,:) Cps(1,1)];
    
    
    [m n] = meshgrid(0:nj+1,0:ni+1);
    
    mi = 0.25*(m(1:end-1,1:end-1) + m(1:end-1,2:end) + m(2:end,1:end-1) + m(2:end,2:end));
    ni = 0.25*(n(1:end-1,1:end-1) + n(1:end-1,2:end) + n(2:end,1:end-1) + n(2:end,2:end));
    
    Body{i}.CPi = interp2(m,n,CPs,mi,ni);
    
    
    
    Body{i}.Xt = eval(['AllBodyPoints_x(BodyMainPointIDS' num2str(i-1) ')' ]);
    Body{i}.Yt = eval(['AllBodyPoints_y(BodyMainPointIDS' num2str(i-1) ')' ]);
    Body{i}.Zt = eval(['AllBodyPoints_z(BodyMainPointIDS' num2str(i-1) ')' ]);
    
    Body{i}.Y = Body{i}.Xt * TRANS_INV(1,1) + Body{i}.Yt * TRANS_INV(2,1) + Body{i}.Zt * TRANS_INV(3,1);
    Body{i}.X = Body{i}.Xt * TRANS_INV(1,2) + Body{i}.Yt * TRANS_INV(2,2) + Body{i}.Zt * TRANS_INV(3,2);
    Body{i}.Z = Body{i}.Xt * TRANS_INV(1,3) + Body{i}.Yt * TRANS_INV(2,3) + Body{i}.Zt * TRANS_INV(3,3);
    
    Body{i}.Cp = eval(['Cp(BodySurface' num2str(i-1) ')']);
    

    
    %surf(Body{i}.X,Body{i}.Y,Body{i}.Z,Body{i}.CPi);
    %surf(Body{i}.X,Body{i}.Y,Body{i}.Z,'FaceColor',[0.95 0.95 0.95], 'LineWidth',0.25);
    
    data = load(infile);
    h = figure();
    Body{i}.CPi(Body{i}.CPi < -2.5) = -2.5;
    contours = linspace(min(Body{i}.CPi(:)) ,1,nc);
    [C] = contourf(Body{i}.CPi,contours,'LineStyle','none');
    
    %colormap(colorGray)
    
    %[C,h] = contour(R,[0.211 0.411 0.611 0.811 0.961],'LineWidth',4);
    
    cax = caxis;
    
    M = blue2red(min(contours),max(contours),nc);
    
    colormap(jet);
    
    axis off
    im = getframe();
    close(h)
    [A,map] = rgb2ind(im.cdata,64);
    s = surface(Body{i}.X,Body{i}.Y,Body{i}.Z,flipdim(im.cdata,1),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    axis equal tight
    
    
    lighting phong
    
    %surface(Body{i}.X,Body{i}.Y,Body{i}.Z,'FaceColor','none','EdgeColor',[.5 .5 .5],'LineWidth',0.1);
    
    Xi = interp2(Body{i}.X,5);
    Yi = interp2(Body{i}.Y,5);
    Zi = interp2(Body{i}.Z,5);
    
    R = Yi;
    R = R/max(R(:));
    
    
    %contourz(Body{i}.X,Body{i}.Y,Body{i}.Z,Body{i}.CPi,linspace(-2.5,1,nc),'k');

    %handles = contourz(Body{i}.X,Body{i}.Y,Body{i}.Z,Body{i}.Y/min(Body{i}.Y(:)),[0.211 0.411 0.611 0.811 0.961],'k');

    
    
    
    Cpi = interp2(Body{i}.Cp,5);
    

    contourz(Body{i}.X,Body{i}.Y,Body{i}.Z,Body{i}.CPi,contours,'k');
    
    axis off
    
    
end

axis equal tight

axis on
box on
grid on
view(3)

axis equal tight;
ax = axis;

voxel([0 ax(3) 0],[1 1 1]./scale,1,'g');

lighting phong; camlight right;

set(gca,'Projection','perspective')

set(gca,'Ztick',linspace(-0.4,0.4,5));
set(gca,'Xtick',linspace(-0.4,0.4,5));

xlabel('$x$')
ylabel('$y$')
zl = zlabel('$z$')
set(zl,'Rotation',0)
caxis([-4 1])
%cb = colorbar;
%set(cb,'Ytick',linspace(-4,1,6))
zoom(1.1);
zoom(1.1);
set(gca,'CameraPosition',[-3.9699   -4.7822    3.7110]);
set(gca,'CameraTarget',[ 0.0308    0.3329   -0.0402]);
set(0,'defaulttextinterpreter','none','defaultaxesposition',[0.10    0.10    .89    .8]);
