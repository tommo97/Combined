function handles = PostPlot(handles)

%%  Begin Postprocessing
for i = 1:size(handles.PostComp.Bodies,1);
    inds = [1:max(handles.PostComp.Bodies{i}.Panels.MainPans(:))]';

    handles.PostComp.Bodies{i}.CPdist = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    
    handles.PostComp.Bodies{i}.Fxdist = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    
    handles.PostComp.Bodies{i}.Fydist = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
        
    handles.PostComp.Bodies{i}.Fzdist = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));

    handles.PostComp.Bodies{i}.GammaDist = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    

    handles.PostComp.Bodies{i}.CPdist(:) = handles.PostComp.Bodies{i}.Faces.Cpress(inds);
    handles.PostComp.Bodies{i}.Fxdist(:) = handles.PostComp.Bodies{i}.Faces.dF(inds,1);
    handles.PostComp.Bodies{i}.Fydist(:) = handles.PostComp.Bodies{i}.Faces.dF(inds,2);
    handles.PostComp.Bodies{i}.Fzdist(:) = handles.PostComp.Bodies{i}.Faces.dF(inds,3);
    handles.PostComp.Bodies{i}.GammaDist(:) = handles.PostComp.Bodies{i}.Faces.Gamma(inds);
    
    handles.PostComp.Bodies{i}.Panels.CP.Body.n = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    handles.PostComp.Bodies{i}.Panels.CP.Body.n(:) = handles.PostComp.Bodies{i}.n(inds);
    handles.PostComp.Bodies{i}.Panels.CP.Body.x = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    handles.PostComp.Bodies{i}.Panels.CP.Body.x(:) = handles.PostComp.Bodies{i}.Faces.CP.Body(inds,1);
    
    handles.PostComp.Bodies{i}.Panels.CP.Body.y = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    handles.PostComp.Bodies{i}.Panels.CP.Body.y(:) = handles.PostComp.Bodies{i}.Faces.CP.Body(inds,2);
    
    handles.PostComp.Bodies{i}.Panels.CP.Body.z = zeros(size(handles.PostComp.Bodies{i}.Panels.MainPans));
    handles.PostComp.Bodies{i}.Panels.CP.Body.z(:) = handles.PostComp.Bodies{i}.Faces.CP.Body(inds,3);
    
    %   Now we can get gradients in a "chordwise" and "spanwise" direction
    slice = str2num(get(handles.slice,'String'));
    CPmid = handles.PostComp.Bodies{i}.CPdist(slice,:);
    
    n1 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c1.Local);
    n2 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c2.Local);
    n3 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c3.Local);
    n4 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c4.Local);

    n = .25*(n1(handles.PostComp.Bodies{1}.Panels.MainPans) +...
        n2(handles.PostComp.Bodies{1}.Panels.MainPans) + ...
        n3(handles.PostComp.Bodies{1}.Panels.MainPans) + ...
        n4(handles.PostComp.Bodies{1}.Panels.MainPans));


    xmid = n(slice,:);
    cmid = (xmid - min(xmid));
    cmid = cmid/max(abs(cmid));

    fx = handles.PostComp.Bodies{i}.Fxdist(slice,:);
    fy = handles.PostComp.Bodies{i}.Fydist(slice,:);
    fz = handles.PostComp.Bodies{i}.Fzdist(slice,:);
    plot(handles.cp_axes,cmid,-CPmid,'--');
    %hold all
    %quiver(x(slice,:),z(slice,:),fx,fz)
    %clf
    
    handles.PostComp.Bodies{i}.GammaDist(slice,:) = 0;
    %surf(handles.PostComp.Bodies{i}.X(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.Y(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.Z(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.GammaDist);
    %hold all
    %scatter3(handles.PostComp.Bodies{i}.Panels.CP.Body.x(:),handles.PostComp.Bodies{i}.Panels.CP.Body.y(:),handles.PostComp.Bodies{i}.Panels.CP.Body.z(:));
end

axis(handles.cp_axes,[-.1 1 -1 4]);
axis(handles.cp_axes,'square');