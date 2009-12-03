function handles = PostPlot(handles,ax)
hold(ax,'off')
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
    
    
    n1 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c1.Local);
    n2 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c2.Local);
    n3 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c3.Local);
    n4 = handles.PostComp.Bodies{i}.n(handles.PostComp.Bodies{i}.Panels.c4.Local);
    
    n = .25*(n1(handles.PostComp.Bodies{1}.Panels.MainPans) +...
        n2(handles.PostComp.Bodies{1}.Panels.MainPans) + ...
        n3(handles.PostComp.Bodies{1}.Panels.MainPans) + ...
        n4(handles.PostComp.Bodies{1}.Panels.MainPans));
    
    if ~handles.PostComp.use_slice
        CPress = interp1(handles.PostComp.Bodies{i}.Panels.CP.Body.y(:,1),handles.PostComp.Bodies{i}.CPdist,handles.PostComp.r_slice,'cubic');
        xPress = interp1(handles.PostComp.Bodies{i}.Panels.CP.Body.y(:,1),n,handles.PostComp.r_slice,'cubic');
    else
        slice = str2num(get(handles.slice,'String'));
        CPress = handles.PostComp.Bodies{i}.CPdist(slice,:);
        xPress = n(slice,:);
        fx = handles.PostComp.Bodies{i}.Fxdist(slice,:);
        fy = handles.PostComp.Bodies{i}.Fydist(slice,:);
        fz = handles.PostComp.Bodies{i}.Fzdist(slice,:);
    end
    %hold all
    %quiver(x(slice,:),z(slice,:),fx,fz)
    %clf
    PressChord = (xPress - min(xPress));
    PressChord = PressChord/max(abs(PressChord));
    PressChord(end+1) = PressChord(1);
    CPress(end+1) = CPress(1);
    plot(ax,PressChord,-CPress);
   
    %surf(handles.PostComp.Bodies{i}.X(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.Y(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.Z(handles.PostComp.Bodies{i}.N.Local),handles.PostComp.Bodies{i}.GammaDist);
    %hold all
    %scatter3(handles.PostComp.Bodies{i}.Panels.CP.Body.x(:),handles.PostComp.Bodies{i}.Panels.CP.Body.y(:),handles.PostComp.Bodies{i}.Panels.CP.Body.z(:));
end

if handles.PostComp.r_R==.30
hold(ax,'on')
data = [0.00000E+0	2.88412
5.19673E-3	2.19367
1.11359E-2	1.85086
1.93022E-2	1.87500
4.00891E-2	1.68670
6.01336E-2	1.48391
9.94803E-2	1.26180
2.00445E-1	1.17006
3.60059E-1	0.99624
5.59020E-1	0.48927
6.80030E-1	0.18991
8.01782E-1	0.02575
1.00074E+0	-0.14807
9.20564E-1	-0.27843
6.80030E-1	-0.15773
4.40980E-1	0.24785
2.80624E-1	0.02575
1.40312E-1	-0.49571
5.93912E-2	-0.84818
1.93022E-2	-1.00751
8.90869E-3	-0.90612
4.45434E-3	-0.54882];
scatter(ax,data(:,1),data(:,2));
elseif handles.PostComp.r_R==.47
hold(ax,'on')
data = [9.65108E-3	2.04677
2.15293E-2	1.96452
4.00891E-2	1.76613
5.93912E-2	1.54355
1.00223E-1	1.39355
1.99703E-1	1.25323
3.60059E-1	1.09839
5.59762E-1	0.61935
6.80030E-1	0.26613
8.01039E-1	0.07258
1.00074E+0	-0.16452
9.20564E-1	-0.34355
6.80030E-1	-0.16935
4.40238E-1	0.29032
2.79881E-1	0.10161
1.39569E-1	-0.39677
6.01336E-2	-0.78871
2.00445E-2	-1.00645
9.65108E-3	-0.95806
4.45434E-3	-0.66290
4.45434E-3	1.88226
-7.42391E-4	1.96935];
scatter(ax,data(:,1),data(:,2));
elseif handles.PostComp.r_R==.63
hold(ax,'on')
data = [0.00776	2.04706
0.02039	1.96486
0.03897	1.76656
0.05979	1.54408
0.09989	1.39912
0.19939	1.25443
0.36051	1.09070
0.55952	0.62231
0.67983	0.28418
0.80012	0.08153
0.99985	-0.15944
0.91969	-0.33886
0.67988	-0.17550
0.43855	0.29753
0.27970	0.09837
0.14093	-0.40553
0.05931	-0.78818
0.01999	-1.01095
0.00959	-0.95294
0.00510	-0.66264
0.00407	1.88736
-0.00114	1.96959];
scatter(ax,data(:,1),data(:,2));
elseif handles.PostComp.r_R==.80
hold(ax,'on')
data = [9.65108E-3	1.76265
1.93022E-2	1.62702
3.93467E-2	1.41227
6.01336E-2	1.31432
9.94803E-2	1.26911
2.00445E-1	1.16738
3.59317E-1	1.10710
5.59020E-1	0.61733
6.79287E-1	0.27072
8.00297E-1	0.07858
9.99258E-1	-0.20398
9.19822E-1	-0.34338
6.79287E-1	-0.11733
4.40238E-1	0.36114
2.80624E-1	0.17277
1.39569E-1	-0.30947
5.86489E-2	-0.72766
2.30141E-2	-0.98385
8.16630E-3	-1.00269
2.96956E-3	-0.80678
0.00000E+0	1.49892
5.19673E-3	1.59688];
scatter(ax,data(:,1),data(:,2));
elseif handles.PostComp.r_R==.95
hold(ax,'on')
data = [2.04987E-2	1.15498
3.97892E-2	1.08725
5.98200E-2	1.00663
9.91584E-2	0.95817
2.00111E-1	0.88695
3.60463E-1	0.86066
5.60088E-1	0.40244
6.79577E-1	0.19260
8.00560E-1	0.03754
9.99477E-1	-0.21767
9.20026E-1	-0.30764
6.80269E-1	-0.10063
4.39836E-1	0.48983
2.80198E-1	0.35178
1.39051E-1	-0.19234
6.02907E-2	-0.58199
1.94041E-2	-0.90409
8.99356E-3	-1.00395
3.06932E-3	-0.91693
-4.01784E-4	0.49125
4.83849E-3	0.74580
1.08250E-2	1.02290];
scatter(ax,data(:,1),data(:,2));
end

axis(ax,[-.1 1.1 -1 4]);
axis(ax,'square');