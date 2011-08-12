function blade = SetBladeData(blade)
% Set current data to the selected data set.
blade.REVERSE = false;
blade.isNREL = false;
blade.TransitionPiece = [];
switch blade.type;
    case 'NREL UAE' % User selects Peaks.
        %%  NRELBlade -- NREL data
        blade.RADIUS=[0.508;0.66;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.305;5.532];
        blade.CHORD=[0.218;0.218;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.328;0.305];
        blade.THETA=[0;0;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-2.191;-2.5];
        blade.THETA = blade.THETA  - blade.THETA(end);
        blade.isNREL = true;
        blade.THICKNESS = [];
    case 'SOTON' % User selects Membrane.
       r = linspace(.2, 1, 17);
       blade.RADIUS = linspace(.08,.4, 17);
       blade.THICKNESS = [24 22.5 20.7 19.5 18.7 18.1 17.6 17.1 16.6 16.1 15.6 15.1 14.6 14.1 13.6 13.1 12.6];
        c = [.125 0.1203 0.1156 0.1109 0.1063 0.1016 0.0969 0.0922 0.0875 0.0828 0.0781 0.0734 0.0688 0.0641 0.0594 0.0547 0.05];
       blade.THETA = [15 12.1 9.5 7.6 6.1 4.9 3.9 3.1 2.4 1.9 1.5 1.2 0.9 0.6 0.4 0.2 0.0];

       blade.CHORD = c.*.4;
    case 'ESRU PoC 1' % User selects Sinc.
        blade.RADIUS=[0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1,1.050000,1.1,1.150000,1.2,1.250000;];
        blade.CHORD=[0.100,0.0965,0.0930,0.0895,0.0860,0.0825,0.0790,0.0755,0.0720,0.0685,0.0650,0.0615,0.0580,0.0545,0.0510,0.0475,0.0440,0.0405,0.0370,0.0335,0.0300;];
        blade.THETA=[35.48354,30.66771,26.60373,23.19068,20.33544,17.95276,15.96524,14.30328,12.90517,11.71700,10.69272,9.794129,8.990847,8.260345,7.587935,6.966768,6.397837,5.889976,5.459857,5.131995,4.938747;];
        blade.THICKNESS = [];
    case 'ESRU PoC 2' % User selects Sinc.
        blade.RADIUS=[0.2500;0.3000;0.3500;0.4000;0.4500;0.5000;0.5500;0.6000;0.6500;0.7000;0.7500;0.8000;0.8500;0.9000;0.9500;1;1.050;1.100;1.150;1.200;1.250;];
        blade.CHORD=0.01*[23.30;20.60;18.05;15.85;13.90;12.30;10.90;10;9.200;8.550;8;7.600;7.250;6.850;6.450;6.050;5.550;5;4.400;3.800;3.200;];
        blade.THETA= [26.40;23.25;20.50;18;15.80;13.75;12;10.55;9.550;8.700;8;7.400;6.850;6.350;5.800;5.250;4.600;4;3.400;2.800;2.200;];
        blade.THICKNESS = [];

    case 'Elliptic'
        blade.RADIUS = 0.1*linspace(-0.5,0.5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = sqrt(0.125*0.125 - linspace(-0.125,0.125));
        blade.THICKNESS = [];
    case 'Straight'
        blade.RADIUS = linspace(-2,2);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = 1*ones(size(blade.RADIUS));
        blade.THICKNESS = [];
        
        
        th = linspace(0,pi/2);
        
        x = cos(th);
        y = sin(th);
        
        x = 1-x;
        
        x = [x linspace(1,4)];
        
        y = [y ones(size(y))];
        
       % x = [x fliplr(4-x)];
       % y = [y fliplr(y)];
       
        blade.RADIUS = linspace(-2,2);
        blade.THETA =   0*ones(size(blade.RADIUS));
        blade.CHORD =  1*ones(size(blade.RADIUS));
        blade.THICKNESS = [];  
end

minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));


if blade.isNREL && ~isempty(blade.y) && (blade.Cutout.Root < 1.257) && (blade.Cutout.Root >= 0.45)
    if (blade.Cutout.Root < 0.660)
        n1 = interp1(blade.y,1:length(blade.y),0.660,'nearest');
        n2 = interp1(blade.y,1:length(blade.y),1.257,'nearest');
        
        blade.y(1:n1) = linspace(blade.Cutout.Root,0.660,n1);
        blade.y(n1:n2) = linspace(max(blade.Cutout.Root,0.660),1.257,1+n2-n1);
        
        
        blade.TransitionPiece = zeros(size(blade.y));
        blade.TransitionPiece(1:n1) = 1;
        if (n2-n1) > 1
        blade.TransitionPiece(n1:n2) = [logspace(0,-1,n2-n1) 0];
        end
    end
end

if ~isempty(blade.y)
    blade.Radius = blade.y;

    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'linear');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'linear');
    blade.Thickness = [];
    if ~isempty(blade.THICKNESS)
        blade.Thickness = interp1(blade.RADIUS,blade.THICKNESS,blade.Radius,'cubic');
    end
end

PlotBlade(blade);

function PlotBlade(blade)
hold(blade.axes,'off')
if ~isempty(blade.y)
scatter(blade.axes,blade.Radius,blade.Chord);
hold(blade.axes,'on');
end
if ~isempty(blade.THICKNESS)
    plot(blade.axes,blade.RADIUS,blade.THICKNESS/100);
    hold(blade.axes,'on');
end
plotyy(blade.axes,blade.RADIUS,blade.CHORD,blade.RADIUS,blade.THETA);

