function blade = SetBladeData(blade)
% Set current data to the selected data set.
blade.MakeCaps = false;
switch blade.type;
    case 'NREL UAE' % User selects Peaks.
        %%  NRELBlade -- NREL data
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.305;5.532];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.328;0.305];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-2.191;-2.5];
        blade.THETA = blade.THETA  - blade.THETA(end);
        blade.THICKNESS = [];
    case 'SOTON' % User selects Membrane.
       r = linspace(.2, 1, 17);
       blade.RADIUS = linspace(.08,.4, 17);
       blade.THICKNESS = [24 22.5 20.7 19.5 18.7 18.1 17.6 17.1 16.6 16.1 15.6 15.1 14.6 14.1 13.6 13.1 12.6];
        c = [.125 0.1203 0.1156 0.1109 0.1063 0.1016 0.0969 0.0922 0.0875 0.0828 0.0781 0.0734 0.0688 0.0641 0.0594 0.0547 0.05];
       blade.THETA = [15 12.1 9.5 7.6 6.1 4.9 3.9 3.1 2.4 1.9 1.5 1.2 0.9 0.6 0.4 0.2 0.0];

       blade.CHORD = c.*.4;
    case 'ESRU PoC 1' % User selects Sinc.
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];
        blade.THICKNESS = [];
    case 'ESRU PoC 2' % User selects Sinc.
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];
        blade.THICKNESS = [];
    case 'Elliptic'
        blade.RADIUS = linspace(-5,5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = sqrt(sin(linspace(0,pi))) + .2;
        blade.THICKNESS = [];
    case 'Straight'
        blade.RADIUS = linspace(-.5,.5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = 1*ones(size(blade.RADIUS));
        blade.THICKNESS = [];
        blade.MakeCaps = false;
        
end

minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));

if ~isempty(blade.y)
    blade.Radius = blade.y;
    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'cubic');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'cubic');
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

