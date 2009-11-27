function blade = SetBladeData(blade)
% Set current data to the selected data set.

switch blade.type;
    case 'NREL UAE' % User selects Peaks.
        %%  NRELBlade -- NREL data
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];
    case 'SOTON' % User selects Membrane.
        blade.RADIUS = [0.0800    0.1200    0.1600    0.2000    0.2400    0.2800    0.3200    0.3600    0.4000];
        blade.CHORD = [0.0500    0.0462    0.0425    0.0388    0.0350    0.0312    0.0275    0.0238    0.0200];
        blade.THETA = [15.0000    9.5000    6.1000    3.9000    2.4000    1.5000    0.9000    0.4000         0];
    case 'ESRU PoC 1' % User selects Sinc.
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];
    case 'ESRU PoC 2' % User selects Sinc.
        blade.RADIUS=[0.45;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.029;5.25];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.666;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.355;0.355];
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-1.815;-1.815];
    case 'Elliptic'
        blade.RADIUS = linspace(-5,5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = sqrt(sin(linspace(0,pi))) + .2;
    case 'Straight'
        blade.RADIUS = linspace(-5,5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = 1*ones(size(blade.RADIUS));
end

minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));

if ~isempty(blade.y)
    blade.Radius = blade.y;
    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'cubic');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'cubic');
end

PlotBlade(blade);

function PlotBlade(blade)
hold(blade.axes,'off')
if ~isempty(blade.y)
scatter(blade.axes,blade.Radius,blade.Chord);
hold(blade.axes,'on');
end
plotyy(blade.axes,blade.RADIUS,blade.CHORD,blade.RADIUS,blade.THETA);

