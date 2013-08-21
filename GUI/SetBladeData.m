function blade = SetBladeData(blade)
% Set current data to the selected data set.
blade.REVERSE = false;
blade.isNREL = false;   %   for transisition piece
blade.isSOTON = false;  %   for transisition piece
blade.isBarltrop = false;
blade.TransitionPiece = [];
blade.SWEEP = [];
switch blade.type;
    case 'NREL UAE' % User selects Peaks.
        %%  NRELBlade -- NREL data
        blade.RADIUS=[0.508;0.66;0.883;1.008;1.067;1.133;1.257;1.343;1.51;1.648;1.9520;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.305;5.532];
        blade.CHORD=[0.218;0.218;0.183;0.349;0.441;0.544;0.737;0.728;0.711;0.697;0.6660;0.636;0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.328;0.305];
               
        blade.THETA=[0;0;0;6.7;9.9;13.4;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-2.191;-2.5];
        %blade.THETA=[20.04;20.04;20.04;20.04;20.04;20.04;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
        %    3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-2.191;-2.5];
          % Set up chord fit and options.
        ChordFit = FitData(blade.RADIUS, blade.CHORD);
        TwistFit = FitData(blade.RADIUS, blade.THETA);
        
        % Fit model to data.
        blade.RADIUS = linspace(min(blade.RADIUS),max(blade.RADIUS),10000);
        blade.CHORD = ChordFit(blade.RADIUS);
        blade.THETA = TwistFit(blade.RADIUS);
        blade.isNREL = true;
        blade.THICKNESS = [];
    case 'SOTON' % User selects Membrane.
        
        bahaj = [20 0.03 15 24
            30 0.03 15 24
            80.0000    0.1250   15.0000   24.0000
            120.0000    0.1156    9.5000   20.7000
            160.0000    0.1063    6.1000   18.7000
            200.0000    0.0969    3.9000   17.6000
            240.0000    0.0875    2.4000   16.6000
            280.0000    0.0781    1.5000   15.6000
            320.0000    0.0688    0.9000   14.6000
            360.0000    0.0594    0.4000   13.6000
            400.0000    0.0500         0   12.6000];
        
        blade.RADIUS = bahaj(:,1)/1000;
        blade.CHORD = bahaj(:,2)*0.4;
        blade.THETA = bahaj(:,3);
        blade.THICKNESS = bahaj(:,4);
        
        
        
        % Set up chord fit and options.
        ChordFit = FitData(blade.RADIUS, blade.CHORD);
        TwistFit = FitData(blade.RADIUS, blade.THETA);
        ThickFit = FitData(blade.RADIUS, blade.THICKNESS);
        
       
        % Fit model to data.
        blade.RADIUS = linspace(min(blade.RADIUS),max(blade.RADIUS),10000);
        blade.CHORD = ChordFit(blade.RADIUS);
        blade.THETA = TwistFit(blade.RADIUS);
        blade.THICKNESS = ThickFit(blade.RADIUS);
        blade.isSOTON = true;
 
        
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
        
        ChordData = [ 0.05 10
            0.075 10
            0.125 20
            0.15	43.910256
            0.200913	42.948718
            0.248858	41.346154
            0.301370	40.064103
            0.351598	38.461538
            0.401826	36.858974
            0.447489	35.576923
            0.500000	34.134615
            0.598174	31.410256
            0.698630	28.525641
            0.801370	25.801282
            0.904110	22.916667
            1.000000	20.192308];
        TwistData = [0.025 32.211
            0.05 32.211
            0.15	32.211538
            0.198630	27.083333
            0.248858	20.032051
            0.301370	15.224359
            0.349315	11.057692
            0.399543	7.532051
            0.449772	5.608974
            0.502283	4.166667
            0.600457	2.243590
            0.698630	0.641026
            0.799087	-0.320513
            0.899543	-1.282051
            1.000	-1.762821];

        ChordFit = FitData(ChordData(:,1), ChordData(:,2));
        TwistFit = FitData(TwistData(:,1), TwistData(:,2));
        
        blade.isBarltrop = true;
        % Fit model to data.
        blade.RADIUS = linspace(0.05,1.0,10000);
        blade.THETA = TwistFit(blade.RADIUS);%interp1(TwistData(:,1),TwistData(:,2),blade.RADIUS);
        blade.CHORD = ChordFit(blade.RADIUS)/1000;%interp1(ChordData(:,1),ChordData(:,2),blade.RADIUS)/1000;
        blade.RADIUS = 0.175*linspace(0.05,1.0,10000);
        blade.THICKNESS = [];
    case 'Straight'
        
        
        US = [0.05 0.018761
            0.251040	0.018761
            0.814147	0.018761
            0.843273	0.020177
            0.889043	0.034336
            0.908460	0.035752
            0.947295	0.020177
            0.968100	0.007434
            1.000000	-0.023717];
        LS = [
            1.000000	-0.097345
            0.990291	-0.087434
            0.911234	-0.053451
            0.251040	-0.053451
            0.05         -0.025];
        rr = linspace(0.05,1);
        uz = interp1(US(:,1),US(:,2),rr,'cubic');
        lz = interp1(LS(:,1),LS(:,2),rr,'cubic');
        
        swp = 0.5*(uz + lz);
%         close all
%         plot(rr,lz)
%         hold all
%         plot(rr,uz)
%         plot(rr,swp);
%         axis equal tight
        


        blade.RADIUS = linspace(-20,20);
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
        
        blade.RADIUS = rr*7.5;%linspace(-10,10);
        blade.THETA =   0*ones(size(blade.RADIUS));
        blade.CHORD =  7.5*(uz-lz);%1*ones(size(blade.RADIUS));
        blade.SWEEP =  -7.5*swp;
        blade.THICKNESS = [];
        
    case 'Wing'
        blade.RADIUS = linspace(-10,10)/5;
        blade.THETA =   [linspace(5,0,50) linspace(0,5,50)];
        
        theta = [4.4  4.3 3.8 3.5 3 2.9 2.5 2 1.9 1.5 0.6 0.1 -0.6 -2 -5.5];
        
        r = linspace(0.1, 0.975,numel(theta))*max(blade.RADIUS);
        theta = -[fliplr(theta) theta];
        r = [fliplr(-r) r];

        blade.THETA = interp1(r,theta,blade.RADIUS,'linear','extrap');
        r = blade.RADIUS;
        rmin = r(20);
        rmax = r(40);
        blade.THETA(boolean(1-((abs(r)>=abs(rmin)) + (abs(r)<=abs(rmax))))) = blade.THETA(boolean(1-((abs(r)>=abs(rmin)) + (abs(r)<=abs(rmax)))))-10;
        
        blade.CHORD =  [linspace(3.983,17.7,50)/17.7 linspace(17.7,3.983,50)/17.7];
        blade.SWEEP =  sind(33.5)*[linspace(10,0,50) linspace(0,10,50)]/5;
end


minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));


% if blade.isNREL || blade.isSOTON || blade.isBarltrop
%     bcR = (blade.RADIUS(find(blade.CHORD==max(blade.CHORD))));
%     if blade.isSOTON
%         shoulder = 0.315;
%     elseif blade.isNREL
%         shoulder = 0.66;
%     else
%         shoulder = 0.024;
%     end
%     if ~isempty(blade.y) && (blade.Cutout.Root < bcR) && (blade.Cutout.Root >= min(blade.RADIUS))
%         if (blade.Cutout.Root <= shoulder)
%             blade.n1 = interp1(blade.y,1:length(blade.y),shoulder,'nearest');
%             blade.n2 = interp1(blade.y,1:length(blade.y),blade.RADIUS(find(blade.CHORD==max(blade.CHORD))),'nearest');
%             disp(blade.RADIUS(find(blade.CHORD==max(blade.CHORD))))
%             blade.y(1:blade.n1) = linspace(blade.Cutout.Root,shoulder,blade.n1);
%             blade.y(blade.n1:blade.n2) = linspace(max(blade.Cutout.Root,shoulder),bcR,1+blade.n2-blade.n1);
%             
%             
%             blade.TransitionPiece = zeros(size(blade.y));
%             blade.TransitionPiece(1:blade.n1) = 1;
%             if (blade.n2-blade.n1) > 1
%                 blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
%             end
%         end
%     end
% end

if blade.isNREL || blade.isSOTON || blade.isBarltrop
    bcR = (blade.RADIUS(find(blade.CHORD==max(blade.CHORD))));
    if blade.isSOTON
        shoulder = 0.0315;
    elseif blade.isNREL
        shoulder = 0.66;
    else
        shoulder = 0.014;
    end
    if ~isempty(blade.y) && (blade.Cutout.Root < bcR) && (blade.Cutout.Root >= min(blade.RADIUS))
        
        if (blade.Cutout.Root < bcR)
            blade.n1 = interp1(blade.y,1:length(blade.y),shoulder,'nearest');
            blade.n2 = interp1(blade.y,1:length(blade.y),blade.RADIUS(find(blade.CHORD==max(blade.CHORD))),'nearest')+1;
            disp(num2str([blade.n1 blade.n2]));
            disp(blade.RADIUS(find(blade.CHORD==max(blade.CHORD))))
            %blade.y(1:blade.n1) = linspace(blade.Cutout.Root,0.05,blade.n1);
            blade.y(1:blade.n2) = linspace(max(blade.Cutout.Root,min(blade.RADIUS)),bcR,blade.n2);
            
            blade.y(blade.n2:end) = linspace(bcR,max(blade.RADIUS),numel(blade.y(blade.n2:end)));
            blade.TransitionPiece = zeros(size(blade.y));
            blade.TransitionPiece(1:blade.n2) = 1;
            blade.CircTrans = 1:blade.n2;
            %     if (blade.n2-blade.n1) > 1
            %     blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
            %     end
        end
    end
end
    
% if blade.isSOTON 
%     bcR = (blade.RADIUS(find(blade.CHORD==max(blade.CHORD))));
%     shoulder = 0.315;
% 
%     if ~isempty(blade.y) && (blade.Cutout.Root < bcR) && (blade.Cutout.Root >= min(blade.RADIUS))
%         
%         if (blade.Cutout.Root < bcR)
%             blade.n1 = interp1(blade.y,1:length(blade.y),0.0315,'nearest');
%             blade.n2 = interp1(blade.y,1:length(blade.y),blade.RADIUS(find(blade.CHORD==max(blade.CHORD))),'nearest');
%             disp(blade.RADIUS(find(blade.CHORD==max(blade.CHORD))))
%             %blade.y(1:blade.n1) = linspace(blade.Cutout.Root,0.05,blade.n1);
%             blade.y(1:blade.n2+1) = linspace(max(blade.Cutout.Root,0.02),bcR,blade.n2+1);
%             
%             blade.y(blade.n2+1:end) = linspace(bcR,0.4,numel(blade.y(blade.n2+1:end)));
%             blade.TransitionPiece = zeros(size(blade.y));
%             blade.TransitionPiece(1:blade.n2) = 1;
%             blade.CircTrans = 1:blade.n2;
%             %     if (blade.n2-blade.n1) > 1
%             %     blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
%             %     end
%         end
%     end
% end
% 
% if blade.isBarltrop
%     bcR = (blade.RADIUS(find(blade.CHORD==max(blade.CHORD))));
%     if ~isempty(blade.y) && (blade.Cutout.Root < bcR) && (blade.Cutout.Root >= 0.0175)
%         if (blade.Cutout.Root < bcR)
%             blade.n1 = interp1(blade.y,1:length(blade.y),0.024,'nearest');
%             blade.n2 = interp1(blade.y,1:length(blade.y),blade.RADIUS(find(blade.CHORD==max(blade.CHORD))),'nearest');
%             disp(blade.RADIUS(find(blade.CHORD==max(blade.CHORD))))
%             
%             %blade.y(1:blade.n1) = linspace(blade.Cutout.Root,0.05,blade.n1);
%             blade.y(1:blade.n2+1) = linspace(max(blade.Cutout.Root,0.0175),bcR,blade.n2+1);
%             
%             blade.y(blade.n2+1:end) = linspace(bcR,0.4,numel(blade.y(blade.n2+1:end)));
%             blade.TransitionPiece = zeros(size(blade.y));
%             blade.TransitionPiece(1:blade.n2) = 1;
%             blade.CircTrans = 1:blade.n2;
%             %     if (blade.n2-blade.n1) > 1
%             %     blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
%             %     end
%         end
%     end
% end


if ~isempty(blade.y)
    blade.Radius = blade.y;
    
    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'linear');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'linear');
    blade.Thickness = [];
    if ~isempty(blade.THICKNESS)
        blade.Thickness = interp1(blade.RADIUS,blade.THICKNESS,blade.Radius,'cubic');
    end
    if ~isempty(blade.SWEEP)
        blade.Sweep = interp1(blade.RADIUS,blade.SWEEP,blade.Radius,'cubic');
    else
            blade.Sweep = zeros(size(blade.Radius));
    end

end

PlotBlade(blade);


function Fit = FitData(xdata,ydata)
 % Set up chord fit and options.
        [xData, yData] = prepareCurveData(xdata, ydata);
        ft = fittype( 'pchipinterp' );
        opts = fitoptions( ft );
        opts.Normalize = 'on';
        Fit = fit( xData, yData, ft, opts );


function PlotBlade(blade)
hold(blade.axes,'off')
if ~isempty(blade.y)
    scatter(blade.axes,blade.Radius,blade.Chord);
    hold(blade.axes,'on');
    if (exist('blade.n1') && exist('blade.n2'))
    plot(blade.axes,blade.Radius(blade.n2),blade.Chord(blade.n2),'o','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','none',...
        'MarkerSize',10)
    
    plot(blade.axes,blade.Radius(blade.n1),blade.Chord(blade.n1),'o','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','none',...
        'MarkerSize',10)
    end
end
if ~isempty(blade.THICKNESS)
    plot(blade.axes,blade.RADIUS,blade.THICKNESS/100);
    
   
    hold(blade.axes,'on');
end
plotyy(blade.axes,blade.RADIUS,blade.CHORD,blade.RADIUS,blade.THETA);

set(blade.axes,'XGrid','on','YGrid','on','box','on')