function blade = SetBladeData(blade)
% Set current data to the selected data set.
blade.REVERSE = false;
blade.isNREL = false;   %   for transisition piece
blade.isSOTON = false;  %   for transisition piece
blade.isBarltrop = false;
blade.isProp = false;
blade.TransitionPiece = [];

disp(blade.type);
blade.LoadFile = false;

if (~strcmp('Load From File',blade.type))
    set(blade.LoadFileBox,'enable','off','string','e.g. data.mat');
end

    set(blade.ScaleFactorNotice,'string',['Suggested Max Scale N/A'],'visible','off');


blade.SKEW = [];
blade.CHORD = [];
blade.RAKE = [];
blade.THETA = [];
blade.THICKNESS = [];
blade.CAMBER = [];
switch blade.type;
    case 'Load From File'
        
        blade.LoadFile = true;
        set(blade.LoadFileBox,'enable','on');
        fname = get(blade.LoadFileBox,'string');
        set(blade.FileNotFound,'String','File not found','visible','off')

        try
            error = false;
            try
            BladeSpec = load(fname);
            catch
                error = true;
            end
            try
                blade.RADIUS = Radius;
            catch
                error = true;
            end
            try
                blade.CHORD = Chord;
            catch
                error = true;
            end
            try
                blade.THETA = Twist;
            catch
                error = true;
            end
            try
                blade.THICKNESS = Thickness;
            end
            try
                blade.SKEW = Skew;
            end
            try
                blade.RAKE = Rake;
            end
            if (~error)
                set(blade.FileNotFound,'visible','off')
            else
                set(blade.FileNotFound,'string','Error in file','visible','on');
            end
        catch
            disp(['File ' fname ' not found']);
            set(blade.FileNotFound,'String','File not found','visible','on')
            
        end
        
    case 'Marine Prop'
        blade.isProp = true;
        blade = BladeFromBladeSpec('../GeomInputMatFiles/DMTP4119Approximation.mat', blade);
        blade.THICKNESS(end) = 0;
    case 'NREL Phase VI' % User selects Peaks.
        blade = BladeFromBladeSpec('../GeomInputMatFiles/NRELPhaseVI.mat', blade);
        blade.isNREL = true;
    case 'Southampton Rotor' % User selects Membrane.
        blade = BladeFromBladeSpec('../GeomInputMatFiles/SOTON.mat', blade);
        blade.isSOTON = true;
    case 'ESRU PoC 1 ''05' % User selects Sinc.
        blade = BladeFromBladeSpec('../GeomInputMatFiles/ESRUPoc1.mat', blade);
    case 'ESRU PoC 2 ''05' % User selects Sinc.
        blade = BladeFromBladeSpec('../GeomInputMatFiles/ESRUPoc2.mat', blade);
    case 'Barltrop ''05'
        blade = BladeFromBladeSpec('../GeomInputMatFiles/Barltrop05.mat', blade);
        blade.isBarltrop = true;
    case 'Elliptic Wing'
        blade.RADIUS = linspace(-5,5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = sqrt(sin(linspace(0,pi))) + .2;
    case 'Straight Wing'    
        blade.RADIUS = linspace(-20,20);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = 1*ones(size(blade.RADIUS));     
    case 'BERP Helicopter'
        blade = BladeFromBladeSpec('../GeomInputMatFiles/BERPapprox.mat', blade);
    case 'A380 Type Planform'
       blade = BladeFromBladeSpec('../GeomInputMatFiles/A380approx.mat', blade);
end


minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));

disp(blade);


if ~isempty(blade.y)
    

    blade.Radius = blade.y;
    blade.SuggestedScaleFactor = 1/max(diff(blade.y));
    set(blade.ScaleFactorNotice,'string',['Suggested Max Scale ' num2str(blade.SuggestedScaleFactor)],'visible','on');
    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'cubic');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'cubic');
    blade.Thickness = [];
    
    
    if (blade.isProp)
        s = cumsum([0 sqrt((diff(blade.Chord).^2) + (diff(blade.Radius).^2))]);
        s = s./max(s);
        sdesired = (blade.y-min(blade.y))/(max(blade.y) - min(blade.y));
        blade.Radius = interp1(s,blade.y,sdesired,'cubic');
        blade.SuggestedScaleFactor = 1/max(diff(blade.y));
        set(blade.ScaleFactorNotice,'string',['Suggested Max Scale ' num2str(blade.SuggestedScaleFactor)],'visible','on');
        blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'cubic');
        blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'cubic');
        blade.Thickness = [];
    end
    
        
    if ~isempty(blade.CAMBER)
        blade.Camber = interp1(blade.RADIUS,blade.CAMBER,blade.Radius,'cubic');
    else
        blade.Camber = [];
    end
    
    if ~isempty(blade.THICKNESS)
        blade.Thickness = interp1(blade.RADIUS,blade.THICKNESS,blade.Radius,'cubic');
    end
    
    if ~isempty(blade.SKEW)
        blade.Skew = interp1(blade.RADIUS,blade.SKEW,blade.Radius,'cubic');
    else
        blade.Skew = zeros(size(blade.Radius));
    end
    
    if ~isempty(blade.RAKE)
        blade.Rake = interp1(blade.RADIUS,blade.RAKE,blade.Radius,'cubic');
    else
        blade.Rake = zeros(size(blade.Radius));
    end
else
        set(blade.ScaleFactorNotice,'visible','off');
end


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




PlotBlade(blade);
%disp(blade);

function blade = BladeFromBladeSpec(infname, blade)
load(infname);

blade.RADIUS = linspace(min(Radius),max(Radius),10000);
blade.Thickness = [];
blade.SKEW = zeros(size(blade.RADIUS));
blade.RAKE = blade.SKEW;
blade.CAMBER = blade.SKEW;

ChordFit = FitData(Radius, Chord);
TwistFit = FitData(Radius, Twist);
blade.CHORD = ChordFit(blade.RADIUS);
blade.THETA = TwistFit(blade.RADIUS);
if ~isempty(Thickness)
    ThickFit = FitData(Radius, Thickness);
    blade.THICKNESS = ThickFit(blade.RADIUS);
end
if ~isempty(Skew)
    SkewFit  = FitData(Radius, Skew);
    blade.SKEW = SkewFit(blade.RADIUS);
end
if ~isempty(Rake)
    RakeFit  = FitData(Radius, Rake);
    blade.RAKE = RakeFit(blade.RADIUS);
end

try 
    C = Camber;
    CamberFit  = FitData(Radius, C);
    blade.CAMBER = CamberFit(blade.RADIUS);
catch
    blade.CAMBER = [];
end



        
        
function Fit = FitData(xdata,ydata)
% Set up chord fit and options.
[xData, yData] = prepareCurveData(xdata, ydata);
ft = fittype( 'pchipinterp' );
opts = fitoptions( ft );
opts.Normalize = 'on';
Fit = fit( xData, yData, ft, opts );


function PlotBlade(blade)
cla(blade.axes,'reset');
hold(blade.axes,'off')

if ~isempty(blade.y)
    plotyy(blade.axes,blade.RADIUS,blade.CHORD,blade.RADIUS,blade.THETA);
    hold(blade.axes,'on');
    scatter(blade.axes,blade.Radius,blade.Chord);
    
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


set(blade.axes,'XGrid','on','YGrid','on','box','on')