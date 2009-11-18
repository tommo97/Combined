function [foil data] = SetFoilData(foil,DistPanel)
type = foil.section;
[Foil Data] = MakeFoil(DistPanel.x,foil.thickness);
set(foil.thickness_slider,'enable','off');
set(foil.thickness_edit_text,'enable','off');
% Set current data to the selected data set.
switch type;
    case 'NREL S809' % User selects Peaks.
        foil.US = Foil.S809.US;
        foil.LS = Foil.S809.LS;
        foil.X = Foil.S809.X;
        data.US = Data.S809.US;
        data.LS = Data.S809.LS;
        data.X = Data.S809.X;
        
    case 'NREL S814' % User selects Membrane.
        foil.US = Foil.S814.US;
        foil.LS = Foil.S814.LS;
        foil.X = Foil.S814.X;
        data.US = Data.S814.US;
        data.LS = Data.S814.LS;
        data.X = Data.S814.X;
    case 'NACA 00xx' % User selects Membrane.
        set(foil.thickness_slider,'enable','on');
        set(foil.thickness_edit_text,'enable','on');
        foil.US = Foil.N00xx.US;
        foil.LS = Foil.N00xx.LS;
        foil.X = Foil.N00xx.X;
        data.US = Data.N00xx.US;
        data.LS = Data.N00xx.LS;
        data.X = Data.N00xx.X;
    otherwise
        foil.US = Foil.N0012.US;
        foil.LS = Foil.N0012.LS;
        foil.X = Foil.N0012.X;
        data.US = Data.N0012.US;
        data.LS = Data.N0012.LS;
        data.X = Data.N0012.X;
end

foil.SectionShape = data;


PlotFoil(foil.X,foil.US,foil.LS,foil.axes);
if ~isempty(DistPanel.x)
    ScatterFoil(data.X,data.US,data.LS,foil.axes);
end


function PlotFoil(ix,ius,ils,iaxes)
hold(iaxes,'off')
plot(iaxes,ix,ius);
hold(iaxes,'on');
plot(iaxes,ix,ils);
axis(iaxes,'equal')
hold(iaxes,'off')

function ScatterFoil(ix,ius,ils,iaxes)
hold(iaxes,'on');
scatter(iaxes,ix,ius);
scatter(iaxes,ix,ils);
axis(iaxes,'equal')
hold(iaxes,'off')