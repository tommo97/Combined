function DistPanel = PanelDistButtonsParam(DistPanel)
minx = DistPanel.minx;
maxx = DistPanel.maxx;
lin = get(DistPanel.lin_button,'Value');
bel = get(DistPanel.bell_button,'Value');
num = str2num(get(DistPanel.NumPanels,'String'));
par = str2num(get(DistPanel.bell_param,'String'));

%   Check which are on, if none, use default values
set(DistPanel.bell_param,'enable','off');
if lin
    set(DistPanel.bell_param,'enable','off');
    DistPanel.x = linspace(minx,maxx,num+1);
elseif bel
    set(DistPanel.bell_param,'enable','on');
    DistPanel.x = BellShape(minx,maxx,num+1,par);
else
    disp('AAArgh');
    DistPanel.x = linspace(minx,maxx,num+1);
end