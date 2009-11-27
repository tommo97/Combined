function TabSelectionFcn(hFig,tag,count,action)
% TABSELECTIONFCN  allows to switch and enable/disable the tabpanels programatically
% Using the TABSELECTIONFCN it is possible to switch, enable and disable
% the tabs (that were created using the TabPanel Constructor) programatically.
%
% usage:
% Select Tabpanel
% ---------------
% TABSELECTIONFCN(<hFig>,<TabTag>,<tabnumber>)
%     <hFig>      the handle of the Figure
%     <TabTag>    the Tag name of the tabpanel
%     <tabnumber> The number of the tabpanel or the tab string
%
% Enable/Disable Tabpanel(s)
% --------------------------
% TABSELECTIONFCN(<hFig>,<TabTag>,<tabnumber(s)>,'on/off')
%     <tabnumber(s)> can be a scalar or vector of indices.
%
% Note: You can not disable the tabpanel that is currently active.
%
% See also TABPANEL.
%
%   Version: v1.1
%      Date: 2008/09/02 00:00:00
%   (c) 2008 By Elmar Tarajan [MCommander@gmx.de]