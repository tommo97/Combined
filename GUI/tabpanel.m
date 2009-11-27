function tabpanel(figname,tag,action)
%TABPANEL  "TabPanel Constructor" offers the easiest way for creating tabpanels in MATLAB
%  Usage:
%  1. Open the figure (where the tab panel needs to be integrated) in
%     GUIDE and create a Text-Object which reserves the place of the future
%     tabpanel.
%  2. Specify and note the unique tag-name of the created text-object and the
%     figure filename.
%  3. Start "TabPanel Constructor" as follows:
%        >> tabpanel('filename.fig','tabpaneltag')
%
%  Options:
%     a. activate "TabPanel Constructor" to edit an existing tabpanel:
%        >> tabpanel('filename.fig','tabpaneltag') 
%     b. remove tabpanel from GUI
%        >> tabpanel('filename.fig','tabpaneltag','delete')
%
%  !!! IMPORTANT !!!
%  The current version fixes the bug (since R2009a) which force the GUI to 
%  restart when the M-File does not contain the TabSelectionFcn.
%  To fix the issue for already existing tabpanel please follow the steps:
%     1. Open the tabpanel created with TPC 2.6 with the current TPC 2.6.1
%     2. Add any new panel
%     3. Remove currently created panel
%     4. Close TPC
%
% See also TABSELECTIONFCN.
%
%   Version: v2.6.1
%      Date: 2009/07/24 00:00:00
%   (c) 2005-2009 By Elmar Tarajan [MCommander@gmx.de]

%   2009/07/24 fixed: Error while closing TPC in R2009a/b
%   2009/07/24 fixed: GUI-Restart while tab switching
%   2009/07/24 some code improvements
%
%   2008/09/10 "tabselectionfcn" for programatically tab switching added
%   2008/08/03 TabSelectionChange_Callback added
%   2008/07/02 supporting of nested tabpanels
%   2008/06/23 works now with R14/SP1/SP2/SP3 versions
%   2008/05/08 many code improvements
%   2008/04/16 improved look - using the mouse on "settings"-button
%   2008/03/28 some code improvements
%   2008/03/17 improved look of tabs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I LOVE MATLAB! You too? :) %%%