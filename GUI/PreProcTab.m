function varargout = PreProcTab(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreProcTab_OpeningFcn, ...
                   'gui_OutputFcn',  @PreProcTab_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PreProcTab is made visible.
function PreProcTab_OpeningFcn(source, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% source    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreProcTab (see VARARGIN)

% Choose default command line output for PreProcTab
handles.output = source;


%%  Set up Simulation panel defaults
set(handles.minmodbeta,'enable','off');
%%  Set up Output panel defaults
set(handles.output_case,'Value',true);
set(handles.output_neu,'Value',true);
set(handles.output_runfile,'Value',true);
set(handles.commit_to_head,'enable','off');
set(handles.output_svn,'enable','off');
set(handles.timestamp,'Value',false);

handles.PostComp.r_R = .47;

handles.PostComp.use_slice = true;
handles.Origin = [0;0;0];
handles.Attitude = [0;0;0];
%%  Clear the axes
cla(handles.root_foil_axes,'reset');
cla(handles.tip_foil_axes,'reset');
cla(handles.goto_axes,'reset');
cla(handles.cp_axes,'reset');
cla(handles.full_turbine_axes,'reset');
cla(handles.blade_surf_axes,'reset');
cla(handles.blade_geom_axes,'reset');
%%  Set up tabs
% hg=uitabgroup;
% ht(1)=uitab(hg,'Title','Tab 1');
% ht(2)=uitab(hg,'Title','Tab 2');
% 
% axes('parent',ht(1));
% plot(1:10,(1:10).^2);
% 
% axes('parent',ht(2));
% plot(1:10,(1:10).^3);


%%  Chord panel and axis


%   Tip
handles.Chord.Tip.thickness = 0;
handles.Chord.Tip.axes = handles.tip_foil_axes;
handles.Chord.Tip.thickness_slider = handles.tip_thickness_slider;
handles.Chord.Tip.thickness_edit_text = handles.tip_thickness_edit_text;
handles.Chord.Tip.section = 'N0012';

%   Root
handles.Chord.Root.thickness = 0;
handles.Chord.Root.axes = handles.root_foil_axes;
handles.Chord.Root.thickness_slider = handles.root_thickness_slider;
handles.Chord.Root.thickness_edit_text = handles.root_thickness_edit_text;
handles.Chord.Root.section = 'N0012';

%   Both
handles.Chord.DistPanel.lin_button = handles.chord_linear;
handles.Chord.DistPanel.bell_button = handles.chord_bell;
handles.Chord.DistPanel.num_panels = handles.chord_pan_count;
handles.Chord.DistPanel.bell_param = handles.chord_bell_param;
handles.Chord.DistPanel.NumPanels = handles.chord_pan_count;
handles.Chord.DistPanel.x = linspace(0,1,16);
handles.Chord.DistPanel.minx = 0;
handles.Chord.DistPanel.maxx = 1;

set(handles.chord_pans_buttongroup,'SelectionChangeFcn',@chord_pans_buttongroup_SelectionChangeFcn);
%%  Span panel and axis
set(handles.span_pans_buttongroup,'SelectionChangeFcn',@span_pans_buttongroup_SelectionChangeFcn);


set(handles.Chord.DistPanel.bell_param,'enable','off')
handles.Span.DistPanel.lin_button = handles.span_linear;
handles.Span.DistPanel.bell_button = handles.span_bell;
handles.Span.DistPanel.num_panels = handles.span_pan_count;
handles.Span.DistPanel.bell_param = handles.span_bell_param;
handles.Span.DistPanel.NumPanels = handles.span_pan_count;


handles.Span.tip_pan_count = handles.tip_pan_count;
handles.Span.RoundTips = false;
handles.Span.Cutout.Root = 0;
handles.Span.Cutout.Tip = 0;
set(handles.Span.DistPanel.bell_param,'enable','off')
handles.Span.DistPanel.x = [];
handles.Span.type = 'NREL Phase VI';
handles.Span.LoadFile = false;
handles.Span.LoadFileBox = handles.load_file_name;

handles.Span.axes = handles.blade_geom_axes;

set(handles.tip_thickness_slider,'enable','off');
set(handles.root_thickness_slider,'enable','off');
set(handles.tip_thickness_edit_text,'enable','off');
set(handles.root_thickness_edit_text,'enable','off');
% Update handles structure
guidata(source, handles);

% UIWAIT makes PreProcTab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PreProcTab_OutputFcn(source, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% source    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit2_Callback(source, eventdata, handles)
function edit2_CreateFcn(source, eventdata, handles)
% source    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function root_menu_Callback(source, eventdata, handles)
        str = get(source, 'String');
        handles.Chord.Root.section = str{get(source,'Value')};
        handles.Chord = UpdateFoil(handles.Chord);
        guidata(source, handles)
function root_menu_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function blade_menu_Callback(source, eventdata, handles)
    str = get(source, 'String');
    handles.Span.type = str{get(source,'Value')};
    handles.Span = UpdateSpan(handles.Span);
    handles = ResetSpanCutouts(handles);
    handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
    handles.Span = UpdateSpan(handles.Span);
    guidata(source, handles);
function blade_menu_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
% --- Executes on selection change in tip_menu.
function tip_menu_Callback(source, eventdata, handles)
        % Determine the selected data set.
        str = get(source, 'String');
        handles.Chord.Tip.section = str{get(source,'Value')};
        handles.Chord = UpdateFoil(handles.Chord);
        guidata(source, handles)
function tip_menu_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function tip_thickness_slider_Callback(source, eventdata, handles)
handles.Chord.Tip = Slider2Text(handles.Chord.Tip);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function tip_thickness_slider_CreateFcn(source, eventdata, handles)
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end
function tip_thickness_edit_text_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function tip_thickness_edit_text_Callback(source, eventdata, handles)
handles.Chord.Tip = Text2Slider(handles.Chord.Tip);
%   UpdateFoils
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function root_thickness_slider_Callback(source, eventdata, handles)
handles.Chord.Root = Slider2Text(handles.Chord.Root);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function root_thickness_slider_CreateFcn(source, eventdata, handles)
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end
function root_thickness_edit_text_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function root_thickness_edit_text_Callback(source, eventdata, handles)
handles.Chord.Root = Text2Slider(handles.Chord.Root);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function chord_pans_buttongroup_SelectionChangeFcn(source, eventdata)
%retrieve GUI data, i.e. the handles structure
handles = guidata(source); 
handles.Chord.DistPanel = PanelDistButtonsParam(handles.Chord.DistPanel);
handles.Chord = UpdateFoil(handles.Chord);
%updates the handles structure
guidata(source, handles);
function chord_pan_count_Callback(source, eventdata, handles)
handles.Chord.DistPanel = PanelDistButtonsParam(handles.Chord.DistPanel);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function chord_pan_count_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function chord_bell_param_Callback(source, eventdata, handles)
handles.Chord.DistPanel = PanelDistButtonsParam(handles.Chord.DistPanel);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
function chord_bell_param_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function span_pans_buttongroup_SelectionChangeFcn(source, eventdata)
%retrieve GUI data, i.e. the handles structure
handles = guidata(source); 
handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
%updates the handles structure
guidata(source, handles);

function span_pan_count_Callback(source, eventdata, handles)
handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);

function span_pan_count_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function span_bell_param_Callback(source, eventdata, handles)
handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
function span_bell_param_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

% --- Executes on button press in MakeBlade.
function MakeBlade_Callback(source, eventdata, handles)
handles.Blade = GeomProcess(handles);
guidata(source, handles);
function psi_Callback(source, eventdata, handles)
function psi_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function theta_Callback(source, eventdata, handles)
function theta_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function phi_Callback(source, eventdata, handles)
function phi_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function Ox_Callback(source, eventdata, handles)
function Ox_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function Oy_Callback(source, eventdata, handles)
function Oy_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function Oz_Callback(source, eventdata, handles)
function Oz_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


function nblades_Callback(source, eventdata, handles)
handles = RotorSpec(handles);
function nblades_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function u_Callback(source, eventdata, handles)
function u_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function v_Callback(source, eventdata, handles)
function v_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function w_Callback(source, eventdata, handles)
function w_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function p_Callback(source, eventdata, handles)
function p_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function q_Callback(source, eventdata, handles)
function q_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function r_Callback(source, eventdata, handles)
function r_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in MakeRotor.
function MakeRotor_Callback(source, eventdata, handles)
handles.Rotor = MakeRotor(handles);
guidata(source, handles);


function pitch_axis_Callback(source, eventdata, handles)
function pitch_axis_CreateFcn(source, eventdata, handles)
% source    handles to pitch_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function theta_hub_Callback(source, eventdata, handles)
handles = RotorSpec(handles);
function theta_hub_CreateFcn(source, eventdata, handles)
% source    handles to theta_hub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function cutout_root_Callback(source, eventdata, handles)
handles.Span.Cutout.Root = str2double(get(source,'String'));
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
function cutout_root_CreateFcn(source, eventdata, handles)
% source    handles to cutout_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function cutout_tip_Callback(source, eventdata, handles)
handles.Span.Cutout.Tip = str2double(get(source,'String'));
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
function cutout_tip_CreateFcn(source, eventdata, handles)
% source    handles to cutout_tip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function cutout_reset_Callback(source, eventdata, handles)
handles = ResetSpanCutouts(handles);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
function handles = ResetSpanCutouts(handles)
handles.Span.Cutout.Root = min(handles.Span.RADIUS);
handles.Span.Cutout.Tip = max(handles.Span.RADIUS);
set(handles.cutout_root,'String', num2str(handles.Span.Cutout.Root));
set(handles.cutout_tip,'String', num2str(handles.Span.Cutout.Tip));
function output_case_Callback(source, eventdata, handles)
function output_neu_Callback(source, eventdata, handles)
function output_runfile_Callback(source, eventdata, handles)
function output_svn_Callback(source, eventdata, handles)
function comment_Callback(source, eventdata, handles)
function timestamp_Callback(source, eventdata, handles)
function commit_to_head_Callback(source, eventdata, handles)
function generate_Callback(source, eventdata, handles)
temp = num2str(round(datevec(now)));
handles.datestamp = temp(~isspace(temp));
handles.fullname = get(handles.name,'String');
if get(handles.timestamp,'Value')
    handles.fullname = [handles.fullname handles.datestamp];
end
set(handles.fullname_box,'String',handles.fullname);
handles = MakeNEU(handles);
%handles = WriteCaseFile(handles);
Bodies = handles.Bodies;
%save(['../mat_files/' handles.fullname '.mat'],'Bodies');
DispMsg(handles);
guidata(source, handles);

function blades_as_bodies_Callback(source, eventdata, handles)
handles = RotorSpec(handles);
function name_Callback(source, eventdata, handles)
function name_CreateFcn(source, eventdata, handles)
% source    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function sftp2hpc_Callback(source, eventdata, handles)

function x_Callback(source, eventdata, handles)
function x_CreateFcn(source, eventdata, handles)
% source    handle to Ox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function rpm_Callback(source, eventdata, handles)
function rpm_CreateFcn(source, eventdata, handles)
% source    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function maxP_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function maxP_CreateFcn(source, eventdata, handles)
% source    handle to maxP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function scale_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function scale_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function fluid_temp_Callback(source, eventdata, handles)
handles = GetFluidProps(handles);
[handles CaseInfo] = SimParam(handles);
guidata(source, handles);
function fluid_temp_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function fluid_menu_Callback(source, eventdata, handles)
contents = get(source,'String');
handles.fluid.type = contents{get(source,'Value')};
handles = GetFluidProps(handles);
[handles CaseInfo] = SimParam(handles);
guidata(source, handles);
function fluid_menu_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function uinf_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function uinf_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function vinf_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function vinf_CreateFcn(source, eventdata, handles);

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function winf_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function winf_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function rotor_axis_x_Callback(source, eventdata, handles)
function rotor_axis_x_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function rotor_axis_y_Callback(source, eventdata, handles)
function rotor_axis_y_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function rotor_axis_z_Callback(source, eventdata, handles)
function rotor_axis_z_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function DispMsg(handles)
str = handles.neumsg;
str = strvcat(str,' ');
n = size(handles.casemsg{1},1);
for i = 1:n
    str =  strvcat(str,handles.casemsg{1}{i});
end
set(handles.terminal_output,'String',str);

function dtinit_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function dtinit_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function dtout_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function dtout_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function tmax_Callback(source, eventdata, handles)
omega = str2num(get(handles.rpm,'String'))/60;
set(handles.nturns,'String',num2str(str2double(get(source,'String'))*omega));
[handles CaseInfo] = SimParam(handles);
function tmax_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function nturns_Callback(source, eventdata, handles)
omega = str2num(get(handles.rpm,'String'))/60;
set(handles.tmax,'String',num2str(str2double(get(source,'String'))/omega));
[handles CaseInfo] = SimParam(handles);
function nturns_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function nss_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function nss_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function cflmax_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function cflmax_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function minmodbeta_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function minmodbeta_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function scheme_Callback(source, eventdata, handles)
contents = get(source,'String');
if strcmp(contents{get(source,'Value')}, 'O(2) MUSCL')
    set(handles.minmodbeta,'enable','on');
    set(handles.limiter,'Value',2,'enable','off');   
elseif strcmp(contents{get(source,'Value')}, 'O(1) Upwind')
    set(handles.limiter,'Value',1,'enable','off');   
else
    set(handles.minmodbeta,'enable','off');
    set(handles.limiter,'Value',1,'enable','on');
end
handles.FVM.scheme = contents{get(source,'Value')};
[handles CaseInfo] = SimParam(handles);
function scheme_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function integrator_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function integrator_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function limiter_Callback(source, eventdata, handles)
[handles CaseInfo] = SimParam(handles);
function limiter_CreateFcn(source, eventdata, handles)

if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on selection change in load_file_box.
function load_file_box_Callback(source, eventdata, handles)
handles = CheckFilesForLoading(handles);
guidata(source, handles);
% source    handle to load_file_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(source,'String') returns load_file_box contents as cell array
%        contents{get(source,'Value')} returns selected item from load_file_box


% --- Executes during object creation, after setting all properties.
function load_file_box_CreateFcn(source, eventdata, handles)
% source    handle to load_file_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(source, eventdata, handles)
handles = BodyDataOut(handles);
guidata(source, handles);

% --- Executes on button press in plot_cp.
function plot_cp_Callback(source, eventdata, handles)
handles = PostPlot(handles);
guidata(source, handles);
% source    handle to plot_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function slice_Callback(source, eventdata, handles)
handles.PostComp.use_slice = true;
slice = str2num(get(source,'String'));
r = handles.PostComp.Bodies{1}.Radius(slice);
R = max(handles.PostComp.Bodies{1}.Radius);
set(handles.r_upon_R,'String',num2str(r/R,3));
guidata(source, handles);
% source    handle to slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of slice as text
%        str2double(get(source,'String')) returns contents of slice as a double


% --- Executes during object creation, after setting all properties.
function slice_CreateFcn(source, eventdata, handles)
% source    handle to slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(source, eventdata, handles)
dirname = pwd;
cd ..
system('rm bin_files/*.bin');
    set(handles.inputFiles_listbox,'Value',1);
    set(handles.inputFiles_listbox,'String','Empty',...
	'Value',get(handles.inputFiles_listbox,'Value'));
[status, result] = system(['export LD_LIBRARY_PATH=/usr/lib64; ./main ' handles.fullname '.cas > dump &']);
cd GUI;
[a b] = system('ps | grep main');
while length(b) > 1
    
    [s r] = system('tail --lines=17 ../dump');
    set(handles.terminal_output,'String',r);
    pause(1);
    [a b] = system('ps | grep main');
    
end
[s r] = system('tail --lines=17 ../dump');
set(handles.terminal_output,'String',r);

% source    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(source, eventdata, handles)
gotonewest(handles);
% source    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function lambda_Callback(source, eventdata, handles)

% source    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of lambda as text
%        str2double(get(source,'String')) returns contents of lambda as a double


% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(source, eventdata, handles)
% source    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function r_upon_R_Callback(source, eventdata, handles)
handles.PostComp.use_slice = false;
r_R = str2num(get(source,'String'));


R = max(handles.PostComp.Bodies{1}.Radius);
handles.PostComp.r_slice = R*r_R;
handles.PostComp.r_R = r_R;
slice = interp1(handles.PostComp.Bodies{1}.Radius,1:numel(handles.PostComp.Bodies{1}.Radius),handles.PostComp.r_slice,'nearest');
set(handles.slice,'String',num2str(round(slice)));

guidata(source, handles);
% source    handle to r_upon_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of r_upon_R as text
%        str2double(get(source,'String')) returns contents of r_upon_R as a double


% --- Executes during object creation, after setting all properties.
function r_upon_R_CreateFcn(source, eventdata, handles)
% source    handle to r_upon_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in rounded_tips.
function rounded_tips_Callback(source, eventdata, handles)
% hObject    handle to rounded_tips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state = get(source,'Value');
handles.Span.RoundTips = state;


if state
    set(handles.Span.tip_pan_count,'enable','on')
    set(handles.Span.tip_pan_count,'string','2')
    handles.Span.num_tip_pans = 2;
    
else
    set(handles.Span.tip_pan_count,'enable','off')
    set(handles.Span.tip_pan_count,'string','Num. Pans.')
    handles.Span.num_tip_pans = [];
end
guidata(source, handles);


function tip_pan_count_Callback(source, eventdata, handles)
% hObject    handle to tip_pan_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tip_pan_count as text
handles.Span.num_tip_pans = str2double(get(source,'String'));
guidata(source, handles);
% --- Executes during object creation, after setting all properties.
function tip_pan_count_CreateFcn(source, eventdata, handles)
% hObject    handle to tip_pan_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
guidata(source, handles);



function CEr_R_Callback(hObject, eventdata, handles)
% hObject    handle to CEr_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CEr_R as text
%        str2double(get(hObject,'String')) returns contents of CEr_R as a double


% --- Executes during object creation, after setting all properties.
function CEr_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CEr_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CEthickness_Callback(hObject, eventdata, handles)
% hObject    handle to CEthickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CEthickness as text
%        str2double(get(hObject,'String')) returns contents of CEthickness as a double


% --- Executes during object creation, after setting all properties.
function CEthickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CEthickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CEsweep_Callback(hObject, eventdata, handles)
% hObject    handle to CEsweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CEsweep as text
%        str2double(get(hObject,'String')) returns contents of CEsweep as a double


% --- Executes during object creation, after setting all properties.
function CEsweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CEsweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CEtwist_Callback(hObject, eventdata, handles)
% hObject    handle to CEtwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CEtwist as text
%        str2double(get(hObject,'String')) returns contents of CEtwist as a double


% --- Executes during object creation, after setting all properties.
function CEtwist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CEtwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CEchord_Callback(hObject, eventdata, handles)
% hObject    handle to CEchord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CEchord as text
%        str2double(get(hObject,'String')) returns contents of CEchord as a double


% --- Executes during object creation, after setting all properties.
function CEchord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CEchord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function load_file_name_Callback(source, eventdata, handles)
cla(handles.blade_surf_axes,'reset');
cla(handles.blade_geom_axes,'reset');
if ~handles.Span.LoadFile
set(handles.load_file_name,'enable','off');
end
handles.Span = UpdateSpan(handles.Span);
% hObject    handle to load_file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of load_file_name as text
%        str2double(get(hObject,'String')) returns contents of load_file_name as a double
guidata(source, handles);

% --- Executes during object creation, after setting all properties.
function load_file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to load_file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% 
% function root_thickness_edit_text_CreateFcn(source, eventdata, handles)
% if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(source,'BackgroundColor','white');
% end
% function root_thickness_edit_text_Callback(source, eventdata, handles)
% handles.Chord.Root = Text2Slider(handles.Chord.Root);
% handles.Chord = UpdateFoil(handles.Chord);
% guidata(source, handles);
