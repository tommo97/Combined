function varargout = PreProcessor(varargin)

% PREPROCESSOR M-file for PreProcessor.fig
%      PREPROCESSOR, by itself, creates a new PREPROCESSOR or raises the existing
%      singleton*.
%
%      H = PREPROCESSOR returns the handles to a new PREPROCESSOR or the handles to
%      the existing singleton*.
%
%      PREPROCESSOR('CALLBACK',source,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSOR.M with the given input arguments.
%
%      PREPROCESSOR('Property','Value',...) creates a new PREPROCESSOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PreProcessor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PreProcessor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PreProcessor

% Last Modified by GUIDE v2.5 19-Nov-2009 18:36:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreProcessor_OpeningFcn, ...
                   'gui_OutputFcn',  @PreProcessor_OutputFcn, ...
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
% --- Executes just before PreProcessor is made visible.
function PreProcessor_OpeningFcn(source, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% source    handles to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreProcessor (see VARARGIN)

% Choose default command line output for PreProcessor
handles.output = source;



%%  Set up Output panel defaults
set(handles.output_case,'Value',true);
set(handles.output_neu,'Value',true);
set(handles.output_runfile,'Value',true);
set(handles.commit_to_head,'enable','off');
set(handles.output_svn,'enable','off');
set(handles.timestamp,'Value',true);


handles.Origin = [0;0;0];
handles.Attitude = [0;0;0];

handles.Chord.Tip.thickness = 0;
handles.Chord.Tip.axes = handles.tip_foil_axes;
handles.Chord.Tip.thickness_slider = handles.tip_thickness_slider;
handles.Chord.Tip.thickness_edit_text = handles.tip_thickness_edit_text;
handles.Chord.Tip.section = 'N0012';

handles.Chord.Root.thickness = 0;
handles.Chord.Root.axes = handles.root_foil_axes;
handles.Chord.Root.thickness_slider = handles.root_thickness_slider;
handles.Chord.Root.thickness_edit_text = handles.root_thickness_edit_text;
handles.Chord.Root.section = 'N0012';

set(handles.chord_pans_buttongroup,'SelectionChangeFcn',@chord_pans_buttongroup_SelectionChangeFcn);

handles.Chord.DistPanel.lin_button = handles.chord_linear;
handles.Chord.DistPanel.bell_button = handles.chord_bell;
handles.Chord.DistPanel.num_panels = handles.chord_pan_count;
handles.Chord.DistPanel.bell_param = handles.chord_bell_param;
handles.Chord.DistPanel.NumPanels = handles.chord_pan_count;
handles.Chord.DistPanel.x = linspace(0,1,16);
handles.Chord.DistPanel.minx = 0;
handles.Chord.DistPanel.maxx = 1;

set(handles.span_pans_buttongroup,'SelectionChangeFcn',@span_pans_buttongroup_SelectionChangeFcn);


set(handles.Chord.DistPanel.bell_param,'enable','off')
handles.Span.DistPanel.lin_button = handles.span_linear;
handles.Span.DistPanel.bell_button = handles.span_bell;
handles.Span.DistPanel.num_panels = handles.span_pan_count;
handles.Span.DistPanel.bell_param = handles.span_bell_param;
handles.Span.DistPanel.NumPanels = handles.span_pan_count;
handles.Span.Cutout.Root = 0;
handles.Span.Cutout.Tip = 0;
set(handles.Span.DistPanel.bell_param,'enable','off')
handles.Span.DistPanel.x = [];
handles.Span.type = 'NREL UAE';


handles.Span.axes = handles.blade_geom_axes;

set(handles.tip_thickness_slider,'enable','off');
set(handles.root_thickness_slider,'enable','off');
set(handles.tip_thickness_edit_text,'enable','off');
set(handles.root_thickness_edit_text,'enable','off');
% Update handles structure
guidata(source, handles);
% UIWAIT makes PreProcessor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = PreProcessor_OutputFcn(source, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% source    handles to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in root_menu.
function root_menu_Callback(source, eventdata, handles)
        str = get(source, 'String');
        handles.Chord.Root.section = str{get(source,'Value')};
        handles.Chord = UpdateFoil(handles.Chord);
        guidata(source, handles)
% --- Executes during object creation, after setting all properties.
function root_menu_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
% --- Executes on selection change in blade_menu.
function blade_menu_Callback(source, eventdata, handles)
    str = get(source, 'String');
    handles.Span.type = str{get(source,'Value')};
    handles.Span = UpdateSpan(handles.Span);
    handles = ResetSpanCutouts(handles);
    handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
    handles.Span = UpdateSpan(handles.Span);
    guidata(source, handles);

% --- Executes during object creation, after setting all properties.
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
% --- Executes during object creation, after setting all properties.
function tip_menu_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on slider movement.
function tip_thickness_slider_Callback(source, eventdata, handles)
handles.Chord.Tip = Slider2Text(handles.Chord.Tip);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
% --- Executes during object creation, after setting all properties.
function tip_thickness_slider_CreateFcn(source, eventdata, handles)
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
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
% --- Executes during object creation, after setting all properties.
function root_thickness_slider_CreateFcn(source, eventdata, handles)
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
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
% --- Executes during object creation, after setting all properties.
function chord_pan_count_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
function chord_bell_param_Callback(source, eventdata, handles)
handles.Chord.DistPanel = PanelDistButtonsParam(handles.Chord.DistPanel);
handles.Chord = UpdateFoil(handles.Chord);
guidata(source, handles);
% --- Executes during object creation, after setting all properties.
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
% --- Executes during object creation, after setting all properties.
function span_pan_count_CreateFcn(source, eventdata, handles)
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function span_bell_param_Callback(source, eventdata, handles)
handles.Span.DistPanel = PanelDistButtonsParam(handles.Span.DistPanel);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
% --- Executes during object creation, after setting all properties.
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
% source    handles to pitch_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of pitch_axis as text
%        str2double(get(source,'String')) returns contents of pitch_axis as a double


% --- Executes during object creation, after setting all properties.
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
% source    handles to theta_hub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of theta_hub as text
%        str2double(get(source,'String')) returns contents of theta_hub as a double


% --- Executes during object creation, after setting all properties.
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
% source    handles to cutout_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of cutout_root as text
%        str2double(get(source,'String')) returns contents of cutout_root as a double


% --- Executes during object creation, after setting all properties.
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
% source    handles to cutout_tip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of cutout_tip as text
%        str2double(get(source,'String')) returns contents of cutout_tip as a double


% --- Executes during object creation, after setting all properties.
function cutout_tip_CreateFcn(source, eventdata, handles)
% source    handles to cutout_tip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in cutout_reset.
function cutout_reset_Callback(source, eventdata, handles)
handles = ResetSpanCutouts(handles);
handles.Span = UpdateSpan(handles.Span);
guidata(source, handles);
% source    handles to cutout_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = ResetSpanCutouts(handles)
handles.Span.Cutout.Root = min(handles.Span.RADIUS);
handles.Span.Cutout.Tip = max(handles.Span.RADIUS);
set(handles.cutout_root,'String', num2str(handles.Span.Cutout.Root));
set(handles.cutout_tip,'String', num2str(handles.Span.Cutout.Tip));


% --- Executes on button press in output_case.
function output_case_Callback(source, eventdata, handles)
% source    handle to output_case (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of output_case


% --- Executes on button press in output_neu.
function output_neu_Callback(source, eventdata, handles)
% source    handle to output_neu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of output_neu


% --- Executes on button press in output_runfile.
function output_runfile_Callback(source, eventdata, handles)
% source    handle to output_runfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of output_runfile


% --- Executes on button press in output_svn.
function output_svn_Callback(source, eventdata, handles)
% source    handle to output_svn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of output_svn


% --- Executes on button press in comment.
function comment_Callback(source, eventdata, handles)
% source    handle to comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of comment


% --- Executes on button press in timestamp.
function timestamp_Callback(source, eventdata, handles)
% source    handle to timestamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of timestamp


% --- Executes on button press in commit_to_head.
function commit_to_head_Callback(source, eventdata, handles)
% source    handle to commit_to_head (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of commit_to_head


% --- Executes on button press in generate.
function generate_Callback(source, eventdata, handles)
MakeNEU(handles.Rotor.Blade,get(handles.name,'String'),handles.Rotor.Split)
% source    handle to generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in blades_as_bodies.
function blades_as_bodies_Callback(source, eventdata, handles)
% source    handle to blades_as_bodies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of blades_as_bodies



function name_Callback(source, eventdata, handles)
% source    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of name as text
%        str2double(get(source,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(source, eventdata, handles)
% source    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on button press in sftp2hpc.
function sftp2hpc_Callback(source, eventdata, handles)
% source    handle to sftp2hpc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(source,'Value') returns toggle state of sftp2hpc



function edit35_Callback(source, eventdata, handles)
% source    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit35 as text
%        str2double(get(source,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(source, eventdata, handles)
% source    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function edit36_Callback(source, eventdata, handles)
% source    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit36 as text
%        str2double(get(source,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(source, eventdata, handles)
% source    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function edit37_Callback(source, eventdata, handles)
% source    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit37 as text
%        str2double(get(source,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(source, eventdata, handles)
% source    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function x_Callback(source, eventdata, handles)
% source    handle to Ox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of Ox as text
%        str2double(get(source,'String')) returns contents of Ox as a double


% --- Executes during object creation, after setting all properties.
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
% source    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of rpm as text
%        str2double(get(source,'String')) returns contents of rpm as a double


% --- Executes during object creation, after setting all properties.
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
% source    handle to maxP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of maxP as text
%        str2double(get(source,'String')) returns contents of maxP as a double


% --- Executes during object creation, after setting all properties.
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
% source    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of scale as text
%        str2double(get(source,'String')) returns contents of scale as a double


% --- Executes during object creation, after setting all properties.
function scale_CreateFcn(source, eventdata, handles)
% source    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function fluid_temp_Callback(source, eventdata, handles)
handles = GetFluidProps(handles);
guidata(source, handles);
% source    handle to fluid_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of fluid_temp as text
%        str2double(get(source,'String')) returns contents of fluid_temp as a double


% --- Executes during object creation, after setting all properties.
function fluid_temp_CreateFcn(source, eventdata, handles)
% source    handle to fluid_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on selection change in fluid_menu.
function fluid_menu_Callback(source, eventdata, handles)
contents = get(source,'String');
handles.fluid.type = contents{get(source,'Value')};
handles = GetFluidProps(handles);
guidata(source, handles);
% source    handle to fluid_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(source,'String') returns fluid_menu contents as cell array
%        contents{get(source,'Value')} returns selected item from fluid_menu


% --- Executes during object creation, after setting all properties.
function fluid_menu_CreateFcn(source, eventdata, handles)
% source    handle to fluid_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function edit45_Callback(source, eventdata, handles)
% source    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit45 as text
%        str2double(get(source,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(source, eventdata, handles)
% source    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function edit46_Callback(source, eventdata, handles)
% source    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit46 as text
%        str2double(get(source,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(source, eventdata, handles)
% source    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function edit47_Callback(source, eventdata, handles)
% source    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of edit47 as text
%        str2double(get(source,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(source, eventdata, handles)
% source    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function rotor_axis_y_Callback(source, eventdata, handles)
% source    handle to rotor_axis_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of rotor_axis_y as text
%        str2double(get(source,'String')) returns contents of rotor_axis_y as a double


% --- Executes during object creation, after setting all properties.
function rotor_axis_y_CreateFcn(source, eventdata, handles)
% source    handle to rotor_axis_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function rotor_axis_z_Callback(source, eventdata, handles)
% source    handle to rotor_axis_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of rotor_axis_z as text
%        str2double(get(source,'String')) returns contents of rotor_axis_z as a double


% --- Executes during object creation, after setting all properties.
function rotor_axis_z_CreateFcn(source, eventdata, handles)
% source    handle to rotor_axis_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end



function rotor_axis_x_Callback(source, eventdata, handles)
% source    handle to rotor_axis_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(source,'String') returns contents of rotor_axis_x as text
%        str2double(get(source,'String')) returns contents of rotor_axis_x as a double


% --- Executes during object creation, after setting all properties.
function rotor_axis_x_CreateFcn(source, eventdata, handles)
% source    handle to rotor_axis_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end
