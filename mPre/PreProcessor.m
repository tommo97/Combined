function varargout = PreProcessor(varargin)

% PREPROCESSOR M-file for PreProcessor.fig
%      PREPROCESSOR, by itself, creates a new PREPROCESSOR or raises the existing
%      singleton*.
%
%      H = PREPROCESSOR returns the handle to a new PREPROCESSOR or the handle to
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

% Last Modified by GUIDE v2.5 17-Nov-2009 19:18:06

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
% source    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreProcessor (see VARARGIN)

% Choose default command line output for PreProcessor
handles.output = source;

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
% source    handle to figure
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
handles = MakeRotor(handles);



function pitch_axis_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pitch_axis as text
%        str2double(get(hObject,'String')) returns contents of pitch_axis as a double


% --- Executes during object creation, after setting all properties.
function pitch_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta_hub_Callback(hObject, eventdata, handles)
% hObject    handle to theta_hub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta_hub as text
%        str2double(get(hObject,'String')) returns contents of theta_hub as a double


% --- Executes during object creation, after setting all properties.
function theta_hub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_hub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
