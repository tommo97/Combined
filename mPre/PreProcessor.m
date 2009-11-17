function varargout = pre(varargin)

% PRE M-file for pre.fig
%      PRE, by itself, creates a new PRE or raises the existing
%      singleton*.
%
%      H = PRE returns the handle to a new PRE or the handle to
%      the existing singleton*.
%
%      PRE('CALLBACK',source,eventData,handles,...) calls the local
%      function named CALLBACK in PRE.M with the given input arguments.
%
%      PRE('Property','Value',...) creates a new PRE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pre_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pre_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pre

% Last Modified by GUIDE v2.5 16-Nov-2009 20:06:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pre_OpeningFcn, ...
                   'gui_OutputFcn',  @pre_OutputFcn, ...
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




% --- Executes just before pre is made visible.
function pre_OpeningFcn(source, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% source    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pre (see VARARGIN)

% Choose default command line output for pre
handles.output = source;
handles.SectionShape.Root.thickness = 0;
handles.SectionShape.Tip.thickness = 0;
handles.SectionShape.Root.axes = handles.root_foil_axes;
handles.SectionShape.Root.thickness_slider = handles.root_thickness_slider;
handles.SectionShape.Root.thickness_edit_text = handles.root_thickness_edit_text;
handles.SectionShape.Tip.axes = handles.tip_foil_axes;
handles.SectionShape.Tip.thickness_slider = handles.tip_thickness_slider;
handles.SectionShape.Tip.thickness_edit_text = handles.tip_thickness_edit_text;
set(handles.tip_thickness_slider,'enable','off');
set(handles.root_thickness_slider,'enable','off');
set(handles.tip_thickness_edit_text,'enable','off');
set(handles.root_thickness_edit_text,'enable','off');
% Update handles structure
guidata(source, handles);

% UIWAIT makes pre wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pre_OutputFcn(source, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% source    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in root_menu.
function root_menu_Callback(source, eventdata, handles)
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        handles.SectionShape.Root = SetFoilData(str{val},handles.SectionShape.Root);


% --- Executes during object creation, after setting all properties.
function root_menu_CreateFcn(source, eventdata, handles)
% source    handle to root_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(source, eventdata, handles)
% source    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(source,'String') returns popupmenu2 contents as cell array
%        contents{get(source,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(source, eventdata, handles)
% source    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on selection change in blade_menu.
function blade_menu_Callback(source, eventdata, handles)
            str = get(source, 'String');
        val = get(source,'Value');
        BladeChordTwist.axes = handles.blade_geom_axes;
        BladeChordTwist = SetBladeData(str{val},BladeChordTwist);


% --- Executes during object creation, after setting all properties.
function blade_menu_CreateFcn(source, eventdata, handles)
% source    handle to blade_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


% --- Executes on selection change in tip_menu.
function tip_menu_Callback(source, eventdata, handles)
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        handles.SectionShape.Tip.axes = handles.tip_foil_axes;
        handles.SectionShape.Tip.thickness_slider = handles.tip_thickness_slider;
        handles.SectionShape.Tip.thickness_edit_text = handles.tip_thickness_edit_text;
        handles.SectionShape.Tip = SetFoilData(str{val},handles.SectionShape.Tip);


% --- Executes during object creation, after setting all properties.
function tip_menu_CreateFcn(source, eventdata, handles)
% source    handle to tip_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end


    function foil = SetFoilData(type,foil)
        foil.x = linspace(0,1);
        [Foil Data] = MakeFoil(foil.x,foil.thickness);
        set(foil.thickness_slider,'enable','off');
        set(foil.thickness_edit_text,'enable','off');
        % Set current data to the selected data set.
        switch type;
            case 'NREL S809' % User selects Peaks.
                foil.US = Data.S809.US;
                foil.LS = Data.S809.LS;
                foil.X = Data.S809.X;
            case 'NREL S814' % User selects Membrane.
                foil.US = Data.S814.US;
                foil.LS = Data.S814.LS;
                foil.X = Data.S814.X;
            case 'NACA 00xx' % User selects Membrane.
                set(foil.thickness_slider,'enable','on');
                set(foil.thickness_edit_text,'enable','on');
                foil.US = Data.N00xx.US;
                foil.LS = Data.N00xx.LS;
                foil.X = Data.N00xx.X;    
            otherwise
                foil.US = Data.N0012.US;
                foil.LS = Data.N0012.LS;
                foil.X = Data.N0012.X;
        end
        PlotFoil(foil.X,foil.US,foil.LS,foil.axes);
        
        
        function blade = SetBladeData(type,blade)
        % Set current data to the selected data set.
        switch type;
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
        end
        PlotBlade(blade.RADIUS,blade.CHORD,blade.THETA,blade.axes);
   
    
        function PlotFoil(ix,ius,ils,iaxes)
        hold(iaxes,'off')
        plot(iaxes,ix,ius);
        hold(iaxes,'on');
        plot(iaxes,ix,ils);
        axis(iaxes,'equal')
        hold(iaxes,'off')
        
        
        function PlotBlade(irad,icrd,ith,iaxes)
        hold(iaxes,'off')
        plotyy(iaxes,irad,icrd,irad,ith);
  


% --- Executes on slider movement.
function tip_thickness_slider_Callback(source, eventdata, handles)
% source    handle to tip_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SectionShape.Tip.thickness = get(source,'Value');
set(handles.tip_thickness_edit_text,'String', num2str(handles.SectionShape.Tip.thickness));
handles.SectionShape.Tip = SetFoilData('NACA 00xx',handles.SectionShape.Tip);
guidata(source, handles);
% Hints: get(source,'Value') returns position of slider
%        get(source,'Min') and get(source,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tip_thickness_slider_CreateFcn(source, eventdata, handles)
% source    handle to tip_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function tip_thickness_edit_text_CreateFcn(source, eventdata, handles)
% source    handle to tip_thickness_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function tip_thickness_edit_text_Callback(source, eventdata, handles)
        % source    handle to tip_thickness_edit_text (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        %get the string for the editText component
        sliderValue = get(handles.tip_thickness_edit_text,'String');
        %convert from string to number if possible, otherwise returns empty
        %if user inputs something is not a number, or if the input is less than 0
        %or greater than 100, then the slider value defaults to 0
        if (isempty(str2num(sliderValue)) || str2num(sliderValue) < 0)
            set(handles.tip_thickness_slider,'Value',0);
            set(handles.tip_thickness_edit_text,'String','0');
            handles.SectionShape.Tip.thickness = 0;
        elseif (str2num(sliderValue) > 1)
            set(handles.tip_thickness_slider,'Value',1);
            set(handles.tip_thickness_edit_text,'String','1');
            handles.SectionShape.Tip.thickness = 1;
        else
            set(handles.tip_thickness_slider,'Value',str2num(sliderValue));
            handles.SectionShape.Tip.thickness = str2num(sliderValue);
        end
        handles.SectionShape.Tip = SetFoilData('NACA 00xx',handles.SectionShape.Tip);
% --- Executes on slider movement.
function root_thickness_slider_Callback(source, eventdata, handles)
% source    handle to root_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SectionShape.Root.thickness = get(source,'Value');
set(handles.root_thickness_edit_text,'String', num2str(handles.SectionShape.Root.thickness));
handles.SectionShape.Root = SetFoilData('NACA 00xx',handles.SectionShape.Root);
guidata(source, handles);
% Hints: get(source,'Value') returns position of slider
%        get(source,'Min') and get(source,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function root_thickness_slider_CreateFcn(source, eventdata, handles)
% source    handle to root_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function root_thickness_edit_text_CreateFcn(source, eventdata, handles)
% source    handle to root_thickness_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(source,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(source,'BackgroundColor','white');
end

function root_thickness_edit_text_Callback(source, eventdata, handles)
% source    handle to root_thickness_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the string for the editText component
sliderValue = get(handles.root_thickness_edit_text,'String');
%convert from string to number if possible, otherwise returns empty
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(str2num(sliderValue)) || str2num(sliderValue) < 0)
    set(handles.root_thickness_slider,'Value',0);
    set(handles.root_thickness_edit_text,'String','0');
    handles.SectionShape.Root.thickness = 0;
elseif (str2num(sliderValue) > 1)
    set(handles.root_thickness_slider,'Value',1);
    set(handles.root_thickness_edit_text,'String','1');
    handles.SectionShape.Root.thickness = 1;
else
    set(handles.root_thickness_slider,'Value',str2num(sliderValue));
    handles.SectionShape.Root.thickness = str2num(sliderValue);
end
handles.SectionShape.Root = SetFoilData('NACA 00xx',handles.SectionShape.Root);