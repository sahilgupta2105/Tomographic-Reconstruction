function varargout = tdlas(varargin)
% TDLAS MATLAB code for tdlas.fig
%      TDLAS, by itself, creates a new TDLAS or raises the existing
%      singleton*.
%
%      H = TDLAS returns the handle to a new TDLAS or the handle to
%      the existing singleton*.
%
%      TDLAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TDLAS.M with the given input arguments.
%
%      TDLAS('Property','Value',...) creates a new TDLAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tdlas_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tdlas_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tdlas

% Last Modified by GUIDE v2.5 10-Jan-2017 22:27:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tdlas_OpeningFcn, ...
                   'gui_OutputFcn',  @tdlas_OutputFcn, ...
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


% --- Executes just before tdlas is made visible.
function tdlas_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tdlas (see VARARGIN)

% Choose default command line output for tdlas
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tdlas wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tdlas_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_pressure_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pressure as text
%        str2double(get(hObject,'String')) returns contents of edit_pressure as a double


% --- Executes during object creation, after setting all properties.
function edit_pressure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_line_strength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_line_strength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_line_strength as text
%        str2double(get(hObject,'String')) returns contents of edit_line_strength as a double


% --- Executes during object creation, after setting all properties.
function edit_line_strength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_line_strength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freq as text
%        str2double(get(hObject,'String')) returns contents of edit_freq as a double


% --- Executes during object creation, after setting all properties.
function edit_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_temperature_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_temperature as text
%        str2double(get(hObject,'String')) returns contents of edit_temperature as a double


% --- Executes during object creation, after setting all properties.
function edit_temperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mm as text
%        str2double(get(hObject,'String')) returns contents of edit_mm as a double


% --- Executes during object creation, after setting all properties.
function edit_mm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bc_air_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bc_air (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bc_air as text
%        str2double(get(hObject,'String')) returns contents of edit_bc_air as a double


% --- Executes during object creation, after setting all properties.
function edit_bc_air_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bc_air (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bc_gas_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bc_gas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bc_gas as text
%        str2double(get(hObject,'String')) returns contents of edit_bc_gas as a double


% --- Executes during object creation, after setting all properties.
function edit_bc_gas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bc_gas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mesh as text
%        str2double(get(hObject,'String')) returns contents of edit_mesh as a double


% --- Executes during object creation, after setting all properties.
function edit_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lasers_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lasers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lasers as text
%        str2double(get(hObject,'String')) returns contents of edit_lasers as a double


% --- Executes during object creation, after setting all properties.
function edit_lasers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lasers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_angles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_angles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_angles as text
%        str2double(get(hObject,'String')) returns contents of edit_angles as a double


% --- Executes during object creation, after setting all properties.
function edit_angles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_angles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inverse.
function inverse_Callback(hObject, eventdata, handles)
% hObject    handle to inverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inverse
cla(handles.fig_inv)
cla(handles.fig_error)
inverse

% --- Executes on key press with focus on inverse and none of its controls.
function inverse_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to inverse (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.fig_mesh)
cla(handles.fig_conc)
cla(handles.fig_od)
cla(handles.fig_inv)
cla(handles.fig_error)


% --- Executes on button press in generate_data.
function generate_data_Callback(hObject, eventdata, handles)
% hObject    handle to generate_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.fig_mesh)
cla(handles.fig_conc)
cla(handles.fig_od)
code



function epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon as text
%        str2double(get(hObject,'String')) returns contents of epsilon as a double


% --- Executes during object creation, after setting all properties.
function epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delta_Callback(hObject, eventdata, handles)
% hObject    handle to delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta as text
%        str2double(get(hObject,'String')) returns contents of delta as a double


% --- Executes during object creation, after setting all properties.
function delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
help

% --- Executes on selection change in popup_conc.
function popup_conc_Callback(hObject, eventdata, handles)
% hObject    handle to popup_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_conc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_conc


% --- Executes during object creation, after setting all properties.
function popup_conc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
