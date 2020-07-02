function varargout = MMM_preferences(varargin)
% MMM_PREFERENCES M-file for MMM_preferences.fig
%      MMM_PREFERENCES, by itself, creates a new MMM_PREFERENCES or raises the existing
%      singleton*.
%
%      H = MMM_PREFERENCES returns the handle to a new MMM_PREFERENCES or the handle to
%      the existing singleton*.
%
%      MMM_PREFERENCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MMM_PREFERENCES.M with the given input arguments.
%
%      MMM_PREFERENCES('Property','Value',...) creates a new MMM_PREFERENCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MMM_preferences_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MMM_preferences_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MMM_preferences

% Last Modified by GUIDE v2.5 24-Jul-2013 15:40:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MMM_preferences_OpeningFcn, ...
                   'gui_OutputFcn',  @MMM_preferences_OutputFcn, ...
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


% --- Executes just before MMM_preferences is made visible.
function MMM_preferences_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MMM_preferences (see VARARGIN)

% Choose default command line output for MMM_preferences
handles.output = hObject;

global web_adr
% global MMM_icon
global browser
global third_party
global general
global ENM_param

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

set(handles.edit_SFX,'String',web_adr.SFX);
set(handles.edit_Modeller_call,'String',third_party.modeller_version);

handles.SFX=web_adr.SFX;

if general.large_windows,
    set(handles.checkbox_large_hierarchy_window,'Value',1);
else
    set(handles.checkbox_large_hierarchy_window,'Value',0);
end;
if strcmp(web_adr.PDB,web_adr.PDB_EU),
    set(handles.radiobutton_PDB_EU,'Value',1),
end;
if strcmp(web_adr.PDB,web_adr.PDB_US),
    set(handles.radiobutton_PDB_USA,'Value',1),
end;
if strcmp(web_adr.PDB,web_adr.PDB_J),
    set(handles.radiobutton_PDB_J,'Value',1),
end;
if strcmpi(browser,'system'),
    set(handles.radiobutton_system_browser,'Value',1),
end;
if strcmpi(browser,'matlab'),
    set(handles.radiobutton_Matlab_browser,'Value',1),
end;
if strcmpi(browser,'mixed'),
    set(handles.radiobutton_mixed,'Value',1),
end;

set(handles.edit_cores,'String',sprintf('%i',general.cpu));
handles.cpu=general.cpu;
 
switch ENM_param.parametrization
    case 'ed-ENM'
        parm=1;
    case 'Jeschke'
        parm=2;
    case 'Hinsen'
        parm=3;
    case 'cutoff10'
        parm=4;
    case 'cutoff13'
        parm=5;
    case 'ed-ENM-p'
        parm=6;
    otherwise
        parm=1;
end;
set(handles.popupmenu_ANM,'Value',parm);

if ENM_param.imANM,
    set(handles.checkbox_imANM,'Value',1);
else
    set(handles.checkbox_imANM,'Value',0);
end;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MMM_preferences wait for user response (see UIRESUME)
% uiwait(handles.my_figure);


% --- Outputs from this function are returned to the command line.
function varargout = MMM_preferences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_SFX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SFX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SFX as text
%        str2double(get(hObject,'String')) returns contents of edit_SFX as a double

% --- Executes during object creation, after setting all properties.
function edit_SFX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SFX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global queries
global general
global browser
global ENM_param
global hMain

load preferences

if get(handles.checkbox_large_hierarchy_window,'Value'),
    general.large_windows = true;
    hMain.large = true;
else
    general.large_windows = false;
    hMain.large = false;
end;
if get(handles.radiobutton_PDB_EU,'Value'),
    web_adr.PDB=web_adr.PDB_EU;
    queries.PDB_structures=queries.PDB_structures_EU;
end;
if get(handles.radiobutton_PDB_USA,'Value'),
    web_adr.PDB=web_adr.PDB_US;
    queries.PDB_structures=queries.PDB_structures_US;
end;
if get(handles.radiobutton_PDB_J,'Value'),
    web_adr.PDB=web_adr.PDB_J;
    queries.PDB_structures=queries.PDB_structures_J;
end;
if get(handles.radiobutton_system_browser,'Value'),
    browser='system';
end;
if get(handles.radiobutton_Matlab_browser,'Value'),
    browser='matlab';
end;
if get(handles.radiobutton_mixed,'Value'),
    browser='mixed';
end;
web_adr.SFX=get(handles.edit_SFX,'String');
third_party.modeller_version=get(handles.edit_Modeller_call,'String');

general.cpu=handles.cpu;
user_preferences.cpu=general.cpu;

parm=get(handles.popupmenu_ANM,'Value');
switch parm
    case 1
        ENM_param.parametrization='ed-ENM';
    case 2
        ENM_param.parametrization='Jeschke';
    case 3
        ENM_param.parametrization='Hinsen';
    case 4
        ENM_param.parametrization='cutoff10';
    case 5
        ENM_param.parametrization='cutoff13';
    case 6
        ENM_param.parametrization='ed-ENM-p';
    otherwise
        ENM_param.parametrization='ed-ENM';        
end;
user_preferences.ANM.parametrization=ENM_param.parametrization;
ENM_param=set_ANM(ENM_param);

imANM_flag=get(handles.checkbox_imANM,'Value');
if imANM_flag,
    ENM_param.imANM=true;
else
    ENM_param.imANM=false;
end;
user_preferences.ANM.imANM=ENM_param.imANM;

ENM_param.mass_weighting=get(handles.checkbox_mass_weighting,'Value');
user_preferences.ANM.mass_weighting=ENM_param.mass_weighting;

user_preferences.SFX=web_adr.SFX;
user_preferences.PDB=web_adr.PDB;
user_preferences.PDB_structures=queries.PDB_structures;
user_preferences.browser=browser;
user_preferences.Modeller_call=third_party.modeller_version;
user_preferences.large = hMain.large;
save([general.rootdir 'preferences.mat'],'virgin','user_preferences');
delete(handles.my_figure);

% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global queries
global browser
global third_party
global general
global ENM_param
global hMain

if get(handles.checkbox_large_hierarchy_window,'Value'),
    general.large_windows = true;
    hMain.large = true;
else
    general.large_windows = false;
    hMain.large = false;
end;
if get(handles.radiobutton_PDB_EU,'Value'),
    web_adr.PDB=web_adr.PDB_EU;
    queries.PDB_structures=queries.PDB_structures_EU;
end;
if get(handles.radiobutton_PDB_USA,'Value'),
    web_adr.PDB=web_adr.PDB_US;
    queries.PDB_structures=queries.PDB_structures_US;
end;
if get(handles.radiobutton_PDB_J,'Value'),
    web_adr.PDB=web_adr.PDB_J;
    queries.PDB_structures=queries.PDB_structures_J;
end;
if get(handles.radiobutton_system_browser,'Value'),
    browser='system';
end;
if get(handles.radiobutton_Matlab_browser,'Value'),
    browser='matlab';
end;
if get(handles.radiobutton_mixed,'Value'),
    browser='mixed';
end;

general.cpu=handles.cpu;

parm=get(handles.popupmenu_ANM,'Value');
switch parm
    case 1
        ENM_param.parametrization='ed-ENM';
    case 2
        ENM_param.parametrization='Jeschke';
    case 3
        ENM_param.parametrization='Hinsen';
    case 4
        ENM_param.parametrization='cutoff10';
    case 5
        ENM_param.parametrization='cutoff13';
    case 6
        ENM_param.parametrization='ed-ENM-p';
    otherwise
        ENM_param.parametrization='ed-ENM';        
end;

ENM_param.mass_weighting=get(handles.checkbox_mass_weighting,'Value');

ENM_param=set_ANM(ENM_param);
imANM_flag=get(handles.checkbox_imANM,'Value');
if imANM_flag,
    ENM_param.imANM=true;
else
    ENM_param.imANM=false;
end;


web_adr.SFX=get(handles.edit_SFX,'String');
third_party.modeller_version=get(handles.edit_Modeller_call,'String');
delete(handles.my_figure);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.my_figure);



function edit_Modeller_call_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Modeller_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Modeller_call as text
%        str2double(get(hObject,'String')) returns contents of edit_Modeller_call as a double


% --- Executes during object creation, after setting all properties.
function edit_Modeller_call_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Modeller_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cores_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cores as text
%        str2double(get(hObject,'String')) returns contents of edit_cores as a double

[v,handback]=edit_update_MMM(handles,hObject,1,100,3,'%i',true);
handles.cpu=v;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_cores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_ANM.
function popupmenu_ANM_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ANM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ANM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ANM


% --- Executes during object creation, after setting all properties.
function popupmenu_ANM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ANM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_imANM.
function checkbox_imANM_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_imANM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_imANM


% --- Executes on button press in checkbox_mass_weighting.
function checkbox_mass_weighting_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mass_weighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mass_weighting


% --- Executes on button press in checkbox_large_hierarchy_window.
function checkbox_large_hierarchy_window_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_large_hierarchy_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_large_hierarchy_window
