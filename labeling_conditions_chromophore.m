function varargout = labeling_conditions_chromophore(varargin)
% LABELING_CONDITIONS M-file for labeling_conditions.fig
%      LABELING_CONDITIONS, by itself, creates a new LABELING_CONDITIONS or raises the existing
%      singleton*.
%
%      H = LABELING_CONDITIONS returns the handle to a new LABELING_CONDITIONS or the handle to
%      the existing singleton*.
%
%      LABELING_CONDITIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABELING_CONDITIONS.M with the given input arguments.
%
%      LABELING_CONDITIONS('Property','Value',...) creates a new LABELING_CONDITIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before labeling_conditions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to labeling_conditions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help labeling_conditions

% Last Modified by GUIDE v2.5 17-Jan-2013 07:58:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labeling_conditions_OpeningFcn, ...
                   'gui_OutputFcn',  @labeling_conditions_OutputFcn, ...
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


% --- Executes just before labeling_conditions is made visible.
function labeling_conditions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labeling_conditions (see VARARGIN)

% Choose default command line output for labeling_conditions

global hMain
% global MMM_icon

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;
handles.eventdata = eventdata;

handles=initialize_popup(handles);
handles=update_selection(handles);

if isfield(hMain,'label_selection')
    title=sprintf('Set labeling conditions for %s',hMain.label_selection);
else
    title='Set labeling conditions for selected residues';
end
set(hObject,'Name',title);

s = load('helpicon.mat');
set(handles.pushbutton_help,'CData',s.cdata);

locked=get(hMain.uitoggletool_lock,'State');
if strcmpi(locked,'on')
    set(handles.checkbox_special,'Enable','off');
    set(handles.popupmenu_library,'Enable','off');
else
    set(handles.checkbox_special,'Enable','on');
    set(handles.popupmenu_library,'Enable','on');
end

% Update handles structure
guidata(hObject, handles);

uiwait(handles.labeling_window);

% UIWAIT makes labeling_conditions wait for user response (see UIRESUME)
% uiwait(handles.labeling_window);


% --- Outputs from this function are returned to the command line.
function varargout = labeling_conditions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on selection change in popupmenu_label.
function popupmenu_label_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_label

handles=update_selection(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

special=get(handles.checkbox_special,'Value');
if special
    item=get(handles.popupmenu_library,'Value');
    liblist=get(handles.popupmenu_library,'String');
    handles.library=liblist{item};
end
hMain.temperature = 298;
hMain.library = handles.library;
delete(handles.labeling_window);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.library = '';
hMain.temperature = [];
delete(handles.labeling_window);


% --- Executes when user attempts to close labeling_window.
function labeling_window_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to labeling_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global hMain

hMain.library = '';
hMain.temperature = [];

delete(hObject);

function handles=initialize_popup(handles)
% initializes the popup menu with list of available labels and select
% default library

global rotamer_libraries

libs=length(rotamer_libraries);

lib_indices = zeros(1,libs);
poi = 0;
for k=1:libs
    % ### here, class must be cahnged to 'chromophore' ###
    if strcmp(rotamer_libraries(k).type,'peptide') && strcmp(rotamer_libraries(k).class,'chromophore')
        poi = poi + 1;
        lib_list{poi}=rotamer_libraries(k).label;
        lib_indices(poi) = k;
    end
end

handles.lib_indices = lib_indices(1:poi);
set(handles.popupmenu_label,'String',lib_list);
set(handles.popupmenu_label,'Value',1);


function handles=update_selection(handles)
% updates the library file name for the currently selected temperature and
% label

global rotamer_libraries

sel=get(handles.popupmenu_label,'Value');
Tvec=rotamer_libraries(handles.lib_indices(sel)).T;

% select the library calibration temperature for which the inverse 
% temperature deviates least from target temperature
Tinv=1/298;
Tinv_vec=ones(size(Tvec))./Tvec;
[~,id]=min(abs(Tinv_vec-Tinv));
handles.library=id2tag(id,rotamer_libraries(handles.lib_indices(sel)).files);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'labeling_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on selection change in popupmenu_library.
function popupmenu_library_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_library contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_library


% --- Executes during object creation, after setting all properties.
function popupmenu_library_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global general

libdir=dir(strcat(general.rootdir,'rotamer_libraries'));

poi=0;
for k=1:length(libdir)
    [~,name,ext] = fileparts(libdir(k).name);
    if strcmpi(ext,'.mat')
        s = load(libdir(k).name);
        % ### class needs to be changed to 'chromophore' here
        if isfield(s,'rot_lib') && isfield(s.rot_lib,'class')...
                && strcmp(s.rot_lib.class,'chromophore')
            poi=poi+1;
            liblist{poi}=name;
        end
    end
end
if poi==0
    liblist{1}='<No appropriate rotamer libraries found>';
end
set(hObject,'String',liblist);


% --- Executes on button press in checkbox_special.
function checkbox_special_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_special (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_special
