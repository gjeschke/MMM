function varargout = labeling_conditions_nucleotides(varargin)
% LABELING_CONDITIONS_NUCLEOTIDES M-file for labeling_conditions_nucleotides.fig
%      LABELING_CONDITIONS_NUCLEOTIDES, by itself, creates a new LABELING_CONDITIONS_NUCLEOTIDES or raises the existing
%      singleton*.
%
%      H = LABELING_CONDITIONS_NUCLEOTIDES returns the handle to a new LABELING_CONDITIONS_NUCLEOTIDES or the handle to
%      the existing singleton*.
%
%      LABELING_CONDITIONS_NUCLEOTIDES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABELING_CONDITIONS_NUCLEOTIDES.M with the given input arguments.
%
%      LABELING_CONDITIONS_NUCLEOTIDES('Property','Value',...) creates a new LABELING_CONDITIONS_NUCLEOTIDES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before labeling_conditions_nucleotides_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to labeling_conditions_nucleotides_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help labeling_conditions_nucleotides

% Last Modified by GUIDE v2.5 15-Jan-2014 18:29:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labeling_conditions_nucleotides_OpeningFcn, ...
                   'gui_OutputFcn',  @labeling_conditions_nucleotides_OutputFcn, ...
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


% --- Executes just before labeling_conditions_nucleotides is made visible.
function labeling_conditions_nucleotides_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labeling_conditions_nucleotides (see VARARGIN)

% Choose default command line output for labeling_conditions_nucleotides

global hMain
global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;

handles.T = 298; % default ambient temperature
handles=initialize_popup(handles);
handles=update_selection(handles);

if isfield(hMain,'label_selection'),
    title=sprintf('Set labeling conditions for %s',hMain.label_selection);
else
    title='Set labeling conditions for nucleotide residues';
end;
set(hObject,'Name',title);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

locked=get(hMain.uitoggletool_lock,'State');
if strcmpi(locked,'on'),
    set(handles.checkbox_special,'Enable','off');
    set(handles.popupmenu_library,'Enable','off');
else
    set(handles.checkbox_special,'Enable','on');
    set(handles.popupmenu_library,'Enable','on');
end;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes labeling_conditions_nucleotides wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = labeling_conditions_nucleotides_OutputFcn(hObject, eventdata, handles) 
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

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_label contents as cell array
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
if special,
    item=get(handles.popupmenu_library,'Value');
    liblist=get(handles.popupmenu_library,'String');
    handles.library=liblist{item};
end;
hMain.temperature=handles.T;
hMain.library=handles.library;
delete(handles.figure1);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.library='';
hMain.temperature=[];
delete(handles.figure1);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'labeling_nucleotides_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on button press in checkbox_special.
function checkbox_special_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_special (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_special


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

global general
global hMain

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

pattern = hMain.residue_pattern;

libdir=dir(strcat(general.rootdir,'rotamer_libraries'));

poi=0;
for k=1:length(libdir),
    [path,name,ext] = fileparts(libdir(k).name);
    if strcmpi(ext,'.mat'),
        load(libdir(k).name);
        if isfield(rot_lib,'attachment'),
            attachment = rot_lib.attachment;
        else
            attachment = 'peptide';
        end;
        relevant = false;
        isA = ~isempty(tag2id('A',attachment));
        isDA = ~isempty(tag2id('DA',attachment));
        isC = ~isempty(tag2id('C',attachment));
        isDC = ~isempty(tag2id('DC',attachment));
        isG = ~isempty(tag2id('G',attachment));
        isDG = ~isempty(tag2id('DG',attachment));
        isU = ~isempty(tag2id('U',attachment));
        isDT = ~isempty(tag2id('DT',attachment));
        for kp = 1:length(pattern),
            switch pattern(kp)
                case 'A'
                    if isA && isDA, relevant = true; end;
                case 'C'
                    if isC && isDC, relevant = true; end;
                case 'G'
                    if isG && isDG, relevant = true; end;
                case 'T'
                    if isDT, relevant = true; end;
                case 'U'
                    if isU, relevant = true; end;
                case 'D'
                    if isDA && isDC && isDG && isDT, relevant = true; end;
                case 'R'
                    if isA && isC && isG && isU, relevant = true; end;
            end;
        end;
        if relevant,
            poi=poi+1;
            liblist{poi}=name;
        end;
    end;
end;
if poi==0,
    liblist{1}='<No relevant rotamer libraries found>';
end;
set(hObject,'String',liblist);

function edit_temperature_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_temperature as text
%        str2double(get(hObject,'String')) returns contents of edit_temperature as a double

[v,handles]=edit_update_MMM(handles,hObject,4,500,175,'%3.0f',0);

handles.T=v;
handles=update_selection(handles);
guidata(hObject,handles);

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

function handles=initialize_popup(handles)
% initializes the popup menu with list of available labels and select
% default library

global rotamer_libraries

libs=length(rotamer_libraries);

lib_list{1} = '<No nucleotide label found>';

lib_indices = zeros(1,libs);
poi = 0;
for k=1:libs,
    if strcmp(rotamer_libraries(k).type,'nucleotide')
        poi = poi + 1;
        lib_list{poi}=rotamer_libraries(k).label;
        lib_indices(poi) = k;
    end;
end;

handles.lib_indices = lib_indices(1:poi);
set(handles.popupmenu_label,'String',lib_list);
set(handles.popupmenu_label,'Value',1);



function handles=update_selection(handles)
% updates the library file name for the currently selected temperature and
% label

global rotamer_libraries

sel=get(handles.popupmenu_label,'Value');
Tvec=rotamer_libraries(handles.lib_indices(sel)).T;

% select the library calibartion temperature for which the inverse 
% temperature deviates least from target temperature
Tinv=1/handles.T;
Tinv_vec=ones(size(Tvec))./Tvec;
[~,id]=min(abs(Tinv_vec-Tinv));
handles.library=id2tag(id,rotamer_libraries(handles.lib_indices(sel)).files);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global hMain

hMain.library='';
hMain.temperature=[];

delete(hObject);


% --- Executes on button press in radiobutton_cryogenic.
function radiobutton_cryogenic_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_cryogenic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_cryogenic

cryo=get(hObject,'Value');
if cryo,
    set(handles.edit_temperature,'String','175');
    handles.T=175;
    set(handles.radiobutton_ambient,'Value',0);
else
    set(handles.radiobutton_ambient,'Value',1);
end;
handles=update_selection(handles);
guidata(hObject,handles);

% --- Executes on button press in radiobutton_ambient.
function radiobutton_ambient_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_ambient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_ambient

ambient=get(hObject,'Value');
if ambient,
    set(handles.radiobutton_cryogenic,'Value',0);
    set(handles.edit_temperature,'String','298');
    handles.T=298;
else
    set(handles.radiobutton_cryogenic,'Value',1);
end;
handles=update_selection(handles);
guidata(hObject,handles);
