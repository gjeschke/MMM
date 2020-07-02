function varargout = site_scan_nucleotides(varargin)
% SITE_SCAN_NUCLEOTIDES M-file for site_scan_nucleotides.fig
%      SITE_SCAN_NUCLEOTIDES, by itself, creates a new SITE_SCAN_NUCLEOTIDES or raises the existing
%      singleton*.
%
%      H = SITE_SCAN_NUCLEOTIDES returns the handle to a new SITE_SCAN_NUCLEOTIDES or the handle to
%      the existing singleton*.
%
%      SITE_SCAN_NUCLEOTIDES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SITE_SCAN_NUCLEOTIDES.M with the given input arguments.
%
%      SITE_SCAN_NUCLEOTIDES('Property','Value',...) creates a new SITE_SCAN_NUCLEOTIDES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before site_scan_nucleotides_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to site_scan_nucleotides_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help site_scan_nucleotides

% Last Modified by GUIDE v2.5 15-Jan-2014 15:19:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @site_scan_nucleotides_OpeningFcn, ...
                   'gui_OutputFcn',  @site_scan_nucleotides_OutputFcn, ...
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


% --- Executes just before site_scan_nucleotides is made visible.
function site_scan_nucleotides_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to site_scan_nucleotides (see VARARGIN)

global hMain
% global MMM_icon

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

% Choose default command line output for site_scan
handles.output = hObject;

hMain.site_scan = 0;
hMain.site_scan_homooligomer=0;
hMain.site_scan_intra=1;
hMain.site_scan_inter=2;
hMain.site_scan_multiplicity=2;
hMain.site_scan_type = 'nucleotide';

handles.multiplicity=2;

if ~hMain.site_scan_residue,
    set(handles.checkbox_no_rot_populations,'Enable','off');
end;

setup_residue_pattern(hObject,handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes site_scan_nucleotides wait for user response (see UIRESUME)
uiwait(handles.site_scan_nuc);


% --- Outputs from this function are returned to the command line.
function varargout = site_scan_nucleotides_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.site_scan=1;
hMain.z_analysis = get(handles.checkbox_z_analysis,'Value');
hMain.statistics = get(handles.checkbox_save_statistics,'Value');
hMain.rotamer_PDB = get(handles.checkbox_PDB_rotamers,'Value');
hMain.no_rot_pop = get(handles.checkbox_no_rot_populations,'Value');
hMain.dynamic_rotamers = 0;
residue_pattern(handles);

delete(handles.site_scan_nuc);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.site_scan=0;

delete(handles.site_scan_nuc);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'site_scan_nucleotides_window.html');
web(entry,'-helpbrowser');


% --- Executes on button press in checkbox_no_rot_populations.
function checkbox_no_rot_populations_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_no_rot_populations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_no_rot_populations

global hMain
no_rot_pop=get(hObject,'Value');
hMain.no_rot_pop=no_rot_pop;

% --- Executes on button press in checkbox_save_statistics.
function checkbox_save_statistics_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_save_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_save_statistics

if ~get(hObject,'Value'),
    set(handles.checkbox_PDB_rotamers,'Value',0);
end;
guidata(hObject,handles);


% --- Executes on button press in checkbox_PDB_rotamers.
function checkbox_PDB_rotamers_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PDB_rotamers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PDB_rotamers

if get(hObject,'Value'),
    ButtonName = questdlg('This is only recommended for methodological work. Do you really need these files?', 'Very storage intensive feature', 'No', 'Yes', 'No');
    if strcmpi(ButtonName,'Yes'),
        set(handles.checkbox_save_statistics,'Value',1);
    else
        set(hObject,'Value',0);
    end;
end;
guidata(hObject,handles);

% --- Executes on button press in radiobutton_intra_all.
function radiobutton_intra_all_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_intra_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_intra_all

global hMain

hMain.site_scan_intra=1;
set(handles.radiobutton_intra_none,'Value',0);
set(handles.radiobutton_intra_all,'Value',1);
guidata(hObject, handles);

% --- Executes on button press in checkbox_homoligomer.
function checkbox_homoligomer_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_homoligomer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_homoligomer

global hMain

hom=get(hObject,'Value');
hMain.site_scan_homooligomer=hom;
if hom,
    hMain.site_scan_multiplicity=handles.multiplicity;
    set(handles.edit_multiplicity,'String',sprintf('%i',handles.multiplicity));
    set(handles.edit_multiplicity,'Enable','on');
    set(handles.text_multiplicity,'Enable','on');
else
    set(handles.edit_multiplicity,'Enable','off');
    set(handles.text_multiplicity,'Enable','off');
end;
guidata(hObject, handles);

function edit_multiplicity_Callback(hObject, eventdata, handles)
% hObject    handle to edit_multiplicity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_multiplicity as text
%        str2double(get(hObject,'String')) returns contents of edit_multiplicity as a double

global hMain

[v,handles]=edit_update_MMM(handles,hObject,2,100,3,'%i',1);
handles.multiplicity=v;
hMain.site_scan_multiplicity=handles.multiplicity;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_multiplicity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_multiplicity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_inter_all.
function radiobutton_inter_all_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_inter_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_inter_all

global hMain

hMain.site_scan_inter=2;
set(handles.radiobutton_inter_equiv,'Value',0);
set(handles.radiobutton_inter_none,'Value',0);
set(handles.radiobutton_inter_all,'Value',1);
guidata(hObject, handles);

% --- Executes on button press in radiobutton_inter_equiv.
function radiobutton_inter_equiv_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_inter_equiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_inter_equiv

global hMain

equiv=get(hObject,'Value');
if equiv,
    hMain.site_scan_inter=1;
    set(handles.radiobutton_inter_none,'Value',0);
    set(handles.radiobutton_inter_all,'Value',0);
else
    hMain.site_scan_inter=2;
    set(handles.radiobutton_inter_equiv,'Value',0);
    set(handles.radiobutton_inter_all,'Value',1);
end;    
guidata(hObject, handles);

% --- Executes on button press in radiobutton_inter_none.
function radiobutton_inter_none_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_inter_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_inter_none

global hMain

none=get(hObject,'Value');
if none,
    hMain.site_scan_inter=0;
    set(handles.radiobutton_inter_equiv,'Value',0);
    set(handles.radiobutton_inter_all,'Value',0);
else
    hMain.site_scan_inter=2;
    set(handles.radiobutton_inter_equiv,'Value',0);
    set(handles.radiobutton_inter_all,'Value',1);
end;    
guidata(hObject, handles);

% --- Executes on button press in radiobutton_intra_none.
function radiobutton_intra_none_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_intra_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_intra_none

global hMain

none=get(hObject,'Value');
if none,
    hMain.site_scan_intra=0;
    set(handles.radiobutton_intra_all,'Value',0);
else
    hMain.site_scan_intra=1;
    set(handles.radiobutton_intra_all,'Value',1);
end;    
guidata(hObject, handles);



% --- Executes on button press in checkbox_z_analysis.
function checkbox_z_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_z_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_z_analysis


% --- Executes on button press in checkbox_adenine.
function checkbox_adenine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_adenine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_adenine


% --- Executes on button press in checkbox_cytosine.
function checkbox_cytosine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cytosine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cytosine


% --- Executes on button press in checkbox_guanine.
function checkbox_guanine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_guanine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_guanine


% --- Executes on button press in checkbox_thymine.
function checkbox_thymine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_thymine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_thymine


% --- Executes on button press in checkbox_uracil.
function checkbox_uracil_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_uracil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_uracil


% --- Executes on button press in checkbox_desoxyribose.
function checkbox_desoxyribose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_desoxyribose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_desoxyribose


% --- Executes on button press in checkbox_ribose.
function checkbox_ribose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ribose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ribose

function setup_residue_pattern(hObject,handles)

global label_defs

set(handles.checkbox_adenine,'Enable','off');
set(handles.checkbox_cytosine,'Enable','off');
set(handles.checkbox_guanine,'Enable','off');
set(handles.checkbox_thymine,'Enable','off');
set(handles.checkbox_uracil,'Enable','off');
set(handles.checkbox_desoxyribose,'Enable','off');
set(handles.checkbox_ribose,'Enable','off');
enabled = '';

for k = 1:length(label_defs.residues),
    attachment = label_defs.residues(k).attachment;
    isA = ~isempty(tag2id('A',attachment));
    isDA = ~isempty(tag2id('DA',attachment));
    if isA && isDA,
        set(handles.checkbox_adenine,'Enable','on');
        if isempty(strfind(enabled,'A')),
            enabled = [enabled 'A'];
        end;
    end;
    isC = ~isempty(tag2id('C',attachment));
    isDC = ~isempty(tag2id('DC',attachment));
    if isC && isDC,
        set(handles.checkbox_cytosine,'Enable','on');
        if isempty(strfind(enabled,'C')),
            enabled = [enabled 'C'];
        end;
    end;
    isG = ~isempty(tag2id('G',attachment));
    isDG = ~isempty(tag2id('DG',attachment));
    if isG && isDG,
        set(handles.checkbox_guanine,'Enable','on');
        if isempty(strfind(enabled,'G')),
            enabled = [enabled 'G'];
        end;
    end;
    isU = ~isempty(tag2id('U',attachment));
    isDT = ~isempty(tag2id('DT',attachment));
    if isDT,
        set(handles.checkbox_thymine,'Enable','on');
        if isempty(strfind(enabled,'T')),
            enabled = [enabled 'T'];
        end;
    end;
    if isU,
        set(handles.checkbox_uracil,'Enable','on');
        if isempty(strfind(enabled,'U')),
            enabled = [enabled 'U'];
        end;
    end;
    if isDA && isDC && isDG && isDT,
        set(handles.checkbox_desoxyribose,'Enable','on');
        if isempty(strfind(enabled,'D')),
            enabled = [enabled 'D'];
        end;
    end;
    if isA && isC && isG && isU,
        set(handles.checkbox_ribose,'Enable','on');
        if isempty(strfind(enabled,'R')),
            enabled = [enabled 'R'];
        end;
    end;
end;

handles.enabled = enabled;

if isempty(enabled),
    set(handles.pushbutton_OK,'Enable','off');
else
    set(handles.pushbutton_OK,'Enable','on');
    for k = 1:length(enabled),
        switch enabled(k),
            case 'A'
                set(handles.checkbox_adenine,'Value',1);
            case 'C'
                set(handles.checkbox_cytosine,'Value',1);
            case 'G'
                set(handles.checkbox_guanine,'Value',1);
            case 'T'
                set(handles.checkbox_thymine,'Value',1);
            case 'U'
                set(handles.checkbox_uracil,'Value',1);
            case 'D'
                set(handles.checkbox_desoxyribose,'Value',1);
            case 'R'
                set(handles.checkbox_ribose,'Value',1);
        end;
    end;
end;

guidata(hObject, handles);

function pattern = residue_pattern(handles)
% store array with residue pattern in hMain and return it
% the pattern is a vector of flags for natural amino acid residues ordered
% alphabetically according to the three-letter code 

global hMain

pattern = '';

if get(handles.checkbox_adenine,'Value')
    pattern = [pattern 'A'];
end;
if get(handles.checkbox_cytosine,'Value')
    pattern = [pattern 'C'];
end;
if get(handles.checkbox_guanine,'Value')
    pattern = [pattern 'G'];
end;
if get(handles.checkbox_thymine,'Value')
    pattern = [pattern 'T'];
end;
if get(handles.checkbox_uracil,'Value')
    pattern = [pattern 'U'];
end;
if get(handles.checkbox_desoxyribose,'Value')
    pattern = [pattern 'D'];
end;
if get(handles.checkbox_ribose,'Value')
    pattern = [pattern 'R'];
end;

hMain.residue_pattern = pattern;
