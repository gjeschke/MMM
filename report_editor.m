function varargout = report_editor(varargin)
% REPORT_EDITOR M-file for report_editor.fig
%      REPORT_EDITOR, by itself, creates a new REPORT_EDITOR or raises the existing
%      singleton*.
%
%      H = REPORT_EDITOR returns the handle to a new REPORT_EDITOR or the handle to
%      the existing singleton*.
%
%      REPORT_EDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REPORT_EDITOR.M with the given input arguments.
%
%      REPORT_EDITOR('Property','Value',...) creates a new REPORT_EDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before report_editor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to report_editor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help report_editor

% Last Modified by GUIDE v2.5 06-Jan-2010 08:06:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @report_editor_OpeningFcn, ...
                   'gui_OutputFcn',  @report_editor_OutputFcn, ...
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


% --- Executes just before report_editor is made visible.
function report_editor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to report_editor (see VARARGIN)

% Choose default command line output for report_editor
handles.output = hObject;

global MMM_icon
global hMain

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

handles.edited=0;
handles.read_only=0;

if ~isempty(hMain.report_file),
    handles.open_file=hMain.report_file;
    handles=load_file(handles);
    hMain.report_file='';
end;

hMain.auxiliary=[hMain.auxiliary hObject];

set(hObject,'Pointer','arrow');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes report_editor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = report_editor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_report_Callback(hObject, eventdata, handles)
% hObject    handle to edit_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_report as text
%        str2double(get(hObject,'String')) returns contents of edit_report as a double

handles.edited=1;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_report_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'report_editor.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.reports);

[fname,pname]=uigetfile({'*.txt';'*.dat';'*.pdb';'*.ent';'*.*'},'Open report text file');
if ~isequal(fname,0) && ~isequal(pname,0)
    reset_user_paths(pname);
    general.reports=pname;
    handles.open_file=fullfile(pname,fname);
    handles=load_file(handles);
end;

if handles.read_only,
    set(handles.pushbutton_save,'Enable','off');
else
    set(handles.pushbutton_save,'Enable','on');
end;
cd(my_path);
guidata(hObject,handles);

function handles=load_file(handles)

fid=fopen(handles.open_file,'r');
if fid==-1,
    handles.open_file='';
    set(handles.edit_report,'String','<file could not be opened>');
    set(handles.figure1,'Name','Report Editor - <no file>');
    set(handles.text_warning,'String','No file.');
else
    handles.read_only=0;
    set(handles.text_warning,'String','OK');
    set(handles.figure1,'Name',sprintf('Report Editor - %s',handles.open_file));
    set(handles.edit_report,'String','Loading...');
    myfig=gcf;
    set(myfig,'Pointer','watch');
    drawnow;
    poi=0;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if length(tline)>80 && ~handles.read_only,
            set(handles.text_warning,'String','More than 80 characters per line. Read only.');
            set(handles.figure1,'Name',sprintf('Report Editor - %s <read only>',handles.open_file));
            handles.read_only=1;
        end;
        poi=poi+1;
        text{poi}=tline;
    end;
    fclose(fid);
    if poi>0,
        set(handles.edit_report,'String',text);
    else
        set(handles.edit_report,'String','<empty file>');
    end;
    set(myfig,'Pointer','arrow');
    drawnow;
end;


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;

if handles.read_only,
    set(handles.text_warning,'String','ERROR: File cannot be saved (read only).');
    guidata(hObject,handles);
    return;
end;

cd(general.reports);
[filename, pathname] = uiputfile({'*.txt';'*.dat';'*.pdb';'*.ent';'*.*'}, 'Save report as...');
if isequal(filename,0) || isequal(pathname,0)
    set(handles.text_warning,'String','WARNING: Report not saved (cancelled).');
else
    reset_user_paths(pathname);
    general.reports=pathname;
    fname=fullfile(pathname, filename);
    set(handles.text_warning,'String',sprintf('Report saved: %s',fname));
    fid=fopen(fname,'wt');
    text=get(handles.edit_report,'String');
    if iscell(text),
        for k=1:numel(text),
            fprintf(fid,'%s\n',text{k});
        end;
    else
        [nl,n]=size(text);
        for k=1:nl
            fprintf(fid,'%s\n',text(k,:));
        end;
    end;
    fclose(fid);
end
cd(my_path);

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.edited && ~handles.read_only,
    cancel=ask_save(handles);
    if cancel,
        set(handles.text_warning,'String','Closing of the window cancelled.');
        guidata(hObject,handles);
        return;
    end;
end;

delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if handles.edited && ~handles.read_only,
    cancel=ask_save(handles);
end;

delete(hObject);


function cancel=ask_save(handles)

global general

my_path=pwd;

cancel=0;

button = questdlg('Do you want to save changes?','Report was possibly edited','Yes');
switch button
    case 'Yes'
        cd(general.reports);
        [filename, pathname] = uiputfile({'*.txt';'*.dat';'*.pdb';'*.ent';'*.*'}, 'Save report as...');
        cd(my_path);
        if isequal(filename,0) || isequal(pathname,0)
            cancel=1;
        else
            reset_user_paths(pname);
            general.reports=pname;
            fname=fullfile(pathname, filename);
            fid=fopen(fname,'wt');
            text=get(handles.edit_report,'String');
            if iscell(text),
                for k=1:numel(text),
                    fprintf(fid,'%s\n',text{k});
                end;
            else
                [nl,n]=size(text);
                for k=1:nl
                    fprintf(fid,'%s\n',text(k,:));
                end;
            end;
            fclose(fid);
        end
    case 'Cancel'
        cancel=1;
end;
