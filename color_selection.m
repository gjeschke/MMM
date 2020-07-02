function varargout = color_selection(varargin)
% COLOR_SELECTION M-file for color_selection.fig
%      COLOR_SELECTION, by itself, creates a new COLOR_SELECTION or raises the existing
%      singleton*.
%
%      H = COLOR_SELECTION returns the handle to a new COLOR_SELECTION or the handle to
%      the existing singleton*.
%
%      COLOR_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLOR_SELECTION.M with the given input arguments.
%
%      COLOR_SELECTION('Property','Value',...) creates a new COLOR_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before color_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to color_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help color_selection

% Last Modified by GUIDE v2.5 16-Sep-2009 14:59:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @color_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @color_selection_OutputFcn, ...
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


% --- Executes just before color_selection is made visible.
function color_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to color_selection (see VARARGIN)

% global MMM_icon
global hMain
global graph_settings

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

% Choose default command line output for graphics_selection
handles.output = hObject;

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.color_selection=[0,0,1];
handles.palette_selection=[0,0,1];
set(handles.pushbutton_palette,'ForegroundColor',handles.palette_selection);
handles.name_selection=10;
handles.named_color=graph_settings.colors(handles.name_selection,:);
set(handles.text_named_color,'ForegroundColor',handles.named_color);

set(handles.listbox_colors,'String',graph_settings.color_names);
set(handles.listbox_colors,'Value',handles.name_selection);



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes color_selection wait for user response (see UIRESUME)
uiwait(handles.my_window);


% --- Outputs from this function are returned to the command line.
function varargout = color_selection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'color_selection.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global graph_settings

if get(handles.radiobutton_name,'Value'),
    sel=get(handles.listbox_colors,'Value');
    hMain.color_selection=graph_settings.colors(sel,:);
end;

if get(handles.radiobutton_palette,'Value'),
    hMain.color_selection=handles.palette_selection;
end;

if get(handles.radiobutton_secondary,'Value'),
    hMain.color_selection='secondary';
end;

if get(handles.radiobutton_chain,'Value'),
    hMain.color_selection='chain';
end;

if get(handles.radiobutton_sequence,'Value'),
    hMain.color_selection='sequence';
end;

if get(handles.radiobutton_Bfactor,'Value'),
    hMain.color_selection='Bfactor';
end;

if get(handles.radiobutton_Bfactor_tight,'Value'),
    hMain.color_selection='Bfactor_tight';
end;

if get(handles.radiobutton_charge,'Value'),
    hMain.color_selection='charge';
end;

if get(handles.radiobutton_hydropathy,'Value'),
    hMain.color_selection='hydropathy';
end;

if get(handles.radiobutton_helix_propensity,'Value'),
    hMain.color_selection='helix_propensity';
end;

if get(handles.radiobutton_ensemble,'Value'),
    hMain.color_selection='ensemble';
end;

close(handles.my_window);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.color_selection=[];

close(handles.my_window);

% --- Executes on selection change in listbox_colors.
function listbox_colors_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_colors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_colors

global graph_settings

sel=get(hObject,'Value');
rgb=graph_settings.colors(sel,:);
set(handles.text_named_color,'ForegroundColor',rgb);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_colors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_palette.
function pushbutton_palette_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_palette (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global graph_settings

rgb=uisetcolor;
handles.palette_selection=rgb;
set(hObject,'ForegroundColor',rgb);
poi=0;
mindist=1e6;
poi2=0;
mindist2=1e6;
[m,n]=size(graph_settings.colors);
for k=1:m,
    dist=norm(rgb-graph_settings.colors(k,:));
    if dist<mindist, poi=k; mindist=dist; end;
end;
set(handles.listbox_colors,'Value',poi);
set(handles.text_named_color,'ForegroundColor',graph_settings.colors(poi,:));

% Update handles structure
guidata(hObject, handles);
