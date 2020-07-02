function varargout = transparency_selection(varargin)
% TRANSPARENCY_SELECTION M-file for transparency_selection.fig
%      TRANSPARENCY_SELECTION, by itself, creates a new TRANSPARENCY_SELECTION or raises the existing
%      singleton*.
%
%      H = TRANSPARENCY_SELECTION returns the handle to a new TRANSPARENCY_SELECTION or the handle to
%      the existing singleton*.
%
%      TRANSPARENCY_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRANSPARENCY_SELECTION.M with the given input arguments.
%
%      TRANSPARENCY_SELECTION('Property','Value',...) creates a new TRANSPARENCY_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before transparency_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to transparency_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help transparency_selection

% Last Modified by GUIDE v2.5 16-Sep-2009 14:00:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @transparency_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @transparency_selection_OutputFcn, ...
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


% --- Executes just before transparency_selection is made visible.
function transparency_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to transparency_selection (see VARARGIN)

% global MMM_icon
global hMain

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

% Choose default command line output for graphics_selection
handles.output = hObject;

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.alpha=1;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes transparency_selection wait for user response (see UIRESUME)
uiwait(handles.my_window);


% --- Outputs from this function are returned to the command line.
function varargout = transparency_selection_OutputFcn(hObject, eventdata, handles) 
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

entry=strcat(help_files,'transparency_selection.html');
webcall(entry,'-helpbrowser');

function edit_value_Callback(hObject, eventdata, handles)
% hObject    handle to edit_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_value as text
%        str2double(get(hObject,'String')) returns contents of edit_value as a double

global hMain

[v,handles]=edit_update_MMM(handles,hObject,0,1,1,'%4.2f',0);

hMain.alpha=v;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_value (see GCBO)
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

global hMain

pop_flag=get(handles.radiobutton_population,'Value');
if pop_flag,
    hMain.alpha=-1;
end;
close(handles.my_window);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.alpha=[];
close(handles.my_window);
