function varargout = nonstandard_replacement(varargin)
% NONSTANDARD_REPLACEMENT M-file for nonstandard_replacement.fig
%      NONSTANDARD_REPLACEMENT, by itself, creates a new NONSTANDARD_REPLACEMENT or raises the existing
%      singleton*.
%
%      H = NONSTANDARD_REPLACEMENT returns the handle to a new NONSTANDARD_REPLACEMENT or the handle to
%      the existing singleton*.
%
%      NONSTANDARD_REPLACEMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NONSTANDARD_REPLACEMENT.M with the given input arguments.
%
%      NONSTANDARD_REPLACEMENT('Property','Value',...) creates a new NONSTANDARD_REPLACEMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nonstandard_replacement_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nonstandard_replacement_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nonstandard_replacement

% Last Modified by GUIDE v2.5 30-Nov-2010 11:44:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nonstandard_replacement_OpeningFcn, ...
                   'gui_OutputFcn',  @nonstandard_replacement_OutputFcn, ...
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


% --- Executes just before nonstandard_replacement is made visible.
function nonstandard_replacement_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nonstandard_replacement (see VARARGIN)

% Choose default command line output for nonstandard_replacement
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nonstandard_replacement wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nonstandard_replacement_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on button press in checkbox_CSE.
function checkbox_CSE_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CSE


% --- Executes on button press in checkbox_MSE.
function checkbox_MSE_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MSE


% --- Executes on button press in checkbox_R1A.
function checkbox_R1A_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_R1A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_R1A


% --- Executes on button press in checkbox_R1B.
function checkbox_R1B_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_R1B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_R1B


% --- Executes on button press in checkbox_R1F.
function checkbox_R1F_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_R1F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_R1F


% --- Executes on button press in checkbox_R7A.
function checkbox_R7A_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_R7A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_R7A


% --- Executes on button press in checkbox_V1A.
function checkbox_V1A_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_V1A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_V1A


% --- Executes on button press in checkbox_IA1.
function checkbox_IA1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_IA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_IA1


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global replacements

if get(handles.radiobutton_selected,'Value'),
    replacements='S:';
else
    replacements='A:';
end;

if get(handles.checkbox_CSE,'Value'),
    replacements=strcat(replacements,'CSE:');
end;
if get(handles.checkbox_MSE,'Value'),
    replacements=strcat(replacements,'MSE:');
end;
if get(handles.checkbox_R1A,'Value'),
    replacements=strcat(replacements,'R1A:');
end;
if get(handles.checkbox_R1B,'Value'),
    replacements=strcat(replacements,'R1B:');
end;
if get(handles.checkbox_R1F,'Value'),
    replacements=strcat(replacements,'R1F:');
end;
if get(handles.checkbox_R7A,'Value'),
    replacements=strcat(replacements,'R7A:');
end;
if get(handles.checkbox_V1A,'Value'),
    replacements=strcat(replacements,'V1A:');
end;
if get(handles.checkbox_IA1,'Value'),
    replacements=strcat(replacements,'IA1:');
end;

delete(handles.figure1);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global replacements

replacements='';

delete(handles.figure1);
