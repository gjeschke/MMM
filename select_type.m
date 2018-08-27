function varargout = select_type(varargin)
% SELECT_TYPE M-file for select_type.fig
%      SELECT_TYPE, by itself, creates a new SELECT_TYPE or raises the existing
%      singleton*.
%
%      H = SELECT_TYPE returns the handle to a new SELECT_TYPE or the handle to
%      the existing singleton*.
%
%      SELECT_TYPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_TYPE.M with the given input arguments.
%
%      SELECT_TYPE('Property','Value',...) creates a new SELECT_TYPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_type_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_type_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_type

% Last Modified by GUIDE v2.5 02-Dec-2010 10:45:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_type_OpeningFcn, ...
                   'gui_OutputFcn',  @select_type_OutputFcn, ...
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


% --- Executes just before select_type is made visible.
function select_type_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_type (see VARARGIN)

% Choose default command line output for select_type
handles.output = hObject;

indices=resolve_address('*');

[m,n]=size(indices);
indices=indices(indices>0);
if m==1,
    [m,n]=size(indices);
end;
if m~=1 || n~=4,
    add_msg_board('ERROR: A single residue must be selected for mutation.');
    delete(hObject);
    return
else
    address=mk_address(indices,true);
    set(handles.text_mutation,'String',sprintf('Residue %s is to be mutated.',address));
end;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_type wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_type_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global exchange_container

hrb=get(handles.uipanel_residues,'SelectedObject');
exchange_container=get(hrb,'String');
delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global exchange_container

exchange_container='';
delete(handles.figure1);
