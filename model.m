function varargout = model(varargin)
% MODEL M-file for model.fig
%      MODEL, by itself, creates a new MODEL or raises the existing
%      singleton*.
%
%      H = MODEL returns the handle to a new MODEL or the handle to
%      the existing singleton*.
%
%      MODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL.M with the given input arguments.
%
%      MODEL('Property','Value',...) creates a new MODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model

% Last Modified by GUIDE v2.5 09-Mar-2009 17:14:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_OpeningFcn, ...
                   'gui_OutputFcn',  @model_OutputFcn, ...
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


% --- Executes just before model is made visible.
function model_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model (see VARARGIN)

% Choose default command line output for model
handles.output = hObject;

% declare global structure variables for all figures

global hModel

hModel.figure=handles.model_plot;
hModel.context_background_black=handles.context_background_black;
hModel.context_background_grey=handles.context_background_grey;
hModel.context_background_white=handles.context_background_white;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model wait for user response (see UIRESUME)
% uiwait(handles.model_plot);


% --- Outputs from this function are returned to the command line.
function varargout = model_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function context_background_Callback(hObject, eventdata, handles)
% hObject    handle to context_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_model_copy_Callback(hObject, eventdata, handles)
% hObject    handle to context_model_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
print(handles.model_plot,'-dbitmap');

% --------------------------------------------------------------------
function context_background_black_Callback(hObject, eventdata, handles)
% hObject    handle to context_background_black (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

hMain.color='black';

set(handles.model_plot,'Color','k');
guidata(hObject,handles);

% --------------------------------------------------------------------
function context_background_grey_Callback(hObject, eventdata, handles)
% hObject    handle to context_background_grey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

hMain.color='grey';

set(handles.model_plot,'Color',[0.941,0.941,0.941]);
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_background_white_Callback(hObject, eventdata, handles)
% hObject    handle to context_background_white (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

hMain.color='white';

set(handles.model_plot,'Color','w');
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_model_Callback(hObject, eventdata, handles)
% hObject    handle to context_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close model_plot.
function model_plot_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain

set(hMain.axes_model,'Position',hMain.oldpos);
set(hMain.axes_model,'Parent',hMain.panel_model);
set(hMain.panel_model,'Title','Model');
switch hMain.color
    case 'black'
        set(hMain.panel_model,'BackgroundColor','k');
        set(hMain.panel_model,'ForegroundColor',[216,41,0]/255);
        set(hMain.panel_model,'HighlightColor',[255,177,100]/255);
        set(hMain.panel_model,'ShadowColor',[222,125,0]/255);
    case 'grey'
        set(hMain.panel_model,'BackgroundColor',[0.941,0.941,0.941]);
        set(hMain.panel_model,'ForegroundColor','k');
        set(hMain.panel_model,'HighlightColor','w');
        set(hMain.panel_model,'ShadowColor',[128,128,128]/255);
    case 'white'
        set(hMain.panel_model,'BackgroundColor','w');
        set(hMain.panel_model,'ForegroundColor','k');
        set(hMain.panel_model,'HighlightColor','w');
        set(hMain.panel_model,'ShadowColor',[128,128,128]/255);
end;
set(hMain.edit_command_line,'String','attach');
hMain.detached=0;
delete(hObject);


% --- Executes when model_plot is resized.
function model_plot_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

newpos=get(handles.model_plot,'Position');
newpos(1)=20;
newpos(2)=20;
newpos(3)=newpos(3)-40;
newpos(4)=newpos(4)-40;
set(hMain.axes_model,'Position',newpos);
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_model_print_Callback(hObject, eventdata, handles)
% hObject    handle to context_model_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.model_plot);
guidata(hObject,handles);

% --------------------------------------------------------------------
function context_model_close_Callback(hObject, eventdata, handles)
% hObject    handle to context_model_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

model_plot_CloseRequestFcn(handles.model_plot, eventdata, handles);
