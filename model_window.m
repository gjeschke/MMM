function varargout = model_window(varargin)
% MODEL_WINDOW M-file for model_window.fig
%      MODEL_WINDOW, by itself, creates a new MODEL_WINDOW or raises the existing
%      singleton*.
%
%      H = MODEL_WINDOW returns the handle to a new MODEL_WINDOW or the handle to
%      the existing singleton*.
%
%      MODEL_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_WINDOW.M with the given input arguments.
%
%      MODEL_WINDOW('Property','Value',...) creates a new MODEL_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_window

% Last Modified by GUIDE v2.5 13-Dec-2010 17:21:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_window_OpeningFcn, ...
                   'gui_OutputFcn',  @model_window_OutputFcn, ...
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


% --- Executes just before model_window is made visible.
function model_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_window (see VARARGIN)

% Choose default command line output for model_window
handles.output = hObject;

% declare global structure variables for all figures

global hModel
% global MMM_icon
global hMain

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

hModel.figure=handles.model_plot;
hModel.context_background_black=handles.context_background_black;
hModel.context_background_grey=handles.context_background_grey;
hModel.context_background_white=handles.context_background_white;
hModel.context_select=handles.context_select;
hModel.context_rotate=handles.context_rotate;
hModel.context_zoom=handles.context_zoom;
hModel.context_pan=handles.context_pan;

set(hMain.uipushtool_copy,'Enable','on');
set(hMain.menu_file_export,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_window wait for user response (see UIRESUME)
% uiwait(handles.model_plot);


% --- Outputs from this function are returned to the command line.
function varargout = model_window_OutputFcn(hObject, eventdata, handles) 
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

global hMain

if hMain.atom_graphics_auto,
    switch_it=true;
    adjust_atom_graphics(false);
else
    switch_it=false;
end;
print(handles.model_plot,'-dbitmap');
if switch_it,
    adjust_atom_graphics(true);
end;

guidata(hObject,handles);

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
set(hMain.uipushtool_copy,'Enable','off');
set(hMain.menu_file_export,'Enable','off');
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

% newpos=get(hObject,'Position');
% newpos(1)=20;
% newpos(2)=20;
% newpos(3)=newpos(3)-40;
% newpos(4)=newpos(4)-40;
% set(hMain.axes_model,'Position',newpos);
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


% --------------------------------------------------------------------
function context_select_Callback(hObject, eventdata, handles)
% hObject    handle to context_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global hModel

if hMain.detached
    fig=hModel.figure;
else
    fig=hMain.figure;
end;

rot_state=get(hMain.uitoggletool_rotate,'State');
if strcmp(rot_state,'on'),
    set(hMain.uitoggletool_rotate,'State','off');
    view3D(fig,'off');
end;
    
zoom_state=get(hMain.uitoggletool_zoom,'State');
if strcmp(zoom_state,'on'),
    set(hMain.uitoggletool_zoom,'State','off');
    view3D(fig,'off');
end;

pan_state=get(hMain.uitoggletool_pan,'State');
if strcmp(pan_state,'on'),
    set(hMain.uitoggletool_pan,'State','off');
    view3D(fig,'off');
end;
    
set(handles.context_rotate,'Checked','off');
set(handles.context_zoom,'Checked','off');
set(handles.context_select,'Checked','on');
set(handles.context_pan,'Checked','off');

guidata(hObject,handles);

% --------------------------------------------------------------------
function context_rotate_Callback(hObject, eventdata, handles)
% hObject    handle to context_rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global hModel

if hMain.detached
    fig=hModel.figure;
else
    fig=hMain.figure;
end;

set(hMain.uitoggletool_rotate,'State','on');
set(hMain.uitoggletool_zoom,'State','off');
set(hMain.uitoggletool_select,'State','off');
set(hMain.uitoggletool_pan,'State','off');
set(handles.context_rotate,'Checked','on');
set(handles.context_zoom,'Checked','off');
set(handles.context_select,'Checked','off');
set(handles.context_pan,'Checked','off');
view3D(fig,'rot');
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function context_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to context_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global hModel

if hMain.detached
    fig=hModel.figure;
else
    fig=hMain.figure;
end;

set(hMain.uitoggletool_rotate,'State','off');
set(hMain.uitoggletool_zoom,'State','on');
set(hMain.uitoggletool_select,'State','off');
set(hMain.uitoggletool_pan,'State','off');
set(handles.context_rotate,'Checked','off');
set(handles.context_zoom,'Checked','on');
set(handles.context_select,'Checked','off');
set(handles.context_pan,'Checked','off');
view3D(fig,'zoom');
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function context_pan_Callback(hObject, eventdata, handles)
% hObject    handle to context_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global hModel

if hMain.detached
    fig=hModel.figure;
else
    fig=hMain.figure;
end;

set(hMain.uitoggletool_rotate,'State','off');
set(hMain.uitoggletool_zoom,'State','off');
set(hMain.uitoggletool_select,'State','off');
set(hMain.uitoggletool_pan,'State','on');
set(handles.context_rotate,'Checked','off');
set(handles.context_zoom,'Checked','off');
set(handles.context_select,'Checked','off');
set(handles.context_pan,'Checked','on');
view3D(fig,'pan');
    
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_color_Callback(hObject, eventdata, handles)
% hObject    handle to context_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

rgb=uisetcolor;
command_line=sprintf('color * %6.3f%6.3f%6.3f',rgb);
hMain=cmd(hMain,command_line);

guidata(hObject,handles);

% --------------------------------------------------------------------
function context_transparency_Callback(hObject, eventdata, handles)
% hObject    handle to context_transparency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

options.Resize='on';
options.WindowStyle='normal';
answer = inputdlg('Transparency value:','Set transparency for selection',1,{'1.0'},options);
command_line=sprintf('transparency * %s',answer{1});
hMain=cmd(hMain,command_line);

guidata(hObject,handles);


% --- Executes on scroll wheel click while the figure is in focus.
function model_plot_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to model_plot (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

view3DScrollFcn(hObject,eventdata);

