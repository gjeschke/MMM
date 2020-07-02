function varargout = inspect_ensemble_fit(varargin)
% INSPECT_ENSEMBLE_FIT M-file for inspect_ensemble_fit.fig
%      INSPECT_ENSEMBLE_FIT, by itself, creates a new INSPECT_ENSEMBLE_FIT or raises the existing
%      singleton*.
%
%      H = INSPECT_ENSEMBLE_FIT returns the handle to a new INSPECT_ENSEMBLE_FIT or the handle to
%      the existing singleton*.
%
%      INSPECT_ENSEMBLE_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSPECT_ENSEMBLE_FIT.M with the given input arguments.
%
%      INSPECT_ENSEMBLE_FIT('Property','Value',...) creates a new INSPECT_ENSEMBLE_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inspect_ensemble_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inspect_ensemble_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inspect_ensemble_fit

% Last Modified by GUIDE v2.5 12-Mar-2012 10:10:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inspect_ensemble_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @inspect_ensemble_fit_OutputFcn, ...
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


% --- Executes just before inspect_ensemble_fit is made visible.
function inspect_ensemble_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inspect_ensemble_fit (see VARARGIN)

global hMain
% global MMM_icon

% Choose default command line output for inspect_ensemble_fit
handles.output = hObject;

if isempty(hMain.fit_diagnostics{1})
    add_msg_board('Warning: No diagnostics information available.');
    delete(hObject);
    return;
end

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

hMain.auxiliary=[hMain.auxiliary hObject];

set(handles.figure1,'Name','Inspect ensemble fit diagnostics');

handles.trial=1;

update(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inspect_ensemble_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inspect_ensemble_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

maxt=length(hMain.fit_diagnostics);
if handles.trial<maxt,
    handles.trial=handles.trial+1;
    update(hObject,handles);
end;


% --- Executes on button press in pushbutton_minus.
function pushbutton_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.trial>1,
    handles.trial=handles.trial-1;
    update(hObject,handles);
end;


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

function update(hObject,handles)

global hMain

k=handles.trial;

maxt=length(hMain.fit_diagnostics);

x=[0 hMain.fit_diagnostics{k}.x];
y=hMain.fit_diagnostics{k}.y;
n=round(length(x)/5);

p1=polyfit(x(1:n),y(1:n),1);
p2=polyfit(x(end-n+1:end),y(end-n+1:end),1);
xp=round((p2(2)-p1(2))/(p1(1)-p2(1)));
y1=polyval(p1,x);
y2=polyval(p2,x);
if xp<2, xp=2; end;
if xp>length(hMain.fit_diagnostics{k}.ys),
    xp=length(hMain.fit_diagnostics{k}.ys);
end;

axes(handles.axes1); cla;
plot([0 hMain.fit_diagnostics{k}.x],hMain.fit_diagnostics{k}.y,'k');
hold on;
plot(x(1:xp),y1(1:xp),'r:');
if ~isempty(y2) && xp>=1 && xp<= length(y2),
    plot(x(xp:end),y2(xp:end),'r:');
    plot(hMain.fit_diagnostics{k}.converged,hMain.fit_diagnostics{k}.y(hMain.fit_diagnostics{k}.converged),'k^');
    plot(xp-1,hMain.fit_diagnostics{k}.y(xp),'ro');
end;
xlabel('Iteration number');
ylabel('r.m.s.d. (?)');
axis tight

if ~isempty(hMain.fit_diagnostics{k}.ys),
    axes(handles.axes3); cla;
    plot(hMain.fit_diagnostics{k}.x,hMain.fit_diagnostics{k}.ys,'b');
    hold on;
    plot(hMain.fit_diagnostics{k}.converged,hMain.fit_diagnostics{k}.ys(hMain.fit_diagnostics{k}.converged-1),'k^');
    plot(xp-1,hMain.fit_diagnostics{k}.ys(xp),'ro');
    xlabel('Iteration number');
    ylabel('r.m.s.d. (?)');
    axis tight
end;

msg=sprintf('Trial %i/%i',k,maxt);
set(handles.text_trial,'String',msg);
messages{1}=sprintf('Converged at iteration %i',hMain.fit_diagnostics{k}.converged);
messages{2}=sprintf('Initial constraint r.m.s.d. %5.2f ?',hMain.fit_diagnostics{k}.drmsd0);
messages{3}=sprintf('Final constraint r.m.s.d.   %5.2f ?',hMain.fit_diagnostics{k}.drmsd);
if ~isempty(hMain.fit_diagnostics{k}.rmsd0) && ~isempty(hMain.fit_diagnostics{k}.rmsd0),
    messages{4}=sprintf('Initial r.m.s.d. to target structure %5.2f ?',hMain.fit_diagnostics{k}.rmsd0);
    messages{5}=sprintf('Final r.m.s.d. to target structure   %5.2f ?',hMain.fit_diagnostics{k}.rmsd);
    messages{6}=sprintf('Fractional coverage %5.3f',(hMain.fit_diagnostics{k}.rmsd0-hMain.fit_diagnostics{k}.rmsd)/hMain.fit_diagnostics{k}.rmsd0);
end;
set(handles.text_info,'String',messages);

guidata(hObject,handles);
