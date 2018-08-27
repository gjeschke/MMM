function varargout = SAXS_fit(varargin)
% SAXS_FIT MATLAB code for SAXS_fit.fig
%      SAXS_FIT, by itself, creates a new SAXS_FIT or raises the existing
%      singleton*.
%
%      H = SAXS_FIT returns the handle to a new SAXS_FIT or the handle to
%      the existing singleton*.
%
%      SAXS_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAXS_FIT.M with the given input arguments.
%
%      SAXS_FIT('Property','Value',...) creates a new SAXS_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SAXS_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SAXS_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SAXS_fit

% Last Modified by GUIDE v2.5 19-Apr-2017 13:07:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SAXS_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @SAXS_fit_OutputFcn, ...
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


% --- Executes just before SAXS_fit is made visible.
function SAXS_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SAXS_fit (see VARARGIN)

global MMM_icon
global hMain
global model

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

handles.SAXS_curve = [];
handles.SAXS_datafile = '';
handles.SAXS_path = '';
handles.SAXS_illres = '';
handles.SAXS_respath = '';
handles.SAXS_fit = [];
handles.sm = 0.30;
handles.edit_max_s.String = sprintf('%5.3f',handles.sm);
handles.no_selection = false;
handles.copy = false;

if ~isfield(model,'selected') || isempty(model.selected)
    model.selected{1} = model.current_structure;
    handles.no_selection = true;
end;

hMain.auxiliary=[hMain.auxiliary hObject];

% Choose default command line output for SAXS_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SAXS_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SAXS_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy = true;
update_plot(handles);


% --- Executes on button press in pushbutton_load_data.
function pushbutton_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load SAXS curve');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('SAXS curve loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name,ext] = fileparts(fname);
    handles.SAXS_datafile = strcat(name,ext);
    handles.SAXS_path = pname;
    handles.SAXS_curve = load_SAXS_curve(fullfile(pname,fname));
    handles.sm = max(handles.SAXS_curve(:,1));
    handles.edit_max_s.String = sprintf('%5.3f',handles.sm);
    handles.SAXS_fit = [];
end;
cd(my_path);
update_plot(handles);



% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.SAXS_curve)
    add_msg_board('Warning. No SAXS curve loaded yet');
    return
end;

pdbfile = sprintf('SAXS_test');
to_be_deleted = sprintf('SAXS_test*.*');
my_dir = pwd;
cd(handles.SAXS_path);
my_fig = gcf;
set(my_fig,'Pointer','watch');
drawnow;
wr_pdb_selected(pdbfile,'SAXS');
[chi2,~,~,result,fit] = fit_SAXS_by_crysol(handles.SAXS_datafile,pdbfile,handles.sm);
if isempty(chi2) || isnan(chi2)
    handles.SAXS_chi2 = 1e6;
    handles.SAXS_fit = [];
    add_msg_board('ERROR: SAXS fitting failed');
    add_msg_board(result);
else
    handles.SAXS_fit = fit;
    handles.SAXS_chi2 = chi2;
end
delete(to_be_deleted);
set(my_fig,'Pointer','arrow');
cd(my_dir);
drawnow;

update_plot(handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.no_selection
    model = rmfield(model,'selected');
end
delete(handles.figure1);

% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function update_plot(handles)

if isempty(handles.SAXS_curve)
    add_msg_board('Warning. No SAXS curve loaded yet');
    return
end;

if handles.copy
    figure;
    handles.copy = false;
else
    axes(handles.axes_plot);
    cla;
end;

if isempty(handles.SAXS_fit)
    plot(handles.SAXS_curve(:,1),handles.SAXS_curve(:,2),'k');
    title(sprintf('SAXS curve %s',handles.SAXS_datafile));
else
    plot(handles.SAXS_fit(:,1),handles.SAXS_fit(:,2),'k');
    hold on
    plot(handles.SAXS_fit(:,1),handles.SAXS_fit(:,3),'Color',[0.75,0,0]);
    title(sprintf('SAXS fit with chi^2 = %4.2f',handles.SAXS_chi2));
end;
xlabel('q');
ylabel('I');

guidata(handles.pushbutton_copy,handles);

function curve = load_SAXS_curve(fname)

fid = fopen(fname);
if fid==-1
    curve = [];
    add_msg_board('Warning. Loading of SAXS curve failed');
    return;
end;
nl=0;
curve = zeros(10000,4);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %         fprintf(1,'%s\n',tline); % echo for debugging
    if nl > 0 % skip first line
        dataset = str2num(tline);
        ncol = length(dataset);
        curve(nl,1:ncol) = dataset;
    end;
    nl = nl + 1;
end;
curve = curve(1:nl-1,:);
curve(:,1) = 10*curve(:,1);
fclose(fid);



function edit_max_s_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_s as text
%        str2double(get(hObject,'String')) returns contents of edit_max_s as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,1,0.3,'%5.3f',0);
if v ~= handles.sm
    handles.SAXS_fit = [];
end
handles.sm = v;
update_plot(handles);

% --- Executes during object creation, after setting all properties.
function edit_max_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
