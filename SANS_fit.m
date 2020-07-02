function varargout = SANS_fit(varargin)
% SANS_FIT MATLAB code for SANS_fit.fig
%      SANS_FIT, by itself, creates a new SANS_FIT or raises the existing
%      singleton*.
%
%      H = SANS_FIT returns the handle to a new SANS_FIT or the handle to
%      the existing singleton*.
%
%      SANS_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SANS_FIT.M with the given input arguments.
%
%      SANS_FIT('Property','Value',...) creates a new SANS_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SANS_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SANS_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SANS_fit

% Last Modified by GUIDE v2.5 05-Apr-2019 12:23:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SANS_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @SANS_fit_OutputFcn, ...
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


% --- Executes just before SANS_fit is made visible.
function SANS_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SANS_fit (see VARARGIN)

% global MMM_icon
global hMain
global model

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

handles.SANS_curve = [];
handles.SANS_datafile = '';
handles.SANS_path = '';
handles.SANS_illres = '';
handles.SANS_respath = '';
handles.SANS_fit = [];
handles.D2O = 0;
handles.D2O_fit = 0;
handles.D2O_variation = 0;
handles.no_selection = false;
handles.copy = false;

if ~isfield(model,'selected') || isempty(model.selected)
    model.selected{1} = model.current_structure;
    handles.no_selection = true;
end;

hMain.auxiliary=[hMain.auxiliary hObject];

% Choose default command line output for SANS_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SANS_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SANS_fit_OutputFcn(hObject, eventdata, handles) 
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

[fname,pname]=uigetfile('*.dat','Load SANS curve');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('SANS curve loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name,ext] = fileparts(fname);
    handles.SANS_datafile = strcat(name,ext);
    handles.SANS_path = pname;
    handles.SANS_curve = load_SANS_curve(fullfile(pname,fname));
    handles.SANS_fit = [];
end;
cd(my_path);
update_plot(handles);



% --- Executes on button press in pushbutton_load_illres.
function pushbutton_load_illres_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_illres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general


my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.res','Load SANS ill.res data');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Resolution parameter loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name,ext] = fileparts(fname);
    handles.SANS_illres = strcat(name,ext);
    if ~strcmpi(pname,handles.SANS_path)
        add_msg_board('ERROR. SANS curve and resolution file must be in the same directory');
        handles.SANS_illres = '';
    end;
end;
cd(my_path);
update_plot(handles);

% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.err = true;
options.lm = 50;
options.fb = 18;
options.D2O = handles.D2O;

if isempty(handles.SANS_curve)
    add_msg_board('Warning. No SANS curve loaded yet');
    return
end

pdbfile = sprintf('SANS_test');
to_be_deleted = sprintf('SANS_test*.*');
my_dir = pwd;
cd(handles.SANS_path);
my_fig = gcf;
set(my_fig,'Pointer','watch');
drawnow;
wr_pdb_selected(pdbfile,'SANS');
if handles.D2O_variation > 0
    handles.SANS_chi2 = 1e6;
    handles.SANS_fit = [];
    min_D2O = handles.D2O - handles.D2O_variation;
    if min_D2O < 0
        min_D2O = 0;
    end
    max_D2O = handles.D2O + handles.D2O_variation;
    if max_D2O > 1
        max_D2O = 1;
    end
    n_D2O = 1 + round((max_D2O-min_D2O)/0.01);
    D2O_vec = linspace(min_D2O,max_D2O,n_D2O);
    chi2_vec = nan(1,n_D2O);
    for kD = 1:n_D2O
        options.D2O = D2O_vec(kD);
        wr_pdb_selected(pdbfile,'SANS');
        [chi2,~,~,result,fit] = fit_SANS_by_cryson(handles.SANS_datafile,pdbfile,handles.SANS_illres,options);
        if isempty(chi2) || isnan(chi2)
            add_msg_board('ERROR: SANS fitting failed');
            add_msg_board(result);
        else
            chi2_vec(kD) = chi2;
            if chi2 < handles.SANS_chi2
                handles.SANS_chi2 = chi2;
                handles.SANS_fit = fit;
                handles.D2O_fit = options.D2O;
            end
        end
        delete(to_be_deleted);
    end
    figure(1000); clf;
    plot(D2O_vec,chi2_vec,'k.');
    xlabel('D_2O content [%]');
    ylabel('\xi ^2_{SANS)');
else
    handles.D2O_fit = options.D2O;
    [chi2,~,~,result,fit] = fit_SANS_by_cryson(handles.SANS_datafile,pdbfile,handles.SANS_illres,options);
    if isempty(chi2) || isnan(chi2)
        handles.SANS_chi2 = 1e6;
        handles.SANS_fit = [];
        add_msg_board('ERROR: SANS fitting failed');
        add_msg_board(result);
    else
        handles.SANS_fit = fit;
        handles.SANS_chi2 = chi2;
    end
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

if isempty(handles.SANS_curve)
    add_msg_board('Warning. No SANS curve loaded yet');
    return
end

if handles.copy
    figure;
    handles.copy = false;
else
    axes(handles.axes_plot);
    cla;
end

if isempty(handles.SANS_fit)
    plot(handles.SANS_curve(:,1),handles.SANS_curve(:,2),'k');
    title(sprintf('SANS curve %s',handles.SANS_datafile));
else
    plot(handles.SANS_fit(:,1),handles.SANS_fit(:,2),'k');
    hold on
    plot(handles.SANS_fit(:,1),handles.SANS_fit(:,3),'Color',[0.75,0,0]);
    title(sprintf('SANS fit with chi^2 = %4.2f (at %4.1f%% D2O)',handles.SANS_chi2,100*handles.D2O_fit));
end
xlabel('q');
ylabel('I');

guidata(handles.pushbutton_copy,handles);

function curve = load_SANS_curve(fname)

fid = fopen(fname);
if fid==-1
    curve = [];
    add_msg_board('Warning. Loading of SANS curve failed');
    return;
end
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
    end
    nl = nl + 1;
end
curve = curve(1:nl-1,:);
fclose(fid);



function edit_percent_D2O_Callback(hObject, eventdata, handles)
% hObject    handle to edit_percent_D2O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_percent_D2O as text
%        str2double(get(hObject,'String')) returns contents of edit_percent_D2O as a double

[v,handles] = edit_update_MMM(handles,hObject,0,100,0,'%4.1f');
handles.D2O = v/100;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_percent_D2O_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_percent_D2O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_variation_D2O_Callback(hObject, eventdata, handles)
% hObject    handle to edit_variation_D2O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_variation_D2O as text
%        str2double(get(hObject,'String')) returns contents of edit_variation_D2O as a double

[v,handles] = edit_update_MMM(handles,hObject,0,50,0,'%4.1f');
handles.D2O_variation = v/100;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_variation_D2O_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_variation_D2O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
