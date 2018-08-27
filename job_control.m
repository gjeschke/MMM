function varargout = job_control(varargin)
% JOB_CONTROL M-file for job_control.fig
%      JOB_CONTROL, by itself, creates a new JOB_CONTROL or raises the existing
%      singleton*.
%
%      H = JOB_CONTROL returns the handle to a new JOB_CONTROL or the handle to
%      the existing singleton*.
%
%      JOB_CONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JOB_CONTROL.M with the given input arguments.
%
%      JOB_CONTROL('Property','Value',...) creates a new JOB_CONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before job_control_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to job_control_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help job_control

% Last Modified by GUIDE v2.5 18-Mar-2011 16:34:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @job_control_OpeningFcn, ...
                   'gui_OutputFcn',  @job_control_OutputFcn, ...
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


% --- Executes just before job_control is made visible.
function job_control_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to job_control (see VARARGIN)

% Choose default command line output for job_control
handles.output = hObject;

handles=set_listbox(handles);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes job_control wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = job_control_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_jobs.
function listbox_jobs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_jobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_jobs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_jobs

global general

lm5='';
lm4='';
lm3='';
lm2='';
lm1='';
lm0='';
if ~isempty(handles.jobs),
    sel=get(handles.listbox_jobs,'Value');
    job=handles.jobs{sel};
    set(handles.text_type,'String',job.type);
    set(handles.text_logfile,'String',job.log);
    runtime=etime(clock,job.started);
    hours=floor(runtime/3600);
    minutes=floor((runtime-3600*hours)/60);
    seconds=round(runtime-3600*hours-60*minutes);
    set(handles.text_runtime,'String',sprintf('%i h %i min % i s.',hours,minutes,seconds));
    fid=fopen(fullfile(general.tmp_files, job.log),'r');
    nl=0;
    while 1
        tline = fgetl(fid);
        nl=nl+1;
        if ~ischar(tline), break, end
        lm5=lm4;
        lm4=lm3;
        lm3=lm2;
        lm2=lm1;
        lm1=lm0;
        lm0=tline;
    end;
    fclose(fid);
    loglines{1}=sprintf('line %i: %s',nl-5,lm5); 
    loglines{2}=sprintf('line %i: %s',nl-4,lm4); 
    loglines{3}=sprintf('line %i: %s',nl-3,lm3); 
    loglines{4}=sprintf('line %i: %s',nl-2,lm2); 
    loglines{5}=sprintf('line %i: %s',nl-1,lm1); 
    loglines{6}=sprintf('line %i: %s',nl,lm0); 
    set(handles.text_logtail,'String',loglines);
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox_jobs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_jobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_logfile.
function pushbutton_logfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_logfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global hMain


sel=get(handles.listbox_jobs,'Value');
if ~isempty(handles.jobs),
    logfile=fullfile(general.tmp_files, handles.jobs{sel}.log);
    hMain.report_file=logfile;
    report_editor;
end;

% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general


sel=get(handles.listbox_jobs,'Value');
if ~isempty(handles.jobs),
    out = timerfindall('Name', handles.jobs{sel}.name);
    stop(out);
    delete(out);
    delete(fullfile(general.tmp_files,[handles.jobs{sel}.name '.mat']));
end;
handles=set_listbox(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

function handles=set_listbox(handles)

global general

poi=0;
tmp_dir=dir(general.tmp_files); % load temporary directory
files=numel(tmp_dir);
for k=1:files, % check for and mark old files
    if ~tmp_dir(k).isdir,
        if strcmpi(tmp_dir(k).name(1:4),'job-'), % this is a job file
            poi=poi+1;
            load(tmp_dir(k).name);
            joblist{poi}=job.name;
            jobs{poi}=job;
            jobfiles{poi}=tmp_dir(k).name;
        end;
    end;
end;

if poi>0,
    handles.jobs=jobs;
    handles.jobfiles=jobfiles;
    set(handles.listbox_jobs,'String',joblist);
    set(handles.listbox_jobs,'Value',1);
    sel=1;
    job=handles.jobs{sel};
    set(handles.text_type,'String',job.type);
    set(handles.text_logfile,'String',job.log);
    runtime=etime(clock,job.started);
    hours=floor(runtime/3600);
    minutes=floor((runtime-3600*hours)/60);
    seconds=round(runtime-3600*hours-60*minutes);
    set(handles.text_runtime,'String',sprintf('%i h %i min % i s.',hours,minutes,seconds));
    fid=fopen(fullfile(general.tmp_files, job.log),'r');
    nl=0;
    lm5='';
    lm4='';
    lm3='';
    lm2='';
    lm1='';
    lm0='';
    while 1
        tline = fgetl(fid);
        nl=nl+1;
        if ~ischar(tline), break, end
        lm5=lm4;
        lm4=lm3;
        lm3=lm2;
        lm2=lm1;
        lm1=lm0;
        lm0=tline;
    end;
    fclose(fid);
    loglines{1}=sprintf('line %i: %s',nl-5,lm5); 
    loglines{2}=sprintf('line %i: %s',nl-4,lm4); 
    loglines{3}=sprintf('line %i: %s',nl-3,lm3); 
    loglines{4}=sprintf('line %i: %s',nl-2,lm2); 
    loglines{5}=sprintf('line %i: %s',nl-1,lm1); 
    loglines{6}=sprintf('line %i: %s',nl,lm0); 
    set(handles.text_logtail,'String',loglines);    
else
    set(handles.listbox_jobs,'String','<none>');
    set(handles.text_type,'String','<none>');
    set(handles.text_logfile,'String','<none>');
    set(handles.text_runtime,'String','<none>');
    set(handles.text_logtail,'String','<empty>');
    set(handles.listbox_jobs,'Value',1);
    handles.jobs=[];
    handles.jobfiles=[];
end;
