function varargout = ensemble_fit(varargin)
% ENSEMBLE_FIT MATLAB code for ensemble_fit.fig
%      ENSEMBLE_FIT, by itself, creates a new ENSEMBLE_FIT or raises the existing
%      singleton*.
%
%      H = ENSEMBLE_FIT returns the handle to a new ENSEMBLE_FIT or the handle to
%      the existing singleton*.
%
%      ENSEMBLE_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENSEMBLE_FIT.M with the given input arguments.
%
%      ENSEMBLE_FIT('Property','Value',...) creates a new ENSEMBLE_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ensemble_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ensemble_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ensemble_fit

% Last Modified by GUIDE v2.5 05-Jul-2019 12:50:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ensemble_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @ensemble_fit_OutputFcn, ...
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


% --- Executes just before ensemble_fit is made visible.
function ensemble_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ensemble_fit (see VARARGIN)

global MMM_icon
global hMain

% Choose default command line output for ensemble_fit
handles.output = hObject;

handles.block_size = 100;
handles.edit_block_size.String = sprintf('%i',handles.block_size); 
handles.threshold = 1;
handles.edit_threshold.String = sprintf('%3.1f',handles.threshold); 
handles.min_mod_depth = 0.005;
handles.edit_min_mod_depth.String = sprintf('%5.3f',handles.min_mod_depth);
handles.max_mod_depth = 0.65;
handles.edit_max_mod_depth.String = sprintf('%5.3f',handles.max_mod_depth);
handles.min_decay_rate = 0;
handles.edit_min_decay_rate.String = sprintf('%5.3f',handles.min_decay_rate);
handles.max_decay_rate = NaN;
handles.edit_max_decay_rate.String = sprintf('%5.3f',handles.max_decay_rate);
handles.cutoff = 0.9;
handles.edit_cutoff.String = sprintf('%5.3f',handles.cutoff);
handles.bckg_dim = 3;
handles.edit_bckg_dim.String = sprintf('%3.1f',handles.bckg_dim);
handles.SAS_max_harmonics = 50;
handles.edit_SAS_max_harmonics.String = sprintf('%i',handles.SAS_max_harmonics);
handles.SAS_Fibonacci = 18;
handles.edit_SAS_Fibonacci.String = sprintf('%i',handles.SAS_Fibonacci);
handles.crosslink_distance = 20;
handles.edit_crosslink_distance.String = sprintf('%4.1f',handles.crosslink_distance);
handles.ensemble_number = 0;
handles.text_ensemble_number.String = sprintf('%i',handles.ensemble_number);
handles.ensemble_size = NaN;
handles.text_ensemble_size.String = sprintf('%i',handles.ensemble_size);
handles.DEER_deficiency = NaN;
handles.text_DEER_deficiency.String = sprintf('%5.3f',handles.DEER_deficiency);
handles.DEER_deficiency_DEER_only = NaN;
handles.text_DEER_deficiency_DEER_only.String = sprintf('%5.3f',handles.DEER_deficiency_DEER_only);
handles.ensemble_size_DEER_only = NaN;
handles.text_ensemble_size_DEER_only.String = sprintf('%i',handles.ensemble_size_DEER_only);
handles.SAS_chi2 = NaN;
handles.text_SAS_chi2.String = sprintf('%5.3f',handles.SAS_chi2);
handles.SAS_chi2_SAS_only = NaN;
handles.text_SAS_chi2_SAS_only.String = sprintf('%5.3f',handles.SAS_chi2_SAS_only);
handles.ensemble_size_SAS_only = NaN;
handles.text_ensemble_size_SAS_only.String = sprintf('%i',handles.ensemble_size_SAS_only);
handles.integrative_loss = NaN;
handles.text_integrative_loss.String = sprintf('%5.3f',handles.integrative_loss);
handles.pdb_ID = 'MMM0';
handles.edit_PDB_ID.String = handles.pdb_ID;
axis(handles.axes_populations);
title('Conformer populations');
xlabel('Conformer number');
ylabel('Population');

handles.restraints = [];
handles.model_list = {};

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ensemble_fit wait for user response (see UIRESUME)
% uiwait(handles.fig_ensemble_fit);


% --- Outputs from this function are returned to the command line.
function varargout = ensemble_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'domain_ensemble.html');
webcall(entry,'-helpbrowser');


function edit_block_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_block_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_block_size as text
%        str2double(get(hObject,'String')) returns contents of edit_block_size as a double

[v,handles]=edit_update_MMM(handles,hObject,10,200,100,'%i',1);
handles.block_size = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_block_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_block_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_mod_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_min_mod_depth as a double

[v,handles]=edit_update_MMM(handles,hObject,0.0005,0.95,0.005,'%5.3f',0);
handles.block_size = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_min_mod_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_mod_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_max_mod_depth as a double

[v,handles]=edit_update_MMM(handles,hObject,0.1,1,0.65,'%5.3f',0);
handles.block_size = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_max_mod_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_decay_rate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_decay_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_decay_rate as text
%        str2double(get(hObject,'String')) returns contents of edit_min_decay_rate as a double

[v,handles]=edit_update_MMM(handles,hObject,0,1,0,'%5.3f',0);
handles.block_size = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_min_decay_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_decay_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_decay_rate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_decay_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_decay_rate as text
%        str2double(get(hObject,'String')) returns contents of edit_max_decay_rate as a double

[v,handles]=edit_update_MMM(handles,hObject,0.0001,5,NaN,'%5.3f',0);
handles.block_size = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_max_decay_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_decay_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_cutoff as a double

[v,handles]=edit_update_MMM(handles,hObject,0.5,1,0.9,'%5.3f',0);
handles.block_size = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAS_max_harmonics_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAS_max_harmonics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAS_max_harmonics as text
%        str2double(get(hObject,'String')) returns contents of edit_SAS_max_harmonics as a double

[v,handles]=edit_update_MMM(handles,hObject,15,50,50,'%i',1);
handles.block_size = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_SAS_max_harmonics_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAS_max_harmonics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAS_Fibonacci_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAS_Fibonacci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAS_Fibonacci as text
%        str2double(get(hObject,'String')) returns contents of edit_SAS_Fibonacci as a double

[v,handles]=edit_update_MMM(handles,hObject,10,18,18,'%i',1);
handles.block_size = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_SAS_Fibonacci_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAS_Fibonacci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_crosslink_distance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crosslink_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crosslink_distance as text
%        str2double(get(hObject,'String')) returns contents of edit_crosslink_distance as a double

[v,handles]=edit_update_MMM(handles,hObject,5,30,20,'%4.1f',0);
handles.block_size = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_crosslink_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crosslink_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_load_restraints.
function pushbutton_load_restraints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_restraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton_load_restraints

global general

my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load restraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Restraint loading cancelled by user');
    cd(my_path);
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name] = fileparts(fname);
    handles.save_path = pname;
    handles.save_name = name;
    
    hfig=gcf;
    set(hfig,'Pointer','watch');

    rtype = get_restraint_type(fullfile(pname,fname));
    unprocessed = true;
    
    switch rtype
        case 'GENERIC'
            restraints = rd_restraints(fullfile(pname,fname));
        case 'FLEX'
            restraints = rd_restraints(fullfile(pname,fname));
        case 'RIGIFLEX'
            restraints = rd_restraints_rigiflex(fullfile(pname,fname),unprocessed);
        otherwise
            add_msg_board(sprintf('ERROR: Restraint type %s is not supported in ensemble fits',rtype));
            restraints = [];
    end

    handles.restraints = restraints;
    
    if ~isempty(restraints)
        handles.pushbutton_load_restraints.Enable = 'off';
        if ~isempty(handles.model_list)
            handles.pushbutton_run_iteration.Enable = 'on';
        end
    end
    cd(my_path);
    set(hfig,'Pointer','arrow');
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_load_conformers.
function pushbutton_load_conformers_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_conformers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load model list from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Model list loading cancelled by user');
    cd(my_path);
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;

    hfig=gcf;
    set(hfig,'Pointer','watch');
    
    model_list = get_file_list(fullfile(pname,fname));
    
    [~,name] = fileparts(fname);
    handles.basname = name;
    handles.baspath = pname;
    tested_file = sprintf('%s_tested.dat',name);
    full_tested_file = fullfile(pname,tested_file);
    handles.tested_file = full_tested_file;
    
    ensemble_file = sprintf('%s_ensemble.dat',name);
    full_ensemble_file = fullfile(pname,ensemble_file);
    handles.ensemble_file = full_ensemble_file;

    tested_list = {};
    if exist(full_tested_file,'file')
        answer = questdlg(sprintf('Open file %s?                                           .',...
            tested_file),'Read file of already tested conformers');
        switch answer
            case 'Yes'
                tested_list = get_file_list(full_tested_file);
            case 'No'
                tested_list = {};
            case 'Cancel'
                add_msg_board('Model list loading cancelled by user');
                cd(my_path);
                set(hfig,'Pointer','arrow');
                return
            otherwise
                add_msg_board('Model list loading cancelled by user');
                cd(my_path);
                set(hfig,'Pointer','arrow');
                return
        end
    end
    
    ensemble_list = {};
    if exist(full_ensemble_file,'file')
        answer = questdlg(sprintf('Open file %s?                                                .',...
            ensemble_file),'Read existing ensemble file');
        switch answer
            case 'Yes'
                [ensemble_list,comments] = get_file_list(full_ensemble_file);
                if ~isempty(comments)
                    prev_iterations = str2double(comments{1});
                    if ~isnan(prev_iterations)
                        handles.ensemble_number = round(prev_iterations);
                        handles.text_ensemble_number.String = sprintf('%i',handles.ensemble_number);
                    end
                end
            case 'No'
                ensemble_list = {};
            case 'Cancel'
                add_msg_board('Model list loading cancelled by user');
                cd(my_path);
                set(hfig,'Pointer','arrow');
                return
            otherwise
                add_msg_board('Model list loading cancelled by user');
                cd(my_path);
                set(hfig,'Pointer','arrow');
                return
        end
    end
    
    handles.model_list = model_list;
    handles.tested_list = tested_list;
    handles.ensemble_list = ensemble_list;

    if ~isempty(handles.model_list)
        handles.pushbutton_load_conformers.Enable = 'off';
        if ~isempty(handles.restraints)
            handles.pushbutton_run_iteration.Enable = 'on';
        end
    end

    if ~isempty(handles.ensemble_list)
        handles.pushbutton_check_ensemble.Enable = 'on';
    end
    
set(hfig,'Pointer','arrow');

end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_run_iteration.
function pushbutton_run_iteration_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run_iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

max_num_conformers = length(handles.model_list) + length(handles.ensemble_list);
conformer_batch = cell(1,max_num_conformers);
newly_tested = cell(1,length(handles.model_list));

hfig=gcf;
set(hfig,'Pointer','watch');
    
poi = 0;
poi2 = 0;
% add conformers from the previous best ensemble
for k = 1:length(handles.ensemble_list)
    poi = poi + 1;
    conformer_batch{poi} = handles.ensemble_list{k};
end
if poi >= handles.block_size
    add_msg_board(sprintf('Warning: Size %i of existing ensemble at block size limit. Only one new model will be tested.',poi));
end
% add conformers from the model list, unless they were already tested
for k = 1:length(handles.model_list)
    to_be_tested = true;
    for k2 = 1:length(handles.tested_list)
        if strcmp(handles.model_list{k},handles.tested_list{k2})
            to_be_tested = false;
            break
        end
    end
    if to_be_tested
        for k2 = 1:length(handles.ensemble_list)
            if strcmp(handles.model_list{k},handles.ensemble_list{k2})
                to_be_tested = false;
                break
            end
        end
    end
    if to_be_tested
        poi = poi + 1;
        conformer_batch{poi} = handles.model_list{k};
        poi2 = poi2 + 1;
        newly_tested{poi2} = handles.model_list{k};
    end
    if poi >= handles.block_size
        break
    end
end
conformer_batch = conformer_batch(1:poi);
newly_tested = newly_tested(1:poi2);

if isempty(conformer_batch) || isempty(newly_tested)
    add_msg_board('Warning: No new conformer found in list.');
    if isempty(handles.ensemble_list)
        add_msg_board('ERROR: No ensemble. Check input data.');
    else
        add_msg_board('Ensemble appears to be fully converged.');
    end
    set(hfig,'Pointer','arrow');
    return
else
    add_msg_board(sprintf('Populations will be fitted for a block of %i conformers',poi));
end

% Generate option structure for fitting/restraint evaluation
options.threshold = handles.threshold/100; % convert percentage to fraction
options.old_ensemble_size = length(handles.ensemble_list);
options.plot_axes = handles.axes_populations; % for display updates during fitting
options.text_DEER_fom0 = handles.text_DEER_deficiency_DEER_only; % for display updates during fitting
options.text_SAS_fom0 = handles.text_SAS_chi2_SAS_only; % for display updates during fitting
options.text_DEER_fom = handles.text_DEER_deficiency; % for display updates during fitting
options.text_SAS_fom = handles.text_SAS_chi2; % for display updates during fitting
options.text_loss = handles.text_integrative_loss;
options.text_DEER_size = handles.text_ensemble_size_DEER_only;
options.text_SAS_size = handles.text_ensemble_size_SAS_only;
options.text_iteration_counter = handles.text_fit_iteration;
options.text_ensemble_size = handles.text_ensemble_size;
options.text_status = handles.text_current_status;
options.text_info_status = handles.text_info_status;
options.depth = [handles.min_mod_depth,handles.max_mod_depth];
options.dim = handles.bckg_dim;
options.kdec = [handles.min_decay_rate, handles.max_decay_rate];
options.cutoff = handles.cutoff;
options.lm = handles.SAS_max_harmonics;
options.fb = handles.SAS_Fibonacci;
options.xlmax = handles.crosslink_distance;

log_file_it = sprintf('%s_ensemble_iteration_%i.log',handles.basname,handles.ensemble_number+1);
options.logfile = fullfile(handles.baspath,log_file_it);
log_fid = fopen(options.logfile,'wt');
if log_fid == -1
    add_msg_board(sprintf('Warning: Log file %s for fitting cannot be written.',options.logfile));
else
    fprintf(log_fid,'*** Log file for MMM integrative ensemble fit iteration %i ***\n\n',handles.ensemble_number+1);
    fprintf(log_fid,'Populations of %i conformers are fitted, of which %i belong to the ensemble from the previous iteration\n\n',...
        length(conformer_batch),options.old_ensemble_size);
    fprintf(log_fid,'Population threshold is %3.1f%% of the maximum population\n\n',100*options.threshold);
    fclose(log_fid);
end

[included,populations,SAS_bsl] = fit_ensemble(handles.restraints,conformer_batch,options);

if ~isempty(included) && ~isempty(populations)
    ensemble_list = cell(1,length(included));
    for kc = 1:length(included)
        ensemble_list{kc} = conformer_batch{included(kc)};
    end
    handles.ensemble_list = ensemble_list;
    handles.populations = populations;
    handles.SAS_bsl = SAS_bsl;
    handles.ensemble_number = handles.ensemble_number + 1;
    handles.text_ensemble_number.String = sprintf('%i',handles.ensemble_number);
    handles.text_ensemble_size.String = sprintf('%i',length(included));
    fname = put_file_list(handles.tested_file,newly_tested,1:length(newly_tested),true);
    if isempty(fname)
        add_msg_board('Warning: Newly tested conformers could not be added to list of already tested conformers.');
    else
        add_msg_board(sprintf('Newly tested conformers added to list %s',fname));
        tested_list = get_file_list(handles.tested_file);
        handles.tested_list = tested_list;
    end
    ensemble_file_it = sprintf('%s_ensemble_iteration_%i.dat',handles.basname,handles.ensemble_number);
    full_ensemble_file_it = fullfile(handles.baspath,ensemble_file_it);
    comment{1} = sprintf('%i',handles.ensemble_number);
    fname = put_file_list(full_ensemble_file_it,conformer_batch,included,false,populations,comment);
    if isempty(fname)
        add_msg_board('Warning: Ensemble file for this iteration could not be written.');
    else
        add_msg_board(sprintf('Ensemble file for this iterations written to %s',fname));
    end
    comment{1} = sprintf('%i',handles.ensemble_number);
    fname = put_file_list(handles.ensemble_file,conformer_batch,included,false,populations,comment);
    if isempty(fname)
        add_msg_board('Warning: File for best ensemble could not be updated.');
    else
        add_msg_board(sprintf('Ensemble file %s was updated',fname));
        hMain.report_file = handles.ensemble_file;
        report_editor;
    end
else
    add_msg_board('ERROR: Ensemble fit failed; ensemble was not updated');
end

if ~isempty(handles.ensemble_list)
    handles.pushbutton_check_ensemble.Enable = 'on';
end

set(hfig,'Pointer','arrow');

guidata(hObject,handles);

% --- Executes on button press in pushbutton_check_ensemble.
function pushbutton_check_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_check_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Generate option structure for fitting/restraint evaluation
options.threshold = handles.threshold/100; % convert percentage to fraction
options.ensemble_size = length(handles.ensemble_list);
options.depth = [handles.min_mod_depth,handles.max_mod_depth];
options.dim = handles.bckg_dim;
options.kdec = [handles.min_decay_rate, handles.max_decay_rate];
options.cutoff = handles.cutoff;
options.lm = handles.SAS_max_harmonics;
options.fb = handles.SAS_Fibonacci;
options.xlmax = handles.crosslink_distance;

options.pdb_file = sprintf('%s_ensemble.pdb',fullfile(handles.baspath,handles.basname));
options.script_file = sprintf('%s_ensemble.mmm',fullfile(handles.baspath,handles.basname));
options.pdb_ID = handles.pdb_ID;

hfig=gcf;
set(hfig,'Pointer','watch');

mk_ensemble(handles.ensemble_list,handles.populations,handles.restraints,options);
[~,snum] = add_pdb(options.pdb_file);
stag = mk_address(snum);
stag = stag(2:end-1);

[sigma,pair_rmsd] = ensemble_comparison(stag,[],handles.populations,[],[],[],options);
ermsd_str = sprintf('Ensemble rmsd is %5.2f Å',sigma);
add_msg_board(ermsd_str);
axes(handles.axes_populations);
cla
[n1,n2] = size(pair_rmsd{1});
plot(1,1,'k.');
plot(n1,n2,'k.');
axis tight
xlabel('Conformer number');
ylabel('Conformer number');
image(pair_rmsd{1},'CDataMapping','scaled');
title(ermsd_str);
colorbar;
axis equal

set(hfig,'Pointer','arrow');

guidata(hObject,handles);

function edit_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_threshold as a double


% --- Executes during object creation, after setting all properties.
function edit_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_ensemble_fit_CloseRequestFcn(handles.fig_ensemble_fit, eventdata, handles);

% --- Executes when user attempts to close fig_ensemble_fit.
function fig_ensemble_fit_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fig_ensemble_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function edit_bckg_dim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg_dim as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg_dim as a double


% --- Executes during object creation, after setting all properties.
function edit_bckg_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PDB_ID_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PDB_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PDB_ID as text
%        str2double(get(hObject,'String')) returns contents of edit_PDB_ID as a double

v = get(hObject,'String');
if length(v) >= 4
    handles.pdb_ID = v(1:4);
else
    handles.pdb_ID = 'XXXX';
    handles.pdb_ID(1:length(v)) = v;
end
hObject.String = handles.pdb_ID;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_PDB_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PDB_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
