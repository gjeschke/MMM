function varargout = domain_ensemble(varargin)
% DOMAIN_ENSEMBLE MATLAB code for domain_ensemble.fig
%      DOMAIN_ENSEMBLE, by itself, creates a new DOMAIN_ENSEMBLE or raises the existing
%      singleton*.
%
%      H = DOMAIN_ENSEMBLE returns the handle to a new DOMAIN_ENSEMBLE or the handle to
%      the existing singleton*.
%
%      DOMAIN_ENSEMBLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOMAIN_ENSEMBLE.M with the given input arguments.
%
%      DOMAIN_ENSEMBLE('Property','Value',...) creates a new DOMAIN_ENSEMBLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before domain_ensemble_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to domain_ensemble_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help domain_ensemble

% Last Modified by GUIDE v2.5 20-Dec-2015 09:58:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @domain_ensemble_OpeningFcn, ...
                   'gui_OutputFcn',  @domain_ensemble_OutputFcn, ...
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


% --- Executes just before domain_ensemble is made visible.
function domain_ensemble_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to domain_ensemble (see VARARGIN)

% global MMM_icon
global hMain
global model

% Choose default command line output for domain_ensemble
handles.output = hObject;

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

if ~exist('model','var') || ~isfield(model,'current_structure'),
    set(hMain.popupmenu_view,'Value',1);
    % initialize display
    axes(hMain.axes_model);
    cla reset;
    axis equal
    axis off
    set(gca,'Clipping','off');
    set(gcf,'Renderer','opengl');
    hold on
    hMain.camlight=camlight;
    guidata(hMain.axes_model,hMain);
    hMain.virgin=0;
    adr = '[none]';
%     add_msg_board('Error: Cannot build a domain model, since no current structure exists.');
%     guidata(hObject,handles);
%     delete(hObject);
%     return;
else
    adr=mk_address(model.current_structure);
end

set(handles.fig_domain_model,'Name',sprintf('Monte Carlo domain ensemble at structure %s',adr));

handles.max_time = 2;
set(handles.edit_max_time,'String',sprintf('%5.2f',handles.max_time));

set(handles.text_time_left,'String','Idle.');
set(handles.text_time_left,'ForegroundColor','b');

handles.n_restraints = 0;
handles.n_monitor = 0;

handles.ensemble = 20;
handles.ensemble_user = handles.ensemble;
set(handles.edit_ensemble,'String',sprintf('%i',handles.ensemble_user));

handles.p_model = 0.5;
handles.p_model_user = handles.p_model;
set(handles.edit_p_model,'String',sprintf('%5.2f',handles.p_model_user));

handles.max_trials = 500000000;
set(handles.edit_max_trials,'String',sprintf('%i',handles.max_trials));

handles.min_approach = 1.2;
set(handles.edit_min_approach,'String',sprintf('%5.2f',handles.min_approach));

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes domain_ensemble wait for user response (see UIRESUME)
% uiwait(handles.fig_domain_model);


% --- Outputs from this function are returned to the command line.
function varargout = domain_ensemble_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'domain_ensemble.html');
webcall(entry,'-helpbrowser');


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global hMain
% global model
% 
% hMain.domain_ensemble_plot=false;
% model.current_structure=handles.current_structure;
% update_current_structure;
delete(handles.fig_domain_model);

% --- Executes on button press in pushbutton_restraints.
function pushbutton_restraints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general
global residue_defs

snum = [];
if isfield(model,'current_structure'),
    snum=model.current_structure;
end

my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load restraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Restraint loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name] = fileparts(fname);
    handles.save_path = pname;
    handles.save_name = name;
    
    hfig=gcf;
    set(hfig,'Pointer','watch');
    restraints=rd_restraints(fullfile(pname,fname));

    handles.ensemble=restraints.ensemble;

    new_template=false;
    if isfield(restraints,'PDB'),
        if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
            id=tag2id(restraints.PDB,model.structure_tags);
            if ~isempty(id),
                snum=id;
                model.current_structure=snum;
            else
                button = questdlg(sprintf('Restraint file specifies structure %s, while current structure is %s. Do you want to load specified template?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
                if strcmp(button,'Yes'),
                    new_template=true;
                    fname=get_pdb_file(restraints.PDB);
                    [message,snum]=add_pdb(fname);
                    handles.template_seq = model.structures{snum}(1).sequence;
                    adr=mk_address(snum);
                    set(handles.fig_domain_model,'Name',sprintf('Monte Carlo domain ensemble at structure %s',adr));
                    if message.error,
                        add_msg_board(sprintf('ERROR: Specified structure PDB file %s could not be retrieved from server',restraints.PDB));
                        add_msg_board(message.text);
                        set(hfig,'Pointer','arrow');
                        return
                    else
                        model.current_structure=snum;
                    end
                end
            end
        end
    end

    restraints.template = snum;
    if ~isempty(snum),
        handles.template_seq = model.structures{snum}(1).sequence;
    else
        handles.template_seq = '';
    end
    
    [restrain,monitor,restraints,number,bnumber,number_monitor,cancelled]=process_domain_restraints(handles,restraints);
    if isfield(model,'current_structure');
        snum = model.current_structure;
    else
        snum = [];
    end
    if cancelled,
        add_msg_board('Processing of restraints cancelled.');
        set(hfig,'Pointer','arrow');
        return
    end
    handles.restrain = restrain;
    handles.monitor = monitor;
    handles.restraints = restraints;
    handles.n_restraints = number;
    handles.b_restraints = bnumber;
    handles.n_monitor = number_monitor;
    add_msg_board(sprintf('%i distribution restraints have been read.',number));
    
    if isfield(restraints,'ensemble'),
        handles.ensemble = restraints.ensemble;
    else
        handles.ensemble = 0;
    end
    if isfield(restraints,'uncertainty'),
        handles.p_model = restraints.uncertainty/restraints.rescale;
    end
    Na_indices = [];
    if isfield(restraints,'Nanchor'),
        [indices,message]=resolve_address(restraints.Nanchor);
        Na_indices = indices;
        if message.error,
            add_msg_board(sprintf('ERROR: N-terminal anchor residue %s does not exist.',restraints.Nanchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        end
        Nname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name; 
        id = tag2id(Nname,upper(residue_defs.restags));
        N_slc = residue_defs.single_letter_code(id);
        indices_p = indices;
        indices_p(4) = indices_p(4) - 1; % index of previous residue before N terminal anchor
        Npname = model.structures{indices_p(1)}(indices_p(2)).residues{indices_p(3)}.info(indices_p(4)).name; 
        id = tag2id(Npname,upper(residue_defs.restags));
        Np_slc = residue_defs.single_letter_code(id);
        restraints.Nanchor_p = mk_address(indices_p);
        anchorNp = zeros(4,3);
        [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor_p),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor N atom does not exist.',restraints.Nanchor_p));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorNp(1,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor_p),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor CA atom does not exist.',restraints.Nanchor_p));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorNp(2,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor_p),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor C atom does not exist.',restraints.Nanchor_p));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorNp(3,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor_p),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor O atom does not exist.',restraints.Nanchor_p));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorNp(4,:) = coor;
        end
        anchorN = zeros(4,3);
        [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s N atom does not exist.',restraints.Nanchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorN(1,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s CA atom does not exist.',restraints.Nanchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorN(2,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s C atom does not exist.',restraints.Nanchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorN(3,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s O atom does not exist.',restraints.Nanchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorN(4,:) = coor;
        end
        restraints.anchorN = anchorN;
        restraints.anchorNp = anchorNp;
        restraints.Nseq = [Np_slc N_slc];
    else
        restraints.anchorN = [];
        restraints.anchorNp = [];
        restraints.Nseq = '';
    end
    Ca_indices = [];
    if isfield(restraints,'Canchor'),
        [indices,message]=resolve_address(restraints.Canchor);
        Ca_indices = indices;
        if message.error,
            add_msg_board(sprintf('ERROR: C-terminal anchor residue %s does not exist.',restraints.Canchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        end
        Cname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name; 
        id = tag2id(Cname,upper(residue_defs.restags));
        C_slc = residue_defs.single_letter_code(id);
        indices_n = indices;
        indices_n(4) = indices_n(4) + 1; % index of next residue after C-terminal anchor
        Cnname = model.structures{indices_n(1)}(indices_n(2)).residues{indices_n(3)}.info(indices_n(4)).name; 
        id = tag2id(Cnname,upper(residue_defs.restags));
        Cn_slc = residue_defs.single_letter_code(id);
        restraints.Canchor_n = mk_address(indices_n);
        anchorCn = zeros(4,3);
        [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor_n),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor N atom does not exist.',restraints.Canchor_n));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorCn(1,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor_n),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor CA atom does not exist.',restraints.Canchor_n));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorCn(2,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor_n),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor C atom does not exist.',restraints.Canchor_n));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorCn(3,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor_n),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor O atom does not exist.',restraints.Canchor_n));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorCn(4,:) = coor;
        end
        anchorC = zeros(4,3);
        [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s N atom does not exist.',restraints.Canchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorC(1,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s CA atom does not exist.',restraints.Canchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorC(2,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s C atom does not exist.',restraints.Canchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorC(3,:) = coor;
        end
        [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor),'coor');
        if msg.error,
            add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s O atom does not exist.',restraints.Canchor));
            add_msg_board(message.text);
            set(hfig,'Pointer','arrow');
            return
        else
            anchorC(4,:) = coor;
        end
        restraints.anchorC = anchorC;
        restraints.anchorCn = anchorCn;
        restraints.Cseq = [C_slc Cn_slc];
    else
        restraints.anchorC = [];
        restraints.anchorCn = [];
        restraints.Cseq = '';
    end
end
restraints.Na_indices = Na_indices;
restraints.Ca_indices = Ca_indices;
if handles.ensemble == 0,
    handles.ensemble = handles.ensemble_user;
else
    handles.ensemble_user = handles.ensemble;
    set(handles.edit_p_model,'String',sprintf('%4.2f',handles.p_model_user));
end

if handles.p_model == 0,
    handles.p_model = handles.p_model_user;
else
    handles.p_model_user = handles.p_model;
    set(handles.edit_ensemble,'String',sprintf('%i',handles.ensemble_user));
end

handles.restraints = restraints;
handles.restrain = restrain;

set(handles.pushbutton_run,'Enable','on');
set(handles.pushbutton_restraints,'Enable','off');
set(hfig,'Pointer','arrow');
cd(my_path);
guidata(hObject,handles);

function [restrain,monitor,restraints,number,bnumber,number_monitor,cancelled]=process_domain_restraints(handles,restraints)

global model
global hMain
global rotamer_libraries
global ligand_libraries

cancelled=false;

restrain = [];
number = 0;
bnumber = 0;
number_monitor = 0;

if ~isfield(restraints,'sequence'),
    add_msg_board('ERROR: Sequence specification is missing.'); 
    restrain=[];
    cancelled = true;
    return;
end

for k = 1:length(restraints.sequence),
    restrain(k).secondary = 0;
    restrain(k).label = [];
    restrain(k).cis = 0;
    restrain(k).r_beacon = [];
    restrain(k).r_intern = [];
    restrain(k).oligomer = [];
    restrain(k).depth = [];
end

monitor = restrain;

if ~isfield(restraints,'domain_start') || ~isfield(restraints,'domain_end'),
    add_msg_board('ERROR: Domain is not specified in restraint file.'); 
    restrain=[];
    cancelled = true;
    return;
else
    poi1 = strfind(restraints.domain_start,'}');
    if ~isempty(poi1);
        poi2 = strfind(restraints.domain_end,'}');
        if ~strcmpi(restraints.domain_start(1:poi1),restraints.domain_end(1:poi2)),
            add_msg_board('ERROR: Addresses of first and last residue of the domain are inconsistent.'); 
            add_msg_board(sprintf('%s does not match %s',restraints.domain_start(1:poi1),restraints.domain_end(1:poi2))); 
            restrain=[];
            cancelled = true;
            return;
        else
            chain_model = restraints.domain_start(1:poi1);
        end
    else
        poi1 = strfind(restraints.domain_start,')');
        if ~isempty(poi1);
            poi2 = strfind(restraints.domain_end,')');
            if ~strcmpi(restraints.domain_start(1:poi1),restraints.domain_end(1:poi2)),
                add_msg_board('ERROR: Addresses of first and last residue of the domain are inconsistent.'); 
                add_msg_board(sprintf('%s does not match %s',restraints.domain_start(1:poi1),restraints.domain_end(1:poi2))); 
                restrain=[];
                cancelled = true;
                return;
            else
                chain_model = restraints.domain_start(1:poi1);
            end
        else
            chain_model = '';
            poi1 = 0;
        end
    end
    res1 = str2double(restraints.domain_start(poi1+1:end));
    rese = str2double(restraints.domain_end(poi1+1:end));
    if isnan(res1) || isnan(rese),
        add_msg_board('ERROR: Residue numbers of domain could not be recognized.'); 
        restrain=[];
        cancelled = true;
        return;
    end
end

restraints.chain_model = chain_model;
restraints.res1 = res1;
restraints.rese = rese;

for k = res1:rese,
    adr = sprintf('%s%i',chain_model,k);
    [indices,message]=resolve_address(adr);
    if message.error ~=2 && message.error ~=13,
        add_msg_board(sprintf('ERROR: Residue %s does exist in the structure. Modelling impossible as this would create a clash.',adr)); 
        address = mk_address(indices,1);
        add_msg_board(address);
        restrain=[];
        cancelled = true;
        return;
    end
end

if isfield(restraints,'DEER'),
    poi = 0;
    llist = cell(0);
    label = cell(0);
    T = zeros(1,length(restraints.DEER));
    for k = 1:length(restraints.DEER),
        sep = strfind(restraints.DEER(k).label,'|');
        if isempty(sep)
            label1 = restraints.DEER(k).label;
            label2 = label1;
        else
            label1 = restraints.DEER(k).label(1:sep-1);
            label2 = restraints.DEER(k).label(sep+1:end);
        end
        [indices,message]=resolve_address(restraints.DEER(k).adr1);
        if message.error ~= 2 && message.error ~=13, % this site exists and is thus a beacon
            if message.error,
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr1)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            llist{poi} = restraints.DEER(k).adr1;
            T(poi) = restraints.DEER(k).T;
            label{poi} = label1;
            restraints.DEER(k).type1 = 1;
            restraints.DEER(k).indices1 = indices;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,label1) || strcmpi(rotamer_libraries(kr).tc,label1),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel1 = NO;
                    restraints.DEER(k).type1 = 0;
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,label1) || strcmpi(ligand_libraries(kr).tc,label1),
                    Tvec = ligand_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,ligand_libraries(kr).files);
                    spin = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                    restraints.DEER(k).NO_rel1 = spin;
                    restraints.DEER(k).type1 = 0;
                end
            end
        end
        [indices,message]=resolve_address(restraints.DEER(k).adr2);
        if message.error ~= 2 && message.error ~= 13, % this site exists and is thus a beacon
            if message.error,
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr2)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            llist{poi} = restraints.DEER(k).adr2;
            T(poi) = restraints.DEER(k).T;
            label{poi} = label2;
            restraints.DEER(k).type2 = 1;
            restraints.DEER(k).indices2 = indices;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,label2) || strcmpi(rotamer_libraries(kr).tc,label2)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel2 = NO;
                    restraints.DEER(k).type2 = 0;
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,label2) || strcmpi(ligand_libraries(kr).tc,label2),
                    Tvec = ligand_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,ligand_libraries(kr).files);
                    spin = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                    restraints.DEER(k).NO_rel2 = spin;
                    restraints.DEER(k).type2 = 0;
                end
            end
        end
    end
    labels = get_labels(llist,label,T,handles);
    for k = 1:length(restraints.DEER),
        clabel = restraints.DEER(k).label;
        cT = restraints.DEER(k).T;
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 0, % internal restraint
            resa = separate_address(restraints.DEER(k).adr1);
            resb = separate_address(restraints.DEER(k).adr2);
            if isempty(resa) || isempty(resb),
                add_msg_board(sprintf('ERROR: Residue address %s or %s is invalid.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO1 = restraints.DEER(k).NO_rel1;
            NO2 = restraints.DEER(k).NO_rel2;
            if resa < resb,
                exch = resa; resa = resb; resb = exch;
                exch = NO1; NO1 = NO2; NO2 = exch;
            end
            if restraints.DEER(k).r ~= 0,
                [restrain,number,bnumber] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT);
            else
                [monitor,number_monitor] = mk_internal_restraint(monitor,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 ==0, % beacon restraint, first residue beacon
            indices = restraints.DEER(k).indices1;
            resb = separate_address(restraints.DEER(k).adr1);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = separate_address(restraints.DEER(k).adr2);
            if isempty(res_loop),
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel2;
            if restraints.DEER(k).r ~= 0,
                [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT,indices,resb);
            else
                [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT,indices,resb);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 1, % beacon restraint, second residue beacon
            indices = restraints.DEER(k).indices2;
            resb = separate_address(restraints.DEER(k).adr2);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = separate_address(restraints.DEER(k).adr1);
            if isempty(res_loop),
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr1)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel1;
            if restraints.DEER(k).r ~= 0,
               [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT,indices,resb);
            else
               [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT,indices,resb);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 == 1, % nonsensical restraint inside defined structure
            add_msg_board(sprintf('Warning: Restraint between residues %s and %s inside defined structure will be ignored.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
        end
    end
end

if isfield(restraints,'depth'),
    for k = 1:length(restraints.depth),
        clabel = restraints.depth(k).label;
        cT = restraints.depth(k).T;
        if strcmpi(restraints.depth(k).label,'CA'),
            NO = [];
        else
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.depth(k).label),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.depth(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
        end
        res = separate_address(restraints.depth(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.depth(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.depth(k).r ~= 0,
            [restrain,number,bnumber] = mk_depth_restraint(restrain,NO,res,restraints.depth(k).r,restraints.depth(k).sigr,res1,number,bnumber,clabel,cT);
        else
            [monitor,number_monitor] = mk_depth_restraint(monitor,NO,res,restraints.depth(k).r,restraints.depth(k).sigr,res1,number_monitor,bnmon,clabel,cT);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'oligomer'),
    for k = 1:length(restraints.oligomer),
        clabel = restraints.oligomer(k).label;
        cT = restraints.oligomer(k).T;
        if strcmpi(restraints.oligomer(k).label,'CA'),
            NO = [];
        else
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.oligomer(k).label),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.oligomer(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,restraints.oligomer(k).label) || strcmpi(ligand_libraries(kr).tc,restraints.oligomer(k).label),
                    libname = id2tag(1,ligand_libraries(kr).files);
                    NO = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                end
            end
        end
        res = separate_address(restraints.oligomer(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.oligomer(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.oligomer(k).r ~= 0,
            [restrain,number,bnumber] = mk_oligomer_restraint(restrain,NO,res,restraints.oligomer(k).multi,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number,bnumber,clabel,cT);
        else
            [monitor,number_monitor] = mk_oligomer_restraint(monitor,NO,res,restraints.oligomer(k).multi,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number_monitor,bnmon,clabel,cT);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'cispeptides'),
    for k = 1:length(restraints.cispeptides),
        kr = restraints.cispeptides(k) - res1 + 1;
        restrain(kr).cis = 1;
    end
end

if isfield(restraints,'aprop'),
    for k = 1:length(restraints.aprop),
        res = separate_address(restraints.aprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.aprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).aprop = restraints.aprop(k).prop;
    end
end

if isfield(restraints,'bprop'),
    for k = 1:length(restraints.bprop),
        res = separate_address(restraints.bprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.bprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).bprop = restraints.bprop(k).prop;
    end
end

if isfield(restraints,'pprop'),
    for k = 1:length(restraints.pprop),
        res = separate_address(restraints.pprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.pprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).pprop = restraints.pprop(k).prop;
    end
end

if isfield(restraints,'helices'),
    for k = 1:length(restraints.helices),
        hpoi = strfind(restraints.helices(k).adr,'-');
        ha = str2double(restraints.helices(k).adr(1:hpoi-1));
        he = str2double(restraints.helices(k).adr(hpoi+1:end));
        if isnan(ha) || isnan(he),
            add_msg_board(sprintf('ERROR: Wrong helix specification %s.',restraints.helices(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence),
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0,
                restrain(ks).secondary = 1;
            end
            if off1 >=2 && off2 >=2,
                restrain(ks).secondary = 3;
            end
        end
    end
end

if isfield(restraints,'strands'),
    for k = 1:length(restraints.strands),
        hpoi = strfind(restraints.strands(k).adr,'-');
        ha = str2double(restraints.strands(k).adr(1:hpoi-1));
        he = str2double(restraints.strands(k).adr(hpoi+1:end));
        if isnan(ha) || isnan(he),
            add_msg_board(sprintf('ERROR: Wrong strand specification %s.',restraints.strands(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence),
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0,
                restrain(ks).secondary = 2;
            end
        end
    end
end

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general
global Ramachandran
global residue_defs
global hMain

parnum = 100; % number of trials performed in the parfor loop

disp_update = 200;

auxiliary = false;

restraints = handles.restraints;
restrain = handles.restrain;

secvec = zeros(1,length(restrain));
for k = 1:length(restrain)
    secvec(k) = restrain(k).secondary;
end

if ~isempty(restraints.Na_indices)
    Nanchor_res = restraints.Na_indices(4);
else
    Nanchor_res = [];
end

if ~isempty(restraints.Ca_indices)
    Canchor_res = restraints.Ca_indices(4);
else
    Canchor_res = [];
end

% Make auxiliary restraints by triangle bound smoothing
% find all beacon residues
beacons = zeros(1,1000);
bref = zeros(1000,2);
bpoi = 0;
for k = 1:length(restrain),
    for kr = 1: length(restrain(k).r_beacon),
        res = restrain(k).r_beacon(kr).resb;
        if isempty(find(beacons == res)),
            bpoi = bpoi+1;
            beacons(bpoi) = res;
            bref(bpoi,1) = k;
            bref(bpoi,2) = kr;
        end
    end
end
ni = length(restrain);
vert = ni + bpoi;
beacons = beacons(1:bpoi);
bref = bref(1:bpoi,:);
lower_bounds = zeros(vert);
upper_bounds = 1e6*ones(vert);
missing = ones(vert);
for k = 1:length(restrain)
    for kr = 1: length(restrain(k).r_beacon)
        if strcmp(restrain(k).r_beacon(kr).type,'Gaussian')
            b = find(beacons == restrain(k).r_beacon(kr).resb) + ni;
            rmean = restrain(k).r_beacon(kr).par1;
            sigr = restrain(k).r_beacon(kr).par2;
            lower_bounds(k,b) = rmean - 2*sigr;
            lower_bounds(b,k) = rmean - 2*sigr;
            upper_bounds(k,b) = rmean + 2*sigr;
            upper_bounds(b,k) = rmean + 2*sigr;
            missing(k,b) = 0;
            missing(b,k) = 0;
        end
    end
    for kr = 1: length(restrain(k).r_intern)
        if strcmp(restrain(k).r_intern(kr).type,'Gaussian')
            b = restrain(k).r_intern(kr).site;
            rmean = restrain(k).r_intern(kr).par1;
            sigr = restrain(k).r_intern(kr).par2;
            lower_bounds(k,b) = rmean - 2*sigr;
            lower_bounds(b,k) = rmean - 2*sigr;
            upper_bounds(k,b) = rmean + 2*sigr;
            upper_bounds(b,k) = rmean + 2*sigr;
            missing(k,b) = 0;
            missing(b,k) = 0;
        end
    end
end
[~,~,err] = triangle_bound_smoothing(lower_bounds,upper_bounds);
if err
    return;
end

if auxiliary
    for k1 = 1:vert-1,
        for k2 = k1+1:vert,
            setbound = false;
            if abs(lb(k1,k2) - lower_bounds(k1,k2)) > eps,
                if missing(k1,k2),
                    if lb(k1,k2) > 2,
                        setbound = true;
                    else
                        setbound = false;
                    end
                end
            end
            if abs(ub(k1,k2) - upper_bounds(k1,k2)) > eps,
                if missing(k1,k2),
                    if ub(k1,k2) < 150,
                        setbound = true;
                    end
                end
            end
            if setbound,
                if k2 > ni,
                    fprintf(1,'Setting auxiliary beacon restraint:');
                    kr = length(restrain(k1).r_beacon)+1;
                    if kr > 1,
                        label_type =  restrain(k).r_beacon(1).label_type;
                        label_T =  restrain(k).r_beacon(1).label_T;
                    else
                        label_type =  restrain(k).r_intern(1).label_type;
                        label_T =  restrain(k).r_intern(1).label_T;
                    end
                    b = k2-ni;
                    kb = bref(b,1);
                    krb = bref(b,2);
                    restrain(k1).r_beacon(kr).xyz = restrain(kb).r_beacon(krb).xyz;
                    restrain(k).r_beacon(kr).label_type = label_type;
                    restrain(k).r_beacon(kr).label_T = label_T;
                    restrain(k).r_beacon(kr).bindices = restrain(kb).r_beacon(krb).bindices;
                    restrain(k).r_beacon(kr).resb = restrain(kb).r_beacon(krb).resb;
                    restrain(k).r_beacon(kr).type = 'bounds';
                    restrain(k).r_beacon(kr).par1 = -lb(k1,k2);
                    restrain(k).r_beacon(kr).par2 = -ub(k1,k2);
                    resa = restraints.res1 + k1 - 1;
                    fprintf(1,'(%i,%i): [%4.2f,%4.2f] ?\n',resa,restrain(k).r_beacon(kr).resb,lb(k1,k2),ub(k1,k2));
                else
                    fprintf(1,'Setting auxiliary internal restraint:');
                    kr = length(restrain(k1).r_intern)+1;
                    if kr > 1,
                        label_type =  restrain(k).r_intern(1).label_type;
                        label_T =  restrain(k).r_intern(1).label_T;
                    else
                        label_type =  restrain(k).r_beacon(1).label_type;
                        label_T =  restrain(k).r_beacon(1).label_T;
                    end
                    restrain(k).r_intern(kr).site = k2;
                    restrain(k).r_intern(kr).label_type = label_type;
                    restrain(k).r_intern(kr).label_T = label_T;
                    restrain(k).r_intern(kr).resb = restraints.res1 + k2 - 1;
                    restrain(k).r_intern(kr).type = 'bounds';
                    restrain(k).r_intern(kr).par1 = -lb(k1,k2);
                    restrain(k).r_intern(kr).par2 = -ub(k1,k2);
                    resa = restraints.res1 + k1 - 1;
                    fprintf(1,'(%i,%i): [%4.2f,%4.2f] ?\n',resa,restrain(k).r_beacon(kr).resb,lb(k1,k2),ub(k1,k2));
                end
            end
        end
    end
end

min_approach = handles.min_approach;

directory = general.tmp_files;
curr_dir = get(handles.checkbox_save_cdir,'Value');
if curr_dir,
    directory = pwd;
end

display_it = get(handles.checkbox_display,'Value');
display_terminus = get(handles.checkbox_display_terminus,'Value');
attach_it = get(handles.checkbox_attach,'Value');

if handles.checkbox_std_seed.Value,
    rng(13);
elseif isfield(handles.restraints,'seed')
    if isnan(handles.restraints.seed) || handles.restraints.seed < 0,
        rng(13);
        add_msg_board(sprintf('Warning: Invalid random number seed value %i specified in restraint file',handles.restraints.seed));
        add_msg_board('Reverting to standard seed value of 13');
    else
        rng(handles.restraints.seed);
        add_msg_board(sprintf('Random number seed value of %i is set.',handles.restraints.seed));
    end
else
    rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs
end

pmodel = handles.p_model;
pthr = exp(-erfinv(pmodel)^2);

max_models = handles.ensemble;

all_p_model = zeros(1,max_models);

snum_vec = zeros(1,max_models);

ntrials = handles.max_trials; % number of Monte Carlo trials
max_seconds = 3600*handles.max_time; % maximum runtime in seconds

time_left = sprintf('%i h %i min estimated run time to completion',floor(max_seconds/3600),round((max_seconds-3600*floor(max_seconds/3600))/60));
set(handles.text_time_left,'String',time_left);
set(handles.text_time_left,'ForegroundColor','r');

Rama_secondary = get(handles.radiobutton_Rama_secondary,'Value');
if Rama_secondary,
    load([general.Ramachandran 'Ramachandran']);
else
    load([general.Ramachandran 'Ramachandran_disordered']);
end

closed_loop = false;
free_standing = false;
set(handles.text_loop_closure_failure,'String','n.a.');
set(handles.text_Rama_mismatch,'String','n.a.');

sequence = restraints.sequence;
if isempty(restraints.anchorN) && ~isempty(restraints.anchorC), % reverse model
    reverse = true;
    template_indices = restraints.Ca_indices(1:3);
%     if isempty(restraints.anchorC), % no anchor, currently not implemented
%         add_msg_board('ERROR: Domain modelling without anchor residues is not yet implemented. No action.');
%         return;
%     end
    sequence = [sequence restraints.Cseq(1)];
elseif ~isempty(restraints.anchorN) % loop anchored at N terminus
    reverse = false;
    template_indices = restraints.Na_indices(1:3);
    sequence = [restraints.Nseq(2) sequence];
    if ~isempty(restraints.anchorC), % closed loop
        sequence = [sequence restraints.Cseq];
        closed_loop = true;
        if attach_it,
            if sum(abs(restraints.Na_indices(1:3)-restraints.Ca_indices(1:3))) > 0,
                add_msg_board('Warning: N-terminal and C-terminal anchors are in different chains or chain models.');
                add_msg_board('Loop models cannot be attached.');
                attach_it = false;
            end
        end
    end
else % free-standing loop
    reverse = false;
    closed_loop = false;
    free_standing = true;
    template_indices = [];
    if attach_it,
        add_msg_board('Warning: No anchor residues defined. Loop is not attached, but models are collected into a single structure.');
    end
%     attach_it = false;
end

if attach_it && ~free_standing, % make template copies in current structure
    for modnum = 2:max_models,
        copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
    end
end


res1 = restraints.res1;
resend = restraints.rese;
terminus = resend;
N_terminus = res1;
if reverse,
    terminus = res1;
end

if ~isempty(restraints.anchorN),
    display_start = res1-1;
else
    display_start = res1;
end

if ~isempty(restraints.anchorC),
    display_end = resend+1;
else
    display_end = resend;
end

if ~free_standing,
    adr=mk_address(restraints.template);
    handles.fname_bas = strcat(adr,'_loop');
else
    handles.fname_bas = 'free_loop';
end

fname = sprintf('%s_%i_%i',handles.fname_bas,res1,resend);

min_prob = pthr^handles.n_restraints;

add_msg_board(sprintf('Cumulative probability threshold for %i restraints is %6.4f.',handles.n_restraints,min_prob));

if ~free_standing,
    [msg,prot_coor] = get_structure(restraints.template,'xyz');
    if msg.error,
        add_msg_board('ERROR: Coordinates of template structure could not be retrieved.');
    end
else
    prot_coor = [];
end

success = 0;
err_count=zeros(1,11);
Ram_fixed = 0;
Ram_fix_clash = 0;
resax = res1:resend;
res_stat = zeros(1,length(resax));
set(gcf,'Pointer','watch');
distributions = cell(handles.n_restraints+handles.b_restraints);
restraint_distr = cell(handles.n_restraints);
descriptors = cell(handles.n_restraints);
monitor_distr = cell(handles.n_monitor);
monitor_descr = cell(handles.n_monitor);
rax = get_distribution;
for k = 1:handles.n_restraints + handles.b_restraints
    distributions{k} = zeros(1,length(rax));
    restraint_distr{k} = zeros(1,length(rax));
end
for k = 1:handles.n_monitor
    monitor_distr{k} = zeros(1,length(rax));
end
for k = 1:handles.n_restraints + handles.n_monitor
    handles.exp_distr{k} = [];
    handles.exp_r{k} = [];
end

fullname = fullfile(handles.save_path,sprintf('%s_%s.log',handles.save_name,datestr(now,30)));
fid_report = fopen(fullname,'wt');

mname_stub = fullfile(handles.save_path,handles.save_name);
add_msg_board(sprintf('Saving models to %s_m#.pdb\n',mname_stub));

p_coor = cell(parnum);
p_restrain = cell(parnum);
p_errcode = zeros(1,parnum);
p_cumprob = zeros(1,parnum);
p_kres = zeros(1,parnum);
drawnow;
tic;
kMC = 1;
p_anchorC = restraints.anchorC;
p_anchorCn = restraints.anchorCn;
p_anchorN = restraints.anchorN;
p_anchorNp = restraints.anchorNp;
n_restraints = handles.n_restraints;

rescodes = zeros(1,length(sequence));
for k = 1:length(sequence),
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end

Rama_res.me = Ramachandran.me;
Rama_res.ephi = Ramachandran.ephi;
Rama_res.epsi = Ramachandran.epsi;
Rama_res.allowed_P = Ramachandran.allowed_P;
Rama_res.allowed_G = Ramachandran.allowed_G;
Rama_res.allowed_gen = Ramachandran.allowed_gen;

while kMC <= ntrials
    parfor kp = 1:parnum % parfor
        if reverse
            [coor,errcode,restrain1,cumprob,kres] = mk_loop_model_reverse(sequence, p_anchorC, p_anchorCn, prot_coor, restrain, Rama_res, rescodes, min_prob,n_restraints);
        else
            [coor,errcode,restrain1,cumprob,kres] = mk_loop_model(sequence, p_anchorN, p_anchorC, p_anchorNp, p_anchorCn, prot_coor, restrain, Rama_res, rescodes, min_prob, n_restraints);
            kres = kres-1;
        end
        p_coor{kp} = coor; 
        p_errcode(kp) = errcode;
        p_restrain{kp} = restrain1;
        p_cumprob(kp) = cumprob;
        p_kres(kp) = kres;
    end
    for kp = 1:parnum
        kMC = kMC + 1;
        coor = p_coor{kp};
        errcode = p_errcode(kp);
        restrain1 = p_restrain{kp};
        cumprob = p_cumprob(kp);
        kres = p_kres(kp);
%         for k = 1:length(restrain),
%             restrain(k).secondary = secvec(k);
%         end
        res_stat(kres) = res_stat(kres)+1;
        runtime = toc;
        if errcode == -1
            Ram_fixed = Ram_fixed + 1;
            errcode = 0;
        end
        if errcode == -4
            Ram_fixed = Ram_fixed + 1;
            Ram_fix_clash = Ram_fix_clash + 1;
            errcode = 4;
        end
        if handles.n_restraints == 0
            p_model = 1;
        else
            p_model = erf(sqrt(-log(cumprob)/handles.n_restraints));
        end
        err_count(errcode+1) = err_count(errcode+1) + 1;

        if ~errcode
            tpm = runtime/err_count(1);
            set(handles.text_time_per_model,'String',sprintf('%8.1f',tpm));
            success = success + 1;
            if success == 1
                bb0 = coor;
            elseif  isempty(restraints.anchorN) && isempty(restraints.anchorC),
                [rms,coor] = rmsd_superimpose(bb0,coor);
                add_msg_board(sprintf('Model superimposes onto first model with rmsd of %4.1f ?',rms));
            end
            loopname = write_pdb_backbone(coor,restraints.sequence,fname,success,res1,directory);
            [pmodel,status,result] = make_SCWRL4_sidegroups(loopname,directory);
    %         fprintf(1,'Model %s: ',pmodel);
    %         if ~status,
    %             fprintf(1,'successful.\n');
    %         else
    %             fprintf(1,'Error in sidegroup creation.\n');
    %         end
            [pclash,iclash] = check_decorated_loop(pmodel,prot_coor,res1,resend,min_approach,directory);

            if pclash
    %             fprintf(1,'Loop with sidegroups clashes with protein\n');
                err_count(10) = err_count(10) + 1;
                success = success - 1;
            elseif iclash
    %             fprintf(1,'Loop with sidegroups clashes with itself\n');
                err_count(11) = err_count(11) + 1;
                success = success - 1; 
            else
                all_p_model(success) = p_model;
                if attach_it,
                    if free_standing && success == 1
                        [msg,snum]=add_pdb(pmodel);
                        template_indices = [snum 1 1];
                        model_indices = template_indices;
                        for modnum = 2:max_models
                            copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
                        end
                    else
                        structure=rd_pdb(pmodel);
                        model_indices = template_indices;
                        model_indices(3) = success;
                        if free_standing
                            replace_model(model_indices,structure);
                        else
                            insert_residues(res1,resend,structure,model_indices,Nanchor_res,Canchor_res);
                        end
                        if isfield(model,'selected')
                            model = rmfield(model,'selected');
                        end
                        model.selected{1} = model_indices;
                        mname = sprintf('%s_m%i',mname_stub,success);
                        wr_pdb_selected(mname,'DENS');
                    end
                    if display_it
                        command=sprintf('show {%i}%i-%i coil',success,display_start,display_end);
                        hMain.store_undo=false;
                        cmd(hMain,command);
                        command=sprintf('color {%i}%i-%i crimson',success,display_start,display_end);
                        hMain.store_undo=false;
                        cmd(hMain,command);
                        command=sprintf('transparency {%i}%i-%i %6.3f',success,display_start,display_end,p_model);
                        hMain.store_undo=false;
                        cmd(hMain,command);
                    end
                    if display_terminus,
                        command=sprintf('show {%i}%i.CA space-filling',success,terminus);
                        hMain.store_undo=false;
                        cmd(hMain,command);
                        command=sprintf('color {%i}%i.CA blue',success,terminus);
                        hMain.store_undo=false;
                        cmd(hMain,command);                    
                        command=sprintf('transparency {%i}%i.CA %6.3f',success,terminus,p_model);
                        hMain.store_undo=false;
                        cmd(hMain,command);      
                        if free_standing,
                            command=sprintf('show {%i}%i.CA space-filling',success,N_terminus);
                            hMain.store_undo=false;
                            cmd(hMain,command);
                            command=sprintf('color {%i}%i.CA darkgreen',success,N_terminus);
                            hMain.store_undo=false;
                            cmd(hMain,command);                    
                            command=sprintf('transparency {%i}%i.CA %6.3f',success,N_terminus,p_model);
                            hMain.store_undo=false;
                            cmd(hMain,command);      
                        end
                    end
                    [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions(fid_report,model_indices,restrain1,handles.monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1);
                elseif display_it || display_terminus,               
                    [msg,snum]=add_pdb(pmodel);
                    snum_vec(success) = snum;
                    model_indices = ones(1,3);
                    model_indices(1) = snum;
                    [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions_free(fid_report,model_indices,restrain1,handles.monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1);
                    if msg.error,
                        add_msg_board(sprintf('ERROR: Conformation %s could not be loaded for display.',pmodel));
                    else
                        if display_it,
                            command=sprintf('show [%i] coil',snum);
                            hMain.store_undo=false;
                            cmd(hMain,command);
                            command=sprintf('color [%i] crimson',snum);
                            hMain.store_undo=false;
                            cmd(hMain,command);
                        end
                        if display_terminus,
                            command=sprintf('show [%i]%i.CA space-filling',snum,terminus);
                            hMain.store_undo=false;
                            cmd(hMain,command);
                            command=sprintf('color [%i]%i.CA blue',snum,terminus);
                            hMain.store_undo=false;
                            cmd(hMain,command);
                            if  free_standing,
                                command=sprintf('show [%i]%i.CA space-filling',snum,N_terminus);
                                hMain.store_undo=false;
                                cmd(hMain,command);
                                command=sprintf('color [%i]%i.CA darkgreen',snum,N_terminus);
                                hMain.store_undo=false;
                                cmd(hMain,command);
                            end
                        end
                        command=sprintf('transparency [%i] %6.3f',snum,p_model);
                        hMain.store_undo=false;
                        cmd(hMain,command);
                    end
                end
            end
        end
        if mod(kMC,disp_update) == 0,
            ftr = (1 - kMC/ntrials)*max_seconds;
            fti = max_seconds - runtime;
            fmo = runtime*(max_models-success)/success;
            time_left = min([fti fmo ftr]);
            hours = floor(time_left/3600);
            minutes = round((time_left-3600*floor(time_left/3600))/60);
            if minutes == 60,
                hours = hours + 1;
                minutes = 0;
            end
            time_left = sprintf('%i h %i min estimated run time to completion',hours,minutes); 
            set(handles.text_time_left,'String',time_left);
            axes(handles.axes_success_length);
            success_distr = (1-success/kMC)*res_stat/sum(res_stat);
            if reverse,
                success_distr = fliplr(success_distr);
            end
            normalize = 1- cumsum(success_distr);
            normalize(normalize==0)=1;
            for ks = 2:length(success_distr),
                success_distr(ks) = success_distr(ks)/normalize(ks-1);
            end
            if reverse,
                success_distr = fliplr(success_distr);
            end
            plot(resax,success_distr);
            set(gca,'FontSize',8);
            xlabel('Residue');
            ylabel('Rejection probability');
            set(handles.text_percent_time,'String',sprintf('%5.2f',100*runtime/max_seconds));
            set(handles.text_max_trials,'String',sprintf('%5.2f',100*kMC/ntrials));
            set(handles.text_success,'String',sprintf('%i',success));
            set(handles.text_backbone_success,'String',sprintf('%i',err_count(1)));
            set(handles.text_restraint_violations,'String',sprintf('%5.2f',100*err_count(6)/kMC));
            set(handles.text_loop_clash,'String',sprintf('%5.2f',100*(err_count(3)+err_count(8))/kMC));
            set(handles.text_protein_clashes,'String',sprintf('%5.2f',100*(err_count(5)+err_count(7))/kMC));
            set(handles.text_sidegroup_clashes_loop,'String',sprintf('%5.2f',100*err_count(11)/err_count(1)));
            set(handles.text_sidegroup_clashes_protein,'String',sprintf('%5.2f',100*err_count(10)/err_count(1)));
            if closed_loop,
                set(handles.text_loop_closure_failure,'String',sprintf('%5.2f',100*err_count(2)/kMC));
                set(handles.text_Rama_mismatch,'String',sprintf('%5.2f',100*err_count(4)/kMC));
            end
            drawnow;
        end
        if success >= max_models,
            break
        end
    end
    if success >= max_models,
        break
    end
%     if mod(kMC,10000)
%         fprintf(1,'%6.3f ms per trial\n',1000*runtime/kMC);
%     end
    if runtime >= max_seconds,
        add_msg_board('Warning: Ensemble generation stopped as maximum allotted time was exceeded.');
        break
    end
end

runtime = toc;
hours = floor(runtime/3600);
minutes = round((runtime-3600*floor(runtime/3600))/60); 

fprintf(fid_report,'\n--- Final statistics ---\n\n');
fprintf(fid_report,'Total runtime: %i h %i min.\n',hours,minutes);
fprintf(fid_report,'%5.2f%% of maximum number of MC trials spent.\n',100*kMC/ntrials);
fprintf(fid_report,'Ensemble with %i models generated.\n',success);
fprintf(fid_report,'%i backbone models were generated.\n',err_count(1));
fprintf(fid_report,'%5.2f%% of all trials had restraint violations.\n',100*err_count(6)/kMC);
fprintf(fid_report,'This is a success rate of %6.1f ppm.\n',1e6*(1-err_count(6)/kMC));
fprintf(fid_report,'%5.2f%% of all trials had internal loop backbone clashes.\n',100*(err_count(3)+err_count(8))/kMC);
fprintf(fid_report,'%5.2f%% of all trials hand backbone clashes with the protein.\n',100*(err_count(5)+err_count(7))/kMC);
fprintf(fid_report,'The sidegroup clash threshold was %4.1f ?.\n',min_approach);
fprintf(fid_report,'%5.2f%% of the backbone models lead to internal sidegroup clashes.\n',100*err_count(11)/err_count(1));
fprintf(fid_report,'%5.2f%% of the backbone models lead to sidegroup clashes with the protein.\n',100*err_count(10)/err_count(1));
if closed_loop,
    fprintf('%5.2f%% of all trials lead to loops that could not be closed.\n',100*err_count(2)/kMC);
    fprintf('%5.2f%% of all trials lead to loops for which the link Ramachandran angles could not be fixed.\n',100*err_count(4)/kMC);
end

hours = floor(max_seconds/3600);
minutes = round((max_seconds-3600*floor(max_seconds/3600))/60);

fprintf(fid_report,'\nThe run had a time limit of: %i h %i min.\n',hours,minutes);
fprintf(fid_report,'The run had a limit of %i Monte Carlo trials.\n',ntrials);
fprintf(fid_report,'%i models had been requested with an ensemble probability of %4.2f.\n',max_models,handles.p_model);
fprintf(fid_report,'The number of Gaussian-type restraints was %i.\n',handles.n_restraints);
Rama_type = 'without';
if Rama_secondary,
    Rama_type = 'with';
end
fprintf(fid_report,'Ramachandran statistics %s secondary structure elements was used.\n\n',Rama_type);

hours = floor(runtime/3600);
minutes = round((runtime-3600*floor(runtime/3600))/60);
 
set(handles.text_time_left,'String',sprintf('Ready after %i h %i min.',hours,minutes));
if success == max_models,
    set(handles.text_time_left,'ForegroundColor',[0 0.5 0]);
else
    if attach_it && exist('model_indices','var'),
        delete_empty_models(model_indices(1),success);
    end
    set(handles.text_time_left,'ForegroundColor',[0.80 0.50 0]);    
end


handles.resax = resax;
success_distr = (1-success/kMC)*res_stat/sum(res_stat);
if reverse,
    success_distr = fliplr(success_distr);
end
normalize = 1- cumsum(success_distr);
normalize(normalize==0)=1;
for ks = 2:length(success_distr),
    success_distr(ks) = success_distr(ks)/normalize(ks-1);
end
if reverse,
    success_distr = fliplr(success_distr);
end
handles.rax = 10*rax;

handles.success_distr = success_distr;
handles.distributions = distributions;
handles.restraint_distr = restraint_distr;
handles.descriptors = descriptors;
handles.monitor_distr = monitor_distr;
handles.monitor_descr = monitor_descr;
if handles.n_restraints > 0 || handles.n_monitor > 0,
    handles.plot_nr = 1;
else
    handles.plot_nr = 0;
end
handles.copy = false; 
handles.snum_vec = snum_vec;

fprintf(fid_report,'--- Distribution shifts and overlaps ---\n\n');
fprintf(fid_report,'Shift [?]\tOverlap\n');
for kd = 1:handles.n_restraints
    [overlap,shift] = get_overlap(handles.rax,handles.distributions{kd},handles.restraint_distr{kd});
    fprintf(fid_report,'%5.2f\t%5.3f\n',shift,overlap);
end
fclose(fid_report);


update_plot(handles);


if kMC == ntrials,
    add_msg_board('Warning: Ensemble generation stopped as maximum allotted number of Monte Carlo trials was exceeded.');
end
set(handles.text_percent_time,'String',sprintf('%5.2f',100*runtime/max_seconds));
set(handles.text_max_trials,'String',sprintf('%5.2f',100*kMC/ntrials));
set(handles.text_success,'String',sprintf('%i',success));
set(handles.text_backbone_success,'String',sprintf('%i',err_count(1)));
set(handles.text_restraint_violations,'String',sprintf('%5.2f',100*err_count(6)/kMC));
set(handles.text_loop_clash,'String',sprintf('%5.2f',100*(err_count(3)+err_count(8))/kMC));
set(handles.text_protein_clashes,'String',sprintf('%5.2f',100*(err_count(5)+err_count(7))/kMC));
set(handles.text_sidegroup_clashes_loop,'String',sprintf('%5.2f',100*err_count(11)/err_count(1)));
set(handles.text_sidegroup_clashes_protein,'String',sprintf('%5.2f',100*err_count(10)/err_count(1)));
if closed_loop,
    set(handles.text_loop_closure_failure,'String',sprintf('%5.2f',100*err_count(2)/kMC));
    set(handles.text_Rama_mismatch,'String',sprintf('%5.2f',100*err_count(4)/kMC));
end
drawnow;
if attach_it && success > 0,
    consolidate_chain(template_indices(1:2));
end

if attach_it,
    set(handles.pushbutton_run,'Enable','off');
end

set(gcf,'Pointer','arrow');
drawnow;

guidata(hObject,handles);

function labels = get_labels(llist,label,T,handles)

global model
global hMain

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end

% check whether sites are already labelled and whether all restraint sites
% do exist
lindices=zeros(length(labels),4);
for k=1:length(labels),
    cindices=labels(k).indices;
    if ~isempty(cindices),
        lindices(k,:)=cindices;
    end
end
poi=0;
to_do_list{1}=' ';
for k=1:length(llist),
    adr1=llist{k};
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
        add_msg_board(sprintf('This site does not exist in current structure %s',mk_address(1)));
    end
    found=false;
    for l=1:length(labels),
        diff=ind1-lindices(l,:);
        if sum(abs(diff))==0,
            found=true;
        end
    end
    if ~found,
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}),
                found=true;
            end
        end
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end
    end
end

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label{k},T(k));
        hMain.store_undo=false;
        hMain.dynamic_rotamers=false;
        cmd(hMain,command);
    end
end

if isfield(model,'sites'),
    labels = label_information(model.sites);
else
    labels = cell(0);
end

function NO = get_relative_label(libname)

load(libname);
pops = rot_lib.calibration.pop;
pops = pops/sum(pops);
if isfield(rot_lib.usefull_atoms,'midNO')
    midNO = rot_lib.usefull_atoms.midNO;
    NO = zeros(1,3);
    for k = 1:length(rot_lib.library),
        coor = rot_lib.library(k).ecoor;
        NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
    end
elseif isfield(rot_lib,'spin_density')
    NO = zeros(1,3);
    for k = 1:length(rot_lib.library),
        coor = rot_lib.library(k).ecoor;
        coor = coor(rot_lib.spin_density(:,1),2:4);
        coor = rot_lib.spin_density(:,2)'*coor;
        NO = NO + pops(k)*coor;
    end
end
NO = NO/(sum(pops));

function spin_coor = get_relative_bilabel(libname,V_ksi)

gas_un = 8.314472;    % universal gas constant in CI (J/(mol*K)       
T = 298;

load(libname);
pops = zeros(1,length(ligand_lib.library));
all_coor = zeros(length(ligand_lib.library),3);
if isfield(ligand_lib,'spin_density')
    spin_density = ligand_lib.spin_density;
else
    spin_density = [1 1];
end
sdp = spin_density(:,2)';
for k = 1:length(ligand_lib.library)
    ksi = pi*ligand_lib.library(k).dihedrals/pi;
    int_en = V_ksi*(1-cos(ksi))/2;
    pops(k) = exp(-int_en/(gas_un*T));
    all_coor(k,:) = sdp*ligand_lib.library(k).ecoor(spin_density(:,1),2:4);
end
spin_coor = pops*all_coor;
spin_coor = spin_coor/(sum(pops));

function labels=label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            labels(poi).indices=sites{k0}(k1).residue(k).indices;
            id=tag2id(sites{k0}(k1).residue(k).label,label_defs.restags);
            labels(poi).name=label_defs.residues(id).short_name;
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd=NOpos_rmsd(NOpos);
        end
    end
end

function [rmsd,xyz]=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
xyz = [xmean,ymean,zmean];
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for ? -> nm

function [res,chain_model] = separate_address(adr)

spoi = strfind(adr,'|');
if isempty(spoi)
    cadr{1} = adr;
    res = 0;
else
    cadr{1} = adr(1:spoi-1);
    cadr{2} = adr(spoi+1:end);
    res = [0,0];
end
for k = 1:length(cadr)
    adr = cadr{k};
    poi = 0;
    poi1 = strfind(adr,']');
    if ~isempty(poi1),
        poi = poi1;
    end
    poi1 = strfind(adr,')');
    if ~isempty(poi1),
        poi = poi1;
    end
    poi1 = strfind(adr,'}');
    if ~isempty(poi1),
        poi = poi1;
    end
    if poi > 0,
        chain_model = adr(1:poi);
    else
        chain_model = '';
    end
    resstr = adr(poi+1:end);
    cres = str2double(resstr);
    if isnan(cres),
        res = [];
        return
    end
    if cres - floor(cres) > eps,
        res = [];
        return
    end
    res(k) = cres;
end

function [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,rmean,sigr,res1,number,bnumber,label,T,bindices,resb)

scale_units = 10;

grace = 0.5; % 5 ? uncertainty of label position

k = res_loop - res1 + 1;
kr = length(restrain(k).r_beacon)+1;
restrain(k).label = NO;
restrain(k).r_beacon(kr).xyz = xyz_beacon;
restrain(k).r_beacon(kr).label_type = label;
restrain(k).r_beacon(kr).label_T = T;
restrain(k).r_beacon(kr).bindices = bindices;
restrain(k).r_beacon(kr).resb = resb;
if rmean > 0 && sigr > 0
    restrain(k).r_beacon(kr).type = 'Gaussian';
    restrain(k).r_beacon(kr).par1 = rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = sqrt(sigr^2 + grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_beacon(kr).type = 'bounds';
    restrain(k).r_beacon(kr).par1 = -rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,rmean,sigr,res1,number,bnumber,label,T)

scale_units = 10;

grace = 0.5; % 5 ? uncertainty of label position

if resa < resb, % restraint must be stored at the later site
    exch = resa;
    resa = resb;
    resb = exch;
    exch = NO1;
    NO1 = NO2;
    NO2 = exch;
end
k = resa - res1 + 1;
k2 = resb - res1 + 1;
kr = length(restrain(k).r_intern)+1;
restrain(k).label = NO1;
restrain(k2).label = NO2;
restrain(k).r_intern(kr).site = k2;
restrain(k).r_intern(kr).label_type = label;
restrain(k).r_intern(kr).label_T = T;
restrain(k).r_intern(kr).resb = resb;
if rmean > 0 && sigr > 0,
    restrain(k).r_intern(kr).type = 'Gaussian';
    restrain(k).r_intern(kr).par1 = rmean*scale_units;
    restrain(k).r_intern(kr).par2 = sqrt(sigr^2 + 2*grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_intern(kr).type = 'bounds';
    restrain(k).r_intern(kr).par1 = -rmean*scale_units;
    restrain(k).r_intern(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_depth_restraint(restrain,NO,res,rmean,sigr,res1,number,bnumber,label,T)

scale_units = 10;

k = res - res1 + 1;
kr = length(restrain(k).depth)+1;
restrain(k).depth(kr).label_type = label;
restrain(k).depth(kr).label_T = T;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).depth(kr).site = 'label';
else
    restrain(k).depth(kr).site = 'CA';
end
if rmean > 0 && sigr > 0,
    restrain(k).depth(kr).type = 'Gaussian';
    restrain(k).depth(kr).par1 = rmean*scale_units;
    restrain(k).depth(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).depth(kr).type = 'bounds';
    restrain(k).depth(kr).par1 = -rmean*scale_units;
    restrain(k).depth(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_oligomer_restraint(restrain,NO,res,n,rmean,sigr,res1,number,bnumber,label,T)

scale_units= 10;

k = res(end) - res1 + 1;
kr = length(restrain(k).oligomer)+1;
restrain(k).oligomer(kr).res0 = res(1)  - res1 + 1;
restrain(k).oligomer(kr).label_type = label;
restrain(k).oligomer(kr).label_T = T;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).oligomer(kr).site = 'label';
else
    restrain(k).oligomer(kr).site = 'CA';
end
if rmean > 0 && sigr > 0,
    restrain(k).oligomer(kr).type = 'Gaussian';
    restrain(k).oligomer(kr).par1 = rmean*scale_units;
    restrain(k).oligomer(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).oligomer(kr).type = 'bounds';
    restrain(k).oligomer(kr).par1 = -rmean*scale_units;
    restrain(k).oligomer(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end
restrain(k).oligomer(kr).n = n;



function edit_max_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_time as text
%        str2double(get(hObject,'String')) returns contents of edit_max_time as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,1000,2,'%5.2f',0);
handles.max_time = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_max_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ensemble as text
%        str2double(get(hObject,'String')) returns contents of edit_ensemble as a double

[v,handles]=edit_update_MMM(handles,hObject,1,10000,20,'%i',1);
handles.ensemble = v;
handles.ensemble_user = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_ensemble_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_trials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_trials as text
%        str2double(get(hObject,'String')) returns contents of edit_max_trials as a double

[v,handles]=edit_update_MMM(handles,hObject,100000,100000000000,5000000000,'%i',1);
handles.max_trials = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_max_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p_model_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p_model as text
%        str2double(get(hObject,'String')) returns contents of edit_p_model as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,0.95,0.5,'%5.2f',0);
handles.p_model = v;
handles.p_model_user = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_p_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_adaptive.
function checkbox_adaptive_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_adaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_adaptive


% --- Executes on button press in checkbox_save_cdir.
function checkbox_save_cdir_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_save_cdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_save_cdir


% --- Executes on button press in checkbox_display_terminus.
function checkbox_display_terminus_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_terminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_terminus


% --- Executes on button press in checkbox_attach.
function checkbox_attach_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_attach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_attach

function loopname = write_pdb_backbone(coor,sequence,fname0,model,res1,directory)

my_dir = pwd;
cd(directory);

residues='ALAARGASNASPCYSGLNGLUGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL';
oneletter='ARNDCQEGHILKMFPSTWYV';

fname = sprintf('%s_m%i',fname0,model);

loopname = [fname '.pdb'];
wfile=fopen(loopname,'w');
for k = 1:length(sequence)
    respoi=strfind(oneletter,sequence(k));
    residue=residues(1+3*(respoi-1):3+3*(respoi-1));
    N = coor(4*k-3,:);
    CA = coor(4*k-2,:);
    C = coor(4*k-1,:);
    O = coor(4*k,:);
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           N\n','ATOM  ',4*k-3,'  N   ',residue,k+res1-1,N(1),N(2),N(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-2,'  CA  ',residue,k+res1-1,CA(1),CA(2),CA(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-1,'  C   ',residue,k+res1-1,C(1),C(2),C(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           O\n','ATOM  ',4*k,'  O   ',residue,k+res1-1,O(1),O(2),O(3));
end
fclose(wfile);

cd(my_dir);

function [outname,status,result] = make_SCWRL4_sidegroups(inname,directory)
%
% Attaches or corrects sidegroups using SCWRL4
% the program SCWRL4 needs to be on the current Matlab path
%
% inname   name of the PDB file to which sidegroups should be attached
% outname  name of the output PDB file with SCWRL4 sidegroups
%

my_dir = pwd;
cd(directory);

poi = strfind(inname,'.pdb');
outname = [inname(1:poi-1) '_SCWRL4.pdb'];

s=which('scwrl4.exe');
cmd=[s ' -i ' inname ' -o ' outname];
[status,result]=dos(cmd);

cd(my_dir);

function [pclash,iclash,approach_prot,approach_loop] = check_decorated_loop(loopname,prot_coor,res1,resend,min_approach,directory)

my_dir = pwd;
cd(directory);

if ~exist('min_approach','var')
    min_approach = 1.2; 
end

approach_prot = -1;
approach_loop = -1;
pclash = 1;
iclash = 1;
loop_coor = zeros(5000,3);
l_res_assign = zeros(5000,3);
fid=fopen(loopname);
if fid==-1
    add_msg_board(sprintf('Warning: Loop structure PDB file %s missing. Rejected.',loopname));
    cd(my_dir);
    return;
end
poi = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if length(tline) >= 6
        record=tline(1:6);
        resnum = str2double(tline(23:26));
        if strcmpi(record,'ATOM  ') || strcmpi(record,'HETATM')
            if length(tline) < 78 || tline(78)~='H'
                if ~strcmpi(tline(14),'H') && resnum ~= res1 && resnum ~= resend
                    poi = poi +1;
                    l_res_assign(poi) = resnum;
                    valstr = [tline(31:38) ' ' tline(39:46) ' ' tline(47:54)];
                    loop_coor(poi,:) = str2num(valstr);
                end
            end
        end
    end
end
fclose(fid);

cd(my_dir);

loop_coor = loop_coor(1:poi,:);
l_res_assign = l_res_assign(1:poi);

pclash = 0;
iclash = 0;

[m1,~] = size(loop_coor); % get sizes of the coordinates arrays
[m2,~] = size(prot_coor);

if m2 > 0
    a2 = repmat(sum(loop_coor.^2,2),1,m2);
    b2 = repmat(sum(prot_coor.^2,2),1,m1).';
    pair_dist = sqrt(abs(a2 + b2 - 2*loop_coor*prot_coor.'));
    min_dist = min(min(pair_dist));

    approach_prot = min_dist;
    if min_dist < min_approach
    %    fprintf(2,'Minimum sidegroup distance to protein is %6.2f ?\n',min_dist);
       pclash = 1;
       cd(my_dir);
       return
    end
end

min_dist = 1e6;
% test for minimum distance within loop
we_clash = zeros(1,2);
for k1 = 1:poi-1
    for k2 = k1+1:poi
        if abs(l_res_assign(k1)-l_res_assign(k2))>1
            approach = norm(loop_coor(k1,:) - loop_coor(k2,:));
            if approach < min_dist
                min_dist = approach;
                we_clash = [k1,k2];
            end
        end
    end
end
approach_loop = min_dist;
if min_dist < min_approach
%    fprintf(2,'Minimum heavy-atom distance: %6.2f ? at (%i,%i)\n',min_dist,we_clash);
    iclash = 1;
end


% --- Executes on button press in checkbox_display.
function checkbox_display_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display



function edit_min_approach_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_approach as text
%        str2double(get(hObject,'String')) returns contents of edit_min_approach as a double

[v,handles]=edit_update_MMM(handles,hObject,0.5,2.5,2.0,'%5.2f',0);
handles.min_approach = v;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_min_approach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function insert_residues(res1,rese,structure,template_indices,Nanchor_res,Canchor_res)

% Bfactor is set to twice the maximum Bfactor of existing residues

global model
global residue_defs


snum = template_indices(1);
cnum = template_indices(2);
mnum = template_indices(3);

if ~isempty(Nanchor_res),
    model.structures{snum}(cnum).residues{mnum}.info(Nanchor_res).terminal = 0;
end

Bmax=max(model.info{snum}.B_range);

newresidues=res1:rese;
[ma,na]=size(model.structures{snum}(cnum).xyz{mnum});
[ma0,maxconn]=size(model.structures{snum}(cnum).conn);
newiso=[model.structures{snum}(cnum).isotopes; zeros(20*length(newresidues),2,'single')];
newxyz=[model.structures{snum}(cnum).xyz{mnum}; zeros(20*length(newresidues),3)];
newBfac=[model.structures{snum}(cnum).Bfactor{mnum}, zeros(1,20*length(newresidues))];
newBtens=[model.structures{snum}(cnum).Btensor{mnum}; zeros(20*length(newresidues),6,'int32')];
sequence=model.structures{snum}(cnum).sequence;
restags=model.structures{snum}(cnum).residues{mnum}.residue_tags;

newrnum=zeros(1,length(newresidues));
nr=length(model.structures{snum}(cnum).residues{mnum}.info);
rpoi=0;
for k=1:length(newresidues),
    rnum0=newresidues(k);
    rnum=0;
    for kk=1:length(structure(1).residues{1}.info),
        if structure(1).residues{1}.info(kk).number == rnum0,
            rnum=kk;
            break;
        end
    end
    if rnum==0,
        add_msg_board(sprintf('Warning: Residue %i not found in modelled loop',rnum0));
        continue;
    end
    tag=structure(1).residues{1}.info(rnum).name;
    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
    sequence(rnum0)=id;
    if mnum==1 && isfield(model.structures{snum}(cnum),'seqexist'),
        model.structures{snum}(cnum).seqexist(rnum0)=1;
    end
    nr=nr+1;
    rpoi=rpoi+1;
    newrnum(rpoi)=nr;
    % model.structures{snum}(cnum).residues{mnum}.info(nr)=structure(1).residues{1}.info(rnum);
    restag=sprintf('%i:',rnum0);
    restags=strcat(restags,restag);
    pointers=structure(1).residues{1}.info(rnum).atom_numbers;
    if isempty(Canchor_res) && k == length(newresidues),
        model.structures{snum}(cnum).residues{mnum}.info(nr).terminal = 1;
    end
    for anum=1:length(pointers), % loop over atoms
        pointer=structure(1).residues{1}.info(rnum).atom_numbers{anum};
        [loc,n]=size(pointer);
        for lnum=1:loc, % loop over locations
            poi=pointer(lnum,1); % actual coordinate set number
            ma=ma+1;
            pointer(lnum,1)=ma;
            newiso(ma,:)=structure(1).isotopes(poi,:);
            newxyz(ma,:)=structure(1).xyz{1}(poi,:);
            newBfac(ma)=2*Bmax;
            newBtens(ma,:)=structure(1).Btensor{1}(poi,:);
        end
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_numbers{anum}=pointer;
        model.structures{snum}(cnum).residues{mnum}.info(nr).name=tag;
        model.structures{snum}(cnum).residues{mnum}.info(nr).type=structure(1).residues{1}.info(rnum).type;
        model.structures{snum}(cnum).residues{mnum}.info(nr).secondary=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).hetflag=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).connected=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).number=structure(1).residues{1}.info(rnum).number;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_tags=structure(1).residues{1}.info(rnum).atom_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).elements=structure(1).residues{1}.info(rnum).elements;
        model.structures{snum}(cnum).residues{mnum}.info(nr).location_tags=structure(1).residues{1}.info(rnum).location_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).insertion_code=structure(1).residues{1}.info(rnum).insertion_code;
    end
end
newiso=newiso(1:ma,:);
newxyz=newxyz(1:ma,:);
newBfac=newBfac(1:ma);
newBtens=newBtens(1:ma,:);

if mnum==1,
    model.structures{snum}(cnum).sequence=sequence;
    model.structures{snum}(cnum).isotopes=newiso;
    model.structures{snum}(cnum).conn=[model.structures{snum}(cnum).conn; zeros(ma-ma0,maxconn)];
end

model.structures{snum}(cnum).xyz{mnum}=newxyz;
model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
model.structures{snum}(cnum).Btensor{mnum}=newBtens;
model.structures{snum}(cnum).residues{mnum}.residue_tags=restags;

% make internal bonds in new residues
for k=1:length(newrnum),
    if newrnum(k) > 0,
        model.structures{snum}(cnum)=mk_internal_bonds(model.structures{snum}(cnum),newrnum(k),residue_defs);
    else
        disp('Aber Hallo!');
    end
end

% sort residues by number
info=model.structures{snum}(cnum).residues{mnum}.info;
numbers=length(info);
for k=1:length(info),
    numbers(k)=info(k).number;
end
info0=info;
[sorted,oldnumbers]=sort(numbers);
tags=':';
for k=1:length(oldnumbers),
    info(k)=info0(oldnumbers(k));
    tag=id2tag(oldnumbers(k),model.structures{snum}(cnum).residues{mnum}.residue_tags);
    tags=[tags tag ':'];
end
model.structures{snum}(cnum).residues{mnum}.residue_tags=tags;
model.structures{snum}(cnum).residues{mnum}.info=info;

function [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions(fid_report,indices,restrain,monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1)

adr = mk_address(indices);

fprintf(fid_report,'\n--- Restraint fulfillment for model %i with weight %5.3f ---\n\n',indices(3),p_model);

poi = 0;
dispoi = 0;
mdispoi = 0;
for k = 1:length(restrain)
    resnum = res1+k-1;
    radr = sprintf('%s%i',adr,resnum);
    fprintf(fid_report,'\nResidue %i\n\n',resnum);
    rindices = resolve_address(radr);
    for kr = 1:length(restrain(k).r_beacon)
        sep = strfind(restrain(k).r_beacon(kr).label_type,'|');
        if isempty(sep)
            label1 = restrain(k).r_beacon(kr).label_type;
            label2 = label1;
        else
            label1 = restrain(k).r_beacon(kr).label_type(1:sep-1);
            label2 = restrain(k).r_beacon(kr).label_type(sep+1:end);
        end
        NO_pos1 = get_NO_pos(rindices,label1,restrain(k).r_beacon(kr).label_T);
        % [rmsd2,xyz1]=NOpos_rmsd(NO_pos1);
        NO_pos2 = get_NO_pos(restrain(k).r_beacon(kr).bindices,label2,restrain(k).r_beacon(kr).label_T);
        % [rmsd2,xyz2]=NOpos_rmsd(NO_pos2);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).r_beacon(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
%                 figure(poi); clf;
                rarg = (rax -restrain(k).r_beacon(kr).par1)/restrain(k).r_beacon(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Beacon restraint %i-%i with <r_{exp}> = %5.2f',resnum,restrain(k).r_beacon(kr).resb,restrain(k).r_beacon(kr).par1);
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                    restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
            case 'bounds'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                exp_distr = (rax>=restrain(k).r_beacon(kr).par1) & (rax<= restrain(k).r_beacon(kr).par2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                if rmean >= restrain(k).r_beacon(kr).par1 && rmean <= restrain(k).r_beacon(kr).par2
                    fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                        restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
                else
                    fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] ? to residue %i violated at <r> = %5.2f\n',...
                        restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
                end
        end
    end

    for kr = 1:length(monitor(k).r_beacon)
        sep = strfind(monitor(k).r_beacon(kr).label_type,'|');
        if isempty(sep)
            label1 = monitor(k).r_beacon(kr).label_type;
            label2 = label1;
        else
            label1 = monitor(k).r_beacon(kr).label_type(1:sep-1);
            label2 = monitor(k).r_beacon(kr).label_type(sep+1:end);
        end
        NO_pos1 = get_NO_pos(rindices,label1,monitor(k).r_beacon(kr).label_T);
        NO_pos2 = get_NO_pos(monitor(k).r_beacon(kr).bindices,label2,monitor(k).r_beacon(kr).label_T);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Beacon monitor %i-%i with <r_{sim}> = %5.2f',resnum,monitor(k).r_beacon(kr).resb,rmean);
        fprintf(fid_report,'Monitored distance to residue %i is <r> = %5.2f\n',monitor(k).r_beacon(kr).resb,rmean);
    end

    for kr = 1:length(restrain(k).r_intern)
        r2adr = sprintf('%s%i',adr,restrain(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        sep = strfind(restrain(k).r_intern(kr).label_type,'|');
        if isempty(sep)
            label1 = restrain(k).r_intern(kr).label_type;
            label2 = label1;
        else
            label1 = restrain(k).r_intern(kr).label_type(1:sep-1);
            label2 = restrain(k).r_intern(kr).label_type(sep+1:end);
        end
        NO_pos1 = get_NO_pos(rindices,label1,restrain(k).r_intern(kr).label_T);
        % [rmsd2,xyz1]=NOpos_rmsd(NO_pos1);
        NO_pos2 = get_NO_pos(r2indices,label2,restrain(k).r_intern(kr).label_T);
        % [rmsd2,xyz2]=NOpos_rmsd(NO_pos2);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).r_intern(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
%                 figure(poi); clf;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Internal restraint %i-%i with <r_{exp}> = %5.2f',resnum,restrain(k).r_intern(kr).resb,restrain(k).r_intern(kr).par1);

%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                    restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
            case 'bounds'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                exp_distr = (rax>=restrain(k).r_intern(kr).par1) & (rax<= restrain(k).r_intern(kr).par2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Internal restraint %i-%i with bounds (%5.2f,%5.2f)',resnum,restrain(k).r_intern(kr).resb,restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2);
                if rmean >= restrain(k).r_intern(kr).par1 && rmean <= restrain(k).r_intern(kr).par2,
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                else
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i violated at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                end
        end
    end
    
    for kr = 1:length(monitor(k).r_intern),
        r2adr = sprintf('%s%i',adr,monitor(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        sep = strfind(monitor(k).r_intern(kr).label_type,'|');
        if isempty(sep)
            label1 = monitor(k).r_intern(kr).label_type;
            label2 = label1;
        else
            label1 = monitor(k).r_intern(kr).label_type(1:sep-1);
            label2 = monitor(k).r_intern(kr).label_type(sep+1:end);
        end
        NO_pos1 = get_NO_pos(rindices,label1,monitor(k).r_intern(kr).label_T);
        NO_pos2 = get_NO_pos(r2indices,label2,monitor(k).r_intern(kr).label_T);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Internal restraint %i-%i',resnum,monitor(k).r_intern(kr).resb);
        fprintf(fid_report,'Monitored distance to residue %i is <r> = %5.2f\n',monitor(k).r_intern(kr).resb,rmean);
    end

    for kr = 1:length(restrain(k).depth),
        NO_pos = get_NO_pos(rindices,restrain(k).depth(kr).label_type,restrain(k).depth(kr).label_T);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).depth(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Depth restraint %i with <r_{exp}> = %5.2f',resnum,restrain(k).r_depth(kr).par1);

%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                    restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
            case 'bounds'
                poi = poi + 1;
                dispoi = dispoi + 1;
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                exp_distr = (rax>=restrain(k).depth(kr).par1) & (rax<= restrain(k).depth(kr).par2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraints-Distr{distpoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Depth restraint %i with bounds = %5.2f,%5.2f',resnum,restrain(k).r_depth(kr).par1,restrain(k).r_depth(kr).par2);
                if rmean >= restrain(k).depth(kr).par1 && rmean <= restrain(k).depth(kr).par2
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                else
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? violated at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                end
        end
    end

    for kr = 1:length(monitor(k).depth),
        NO_pos = get_NO_pos(rindices,monitor(k).depth(kr).label_type,monitor(k).depth(kr).label_T);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored depth %i is <r_{sim}> = %5.2f',resnum,rmean);
        fprintf(fid_report,'Monitored depth is <r> = %5.2f\n',rmean);
    end

    for kr = 1:length(restrain(k).oligomer),
        NO_pos = get_NO_pos(rindices,restrain(k).oligomer(kr).label_type,restrain(k).oligomer(kr).label_T);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,restrain(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).oligomer(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).oligomer(kr).par1)/restrain(k).oligomer(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Oligomer restraint %i (n = %i) with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).oligomer(kr).n,restrain(k).oligomer(kr).par1,rmean);
%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                    restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
            case 'bounds'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                exp_distr = (rax>=restrain(k).oligomer(kr).par1) & (rx<= restrain(k).oligomer(kr).par2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Oligomer restraint %i (n = %i) with <r_{sim}> = %5.2f',resnum,restrain(k).oligomer(kr).n,rmean);
                if rmean >= restrain(k).oligomer(kr).par1 && rmean <= restrain(k).oligomer(kr).par2,
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                else
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? violated at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                end
        end
    end

    for kr = 1:length(monitor(k).oligomer),
        NO_pos = get_NO_pos(rindices,monitor(k).oligomer(kr).label_type,monitor(k).oligomer(kr).label_T);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,monitor(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored oligomer distance %i (n = %i) is <r_{sim}> = %5.2f',resnum,monitor(k).oligomer(kr).n,rmean);
        fprintf(fid_report,'Monitored oligomer distance is <r> = %5.2f\n',rmean);
    end

end

function [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions_free(fid_report,indices,restrain,monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1)

rindices = indices;
rindices(3) = 1;
adr = mk_address(rindices);

fprintf(fid_report,'\n--- Restraint fulfillment for model %i with weight %5.3f ---\n\n',indices(3),p_model);

poi = 0;
dispoi = 0;
mdispoi = 0;
for k = 1:length(restrain),
    resnum = res1+k-1;
    radr = sprintf('%s%i',adr,resnum);
    fprintf(fid_report,'\nResidue %i\n\n',resnum);
    rindices = resolve_address(radr);

    for kr = 1:length(restrain(k).r_intern),
        r2adr = sprintf('%s%i',adr,restrain(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        NO_pos1 = get_NO_pos(rindices,restrain(k).r_intern(kr).label_type,restrain(k).r_intern(kr).label_T);
        % [rmsd2,xyz1]=NOpos_rmsd(NO_pos1);
        NO_pos2 = get_NO_pos(r2indices,restrain(k).r_intern(kr).label_type,restrain(k).r_intern(kr).label_T);
        % [rmsd2,xyz2]=NOpos_rmsd(NO_pos2);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).r_intern(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
%                 figure(poi); clf;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Internal restraint %i-%i with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).r_intern(kr).resb,restrain(k).r_intern(kr).par1,rmean);

%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                    restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
            case 'bounds'
                if rmean >= restrain(k).r_intern(kr).par1 && rmean <= restrain(k).r_intern(kr).par2,
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i fulfilled at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                else
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] ? to residue %i violated at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                end
        end
    end
    
    for kr = 1:length(monitor(k).r_intern),
        r2adr = sprintf('%s%i',adr,monitor(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        NO_pos1 = get_NO_pos(rindices,monitor(k).r_intern(kr).label_type,monitor(k).r_intern(kr).label_T);
        NO_pos2 = get_NO_pos(r2indices,monitor(k).r_intern(kr).label_type,monitor(k).r_intern(kr).label_T);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Internal restraint %i-%i with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).r_intern(kr).resb,restrain(k).r_intern(kr).par1,rmean);
        fprintf(fid_report,'Monitored distance to residue %i is <r> = %5.2f\n',monitor(k).r_intern(kr).resb,rmean);
    end

    for kr = 1:length(restrain(k).depth),
        NO_pos = get_NO_pos(rindices,restrain(k).depth(kr).label_type,restrain(k).depth(kr).label_T);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).depth(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Depth restraint %i with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).r_depth(kr).par1,rmean);

%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                    restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
            case 'bounds'
                if rmean >= restrain(k).depth(kr).par1 && rmean <= restrain(k).depth(kr).par2,
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                else
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] ? violated at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                end
        end
    end

    for kr = 1:length(monitor(k).depth),
        NO_pos = get_NO_pos(rindices,monitor(k).depth(kr).label_type,monitor(k).depth(kr).label_T);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored depth %i is <r_{sim}> = %5.2f',resnum,rmean);
        fprintf(fid_report,'Monitored depth is <r> = %5.2f\n',rmean);
    end

    for kr = 1:length(restrain(k).oligomer),
        NO_pos = get_NO_pos(rindices,restrain(k).oligomer(kr).label_type,restrain(k).oligomer(kr).label_T);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,restrain(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).oligomer(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).oligomer(kr).par1)/restrain(k).oligomer(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Oligomer restraint %i (n = %i) with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).oligomer(kr).n,restrain(k).oligomer(kr).par1,rmean);
%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                    restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
            case 'bounds'
                if rmean >= restrain(k).oligomer(kr).par1 && rmean <= restrain(k).oligomer(kr).par2,
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? fulfilled at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                else
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] ? violated at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                end
        end
    end

    for kr = 1:length(monitor(k).oligomer),
        NO_pos = get_NO_pos(rindices,monitor(k).oligomer(kr).label_type,monitor(k).oligomer(kr).label_T);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,monitor(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored oligomer distance %i (n = %i) is <r_{sim}> = %5.2f',resnum,monitor(k).oligomer(kr).n,rmean);
        fprintf(fid_report,'Monitored oligomer distance is <r> = %5.2f\n',rmean);
    end

end


function NO_pos = get_NO_pos(indices,label,T)

global model
global label_defs
global hMain

if strcmpi(label,'CA'),
    adr = sprintf('%s.CA',mk_address(indices));
    [~,xyz] = get_object(adr,'coor');
    NO_pos = [xyz 1];
    return
end

NO_pos = [];
if isfield(model,'sites'),
    for k0=1:length(model.sites),
        for k1=1:length(model.sites{k0}),
            for k=1:length(model.sites{k0}(k1).residue),
                if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                    id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                    if strcmpi(label,label_defs.residues(id).short_name) || strcmpi(label,label_defs.residues(id).tc),
                        if T == model.sites{k0}(k1).residue(k).T,
                            NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                        end
                    end
                end
            end
        end
    end
end

if isempty(NO_pos),
    adr = mk_address(indices);
    command=sprintf('rotamers %s %s %i',adr,label,T);
    hMain.store_undo=false;
    hMain.dynamic_rotamers=false;
    cmd(hMain,command);
end

for k0=1:length(model.sites),
    for k1=1:length(model.sites{k0}),
        for k=1:length(model.sites{k0}(k1).residue),
            if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                if strcmpi(label,label_defs.residues(id).short_name) ||  strcmpi(label,label_defs.residues(id).tc)
                    if T == model.sites{k0}(k1).residue(k).T,
                        NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                    end
                end
            end
        end
    end
end


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy = true;
update_plot(handles);

% --- Executes on button press in pushbutton_forward.
function pushbutton_forward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plot_nr < handles.n_restraints + handles.n_monitor,
    handles.plot_nr = handles.plot_nr + 1;
    update_plot(handles);
else
    add_msg_board('Warning: This was the last plot available.');
end

% --- Executes on button press in pushbutton_backward.
function pushbutton_backward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plot_nr > 0,
    handles.plot_nr = handles.plot_nr - 1;
    update_plot(handles);
else
    add_msg_board('Warning: Use forward (>) button to access distribution plots.');
end

function update_plot(handles)

if handles.copy
    figure;
    handles.copy = false;
else
    axes(handles.axes_success_length);
    cla;
end
sc = 1;
if handles.plot_nr == 0
    plot(handles.resax,handles.success_distr,'k');
    set(gca,'FontSize',8);
    xlabel('Residue');
    ylabel('Rejection probability');
    title('');
    set(handles.text_axes_title,'String','Rejection distribution along the loop');
elseif handles.plot_nr <= handles.n_restraints
    plot(handles.rax,handles.distributions{handles.plot_nr},'k');
    hold on;
    plot(handles.rax,handles.restraint_distr{handles.plot_nr},'r');
    [overlap,shift] = get_overlap(handles.rax,handles.distributions{handles.plot_nr},handles.restraint_distr{handles.plot_nr});
    set(handles.text_overlap_distr,'String',sprintf('%5.3f',overlap));
    set(handles.text_shift_distr,'String',sprintf('%5.2f',shift));
    set(gca,'FontSize',8);
    xlabel('Distance [?]');
    ylabel('Probability density');
    sc = max(handles.restraint_distr{handles.plot_nr});
    title(handles.descriptors{handles.plot_nr});
    set(handles.text_axes_title,'String','Restraint fulfillment by the ensemble');
else
    plot(handles.rax,handles.monitor_distr{handles.plot_nr-handles.n_restraints},'r');
    sc = max(handles.monitor_distr{handles.plot_nr-handles.n_restraints});
    hold on;
    set(gca,'FontSize',8);
    xlabel('Distance [?]');
    ylabel('Probability density');
    title(handles.monitor_descr{handles.plot_nr-handles.n_restraints});
    set(handles.text_axes_title,'String','Distance distribution in the ensemble');
end
if handles.plot_nr > 0 && ~isempty(handles.exp_distr{handles.plot_nr}) && ~isempty(handles.exp_r{handles.plot_nr})
    plot(handles.exp_r{handles.plot_nr},sc*handles.exp_distr{handles.plot_nr}/max(handles.exp_distr{handles.plot_nr}),'b');
end
guidata(handles.text_axes_title,handles);


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plot_nr > 0,
    handles = load_DeerAnalysis_MMM(handles,true);
    handles.exp_distr{handles.plot_nr} = handles.dexp/sum(handles.dexp);
    handles.exp_r{handles.plot_nr} = 10*handles.rexp;
else
    add_msg_board('Warning: Distance distribution cannot be loaded for this plot.');
end
update_plot(handles)
guidata(hObject,handles);

function delete_empty_models(snum,success)

global model

chains=length(model.structures{snum});
for cnum=1:chains,
    for k=1:success,
        atoms{k} = model.structures{snum}(cnum).atoms{k};
    end
    model.structures{snum}(cnum).atoms = atoms;
    for k=1:success,
        residues{k} = model.structures{snum}(cnum).residues{k};
    end
    model.structures{snum}(cnum).residues = residues;    
    for k=1:success,
        xyz{k} = model.structures{snum}(cnum).xyz{k};
    end
    model.structures{snum}(cnum).xyz = xyz;
    for k=1:success,
        Bfactor{k} = model.structures{snum}(cnum).Bfactor{k};
    end
    model.structures{snum}(cnum).Bfactor = Bfactor;
    for k=1:success,
        Btensor{k} = model.structures{snum}(cnum).Btensor{k};
    end
    model.structures{snum}(cnum).Btensor = Btensor;
end

function replace_model(indices,structure)

global model

snum = indices(1);
cnum = indices(2);
modnum = indices(3);

model.structures{snum}(cnum).atoms{modnum} = structure(1).atoms{1};
model.structures{snum}(cnum).residues{modnum} = structure(1).residues{1};
model.structures{snum}(cnum).xyz{modnum} = structure(1).xyz{1};
model.structures{snum}(cnum).Bfactor{modnum} = structure(1).Bfactor{1};
model.structures{snum}(cnum).Btensor{modnum} = structure(1).Btensor{1};
        


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sc = 1;
if handles.plot_nr == 0,
    data = [handles.resax',handles.success_distr'];
    comment = '% Rejection distribution along the loop';
    suggestion = 'rejection_distribution.dat';
    save_data(data,comment,suggestion);
elseif handles.plot_nr <= handles.n_restraints,
    data= [handles.rax',handles.distributions{handles.plot_nr}',handles.restraint_distr{handles.plot_nr}'];
    comment = ['% ' handles.descriptors{handles.plot_nr}];
    suggestion = sprintf('restraint_fullfillment_%i.dat',handles.plot_nr);
    save_data(data,comment,suggestion);
else
    data = [handles.rax',handles.monitor_distr{handles.plot_nr-handles.n_restraints}'];
    comment = ['% ' handles.monitor_descr{handles.plot_nr-handles.n_restraints}];
    suggestion = sprintf('monitored_distribution_%i.dat',handles.plot_nr);
    save_data(data,comment,suggestion);
end
guidata(handles.pushbutton_save,handles);

function save_data(data,comment,suggestion)

global general

my_path=pwd;
cd(general.DEER_files);

[fname,pname]=uiputfile('*.txt','Save plotted data',suggestion);
if isequal(fname,0) || isequal(pname,0),
    add_msg_board('Saving of distance distribution and restraint canceled by user.');
    return;
end
reset_user_paths(pname);
general.DEER_files=pname;

wfid = fopen(fullfile(pname,fname),'wt');
fprintf(wfid,'%s\n',comment);
[m,n] = size(data);
for k = 1:m,
    for k2  = 1:n,
        fprintf(wfid,'%12.8f',data(k,k2));
        if k2 < n,
            fprintf(wfid,'\t');
        else
            fprintf(wfid,'\n');
        end
    end
end
    
fclose(wfid);

cd(my_path);


% --- Executes on button press in checkbox_std_seed.
function checkbox_std_seed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_std_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_std_seed

function [overlap,shift] = get_overlap(rax,distr1,distr2)

distr1 = distr1/sum(distr1);
distr2 = distr2/sum(distr2);
shift = sum(rax.*distr1) - sum(rax.*distr2); 
% fprintf(1,'Sim. mean distance: %5.2f\n',sum(rax.*distr1));
% fprintf(1,'Restraint mean distance: %5.2f\n',sum(rax.*distr2));
overlap = sum(min([distr1;distr2]));
    

function spin_coor = get_bilabel_pos(rindices,label)
% provides cooordinates and populations of all rotamers/conformations of a
% label attached at two sites
% rindices  [2x4] array of two residue addresses
% label     label type

global model

spin_coor = [];

% check the existing bilabels
if isfield(model,'bisites')
    for kscan = 1:length(model.bisites)
        for ksite = 1:length(model.bisites{kscan}.sites)
            sindices = model.bisites{kscan}.sites{ksites}.indices;
            if sum(abs(sindices-rindices)) == 0 || sum(abs(sindices-flipud(rindices))) == 0
                if strcmpi(model.bisites{kscan}.sites{ksites}.label,label)
                    spin_coor = model.bisites{kscan}.sites{ksites}.coor(:,1:4);
                end
            end
        end
    end
end

if isempty(spin_coor) % the requested bilabel was not found
    adr1 = mk_address(rindices(1,:));
    adr2 = mk_address(rindices(2,:));
    [msg,coor1] = get_object(ind1(1:3),'xyz_heavy');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates of the labelled chain could not be obtained (%s)',msg.text));
        return
    end
    [msg,elements] = get_object(ind1(1:3),'elements_heavy');
    if msg.error
        add_msg_board(sprintf('ERROR: Element specifiers of the labelled chain could not be obtained (%s)',msg.text));
        return
    end
    ecoor = [elements,coor1];
    spin_coor = two_pronged_label(adr1,adr2,label,ecoor);
end