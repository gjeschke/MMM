function varargout = compare_ensembles(varargin)
% COMPARE_ENSEMBLES MATLAB code for compare_ensembles.fig
%      COMPARE_ENSEMBLES, by itself, creates a new COMPARE_ENSEMBLES or raises the existing
%      singleton*.
%
%      H = COMPARE_ENSEMBLES returns the handle to a new COMPARE_ENSEMBLES or the handle to
%      the existing singleton*.
%
%      COMPARE_ENSEMBLES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPARE_ENSEMBLES.M with the given input arguments.
%
%      COMPARE_ENSEMBLES('Property','Value',...) creates a new COMPARE_ENSEMBLES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before compare_ensembles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to compare_ensembles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help compare_ensembles

% Last Modified by GUIDE v2.5 04-Feb-2020 10:55:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @compare_ensembles_OpeningFcn, ...
                   'gui_OutputFcn',  @compare_ensembles_OutputFcn, ...
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


% --- Executes just before compare_ensembles is made visible.
function compare_ensembles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to compare_ensembles (see VARARGIN)

global MMM_icon
global hMain

% Choose default command line output for compare_ensembles
handles.output = hObject;

handles.ensemble1.exist = false;
handles.ensemble2.exist = false;

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

hMain.auxiliary=[hMain.auxiliary hObject];

handles.pushbutton_compare.Enable = 'off';
handles.pushbutton_correlation1.Enable = 'off';
handles.pushbutton_correlation2.Enable = 'off';
handles.pushbutton_restraint_plots.Enable = 'off';
handles.ensemble1.pair_rmsd = [];
handles.ensemble2.pair_rmsd = [];
handles.cross_std = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes compare_ensembles wait for user response (see UIRESUME)
% uiwait(handles.figure_ensemble_comparison);


% --- Outputs from this function are returned to the command line.
function varargout = compare_ensembles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes_rmsd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_rmsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_rmsd

axis(hObject);
xlabel('Conformer');
ylabel('Conformer');


% --- Executes on button press in pushbutton_load_ensemble1.
function pushbutton_load_ensemble1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_ensemble1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[fname,pname]=uigetfile('*.dat','Load ensemble description from file');

poi = strfind(fname,'.dat');
basname = fname(1:poi-1);
ensemble_pdb = fullfile(pname,strcat(basname,'.pdb'));

hfig = gcf;

if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Ensemble loading cancelled by user');
    return
elseif ~exist(ensemble_pdb,'file')
    add_msg_board('PDB file of ensemble not found');
    return    
else
    set(hfig,'Pointer','watch');
    [file_list,pop] = rd_ensemble_description(fullfile(pname,fname));
    if isempty(file_list)
        add_msg_board('ERROR: Reading of ensemble file failed.');
        return
    end
    [message,snum] = add_pdb(ensemble_pdb);     
    if message.error
        add_msg_board(sprintf('Ensemble PDB file could not be read (%s)',message.text));
        return
    end
    handles.ensemble1.name = basname;
    handles.ensemble1.file_list = file_list;
    handles.ensemble1.pop = pop;
    handles.ensemble1.snum = snum;
    handles = setup_ensemble(handles,1);
    handles.pushbutton_correlation1.Enable = 'on';
    handles.pushbutton_restraint_plots.Enable = 'on';
    handles.cross_std = [];
    set(hfig,'Pointer','arrow');
end

handles.text_ensemble1.String = basname;
guidata(hObject,handles);

% --- Executes on button press in pushbutton_load_ensemble2.
function pushbutton_load_ensemble2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_ensemble2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname,pname]=uigetfile('*.dat','Load ensemble description from file');

poi = strfind(fname,'.dat');
basname = fname(1:poi-1);
ensemble_pdb = fullfile(pname,strcat(basname,'.pdb'));

hfig = gcf;

if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Ensemble loading cancelled by user');
    return
elseif ~exist(ensemble_pdb,'file')
    add_msg_board('PDB file of ensemble not found');
    return    
else
    set(hfig,'Pointer','watch');
    [file_list,pop] = rd_ensemble_description(fullfile(pname,fname));
    if isempty(file_list)
        add_msg_board('ERROR: Reading of ensemble file failed');
        return
    end
    [message,snum] = add_pdb(ensemble_pdb);     
    if message.error
        add_msg_board(sprintf('Ensemble PDB file could not be read (%s)',message.text));
        return
    end
    handles.ensemble2.name = basname;
    handles.ensemble2.file_list = file_list;
    handles.ensemble2.pop = pop;
    handles.ensemble2.snum = snum;
    handles = setup_ensemble(handles,2);
    handles.cross_std = [];
    handles.pushbutton_correlation2.Enable = 'on';
    set(hfig,'Pointer','arrow'); 
end

handles.text_ensemble2.String = basname;
guidata(hObject,handles);


% --- Executes on selection change in popupmenu_chains1.
function popupmenu_chains1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_chains1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chains1


% --- Executes during object creation, after setting all properties.
function popupmenu_chains1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_chains2.
function popupmenu_chains2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_chains2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chains2


% --- Executes during object creation, after setting all properties.
function popupmenu_chains2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_assignment.
function popupmenu_assignment_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_assignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_assignment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_assignment


% --- Executes during object creation, after setting all properties.
function popupmenu_assignment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_assignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_remove_pair.
function pushbutton_remove_pair_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.ensemble1.exist && handles.ensemble2.exist
    sel1 = handles.popupmenu_chains1.Value;
    sel2 = handles.popupmenu_chains2.Value;
    if handles.pairing(sel1) == sel2
        handles.pairing(sel1) = 0;
        handles = update_pairlist(handles);
    else
        add_msg_board('Warning: This chain pair does not exist in chain assignment.');
    end
else
    add_msg_board('Warning: Chain assignment requires that two valid ensemble have been loaded.');
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_add_pair.
function pushbutton_add_pair_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ensemble1.exist && handles.ensemble2.exist
    sel1 = handles.popupmenu_chains1.Value;
    sel2 = handles.popupmenu_chains2.Value;
    handles.pairing(sel1) = sel2;
    handles = update_pairlist(handles);
else
    add_msg_board('Warning: Chain assignment requires that two valid ensemble have been loaded.');
end
guidata(hObject,handles);

% --- Executes on selection change in popupmenu_display.
function popupmenu_display_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_display

sel = hObject.Value;
switch sel
    case 1
        if ~isempty(handles.ensemble1.pair_rmsd)
            options.xlabel = 'Conformer number';
            options.ylabel = 'Conformer number';
            handles = plot_pair_rmsd(handles,handles.ensemble1.pair_rmsd,options);
        end
    case 2
        if ~isempty(handles.ensemble2.pair_rmsd)
            options.xlabel = 'Conformer number';
            options.ylabel = 'Conformer number';
            handles = plot_pair_rmsd(handles,handles.ensemble2.pair_rmsd,options);
        end
    case 3
        if isfield(handles,'cross_std') && ~isempty(handles.cross_std)
            options.xlabel = 'Conformer number';
            options.ylabel = 'Conformer number';
            handles = plot_pair_rmsd(handles,handles.cross_std,options);
        end
    case 4
        if isfield(handles.ensemble1,'corr_isotropic') && ~isempty(handles.ensemble1.corr_isotropic)
            options.xlabel = 'Residue number';
            options.ylabel = 'Residue number';
            handles = plot_pair_rmsd(handles,handles.ensemble1.corr_isotropic,options);
        end
    case 5
        if isfield(handles.ensemble2,'corr_isotropic') && ~isempty(handles.ensemble2.corr_isotropic)
            options.xlabel = 'Residue number';
            options.ylabel = 'Residue number';
            handles = plot_pair_rmsd(handles,handles.ensemble2.corr_isotropic,options);
        end
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_compare.
function pushbutton_compare_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

stag1 = id2tag(handles.ensemble1.snum,model.structure_tags);
pop1 = handles.ensemble1.pop;
stag2 = id2tag(handles.ensemble2.snum,model.structure_tags);
pop2 = handles.ensemble2.pop;
chains1 = '';
chains2 = '';
for k = 1:length(handles.pairing)
    if handles.pairing(k) ~= 0
        chains1 = [chains1 handles.ensemble1.chains{k}]; 
        chains2 = [chains2 handles.ensemble1.chains{handles.pairing(k)}]; 
    end
end
options.mode = 'backbone';
[sigma,pair_rmsd,ordering] = ensemble_comparison(stag1,chains1,pop1,stag2,chains2,pop2,options);

handles.text_cross_std.String = sprintf('%5.2f Å',sigma(3));
ensemble1_stdv = pair_rmsd(1);
ensemble2_stdv = pair_rmsd(2);
cross_stdv = pair_rmsd{3};
options.xlabel = 'Conformer number';
options.ylabel = 'Conformer number';
handles = plot_pair_rmsd(handles,cross_stdv,options);
handles.cross_std = cross_stdv;
handles.cross_xaxis = ordering{3};
handles.cross_yaxis = ordering{4};
handles.all_sigma = sigma;

oname = sprintf('%s_%s_cross_stdv.mat',handles.ensemble1.name,handles.ensemble2.name);
save(oname,'sigma','cross_stdv','ensemble1_stdv','ensemble2_stdv');

handles.popupmenu_display.Value = 3;

guidata(hObject,handles);


function [file_list,pop] = rd_ensemble_description(fname)

file_list = cell(1,1000);
pop = zeros(1,1000);
poi = 0;

fid=fopen(fname);
if fid == -1
    add_msg_board('Warning: File list does not exist');
    file_list = {};
    return;
end

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        if ~contains(tline,'%')
            myline = textscan(tline,'%s');
            args = myline{1};
            arg1 = char(args(1));
            if ~isempty(arg1)
                poi = poi + 1;
                file_list{poi} = arg1;
                pop(poi) = str2double(char(args(2)));
            end
        end
    end
end

fclose(fid);

file_list = file_list(1:poi);
pop = pop(1:poi);

function handles = setup_ensemble(handles,num)

global model

switch num
    case 1
        basname = handles.ensemble1.name;
        snum = handles.ensemble1.snum;
        pop = handles.ensemble1.pop;
        text_rmsd = handles.text_rmsd1;
        text_nconf = handles.text_nconformers1;
        popup_chains = handles.popupmenu_chains1;
    case 2
        basname = handles.ensemble2.name;
        snum = handles.ensemble2.snum;
        pop = handles.ensemble2.pop;
        text_rmsd = handles.text_rmsd2;
        text_nconf = handles.text_nconformers2;
        popup_chains = handles.popupmenu_chains2;
    otherwise
        return
end

nch = length(model.structures{snum});
nconf = length(model.structures{snum}(1).xyz);
text_nconf.String = sprintf('%i',nconf);
resnums = zeros(1,nch);
seqtypes = zeros(1,nch);
for k = 1:nch
    chain_list{k} = model.structures{snum}(k).name;
    resnums(k) = length(model.structures{snum}(k).residues{1}.info);
    seqtypes(k) = model.structures{snum}(k).seqtype;
end
popup_chains.String = chain_list;
popup_chains.String = popup_chains.String(1:nch);

stag = id2tag(snum,model.structure_tags);

options.mode = 'backbone';
options.cluster_sorting = true;

[sigma,pair_rmsd,ordering] = ensemble_comparison(stag,[],pop,[],[],[],options);
text_rmsd.String = sprintf('%5.2f Å',sigma);
populations = pop(ordering{1});

pair_rmsd = pair_rmsd{1};
options.xlabel = 'Conformer number';
options.ylabel = 'Conformer number';
handles = plot_pair_rmsd(handles,pair_rmsd,options);
oname = strcat(basname,'_pair_rmsd.mat');
save(oname,'sigma','pair_rmsd','resnums','seqtypes','ordering','populations');

switch num
    case 1
        handles.ensemble1.exist = true;
        handles.ensemble1.sigma = sigma;
        handles.ensemble1.pair_rmsd = pair_rmsd;
        handles.ensemble1.pr_axis = ordering{1};
        handles.ensemble1.resnums = resnums;
        handles.ensemble1.seqtypes = seqtypes;
        handles.ensemble1.chains = chain_list;
        handles.ensemble1.ordering = ordering{1};
    case 2
        handles.ensemble2.exist = true;
        handles.ensemble2.sigma = sigma;
        handles.ensemble2.pair_rmsd = pair_rmsd;
        handles.ensemble2.pr_axis = ordering{1};
        handles.ensemble2.resnums = resnums;
        handles.ensemble2.seqtypes = seqtypes;
        handles.ensemble2.chains = chain_list;
        handles.ensemble2.ordering = ordering{2};
end

handles.popupmenu_display.Value = num;

if handles.ensemble1.exist && handles.ensemble2.exist
    handles = assign_pairs(handles);
    if sum(handles.pairing) > 0
        handles.pushbutton_compare.Enable = 'on';
    end
end

function handles = plot_pair_rmsd(handles,pair_rmsd,options)

if ~exist('options','var') || ~isfield(options,'figure') || isempty(options.figure)
    axes(handles.axes_rmsd);
    cla
else
    figure(options.figure);
end
[n1,n2] = size(pair_rmsd);
plot(1,1,'k.');
plot(n1,n2,'k.');
h = image(pair_rmsd,'CDataMapping','scaled','ButtonDownFcn',@(hObject,eventdata,handles)axes_rmsd_ButtonDownFcn);
curr_axis = gca;
curr_axis.YDir = 'normal';
colorbar;
axis tight
xlabel('Conformer number');
ylabel('Conformer number');
axis equal
if exist('options','var') && isfield(options,'title') && ~isempty(options.title)
    title(options.title);
end
if exist('options','var') && isfield(options,'xlabel') && ~isempty(options.xlabel)
    xlabel(options.xlabel);
end
if exist('options','var') && isfield(options,'xlabel') && ~isempty(options.xlabel)
    ylabel(options.ylabel);
end

% --- Executes on button press in pushbutton_autoassign_pairs.
function pushbutton_autoassign_pairs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autoassign_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = assign_pairs(handles);
guidata(hObject,handles);

function handles = assign_pairs(handles)

n1 = length(handles.ensemble1.chains);
n2 = length(handles.ensemble2.chains);
chains1 = handles.ensemble1.chains;
chains2 = handles.ensemble2.chains;

pairing = zeros(1,n1);
assigned = zeros(1,n2);
pairlist{1} = '<no chain pairs>';
ppoi = 0;
for k1 = 1:n1
    for k2 = 1:n2
        if ~assigned(k2)
            if handles.ensemble1.resnums(k1) == handles.ensemble2.resnums(k2)...
                    && handles.ensemble1.seqtypes(k1) == handles.ensemble2.seqtypes(k2)
                ppoi = ppoi+ 1;
                pairlist{ppoi} = sprintf('%s %s',chains1{k1},chains2{k2});
                pairing(k1) = k2;
            end
        end
    end
end

handles.pairing = pairing;
handles.popupmenu_assignment.String = pairlist;

function handles = update_pairlist(handles)

ppoi = 0;
pairlist{1} = '<no pairs selected>';
for k1 = 1:length(handles.pairing)
    if handles.pairing(k1) ~= 0
        ppoi = ppoi+ 1;
        pairlist{ppoi} = sprintf('%s %s',handles.ensemble1.chains{k1},handles.ensemble2.chains{handles.pairing(k1)});
    end
end

if ppoi == 0
    ppoi = 1;
end
pairlist = pairlist(1:ppoi);

handles.popupmenu_assignment.String = pairlist;


% --- Executes on button press in pushbutton_correlation1.
function pushbutton_correlation1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_correlation1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ensemble1.exist
    hfig = gcf;
    set(hfig,'Pointer','watch');
    basname = handles.ensemble1.name;
    snum = handles.ensemble1.snum;
    pop = handles.ensemble1.pop;
    [handles,corr_isotropic,corr_axis,dmat] = get_correlations(handles,basname,snum,pop);
    handles.ensemble1.corr_isotropic = corr_isotropic;
    handles.ensemble1.dmat = dmat;
    handles.ensemble1.corr_axis = corr_axis;
    options.xlabel = 'Residue number';
    options.ylabel = 'Residue number';
    handles = plot_pair_rmsd(handles,handles.ensemble1.corr_isotropic,options);
    handles.popupmenu_display.Value = 4;
    set(hfig,'Pointer','arrow');
 else
    add_msg_board('Warning: No valid ensemble 1 loaded. Aborting.'); 
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_correlation2.
function pushbutton_correlation2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_correlation2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ensemble2.exist
    hfig = gcf;
    set(hfig,'Pointer','watch');
    basname = handles.ensemble2.name;
    snum = handles.ensemble2.snum;
    pop = handles.ensemble2.pop;
    [handles,corr_isotropic,corr_axis,dmat] = get_correlations(handles,basname,snum,pop);
    handles.ensemble2.corr_isotropic = corr_isotropic;
    handles.ensemble2.dmat = dmat;
    handles.ensemble2.corr_axis = corr_axis;
    options.xlabel = 'Residue number';
    options.ylabel = 'Residue number';
    handles = plot_pair_rmsd(handles,handles.ensemble2.corr_isotropic,options);
    handles.popupmenu_display.Value = 5;
    set(hfig,'Pointer','arrow');
else
    add_msg_board('Warning: No valid ensemble 2 loaded. Aborting.'); 
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_ensemble_comparison_CloseRequestFcn(handles.figure_ensemble_comparison, eventdata, handles);


% --- Executes when user attempts to close figure_ensemble_comparison.
function figure_ensemble_comparison_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_ensemble_comparison (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

function [handles,corr_isotropic,corr_axis,dmat] = get_correlations(handles,basname,snum,pop)

global model

nch = length(model.structures{snum});
nconf = length(model.structures{snum}(1).xyz);
resnums = zeros(1,nch);
for k = 1:nch
    resnums(k) = length(model.structures{snum}(k).residues{1}.info);
end
resnum = sum(resnums);
corr_isotropic = zeros(resnum);
dmat = zeros(resnum);
corr_axis = cell(1,resnum);
poi1 = 0;
for kc = 1:nch
    for kr = 1:resnums(kc)
        poi1 = poi1 + 1;
        corr_axis{poi1} = sprintf('%s.%i',model.structures{snum}(kc).name,model.structures{snum}(kc).residues{1}.info(kr).number);
        poi2 = 0;
        for kc2 = 1:nch
            for kr2 = 1:resnums(kc2)
                poi2 = poi2 + 1;
                rvec = zeros(1,nconf);
                for km = 1:nconf
                    xyz1 = model.structures{snum}(kc).xyz{km};
                    xyz2 = model.structures{snum}(kc2).xyz{km};
                    tags = model.structures{snum}(kc).residues{km}.info(kr).atom_tags;
                    id = tag2id('CA',tags);
                    if isempty(id)
                        id = tag2id('C4''',tags);
                    end
                    if ~isempty(id)
                        atnum = model.structures{snum}(kc).residues{km}.info(kr).atom_numbers{id};
                        if length(atnum) == 1
                            corr1 = xyz1(atnum,:);
                        else
                            [m,~] = size(atnum);
                            corr1 = zeros(1,3);
                            for kl = 1:m
                                corr1 = corr1 + atnum(kl,2)*xyz1(atnum(kl,1),:);
                            end
                        end
                        tags = model.structures{snum}(kc2).residues{km}.info(kr2).atom_tags;
                        id2 = tag2id('CA',tags);
                        if isempty(id2)
                            id2 = tag2id('C4''',tags);
                        end
                        if ~isempty(id2)
                            atnum2 = model.structures{snum}(kc2).residues{km}.info(kr2).atom_numbers{id2};
                            if length(atnum2) == 1
                                corr2 = xyz2(atnum2,:);
                            else
                                [m,~] = size(atnum2);
                                corr2 = zeros(1,3);
                                for k2 = 1:m
                                    corr2 = corr2 + atnum2(kl,2)*xyz2(atnum2(kl,1),:);
                                end
                            end
                            rvec(km) = norm(corr1-corr2);
                        end
                    end
                end
                rmean = mean(rvec);
                rvec = rvec - rmean;
                corr_isotropic(poi1,poi2) = sqrt(sum(pop.*rvec.^2));
                corr_isotropic(poi2,poi1) = corr_isotropic(poi1,poi2);
                dmat(poi1,poi2) = rmean;
                dmat(poi2,poi1) = rmean;
            end
        end
    end
end

oname = sprintf('%s_correlations.mat',basname);
save(oname,'corr_isotropic');


% --- Executes on mouse press over axes background.
function axes_rmsd_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_rmsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(gcf);
cp = get(handles.axes_rmsd,'CurrentPoint');
x = round(cp(1,1));
y = round(cp(1,2));
sel = handles.popupmenu_display.Value;
switch sel
    case 1
        my_data = handles.ensemble1.pair_rmsd;
        xaxis = handles.ensemble1.pr_axis;
        yaxis = xaxis;
        mode = 'integer';
    case 2
        my_data = handles.ensemble2.pair_rmsd;
        xaxis = handles.ensemble2.pr_axis;
        yaxis = xaxis;
        mode = 'integer';
    case 3
        my_data = handles.cross_std;
        xaxis = handles.cross_xaxis;
        yaxis = handles.cross_yaxis;
        mode = 'integer';
    case 4
        my_data = handles.ensemble1.corr_isotropic;
        xaxis = handles.ensemble1.corr_axis;
        yaxis = xaxis;
        mode = 'string';
    case 5
        my_data = handles.ensemble2.corr_isotropic;
        xaxis = handles.ensemble2.corr_axis;
        yaxis = xaxis;
        mode = 'string';
end
assignment = '???';
if strcmpi(mode,'integer')
    assignment = sprintf('Conformers (%i,%i)',xaxis(x),yaxis(y));
elseif strcmpi(mode,'string')
    assignment = sprintf('Residues (%s,%s)',xaxis{x},yaxis{y});
end
handles.text_selection.String = sprintf('%s: %5.2f Å',assignment,my_data(y,x));


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.figure = 1;

handles = update_display(handles,options);

guidata(hObject,handles);

function handles = update_display(handles,options)

normalize = 0;
if handles.radiobutton_normalization_r.Value
    normalize = 1;
elseif handles.radiobutton_normalization_sqrt_r.Value
    normalize = 2;
end
    

sel = handles.popupmenu_display.Value;
switch sel
    case 1
        my_data = handles.ensemble1.pair_rmsd;
        options.xlabel = 'Conformer number';
        options.ylabel = 'Conformer number';
        options.title = 'Pair rmsd for ensemble 1 [Å]';
    case 2
        my_data = handles.ensemble2.pair_rmsd;
        options.xlabel = 'Conformer number';
        options.ylabel = 'Conformer number';
        options.title = 'Pair rmsd for ensemble 2 [Å]';
    case 3
        my_data = handles.cross_std;
        options.xlabel = 'Conformer number';
        options.ylabel = 'Conformer number';
        options.title = 'Cross standard deviation [Å]';
    case 4
        my_data = handles.ensemble1.corr_isotropic;
        dist_matrix = handles.ensemble1.dmat;
        [m,~] = size(dist_matrix);
        dist_matrix = dist_matrix + eye(m);
        switch normalize
            case 1
                my_data = my_data./dist_matrix;
            case 2
                my_data = my_data./sqrt(dist_matrix);
        end
        options.xlabel = 'Residue number';
        options.ylabel = 'Residue number';
        options.title = 'Residue-residue correlation (ensemble 1)';
    case 5
        my_data = handles.ensemble2.corr_isotropic;
        dist_matrix = handles.ensemble2.dmat;
        [m,~] = size(dist_matrix);
        dist_matrix = dist_matrix + eye(m);
        switch normalize
            case 1
                my_data = my_data./dist_matrix;
            case 2
                my_data = my_data./sqrt(dist_matrix);
        end
        options.xlabel = 'Residue number';
        options.ylabel = 'Residue number';
        options.title = 'Residue-residue correlation (ensemble 2)';
end

handles = plot_pair_rmsd(handles,my_data,options);


% --- Executes on button press in radiobutton_normalization_none.
function radiobutton_normalization_none_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_normalization_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_normalization_none

options.figure = [];
handles = update_display(handles,options);
guidata(hObject,handles);

% --- Executes on button press in radiobutton_normalization_r.
function radiobutton_normalization_r_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_normalization_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_normalization_r

options.figure = [];
handles = update_display(handles,options);
guidata(hObject,handles);

% --- Executes on button press in radiobutton_normalization_sqrt_r.
function radiobutton_normalization_sqrt_r_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_normalization_sqrt_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_normalization_sqrt_r

options.figure = [];
handles = update_display(handles,options);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_save_ensemble1.
function pushbutton_save_ensemble1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_ensemble1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general

model.current_structure = handles.ensemble1.snum;
idCode = model.info{handles.ensemble1.snum}.idCode;
model_sequence = handles.ensemble1.ordering;

default_name = strcat(handles.ensemble1.name,'_ordered.pdb');

[filename, pathname] = uiputfile(default_name, 'Save ordered ensemble 1 as PDB');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Save as PDB cancelled by user');
else
    reset_user_paths(pathname);
    general.pdb_files=pathname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Ordered ensemble 1 saved as PDB file: %s',fname);
    add_msg_board(msg);
    set(gcf,'Pointer','watch');
    wr_pdb(fname,idCode,model_sequence);
end

set(gcf,'Pointer','arrow');
guidata(hObject,handles);


% --- Executes on button press in pushbutton_save_ensemble2.
function pushbutton_save_ensemble2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_ensemble2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general

model.current_structure = handles.ensemble2.snum;
idCode = model.info{handles.ensemble2.snum}.idCode;
model_sequence = handles.ensemble2.ordering;

default_name = strcat(handles.ensemble1.name,'_ordered.pdb');

[filename, pathname] = uiputfile(default_name, 'Save ordered ensemble 2 as PDB');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Save as PDB cancelled by user');
else
    reset_user_paths(pathname);
    general.pdb_files=pathname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Ordered ensemble 2 saved as PDB file: %s',fname);
    add_msg_board(msg);
    set(gcf,'Pointer','watch');
    wr_pdb(fname,idCode,model_sequence);
end

set(gcf,'Pointer','arrow');
guidata(hObject,handles);

% --- Executes on button press in pushbutton_restraint_plots.
function pushbutton_restraint_plots_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraint_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

% Generate option structure for fitting/restraint evaluation
options.threshold = 0.01; % convert percentage to fraction
options.ensemble_size = length(handles.ensemble1.file_list);
options.depth = [0.005,0.65];
options.dim = 3;
options.kdec = [0, NaN];
options.cutoff = 0.9;
options.lm = 50;
options.fb = 18;
options.xlmax = 20;

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
    plot_groups = get_plot_groups(fullfile(pname,fname));
    unprocessed = true;
    
    switch rtype
        case 'GENERIC'
            restraints = rd_restraints(fullfile(pname,fname));
            if isfield(restraints,'DEER') && ~isempty(restraints.DEER)
                for k = 1:length(restraints.DEER)
                    restraints.DEER(k).r = 10*restraints.DEER(k).r;
                    restraints.DEER(k).sigr = 10*restraints.DEER(k).sigr;
                end
            end
        case 'FLEX'
            restraints = rd_restraints(fullfile(pname,fname));
            if isfield(restraints,'DEER') && ~isempty(restraints.DEER)
                for k = 1:length(restraints.DEER)
                    restraints.DEER(k).r = 10*restraints.DEER(k).r;
                    restraints.DEER(k).sigr = 10*restraints.DEER(k).sigr;
                end
            end
        case 'RIGIFLEX'
            restraints = rd_restraints_rigiflex(fullfile(pname,fname),unprocessed);
        otherwise
            add_msg_board(sprintf('ERROR: Restraint type %s is not supported in ensemble fits',rtype));
            restraints = [];
    end
    cd(my_path);
    options.script = false;
    options.write_pdb = false;
    options.ordering = handles.ensemble1.ordering;
    options.individual = false;
    options.groups = plot_groups;
    [score_DEER,sum_chi2] = mk_ensemble(handles.ensemble1.file_list,handles.ensemble1.pop,restraints,options);
    nc = length(handles.ensemble1.pop);
    fprintf(1,'--- Fitted ensemble with %i conformers---\n',nc);
    fprintf(1,'Total SAS chi^2 value: %6.3f\n',sum_chi2);
    fprintf(1,'Total DEER score : %5.3f\n',score_DEER);

    if handles.checkbox_null_test.Value
        [fname_list,pname_list]=uigetfile('*.dat','Load conformer file list');
        if isequal(fname,0) || isequal(pname,0)
            add_msg_board('Conformer list loading cancelled by user');
            cd(my_path);
            set(hfig,'Pointer','arrow');
            guidata(hObject,handles);
            return
        end
        random_file_list = handles.ensemble1.file_list;
        all_conformers = get_file_list(fullfile(pname_list,fname_list));
        nac = length(all_conformers);
        add_msg_board(sprintf('Evaluating fit quality for %i random ensembles',restraints.random_ensembles));
        fprintf(1,'--- Random ensembles ---\n');
        fprintf(1,'#   SAS chi^2  DEER score\n');
        options.plot = false;
        uniform_populations = ones(1,nc)/nc;
        for krnd = 1:restraints.random_ensembles
            sequencer = rand(1,nac);
            [~,sequence] = sort(sequencer);
            for k = 1:nc
                random_file_list{k} = all_conformers{sequence(k)};
            end
            [score_DEER,sum_chi2] = mk_ensemble(random_file_list,uniform_populations,restraints,options);
            fprintf(1,'%3i%10.3f%12.3f\n',krnd,sum_chi2,score_DEER);
        end
    end
    set(hfig,'Pointer','arrow');
end


guidata(hObject,handles);


% --- Executes on button press in checkbox_null_test.
function checkbox_null_test_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_null_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_null_test
