function varargout = ensemble_flexibility(varargin)
% ENSEMBLE_FLEXIBILITY MATLAB code for ensemble_flexibility.fig
%      ENSEMBLE_FLEXIBILITY, by itself, creates a new ENSEMBLE_FLEXIBILITY or raises the existing
%      singleton*.
%
%      H = ENSEMBLE_FLEXIBILITY returns the handle to a new ENSEMBLE_FLEXIBILITY or the handle to
%      the existing singleton*.
%
%      ENSEMBLE_FLEXIBILITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENSEMBLE_FLEXIBILITY.M with the given input arguments.
%
%      ENSEMBLE_FLEXIBILITY('Property','Value',...) creates a new ENSEMBLE_FLEXIBILITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ensemble_flexibility_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ensemble_flexibility_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ensemble_flexibility

% Last Modified by GUIDE v2.5 12-Jan-2018 12:45:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ensemble_flexibility_OpeningFcn, ...
                   'gui_OutputFcn',  @ensemble_flexibility_OutputFcn, ...
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


% --- Executes just before ensemble_flexibility is made visible.
function ensemble_flexibility_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ensemble_flexibility (see VARARGIN)

% Choose default command line output for ensemble_flexibility
handles.output = hObject;

global model
% global MMM_icon

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

CData = get_detach_icon;
set(handles.pushbutton_detach_plot,'CData',CData/255);
handles.detach_icon = CData/255;
CData = get_attach_icon;
handles.attach_icon = CData/255;
handles.plot_attached = true;
handles.margin = 30; % margin for detached plots in pixels
handles.add_margin = 30; % additional margin for axes labels
handles.plot_colors = [0.5,0,0;0,0.5,0;0,0,0.5;0.8,0.4,0.11;0.8,0,0.8;0.51,0.56,0.87;1/3,1/3,1/3;0.6,0.8,0.2;0.94,0.5,0.5;0.27,0.51,0.71];

handles.rigid_core = 2; % factor with respect to smallest local rmsd that defines rigid-core residues
handles.rigid_core_min = 0.5; % minimum rmsd threshold for rigid-core residues
handles.resnum_offset = 0;

cnum=tag2id(model.current_chain,model.chain_tags{model.current_structure});
idCode=model.info{model.current_structure}.idCode;
type=-1;
sequence=model.structures{model.current_structure}(cnum).sequence;
if isfield(model.structures{model.current_structure}(cnum),'seqtype'), % determine sequence type, information is present
    type = model.structures{model.current_structure}(cnum).seqtype;
    seqtest = double(sequence)>=double('a');
    if type == 2 && sum(seqtest),
        type=3;
    end;
end;
switch type
    case 1
        ctype = 'peptide';
    case 2
        ctype = 'DNA';
    case 3
        ctype = 'RNA';
    otherwise
        add_msg_board('ERROR: Current chain must be a peptide or nucleic acid chain.');
        figure1_CloseRequestFcn(handles.figure1, eventdata, handles);
        return
end
handles.figure1.Name = sprintf('Flexibility analysis for %s chains',ctype);
handles.type = type;
        
chain1 = sprintf('[%s](%s)',idCode,model.current_chain);
handles.chain_address = chain1;
chain_addresses{1} = chain1;
handles.edit_chain_address.String = chain1;
handles.popupmenu_chains.String = chain_addresses;
numrep = length(model.structures{model.current_structure}(cnum).residues);
% if numrep < 2
%     add_msg_board(sprintf('ERROR: Current chain %s is not an ensemble of conformations',chain1));
%     figure1_CloseRequestFcn(handles.figure1, eventdata, handles); 
%     return
% end;
handles.text_curr_chain_models.String = sprintf('%i',numrep);

my_chain = get_chain(chain1,handles.type);
if isempty(my_chain.coor)
    add_msg_board(sprintf('ERROR: No backbone atoms found for %s chain %s',ctype,chain1));
    figure1_CloseRequestFcn(handles.figure1, eventdata, handles); 
    return
end;
handles.chains{1} = my_chain;
info = analyze_rmsd(handles.chains,1,handles.type);
handles.text_central_model.String = sprintf('%i',info.central);
handles.text_rmsd_central.String = sprintf('%5.2f',info.central_rmsd);
handles.info{1} = info;
plot_flexibility(handles,1,true);
update_list(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ensemble_flexibility wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ensemble_flexibility_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if exist('handles','var') && isfield(handles,'output'),
    varargout{1} = handles.output;
else
    varargout{1} = 'cancelled';
end;


% --- Executes on selection change in popupmenu_chains.
function popupmenu_chains_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_chains contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chains

info = handles.info{hObject.Value};
my_chain = handles.chains{hObject.Value};
handles.text_central_model.String = sprintf('%i',info.central);
handles.text_rmsd_central.String = sprintf('%5.2f',info.central_rmsd);
handles.text_curr_chain_models.String = sprintf('%i',my_chain.numrep);
if my_chain.numrep == 1
    set(handles.figure1,'Pointer','watch');
    handles = analyze_rmsd_to_xray(handles,hObject.Value);
    set(handles.figure1,'Pointer','arrow');
    handles.checkbox_single_structure_reference.Enable = 'on';
    handles.checkbox_single_structure_reference.Value = 1;
else
    handles.checkbox_single_structure_reference.Enable = 'off';
    handles.checkbox_single_structure_reference.Value = 0;
end
plot_flexibility(handles,1,1);
for k = 2:length(handles.info),
    plot_flexibility(handles,k,0);
end;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_chains_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add_chain.
function pushbutton_add_chain_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_chain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Pointer','watch');
drawnow;
numchains = length(handles.popupmenu_chains.String);
if numchains >= 10,
    add_msg_board('ERROR: No more than 10 structure (ensembles) can be compared.');
    return
end;
failed = false;
for k = 1:numchains
    if strcmp(handles.popupmenu_chains.String{k},handles.chain_address)
        add_msg_board(sprintf('ERROR: The same chain %s cannot be added a second time.',handles.popupmenu_chains.String{k}));
        failed = true;
    end;
end
if ~failed
    handles.popupmenu_chains.String{numchains+1} = handles.chain_address;
    my_chain = get_chain(handles.chain_address,handles.type);
    handles.chains{numchains+1} = my_chain;
    handles.info{numchains+1} = analyze_rmsd(handles.chains,numchains+1,handles.type);
    
    if length(handles.chains) > 1
        seqs = cell(1,length(handles.chains));
        sindices = zeros(length(handles.chains),2);
        for k = 1: length(handles.chains)
            indices = resolve_address(handles.popupmenu_chains.String{k});
            sindices(k,:) = indices;
            my_chain = handles.chains{k};
            seqs{k} = my_chain.seq;
            %         if handles.type > 1,
            %            seqs{k} = cheat_muscle(seqs{k});
            %         end;
        end
        [message,outfile] = align_sequences(seqs,[],~handles.checkbox_show_alignment.Value);
        if message.error
            add_msg_board(sprintf('%s',message.text));
        end
        alignment = get_multiple_clustal(outfile);
        seqorder = zeros(1,length(alignment));
        for k = 1:length(alignment)
            seqorder(k) = str2double(alignment(k).name(8:end));
        end;
        for k = 1:length(alignment)
            refseq = alignment(k).sequence;
            refstart = 0;
            start = false;
            while ~start && refstart < length(refseq)
                refstart = refstart + 1;
                if refseq(refstart)~='X' && refseq(refstart)~='-'
                    start = true;
                end;
            end;
            my_chain = handles.chains{seqorder(k)};
            my_chain.align_point = refstart;
            handles.chains{seqorder(k)} = my_chain;
        end;
    end;
    update_list(handles);
    set(handles.figure1,'Pointer','arrow');
    handles.popupmenu_chains.Value = numchains+1;
    info = handles.info{numchains+1};
    my_chain = handles.chains{numchains+1};
    handles.text_central_model.String = sprintf('%i',info.central);
    handles.text_rmsd_central.String = sprintf('%5.2f',info.central_rmsd);
    handles.text_curr_chain_models.String = sprintf('%i',my_chain.numrep);
    if my_chain.numrep == 1
        set(handles.figure1,'Pointer','watch');
        handles = analyze_rmsd_to_xray(handles,numchains+1);
        set(handles.figure1,'Pointer','arrow');
        handles.checkbox_single_structure_reference.Enable = 'on';
        handles.checkbox_single_structure_reference.Value = 1;
    else
        handles.checkbox_single_structure_reference.Enable = 'off';
        handles.checkbox_single_structure_reference.Value = 0;
    end
    plot_flexibility(handles,1,1);
    for k = 2:length(handles.info),
        plot_flexibility(handles,k,0);
    end;
end;
guidata(hObject, handles);


function edit_chain_address_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chain_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chain_address as text
%        str2double(get(hObject,'String')) returns contents of edit_chain_address as a double

[indices,message]=resolve_address(hObject.String);
if message.error || length(indices) ~= 2
    add_msg_board(sprintf('Warning: %s is not a chain address. Reverting to current selection',hObject.String));
    sel = handles.popupmenu_chains.Value;
    adr = handles.popupmenu_chains.String{sel};
    hObject.String = adr;
    handles.chain_address = adr;
else
    [stag,ctag] = mk_address_parts(indices);
    handles.chain_address = sprintf('[%s](%s)',stag,ctag);
end;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_chain_address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chain_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

to_be_removed = handles.popupmenu_chains.Value;
chains = length(handles.popupmenu_chains.String);
if chains == 1,
    add_msg_board('ERROR: At least one chain must remain.');
    guidata(hObject, handles);
    return
end
switch to_be_removed
    case 1
        keep = 2:chains;
    case chains
        keep = 1:chains-1;
    otherwise
        keep = [1:to_be_removed-1 to_be_removed+1:chains];
end
handles.chains = handles.chains(keep);
handles.info = handles.info(keep);
handles.popupmenu_chains.Value = 1;
handles.popupmenu_chains.String = handles.popupmenu_chains.String(keep);
info = handles.info{1};
my_chain = handles.chains{1};
handles.text_central_model.String = sprintf('%i',info.central);
handles.text_rmsd_central.String = sprintf('%5.2f',info.central_rmsd);
handles.text_curr_chain_models.String = sprintf('%i',my_chain.numrep);
if length(handles.chains) > 1
    seqs = cell(1,length(handles.chains));
    sindices = zeros(length(handles.chains),2);
    for k = 1: length(handles.chains)
        indices = resolve_address(handles.popupmenu_chains.String{k});
        sindices(k,:) = indices;
        my_chain = handles.chains{k};
        seqs{k} = my_chain.seq;
    end
    [message,outfile] = align_sequences(seqs,[],~handles.checkbox_show_alignment.Value);
    if message.error
        add_msg_board(sprintf('%s',message.text));
    end
    alignment = get_multiple_clustal(outfile);
    seqorder = zeros(1,length(alignment));
    for k = 1:length(alignment)
        seqorder(k) = str2double(alignment(k).name(8:end));
    end;
    for k = 1:length(alignment)
        refseq = alignment(k).sequence;
        refstart = 0;
        start = false;
        while ~start && refstart < length(refseq)
            refstart = refstart + 1;
            if refseq(refstart)~='X' && refseq(refstart)~='-'
                start = true;
            end;
        end;
        my_chain = handles.chains{seqorder(k)};
        my_chain.align_point = refstart;
        handles.chains{seqorder(k)} = my_chain;
    end;
end;
plot_flexibility(handles,1,1);
for k = 2:length(handles.info),
    plot_flexibility(handles,k,0);
end;
update_list(handles);
set(handles.figure1,'Pointer','arrow');
guidata(hObject, handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

basname = handles.popupmenu_chains.String{1};
suggestion=[basname '_flexibility.mat'];
[fname,pname]=uiputfile('*.mat','Save results',suggestion);
if isequal(fname,0) || isequal(pname,0),
    add_msg_board('Saving of distance distribution and DEER simulation canceled by user.');
    return;
end;
save_flexibility(fullfile(pname,fname),handles);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function my_chain = get_chain(address,type)

global model
global residue_defs

switch type
    case 1
        bbapr = 3; % backbone atoms per residue
        bbat{1} = 'N';
        bbat{2} = 'CA';
        bbat{3} = 'C';
    otherwise
        bbapr = 6;
        bbat{1} = 'P';
        bbat{2} = 'O5''';
        bbat{3} = 'C5''';
        bbat{4} = 'C4''';
        bbat{5} = 'C3''';
        bbat{6} = 'O3''';
end;
my_chain.seq = '';
my_chain.coor = {};
my_chain.numrep = 0;
my_chain.resnr = 0;
my_chain.address = address;

[indices,message]=resolve_address(address);
if message.error,
    add_msg_board(message.text);
    return
end;

my_chain.numrep = length(model.structures{indices(1)}(indices(2)).residues);
my_chain.resnr = length(model.structures{indices(1)}(indices(2)).residues{1}.info);

my_chain.coor = cell(1,my_chain.numrep);
my_chain.mask = ones(bbapr*my_chain.resnr,1);
my_chain.resmask = ones(my_chain.resnr,1);
my_chain.resnums = zeros(1,my_chain.resnr);
my_chain.seq = char(1,my_chain.resnr);

for m = 1:my_chain.numrep,
    if length(model.structures{indices(1)}(indices(2)).residues{m}.info) ~= my_chain.resnr
        add_msg_board(sprintf('ERROR: Inconsistent numbers of residues in different chain models'));
        my_chain.seq = '';
        my_chain.coor = {};
        return
    end;
    coor = zeros(3*my_chain.resnr,3);
    maxres = 0;
    for r = 1:my_chain.resnr
        baspoi = bbapr*(r-1);
        resname = model.structures{indices(1)}(indices(2)).residues{m}.info(r).name;
        if ~strcmpi(resname,'HOH'),
            maxres = r;
            resnum = model.structures{indices(1)}(indices(2)).residues{m}.info(r).number;
            my_chain.resnums(r) = resnum;
            if type < 2,
                slc = tag2id(resname,upper(residue_defs.restags),residue_defs.single_letter_code);
            else
                slc = tag2id(resname,upper(residue_defs.nucleotide_tags),residue_defs.nucleotide_slc);
            end;
            if ~isempty(slc),
                my_chain.seq(r) = slc;
            else
                my_chain.seq(r) = '?';
            end;
            for kat = 1:bbapr,
                adr = sprintf('%s{%i}%i.%s',address,m,resnum,bbat{kat});
                [msg,atcoor] = get_object(adr,'coor');
                if msg.error,
                    my_chain.mask(baspoi+1) = 0;
                    my_chain.resmask(r) = 0;
                else
                    coor(baspoi+kat,:) = atcoor;
                end;
            end;
        end;
    end;
    my_chain.resnr = maxres;
    my_chain.resnums = my_chain.resnums(1:maxres);
    my_chain.mask = my_chain.mask(1:bbapr*my_chain.resnr,1);
    my_chain.resmask = my_chain.resmask(1:my_chain.resnr,1);
    my_chain.seq = my_chain.seq(1:my_chain.resnr);
    coor = coor(1:bbapr*maxres,:);
    my_chain.coor{m} = coor;
end;
my_chain.align_point = my_chain.resnums(1);

if sum(my_chain.mask) == 0
    my_chain.seq = '';
    my_chain.coor = {};
end;

function info = analyze_rmsd(all_chains,chain,type)
% current version for only a single chain

switch type
    case 1
        bbapr = 3; % backbone atoms per residue
    otherwise
        bbapr = 6;
end;

my_chain = all_chains{chain};
info.resnums = my_chain.resnums;
info.mask = my_chain.mask;

rmsd_tab = zeros(my_chain.numrep);
weights = ones(1,bbapr*my_chain.resnr);

if my_chain.numrep == 1
    info.central = 1;
    info.central_rmsd = 0;
    info.rmsd_vec = zeros(1,my_chain.resnr);
    return
end;

for m1 = 1:my_chain.numrep-1
    coor1 = my_chain.coor{m1};
    coor1 = coor1(my_chain.mask == 1,:);
    for m2 = m1+1:my_chain.numrep
        coor2 = my_chain.coor{m2};
        coor2 = coor2(my_chain.mask == 1,:);
        rms = rmsd_superimpose(coor1,coor2,weights);
        rmsd_tab(m1,m2) = rms;
        rmsd_tab(m2,m1) = rms;
    end
end;

rmsd_min = 1e19;
best = 0;
for m = 1:my_chain.numrep
    rmsd = sqrt(sum(rmsd_tab(m,:).^2)/(my_chain.numrep-1));
    if rmsd < rmsd_min
        best = m;
        rmsd_min = rmsd;
    end
end

info.central = best;
info.central_rmsd = rmsd_min;

coor1 = my_chain.coor{best}; % central chain as reference
coor1 = coor1(my_chain.mask == 1,:);
weights0 = weights;
iter = 0;
det = 1;
resassign = zeros(bbapr*my_chain.resnr,1);
poi = 0;
for kr = 1:my_chain.resnr
    for kat = 1:bbapr
        poi = poi+1;
        resassign(poi) = kr;
    end
end
resassign = resassign(my_chain.mask == 1);
while det > 1e-4 && iter < 1000,
    iter = iter + 1;
    rmsd_vec = zeros(1,my_chain.resnr);
    msd_vec_all = zeros(1,bbapr*my_chain.resnr);
    for m = 1:my_chain.numrep
        if m ~= best
            coor2 = my_chain.coor{m};
            coor2 = coor2(my_chain.mask == 1,:);
            [~,coor2b] = rmsd_superimpose(coor1,coor2,weights);
            for r = 1:my_chain.resnr
                rcoor1 = coor1(resassign == r,:);
                rcoor2 = coor2b(resassign == r,:);
                [nat,~] = size(rcoor1);
                if nat > 0
                    rmsd = sum(sum((rcoor1-rcoor2).^2))/nat;
                    rmsd_vec(r) = rmsd_vec(r) + rmsd;
                    bas = bbapr*(r-1);
                    weights(bas+1:bas+bbapr) = rmsd;
                end;
            end;
        end;
    end
    weights = msd_vec_all;
    weights(weights == 0) = 0.025; % avoid numerical problems for coinciding atoms
    weights = 1./weights;
    % weights = exp(-weights.^2);
    weights(weights0 == 0) = 0;
    weights = weights/max(weights);
    det = sqrt(sum((weights0-weights).^2)/sum(weights.^2));
    weights0 = weights;
end;
% fprintf(1,'Iterations needed: %i\n',iter);
% fprintf(1,'Final rms weight change: %g\n',det);

info.rmsd_vec = sqrt(rmsd_vec/(my_chain.numrep-1));

function CData = get_detach_icon

CData(:,:,1) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,2) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,3) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

function CData = get_attach_icon

CData(:,:,1) =[    10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,2) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,3) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];


% --- Executes on button press in pushbutton_detach_plot.
function pushbutton_detach_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plot_attached,
    handles.plot_attached = false;
    set(hObject,'TooltipString','Attach original data plot');
    handles.plot_fig = figure('NumberTitle','off','Name','Flexibility plot','Units','pixels');
    set(hObject,'CData',handles.attach_icon);
    figpos = get(handles.plot_fig,'Position');
    pos = get(handles.axes_flexibility,'Position');
    handles.plot_pos = pos;
    set(handles.axes_flexibility,'Parent', handles.plot_fig);
    set(handles.axes_flexibility,'Units', 'pixels');
    set(handles.plot_fig,'CloseRequestFcn',@attach_plot);
    pos(1)= handles.margin+handles.add_margin;
    pos(2)= handles.margin+handles.add_margin;
    pos(3)= figpos(3)-2*handles.margin-handles.add_margin;
    pos(4)= figpos(4)-2*handles.margin-handles.add_margin;
    set(handles.axes_flexibility,'Position',pos);
    set(handles.plot_fig,'ResizeFcn',@plot_ResizeFcn);
    handles.plot_fig.UserData = hObject;
else
    handles.plot_attached = true;
    set(hObject,'TooltipString','Detach original data plot');
    set(hObject,'CData',handles.detach_icon);
    set(handles.axes_flexibility,'Parent', handles.figure1);
    set(handles.axes_flexibility,'Units', 'characters');
    set(handles.axes_flexibility,'Position', handles.plot_pos);
    delete(handles.plot_fig);
end;    
% Update handles structure
guidata(hObject, handles);

function attach_plot(hObject,callbackdata)

hMain = guidata(hObject.UserData);

hMain.plot_attached = true;
if ishandle(hMain.pushbutton_detach_plot),
    set(hMain.pushbutton_detach_plot,'CData',hMain.detach_icon);
    set(hMain.pushbutton_detach_plot,'TooltipString','Detach original data plot');
    set(hMain.axes_flexibility,'Parent', hMain.figure1);
    set(hMain.axes_flexibility,'Units', 'characters');
    set(hMain.axes_flexibility,'Position', hMain.plot_pos);
    guidata(hMain.pushbutton_detach_plot, hMain);
end;
delete(hObject);

% --- Executes when model_plot is resized.
function plot_ResizeFcn(hObject, callbackdata)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMain = guidata(hObject.UserData);

newpos=get(hObject,'Position');
newpos(1)= hMain.margin + hMain.add_margin;
newpos(2)= hMain.margin + hMain.add_margin;
newpos(3)= newpos(3)-2*hMain.margin - hMain.add_margin;
newpos(4)= newpos(4)-2*hMain.margin - hMain.add_margin;
set(hMain.axes_flexibility,'Position',newpos);



function edit_residue_number_offset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_residue_number_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_residue_number_offset as text
%        str2double(get(hObject,'String')) returns contents of edit_residue_number_offset as a double

[v,handles]=edit_update_MMM(handles,hObject,-10000,10000,0,'%i',1);
handles.resnum_offset = v;
plot_flexibility(handles,1,true);
for k = 1:length(handles.info)
    plot_flexibility(handles,k,false);
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_residue_number_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_residue_number_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_flexibility(handles,chain,initialize)

ref_mode = handles.checkbox_single_structure_reference.Value;
xray_ref = handles.popupmenu_chains.Value;

if ref_mode && chain == xray_ref
    return
end;

my_chain = handles.chains{1};
info = handles.info{1};
ref_align = my_chain.align_point;
ref_resnum = info.resnums(1);

my_chain = handles.chains{chain};
info = handles.info{chain};
align_offset = my_chain.align_point - ref_align;
resnum_offset = ref_resnum - info.resnums(1) + align_offset;
offset = resnum_offset + handles.resnum_offset;

rmsd_vec = info.rmsd_vec;
if ref_mode
    rmsd_vec = info.rmsd2xray;
end;

axes(handles.axes_flexibility); 
if initialize
    cla;
    hold on;
end
for k = 1:length(my_chain.resmask)-1
    if my_chain.resmask(k)*my_chain.resmask(k+1)
        plot(info.resnums(k:k+1)+offset,rmsd_vec(k:k+1),'Color',handles.plot_colors(chain,:));
    end;
end;
% plot(info.resnums(my_chain.resmask == 1)+offset,rmsd_vec(my_chain.resmask == 1),'k.');
% core_level = handles.rigid_core*min(rmsd_vec(info.rmsd_vec ~= 0));
% if core_level < handles.rigid_core_min
%     core_level = handles.rigid_core_min;
% end
% if ~isempty(core_level)
%     plot([info.resnums(1)+offset,info.resnums(end)+offset],[core_level,core_level],':','Color',handles.plot_colors(chain,:));
% end;
xlabel('Residue number');
ylabel('Local r.m.s.d. [?]');
guidata(handles.axes_flexibility,handles);


function save_flexibility(fname,handles)

my_chain = handles.chains{1};
info = handles.info{1};
ref_align = my_chain.align_point;
ref_resnum = info.resnums(1);

selected = handles.popupmenu_chains.String{handles.popupmenu_chains.Value};

for chain = 1:length(handles.chains),

    my_chain = handles.chains{chain};
    info = handles.info{chain};
    align_offset = my_chain.align_point - ref_align;
    resnum_offset = ref_resnum - info.resnums(1) + align_offset;
    offset = resnum_offset + handles.resnum_offset;
    
    flexibility(chain).name = my_chain.address;
    flexibility(chain).resnums = info.resnums + offset;
    flexibility(chain).rmsd = info.rmsd_vec;
    if isfield(info,'rmsd2xray')
        flexibility(chain).rmsd2xray = info.rmsd2xray;
    end;
end;

save(fname,'flexibility','selected');


% --- Executes on button press in checkbox_show_alignment.
function checkbox_show_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_alignment


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if isfield(handles,'plot_fig') && ishandle(handles.plot_fig),
    delete(handles.plot_fig);
end;
delete(handles.figure1);

function update_list(handles)

items = length(handles.popupmenu_chains.String);
for k = 1:items
    my_chain = handles.chains{k};
    tfield = sprintf('text_trace%i',k);
    handles.(tfield).String = handles.popupmenu_chains.String{k};
    if my_chain.numrep == 1
        handles.(tfield).FontAngle = 'italic';
        handles.(tfield).ForegroundColor = ([1,1,1] + handles.plot_colors(k,:))/2;
    else
        handles.(tfield).FontAngle = 'normal';
        handles.(tfield).ForegroundColor = handles.plot_colors(k,:);
    end
end;
for k = items+1:10
    tfield = sprintf('text_trace%i',k);
    handles.(tfield).String = '';
    handles.(tfield).ForegroundColor = handles.plot_colors(k,:);
end;
    
function handles = analyze_rmsd_to_xray(handles,ref)

switch handles.type
    case 1
        bbapr = 3; % backbone atoms per residue
    otherwise
        bbapr = 6;
end;


ref_chain = handles.chains{ref};
refcoor = ref_chain.coor{1};
ref_align = ref_chain.align_point;
ref_mask = ref_chain.mask;

refinfo = handles.info{ref};
numref = length(ref_chain.resnums);
refinfo.rmsd2xray = zeros(1,numref);
handles.info{ref} = refinfo;

for c = 1:length(handles.chains)
    if c~= ref
        info = handles.info{c};
        my_chain = handles.chains{c};
        my_align = my_chain.align_point;
        my_mask = my_chain.mask;
        nummy = length(my_chain.resnums);
        offset = ref_align - my_align;
        rmsd2xray = zeros(1,nummy);
        msdcount = zeros(1,nummy);
        ref_mask2 = zeros(size(ref_mask));
        my_mask2 = zeros(size(my_mask));
        aligned = 0;
        myindices = zeros(1,nummy);
        for rref = 1:numref
            basref = bbapr*(rref-1);
            rmy = rref + offset;
            if rmy >=1 && rmy <= nummy
                basmy = bbapr*(rmy-1);
                for at = 1:bbapr
                    if ref_mask(basref+at) && my_mask(basmy+at)
                        ref_mask2(basref+at) = 1;
                        my_mask2(basmy+at) = 1;
                        aligned = aligned + 1;
                        myindices(aligned) = basmy+at;
                    end;
                end;
            end;
        end;
        refcoor2 = refcoor((ref_mask2 == 1),:);
        for m = 1:my_chain.numrep
            mycoor = my_chain.coor{m};
            mycoor2 = mycoor((my_mask2 == 1),:);
            [~,mycoor2b] = rmsd_superimpose(refcoor2,mycoor2);
            tcoor = nan(size(mycoor));
            rcoor = nan(size(mycoor));
            for ka = 1:aligned,
                tcoor(myindices(ka),:) = mycoor2b(ka,:);
                rcoor(myindices(ka),:) = refcoor2(ka,:);
            end;
            msq = sum((tcoor-rcoor).^2,2);
            for rmy = 1:nummy
                basmy = bbapr*(rmy-1);
                msd = 0;
                for at = 1:bbapr
                    msd = msd + msq(basmy+at);
                end;
                if ~isnan(msd)
                    rmsd2xray(rmy) = rmsd2xray(rmy) + msd;
                    msdcount(rmy) = msdcount(rmy) + bbapr;
                end;
            end;
        end;
        info.rmsd2xray = sqrt(rmsd2xray./msdcount);
        handles.info{c} = info;
    end;
end;


% --- Executes on button press in checkbox_single_structure_reference.
function checkbox_single_structure_reference_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_single_structure_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_single_structure_reference

plot_flexibility(handles,1,1);
for k = 2:length(handles.info),
    plot_flexibility(handles,k,0);
end;
guidata(hObject, handles);
