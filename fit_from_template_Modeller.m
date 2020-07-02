function varargout = fit_from_template_Modeller(varargin)
% fit_from_template_Modeller M-file for fit_from_template_Modeller.fig
%      fit_from_template_Modeller, by itself, creates a new fit_from_template_Modeller or raises the existing
%      singleton*.
%
%      H = fit_from_template_Modeller returns the handle to a new fit_from_template_Modeller or the handle to
%      the existing singleton*.
%
%      fit_from_template_Modeller('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in fit_from_template_Modeller.M with the given input arguments.
%
%      fit_from_template_Modeller('Property','Value',...) creates a new fit_from_template_Modeller or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fit_from_template_Modeller_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fit_from_template_Modeller_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fit_from_template_Modeller

% Last Modified by GUIDE v2.5 09-Jul-2012 07:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fit_from_template_Modeller_OpeningFcn, ...
                   'gui_OutputFcn',  @fit_from_template_Modeller_OutputFcn, ...
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


% --- Executes just before fit_from_template_Modeller is made visible.
function fit_from_template_Modeller_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fit_from_template_Modeller (see VARARGIN)

% global MMM_icon
global hMain
global model

handles.darken=0.75; % darken chain colorscheme colors for better visibility and match with lighted ribbon model

handles.min_GA341=0.75;
handles.targ_ID='targ';
handles.temp_ID='temp';

handles.ensemble_size=10;
handles.models=20;

set(handles.edit_models,'String',sprintf('%i',handles.models));
set(handles.edit_ensemble_size,'String',sprintf('%i',handles.ensemble_size));

handles.DEER=[];
handles.direct=[];
handles.alignment=[];
handles.template_info.alg=0;
handles.template_info.name='';
handles.template_info.mode=-1;
handles.template_info.identity=1;
handles.target_info.alg=0;
handles.target_info.identity=1;
handles.translation=[];
handles.target = [];

% Choose default command line output for fit_from_template_Modeller
handles.output = hObject;

set(hMain.figure,'Pointer','watch');

drawnow;

h = msgbox('Please be patient. This can take some time.','Coarse-grained model is computed');

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

if ~exist('model','var') || ~isfield(model,'current_structure'),
    add_msg_board('Error: Cannot build from template, since no current structure exists.');
    guidata(hObject,handles);
    if ishandle(h),
        delete(h);
    end;
    set(hMain.figure,'Pointer','arrow');
    delete(hObject);
    return;
end;

handles.current_structure=model.current_structure;

hMain.auxiliary=[hMain.auxiliary hObject];
hMain.fit_plot=true;
hMain.fit_axes=handles.axes_model;

snum=model.current_structure;
adr=mk_address(snum);
set(handles.figure1,'Name',sprintf('Model transition from template structure %s',adr));

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
    if isempty(Ca_coor) || length(masses)<2,
        set(hMain.figure,'Pointer','arrow');
        add_msg_board('Error: Cannot build from template, since current structure has less than two amino acid residues.');
        if ishandle(h),
            delete(h);
        end;        
        guidata(hObject,handles);
        delete(hObject);
        return
    end;
    model.coarse(snum).Ca_coor=Ca_coor;
    model.coarse(snum).indices=rindices;
    model.coarse(snum).masses=masses;
    model.coarse(snum).Bfactors=Bfactors;
    model.coarse(snum).restypes=restypes;
    % make and store chain assignment
    [mr,nr]=size(rindices); 
    poi=1;
    chains=zeros(50,3);
    currchain=rindices(1,2);
    chains(1,1)=currchain;
    chains(1,2)=1;
    for k=2:mr,
        if rindices(k,2)~=currchain,
            currchain=rindices(k,2);
            chains(poi,3)=k-1;
            poi=poi+1;
            chains(poi,1)=currchain;
            chains(poi,2)=k;
        end;
    end;
    chains(poi,3)=mr;
    chains=chains(1:poi,:);
    model.coarse(snum).chains=chains;
else
    Ca_coor=model.coarse(snum).Ca_coor;
    rindices=model.coarse(snum).indices;
    restypes=model.coarse(snum).restypes;
    Bfactors=model.coarse(snum).Bfactors;
end;

axes(handles.axes_plot);
cla;
xlabel('Constraint');
ylabel('Distance (nm)');

axes(handles.axes_model);
hold on;

% determine residue ranges for chains
handles.chains=zeros(100,3);
poi=1;
handles.chains(1,1)=rindices(1,2);
handles.chains(1,2)=1;
[mm,nn]=size(rindices);
for k=1:mm,
    if rindices(k,2)~=handles.chains(poi,1),
        handles.chains(poi,3)=k-1;
        poi=poi+1;
        handles.chains(poi,1)=rindices(k,2);
        handles.chains(poi,2)=k;
    end;
end;
handles.chains(poi,3)=mm;
handles.chains=handles.chains(1:poi,:);
handles.wire=zeros(1,poi);

for k=1:poi,
    col=handles.darken*color_grade(handles.chains(k,1),poi);
    x=Ca_coor(handles.chains(k,2):handles.chains(k,3),1);
    y=Ca_coor(handles.chains(k,2):handles.chains(k,3),2);
    z=Ca_coor(handles.chains(k,2):handles.chains(k,3),3);
    handles.wire(k)=line(x,y,z,'color',col,'LineWidth',2);
end;

axis equal
axis off
hold on;

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
cam_up=get(hMain.axes_model,'CameraUpVector');
set(gca,'CameraPosition',cam_pos);
set(gca,'CameraTarget',cam_tar);
set(gca,'CameraUpVector',cam_up);
camlookat(handles.axes_model);

if ishandle(h),
    delete(h);
end;

set(hMain.figure,'Pointer','arrow');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fit_from_template_Modeller wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fit_from_template_Modeller_OutputFcn(hObject, eventdata, handles) 
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

entry=strcat(help_files,'fit_from_template_Modeller.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

hMain.fit_plot=false;
model.current_structure=handles.current_structure;
update_current_structure;
delete(handles.figure1);




% --------------------------------------------------------------------
function uitoggletool_chain_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_chain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

c=length(handles.wire);
for k=1:c,
    col=handles.darken*color_grade(handles.chains(k,1),c);
    set(handles.wire(k),'Color',col);
end;
guidata(hObject,handles);



% --------------------------------------------------------------------
function uipushtool_grey_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_grey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

c=length(handles.wire);
for k=1:c,
    set(handles.wire(k),'Color',[0.25,0.25,0.25]);
end;
guidata(hObject,handles);
 

% --- Executes on button press in pushbutton_decrease_GA341.
function pushbutton_decrease_GA341_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_decrease_GA341 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.min_GA341>0.5,
    handles.min_GA341=handles.min_GA341-0.05;
end;
set(handles.edit_min_GA341,'String',sprintf('%5.2f',handles.min_GA341));
guidata(hObject,handles);


% --- Executes on button press in pushbutton_increase_GA341.
function pushbutton_increase_GA341_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_increase_GA341 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.min_GA341<1,
    handles.min_GA341=handles.min_GA341+0.05;
    if handles.min_GA341>1,
        handles.min_GA341=1;
    end;
end;
set(handles.edit_min_GA341,'String',sprintf('%5.2f',handles.min_GA341));
guidata(hObject,handles);



function edit_min_GA341_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_GA341 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_GA341 as text
%        str2double(get(hObject,'String')) returns contents of edit_min_GA341 as a double

[v,handles]=edit_update_MMM(handles,hObject,0.3,1.00,0.75,'%5.2f');
handles.min_GA341=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_min_GA341_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_GA341 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkbox_background.
function checkbox_background_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_background


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% in generation of an ensemble, all direct C_alpha-C_alpha constraints are
% always used, only DEER constraints can be ommitted, a warning is issued
% if the number of DEER constraints is insufficient to generate the
% specified ensemble

global model
global general
global third_party
global help_files

if isempty(handles.DEER) && isempty(handles.direct),
    if isempty(handles.alignment) || handles.template_info.alg==handles.target_info.alg,
        add_msg_board('ERROR: No experimental restraints and target sequence identical to template sequence.');
        add_msg_board('Modeling is aborted.');
        return
    end;
end;

if handles.target_info.identity<0.2,
    add_msg_board('### Strong warning ###: Insufficient identity of target and template sequence.');
    add_msg_board('Homology modeling is most likely unreliable.');
elseif handles.target_info.identity<0.4,
    add_msg_board('Warning: Identity of target and template sequence lower than 40%%.');
    add_msg_board('Homology modeling may be unreliable.');
end;

entry=strcat(help_files,'third_party.html#Modeller');

dospath=which([third_party.modeller_version  '.exe']);
if isempty(dospath),
    message.error=2;
    message.text='Modeller software not found on Matlab path.';
    add_msg_board('This feature requires Modeller from the Sali lab');
    add_msg_board('ERROR: Modeller could not be found on the Matlab path.');
    add_msg_board('Please check whether Modeller is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);

handles.test_mode=false;
if ~isempty(handles.target),
    snum=model.current_structure;
    [coarse_correspondence,coarse_target]=align_template_target(mk_address(snum),handles.target);
    if ~isempty(coarse_correspondence),
        handles.test_mode=true;
        if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
            [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
            if isempty(Ca_coor) || length(masses)<2,
                add_msg_board('Error: Cannot generate test information, since current structure has less than two amino acid residues.');
                handles.test_mode=false;
            else
                handles.coarse_template=Ca_coor;
                model.coarse(snum).Ca_coor=Ca_coor;
                model.coarse(snum).indices=rindices;
                model.coarse(snum).masses=masses;
                model.coarse(snum).Bfactors=Bfactors;
                model.coarse(snum).restypes=restypes;
            end;
        else
            handles.coarse_template=model.coarse(snum).Ca_coor;
            handles.coarse_target=coarse_target;
            handles.coarse_correspondence=coarse_correspondence;
        end;
    end;
end;

hwin=gcf;
set(hwin,'Pointer','watch');

idCode=model.info{model.current_structure}.idCode;
if isempty(idCode), idCode='AMMM'; end;
idCode(1)=char(idCode(1)+16);
filename=sprintf('%s.pdb',handles.temp_ID);

fname=fullfile(general.tmp_files, filename);
msg=sprintf('Template structure saved as PDB file: %s',fname);
add_msg_board(msg);
if handles.template_info.mode>0,
    [message,info]=wr_pdb_paradigm(fname,idCode,1,handles.template_info.mode);
else
    [message,info]=wr_pdb_paradigm(fname,idCode);
end;
if message.error,
    add_msg_board(message.text);
    set(hwin,'Pointer','arrow');
    return
end;
for k=1:length(info.offsets),
    handles.translation(k).targ_offset=info.offsets(k);
    if ~isempty(handles.alignment) && isfield(handles.alignment(1),'first_res'),
        handles.translation(k).targ_offset=...
            handles.translation(k).targ_offset+str2double(handles.alignment(1).first_res)-1;
    end;
    if isfield(handles,'restraints') && isfield(handles.restraints,'chains'),
        targ_tag=id2tag(k,handles.restraints.chains);
        if isempty(targ_tag),
            add_msg_board('ERROR: More chains in template file than specified in restraint file.');
            return
        end;
    else
        targ_tag=char(double('A')-1+k);
    end;    
    handles.translation(k).tchain=targ_tag;
end;
if length(info.offsets)==1,
    handles.translation(1).tchain='';
end;
temp_seq=info.seq;
targ_seq='';
arestarg=[];
if length(handles.alignment)>1,
    ali(1)=handles.alignment(handles.template_info.alg);
    ali(2)=handles.alignment(handles.target_info.alg);
    if isfield(ali,'first_res'),
        offset=str2double(ali(1).first_res)-info.first;
        if offset<0,
            add_msg_board(sprintf('Warning: Eliminating from alignment %i N-terminal residues that are missing in template',-offset));
            ali(1).first_res=sprintf('%i',info.first);
        else
            info.first=info.first+offset;
        end;
    end;
    if isfield(ali,'first_res') && str2double(ali(2).first_res)>0,
        arestarg=str2double(ali(2).first_res);
    else
        arestarg=1;
    end;
    [temp_seq,targ_seq,darestarg,clean_targ_seq]=adjust_alignment(ali,temp_seq);
    if isfield(ali,'first_res')
        info.first=str2double(ali(1).first_res);
        if isfield(ali,'last_res')
            info.last=str2double(ali(1).last_res);
        end;
    end;
    arestarg=arestarg+darestarg;
    handles.translation(k).targ_offset=1-arestarg;
    handles.alignment(handles.target_info.alg).first_res=sprintf('%i',arestarg);
else
    if ~isempty(handles.alignment),
        clean_targ_seq=handles.alignment(handles.template_info.alg).sequence;
    else
        clean_targ_seq=handles.clean_targ_seq(1-handles.translation(k).targ_offset:end);
        if length(clean_targ_seq)>length(temp_seq),
            clean_targ_seq=clean_targ_seq(1:length(temp_seq));
        end;
        targ_seq=clean_targ_seq;
        for k=1:length(targ_seq),
            if char(temp_seq(k))=='!' || char(temp_seq(k))==';',
                targ_seq(k)=temp_seq(k);
            end;
        end;
        arestarg=info.first;
        handles.translation(1).targ_offset=1-arestarg;
        handles.translation(1).chain=info.first_chain;
        % arestarg=1;
%        handles.alignment(handles.target_info.alg).first_res=sprintf('%i',arestarg);
    end;
end;
handles.ares=arestarg;
handles.clean_targ_seq=clean_targ_seq;
basname=sprintf('%s_to_%s',handles.temp_ID,handles.targ_ID); % ###
[alignfile,message]=seq2pir(basname,handles.temp_ID,handles.targ_ID,temp_seq,info.first,info.first_chain,info.last,info.last_chain,targ_seq,arestarg);

if message.error,
    add_msg_board(message.text);
    add_message_board('ERROR: Modelling aborted.');
    set(hwin,'Pointer','arrow');
    return
end;

add_msg_board(sprintf('Alignment file saved as: %s',alignfile));
[algpath, algfile, algext] = fileparts(alignfile);
alignfile=strcat(algfile,algext);

runfilename=sprintf('model_%s_to_%s.py',handles.temp_ID,handles.targ_ID);
runfile=fullfile(general.tmp_files, runfilename);
if ~isfield(handles,'restraints'),
    handles.restraints=[];
end;
message=mk_modeller_input(runfile,alignfile,handles.temp_ID,handles.targ_ID,handles.restraints,handles.translation,handles.models);

if message.error,
    add_msg_board(message.text);
    add_message_board('ERROR: Modelling aborted.');
    set(hwin,'Pointer','arrow');
    return
end;

msg='Running Modeller. Please be patient...';
is_bckg=false;
if get(handles.checkbox_background,'Value'),
    is_bckg=true;
    msg='Started Modeller run as a background job.';
end;


[batcmd,message]=mk_modeller_bat(modpath,modcmd,runfile,is_bckg);
if message.error,
    add_msg_board('ERROR: Modeller could not be initialized');
    add_msg_board(message.text);
    set(hwin,'Pointer','arrow');
    return
end;

set(handles.text_info,'String',msg);
add_msg_board(msg);
drawnow;
my_dir=pwd;
cd(modpath);
[s, w] = dos(batcmd);
if is_bckg,
    if length(handles.alignment)>1,
        ali=handles.alignment(handles.target_info.alg);
        if isfield(ali,'first_res') && str2double(ali.first_res)>0,
            ares=str2double(ali.first_res);
        end;
    else
        ares=1;
    end;
    if isfield(handles.restraints,'chains'),
        chain_tags=handles.restraints.chains;
    else
        chain_tags='';
    end;
    job.type='Modeller';
    job.started=clock;
    job.name=sprintf('job-%s',datestr(job.started,'yyyy-mm-dd_HH-MM-SS'));
    job.log=sprintf('model_%s_to_%s.log',handles.temp_ID,handles.targ_ID);
    job.finished='Total CPU time [seconds]';
    job.min_GA341=handles.min_GA341;
    job.ensemble_size=handles.ensemble_size;
    job.targ_ID=handles.targ_ID;
    job.ares=ares;
    job.chain_tags=chain_tags;
    job.clean_targ_seq=clean_targ_seq;
    t = timer('TimerFcn',@check_for_Modeller_completion, 'StartDelay', 120.0, 'Period', 120.0, 'ExecutionMode', 'fixedSpacing' ,'UserData', job, 'Tag', 'MMM', 'Name', job.name);
    job.timer=t;
    k=length(general.timers);
    general.timers{k+1}=t;
    save(fullfile(general.tmp_files,job.name),'job');
    start(t);
    set(hwin,'Pointer','arrow');
    cd(my_dir);
    return
end;
if s~=0,
    rem=w;
    while ~isempty(rem),
        [token,rem]=strtok(rem,char(10));
        if ~isempty(token),
            add_msg_board(token);
        end;
    end;
    message.error=2;
    message.text='Modeller error.';
    add_msg_board('ERROR: Modeller did not run successfully.');
    set(hwin,'Pointer','arrow');
    cd(my_dir);
    return
else
    logname=sprintf('model_%s_to_%s.log',handles.temp_ID,handles.targ_ID);
    msg=sprintf('Modeller job %s_to_%s completed.',handles.temp_ID,handles.targ_ID);
    add_msg_board(msg);
    set(handles.text_info,'String',msg);
    analyze_modeller_log(logname,true);
    set(handles.pushbutton_import,'UserData',logname);
end;
cd(my_dir);


set(hwin,'Pointer','arrow');

set(handles.pushbutton_import,'Enable','on');
guidata(hObject,handles);


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global ENM_param
global general
global hMain

snum=model.current_structure;

handles.clean_targ_seq=model.structures{snum}(1).sequence;
axes(handles.axes_model);

my_path=pwd;
cd(general.restraint_files);

[fname,pname,findex]=uigetfile('*.dat','Load restraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Constraint loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    hfig=gcf;
    set(hfig,'Pointer','watch');
    handles.basis=ENM_param.fit_basis;
    restraints=rd_restraints(fullfile(pname,fname));

    handles.ensemble=restraints.ensemble;
    handles.uncertainty=restraints.uncertainty;
    handles.exclude=restraints.exclude;
    handles.target=restraints.target;
    if isfield(restraints,'basis'),
        handles.basis=restraints.basis;
    end;

    new_template=false;
    if isfield(restraints,'PDB'),
        if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
            id=tag2id(restraints.PDB,model.structure_tags);
            if ~isempty(id),
                snum=id;
                model.current_structure=snum;
            else
                button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to load specified template?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
                if strcmp(button,'Yes'),
                    new_template=true;
                    fname=get_pdb_file(restraints.PDB);
                    [message,snum]=add_pdb(fname);
                    adr=mk_address(snum);
                    set(handles.figure1,'Name',sprintf('Model transition from template structure %s',adr));
                    if message.error,
                        add_msg_board(sprintf('ERROR: Specified template PDB file %s could not be retrieved from server',restraints.PDB));
                        add_msg_board(message.text);
                        set(hfig,'Pointer','arrow');
                        return
                    else
                        model.current_structure=snum;
                    end;
                end;
            end;
        end;
    end;

    if isfield(restraints,'template') && length(restraints.template)>1,
        for k=1:length(restraints.template),
            idCode=restraints.template(k).pdbid;
            id=tag2id(upper(idCode),upper(model.structure_tags));
            if isempty(id),
                new_template=true;
                add_msg_board(sprintf('Retrieving template PDB file %s from server',idCode));
                fname=get_pdb_file(idCode);
                [message,snum]=add_pdb(fname);
                if message.error,
                    add_msg_board(sprintf('ERROR: Specified template PDB file %s could not be retrieved from server',idCode));
                    add_msg_board(message.text);
                    set(hfig,'Pointer','arrow');
                    return
                elseif k==1,
                    snum1=snum;
                    model.current_structure=snum;
                    adr=mk_address(snum);
                    set(handles.figure1,'Name',sprintf('Model transition from template structure %s',adr));
                end;
            elseif k==1,
                snum1=id;
            end;
        end;
        model.current_structure=snum1;
    end;

    if isfield(restraints,'alignment'),
        [pstr,fname,ext]=fileparts(restraints.alignment);
        switch lower(ext)
            case {'.afa','.fa','.fasta','.ebi','.txt'}
                alignment=get_multiple_fasta(restraints.alignment);
            case {'.ali','.pir'}
                alignment=get_multiple_pir(restraints.alignment);
            case {'.clw','.aln','.clustalw'}
                alignment=get_multiple_clustal(restraints.alignment);
            otherwise
                add_msg_board('Warning: Unknown sequence format. Assuming FASTA.');
                alignment=get_multiple_fasta(restraints.alignment);
        end;
        if isempty(alignment),
            add_msg_board(sprintf('ERROR: Alignment file %s was specified in restraint list',restraints.alignment));
            add_msg_board('This file could not be read or alignment format was not recognized.');
            set(hfig,'Pointer','arrow');
            return
        end;
        alignment=test_repair_alignment(alignment);
        set(gcf,'Pointer','watch');
        drawnow;
        if ~isempty(alignment),
            handles.alignment=alignment;
            [best_match,name,mode,identity]=match_alignment(alignment);
            if identity<0.95,
                add_msg_board('ERROR: Template structure does not seem to exist in alignment file.');
                set(gcf,'Pointer','arrow');
                return
            end;
            if length(alignment)~=2,
                add_msg_board('ERROR: This mode requires exactly two sequences in the alignment file.');
                set(gcf,'Pointer','arrow');
                return
            end;
            handles.template_info.alg=best_match;
            handles.template_info.name=name;
            handles.template_info.mode=mode;
            handles.template_info.identity=identity;
            add_msg_board(sprintf('Sequence %i of %s in the alignment array best matches template structure',best_match,name));
            if mode==0,
                add_msg_board('This sequence matches all chains of the template');
            else
                ctag=id2tag(mode,char(model.chain_tags{model.current_structure}));
                add_msg_board(sprintf('This sequence matches chain %s of the template',ctag));
            end;
            add_msg_board(sprintf('Sequence identity of match with template is %4.1f%%',100*identity));
            if best_match==1,
                handles.target_info.alg=2;
            else
                handles.target_info.alg=1;
            end;
            ali(1)=alignment(handles.target_info.alg);
            ali(2)=alignment(handles.template_info.alg);
            handles.target_info.identity=get_identity(ali);
            msg=sprintf('Selected target sequence %s with %4.1f%% identity to template.',...
                alignment(handles.target_info.alg).name,...
                100*handles.target_info.identity);
            add_msg_board(msg);
            ident=handles.target_info.identity;
            if ident>0.4,
                col=[0,1,0];
            elseif ident<0.2,
                col=[1,0,0];
            else
                twilight=(ident-0.2)/0.2;
                col=twilight*[0,1,0]+(1-twilight)*[1,0,0];
            end;
            set(handles.pushbutton_alignment,'ForegroundColor',col);
        else
            add_msg_board('ERROR: Alignment file is invalid. No sequence found.');
            return
        end;
        set(gcf,'Pointer','arrow');
        drawnow;
        ali(1)=alignment(handles.template_info.alg);
        ali(2)=alignment(handles.target_info.alg);
        handles.alignment=ali;
        handles.template_info.alg=1;
        handles.target_info.alg=2;
    end;
    
    if ~isempty(handles.alignment),
        if restraints.realign,
            handles.alignment=realign(handles.alignment,restraints);
        end;
        handles.translation=trans_table(handles.alignment,restraints);
    else
        if restraints.realign,
            add_msg_board(sprintf('Warning: Constraint file requests sequence realignment, but no alignment was loaded or specified.'));
        end;
        handles.translation=[];
    end;
    [DEER,cancelled]=process_DEER_restraints(handles,restraints,handles.alignment);
    snum=model.current_structure;
    if cancelled,
        add_msg_board('Processing of DEER constraints cancelled.');
        set(hfig,'Pointer','arrow');
        return
    end;
    handles.DEER=DEER;
    restraints.DEER=DEER;
    [direct,cancelled]=process_direct_restraints(restraints,handles.alignment);
    if cancelled,
        add_msg_board('Processing of direct constraints cancelled.');
        set(hfig,'Pointer','arrow');
        return
    end;
    handles.direct=direct;
    handles.restraints=restraints;
    
    c=length(handles.wire);
    for k=1:c,
        set(handles.wire(k),'Color',[0.25,0.25,0.25]);
    end;
    handles.restraint_graphics=hgtransform;
    if ~isempty(handles.DEER),
        dvec=zeros(length(handles.DEER),3);
        for k=1:length(handles.DEER),
            handles.DEER(k).l1=line(handles.DEER(k).xyz1(1),handles.DEER(k).xyz1(2),handles.DEER(k).xyz1(3),...
                'Parent',handles.restraint_graphics,'Marker','.','Color','b');
            handles.DEER(k).l2=line(handles.DEER(k).xyz2(1),handles.DEER(k).xyz2(2),handles.DEER(k).xyz2(3),...
                'Parent',handles.restraint_graphics,'Marker','.','Color','b');
            x=[handles.DEER(k).xyz1(1) handles.DEER(k).xyz2(1)];
            y=[handles.DEER(k).xyz1(2) handles.DEER(k).xyz2(2)];
            z=[handles.DEER(k).xyz1(3) handles.DEER(k).xyz2(3)];
            r0=norm(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/10;
            dvec(k,:)=(handles.DEER(k).r-r0)*(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/(10*r0);
            det=abs(r0-handles.DEER(k).r)/handles.DEER(k).sigr;
            if det>2,
                col='r';
            elseif det>1,
                col=[255,190,0]/255; % yellow-orange
            else
                col='g';
            end;
            handles.DEER(k).ll=line(x,y,z,'Parent',handles.restraint_graphics,'Color',col,'LineWidth',1.5,'LineStyle',':');
            cindices=handles.DEER(k).indices;
            f1=false;
            f2=false;
            for l=1:length(model.coarse(snum).indices),
                diff=cindices(1,:)-model.coarse(snum).indices(l,:);
                if sum(abs(diff))==0,
                    handles.DEER(k).res1=l;
                    f1=true;
                    x=[handles.DEER(k).xyz1(1) model.coarse(snum).Ca_coor(l,1)];
                    y=[handles.DEER(k).xyz1(2) model.coarse(snum).Ca_coor(l,2)];
                    z=[handles.DEER(k).xyz1(3) model.coarse(snum).Ca_coor(l,3)];
                    handles.DEER(k).rl1=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
                end;
                diff=cindices(2,:)-model.coarse(snum).indices(l,:);
                if sum(abs(diff))==0,
                    handles.DEER(k).res2=l;
                    f2=true;
                    x=[handles.DEER(k).xyz2(1) model.coarse(snum).Ca_coor(l,1)];
                    y=[handles.DEER(k).xyz2(2) model.coarse(snum).Ca_coor(l,2)];
                    z=[handles.DEER(k).xyz2(3) model.coarse(snum).Ca_coor(l,3)];
                    handles.DEER(k).rl2=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
                end;
            end;
            if ~f1,
                add_msg_board(sprintf('Warning: Residue for first label of restraint %i not in coarse-grained model.',k));
            end;
            if ~f2,
                add_msg_board(sprintf('Warning: Residue for second label of restraint %i not in coarse-grained model.',k));
            end;
        end;
        dvec=dvec/sqrt(sum(sum(dvec.^2)));
    end;

    if new_template,
        handles = update_3D_model(handles,snum,handles.DEER);
    end;
end;
if ~isempty(handles.DEER),
    filename=sprintf('%s_restraint_matching.txt',id2tag(snum,model.structure_tags));
    fname=fullfile(general.tmp_files, filename);
    fid=fopen(fname,'wt');
    if fid==-1,
        add_msg_board('Constraint matching report file could not be written');
    else
        fprintf(fid,'Constraint matching report for structure %s\n',id2tag(snum,model.structure_tags));
        fprintf(fid,'(all distances in nm)\n\n');
    end;
    axes(handles.axes_plot);
    cla;
    set(handles.text_auxiliary,'String','Constraint matching');
    ma=0;
    rmsd=0;
    fom=0; % figure of merit
    if fid~=-1,
        fprintf(fid,'--- Constraint matching for template structure ---\n\n');
        fprintf(fid,'Resid. 1  Resid. 2   r(exp)  sr(exp) r(mod)  |dr/sr|  Matching\n');
    end;
    for k=1:length(handles.DEER),
        r0=norm(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/10;
        rmsd=rmsd+(r0-handles.DEER(k).r)^2;
        det=abs(r0-handles.DEER(k).r)/handles.DEER(k).sigr;
        fom=fom+det^2;
        if det>2,
            col='r';
            match_str='### not matched ###';
        elseif det>1,
            col=[255,190,0]/255; % yellow-orange
            match_str='*** poorly matched ***';
        else
            col='g';
            match_str='matched';
        end;
        if fid~=-1,
            fprintf(fid,'%8s  %8s%8.2f%8.2f%8.2f %8.2f   %s\n',DEER(k).adr1,DEER(k).adr2,DEER(k).r,DEER(k).sigr,r0,det,match_str);
        end;
        errorbar(k,handles.DEER(k).r,handles.DEER(k).sigr,'k');
        if handles.DEER(k).r+handles.DEER(k).sigr>ma,
            ma=handles.DEER(k).r+handles.DEER(k).sigr;
        end;
        if r0>ma,
            ma=r0;
        end;
        hold on;
        plot(k,r0,'.','Color',col);
    end;
    handles.rmsd=sqrt(rmsd/length(handles.DEER));
    handles.fom=sqrt(fom/length(handles.DEER));
    axis([0,length(handles.DEER)+1,0,1.05*ma]);
    xlabel('Constraint number');
    ylabel('Distance (nm)');
    set(handles.text_info,'String',sprintf('Loaded %i DEER restraints with rmsd of %5.2f nm to template',length(handles.DEER),handles.rmsd));
    set(handles.text_auxiliary_msg,'String',sprintf('Figure of merit: %6.4f',handles.fom));
    if fid~=-1,
        fclose(fid);
        hMain.report_file=fname;
        report_editor;
    end;
end;
set(handles.pushbutton_run,'Enable','on');
set(handles.pushbutton_alignment,'Enable','off');
set(hfig,'Pointer','arrow');
cd(my_path);
guidata(hObject,handles);

function [DEER,cancelled]=process_DEER_restraints(handles,restraints,alignment)

global model
global hMain

cancelled=false;

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

snum=model.current_structure;

if get(handles.checkbox_dynamic_rotamers,'Value'),
    hMain.dynamic_rotamers=true;
else
    hMain.dynamic_rotamers=false;
end;


restraints=label_all_sites(restraints,handles.translation,true);

labels=label_information(model.sites);

for k=1:length(restraints.DEER),
    adr1=restraints.DEER(k).adr1;
    tempadr1=translate_address(adr1,handles.translation);
    ind1=resolve_address(tempadr1);
    adr2=restraints.DEER(k).adr2;
    tempadr2=translate_address(adr2,handles.translation);
    ind2=resolve_address(tempadr2);
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).indices=[ind1;ind2];
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
    DEER(k).label=restraints.DEER(k).label;
    DEER(k).T=restraints.DEER(k).T;
    f1=false;
    f2=false;
    for l=1:length(labels),
        diff1=ind1-labels(l).indices;
        if sum(abs(diff1))==0,
            f1=true;
            DEER(k).xyz1=labels(l).xyz;
            DEER(k).rmsd1=labels(l).rmsd;
        end;
        diff2=ind2-labels(l).indices;
        if sum(abs(diff2))==0,
            f2=true;
            DEER(k).xyz2=labels(l).xyz;
            DEER(k).rmsd2=labels(l).rmsd;
        end;
    end;
    if ~f1 || ~f2,
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
        cancelled=true;
        DEER=[];
        return;
    end;
end;


function [direct,cancelled]=process_direct_restraints(restraints,alignment)

global model

cancelled=false;
if ~isfield(restraints,'direct'),
    direct=[];
    return;
end;

snum=model.current_structure;

if ~isempty(alignment),
end;

md=length(restraints.direct);
direct=zeros(md,5);

cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);
poi=0;
add_msg_board(sprintf('Checking and processing %i direct C_alpha-C_alpha restraints',md));
for k=1:md,
    possible=true;
    adr1=restraints.direct(k).adr1;
    ind1=resolve_address(adr1);
    if isempty(ind1),
        possible=false;
        add_msg_board(sprintf('Warning: Specified residue %s does not exist in template structure.',adr1));
        continue
    end;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    adr2=restraints.direct(k).adr2;
    ind2=resolve_address(adr2);
    if isempty(ind2),
        possible=false;
        add_msg_board(sprintf('Warning: Specified residue %s does not exist in template structure.',adr2));
        continue;
    end;
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    if net1==0 || net2==0,
        possible=false;
        add_msg_board(sprintf('Warning: Constrained C_alpha-C_alpha distance %s to %s does not exist in template structure.',adr1,adr2));
    end;
    if possible,
        poi=poi+1;
        direct(poi,1)=net1;
        direct(poi,2)=net2;    
        direct(poi,4)=restraints.direct(k).r;
        direct(poi,5)=restraints.direct(k).sigr;
    end;
end;
if poi<md,
    add_msg_board(sprintf('Warning: %i out of %i direct restraints had to be removed.',md-poi,md));
else
    add_msg_board('All direct restraints can be kept.');
end;
direct=direct(1:poi,:);

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
        end;
    end;
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for ? -> nm


function [handles, all_constraint_rmsd] = update_3D_model(handles,snum,DEER)

global model 
global general
global hMain

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(mk_address(snum));
    model.coarse(snum).Ca_coor=Ca_coor;
    model.coarse(snum).indices=rindices;
    model.coarse(snum).masses=masses;
    model.coarse(snum).Bfactors=Bfactors;
    model.coarse(snum).restypes=restypes;
end;

modnum=length(model.structures{snum}(1).xyz);

axes(handles.axes_model);

[network,masses,rindices]=coarse_residues(mk_address(snum));

% determine residue ranges for chains
handles.chains=zeros(100,3);
poi=1;
handles.chains(1,1)=rindices(1,2);
handles.chains(1,2)=1;
[mm,nn]=size(rindices);
for k=1:mm,
    if rindices(k,2)~=handles.chains(poi,1),
        handles.chains(poi,3)=k-1;
        poi=poi+1;
        handles.chains(poi,1)=rindices(k,2);
        handles.chains(poi,2)=k;
    end;
end;
handles.chains(poi,3)=mm;
handles.chains=handles.chains(1:poi,:);
handles.wire=zeros(1,poi);

for k=1:poi,
    col=handles.darken*color_grade(handles.chains(k,1),poi);
    x=network(handles.chains(k,2):handles.chains(k,3),1);
    y=network(handles.chains(k,2):handles.chains(k,3),2);
    z=network(handles.chains(k,2):handles.chains(k,3),3);
    handles.wire(k)=line(x,y,z,'color',col,'LineWidth',2);
end;


% update N-O midpoint coordinates for all label pairs
for k=1:length(DEER),
    lindices=DEER(k).indices;
    res1adr=mk_address([snum lindices(1,2) 1 lindices(1,4)]);
    [msg,N1coor]=get_object([res1adr '.N1'],'coor');
    [msg,O1coor]=get_object([res1adr '.O1'],'coor');
    DEER(k).xyz1=(N1coor+O1coor)/2;
    res2adr=mk_address([snum lindices(2,2) 1 lindices(2,4)]);
    [msg,N1coor]=get_object([res2adr '.N1'],'coor');
    [msg,O1coor]=get_object([res2adr '.O1'],'coor');
    DEER(k).xyz2=(N1coor+O1coor)/2;
end;

chains=length(handles.wire);
for c=1:chains,
    ia=handles.chains(c,2);
    ie=handles.chains(c,3);
    x=network(ia:ie,1);
    y=network(ia:ie,2);
    z=network(ia:ie,3);
    set(handles.wire(c),'XData',x,'YData',y,'ZData',z);
end;
for k=1:length(DEER),
    set(DEER(k).l1,'XData',DEER(k).xyz1(1),'YData',DEER(k).xyz1(2),'ZData',DEER(k).xyz1(3));
    set(DEER(k).l2,'XData',DEER(k).xyz2(1),'YData',DEER(k).xyz2(2),'ZData',DEER(k).xyz2(3));
    x=[DEER(k).xyz1(1) DEER(k).xyz2(1)];
    y=[DEER(k).xyz1(2) DEER(k).xyz2(2)];
    z=[DEER(k).xyz1(3) DEER(k).xyz2(3)];
    r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
    det=abs(r0-DEER(k).r)/DEER(k).sigr;
    if det>2,
        col='r';
    elseif det>1,
        col=[255,190,0]/255; % yellow-orange
    else
        col='g';
    end;
    set(DEER(k).ll,'XData',x,'YData',y,'ZData',z,'Color',col);
    cindices=DEER(k).indices;
    cindices(1,1)=snum;
    cindices(2,1)=snum;
    cindices(1,3)=1;
    cindices(2,3)=1;
    f1=false;
    f2=false;
    for l=1:length(model.coarse(snum).indices),
        diff=cindices(1,:)-model.coarse(snum).indices(l,:);
        if sum(abs(diff))==0,
            f1=true;
            x=[DEER(k).xyz1(1) network(l,1)];
            y=[DEER(k).xyz1(2) network(l,2)];
            z=[DEER(k).xyz1(3) network(l,3)];
            set(DEER(k).rl1,'XData',x,'YData',y,'ZData',z);
        end;
        diff=cindices(2,:)-model.coarse(snum).indices(l,:);
        if sum(abs(diff))==0,
            f2=true;
            x=[DEER(k).xyz2(1) network(l,1)];
            y=[DEER(k).xyz2(2) network(l,2)];
            z=[DEER(k).xyz2(3) network(l,3)];
            set(handles.DEER(k).rl2,'XData',x,'YData',y,'ZData',z);
        end;
    end;
    if ~f1,
        add_msg_board(sprintf('Warning: Residue for first label of restraint %i not in network model.',k));
    end;
    if ~f2,
        add_msg_board(sprintf('Warning: Residue for second label of restraint %i not in network model.',k));
    end;
end;
camlookat(handles.axes_model);
axes(handles.axes_plot);
cla;

filename=sprintf('%s_restraint_matching.txt',id2tag(snum,model.structure_tags));
fname=fullfile(general.tmp_files, filename);
fid=fopen(fname,'wt');
if fid==-1,
    add_msg_board('Constraint matching report file could not be written');
else
    fprintf(fid,'Constraint matching report for structure %s\n',id2tag(snum,model.structure_tags));
    fprintf(fid,'(all distances in nm)\n');
end;

set(handles.text_auxiliary,'String','Constraint matching');
full_rmsd=0;
full_num=0;
all_constraint_rmsd=zeros(1,modnum);
for km=1:modnum,
    if fid~=-1,
        fprintf(fid,'\n--- Constraint matching for model %i ---\n\n',km);
        fprintf(fid,'Resid. 1  Resid. 2   r(exp)  sr(exp) r(mod)  |dr/sr|  Matching\n');
    end;
    ma=0;
    rmsd=0;
    num=0;
    fom=0; % figure of merit
    if ~isempty(DEER),
        for k=1:length(DEER),
            lindices=DEER(k).indices;
            res1adr=mk_address([snum lindices(1,2) km lindices(1,4)]);
            [stag,ctag1,modelnum,resnum1]=mk_address_parts([snum lindices(1,2) km lindices(1,4)]);
            r1adrshort=sprintf('(%s)%i',ctag1,resnum1);
            [msg,N1coor]=get_object([res1adr '.N1'],'coor');
            [msg,O1coor]=get_object([res1adr '.O1'],'coor');
            DEER(k).xyz1=(N1coor+O1coor)/2;
            res2adr=mk_address([snum lindices(2,2) km lindices(2,4)]);
            [stag,ctag2,modelnum,resnum2]=mk_address_parts([snum lindices(2,2) km lindices(2,4)]);
            r2adrshort=sprintf('(%s)%i',ctag2,resnum2);
            [msg,N1coor]=get_object([res2adr '.N1'],'coor');
            [msg,O1coor]=get_object([res2adr '.O1'],'coor');
            DEER(k).xyz2=(N1coor+O1coor)/2;
            r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
            rmsd=rmsd+(r0-DEER(k).r)^2;
            num=num+1;
            full_num=full_num+1;
            det=abs(r0-DEER(k).r)/DEER(k).sigr;
            fom=fom+det^2;
            if det>2,
                col='r';
                match_str='### not matched ###';
            elseif det>1,
                col=[255,190,0]/255; % yellow-orange
                match_str='*** poorly matched ***';
            else
                col='g';
                match_str='well matched';
            end;
            if fid~=-1,
                fprintf(fid,'%10s%10s%8.2f%8.2f%8.2f %8.2f    %s\n',r1adrshort,r2adrshort,DEER(k).r,DEER(k).sigr,r0,det,match_str);
            end;
            errorbar(k,DEER(k).r,DEER(k).sigr,'k');
            if DEER(k).r+DEER(k).sigr>ma,
                ma=DEER(k).r+DEER(k).sigr;
            end;
            if r0>ma,
                ma=r0;
            end;
            hold on;
            plot(k,r0,'.','Color',col);
        end;
        if num==0, num=1; end;
        fprintf(fid,'\nMean r.m.s.d. of all restraints for model %i: %5.2f nm\n',km,sqrt(rmsd/num));
        all_constraint_rmsd(km)=sqrt(rmsd/num);
    end;
    if ma>0 && length(DEER)>0,
        axis([0,length(DEER)+1,0,1.05*ma]);
    end;
    xlabel('Constraint number');
    ylabel('Distance (nm)');
    full_rmsd=full_rmsd+rmsd;
end;
if fid~=-1,
    full_rmsd=sqrt(full_rmsd/full_num);
    fprintf(fid,'\nMean r.m.s.d. of all restraints for all models: %5.2f nm\n',full_rmsd);
    fclose(fid);
    hMain.report_file=fname;
    curr_fig=gcf;
    report_editor;
    figure(curr_fig);
end;

% --- Executes on button press in pushbutton_import.
function pushbutton_import_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global model

pname=general.tmp_files;
idCode=handles.targ_ID;
logname=get(handles.pushbutton_import,'UserData');
ensemble=assess_ensemble(logname,'',handles.ensemble_size,handles.min_GA341);
if isempty(ensemble),
    msg='ERROR: No ensemble model could be generated.';
    set(handles.text_info,'String',msg);
    add_msg_board(msg);
    return
end;
ensemble_size=length(ensemble);
if ensemble_size < handles.ensemble_size,
    add_msg_board(sprintf('Warning: only %i models have sufficient GA341 score, while %i models were requested.',ensemble_size,handles.ensemble_size));
    add_msg_board('Importing smaller ensemble. Consider increasing the number of models.');
end;
[message,snum,stag,models]=rd_pdb_ensemble(ensemble,idCode,pname);
ares=[];
if length(handles.alignment)>1,
    ali=handles.alignment(handles.target_info.alg);
    if isfield(ali,'first_res') && str2double(ali.first_res)>0,
        ares=str2double(ali.first_res);
    end;
elseif isfield(handles,'ares'),
    ares=handles.ares;
end;
if isfield(handles.restraints,'chains'),
    chain_tags=handles.restraints.chains;
else
    chain_tags='';
end;
correct_Modeller_ensemble(snum,ares,chain_tags);
if message.error,
    add_msg_board(message.text);
else
    add_msg_board(sprintf('Ensemble read into structure %i with tag [%s].',snum,stag));
    add_msg_board(sprintf('This ensemble has %i chain models',models));
end;
DEER=handles.DEER;
for k=1:length(DEER),
    [targadr1,cid,ctag1,resnum1]=target_address_to_modeller(DEER(k).adr1,handles.translation);
    [targadr2,cid,ctag2,resnum2]=target_address_to_modeller(DEER(k).adr2,handles.translation);
    targadr1=sprintf('[%s]{1}%i',stag,resnum1);
    targadr2=sprintf('[%s]{1}%i',stag,resnum2);
%     ind1=resolve_address(DEER(k).adr1);
%     ind2=resolve_address(DEER(k).adr2);
    ind1=resolve_address(targadr1);
    ind2=resolve_address(targadr2);
    DEER(k).indices=[ind1;ind2];
end;
handles.DEER=DEER;
[handles,all_distance_rmsd]=update_3D_model(handles,snum,DEER);
id=[stag 'r'];
exists=tag2id(id,model.structure_tags);
while ~isempty(exists),
    id=[id 'a'];
    exists=tag2id(id,model.structure_tags);
end;
[snum1,message]=copy_structure(snum,id);
add_msg_board('Now removing spin labels...');
[repnum,msg]=replace(snum1,':R1A:IA1:',false,false);
add_msg_board(sprintf('%i spin labels were removed.',repnum));
add_msg_board('Now repacking side chains...');
for modnum=1:models,
    message=repack(snum1,modnum,0,handles.clean_targ_seq);
end;
msg='Modeller ensemble was succesfully imported.';
set(handles.text_info,'String',msg);
add_msg_board(msg);
if handles.test_mode,
    ensemble_address=mk_address(snum1);
    model_address=sprintf('%s{1}',ensemble_address);
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(model_address);
    model.coarse(snum1).Ca_coor=Ca_coor;
    model.coarse(snum1).indices=rindices;
    model.coarse(snum1).masses=masses;
    model.coarse(snum1).Bfactors=Bfactors;
    model.coarse(snum1).restypes=restypes;
    [coarse_correspondence,coarse_target]=align_template_target(mk_address(snum1),handles.target);
    add_msg_board(sprintf('Generating test information for ensemble %s with respect to target %s...',ensemble_address,handles.target));
    rms0=rmsd_superimpose(handles.coarse_target(handles.coarse_correspondence(2,:),:),handles.coarse_template(handles.coarse_correspondence(1,:),:));
    logname=sprintf('%s_Modeller_ensemble_test_%i_constraints.log',handles.target,length(handles.DEER));
    fid=fopen(logname,'w');
    fprintf(fid,'Ensemble logfile for target %s\n\n',handles.target);
    all_coverage=zeros(1,models);
    all_rmsd=zeros(1,models);
    for k=1:models,
        model_address=sprintf('%s{%i}',ensemble_address,k);
        [network]=coarse_residues(model_address);
        rms=rmsd_superimpose(coarse_target(coarse_correspondence(2,:),:),network(coarse_correspondence(1,:),:));
        all_rmsd(k)=rms;
        coverage=(rms0-rms)/rms0;
        fprintf(fid,'\n*** Trial %i ***\n\n',k);
        add_msg_board(sprintf('Trial %i: Initial r.m.s.d. w.r.t. target structure: %5.2f ?',k,rms0));
        add_msg_board(sprintf('Trial %i: Final   r.m.s.d. w.r.t. target structure: %5.2f ?',k,rms));
        add_msg_board(sprintf('Trial %i: Coverage of the structural change: %5.3f\n',k,coverage));
        fprintf(fid,'\nFinal   r.m.s.d. w.r.t. target structure: %5.2f ?\n',rms);
        fprintf(fid,'Distance constraint r.m.s.d.     : %5.2f\n',all_distance_rmsd(k));
        fprintf(fid,'Coverage of the structural change: %5.3f\n',coverage);
        all_coverage(k)=coverage;
    end;
    fprintf(fid,'\n\nMean final   r.m.s.d. w.r.t. target structure  : %5.2f ?\n',mean(all_rmsd));
    fprintf(fid,'Std. dev. of r.m.s.d. w.r.t. target structure  : %5.2f ?\n',std(all_rmsd));
    fprintf(fid,'Mean distance constraint r.m.s.d.              : %5.2f ?\n',mean(all_distance_rmsd));
    fprintf(fid,'Std. dev. of distance constraint r.m.s.d.      : %5.2f ?\n',std(all_distance_rmsd));
    fprintf(fid,'Mean coverage of the structural change         : %5.3f\n',mean(all_coverage));
    fprintf(fid,'Std. dev. of coverage of the structural change : %5.3f\n',std(all_coverage));
    fprintf(fid,'Worst        coverage of the structural change : %5.3f\n',min(all_coverage));
    fprintf(fid,'Best         coverage of the structural change : %5.3f\n',max(all_coverage));
    diff=all_coverage-mean(all_coverage);
    [~,poi]=min(abs(diff));
    fprintf(fid,'The most typical result was obtained in trial %i.\n',poi);
    fprintf(fid,'Deviation of this result from mean coverage is: %5.3f\n\n',diff(poi));
    fclose(fid);
    add_msg_board('Test information generated. Import is complete.');
end;
set(gcf,'Pointer','arrow');
set(handles.figure1,'Pointer','arrow');
guidata(hObject,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain
global model

hMain.fit_plot=false;
model.current_structure=handles.current_structure;
update_current_structure;
delete(hObject);


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

[c,n]=size(handles.chains);

figure; clf; hold on

snum=model.current_structure;

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(mk_address(snum));
    model.coarse(snum).Ca_coor=Ca_coor;
    model.coarse(snum).indices=rindices;
    model.coarse(snum).masses=masses;
    model.coarse(snum).Bfactors=Bfactors;
    model.coarse(snum).restypes=restypes;
end;

Ca_coor=model.coarse(snum).Ca_coor;

wires=zeros(1,c);

for k=1:c,
    col=handles.darken*color_grade(handles.chains(k,1),c);
    x=Ca_coor(handles.chains(k,2):handles.chains(k,3),1);
    y=Ca_coor(handles.chains(k,2):handles.chains(k,3),2);
    z=Ca_coor(handles.chains(k,2):handles.chains(k,3),3);
    wires(k)=line(x,y,z,'color',col,'LineWidth',2);
end;

if isfield(handles,'DEER'),
    DEER=handles.DEER;
    for k=1:c,
        set(wires(k),'Color',[0.25,0.25,0.25]);
    end;
    for k=1:length(DEER),
        lindices=DEER(k).indices;
        res1adr=mk_address([snum lindices(1,2) 1 lindices(1,4)]);
        [msg,N1coor]=get_object([res1adr '.N1'],'coor');
        [msg,O1coor]=get_object([res1adr '.O1'],'coor');
        DEER(k).xyz1=(N1coor+O1coor)/2;
        res2adr=mk_address([snum lindices(2,2) 1 lindices(2,4)]);
        [msg,N1coor]=get_object([res2adr '.N1'],'coor');
        [msg,O1coor]=get_object([res2adr '.O1'],'coor');
        DEER(k).xyz2=(N1coor+O1coor)/2;
        line(DEER(k).xyz1(1),DEER(k).xyz1(2),DEER(k).xyz1(3),...
            'Marker','.','Color','b');
        line(DEER(k).xyz2(1),DEER(k).xyz2(2),DEER(k).xyz2(3),...
            'Marker','.','Color','b');
        x=[DEER(k).xyz1(1) DEER(k).xyz2(1)];
        y=[DEER(k).xyz1(2) DEER(k).xyz2(2)];
        z=[DEER(k).xyz1(3) DEER(k).xyz2(3)];
        r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
        det=abs(r0-DEER(k).r)/DEER(k).sigr;
        if det>2,
            col='r';
        elseif det>1,
            col=[255,190,0]/255; % yellow-orange
        else
            col='g';
        end;
        line(x,y,z,'Color',col,'LineWidth',1.5,'LineStyle',':');
        cindices=DEER(k).indices;
        cindices(1,1)=snum;
        cindices(2,1)=snum;
        cindices(1,3)=1;
        cindices(2,3)=1;
        f1=false;
        f2=false;
        for l=1:length(model.coarse(snum).indices),
            diff=cindices(1,:)-model.coarse(snum).indices(l,:);
            if sum(abs(diff))==0,
                DEER(k).res1=l;
                f1=true;
                x=[DEER(k).xyz1(1) Ca_coor(l,1)];
                y=[DEER(k).xyz1(2) Ca_coor(l,2)];
                z=[DEER(k).xyz1(3) Ca_coor(l,3)];
                line(x,y,z,'Color','b','LineWidth',1);
            end;
            diff=cindices(2,:)-model.coarse(snum).indices(l,:);
            if sum(abs(diff))==0,
                DEER(k).res2=l;
                f2=true;
                x=[DEER(k).xyz2(1) Ca_coor(l,1)];
                y=[DEER(k).xyz2(2) Ca_coor(l,2)];
                z=[DEER(k).xyz2(3) Ca_coor(l,3)];
                DEER(k).rl2=line(x,y,z,'Color','b','LineWidth',1);
            end;
        end;
    end;
end;

axis equal
axis off
hold on;

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
cam_up=get(hMain.axes_model,'CameraUpVector');
set(gca,'CameraPosition',cam_pos);
set(gca,'CameraTarget',cam_tar);
set(gca,'CameraUpVector',cam_up);
camlookat(gca);


% --- Executes on button press in pushbutton_copy_2.
function pushbutton_copy_2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

snum=model.current_structure;
modnum=length(model.structures{snum}(1).xyz);

if isfield(handles,'DEER'),
    DEER=handles.DEER;
    figure; clf;
    hold on;
    set(gca,'FontSize',16);
    for km=1:modnum,
        ma=0;
        rmsd=0;
        fom=0; % figure of merit
        if ~isempty(DEER),
            for k=1:length(DEER),
                lindices=DEER(k).indices;
                res1adr=mk_address([snum lindices(1,2) km lindices(1,4)]);
                [msg,N1coor]=get_object([res1adr '.N1'],'coor');
                [msg,O1coor]=get_object([res1adr '.O1'],'coor');
                DEER(k).xyz1=(N1coor+O1coor)/2;
                res2adr=mk_address([snum lindices(2,2) km lindices(2,4)]);
                [msg,N1coor]=get_object([res2adr '.N1'],'coor');
                [msg,O1coor]=get_object([res2adr '.O1'],'coor');
                DEER(k).xyz2=(N1coor+O1coor)/2;
                r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
                rmsd=rmsd+(r0-DEER(k).r)^2;
                det=abs(r0-DEER(k).r)/DEER(k).sigr;
                fom=fom+det^2;
                if det>2,
                    col='r';
                elseif det>1,
                    col=[255,190,0]/255; % yellow-orange
                else
                    col='g';
                end;
                errorbar(k,DEER(k).r,DEER(k).sigr,'k');
                if DEER(k).r+DEER(k).sigr>ma,
                    ma=DEER(k).r+DEER(k).sigr;
                end;
                if r0>ma,
                    ma=r0;
                end;
                hold on;
                plot(k,r0,'.','Color',col);
            end;
            axis([0,length(DEER)+1,0,1.05*ma]);
            xlabel('Constraint number');
            ylabel('Distance (nm)');
        end;
    end;
end;



function edit_models_Callback(hObject, eventdata, handles)
% hObject    handle to edit_models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_models as text
%        str2double(get(hObject,'String')) returns contents of edit_models as a double

[v,handles]=edit_update_MMM(handles,hObject,1,500,200,'%i',1);
handles.models=v;
if v<handles.ensemble_size,
    handles.ensemble_size=v;
    set(handles.edit_ensemble_size,'String',sprintf('%i',v));
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_models_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ensemble_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ensemble_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ensemble_size as text
%        str2double(get(hObject,'String')) returns contents of edit_ensemble_size as a double

[v,handles]=edit_update_MMM(handles,hObject,1,100,20,'%i',1);
if v>handles.models,
    handles.models=v;
    set(handles.edit_models,'String',sprintf('%i',v));
end;
handles.ensemble_size=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_ensemble_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ensemble_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_template_ID_Callback(hObject, eventdata, handles)
% hObject    handle to edit_template_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_template_ID as text
%        str2double(get(hObject,'String')) returns contents of edit_template_ID as a double

tempID=get(hObject,'String');
if isempty(tempID), 
    tempID='temp';
    set(hObject,'String',tempID);
end;
handles.temp_ID=tempID;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_template_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_template_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_target_ID_Callback(hObject, eventdata, handles)
% hObject    handle to edit_target_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_target_ID as text
%        str2double(get(hObject,'String')) returns contents of edit_target_ID as a double

targID=get(hObject,'String');
if isempty(targID), 
    targID='targ'; 
    set(hObject,'String',targID);
end;
handles.targ_ID=targID;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_target_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_target_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function message=mk_modeller_input(runfile,alignfile,temp_ID,targ_ID,restraints,translation,models,weighting)

tightness=0.04; % tightness of DEER distance restraint fitting, 0.040 is optimized
tightness_Ca=0.04; % tightness of Ca distance restraint fitting, 0.040 is optimized

% redefine relative weighting of distance constraints from optional client
% input
if exist('weighting','var'),
    tightness=weighting.DEER;
    tightness_Ca=weighting.Ca;
end;

%tightness=tightness*sqrt(2);

message.error=0;
message.text='No error.';

single_chain=true; % currently, MMM supports only single-chain homology modeling

cali=get_multiple_pir(alignfile);

if isfield(restraints,'DEER'),
    DEER=restraints.DEER;
else
    DEER=[];
end;
if isfield(restraints,'direct'),
    direct=restraints.direct;
else
    direct=[];
end;
if isfield(restraints,'helices'),
    helices=restraints.helices;
else
    helices=[];
end;
if isfield(restraints,'strands'),
    strands=restraints.strands;
else
    strands=[];
end;
if isfield(restraints,'sheets'),
    sheets=restraints.sheets;
else
    sheets=[];
end;
if ~isempty(translation),
    chains=length(translation);
    offsets=zeros(1,chains);
    new_targ_tags=':';
    for k=1:chains,
        new_targ_tags=[new_targ_tags char(double('A')-1+k) ':'];
        offsets(k)=translation(k).targ_offset;
    end;
    chain_tags=':';
    if isfield(translation,'chain'),
        for k=1:length(translation),
            chain_tags=[chain_tags translation(k).chain ':'];
        end;
    else
        offsets=zeros(1,100);
        chain_tags=':A:';
    end;
else
    offsets=zeros(1,100);
    chain_tags=':A:';
end;


% fid_trans=fopen('PutP_to_vSGLT_Olkhova_alignment.txt','wt');
% fprintf(fid_trans,'PutP res1  PutP res2  vSGLT res1  vSGLT res2\n');


fid=fopen(runfile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller input file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'from modeller import *\n');
fprintf(fid,'from modeller.automodel import *\n\n');
fprintf(fid,'log.verbose()\n');
fprintf(fid,'env = environ(rand_seed=-8123, restyp_lib_file=''$(LIB)/restyp_EPR.lib'', copy=None)\n\n');

fprintf(fid,'env.io.atom_files_directory = [''.'', ''../atomfiles'']\n\n');

fprintf(fid,'env.io.hetatm = True\n\n');

fprintf(fid,'class MyModel(automodel):\n');
fprintf(fid,'    def special_restraints(self, aln):\n');
fprintf(fid,'        rsr = self.restraints\n');
fprintf(fid,'        at  = self.atoms\n');

if ~isempty(DEER),
    % make list of labeled residues
    reslist=zeros(1,400);
    pairlist=zeros(length(DEER),2);
    poi=0;
    for k=1:length(DEER),
        [targadr,cid,ctag1,resnum1]=target_address_to_modeller(DEER(k).adr1,translation);
        % [stag,ctag,modelnum,resnum1]=mk_address_parts(DEER(k).indices(1,:));
        [targadr,cid,ctag2,resnum2]=target_address_to_modeller(DEER(k).adr2,translation);
        
        % fprintf(fid_trans,'   %4s       %4s        %4i        %4i\n',DEER(k).adr1,DEER(k).adr2,resnum1,resnum2);
        % [stag,ctag,modelnum,resnum2]=mk_address_parts(DEER(k).indices(2,:));
        resnum1=resnum1+offsets;
        resnum2=resnum2+offsets;
        pairlist(k,:)=[resnum1 resnum2];
        if isempty(find(reslist(1:poi)==resnum1,1)),
            poi=poi+1;
            reslist(poi)=resnum1;
            chain_assign{poi}=ctag1;
        end;
        if isempty(find(reslist(1:poi)==resnum2,1)),
            poi=poi+1;
            reslist(poi)=resnum2;
            chain_assign{poi}=ctag2;
        end;
    end;
    reslist=sort(reslist(1:poi));
    for k=1:length(reslist);
        if isempty(chain_assign{k}) || single_chain,
            fprintf(fid,'        NO%i = at[''N1:%i''], at[''O1:%i'']\n',reslist(k),reslist(k),reslist(k));
        else
            fprintf(fid,'        NO%i = at[''N1:%i:%s''], at[''O1:%i:%s'']\n',reslist(k),reslist(k),chain_assign{k},reslist(k),chain_assign{k});
        end;
        fprintf(fid,'        label%i = pseudo_atom.gravity_center(NO%i)\n',reslist(k),reslist(k));
        fprintf(fid,'        rsr.pseudo_atoms.append(label%i)\n',reslist(k));
    end;
    fprintf(fid,'\n');
    
    for k=1:length(reslist);
        if isempty(chain_assign{k}) || single_chain,
            fprintf(fid,'	r%i = rigid_body(self.residue_range(''%i'', ''%i''))\n',reslist(k),reslist(k),reslist(k));
        else
            fprintf(fid,'	r%i = rigid_body(self.residue_range(''%i:%s'', ''%i:%s''))\n',reslist(k),reslist(k),chain_assign{k},reslist(k),chain_assign{k});
        end;
        fprintf(fid,'	rsr.rigid_bodies.append(r%i)\n',reslist(k));
    end;
    fprintf(fid,'\n');

    for k=1:length(DEER);
        fprintf(fid,'	rsr.add(forms.gaussian(group=physical.xy_distance,\n');
        fprintf(fid,'                               feature=features.distance(label%i,\n',pairlist(k,1));
        fprintf(fid,'                                                         label%i),\n',pairlist(k,2));
        fprintf(fid,'                               mean=%6.3f, stdev=%5.3f))\n',10*DEER(k).r,tightness);
    end;
    fprintf(fid,'\n');
end;

% fclose(fid_trans);
if ~isempty(helices),
    for k=1:length(helices);
        [cid,resnum1,resnum2]=dissect_address(helices(k).adr,chain_tags);
        resnum1=resnum1+offsets;
        resnum2=resnum2+offsets;
        if ~isempty(cid),
            adr1=sprintf('(%s)%i',id2tag(cid,chain_tags),resnum1);
            adr2=sprintf('(%s)%i',id2tag(cid,chain_tags),resnum2);
        else
            adr1=sprintf('%i',resnum1);
            adr2=sprintf('%i',resnum2);
        end;
        [targadr,cid,ctag1,resnum1]=target_address_to_modeller(adr1,translation);
        [targadr,cid,ctag2,resnum2]=target_address_to_modeller(adr2,translation);
        fprintf(fid,'	rsr.add(secondary_structure.alpha(self.residue_range(');
        if isempty(ctag1) || isempty(ctag2) || single_chain,
            fprintf(fid,'''%i:'',''%i:'')))\n',resnum1,resnum2);
        else
            fprintf(fid,'''%i:%s'',''%i:%s'')))\n',resnum1,ctag1,resnum2,ctag2);
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(strands),
    for k=1:length(strands);
        [cid,resnum1,resnum2]=dissect_address(strands(k).adr,chain_tags);
        resnum1=resnum1+offsets;
        resnum2=resnum2+offsets;
        if ~isempty(cid),
            adr1=sprintf('(%s)%i',id2tag(cid,chain_tags),resnum1);
            adr2=sprintf('(%s)%i',id2tag(cid,chain_tags),resnum2);
        else
            adr1=sprintf('%i',resnum1);
            adr2=sprintf('%i',resnum2);
        end;
        [targadr,cid,ctag1,resnum1]=target_address_to_modeller(adr1,translation);
        [targadr,cid,ctag2,resnum2]=target_address_to_modeller(adr2,translation);
        fprintf(fid,'	rsr.add(secondary_structure.strand(self.residue_range(');
        if isempty(ctag1) || isempty(ctag2) || single_chain,
            fprintf(fid,'''%i:'',''%i:'')))\n',resnum1,resnum2);
        else
            fprintf(fid,'''%i:%s'',''%i:%s'')))\n',resnum1,ctag1,resnum2,ctag2);
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(sheets),
    for k=1:length(sheets);
        [targadr,cid,ctag1,resnum1]=target_address_to_modeller(sheets(k).adr1,translation);
        [targadr,cid,ctag2,resnum2]=target_address_to_modeller(sheets(k).adr2,translation);
        resnum1=resnum1+offsets;
        resnum2=resnum2+offsets;
        fprintf(fid,'	rsr.add(secondary_structure.sheet(at[,\n');
        if isempty(ctag1) || isempty(ctag2) || single_chain,
            fprintf(fid,'''N:%i''], at[''O:%i''],\n',resnum1,resnum2);
            fprintf(fid,'                               sheet_h_bonds=%i))\n',str2double(sheets(k).length));
        else
            fprintf(fid,'''N:%i:%s''], at[''O:%i:%s''],\n',resnum1,ctag1,resnum2,ctag2);
            fprintf(fid,'                               sheet_h_bonds=%i))\n',str2double(sheets(k).length));
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(direct),
    for k=1:length(direct);
        [targadr,cid,ctag1,resnum1]=target_address_to_modeller(direct(k).adr1,translation);
        [targadr,cid,ctag2,resnum2]=target_address_to_modeller(direct(k).adr2,translation);
        resnum1=resnum1+offsets;
        resnum2=resnum2+offsets;
        fprintf(fid,'	rsr.add(forms.gaussian(group=physical.xy_distance,\n');
        if isempty(ctag1) || isempty(ctag2) || single_chain,
            fprintf(fid,'                               feature=features.distance(at[''CA:%i''],\n',resnum1);
            fprintf(fid,'                                                         at[''CA:%i'']),\n',resnum2);
        else
            fprintf(fid,'                               feature=features.distance(at[''CA:%i:%s''],\n',resnum1,ctag1);
            fprintf(fid,'                                                         at[''CA:%i:%s'']),\n',resnum2,ctag2);
        end;
        fprintf(fid,'                               mean=%6.3f, stdev=%5.3f))\n',10*direct(k).r,tightness_Ca);
    end;
    fprintf(fid,'\n');
end;

fprintf(fid,'a = MyModel(env, alnfile = ''%s'',\n',alignfile);
fprintf(fid,'            knowns = ''%s'', sequence = ''%s'',\n',temp_ID,targ_ID);
fprintf(fid,'            assess_methods=(assess.DOPE, assess.GA341))\n');

fprintf(fid,'a.starting_model= 1\n');
fprintf(fid,'a.ending_model  = %i\n\n',models);

fprintf(fid,'a.make()\n');

fclose(fid);

function [batname,message]=mk_modeller_bat(modpath,modcmd,jobfile,is_bckg)

global general

if nargin<4,
    is_bckg=false;
end;

batname='run_modeller.bat';
envfile= which('modenv.bat');
batfile=[modpath filesep  batname];

message.error=0;
message.text='No error.';

fid=fopen(batfile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller batch file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

rfid=fopen(envfile,'rt');
if rfid==-1,
    message.error=3;
    message.text='Modeller environment file could not be read.';
    add_msg_board('ERROR: Modeller environment file could not be read.');
    add_msg_board(message.text);
    return;
end;

while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    fprintf(fid,'%s\n',tline);
end;

fprintf(fid,'echo modeller installed\n');

k=strfind(general.tmp_files,':');
if ~isempty(k),
    fprintf(fid,'%s\n',general.tmp_files(1:k));
end;

fprintf(fid,'cd %s\n',general.tmp_files);
if is_bckg,
    fprintf(fid,'start %s %s\n',modcmd,jobfile);
else
    fprintf(fid,'%s %s\n',modcmd,jobfile);
end;

fclose(fid);


% --- Executes on button press in pushbutton_alignment.
function pushbutton_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general
global hMain

my_path=pwd;
cd(general.restraint_files);

[fname,pname,findex]=uigetfile(...
    {'*.ali;*.pir', 'PIR files (*.ali, *.pir)'; ...
     '*.afa;*.ebi;*.fa;*.fasta;*.txt','FASTA files (*.afa, *.ebi, *.fa, *.fasta, *.txt)'; ...
     '*.aln;*.clw;*.clustalw','Clustal W files (*.aln, *.clw, *.clustalw)'; ...
     '*.*', 'All files'}, ...
     'Load alignment file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Loading of alignment cancelled by user');
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [pathstr,name,ext]=fileparts(fname);
    switch ext
        case {'.afa','.fa','.fasta','.ebi','.txt'}
            alignment=get_multiple_fasta(fullfile(pname,fname));
        case {'.ali','.pir'}
            alignment=get_multiple_pir(fullfile(pname,fname));
        case {'.clw','.aln','.clustalw'}
            alignment=get_multiple_clustal(fullfile(pname,fname));
        otherwise
            add_msg_board('Warning: Unknown alignment file format. Assuming FASTA.');
            alignment=get_multiple_fasta(fullfile(pname,fname));
    end;
    alignment=test_repair_alignment(alignment);
    set(gcf,'Pointer','watch');
    drawnow;
    if ~isempty(alignment),
        handles.alignment=alignment;
        [best_match,name,mode,identity]=match_alignment(alignment);
        if identity<0.95,
            add_msg_board('Template structure does not seem to exist in alignment file.');
            add_msg_board('Adding current chain of template to alignment array and realigning by MUSCLE.');
            n=length(alignment);
            alignment(n+1).db='MMM';
            alignment(n+1).name=[model.info{model.current_structure}.idCode model.current_chain];
            cnum=tag2id(model.current_chain,model.chain_tags{model.current_structure});
            alignment(n+1).sequence=model.structures{model.current_structure}(cnum).sequence;
            alignment0=alignment;
            alignment=test_repair_alignment(alignment);
            if ~isempty(alignment),
                [best_match,name,mode,identity]=match_alignment(alignment);
            else
                alignment=alignment0;
            end;
        end;
        handles.template_info.alg=best_match;
        handles.template_info.name=name;
        handles.template_info.mode=mode;
        handles.template_info.identity=identity;
        add_msg_board(sprintf('Sequence %i of %s in the alignment array best matches template structure',best_match,name));
        if mode==0,
            add_msg_board('This sequence matches all chains of the template');
        else
            ctag=id2tag(mode,char(model.chain_tags{model.current_structure}));
            add_msg_board(sprintf('This sequence matches chain %s of the template',ctag));
        end;
        add_msg_board(sprintf('Sequence identity of match with template is %4.1f%%',100*identity));
        if length(alignment)==1,
            add_msg_board('Warning: Only one sequence found. Template and target sequences are assumed to be identical.');
            handles.target_info.alg=1;
            handles.target_info.identity=1;
        elseif length(alignment)==2,
            if best_match==1,
                handles.target_info.alg=2;
            else
                handles.target_info.alg=1;
            end;
            ali(1)=alignment(handles.target_info.alg);
            ali(2)=alignment(handles.template_info.alg);
            handles.target_info.identity=get_identity(ali);
            msg=sprintf('Selected target sequence %s with %4.1f%% identity to template.',...
                alignment(handles.target_info.alg).name,...
                100*handles.target_info.identity);
            add_msg_board(msg);
        else
            if isfield(hMain,'container'),
                hMain=rmfield(hMain,'container');
            end;
            hMain.container.alignment=alignment;
            hMain.container.template_id=best_match;
            select_target;
            if isempty(hMain.container.alg) || isempty(hMain.container.identity),
                add_msg_board('Warning: Target sequence selection canceled.');
                return
            end;
            handles.target_info.alg=hMain.container.alg;
            handles.target_info.identity=hMain.container.identity;
            msg=sprintf('Selected target sequence %s with %4.1f%% identity to template.',...
                alignment(handles.target_info.alg).name,...
                100*handles.target_info.identity);
            add_msg_board(msg);
        end;
        ident=handles.target_info.identity;
        if ident>0.4,
            col=[0,1,0];
        elseif ident<0.2,
            col=[1,0,0];
        else
            twilight=(ident-0.2)/0.2;
            col=twilight*[0,1,0]+(1-twilight)*[1,0,0];
        end;
        set(handles.pushbutton_alignment,'ForegroundColor',col);
    else
        add_msg_board('ERROR: Alignment file is invalid. No sequence found.');
        return
    end;
    set(gcf,'Pointer','arrow');
    drawnow;
end;

if handles.target_info.identity>0.15,
    if length(alignment)>1,
        set(handles.pushbutton_run,'Enable','on');
        if handles.target_info.identity>0.4,
            set(handles.text_info,'String','Valid alignment for homology modeling.');
        else
            set(handles.text_info,'String','Twilight zone alignment for homology modeling.');
        end;
    else
        add_msg_board('Run button remains disabled.');
        set(handles.pushbutton_alignment,'ForegroundColor',[0,0,0]);
        set(handles.text_info,'String','Template and target sequence identical. No homology modeling.');
    end;        
else
    set(handles.pushbutton_run,'Enable','off');
    add_msg_board('Warning: Sequence identity is in midnight zone.');
    add_msg_board('Run button remains disabled.');
    set(handles.text_info,'String','Midnight zone alignment. Homology modeling disabled.');
end;
cd(my_path);
guidata(hObject,handles);

function [best_match,name,mode,identity]=match_alignment(alignment,snum)
% function [best_match,name,mode,identity]=match_alignment(alignment,snum)
%
% finds the best matching record in an alignment array for the current
% structure or structure with number snum
% either a single chain or all chains can be matched, for all-chain
% matches, the aligned sequences must have '/' characters for chain breaks
%
% alignment     MMM alignment array, each alignemnt must have the field
%               .sequence and the field .name
% snum          (optional) structure number, if missing, the current
%               structure is used
%
% best_match    index into alignment for the best matching record
% name          name of the best matchin sequence
% mode          is 0 for multiple-chain, else the MMM index of the matching
%               chain
% identity      sequence identity of the best match
%
% Remarks:
% - identity should ideally be 1.0, but PDB files do not necessarily match
%   the database sequences to which they refer
% - note that alignment arrays are a rather flexible data structure, in
%   which sequences may not yet be aligned (and may thus have different
%   length)
% - best_match is empty if the alignment array is empty, identity is 0 in
%   that case, as is mode
% - matching of multiple chains requires that the number of chains in the
%   alignment record and in the query structure agree
% - multiple-chain matches are always preferred to single-chain matches, if
%   they exist; don't provide multiple-chain alignment records if you want
%   to force single-chain matching
% - multiple-chain alignment records that do not match the chain number in
%   the query structure are ignored
%
% G. Jeschke, 2011

global model

if nargin<2,
    snum=model.current_structure;
end;

chains=length(model.structures{snum});

best_match=[];
mode=0;
identity=0;
multiple=false;

if isempty(alignment),
    return;
end;

for k=1:length(alignment),
    seq=alignment(k).sequence;
    poi=0;
    for kk=1:length(seq), % remove gaps
        if char(seq(kk))~='-',
            poi=poi+1;
            seq(poi)=seq(kk);
        end;
    end;
    seq=seq(1:poi);
    nonsense=textscan(seq,'%s','Delimiter','/'); % separate PIR format into single chain sequences
    seqs=nonsense{1};
    if length(seqs)>1 && length(seqs)==chains, % number of chains in alignment and query structure matches
        cmode=0;
        if ~multiple,
            identity=0; % cancel all previous single-chain matches
        end;
        multiple=true;
    elseif length(seqs)>1,
        continue;
    end;
    if ~multiple, % only if no multiple-chain match was ever done
        cseqs{1}=char(seqs(1));
        bident=0;
        for kk=1:chains,
            cseqs{2}=model.structures{snum}(kk).sequence;
            [message,inname]=align_sequences(cseqs,[],true);
            if message.error,
                continue;
            else
                ali=get_multiple_clustal(inname);
                cident=get_identity(ali);
                if cident>bident,
                    bident=cident;
                    cmode=kk;
                end;
            end;
        end;
        if bident>identity,
            identity=bident;
            mode=cmode;
            best_match=k;
            name=alignment(k).name;
        end;
    else % multiple-chain match
        bident=0;
        totlen=0;
        for kk=1:chains,
            cseqs{1}=char(seqs(kk));
            cseqs{2}=model.structures{snum}(kk).sequence;
            [message,inname]=align_sequences(cseqs,[],true);
            if message.error,
                continue;
            else
                ali=get_multiple_clustal(inname);
                [cident,len]=get_identity(ali);
                totlen=totlen+len;
                bident=bident+len*cident;
            end;
        end;
        bident=bident/totlen;
        if bident>identity,
            identity=bident;
            mode=cmode;
            best_match=k;
            name=alignment(k).name;
        end;
    end;
end;

function [ident,len]=get_identity(ali)
% in an alignment array with exactly two sequences, get_identity
% determines relative sequence identity between the sequences normalized to
% the length of the first sequence
% for an alignment of more than two sequences, identity of zero is returned
% identity of zero is also returned if sequence lengths do not match (not
% aligned) or sequences are empty
% the length len of the first sequence is also returned

ident=0;
if length(ali)~=2,
    return
end;
seq1=ali(1).sequence;
seq2=ali(2).sequence;
if length(seq1)~=length(seq2),
    return
end;
len=0;
matches=0;
for k=1:length(seq1),
    if char(seq1(k))~='-', % no gap in sequence 1
        len=len+1; 
        if char(seq1(k))==char(seq2(k)),
            matches=matches+1;
        end;
    end;
end;
if len>0,
    ident=matches/len;
end;

function alignment=test_repair_alignment(alignment,force)
% tests whether all sequences in an alignment array have the same length,
% if not, the alignment is repaired by calling MUSCLE
% if alignment by MUSCLE fails, an empty alignment record is returned
%
% force     flag that forces realignment even if sequences have the same
%           length

global model
global third_party
global help_files
global general

if nargin<2,
    force=false;
end;

if length(alignment)<2,
    return
end;

aligned=false;
if ~force,
    len=length(alignment(1).sequence),
    aligned=true;
    for k=2:length(alignment),
        len2 = length(alignment(k).sequence),
        if length(alignment(k).sequence)~=len,
            aligned=false;
        end;
    end;
end;

if ~aligned,
    add_msg_board('Warning: Sequence length mismatch in imported alignment.');
    add_msg_board('Alignment by MUSCLE will be attempted.');
    dospath=which('muscle.exe');
    entry=strcat(help_files,'third_party.html#MUSCLE');
    if isempty(dospath),
        add_msg_board('ERROR: MUSCLE could not be found on the Matlab path.');
        add_msg_board('Please check whether MUSCLE is installed and the path set.');
        add_msg_board('(see also help browser)');
        webcall(entry,'-helpbrowser');
        alignment=[];
        return
    end;
    % add the reference, if it does not yet exist
    muscle_ref=true;
    id=tag2id('Edgar:2004_muscle',third_party.tags,[],'|');
    if ~isempty(id),
        if isfield(model,'auto_references'),
            if ~isempty(find(id==model.auto_references, 1)),
                muscle_ref=false;
            end;
        else
            model.auto_references=[];
        end;
        if muscle_ref,
            if ~isfield(model,'references'),
                model.references(1)=third_party.references(id);
            elseif isempty(model.references)
                model=rmfield(model,'references');
                model.references(1)=third_party.references(id);
            else
                model.references(end+1)=third_party.references(id);
            end;
            model.auto_references(end+1)=id;
        end;
    end;
    infile=fullfile(general.tmp_files,'doalign.fa');
    wr_multiple_fasta(infile,alignment);
    outfile=[general.tmp_files 'isaligned.afa'];
    comd=[dospath ' -in ' infile ' -out ' outfile];
    [s,w] = dos(comd);
    if s~=0,
        message.error=3;
        message.text='MUSCLE error.';
        add_msg_board('ERROR: MUSCLE did not run successfully.');
        alignment=[];
        return
    else
        add_msg_board('Now importing aligned sequences.');
        alignment=get_multiple_fasta(outfile);
    end;
end;

function [temp_seq,targ_seq,darestarg,clean_targ_seq]=adjust_alignment(ali,temp_seq)
% adjusts a given alignment to the actual template sequence as determined
% by wr_pdb_paradigm, 
% this takes care of mutations present in the PDB file, but not the
% alignment file
% spin labels in the template are also defined in the target
% darestarg is the offset of the first residue of the remaining target
% sequence with respect to the original target sequence

global Modeller_info

if isfield(ali,'first_res')
    if isfield(ali,'last_res'),
        templen=str2double(ali(1).last_res)-str2double(ali(1).first_res)+1;
        temp_seq=temp_seq(1:templen);
    end;
end;


label_codes=Modeller_info.label_codes;
het_codes=Modeller_info.het_codes;

labels=0;
cofactors=0;
blocks=0;
gaps=0;

% clean up Modeller format sequence for MUSCLE alignment
temp_seq_rep=temp_seq;
for k=1:length(temp_seq),
    if char(temp_seq(k))=='-',
        temp_seq_rep(k)='X';
        gaps=gaps+1;
    end;
    if char(temp_seq(k))=='.',
        temp_seq_rep(k)='X';
        blocks=blocks+1;
    end;
    if ~isempty(strfind(label_codes,temp_seq(k))),
        temp_seq_rep(k)='X';
        labels=labels+1;
    end;        
    if ~isempty(strfind(het_codes,temp_seq(k))),
        temp_seq_rep(k)='X';
        cofactors=cofactors+1;
    end;        
end;
add_msg_board(sprintf('Template %s has %i spin labels, %i parametrized cofactors, %i rigid cofactors, and %i gap residues.',...
    ali(1).name,labels,cofactors,blocks));


alg_seq_0=ali(1).sequence;
targ_seq_0=ali(2).sequence;

ali2=ali;
ali2(1).name='alg_template';
ali2(2).db='MMM';
ali2(2).sequence=temp_seq_rep;
ali2(2).name='pdb_template';
ali2=test_repair_alignment(ali2,true);
ident=get_identity(ali2);
add_msg_board(sprintf('Template PDB sequence has %4.1f%% identity with template sequence used in alignment.',100*ident));
% find record of pdb template sequence and template sequence from alignment
% necessary since MUSCLE may change order of sequences
atemp=0;
ptemp=0;
for k=1:2,
    if strcmpi(ali2(k).name,'pdb_template'),
        ptemp=k;
    end;
    if strcmpi(ali2(k).name,'alg_template'),
        atemp=k;
    end;
end;
pdb_seq=ali2(ptemp).sequence;
alg_seq=ali2(atemp).sequence;
% make translation table pdb template to intermediate alignment
poi=0;
trans1=-ones(1,length(temp_seq));
for k=1:length(pdb_seq),
    if char(pdb_seq(k))~='-',
        poi=poi+1;
        while char(temp_seq_rep(poi))=='-',
            poi=poi+1;
            if poi>length(temp_seq), break; end;
        end;
        test=char(pdb_seq(k))==char(temp_seq_rep(poi));
        if ~test,
            error('MMM:fit_from_template_Modeller_adjust_sequence',...
                'Inexplicable sequence mismatch in intermediate alignment 1 at position %i (%i)', k, poi);
        end;
        trans1(poi)=k;
    end;
end;
% make translation table from original alignment template to intermediate
% alignment
trans2=zeros(1,length(alg_seq_0));
poi=0;
for k=1:length(alg_seq),
    poi=poi+1;
    if poi>length(alg_seq_0), break; end;
    while char(alg_seq_0(poi))=='-',
        poi=poi+1;
        if poi>length(alg_seq_0), break; end;
    end;
    if poi>length(alg_seq_0), break; end;
    if char(alg_seq(k))==char(alg_seq_0(poi)),
        trans2(poi)=k;
        test=char(alg_seq_0(poi))==char(alg_seq(k));
        if ~test,
            error('MMM:fit_from_template_Modeller_adjust_sequence',...
                'Inexplicable sequence mismatch in intermediate alignment 2 at position %i (%i)', k, poi);
        end;
    end;
end;
% make translation table from original alignment template to pdb template
trans3=trans1;
for k=1:length(temp_seq),
    poi=find(trans1(k)==trans2);
    if ~isempty(poi),
        trans3(k)=poi;
    else
        trans3(k)=-1;
    end;
    ccode=char(temp_seq_rep(k));
    test=ccode==char(alg_seq_0(poi));
    if ~test,
        if ~isempty(strfind(label_codes,ccode)) && ~isempty(strfind(het_codes,ccode)) && ccode~='.',
            error('MMM:fit_from_template_Modeller_adjust_sequence',...
                'Inexplicable sequence mismatch in final alignment at position %i (%i)', k, poi);
        end;
    end;
end;

% test that all unknown cofactors are at the end
failure=false;
cofac=false;
for k=1:length(trans3),
    if trans3(k)==-1 && char(temp_seq(k))~='-',
        cofac=true;
    end;
    if cofac && trans3(k)>0,
        failure=true;
        cofac=k;
    end;
end;

if failure,
    error('MMM:fit_from_template_Modeller_adjust_sequence',...
        'Program cannot handle inserted cofactors (position %i)', cofac);
end;



temp_seq_0=temp_seq;
temp_seq=zeros(1,length(alg_seq));
targ_seq=zeros(1,length(alg_seq));

ind=0;
for k=1:length(alg_seq_0),
    poi=find(trans3==k);
    if isempty(poi),
        if char(targ_seq_0(k))~='-',
            ind=ind+1;
            temp_seq(ind)='-';
            targ_seq(ind)=targ_seq_0(k);
        end;
    else
        ind=ind+1;
        temp_seq(ind)=temp_seq_0(poi);
        targ_seq(ind)=targ_seq_0(k);
    end;
end;
poi=find(trans3==-1);
for k=poi,
    ind=ind+1;
    temp_seq(ind)=temp_seq_0(k);
    targ_seq(ind)=temp_seq_0(k);
end;
k=1;
while k<=length(temp_seq) && char(temp_seq(k))=='-', % cutoff of template gap at the beginning
    k=k+1;
end;
while ind>1 && char(temp_seq(ind))=='-',
    ind=ind-1;
end;
ind1=ind;
while ind1>1 && char(temp_seq(ind))=='.' && char(targ_seq(ind))~='.',
    targ_seq(ind)='-';
    ind1=ind1-1;
end;

darestarg=0;
for kk=1:k-1,
    if char(targ_seq(kk))~='-',
        darestarg=darestarg+1;
    end;
end;

temp_seq=char(temp_seq(k:ind));
targ_seq=char(targ_seq(k:ind));

clean_targ_seq=targ_seq;
poi=0;
for k=1:length(targ_seq),
    if targ_seq(k)~='-',
        poi=poi+1;
        clean_targ_seq(poi)=targ_seq(k);
    end;
end;
clean_targ_seq=clean_targ_seq(1:poi);
% copy labels
for k=1:length(temp_seq),
    if ~isempty(strfind(label_codes,temp_seq(k))) && isempty(strfind('-.',targ_seq(k))),
        targ_seq(k)=temp_seq(k);
    end;        
end;

function restraints=label_all_sites(restraints,translation,do_attach)
% Given a restraint array, all spin labels are computed that are required
% for specifying DEER restraints
% restraints for which the labeling sites do not exist are removed
% labels are attached, if requested
%
% restraints    array of existing restraints
% translation   translation tables from target to template addresses
% do_attach     (optional) flag that determines whether labels are
%               attached, true: labels are attached, defaults to false
%

global model
global hMain

snum=model.current_structure;

if isempty(restraints.DEER),
    return
else
    DEER=restraints.DEER;
    poi=0;
end;

if nargin<2,
    do_attach=false;
end;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

if ~isempty(DEER),
    add_msg_board(sprintf('Generating labels for %i DEER restraints',length(DEER)));
end;
if ~isempty(labels),
    add_msg_board(sprintf('Structure already contains %i labels',length(labels)));
end;

allowed=0;
% check whether sites are already labelled and whether all restraint sites
% do exist
% identity of the label is checked
% labeling temperature is NOT checked
T_list=zeros(1,200);
if ~isempty(labels),
    lindices=zeros(length(labels),4);
    for k=1:length(labels),
        cindices=labels(k).indices;
        if ~isempty(cindices),
            lindices(k,:)=cindices;
        end;
    end;
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        tempadr1=translate_address(adr1,translation);
        ind1=resolve_address(tempadr1);
        if isempty(ind1),
            add_msg_board(sprintf('Warning: Constraint %i removed as it has first label at site %s',k,tempadr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        end;
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(tempadr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=tempadr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,tempadr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        tempadr2=translate_address(adr2,translation);
        ind2=resolve_address(tempadr2);
        if isempty(ind2),
            add_msg_board(sprintf('Warning: Constraint %i removed as it has second label at site %s',k,tempadr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        end;
        found=false;
        for l=1:length(labels),
            diff=ind2-lindices(l,:);
            if sum(abs(diff))==0  && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(tempadr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=tempadr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,tempadr2));
            end;
        end;
        allowed=allowed+1;
        DEER(allowed)=restraints.DEER(k);
    end;
    if ~isempty(to_do_list),
        add_msg_board('Warning: New labels are computed in a structure that may already contain attached spin labels.');
        add_msg_board('This may lead to erroneous rotamer distributions.');
        add_msg_board('Consider starting from an unlabeled template or from a template where rotamers are computed but labels are not attached.');
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        tempadr1=translate_address(adr1,translation);
        ind1=resolve_address(tempadr1);
        if isempty(ind1),
            add_msg_board(sprintf('Warning: Constraint %i removed as it has first label at site %s',k,tempadr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        else
            found=false;
            for l=1:length(to_do_list),
                if strcmp(tempadr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=tempadr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label at site %s will be generated.',tempadr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        tempadr2=translate_address(adr2,translation);
        ind2=resolve_address(tempadr2);
        if isempty(ind2),
            add_msg_board(sprintf('Warning: Constraint %i removed as it has second label at site %s',k,tempadr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        else
            found=false;
            for l=1:length(to_do_list),
                if strcmp(tempadr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=tempadr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label at site %s will be generated.',tempadr2));
            end;
        end;
        allowed=allowed+1;
        DEER(allowed)=restraints.DEER(k);
    end;
end;
if allowed<length(restraints.DEER),
    add_msg_board(sprintf('Warning: %i out of %i DEER restraints had to be removed.',...
        length(restraints.DEER)-allowed,length(restraints.DEER)));
else
    add_msg_board('All DEER restraints can be kept.');
end;

restraints.DEER=DEER(1:allowed);

if ~isempty(to_do_list),
    stag=model.info{snum}.idCode;
    k=1;
    newtag=sprintf('%s_%i',stag,k);
    while ~isempty(tag2id(newtag,model.structure_tags)),
        k=1+1;
        newtag=sprintf('%s_%i',stag,k);
    end;
    snum=copy_structure(model.current_structure,newtag);
    model.current_structure=snum;
    for k=1:length(restraints.DEER),
        restraints.DEER(k).indices(1,1)=snum;
        restraints.DEER(k).indices(2,1)=snum;
    end;
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
    model.coarse(snum).Ca_coor=Ca_coor;
    model.coarse(snum).indices=rindices;
    model.coarse(snum).masses=masses;
    model.coarse(snum).Bfactors=Bfactors;
    model.coarse(snum).restypes=restypes;
end;


for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

if do_attach,
    for k=1:length(to_do_list),
        if ~strcmp(to_do_list{k},' '),
            command=sprintf('label %s %s %i',to_do_list{k},label_list{k},T_list(k));
            hMain.store_undo=false;
            cmd(hMain,command);
        end;
    end;
end;

function translation=trans_table(alignment,restraints)
% returns a chainwise residue translation table from target residue numbers
% and chain tags to template residue numbers and chain tags
% an empty table is returned if the alignment or restraints record is not
% self-consistent

global model

translation=[];

if isempty(alignment); return; end;

if length(alignment)<2,
    translation(1).table=[];
    translation(1).chain='';
    translation(1).tchain='';
    return;
end;

chain_tags=model.chain_tags{model.current_structure};
tchains=length(model.structures{model.current_structure});
tempseq=alignment(1).sequence;
nonsense=textscan(tempseq,'%s','Delimiter','/');
tempchains=nonsense{1};
chains=length(tempchains);
targseq=alignment(2).sequence;
nonsense=textscan(targseq,'%s','Delimiter','/');
targchains=nonsense{1};
chains2=length(targchains);
if length(tempseq)~=length(targseq),
    add_msg_board('ERROR: Different lengths of template and target sequence');
    return
end;
if chains2~=chains,
    add_msg_board('ERROR: Different number of chains in template and target sequence');
    return
end;
if chains>tchains,
    add_msg_board('ERROR: More chains in alignment than in template structure');
    return
end;
if isfield(restraints,'chains'),
    tag=id2tag(chains,restraints.chains);
    if isempty(tag),
        add_msg_board('ERROR: More chains in alignment than specified in restraint file.');
        return
    end;
    targ_tags=restraints.chains;
else
    targ_tags=':';
    for k=1:chains,
        targ_tags=[targ_tags char(double('A')-1+k) ':'];
    end;
end;

if strcmp(alignment(1).db,'Modeller') && isfield(alignment(1),'first_res'),
    offset=str2double(alignment(1).first_res)-1;
else
    offset=0;
end;
translation(1).offset=offset;
if strcmp(alignment(2).db,'Modeller') && isfield(alignment(2),'first_res'),
    offset2=str2double(alignment(2).first_res)-1;
else
    offset2=0;
end;
translation(1).targ_offset=-offset2;
for k=1:chains,
    chain_tag=id2tag(k,chain_tags);
    translation(k).chain=chain_tag;
    targ_tag=id2tag(k,targ_tags);
    translation(k).tchain=targ_tag;
    ctempseq=char(tempchains(k));
    ctargseq=char(targchains(k));
    if length(ctempseq)~=length(ctargseq),
        add_msg_board(sprintf('ERROR: Different sequence lengths in alignment for chain %i.',k));
        translation=[];
        return
    end;
    ptemp=0;
    ptarg=0;
    if k==1,
        ptemp=ptemp+offset;
        ptarg=ptarg+offset2;
    end;
    transtab=zeros(1,2000);
    for kk=1:length(ctempseq),
        if char(ctempseq(kk))~='-',
            ptemp=ptemp+1;
            present=true;
        else
            present=false;
        end;
        if char(ctargseq(kk))~='-',
            ptarg=ptarg+1;
            if present,
                transtab(ptarg)=ptemp;
            end;
        end;
    end;
    translation(k).table=transtab(1:ptarg);
end;

function tempadr=translate_address(adr,translation)
% translates an target address to a template address
% if no chain is specified, the target address refers to the first chain
% if the chain identifier or the residue number is invalid, an empty string
% is returned

if isempty(translation),
    tempadr=adr;
    return
end;

chain_tags=':';
for k=1:length(translation),
    chain_tags=[chain_tags translation(k).chain ':'];
end;
cprefix=strfind(adr,'(');
csuffix=strfind(adr,')');
if ~isempty(csuffix),
    if csuffix<length(adr),
        resnum=str2double(adr(csuffix+1:end));
    else
        return;
    end;
else
    resnum=str2double(adr);
end;
if isnan(resnum),
    return;
end;

if ~isempty(cprefix) && ~isempty(csuffix),
    ctag=adr(cprefix+1:csuffix-1);
    cid=tag2id(ctag,translation.targ_tags);
else
    cid=1;
end;
ctag=id2tag(cid,chain_tags);

if resnum<=length(translation(cid).table),
    newresnum=translation(cid).table(resnum);
    if newresnum>0,
        if isempty(ctag),
            tempadr=sprintf('%i',newresnum);
        else
            tempadr=sprintf('(%s)%i',ctag,newresnum);
        end;
    else
        tempadr='';
    end;
else
    tempadr='';
end;
% fprintf(1,'targ2temp: |%s| -> |%s|\n',adr,tempadr);

function [targadr,cid,ctag,resnum]=target_address_to_modeller(adr,translation)
% translates a target address to a Modeller PDB output address
% if no chain is specified, the target address refers to the first chain
% if the chain identifier or the residue number is invalid, an empty string
% is returned


if isempty(translation),
    targadr=adr;
end;

targadr='';

targ_tags=':';
chains=length(translation);
for k=1:length(translation),
    targ_tags=[targ_tags translation(k).tchain ':'];
end;
cprefix=strfind(adr,'(');
csuffix=strfind(adr,')');
if ~isempty(csuffix),
    if csuffix<length(adr),
        resnum=str2double(adr(csuffix+1:end));
    else
        return;
    end;
else
    resnum=str2double(adr);
end;
if isnan(resnum),
    return;
end;

% ### the following is a translation to template residue numbers ###
% ### don't activate this, except for tests ###
% if isfield(translation,'table'),
%     resnum=translation.table(resnum);
% end;

if ~isempty(cprefix) && ~isempty(csuffix),
    specify_chain=true;
    ctag=adr(cprefix+1:csuffix-1);
    cid=tag2id(ctag,targ_tags);
    if isempty(cid),
        return;
    end;
else
    specify_chain=false;
    cid=1;
end;

offsets=zeros(1,chains);
if specify_chain,
    new_targ_tags=':';
    for k=1:chains,
        new_targ_tags=[new_targ_tags char(double('A')-1+k) ':'];
        if isfield(translation,'targ_offset'),
            offsets(k)=translation(k).targ_offset;
        else
            offsets(k)=0;
        end;
    end;
    ctag=id2tag(cid,new_targ_tags);
else
    if isfield(translation,'targ_offset'),
        offsets=translation(1).targ_offset;
    else
        offsets=0;
    end;
    ctag='';
end;

% resnum=resnum+offsets(cid);
if cid==1 && isfield(translation,'start_targ'),
    resnum=resnum-translation(1).start_targ;
end;
if isempty(ctag),
    targadr=sprintf('%i',resnum);
else
    targadr=sprintf('(%s)%i',ctag,resnum);
end;
%fprintf(1,'targ2mod: |%s| -> |%s|\n',adr,targadr);


function alignment=realign(alignment,restraints)
% Realigns two sequences given in an MMM alignment record considering an
% optional list of alignment restraints (matching residues)

global model

chain_tags=model.chain_tags{model.current_structure};
tchains=length(model.structures{model.current_structure});
tempseq=alignment(1).sequence;
nonsense=textscan(tempseq,'%s','Delimiter','/');
tempchains=nonsense{1};
chains1=length(tempchains);
targseq=alignment(2).sequence;
nonsense=textscan(targseq,'%s','Delimiter','/');
targchains=nonsense{1};
chains2=length(targchains);

if tchains~=chains1,
    add_msg_board('Warning: Realignment failed since number of chains in alignment file does not match number of chains in template');
    return;
end;
if chains2~=chains1,
    add_msg_board('Warning: Realignment failed since number of chains in target does not match number of chains in template');
    return;
end;

aligned_pairs=zeros(500,3);
apoi=0;
if isfield(restraints,'aligned'),
    if isfield(restraints,'chains'),
        targ_tags=restraints.chains;
    else
        targ_tags=':';
        for k=1:chains2,
            targ_tags=[targ_tags char(double('A')-1+k) ':'];
        end;
    end;
    for k=1:length(restraints.aligned),
        [cid1,resnum1]=dissect_address(restraints.aligned(k).adr1,chain_tags);
        if isempty(cid1) || isempty(resnum1),
            add_msg_board('Warning: Realignment failed due to incorrect specification of aligned template residue.');
            add_msg_board(sprintf('Incorrect address: %s',restraints.aligned(k).adr1));
            return;
        end;
        [cid2,resnum2]=dissect_address(restraints.aligned(k).adr2,targ_tags);
        if isempty(cid2) || isempty(resnum2),
            add_msg_board('Warning: Realignment failed due to incorrect specification of aligned target residue.');
            add_msg_board(sprintf('Incorrect address: %s',restraints.aligned(k).adr2));
            return;
        end;
        if cid1~=cid2,
            add_msg_board('Warning: Realignment failed due to aligned residues in non-matching chains.');
            add_msg_board(sprintf('Incorrect addresses: %s, %s',restraints.aligned(k).adr1,restraints.aligned(k).adr2));
            return;
        end;
        apoi=apoi+1;
        aligned_pairs(apoi,:)=[cid1,resnum1,resnum2];
    end;
end;
aligned_pairs=aligned_pairs(1:apoi,:);

tempseq='';
targseq='';
for k=1:chains1,
    apoic=0;
    curr_align=zeros(500,2);
    ctempchain0=char(tempchains(k));
    ctargchain0=char(targchains(k));
    ctempchain='';
    ctargchain='';
    for kk=1:length(ctempchain0),
        if ctempchain0(kk)~='-',
            ctempchain=[ctempchain ctempchain0(kk)];
        end;
    end;
    for kk=1:length(ctargchain0),
        if ctargchain0(kk)~='-',
            ctargchain=[ctargchain ctargchain0(kk)];
        end;
    end;
    for kk=1:apoi,
        if aligned_pairs(kk,1)==k,
            apoic=apoic+1;
            curr_align(apoic,:)=aligned_pairs(kk,2:3);
        end;
    end;
    curr_align=curr_align(1:apoic,:);
    if isempty(curr_align),
        seqs{1}=ctempchain;
        seqs{2}=ctargchain;
        [message,outfile]=align_sequences(seqs,'',true);
        if message.error,
            add_msg_board(sprintf('Warning: Realignment failed due to MUSCLE alignment error: %s.',message.text));
        end;
        cali=get_multiple_clustal(outfile);
        if strcmpi(cali(1).name,'MMM|seq1'),
            tempseq=[tempseq cali(1).sequence '/'];
            targseq=[targseq cali(2).sequence '/'];
        else
            tempseq=[tempseq cali(2).sequence '/'];
            targseq=[targseq cali(1).sequence '/'];
        end;
    else
        temppoi=0;
        targpoi=0;
        for kk=1:apoic,
            al1=curr_align(kk,1);
            al2=curr_align(kk,2);
            if al1>temppoi+1 && al2>targpoi+1,
                seqs{1}=ctempchain(temppoi+1:al1-1);
                seqs{2}=ctargchain(targpoi+1:al2-1);
                [message,outfile]=align_sequences(seqs,'',true);
                if message.error,
                    add_msg_board(sprintf('Warning: Realignment failed due to MUSCLE alignment error: %s.',message.text));
                end;
                cali=get_multiple_clustal(outfile);
                if strcmpi(cali(1).name,'MMM|seq1'),
                    tempseq=[tempseq cali(1).sequence];
                    targseq=[targseq cali(2).sequence];
                else
                    tempseq=[tempseq cali(2).sequence];
                    targseq=[targseq cali(1).sequence];
                end;
            elseif al1>temppoi+1,
                shelp='';
                ah=al2;
                while ah<targpoi+1,
                    ah=ah+1;
                    shelp=[shelp '-'];
                end;
                tempseq=[tempseq ctempchain(temppoi+1:al1-1)];
                targseq=[targseq shelp];
            elseif al2>targpoi+1,
                shelp='';
                ah=al1;
                while ah<=temppoi+1,
                    ah=ah+1;
                    shelp=[shelp '-'];
                end;
                tempseq=[tempseq shelp];
                targseq=[targseq ctargchain(targpoi+1:al2-1)];
            end;
            tempseq=[tempseq ctempchain(al1)];
            targseq=[targseq ctargchain(al2)];
            temppoi=al1;
            targpoi=al2;
        end;
        if temppoi<length(ctempchain) && targpoi<length(ctargchain),
            seqs{1}=ctempchain(temppoi+1:end);
            seqs{2}=ctargchain(targpoi+1:end);
            [message,outfile]=align_sequences(seqs,'',true);
            if message.error,
                add_msg_board(sprintf('Warning: Realignment failed due to MUSCLE alignment error: %s.',message.text));
            end;
            cali=get_multiple_clustal(outfile);
            if strcmpi(cali(1).name,'MMM|seq1'),
                tempseq=[tempseq cali(1).sequence '/'];
                targseq=[targseq cali(2).sequence '/'];
            else
                tempseq=[tempseq cali(2).sequence '/'];
                targseq=[targseq cali(1).sequence '/'];
            end;
        elseif temppoi<length(ctempchain),
            shelp='';
            tempc=temppoi;
            while tempc<length(ctempchain),
                tempc=tempc+1;
                shelp=[shelp '-'];
            end;
            tempseq=[tempseq ctempchain(temppoi:end) '/'];
            targseq=[targseq shelp '/'];
        else
            shelp='';
            targc=targpoi;
            while targc<length(ctargchain),
                targc=targc+1;
                shelp=[shelp '-'];
            end;
            tempseq=[tempseq shelp '/'];
            targseq=[targseq ctargchain(targpoi:end) '/'];
        end;
    end;
end;
alignment(1).sequence=tempseq(1:end-1);
alignment(2).sequence=targseq(1:end-1);

function [cid,resnum,resnum2]=dissect_address(adr,chain_tags)

resnum=[];
resnum2=[];
cid=[];
cprefix=strfind(adr,'(');
csuffix=strfind(adr,')');
hyphen=strfind(adr,'-');
if ~isempty(csuffix),
    if csuffix<length(adr),
        if isempty(hyphen),
            resnum=str2double(adr(csuffix+1:end));
        else
            resnum=str2double(adr(csuffix+1:hyphen-1));
            if hyphen<length(adr),
                resnum2=str2double(adr(hyphen+1:end));
            end;
        end;
    else
        return;
    end;
else
    if isempty(hyphen),
        resnum=str2double(adr);
    else
        resnum=str2double(adr(1:hyphen-1));
        if hyphen<length(adr),
            resnum2=str2double(adr(hyphen+1:end));
        end;
    end;
end;
if isnan(resnum),
    resnum=[];
    return;
end;
if isnan(resnum2),
    resnum2=[];
    return;
end;

if ~isempty(cprefix) && ~isempty(csuffix),
    ctag=adr(cprefix+1:csuffix-1);
    cid=tag2id(ctag,chain_tags);
end;


% --- Executes on button press in checkbox_dynamic_rotamers.
function checkbox_dynamic_rotamers_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dynamic_rotamers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_dynamic_rotamers

function [correspondence,Ca_coor_target]=align_template_target(id_template,id_target)
% returns Calpha coordinates for a target structure and a correspondence
% table between Ca coordinate arrays of template and target structure based
% on sequence alignment with MUSCLE

global model

correspondence=[];
Ca_coor_target=[];

snum_template=resolve_address(id_template);
snum_target=resolve_address([ '[' id_target ']']);
if isempty(snum_target),
    add_msg_board('Warning! Target structure was not loaded.');
    add_msg_board('Comparison with target structure will be skipped.');
    return
end;
rindices1=model.coarse(snum_template).indices;
[Ca_coor_target,masses2,rindices2]=coarse_residues([ '[' id_target ']']);

cid1=model.chain_ids(snum_template);
cid1=cid1{1};
cid2=model.chain_ids(snum_target);
cid2=cid2{1};
if length(cid1)~=length(cid2),
    add_msg_board('ERROR: Different number of chains in template and moving structure.');
    add_msg_board('Deactivate "whole structure".');
    return
end;
sel1=zeros(2000,4);
sel2=sel1;
psel=0;
for k=1:length(cid1),
    seqs{1}=model.structures{snum_template}(cid1(k)).sequence;
    seqs{2}=model.structures{snum_target}(cid2(k)).sequence;
    sindices=[snum_template,cid1(k);snum_target,cid2(k)];
    [message,inname]=align_sequences(seqs,sindices,true);
    if message.error,
        add_msg_board('ERROR: MUSCLE sequence alignment failed.');
        add_msg_board(message.text);
        add_msg_board('Deactivate "whole structure".');
        return
    end;
    alignment=get_multiple_clustal(inname);
    [asel1,asel2]=select_aligned(alignment,sindices);
    [msel,nsel]=size(asel1);
    sel1(psel+1:psel+msel,:)=asel1;
    sel2(psel+1:psel+msel,:)=asel2;
    psel=psel+msel;
end;
sel1=sel1(1:psel,:);
sel2=sel2(1:psel,:);

[mc1,nc1]=size(rindices1);
[mc2,nc2]=size(rindices2);
correspondence=zeros(2,mc1);
poi=0;
for k=1:psel,
    cindices=sel1(k,:);
    diff=repmat(cindices,mc1,1)-rindices1; % check if residue is in Calpha coordinate array of template structure
    adiff=sum(abs(diff),2);
    cpoi1=find(adiff==0);
    if ~isempty(cpoi1), % if index into Calpha coordinate array was found for template
        cindices=sel2(k,:);
        diff=repmat(cindices,mc2,1)-rindices2; % check if residue is in Calpha coordinate array of target structure
        adiff=sum(abs(diff),2);
        cpoi2=find(adiff==0);
        if ~isempty(cpoi2), % corresponding residues exist in both Calpha coordinate arrays
            poi=poi+1;
            correspondence(1,poi)=cpoi1;
            correspondence(2,poi)=cpoi2;
        end;
    end;
end;

correspondence=correspondence(:,1:poi);

function [sel1,sel2]=select_aligned(alignment,sindices)
% select only residues that are aligned

global model

seq1=alignment(1).sequence;
seq2=alignment(2).sequence;
sel1=zeros(length(seq1),4);
sel2=sel1;
poi1=0;
poi2=0;
pois=0;
for k=1:length(seq1),
    match1=false;
    match2=false;
    if char(seq1(k))~='-',
        poi1=poi1+1;
        match1=true;
    end;
    if char(seq2(k))~='-',
        poi2=poi2+1;
        match2=true;
    end;
    if match1 && match2
        tag1=sprintf('%i',poi1);
        id1=tag2id(tag1,model.structures{sindices(1,1)}(sindices(1,2)).residues{1}.residue_tags);
        tag2=sprintf('%i',poi2);
        id2=tag2id(tag2,model.structures{sindices(2,1)}(sindices(2,2)).residues{1}.residue_tags);
        if ~isempty(id1) && ~isempty(id2),
            pois=pois+1;
            sel1(pois,:)=[sindices(1,:) 1 id1];
            sel2(pois,:)=[sindices(2,:) 1 id2];
        end;
    end;
end;
sel1=sel1(1:pois,:);
sel2=sel2(1:pois,:);


% --- Executes on button press in pushbutton_parametrize.
function pushbutton_parametrize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parametrize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general
global third_party
global help_files
global hMain

parametrization=false;

tightness=0.042; % tightness of DEER distance restraint fitting, 0.048 is optimized for uncertainity 0.3 nm
tightness_Ca=0.040; % tightness of Ca distance restraint fitting, 0.040 is optimized

if parametrization,
    opt_Ca=0.028/4;
    for k=1:9, 
        weightings.Ca(k)=opt_Ca; 
        opt_Ca=sqrt(2)*opt_Ca; 
    end;

    opt_DEER=0.060/4;
    for k=1:9, 
        weightings.DEER(k)=opt_DEER; 
        opt_DEER=sqrt(2)*opt_DEER; 
    end;
else
    weightings.Ca=tightness_Ca;
    weightings.DEER=tightness;
end;

if isempty(handles.DEER) && isempty(handles.direct),
    if isempty(handles.alignment) || handles.template_info.alg==handles.target_info.alg,
        add_msg_board('ERROR: No experimental restraints and target sequence identical to template sequence.');
        add_msg_board('Modeling is aborted.');
        return
    end;
end;

if handles.target_info.identity<0.2,
    add_msg_board('### Strong warning ###: Insufficient identity of target and template sequence.');
    add_msg_board('Homology modeling is most likely unreliable.');
elseif handles.target_info.identity<0.4,
    add_msg_board('Warning: Identity of target and template sequence lower than 40%%.');
    add_msg_board('Homology modeling may be unreliable.');
end;

entry=strcat(help_files,'third_party.html#Modeller');

dospath=which([third_party.modeller_version  '.exe']);
if isempty(dospath),
    message.error=2;
    message.text='Modeller software not found on Matlab path.';
    add_msg_board('This feature requires Modeller from the Sali lab');
    add_msg_board('ERROR: Modeller could not be found on the Matlab path.');
    add_msg_board('Please check whether Modeller is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);

handles.test_mode=false;
if ~isempty(handles.target),
    snum=model.current_structure;
    [coarse_correspondence,coarse_target]=align_template_target(mk_address(snum),handles.target);
    if ~isempty(coarse_correspondence),
        handles.test_mode=true;
        if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
            [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
            if isempty(Ca_coor) || length(masses)<2,
                add_msg_board('Error: Cannot generate test information, since current structure has less than two amino acid residues.');
                handles.test_mode=false;
            else
                handles.coarse_template=Ca_coor;
                model.coarse(snum).Ca_coor=Ca_coor;
                model.coarse(snum).indices=rindices;
                model.coarse(snum).masses=masses;
                model.coarse(snum).Bfactors=Bfactors;
                model.coarse(snum).restypes=restypes;
            end;
        else
            handles.coarse_template=model.coarse(snum).Ca_coor;
            handles.coarse_target=coarse_target;
            handles.coarse_correspondence=coarse_correspondence;
        end;
    end;
end;

hwin=gcf;
set(hwin,'Pointer','watch');

idCode=model.info{model.current_structure}.idCode;
if isempty(idCode), idCode='AMMM'; end;
idCode(1)=char(idCode(1)+16);
filename=sprintf('%s.pdb',handles.temp_ID);

fname=fullfile(general.tmp_files, filename);
msg=sprintf('Template structure saved as PDB file: %s',fname);
add_msg_board(msg);
if handles.template_info.mode>0,
    [message,info]=wr_pdb_paradigm(fname,idCode,1,handles.template_info.mode);
else
    [message,info]=wr_pdb_paradigm(fname,idCode);
end;
if message.error,
    add_msg_board(message.text);
    set(hwin,'Pointer','arrow');
    return
end;
for k=1:length(info.offsets),
    handles.translation(k).targ_offset=info.offsets(k);
    if ~isempty(handles.alignment) && isfield(handles.alignment(1),'first_res'),
        handles.translation(k).targ_offset=...
            handles.translation(k).targ_offset+str2double(handles.alignment(1).first_res)-1;
    end;
    if isfield(handles,'restraints') && isfield(handles.restraints,'chains'),
        targ_tag=id2tag(k,handles.restraints.chains);
        if isempty(targ_tag),
            add_msg_board('ERROR: More chains in template file than specified in restraint file.');
            return
        end;
    else
        targ_tag=char(double('A')-1+k);
    end;    
    handles.translation(k).tchain=targ_tag;
end;
if length(info.offsets)==1,
    handles.translation(1).tchain='';
end;
temp_seq=info.seq;
targ_seq='';
arestarg=[];
if length(handles.alignment)>1,
    ali(1)=handles.alignment(handles.template_info.alg);
    ali(2)=handles.alignment(handles.target_info.alg);
    if isfield(ali,'first_res'),
        offset=str2double(ali(1).first_res)-info.first;
        if offset<0,
            add_msg_board(sprintf('Warning: Eliminating from alignment %i N-terminal residues that are missing in template',-offset));
            ali(1).first_res=sprintf('%i',info.first);
        else
            info.first=info.first+offset;
        end;
    end;
    if isfield(ali,'first_res') && str2double(ali(2).first_res)>0,
        arestarg=str2double(ali(2).first_res);
    else
        arestarg=1;
    end;
    [temp_seq,targ_seq,darestarg,clean_targ_seq]=adjust_alignment(ali,temp_seq);
    if isfield(ali,'first_res')
        info.first=str2double(ali(1).first_res);
        if isfield(ali,'last_res')
            info.last=str2double(ali(1).last_res);
        end;
    end;
    arestarg=arestarg+darestarg;
    handles.translation(k).targ_offset=1-arestarg;
    handles.alignment(handles.target_info.alg).first_res=sprintf('%i',arestarg);
else
    if ~isempty(handles.alignment),
        clean_targ_seq=handles.alignment(handles.template_info.alg).sequence;
    else
        clean_targ_seq=handles.clean_targ_seq(1-handles.translation(k).targ_offset:end);
        if length(clean_targ_seq)>length(temp_seq),
            clean_targ_seq=clean_targ_seq(1:length(temp_seq));
        end;
        targ_seq=clean_targ_seq;
        for k=1:length(targ_seq),
            if char(temp_seq(k))=='!' || char(temp_seq(k))==';',
                targ_seq(k)=temp_seq(k);
            end;
        end;
        arestarg=info.first;
        handles.translation(1).targ_offset=1-arestarg;
        handles.translation(1).chain=info.first_chain;
        % arestarg=1;
%        handles.alignment(handles.target_info.alg).first_res=sprintf('%i',arestarg);
    end;
end;
handles.ares=arestarg;
handles.clean_targ_seq=clean_targ_seq;
basname=sprintf('%s_to_%s',handles.temp_ID,handles.targ_ID); % ###
[alignfile,message]=seq2pir(basname,handles.temp_ID,handles.targ_ID,temp_seq,info.first,info.first_chain,info.last,info.last_chain,targ_seq,arestarg);

if message.error,
    add_msg_board(message.text);
    add_message_board('ERROR: Modelling aborted.');
    set(hwin,'Pointer','arrow');
    return
end;

add_msg_board(sprintf('Alignment file saved as: %s',alignfile));
[algpath, algfile, algext] = fileparts(alignfile);
alignfile=strcat(algfile,algext);

runfilename=sprintf('model_%s_to_%s.py',handles.temp_ID,handles.targ_ID);
runfile=fullfile(general.tmp_files, runfilename);
if ~isfield(handles,'restraints'),
    handles.restraints=[];
end;

uncertainty=handles.uncertainty;
esize=handles.ensemble;
DEER0=handles.restraints.DEER;
coverages=zeros(esize,length(weightings.DEER));

par_logname=sprintf('%s_Modeller_parametrization_50.log',handles.target);
my_dir=pwd;
fid=fopen(par_logname,'w');
fprintf(fid,'Parametrization file for target %s\n',handles.target);

for kp=1:length(weightings.DEER),
    weighting.DEER=weightings.DEER(kp);
    weighting.Ca=weightings.Ca(kp);
    for ke=1:esize,
        fprintf(fid,'\nConstraint set #%i\n',ke);
        DEER=DEER0;
        for kd=1:length(DEER),
            DEER(kd).r=DEER0(kd).r+uncertainty*randn;
            fprintf(fid,'%12s %12s %5.2f (original: %5.2f)\n',DEER(kd).adr1,DEER(kd).adr2,DEER(kd).r,DEER0(kd).r); 
        end;
        handles.restraints.DEER=DEER;
        message=mk_modeller_input(runfile,alignfile,handles.temp_ID,handles.targ_ID,handles.restraints,handles.translation,1,weighting);

        if message.error,
            add_msg_board(message.text);
            add_message_board('ERROR: Modelling aborted.');
            set(hwin,'Pointer','arrow');
            return
        end;

        msg='Running Modeller. Please be patient...';

        [batcmd,message]=mk_modeller_bat(modpath,modcmd,runfile,false);
        if message.error,
            add_msg_board('ERROR: Modeller could not be initialized');
            add_msg_board(message.text);
            set(hwin,'Pointer','arrow');
            return
        end;

        set(handles.text_info,'String',msg);
        add_msg_board(msg);
        drawnow;
        cd(modpath);
        [s, w] = dos(batcmd);
        if s~=0,
            rem=w;
            while ~isempty(rem),
                [token,rem]=strtok(rem,char(10));
                if ~isempty(token),
                    add_msg_board(token);
                end;
            end;
            message.error=2;
            message.text='Modeller error.';
            add_msg_board('ERROR: Modeller did not run successfully.');
            set(hwin,'Pointer','arrow');
            cd(my_dir);
            return
        else
            logname=sprintf('model_%s_to_%s.log',handles.temp_ID,handles.targ_ID);
            msg=sprintf('Modeller job %s_to_%s completed.',handles.temp_ID,handles.targ_ID);
            add_msg_board(msg);
            set(handles.text_info,'String',msg);
            analyze_modeller_log(logname,true);
            coverages(ke,kp)=get_coverage(handles,logname);
        end;
    end;
end;
fprintf(fid,'\nParametrization result for ensemble size %i with uncertainty %5.2f ?\n\n',esize,10*uncertainty);
n_direct=length(handles.direct);
n_DEER=length(handles.DEER);
fprintf(fid,'%i Calpha constraints\n',n_direct);
fprintf(fid,'%i DEER constraints\n\n',n_DEER);
[ne,np]=size(coverages);
for kp=1:np,
    if n_direct>0,
        fprintf(fid,'w_Ca = %5.3f, coverage(s): ',weightings.Ca(kp));
        for ke=1:ne,
            fprintf(fid,'%7.3f',coverages(ke,kp));
        end;
        fprintf(fid,'\n');
    end;
    if n_DEER>0,
        fprintf(fid,'w_DEER = %5.3f, coverage(s): ',weightings.DEER(kp));
        for ke=1:ne,
            fprintf(fid,'%7.3f',coverages(ke,kp));
        end;
        fprintf(fid,'\n');
    end;
end;
fprintf(fid,'\n\nMean coverages for the ensemble:\n');
mcon=mean(coverages,1);
for kp=1:length(mcon),
    fprintf(fid,'%6.3f\t',mcon(kp));
end;
fprintf(fid,'\n');
fclose(fid);
hMain.report_file=par_logname;
cd(my_dir);
report_editor;


set(hwin,'Pointer','arrow');

guidata(hObject,handles);

function coverage=get_coverage(handles,logname)

global general
global model

pname=general.tmp_files;
idCode=handles.targ_ID;
ensemble=assess_ensemble(logname,'',handles.ensemble_size,0.70);
if isempty(ensemble),
    msg='Warning: Model has GA341 lower than 0.70.';
    set(handles.text_info,'String',msg);
    add_msg_board(msg);
    ensemble=assess_ensemble(logname,'',handles.ensemble_size,0.30);
end;
ensemble_size=length(ensemble);
if ensemble_size < 1,
    add_msg_board(sprintf('Warning: only %i models have sufficient GA341 score, while %i models were requested.',ensemble_size,handles.ensemble_size));
    add_msg_board('Importing smaller ensemble. Consider increasing the number of models.');
end;
[message,snum,stag,models]=rd_pdb_ensemble(ensemble,idCode,pname);
ares=[];
if length(handles.alignment)>1,
    ali=handles.alignment(handles.target_info.alg);
    if isfield(ali,'first_res') && str2double(ali.first_res)>0,
        ares=str2double(ali.first_res);
    end;
elseif isfield(handles,'ares'),
    ares=handles.ares;
end;
if isfield(handles.restraints,'chains'),
    chain_tags=handles.restraints.chains;
else
    chain_tags='';
end;
correct_Modeller_ensemble(snum,ares,chain_tags);
if message.error,
    add_msg_board(message.text);
else
    add_msg_board(sprintf('Ensemble read into structure %i with tag [%s].',snum,stag));
    add_msg_board(sprintf('This ensemble has %i chain models',models));
end;
ensemble_address=mk_address(snum);
model_address=sprintf('%s{1}',ensemble_address);
[Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(model_address);
model.coarse(snum).Ca_coor=Ca_coor;
model.coarse(snum).indices=rindices;
model.coarse(snum).masses=masses;
model.coarse(snum).Bfactors=Bfactors;
model.coarse(snum).restypes=restypes;
[coarse_correspondence,coarse_target]=align_template_target(mk_address(snum),handles.target);
add_msg_board(sprintf('Generating test information for ensemble %s with respect to target %s...',ensemble_address,handles.target));
rms0=rmsd_superimpose(handles.coarse_target(handles.coarse_correspondence(2,:),:),handles.coarse_template(handles.coarse_correspondence(1,:),:));
model_address=sprintf('%s{1}',ensemble_address);
[network]=coarse_residues(model_address);
rms=rmsd_superimpose(coarse_target(coarse_correspondence(2,:),:),network(coarse_correspondence(1,:),:));
coverage=(rms0-rms)/rms0;
