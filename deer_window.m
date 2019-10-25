function varargout = deer_window(varargin)
% DEER_WINDOW M-file for deer_window.fig
%      DEER_WINDOW, by itself, creates a new DEER_WINDOW or raises the existing
%      singleton*.
%
%      H = DEER_WINDOW returns the handle to a new DEER_WINDOW or the handle to
%      the existing singleton*.
%
%      DEER_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEER_WINDOW.M with the given input arguments.
%
%      DEER_WINDOW('Property','Value',...) creates a new DEER_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deer_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deer_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deer_window

% Last Modified by GUIDE v2.5 13-Nov-2015 14:28:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deer_window_OpeningFcn, ...
                   'gui_OutputFcn',  @deer_window_OutputFcn, ...
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


% --- Executes just before deer_window is made visible.
function deer_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deer_window (see VARARGIN)

% Choose default command line output for deer_window
handles.output = hObject;

global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.updated=0;
handles.bckg_dim=3;
set(handles.edit_bckg,'String',sprintf('%4.2f',handles.bckg_dim));
handles.mod_depth=0.4;
handles.exp_depth=0.4;
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.zero_time=0;
set(handles.edit_zero_time,'String',sprintf('%5i',handles.zero_time));
handles.sel_distr=1;
handles.range=[2,8];
handles.spin_system=[];
handles.labels=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
handles.expanded=false;
handles.flex_color = [0,0.5,0.5];
handles.fit_color = [0,0.5,0];

set(handles.text_mean,'String','no distribution');
set(handles.text_stddev,'String',' ');
set(handles.text_coor,'String',' ');

% Experimental data set defaults
handles.t_orig=[];
handles.v_orig=[];
handles.dt=0.008;
handles.bas_name='';
handles.texp=[];
handles.texp_fit=[];
handles.vexp=[];
handles.vexp_fit=[];
handles.vb={};
handles.rexp=1.5:0.05:10;
handles.dexp=[];
handles.tdip=[];
handles.tsim=[];
handles.cluster=[];
handles.ff_multi=[];
handles.ff=[];
handles.rsim=[];
handles.dexp_fit=[];
handles.dsim=[];
handles.bckg=[];
handles.bckg_k=0;
handles.tmax_min=0;
handles.tmax_opt=0;
handles.rmsd=0;
handles.flex_move = 0;
handles.dr = 0;

% kernel for fast Deer simulations
load('pake_base40_MMM.mat');
kernel=base-ones(size(base)); % kernel format for pcf2deer
handles.Pake_kernel=kernel;
handles.Pake_t=t;
handles.Pake_r=r;
handles.Pake_wd=wd;
handles.texp=t';
handles.tdip=t';

% copy control flags
handles.copy_DEER=0;
handles.copy_distr=0;

% experimental data names

handles.project_dir='';
handles.bas_name='';

[handles,success] = mk_label_list(handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

set(handles.edit_zero_time,'Enable','off');

if ~success,
    prompt = 'Labels per trajectory frame:';
    dlg_title = 'No labels. Trajectory must be loaded.';
    num_lines = 1;
    def = {'2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer),
        lpf = str2double(answer{1});
        if ~isnan(lpf) && lpf>=2,
            lpf=round(lpf);
            [filename,pathname] = uigetfile('.dat','Load trajectory file');
            try
                data = load(fullfile(pathname, filename));
            catch exception
                data=[];
            end;
            [ln,n]=size(data);
            if mod(ln,lpf)~=0,
                add_msg_board(sprintf('Number of data in trajectory file %s is not a multiple of labels per frame',fullfile(pathname, filename)));
            elseif n~=3 && n~=4,
                add_msg_board(sprintf('Trajectory file %s contains no valid coordinates',fullfile(pathname, filename)));
            else
                success=true;
                for k=1:lpf,
                    handles.trajectory{k}=sprintf('#traj%i',k);
                    handles.frames{k}=zeros(ln/lpf,4);
                end;
                for kk=1:ln/lpf,
                    bas=(kk-1)*lpf;
                    for k=1:lpf,
                        coor=ones(1,4)/(sqrt(ln/lpf));
                        coor(1:n)=data(bas+k,:);
                        handles.frames{k}(kk,:)=coor;
                    end;
                end;
            end;
        end;
    end;
    [handles,success]=mk_label_list(handles);
    if ~success,
       msgbox('Model must feature at least two spin labels or selected atoms or a valid trajectory','DEER simulation impossible','error');
       DEER_CloseRequestFcn(handles.DEER, eventdata, handles); 
    else
        guidata(hObject, handles);
    end;
else
    % Update handles structure
    guidata(hObject, handles);
end;

% UIWAIT makes deer_window wait for user response (see UIRESUME)
% uiwait(handles.DEER);


% --- Outputs from this function are returned to the command line.
function varargout = deer_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in listbox_label.
function listbox_label_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_label

global model

if isfield(handles,'trajectory'),
    trj_labels=length(handles.trajectory);
else
    trj_labels=0;
end;

bilabels = 0;
if isfield(handles,'bilabels')
    bilabels = length(handles.bilabels);
end
sel=get(hObject,'Value');
if trj_labels && sel<=trj_labels,
    NOpos=handles.frames{sel};
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name=handles.trajectory{sel};
    msg{1}=name;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
elseif isfield(model,'bisites') && sel <= length(handles.bilabels)+trj_labels,
    NOpos=handles.bilabels(sel-trj_labels).popcoor;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    msg{1}=handles.bilabels(sel-trj_labels).adr;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
elseif isfield(model,'labels') && sel<=length(handles.bilabels) + length(model.labels)+trj_labels,
    NOpos=model.labels(sel-trj_labels).NOpos;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    msg{1}=model.labels(sel-trj_labels-bilabels).adr;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
else
    if isfield(model,'labels'),
        asel=sel-length(model.labels)-bilabels-trj_labels;
    else
        asel=sel;
    end;
    indices=handles.atoms(asel,:);
    [msg0,xyz]=get_atom(indices,'xyz');
    [msg0,pop]=get_atom(indices,'populations');
    pop=pop/sum(pop);
    xyz=pop'*xyz;
    adr=mk_address(indices,true);
    msg{1}=adr;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',xyz(1),xyz(2),xyz(3));
end;
set(handles.text_coor,'String',msg);
addresses=get(hObject,'String');
[msg,absresnum]=get_object(sprintf('%s',addresses{sel}),'absresnum');
if msg.error
    absresnum = 0;
end
set(handles.text_abs_number,'String',sprintf('%i',absresnum));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_label.
function pushbutton_select_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

which_one=get(handles.listbox_label,'Value');
handles.labels = which_one;
handles.ff_multi=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
set(handles.text_selected_distribution,'String','*** none ***');
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_add_label.
function pushbutton_add_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

which_one=get(handles.listbox_label,'Value');
if isempty(handles.labels)
    handles.labels= which_one;
    handles.ff_multi=[];
    handles.tweak_distr=[];
    handles.tweak_rax=[];
else
    curr_labels = handles.labels;
    double_sel = false;
    if min(abs(curr_labels-which_one)) == 0
        double_sel = true;
    end
    if ~double_sel
        handles.labels = [curr_labels which_one];
        handles.ff_multi=[];
        handles.tweak_distr=[];
        handles.tweak_rax=[];
    else
        msgbox('The same label cannot be selected twice.','Double selection ignored','warn');
    end
    poi=length(handles.labels);
    if poi == 2
        if handles.site_list(handles.labels(1)).chain == handles.site_list(handles.labels(2)).chain, % same chain
            set(handles.checkbox_coupled_ensemble,'Value',1);
        else
            set(handles.checkbox_coupled_ensemble,'Value',0);
        end
    end
end
set(handles.text_selected_distribution,'String','*** none ***');
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_select_all.
function pushbutton_select_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addresses=get(handles.listbox_label,'String');
handles.labels=1:length(addresses);
handles.ff_multi=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
poi=length(handles.labels);
if poi == 2
    ind1 = resolve_address(addresses{handles.labels(1)});
    ind2 = resolve_address(addresses{handles.labels(2)});
    if abs(sum(ind1(1:2)-ind2(1:2))) == 0 % same chain
        set(handles.checkbox_coupled_ensemble,'Value',1);
    else
        set(handles.checkbox_coupled_ensemble,'Value',0);
    end
end
set(handles.text_selected_distribution,'String','*** none ***');
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_remove_label.
function pushbutton_remove_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_labels= zeros(1,10);;
which_one=get(handles.listbox_label,'Value');
addresses=get(handles.listbox_label,'String');
adr=addresses{which_one};
npoi=0;
removal=0;
if ~isempty(handles.labels)
    poi=length(handles.labels);
    for k=1:poi
        if strcmp(handles.site_list(handles.labels(k)).adr,adr)
            removal=1;
        else
            npoi=npoi+1;
            new_labels(npoi)=handles.labels(k);
        end % check for double selection of an existing label
    end
end
if removal
    handles.labels=new_labels(1:npoi);
    handles.ff_multi=[];
    handles.tweak_distr=[];
    handles.tweak_rax=[];
    if npoi == 2
        ind1 = resolve_address(handles.site_list(handles.labels(1)).adr);
        ind2 = resolve_address(handles.site_list(handles.labels(2)).adr);
        if abs(sum(ind1(1:2)-ind2(1:2))) == 0 % same chain
            set(handles.checkbox_coupled_ensemble,'Value',1);
        else
            set(handles.checkbox_coupled_ensemble,'Value',0);
        end
    end
else
    msgbox('Only a selected label can be deselected.','Deselection ignored','warn');
end
set(handles.text_selected_distribution,'String','*** none ***');
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_load_Xepr.
function pushbutton_load_Xepr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_Xepr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=get_dataset_MMM(handles);
set(handles.checkbox_form_factor,'Value',0);
set(handles.edit_zero_time,'Enable','on');
handles.ff_multi=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_load_DeerAnalysis.
function pushbutton_load_DeerAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_DeerAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = load_DeerAnalysis_MMM(handles);
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.exp_depth));
set(handles.edit_zero_time,'Enable','off');
set(handles.pushbutton_any_rotamers,'Enable','on');
handles.ff_multi=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
handles=update(handles);
guidata(hObject,handles);

% --- Executes on selection change in popupmenu_bckg.
function popupmenu_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_bckg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bckg

% items:    1 is 3D
%           2 is 2D
%           3 is fractal
%           4 is none

contents = get(hObject,'String');
mode=contents{get(hObject,'Value')};

switch mode
    case '3D'
        set_bckg_edit(hObject,handles,3);
    case '2D'
        set_bckg_edit(hObject,handles,2);
    case 'fractal'
        set_bckg_edit(hObject,handles,-1);
    case 'none'
        set_bckg_edit(hObject,handles,NaN);
end;
handles=update(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_bckg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.DEER_files);

label_list = get(handles.listbox_label,'String');
if isfield(handles,'rsim')
    labstring='';
    for k=1:length(handles.labels)
        labstring=sprintf('%s%s_',labstring,label_list{handles.labels(k)});
    end
    suggestion=[labstring 'res.txt']; 
    [fname,pname]=uiputfile('*.txt','Save results',suggestion);
    if isequal(fname,0) || isequal(pname,0)
        add_msg_board('Saving of distance distribution and DEER simulation canceled by user.');
        return;
    end
    reset_user_paths(pname);
    general.DEER_files=pname;
    % Remove (last) extension, if any
    s=strfind(fname,'.');
    if ~isempty(s)
        fname=fname(1:s(length(s))-1);
    end;
    % Remove suffix '_res', if present
    s=strfind(fname,'_res');
    if ~isempty(s),
        fname=fname(1:s(length(s))-1);
    end;
    handles=save_result_MMM(handles,fname,pname);

else
    msgbox('Select spin system with at least two labels before saving.','Nothing to save','warn')
end;

cd(my_path);
guidata(hObject,handles);

% --- Executes on button press in checkbox_mod_depth.
function checkbox_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mod_depth

handles=update(handles);
guidata(hObject,handles);

function edit_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg as a double
[v,handles]=edit_update_MMM(handles,hObject,0.0,10.0,3.0,'%4.2f',0);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_bckg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mod_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_mod_depth as a double

v0=handles.exp_depth;
[v,handles]=edit_update_MMM(handles,hObject,0.0,1.0,0.4,'%5.3f',0);
handles.mod_depth=v;
handles.exp_depth=v;
if abs(v-v0)>5e-4,
    handles.ff_multi=[];
    handles.tweak_distr=[];
    handles.tweak_rax=[];
    handles=update(handles);
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_mod_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_labeling_efficiencies.
function checkbox_labeling_efficiencies_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_labeling_efficiencies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_labeling_efficiencies


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_detach_distr.
function pushbutton_detach_distr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_distr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy_distr=1;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DEER_CloseRequestFcn(handles.DEER, eventdata, handles);


function set_bckg_edit(hObject,handles,dim)
% updates edit_bckg after popupmenu selection of background function

if isnan(dim),
    set(handles.edit_bckg,'String','n.a.');
else
    if dim<0, dim=handles.bckg_dim; end;
    set(handles.edit_bckg,'String',sprintf('%4.2f',dim));
    handles.bckg_dim=dim;
end;

% Update handles structure
guidata(hObject, handles);


% --- Executes when user attempts to close DEER.
function DEER_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to DEER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function edit_zero_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zero_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zero_time as text
%        str2double(get(hObject,'String')) returns contents of edit_zero_time as a double

[zt0,handles]=edit_update_MMM(handles,hObject,0,max(handles.t_orig),handles.zero_time,'%5i',1);
[texp,vexp,zt,dt]=pre_process_MMM(handles.t_orig,handles.v_orig,zt0);
handles.texp=texp'/1000;
handles.vexp=vexp';
handles.dt=dt;
handles.zero_time=zt;
handles.ff_multi=[];
handles.tweak_distr=[];
handles.tweak_rax=[];
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_zero_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zero_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=update(handles)
% Display update

global model

col1=[0,1,0];
col2=[0,0,1];

if handles.copy_distr
    figure(1); clf;
    set(gca,'FontSize',16);
else
    axes(handles.axes_distribution);
    cla;
    if isfield(handles,'left_crsr')
        handles=rmfield(handles,'left_crsr');
    end
    if isfield(handles,'right_crsr')
        handles=rmfield(handles,'right_crsr');
    end
end
hold on;

if ~isempty(handles.labels)
    lab_string = handles.site_list(handles.labels(1)).adr;
    if length(handles.labels)>1
        for k=2:length(handles.labels)
            lab_string=sprintf('%s; %s',lab_string,handles.site_list(handles.labels(k)).adr);
        end
    end
    set(handles.text_spin_system,'String',lab_string);
else
    set(handles.text_spin_system,'String','');
end    

simdistr=0;
numlab=length(handles.labels);
% for k=1:numlab,
%     [msg,absresnum]=get_object(sprintf('%s',handles.labels{k}),'absresnum');
%     if msg.error,
%         fprintf(1,'Error %i:%s\n',msg.error,message.text);
%     else
%         fprintf(1,'Label at %s has absolute residue number %i\n',handles.labels{k},absresnum);
%     end;
% end;

do_ensemble = get(handles.checkbox_ensemble,'Value');
if numlab > 2 && do_ensemble
    add_msg_board('Warning: Ensemble computation not implemented for multi-spin effects.');
    do_ensemble = false;
    set(handles.checkbox_ensemble,'Value',0);
    set(handles.checkbox_coupled_ensemble,'Enable','off');
end

rmin0=1e6;
rmax0=0;
% Display distance distributions
if numlab>1
    poi=0;
    pairs=numlab*(numlab-1)/2; % number of spin pairs
    for k=1:numlab-1
        currlabel = handles.labels(k);
        % find the Calpha coordinate
        switch handles.site_list(currlabel).type
            case {'label','atom','trajectory'}
                adr1 = handles.site_list(currlabel).adr;
                dadr1 = adr1;
                if do_ensemble
                    ind1 = resolve_address(adr1);
                    ne1 = length(model.structures{ind1(1)}(ind1(2)).xyz);
                    ind1e = ind1;
                    CA1 = [0,0,0];
                    found = false;
                    for ke = 1:ne1
                        ind1e(3) = ke;
                        adr1e = mk_address(ind1e);
                        [~,CA01]=get_object(sprintf('%s.CA',adr1e),'coor');
                        if isempty(CA01) % for nucleotides
                            [~,CA01]=get_object(sprintf('%s.C1''',adr1e),'coor');
                        end
                        if ~isempty(CA01)
                            CA1 = CA1 + CA01;
                            found = true;
                        end
                    end
                    if found
                        CA1 = CA1/ne1;
                    else
                        CA1 = [];
                    end
                else
                    [~,CA1]=get_object(sprintf('%s.CA',adr1),'coor');
                    if isempty(CA1) % for nucleotides
                        [~,CA1]=get_object(sprintf('%s.C1''',adr1),'coor');
                    end
                end
            case {'bilabel'}
                addresses = handles.site_list(currlabel).adr;
                dadr1 = addresses;
                seppoi = strfind(addresses,'|');
                adr1 = addresses(1:seppoi-1);
                adr2 = addresses(seppoi+1:end);
                if do_ensemble
                    ind1 = resolve_address(adr1);
                    ne1 = length(model.structures{ind1(1)}(ind1(2)).xyz);
                    ind1e = ind1;
                    CA1a = [0,0,0];
                    found = false;
                    for ke = 1:ne1
                        ind1e(3) = ke;
                        adr1e = mk_address(ind1e);
                        [~,CA01]=get_object(sprintf('%s.CA',adr1e),'coor');
                        if isempty(CA01) % for nucleotides
                            [~,CA01]=get_object(sprintf('%s.C1''',adr1e),'coor');
                        end
                        if ~isempty(CA01)
                            CA1a = CA1a + CA01;
                            found = true;
                        end
                    end
                    if found
                        CA1a = CA1a/ne1;
                    else
                        CA1a = [];
                    end
                    ind2 = resolve_address(adr2);
                    ne2 = length(model.structures{ind2(1)}(ind2(2)).xyz);
                    ind2e = ind2;
                    CA1b = [0,0,0];
                    found = false;
                    for ke = 1:ne2
                        ind2e(3) = ke;
                        adr2e = mk_address(ind2e);
                        [~,CA01]=get_object(sprintf('%s.CA',adr2e),'coor');
                        if isempty(CA01) % for nucleotides
                            [~,CA01]=get_object(sprintf('%s.C1''',adr2e),'coor');
                        end
                        if ~isempty(CA01)
                            CA1b = CA1b + CA01;
                            found = true;
                        end
                    end
                    if found
                        CA1b = CA1b/ne2;
                    else
                        CA1b = [];
                    end
                    CA1 = (CA1a+CA1b)/2;
                else
                    [~,CA1a]=get_object(sprintf('%s.CA',adr1),'coor');
                    if isempty(CA1a) % for nucleotides
                        [~,CA1a]=get_object(sprintf('%s.C1''',adr1),'coor');
                    end
                    [~,CA1b]=get_object(sprintf('%s.CA',adr2),'coor');
                    if isempty(CA1b) % for nucleotides
                        [~,CA1b]=get_object(sprintf('%s.C1''',adr2),'coor');
                    end
                    CA1 = (CA1a+CA1b)/2;
                end
        end
        for kk=k+1:numlab
            poi=poi+1;
            if pairs>1
                cf=(poi-1)/(pairs-1);
            else
                cf=0;
            end
            col=(1-cf)*col1+cf*col2;
            col=col/max(col);
            currlabel2 = handles.labels(kk);
            % find the Calpha coordinate
            switch handles.site_list(currlabel2).type
                case {'label','atom','trajectory'}
                    adr2= handles.site_list(handles.labels(kk)).adr;
                    dadr2 = adr2;
                    if do_ensemble
                        ind2 = resolve_address(adr1);
                        ne2 = length(model.structures{ind2(1)}(ind2(2)).xyz);
                        ind2e = ind2;
                        CA2 = [0,0,0];
                        found = false;
                        for ke = 1:ne2
                            ind2e(3) = ke;
                            adr2e = mk_address(ind2e);
                            [~,CA02]=get_object(sprintf('%s.CA',adr2e),'coor');
                            if isempty(CA02) % for nucleotides
                                [~,CA02]=get_object(sprintf('%s.C1''',adr2e),'coor');
                            end
                            if ~isempty(CA02)
                                CA2 = CA2 + CA02;
                                found = true;
                            end
                        end
                        if found
                            CA2 = CA2/ne2;
                        else
                            CA2 = [];
                        end
                    else
                        [~,CA2]=get_object(sprintf('%s.CA',adr2),'coor');
                        if isempty(CA2) % for nucleotides
                            [~,CA2]=get_object(sprintf('%s.C1''',adr2),'coor');
                        end
                    end
                case 'bilabel'
                    addresses = handles.site_list(currlabel2).adr;
                    dadr2 = addresses;
                    seppoi = strfind(addresses,'|');
                    adr1 = addresses(1:seppoi-1);
                    adr2 = addresses(seppoi+1:end);
                    if do_ensemble
                        ind1 = resolve_address(adr1);
                        ne1 = length(model.structures{ind1(1)}(ind1(2)).xyz);
                        ind1e = ind1;
                        CA2a = [0,0,0];
                        found = false;
                        for ke = 1:ne1
                            ind1e(3) = ke;
                            adr1e = mk_address(ind1e);
                            [~,CA02]=get_object(sprintf('%s.CA',adr1e),'coor');
                            if isempty(CA02) % for nucleotides
                                [~,CA02]=get_object(sprintf('%s.C1''',adr1e),'coor');
                            end
                            if ~isempty(CA02)
                                CA2a = CA2a + CA02;
                                found = true;
                            end
                        end
                        if found
                            CA2a = CA2a/ne1;
                        else
                            CA2a = [];
                        end
                        ind2 = resolve_address(adr2);
                        ne2 = length(model.structures{ind2(1)}(ind2(2)).xyz);
                        ind2e = ind2;
                        CA2b = [0,0,0];
                        found = false;
                        for ke = 1:ne2
                            ind2e(3) = ke;
                            adr2e = mk_address(ind2e);
                            [~,CA02]=get_object(sprintf('%s.CA',adr2e),'coor');
                            if isempty(CA02) % for nucleotides
                                [~,CA02]=get_object(sprintf('%s.C1''',adr2e),'coor');
                            end
                            if ~isempty(CA02)
                                CA2b = CA2b + CA02;
                                found = true;
                            end
                        end
                        if found
                            CA2b = CA2b/ne2;
                        else
                            CA2b = [];
                        end
                        CA2 = (CA2a+CA2b)/2;
                    else
                        [~,CA2a]=get_object(sprintf('%s.CA',adr1),'coor');
                        if isempty(CA2a) % for nucleotides
                            [~,CA2a]=get_object(sprintf('%s.C1''',adr1),'coor');
                        end
                        [~,CA2b]=get_object(sprintf('%s.CA',adr2),'coor');
                        if isempty(CA2b) % for nucleotides
                            [~,CA2b]=get_object(sprintf('%s.C1''',adr2),'coor');
                        end
                        CA2 = (CA2a+CA2b)/2;
                    end
            end
            if ~isempty(CA1) && ~isempty(CA2)
                r0=norm(CA1-CA2)/10; % C-alpha/C-alpha distance in nm
                if r0-0.1<rmin0, rmin0=r0-0.1; end
                if r0+0.1>rmax0, rmax0=r0+0.1; end
            else
                r0=[];
            end
            adr1 = dadr1;
            adr2 = dadr2;
            handles.pairs{poi}=[adr1 '-' adr2];
            set(handles.DEER,'Pointer','watch');
            if do_ensemble
               [rax,act_distr] = pair_distribution(handles,adr1,adr2);
               act_distr = zeros(size(act_distr));
               coupled = get(handles.checkbox_coupled_ensemble,'Value');
               if coupled % use same coordinate set numbers for both labels
                   for ke = 1:ne1
                       [~,act_distr0]=pair_distribution(handles,adr1,adr2,ke,ke);
                       if ~isempty(act_distr0)
                        act_distr = act_distr + act_distr0;
                       end
                   end
                   act_distr = act_distr/ne1;
               else % combine labels from all pairs of coordinate set numbers
                   ind1e = ind1;
                   for ke = 1:ne1
                       for ke2 = 1:ne2
                           [~,act_distr0]=pair_distribution(handles,adr1,adr2,ke,ke2);
                           if ~isempty(act_distr0)
                            act_distr = act_distr + act_distr0;
                           end
                       end
                   end
                   act_distr = act_distr/(ne1*ne2);
               end
            else
                [rax,act_distr]=pair_distribution(handles,adr1,adr2);
            end
            set(handles.DEER,'Pointer','arrow');
            h=plot(rax,act_distr,'Color',col);
            set(h,'ButtonDownFcn',{@pair_ButtonDownFcn,poi});
            if isempty(r0)
                ndis=act_distr/sum(act_distr);
                r0=sum(rax.*ndis);
            end
            handles.pair_plots{poi}=h;
            handles.pair_distributions{poi}=act_distr;
            handles.pair_r0{poi}=r0;
            h=plot([r0,r0],[0,max(act_distr)],':','Color',col);
            handles.pair_Ca_plots{poi}=h;
            if poi==1
                full_distr=act_distr;
                rmin=min(rax);
                rmax=max(rax);
            else
                full_distr=full_distr+act_distr;
            end
        end
    end
    handles.rsim=rax;
    handles.dsim=full_distr;
    ndistr=full_distr/sum(full_distr);
    rmean=sum(ndistr.*rax);
    set(handles.text_mean,'String',sprintf('<r> = %4.2f nm',rmean));
    dr=rax-rmean;
    stddev=sqrt(sum(dr.^2.*ndistr));
    set(handles.text_stddev,'String',sprintf('sr = %4.2f nm',stddev));
    tmax_min=2*((rmean+stddev)/5)^3;
    handles.tmax_min=tmax_min;
    tmax_opt=2*((rmean+stddev)/3.2)^3;
    handles.tmax_opt=tmax_opt;
    set(handles.text_minimum_tmax,'String',sprintf('%3.1f us',tmax_min));
    set(handles.text_optimum_tmax,'String',sprintf('%3.1f us',tmax_opt));
    rfilled=rax(full_distr>0.01*max(full_distr));
    rmin=min([min(rfilled) rmin0]);
    rmax=max([max(rfilled) rmax0]);
    bsl=0.1*(rmax-rmin);
    rmin=rmin-bsl;
    rmax=rmax+bsl;
    dmax=1.1*max(full_distr);
    plot(rax,full_distr,'r','LineWidth',1);
    if ~isempty(handles.tweak_distr) && ~isempty(handles.tweak_rax)
        plot(handles.tweak_rax,handles.tweak_distr,'r:','LineWidth',2);
        dmax=1.1*max([max(full_distr) max(handles.tweak_distr)]);
    end
    dmin=-0.09*dmax;
    simdistr=1;
    if ~handles.expanded
        if rmin> handles.range(1)
            handles.range(1)=rmin;
        end
        if rmax< handles.range(2)
            handles.range(2)=rmax;
        end
    end
    handles=range_update(handles);
    if handles.expanded
        axis([handles.range(1),handles.range(2),dmin,dmax]);
    else
        axis([rmin,rmax,dmin,dmax]);
    end
    if isfield(handles,'MtsslWizard_rax') && isfield(handles,'MtsslWizard_distr')
        sc=max(full_distr)/max(handles.MtsslWizard_distr);
        plot(handles.MtsslWizard_rax,sc*handles.MtsslWizard_distr,':','Color',[0 0.75 0.25],'LineWidth',1.5);
    end
    if isfield(handles,'PRONOX_rax') && isfield(handles,'PRONOX_distr')
        sc=max(full_distr)/max(handles.PRONOX_distr);
        plot(handles.PRONOX_rax,sc*handles.PRONOX_distr,':','Color',[0.25 0.0 0.75],'LineWidth',1.5);
    end
end

ff_flag=get(handles.checkbox_form_factor,'Value');

flexible = get(handles.checkbox_flexible,'Value');
moved = false;

if handles.copy_DEER
    figure(2); clf;
    set(gca,'FontSize',16);
    hold on;
else
    axes(handles.axes_DEER);
    cla;
    hold on;
end

if ~isempty(handles.texp) && ~isempty(handles.vexp)
    vsim_M=[];
    vsim_P=[];
    vsim_F=[];
    dr = [];
    if ff_flag
        set(handles.edit_zero_time,'String',sprintf('%5i',0));
        plot(handles.tdip,handles.cluster,'k');
        full=max(handles.cluster)-min(handles.cluster);
        axis([min(handles.tdip),max(handles.tdip),min(handles.cluster)-0.1*full,max(handles.cluster)+0.1*full]);
    else
        set(handles.edit_zero_time,'String',sprintf('%5i',handles.zero_time));
        plot(handles.texp,handles.vexp,'k');
        full=max(handles.vexp)-min(handles.vexp);
        axis([min(handles.texp),max(handles.texp),min(handles.vexp)-0.1*full,max(handles.vexp)+0.1*full]);
    end
    hold on;
    if simdistr
        if ff_flag
            if pairs==1
                [vsim,rmsmin,handles]=fit_formfactor(handles,rax,full_distr,handles.tdip,handles.cluster);
                if isfield(handles,'MtsslWizard_rax') && isfield(handles,'MtsslWizard_distr')
                    [vsim_M,rmsmin_M]=fit_formfactor(handles,handles.MtsslWizard_rax,handles.MtsslWizard_distr,handles.tdip,handles.cluster);
                end
                if isfield(handles,'PRONOX_rax') && isfield(handles,'PRONOX_distr'),
                    [vsim_P,rmsmin_P]=fit_formfactor(handles,handles.PRONOX_rax,handles.PRONOX_distr,handles.tdip,handles.cluster);
                end
                vsim_F = [];
                if flexible
                    set(handles.DEER,'Pointer','watch');
                    drawnow
                    [vsim_F,rmsmin_F,rax_F,distr_F,dr]=fit_formfactor_flexible(handles,handles.tdip,handles.cluster);
                    set(handles.DEER,'Pointer','arrow');
                    set(handles.edit_flex_move,'String',sprintf('%4.2f',dr/10));
                    set(handles.edit_flex_move,'ForegroundColor',handles.fit_color);
                elseif handles.flex_move ~=0
                    moved = true;
                    [vsim_D,rmsmin_D,rax_D,distr_D]=fit_formfactor_shifted(handles,handles.tdip,handles.cluster);
                end
            else
                if flexible
                    add_msg_board('Warning: Flexible fit is defined only for site pairs.');
                    set(handles.checkbox_flexible,'Value',0);
                end
                [vsim,rmsmin,handles]=fit_formfactor_multi(handles,rax,full_distr,handles.tdip,handles.cluster);
                handles.ff_multi=vsim;
                vsim_M=[];
                vsim_F = [];
                vsim_P=[];
            end
            handles.tsim=handles.tdip;
            handles.vsim=vsim;
            handles.vsim_M=vsim_M;
            handles.vsim_F=vsim_F;
            handles.vsim_P=vsim_P;
            handles.dr = dr;
        else
            [mimi,poi]=min(abs(handles.texp));
            texp=handles.texp(poi:end);
            vexp=handles.vexp(poi:end);
            if pairs==1
                [vsim,rmsmin,handles]=fit_DEER(handles,rax,full_distr,texp,vexp);
                if isfield(handles,'MtsslWizard_rax') && isfield(handles,'MtsslWizard_distr')
                    [vsim_M,rmsmin_M]=fit_DEER(handles,handles.MtsslWizard_rax,handles.MtsslWizard_distr,texp,vexp);
                end
                if isfield(handles,'PRONOX_rax') && isfield(handles,'PRONOX_distr')
                    [vsim_P,rmsmin_P]=fit_DEER(handles,handles.PRONOX_rax,handles.PRONOX_distr,texp,vexp);
                end
                if flexible
                    set(handles.DEER,'Pointer','watch');
                    drawnow
                    [vsim_F,rmsmin_F,rax_F,distr_F,dr]=fit_DEER_flexible(handles,texp,vexp);
                    set(handles.DEER,'Pointer','arrow');
                    set(handles.edit_flex_move,'String',sprintf('%4.2f',dr/10));
                    set(handles.edit_flex_move,'ForegroundColor',handles.fit_color);
                elseif handles.flex_move ~=0
                    moved = true;
                    [vsim_D,rmsmin_D,rax_D,distr_D]=fit_DEER_shifted(handles,texp,vexp);
                end
            else
                [vsim,rmsmin,handles]=fit_DEER_multi(handles,rax,full_distr,texp,vexp);
                if flexible
                    add_msg_board('Warning: Flexible fit is defined only for site pairs.');
                    set(handles.checkbox_flexible,'Value',0);
                end
            end
            handles.tsim=texp;
            handles.vsim=vsim;
            handles.vsim_M=vsim_M;
            handles.vsim_P=vsim_P;
            handles.vsim_F=vsim_F;
            handles.dr = dr;
            handles.vexp_fit=vexp;
            handles.texp_fit=texp;
        end
        handles.rmsd=rmsmin;
        if handles.copy_DEER
            figure(2);
        end
        plot(handles.tsim,handles.vsim,'r');
        if ~isempty(handles.vsim_M)
            plot(handles.tsim,handles.vsim_M,':','Color',[0 0.75 0.25],'LineWidth',1.5);
        end
        if ~isempty(handles.vsim_P)
            plot(handles.tsim,handles.vsim_P,':','Color',[0.25 0.0 0.75],'LineWidth',1.5);
        end
        if flexible
            plot(handles.tsim,vsim_F,'Color',handles.fit_color);
            set(handles.text_rmsd_flex,'String',sprintf('%8.5f',rmsmin_F));
            set(handles.text_rmsd_flex,'ForegroundColor',handles.fit_color);
%            fprintf(1,'shift: %4.1f Å, r.m.s.d.: %8.5f, r.m.s.d. flex: %8.5f\n',dr,rmsmin,rmsmin_F);
        elseif moved
            plot(handles.tsim,vsim_D,'Color',handles.flex_color);
            set(handles.text_rmsd_flex,'String',sprintf('%8.5f',rmsmin_D));
            set(handles.text_rmsd_flex,'ForegroundColor',handles.flex_color);
        else
            set(handles.text_rmsd_flex,'String','n.a.');
            set(handles.text_rmsd_flex,'ForegroundColor',[0,0,0]);
        end
        set(handles.text_rmsd,'String',sprintf('%8.5f',rmsmin));
    else
        set(handles.text_rmsd,'String','n.a.');
    end
else
    if simdistr
        if tmax_opt>8, tmax_opt=8; end
        n=round(tmax_opt/0.008);
        tmax_opt=n*0.008;
        tsim=linspace(0,tmax_opt,n+1);
        handles.tsim=tsim';
        if pairs==1
            pcf=get_std_distr_MMM(rax,full_distr,handles.Pake_r);
            ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
            ff=ff*handles.mod_depth+ones(size(ff))*(1-handles.mod_depth);
            ffdisplay=interp1(handles.Pake_t,ff,tsim,'pchip');
        else
            pcf=get_std_distr_MMM(rax,full_distr,handles.Pake_r);
            ff0=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
            ff0=ff0*handles.mod_depth+ones(size(ff0))*(1-handles.mod_depth);
            if tmax_opt>8, tmax_opt=8; end
            n=round(tmax_opt/0.008);
            tmax_opt=n*0.008;
            tsim=linspace(0,tmax_opt,n+1);
            handles.tsim=tsim';
            ffdisplay=interp1(handles.Pake_t,ff0,tsim,'pchip');
            plot(tsim,ffdisplay,'m:');
            hold on
            add_msg_board('Computing higher-order correlations');
            set(handles.DEER,'Pointer','watch');
            drawnow;
            addresses=get(handles.listbox_label,'String');
            ffdisplay=multi_formfactor(handles,addresses,tsim,handles.exp_depth);
            handles.ff_multi=ffdisplay;
            if ~isempty(ffdisplay)
                handles.ff=ffdisplay;
            end
            set(handles.DEER,'Pointer','arrow');
        end
        handles.vsim=ffdisplay;
        handles.ff=ffdisplay';
        % handles.tdip=handles.tsim;
        if handles.copy_DEER
            figure(2);
        end
        if ~isempty(handles.vsim)
            plot(handles.tsim,handles.vsim,'r');
        end
    end
end

if handles.copy_distr
    figure(1);
else
    axes(handles.axes_distribution);
end

if flexible && simdistr && ~isempty(handles.texp) && ~isempty(handles.vexp)
    plot(rax_F,distr_F,'Color',handles.fit_color);
end

if moved && simdistr && ~isempty(handles.texp) && ~isempty(handles.vexp)
    plot(rax_D,distr_D,'Color',handles.flex_color);
end

if ~isempty(handles.rexp) && ~isempty(handles.dexp)
    dexp=handles.dexp;
    if simdistr
        dexp0=interp1(handles.rexp,dexp,rax,'pchip',0);
        sc0=sum(full_distr.*dexp0)/sum(dexp0.*dexp0);
        sc=max(full_distr)/max(dexp);
        if ~isempty(handles.tweak_distr)
            sc=max(handles.tweak_distr)/max(dexp);
        end
        dexp=sc*dexp;
        handles.dexp_fit=sc0*dexp0;
    else
        rmin=50;
        rmax=0;
        dmin=0;
        dmax=0;
    end
    plot(handles.rexp,dexp,'k','LineWidth',1);
    full=max(handles.dexp)-min(handles.dexp);
    rmin=min([rmin min(handles.rexp)]);
    rmax=max([rmax max(handles.rexp)]);
    dmin=min([dmin min(dexp)-0.1*full]);
    dmax=max([dmax max(dexp)+0.1*full]);
    if handles.expanded
        axis([handles.range(1),handles.range(2),dmin,dmax]);
    else
        axis([rmin,rmax,dmin,dmax]);
    end
    hold on;
end
    
handles.copy_distr=0;
handles.copy_DEER=0;


% --- Executes on button press in checkbox_form_factor.
function checkbox_form_factor_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_form_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_form_factor

handles=update(handles);
% Update handles structure
guidata(hObject, handles);

function [handles,success]=mk_label_list(handles)
% Fills listbox_label with the list of available labels
% returns a success flag success=1 if there are at least two labels,
% success=0 otherwise

global model

success=0;
sites=0;
sitelistpoi = 0;

if isfield(handles,'trajectory')
    success=1;
    sites=length(handles.trajectory);
    trj_labels=sites;
    for k=1:length(handles.trajectory)
        label_list{k} = handles.trajectory{k};
        sitelistpoi = sitelistpoi + 1;
        site_list(sitelistpoi).type = 'trajectory';
        site_list(sitelistpoi).item = k;
        site_list(sitelistpoi).chain = 0;
        site_list(sitelistpoi).adr = handles.trajectory{k};
    end
    sites = length(handles.trajectory);
end;

bpoi = 0;
bilabels = 0;
if isfield(model,'bisites') && ~isempty(model.bisites)
    bilabels = 0;
    for k=1:length(model.bisites),
        sites = sites+length(model.bisites{k}.sites);
        bilabels = bilabels + sites;
        for r = 1:length(model.bisites{k}.sites)
            bpoi = bpoi + 1;
            adr1 = mk_address(model.bisites{k}.sites{r}.indices(1,:));
            adr2 = mk_address(model.bisites{k}.sites{r}.indices(2,:));
            handles.bilabels(bpoi).popcoor = model.bisites{k}.sites{r}.coor(:,1:4);
            handles.bilabels(bpoi).adr = sprintf('%s|%s',adr1,adr2);
            label_list{bpoi}= sprintf('%s|%s',adr1,adr2);
            sitelistpoi = sitelistpoi + 1;
            site_list(sitelistpoi).type = 'bilabel';
            site_list(sitelistpoi).item = bpoi;
            site_list(sitelistpoi).chain = model.bisites{k}.sites{r}.indices(1,2);
            site_list(sitelistpoi).adr =  sprintf('%s|%s',adr1,adr2);
            success = 1;
        end
    end;
end;
if bilabels == 0
    handles.bilabels = [];
end;

indices=resolve_address('*');

if ~isfield(model,'labels')
    labels=0;
    if isempty(indices) && ~success,
        return
    end;
end;


if isfield(model,'labels') && ~isempty(model.labels),
    sites=sites+length(model.labels);
    labels=sites;
    for k=1:length(model.labels),
        sitelistpoi = sitelistpoi + 1;
        label_list{sitelistpoi}=model.labels(k).adr;
        site_list(sitelistpoi).type = 'label';
        site_list(sitelistpoi).item = k;
        indices = resolve_address(model.labels(k).adr);
        site_list(sitelistpoi).chain = indices(2);
        site_list(sitelistpoi).adr = model.labels(k).adr;
    end;
end;

indices=resolve_address('*');

% make list of selected atoms
if ~isempty(indices),
    [m,n]=size(indices);
    poi=0;
    handles.atoms=zeros(m,5);
    for ko=1:m, % loop over all objects
        idepth=length(find(indices(ko,:)>0)); % determine type of current object
        if idepth==5 || idepth == 6,
            poi=poi+1;
            cindices=indices(ko,1:5);
            address=mk_address(cindices,true);
            label_list{poi+sites}=address;
            handles.atoms(poi,:)=cindices;
            handles.atom_adr{poi}=address;
            sitelistpoi = sitelistpoi + 1;
            site_list(sitelistpoi).type = 'atom';
            site_list(sitelistpoi).item = poi;
            site_list(sitelistpoi).chain = cindices(2);
            site_list(sitelistpoi).adr = address;
        end;
    end;
    if poi>0,
        handles.atoms=handles.atoms(1:poi,:);
        sites=sites+poi;
    else
        handles.atoms=[];
    end;
end;

if sites>=2,
    success=1;
else
    set(handles.listbox_label,'String',' ');
    set(handles.listbox_label,'Value',1);
    return;
end;

handles.site_list = site_list;

set(handles.listbox_label,'String',label_list);
set(handles.listbox_label,'Value',1);

if exist('trj_labels','var') && trj_labels>0,
    NOpos=handles.frames{1};
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name=handles.trajectory{1};
    msg{1}=name;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
elseif bilabels > 0
    NOpos = handles.bilabels(1).popcoor;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name = handles.bilabels(1).adr;
    msg{1}=name;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
elseif labels>0,
    NOpos=model.labels(1).NOpos;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name=model.labels(1).adr;
    msg{1}=name;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',x,y,z);
else
    indices=handles.atoms(1,:);
    [msg0,xyz]=get_atom(indices,'xyz');
    [msg0,pop]=get_atom(indices,'populations');
    pop=pop/sum(pop);
    xyz=pop'*xyz;
    adr=mk_address(indices,true);
    msg{1}=adr;
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] Å',xyz(1),xyz(2),xyz(3));    
end;
set(handles.text_coor,'String',msg);

function [rax,distr]=pair_distribution(handles,adr1,adr2,km1,km2)
% Distance distribution function for a pair of spin labels
% km1, km2 replace the model index, if present

global model

seppoi = strfind(adr1,'|');
if ~isempty(seppoi)
    bilabel1 = true;
    adr1a = adr1(1:seppoi-1);
    if exist('km1','var')
        ind1 = resolve_address(adr1a);
        ind1(3) = km1;
        adr1a = mk_adress(ind1);
    end
    adr1b = adr1(seppoi+1:end);
    if exist('km1','var')
        ind1 = resolve_address(adr1b);
        ind1(3) = km1;
        adr1b = mk_adress(ind1);
    end
    adr1 = sprintf('%s|%s',adr1a,adr1b);
else
    bilabel1 = false;
    if exist('km1','var')
        ind1 = resolve_address(adr1);
        ind1(3) = km1;
        adr1 = mk_address(ind1);
    end
end

seppoi = strfind(adr2,'|');
if ~isempty(seppoi)
    bilabel2 = true;
    adr2a = adr2(1:seppoi-1);
    if exist('km2','var')
        ind2 = resolve_address(adr2a);
        ind2(3) = km2;
        adr2a = mk_adress(ind2);
    end
    adr2b = adr2(seppoi+1:end);
    if exist('km2','var')
        ind2 = resolve_address(adr2b);
        ind2(3) = km2;
        adr2b = mk_adress(ind2);
    end
    adr2 = sprintf('%s|%s',adr2a,adr2b);
else
    bilabel2 = false;
    if exist('km2','var')
        ind2 = resolve_address(adr2);
        ind2(3) = km2;
        adr2 = mk_address(ind2);
    end
end

rax = [];
distr = [];

sig=0.1; % Gaussian broadening of the distance distribution

options.trajectory=[0,0];

NOpos1 = [];
NOpos2 = [];
if isfield(model,'labels'),
    for k=1:length(model.labels),
        if strcmp(model.labels(k).adr,adr1),
            NOpos1=model.labels(k).NOpos;
        end;
        if strcmp(model.labels(k).adr,adr2),
            NOpos2=model.labels(k).NOpos;
        end;
    end;
end;
if isfield(handles,'atom_adr'),
    for k=1:length(handles.atom_adr),
        if ~isempty(strfind(handles.atom_adr{k},adr1)),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            if isempty(NOpos1)
                NOpos1=[xyz pop];
            else
                NOpos1 = [NOpos1; xyz pop];
            end
        end;
        if ~isempty(strfind(handles.atom_adr{k},adr2)),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            if isempty(NOpos2)
                NOpos2=[xyz pop];
            else
                NOpos2 = [NOpos2; xyz pop];
            end
        end;
    end;
end;
if isfield(handles,'trajectory'),
    for k=1:length(handles.trajectory),
        if strcmp(handles.trajectory{k},adr1),
            NOpos1=handles.frames{k}(k,:);
            options.trajectory(1)=1;
        end;
        if strcmp(handles.trajectory{k},adr2),
            NOpos2=handles.frames{k}(k,:);
            options.trajectory(2)=1;
        end;
    end;
end;
if isfield(handles,'bilabels'),
    for k=1:length(handles.bilabels),
        if strcmp(handles.bilabels(k).adr,adr1),
            NOpos1 = handles.bilabels(k).popcoor;
        end;
        if strcmp(handles.bilabels(k).adr,adr2),
            NOpos2 = handles.bilabels(k).popcoor;
        end;
    end;
end;

if exist('NOpos1','var') && exist('NOpos2','var')
    [rax,distr]=get_distribution(NOpos1,NOpos2,sig,[],[],[],options);
end

function [NOpos1,NOpos2,CACA,options]=get_NOpos(handles,adr1,adr2)
% NO positions for a pair of label sites

global model


options.trajectory=[0,0];

if isfield(model,'labels'),
    for k=1:length(model.labels),
        if strcmp(model.labels(k).adr,adr1),
            NOpos1=model.labels(k).NOpos;
            [~,CA1]=get_object(sprintf('%s.CA',adr1),'coor');
        end;
        if strcmp(model.labels(k).adr,adr2),
            NOpos2=model.labels(k).NOpos;
            [~,CA2]=get_object(sprintf('%s.CA',adr2),'coor');
        end;
    end;
end;
if isfield(handles,'atom_adr'),
    for k=1:length(handles.atom_adr),
        if strcmp(handles.atom_adr{k},adr1),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            NOpos1=[xyz pop];
            CA1 = mean(NOpos1);
        end;
        if strcmp(handles.atom_adr{k},adr2),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            NOpos2=[xyz pop];
            CA2 = mean(NOpos2);
        end;
    end;
end;
if isfield(handles,'trajectory'),
    for k=1:length(handles.trajectory),
        if strcmp(handles.trajectory{k},adr1),
            NOpos1=handles.frames{k}(k,:);
            CA1 = mean(NOpos1);
            options.trajectory(1)=1;
        end;
        if strcmp(handles.trajectory{k},adr2),
            NOpos2=handles.frames{k}(k,:);
            CA2 = mean(NOpos2);
            options.trajectory(2)=1;
        end;
    end;
end;

CACA = CA2 - CA1;


function [rax,distr]=pair_distribution_by_pos(NOpos1,NOpos2,options)
% Distance distribution function for a pair of spin labels

sig=0.1; % Gaussian broadening of the distance distribution

[rax,distr]=get_distribution(NOpos1,NOpos2,sig,[],[],[],options);


function pair_ButtonDownFcn(hObject,eventdata,mynum)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
set(handles.text_selected_distribution,'String',handles.pairs{mynum},'FontSize',8);
for k=1:length(handles.pair_plots),
    if handles.pair_plots{k}>0,
        set(handles.pair_plots{k},'LineWidth',0.5);
        set(handles.pair_Ca_plots{k},'LineWidth',0.5);
    end;
end;
set(handles.pair_plots{mynum},'LineWidth',2);
set(handles.pair_Ca_plots{mynum},'LineWidth',2);

function [deer,rmsd,handles]=fit_DEER(handles,rax,distr,texp,vexp)

pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;
[deer,bckg,bckg_k,mod_depth,dim,rmsd]=fit_primary(handles,texp,vexp,ff1);
handles.ff=(deer./bckg);
handles.ff=handles.ff/max(handles.ff);
handles.bckg=bckg;
handles.mod_depth=mod_depth;
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.rmsd=rmsd;
set(handles.text_rmsd,'String',sprintf('%7.5f',handles.rmsd));
handles.bckg_dim=dim;
set(handles.edit_bckg,'String',sprintf('%5.3f',handles.bckg_dim));
handles.bckg_k=bckg_k;
set(handles.text_bckg_k,'String',sprintf('%5.3f',handles.bckg_k));

function [deer,rmsd,handles]=fit_DEER_multi(handles,rax,distr,texp,vexp)

addresses=get(handles.listbox_label,'String');
ff=multi_formfactor(handles,addresses,handles.Pake_t,handles.exp_depth);
handles.ff_multi=ff;
handles.tdip=texp;
if isempty(ff),
    add_msg_board('### Warning ###: Higher-order correlations neglected');
    pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
    ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
end;
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
if ~isempty(ff),
    ff2=ff1-(1-handles.exp_depth);
    ff2=ff2/max(ff2);
else
    ff2=ff1;
end;
[deer,bckg,bckg_k,mod_depth,dim,rmsd]=fit_primary(handles,texp,vexp,ff2);
handles.ff=(deer./bckg);
handles.ff=handles.ff/max(handles.ff);
handles.mod_depth=mod_depth;
handles.bckg=(1-handles.mod_depth)*bckg/max(bckg);
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.rmsd=rmsd;
set(handles.text_rmsd,'String',sprintf('%7.5f',handles.rmsd));
handles.bckg_dim=dim;
set(handles.edit_bckg,'String',sprintf('%5.3f',handles.bckg_dim));
handles.bckg_dim=dim;
handles.bckg_k=bckg_k;
set(handles.text_bckg_k,'String',sprintf('%5.3f',handles.bckg_k));

function [vsim_F,rmsmin_F,rax_F,distr_F,dr] = fit_DEER_flexible(handles,texp,vexp)

label_list = get(handles.listbox_label,'String');
adr1=label_list{handles.labels(1)};
adr2=label_list{handles.labels(2)};
[NOpos1,NOpos2,CACA,options]=get_NOpos(handles,adr1,adr2);

dr = fminsearch(@fit_shift,1,[],handles,NOpos1,NOpos2,CACA,options,texp,vexp);

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;
[vsim_F,~,~,~,~,rmsmin_F]=fit_primary(handles,texp,vexp,ff1);
rax_F = rax;
distr_F = distr;

function rmsd = fit_shift(dr,handles,NOpos1,NOpos2,CACA,options,texp,vexp)

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;
[~,~,~,~,~,rmsd]=fit_primary(handles,texp,vexp,ff1);
% fprintf(1,'dr= %4.2f Å, r.m.s.d.: %8.5f\n',dr,rmsd);

function [vsim_D,rmsmin_D,rax_D,distr_D] = fit_DEER_shifted(handles,texp,vexp)

label_list = get(handles.listbox_label,'String');
adr1=label_list{handles.labels(1)};
adr2=label_list(handles.labels{2});
[NOpos1,NOpos2,CACA,options]=get_NOpos(handles,adr1,adr2);

dr = 10*handles.flex_move;

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;
[vsim_D,~,~,~,~,rmsmin_D]=fit_primary(handles,texp,vexp,ff1);
rax_D = rax;
distr_D = distr;


function [deer,rmsd,handles]=fit_formfactor(handles,rax,distr,texp,vexp)

pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
handles.ff=ff;
deer=interp1(handles.Pake_t,ff,texp,'pchip');
deer=deer/max(deer);
handles.ff=deer;
fit_depth=get(handles.checkbox_mod_depth,'Value');
depth0=handles.mod_depth; % starting value for modulation depth
if ~isempty(handles.exp_depth),
    depth0=handles.exp_depth;
end;
if fit_depth,
    [mod_depth,rmsd]=fminsearch(@rmsd_depth,depth0,[],vexp,deer);
else
    mod_depth=depth0;
    sim=deer*mod_depth+ones(size(deer))*(1-mod_depth);
    diff=sim-vexp;
    rmsd=sqrt(diff.^2/length(diff));
end;
deer=deer*mod_depth+ones(size(deer))*(1-mod_depth);
handles.ff=deer;

% [bsl,rmsd]=fminbnd(@rms_cluster,0,100,[],deer,vexp);
% deer=deer+bsl;
% handles.mod_depth=1/(1+bsl);
% set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
% sc=sum(vexp.*vexp)/sum(vexp.*deer);
% deer=sc*deer;
handles.bckg=[];
handles.mod_depth=mod_depth;
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.rmsd=rmsd;
set(handles.text_rmsd,'String',sprintf('%7.5f',handles.rmsd));
handles.bckg_k=0;
set(handles.text_bckg_k,'String',sprintf('%5.3f',handles.bckg_k));

function [deer,rmsd,handles]=fit_formfactor_multi(handles,rax,distr,texp,vexp)

addresses=get(handles.listbox_label,'String');
deer=multi_formfactor(handles,addresses,texp,handles.exp_depth);
handles.ff_multi=deer;
deer=deer';
if isempty(deer),
    add_msg_board('### Warning ###: Higher-order correlations neglected');
    pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
    ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
    deer=interp1(handles.Pake_t,ff,texp,'pchip');
    deer=deer/max(deer);
end;
fit_depth=get(handles.checkbox_mod_depth,'Value');
depth0=handles.mod_depth; % starting value for modulation depth
if ~isempty(handles.exp_depth),
    depth0=handles.exp_depth;
end;
if fit_depth,
    [mod_depth,rmsd]=fminsearch(@rmsd_depth,depth0,[],vexp,deer);
else
    mod_depth=depth0;
    sim=deer*mod_depth+ones(size(deer))*(1-mod_depth);
    diff=sim-vexp;
    rmsd=sqrt(diff.^2/length(diff));
end;
% [bsl,rmsd]=fminbnd(@rms_cluster,0,100,[],deer,vexp);
% deer=deer+bsl;
% handles.mod_depth=handles.exp_depth/(1+bsl);
% set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
% sc=sum(vexp.*vexp)/sum(vexp.*deer);
% deer=sc*deer;
handles.bckg=[];
handles.mod_depth=mod_depth;
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.rmsd=rmsd;
set(handles.text_rmsd,'String',sprintf('%7.5f',handles.rmsd));
handles.bckg_k=0;
set(handles.text_bckg_k,'String',sprintf('%5.3f',handles.bckg_k));
deer=deer*mod_depth+ones(size(deer))*(1-mod_depth);
handles.ff=deer;


function [vsim_F,rmsmin_F,rax_F,distr_F,dr]=fit_formfactor_flexible(handles,texp,vexp)

label_list=get(handles.listbox_label,'String');
adr1=label_list{handles.labels(1)};
adr2=label_list{handles.labels(2)};
[NOpos1,NOpos2,CACA,options]=get_NOpos(handles,adr1,adr2);

fit_depth=get(handles.checkbox_mod_depth,'Value');
depth0=handles.mod_depth; % starting value for modulation depth
if ~isempty(handles.exp_depth),
    depth0=handles.exp_depth;
end;

dr = fminsearch(@fit_shift_ff,1,[],handles,NOpos1,NOpos2,CACA,options,texp,vexp,fit_depth,depth0);

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;

if fit_depth,
    [mod_depth,rmsmin_F]=fminsearch(@rmsd_depth,depth0,[],vexp,ff1);
    vsim_F=ff1*mod_depth+ones(size(ff1))*(1-mod_depth);
else
    mod_depth=depth0;
    vsim_F=ff1*mod_depth+ones(size(ff1))*(1-mod_depth);
    diff=sim-vexp;
    rmsmin_F=sqrt(diff.^2/length(diff));
end;

rax_F = rax;
distr_F = distr;

function rmsd = fit_shift_ff(dr,handles,NOpos1,NOpos2,CACA,options,texp,vexp,fit_depth,depth0)

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
if fit_depth,
    [~,rmsd]=fminsearch(@rmsd_depth,depth0,[],vexp,ff1);
else
    mod_depth=depth0;
    sim=ff1*mod_depth+ones(size(ff1))*(1-mod_depth);
    diff=sim-vexp;
    rmsd=sqrt(diff.^2/length(diff));
end;
% fprintf(1,'dr= %4.2f Å, r.m.s.d.: %8.5f\n',dr,rmsd);

function [vsim_D,rmsmin_D,rax_D,distr_D]=fit_formfactor_shifted(handles,texp,vexp)

label_list = get(handles.listbox_label,'String');
adr1=label_list{handles.labels(1)};
adr2=label_list{handles.labels(2)};
[NOpos1,NOpos2,CACA,options]=get_NOpos(handles,adr1,adr2);

fit_depth=get(handles.checkbox_mod_depth,'Value');
depth0=handles.mod_depth; % starting value for modulation depth
if ~isempty(handles.exp_depth),
    depth0=handles.exp_depth;
end;

dr = 10*handles.flex_move;

rCACA = norm(CACA);
shift = CACA*dr/(2*rCACA);
[m1,~] = size(NOpos1);
[m2,~] = size(NOpos2);
[rax,distr]=pair_distribution_by_pos(NOpos1-reprowvector([shift 0],m1),NOpos2+reprowvector([shift 0],m2),options);
pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
ff1=interp1(handles.Pake_t,ff,texp,'pchip');
ff1=ff1/max(ff1);
handles.ff=ff1;
handles.tsim=texp;

if fit_depth,
    [mod_depth,rmsmin_D]=fminsearch(@rmsd_depth,depth0,[],vexp,ff1);
    vsim_D=ff1*mod_depth+ones(size(ff1))*(1-mod_depth);
else
    mod_depth=depth0;
    vsim_D=ff1*mod_depth+ones(size(ff1))*(1-mod_depth);
    diff=sim-vexp;
    rmsmin_D=sqrt(diff.^2/length(diff));
end;

rax_D = rax;
distr_D = distr;

% --- Executes on button press in pushbutton_detach_DEER.
function pushbutton_detach_DEER_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_DEER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy_DEER=1;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'deer_window.html');
webcall(entry,'-helpbrowser');

function ff = multi_formfactor(handles,addresses,t,depth)
% Distance distribution function for a pair of spin labels

global model

if ~isempty(handles.ff_multi), % avoid computing the same thing twice
    ff=handles.ff_multi;
    return
end;

if isfield(model,'labels'),
    for k=1:length(model.labels),
        for kk=1:length(addresses),
            if strcmp(model.labels(k).adr,addresses{kk}),
                NOpos{kk}=model.labels(k).NOpos;
            end;
        end;
    end;
end;
if isfield(handles,'atom_adr'),
    for k=1:length(handles.atom_adr),
        for kk=1:length(addresses),
            if strcmp(handles.atom_adr{k},addresses{kk}),
                indices=handles.atoms(k,:);
                [msg0,xyz]=get_atom(indices,'xyz');
                [msg0,pop]=get_atom(indices,'populations');
                NOpos{kk}=[xyz pop];
            end;
        end;
    end;
end;

ff=get_ff_higher_order(NOpos,t,depth,'medium');


% --- Executes on button press in pushbutton_any_rotamers.
function pushbutton_any_rotamers_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_any_rotamers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

numlab=length(handles.labels);

if numlab~=2
    add_msg_board('Error: "Any rotamers?" requires exactly two labels.');
    return
else
    label_list = get(handles.listbox_label,'String');
    adr1=label_list{handles.labels(1)};
    adr2=label_list{handles.labels(2)};
    set(handles.DEER,'Pointer','watch');
    [rax,distr]=tweak_rotamer_populations(adr1,adr2,handles.rexp,handles.dexp);
    handles.tweak_rax=rax;
    handles.tweak_distr=distr;
    set(handles.DEER,'Pointer','arrow');
    handles=update(handles);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_range_analysis.
function pushbutton_range_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_range_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

numlab=length(handles.labels);

if numlab~=2,
    add_msg_board('Error: Range analysis requires exactly two labels.');
    return
else
    label_list = get(handles.listbox_label,'String');
    adr1=label_list{handles.labels(1)};
    adr2=label_list{handles.labels(2)};
    set(handles.DEER,'Pointer','watch');
    [pairs,fname]=analyze_range(handles,adr1,adr2);
    set(handles.DEER,'Pointer','arrow');
    if ~isempty(fname),
        hMain.report_file=fname;
        report_editor;
    end;
end;
% Update handles structure
guidata(hObject, handles);

function edit_lower_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower as text
%        str2double(get(hObject,'String')) returns contents of edit_lower as a double

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles]=edit_update_MMM(handles,hObject,min(handles.rsim),handles.range(2)-0.05,handles.range(1),'%4.2f',0);
    handles.range=[v,handles.range(2)];
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes during object creation, after setting all properties.
function edit_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_lower_left.
function pushbutton_lower_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lower_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)-0.05>=min(handles.rsim)
    handles.range(1)=handles.range(1)-0.05;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;


% --- Executes on button press in pushbutton_lower_right.
function pushbutton_lower_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lower_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)+0.05<=handles.range(2)-0.05,
    handles.range(1)=handles.range(1)+0.05;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;


function edit_upper_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upper as text
%        str2double(get(hObject,'String')) returns contents of edit_upper as a double

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles]=edit_update_MMM(handles,hObject,handles.range(1)+0.05,max(handles.rsim),handles.range(2),'%4.2f',0);
    handles.range=[handles.range(1),v];
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes during object creation, after setting all properties.
function edit_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_upper_left.
function pushbutton_upper_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_upper_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)-0.05>=handles.range(1)+0.05
    handles.range(2)=handles.range(2)-0.05;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes on button press in pushbutton_upper_right.
function pushbutton_upper_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_upper_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)+0.05<=max(handles.rsim)
    handles.range(2)=handles.range(2)+0.05;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

function handles=range_update(handles)

if isfield(handles,'left_crsr') && ~ishandle(handles.left_crsr)
    set(handles.left_crsr,'XData',[handles.range(1) handles.range(1)]);
else
    axes(handles.axes_distribution);
    handles.left_crsr=plot([handles.range(1) handles.range(1)],[min(handles.dsim),max(handles.dsim)],'b','LineWidth',0.5);
end
set(handles.edit_lower,'String',sprintf('%4.2f',handles.range(1)));

if isfield(handles,'right_crsr') && ~ishandle(handles.right_crsr)
    set(handles.right_crsr,'XData',[handles.range(2) handles.range(2)]);
else
    axes(handles.axes_distribution);
    handles.right_crsr=plot([handles.range(2) handles.range(2)],[min(handles.dsim),max(handles.dsim)],'m','LineWidth',0.5);
end
set(handles.edit_upper,'String',sprintf('%4.2f',handles.range(2)));


function [pairs,fname]=analyze_range(handles,adr1,adr2)
% Analyze which rotamer pairs contribute to distance distribution in the
% range speciefied by handles.range

global model
global general

pairs=[];
fname='';

rotamers1=0;
rotamers2=0;

library1='';
library2='';
if isfield(model,'labels'),
    for k=1:length(model.labels),
        if strcmp(model.labels(k).adr,adr1),
            NOpos1=model.labels(k).NOpos;
            [significant1,n1]=size(NOpos1);
            rotamers1=true;
            [library1,NOpos1]=find_library_NOpos(adr1);
        end;
        if strcmp(model.labels(k).adr,adr2),
            NOpos2=model.labels(k).NOpos;
            [significant2,n2]=size(NOpos2);
            rotamers2=true;
            [library2,NOpos2]=find_library_NOpos(adr2);
        end;
    end;
end;
if ~isempty(library1),
    char_table1=rotamer_char_table(library1);
end;
if ~isempty(library2),
    char_table2=rotamer_char_table(library2);
end;
if isfield(handles,'atom_adr'),
    for k=1:length(handles.atom_adr),
        if strcmp(handles.atom_adr{k},adr1),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            indy=1:length(pop);
            NOpos1=[xyz pop indy'];
            rotamers1=false;
        end;
        if strcmp(handles.atom_adr{k},adr2),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            indy=1:length(pop);
            NOpos2=[xyz pop indy'];
            rotamers2=false;
        end;
    end;
end;

if ~rotamers1 && ~rotamers2,
    add_msg_board('Error: For range analysis, at least one site must be a spin label with rotamers');
    return
end;

[pairs,popsum]=pairs_in_range(handles.range,NOpos1,NOpos2); % range must be supplied in nm

if isempty(pairs),
    add_msg_board('No rotamer pairs contribute to distance distribution in selected range.');
else
    frac=sum(pairs(:,3))/popsum;
    add_msg_board(sprintf('Selected range corresponds to %6.1f%% of total pair population',100*frac));
    [sortpop,poppoi]=sort(pairs(:,3),'descend');
    sortpop=sortpop/sum(sortpop);
    [mi,significant]=min(abs(cumsum(sortpop)-0.99));
    add_msg_board(sprintf('%i pairs contribute in this range, %i of those contribute 99%% population of the range',length(poppoi),significant));
    spairs=pairs(poppoi(1:significant),:); % leading 99% of pairs sorted by descending population
    poi1=0;
    [m1,n1]=size(NOpos1);
    rot1=zeros(m1,2);
    for k=1:m1, % analyze site 1
        cpairs=find(spairs(:,1)==NOpos1(k,5));
        if ~isempty(cpairs),
            pops=sum(spairs(cpairs,3));
            poi1=poi1+1;
            rot1(poi1,1)=NOpos1(k,5);
            rot1(poi1,2)=pops;
        end;
    end;
    rot1=rot1(1:poi1,:);
    add_msg_board(sprintf('At site 1: %i out of %i rotamers contribute in selected range.',poi1,m1));
    poi2=0;
    [m2,n2]=size(NOpos2);
    rot2=zeros(m2,2);
    for k=1:m2, % analyze site 2
        cpairs=find(spairs(:,2)==NOpos2(k,5));
        if ~isempty(cpairs),
            pops=sum(spairs(cpairs,3));
            poi2=poi2+1;
            rot2(poi2,1)=NOpos2(k,5);
            rot2(poi2,2)=pops;
        end;
    end;
    rot2=rot2(1:poi2,:);
    add_msg_board(sprintf('At site 2: %i out of %i rotamers contribute in selected range.',poi2,m2));
    fname=[general.tmp_files 'chi1_chi2_analysis.dat'];
    fid=fopen(fname,'w');
    fprintf(fid,'--- MMM chi1/chi2 analysis of rotamers ---\nSite 1: %s, Site 2: %s\ncontributions to the distance range [%4.2f, %4.2f] nm\n\n',adr1,adr2,handles.range(1),handles.range(2));
    fprintf(fid,'Selected range corresponds to %6.1f%% of total pair population\n',100*frac);
    fprintf(fid,'%i pairs contribute in this range,\n%i of those contribute 99%% population of the range\n\n',length(poppoi),significant);
    fprintf(fid,'At %s: %i out of %i rotamers contribute in selected range.\n',adr1,poi1,m1);
    if ~isempty(library1),
        fprintf(fid,'The chi1/chi2 analysis of rotamers at %s follows:\n\n',adr1);
        fprintf(fid,'rotamer   total population (%%)   pop. in range (%%)\n');
        [pops,chi1_chi2]=chi1_chi2_analysis(NOpos1,rot1,char_table1);
        [m,n]=size(pops);
        for k=1:m,
            fprintf(fid,'(%s) %15.1f%20.1f\n',chi1_chi2{k},100*pops(k,1)/sum(pops(:,1)),100*pops(k,2)/sum(pops(:,2)));
        end;
    end;
    fprintf(fid,'\nAt %s: %i out of %i rotamers contribute in selected range.\n',adr2,poi2,m2);
    if ~isempty(library2),
        [pops,chi1_chi2]=chi1_chi2_analysis(NOpos2,rot2,char_table2);
        fprintf(fid,'The chi1/chi2 analysis of rotamers at %s follows:\n\n',adr2);
        fprintf(fid,'rotamer   total population (%%)   pop. in range (%%)\n');
        [m,n]=size(pops);
        for k=1:m,
            fprintf(fid,'(%s) %15.1f%20.1f\n',chi1_chi2{k},100*pops(k,1)/sum(pops(:,1)),100*pops(k,2)/sum(pops(:,2)));
        end;
    end;
    fclose(fid);
    % Now color rotamers in range
    rotlist=NOpos1(:,5);
    poplist=NOpos1(:,4);
    [sortpop,poppoi]=sort(poplist,'descend');
    attached1=rotlist(poppoi(1:significant1)); % list of numbers of significant rotamers
    for k=1:length(attached1),
        rotamer=attached1(k);
        in_range=sum(find(spairs(:,1)==rotamer));
        if in_range,
            loctag=get_loctag(k);
            locadr=sprintf('%s.:%s',adr1,loctag);
            set_object(locadr,'color',{[0,0,0]});
        end;
    end;
    rotlist=NOpos2(:,5);
    poplist=NOpos2(:,4);
    [sortpop,poppoi]=sort(poplist,'descend');
    attached2=rotlist(poppoi(1:significant2)); % list of numbers of significant rotamers
    for k=1:length(attached2),
        rotamer=attached2(k);
        in_range=sum(find(spairs(:,2)==rotamer));
        if in_range,
            loctag=get_loctag(k);
            locadr=sprintf('%s.:%s',adr2,loctag);
            set_object(locadr,'color',{[0,0,0]});
        end;
    end;
end;

function [library,NOpos]=find_library_NOpos(adr)
% returns the rotamer library for a labeled residue at a given address
% empty string is returned, if residue is not found or library was not
% stored

global model

library='';
NOpos=[];
indices=resolve_address(adr);
if isempty(indices) || length(indices)~=4,
    return;
end;


nscan=length(model.sites);
for k=1:nscan,
    if isfield(model.sites{k},'library'),
        currlib=model.sites{k}.library;
    else
        currlib='';
    end;    
    for kc=1:length(model.sites{k}),
        for kk=1:length(model.sites{k}(kc).residue),
            diff=sum(abs(indices-model.sites{k}(kc).residue(kk).indices));
            if diff==0,
                library=currlib;
                NOpos=model.sites{k}(kc).residue(kk).NOpos;
            end;
        end;
    end;
end;

function [pops,chi1_chi2]=chi1_chi2_analysis(NOpos,rot,char_table)
% Populations of chi1-chi2 rotamer groups in full distribution NOpos and
% range-selected distribution rot, the character table char_table for the
% corresponding rotamer library must be provided

chi1_chi2{1}='m,m';
chi1_chi2{2}='m,t';
chi1_chi2{3}='m,p';
chi1_chi2{4}='t,m';
chi1_chi2{5}='t,t';
chi1_chi2{6}='t,p';
chi1_chi2{7}='p,m';
chi1_chi2{8}='p,t';
chi1_chi2{9}='p,p';

pops=zeros(9,2);

[m,n]=size(NOpos);

for k=1:m,
    indy=NOpos(k,5);
    character=char_table{indy};
    short_char=character(1:3);
    for kk=1:9,
        if strcmp(short_char,chi1_chi2{kk}),
            pops(kk,1)=pops(kk,1)+NOpos(k,4);
        end;
    end;
end;

[mr,nr]=size(rot);

for k=1:mr,
    indy=rot(k,1);
    character=char_table{indy};
    short_char=character(1:3);
    for kk=1:9,
        if strcmp(short_char,chi1_chi2{kk}),
            pops(kk,2)=pops(kk,2)+rot(k,2);
        end;
    end;
end;

function loctag=get_loctag(rotamer)
% return the location tag for a numbered rotamer

loctags0=':A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:';

locid=mod(rotamer,26);
if locid==0, locid=26; end;
locmid=floor((rotamer-1)/26);
if locmid>0,
    loctag=strcat(id2tag(locmid,loctags0),id2tag(locid,loctags0));
else
    loctag=id2tag(locid,loctags0);
end;


% --- Executes on button press in togglebutton_expand.
function togglebutton_expand_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state=get(hObject,'Value');
handles.expanded=state;
if state,
    set(hObject,'TooltipString','Return to full view');
    set(hObject,'String','Full view');
else
    set(hObject,'TooltipString','Expand or shrink view to selected range');
    set(hObject,'String','Expand');
end;

handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_uniform.
function pushbutton_uniform_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_uniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

numlab=length(handles.labels);

if numlab~=2
    add_msg_board('Error: "Uniform!" requires exactly two labels.');
    return
else
    label_list = get(handles.listbox_label,'String');
    adr1=label_list{handles.labels(1)};
    adr2=label_list{handles.labels(2)};
    set(handles.DEER,'Pointer','watch');
    [rax,distr]=uniform_rotamers(adr1,adr2,handles);
    handles.tweak_rax=rax;
    handles.tweak_distr=distr;
    set(handles.DEER,'Pointer','arrow');
    handles=update(handles);
end;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_PRONOX.
function pushbutton_PRONOX_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PRONOX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('.add','Load PRONOX conformation file');
fname = fullfile(pathname, filename);
set(handles.DEER,'Pointer','watch');
drawnow
[rax,distr]=rd_PRONOX(fname);
set(handles.DEER,'Pointer','arrow');
handles.PRONOX_rax=rax;
handles.PRONOX_distr=distr;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in pushbutton_MtsslWizard.
function pushbutton_MtsslWizard_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MtsslWizard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('.txt','Load MtsslWizard distance distribution file');
fname = fullfile(pathname, filename);
set(handles.DEER,'Pointer','watch');
drawnow
[rax,distr]=rd_MtsslWizard(fname);
set(handles.DEER,'Pointer','arrow');
handles.MtsslWizard_rax=rax;
handles.MtsslWizard_distr=distr;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_call_PRONOX.
function pushbutton_call_PRONOX_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_call_PRONOX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global web_adr

url=web_adr.PRONOX;
webcall(url);


% --- Executes on button press in checkbox_flexible.
function checkbox_flexible_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_flexible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_flexible

if ~get(hObject,'Value'),
    set(handles.edit_flex_move,'String',sprintf('%5.2f',handles.flex_move));
    set(handles.edit_flex_move,'ForegroundColor',handles.flex_color);
end;
handles=update(handles);
guidata(hObject,handles);



function edit_flex_move_Callback(hObject, eventdata, handles)
% hObject    handle to edit_flex_move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_flex_move as text
%        str2double(get(hObject,'String')) returns contents of edit_flex_move as a double

[v,handles]=edit_update_MMM(handles,hObject,-4.0,4.0,0.0,'%5.2f',0);
handles.flex_move = v;
set(handles.edit_flex_move,'ForegroundColor',handles.flex_color);
handles=update(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_flex_move_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_flex_move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ensemble.
function checkbox_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ensemble

if get(hObject,'Value'),
    set(handles.checkbox_coupled_ensemble,'Enable','on');
else
    set(handles.checkbox_coupled_ensemble,'Enable','off');
end;
handles=update(handles);
guidata(hObject,handles);


% --- Executes on button press in checkbox_coupled_ensemble.
function checkbox_coupled_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_coupled_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_coupled_ensemble

handles=update(handles);
guidata(hObject,handles);
