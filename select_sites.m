function varargout = select_sites(varargin)
% SELECT_SITES M-file for select_sites.fig
%      SELECT_SITES, by itself, creates a new SELECT_SITES or raises the existing
%      singleton*.
%
%      H = SELECT_SITES returns the handle to a new SELECT_SITES or the handle to
%      the existing singleton*.
%
%      SELECT_SITES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_SITES.M with the given input arguments.
%
%      SELECT_SITES('Property','Value',...) creates a new SELECT_SITES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_sites_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_sites_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_sites

% Last Modified by GUIDE v2.5 19-Oct-2012 13:51:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_sites_OpeningFcn, ...
                   'gui_OutputFcn',  @select_sites_OutputFcn, ...
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


% --- Executes just before select_sites is made visible.
function select_sites_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_sites (see VARARGIN)

% Choose default command line output for select_sites
handles.output = hObject;

% global MMM_icon
global hMain
global model

if ~isfield(model,'sites') || isempty(model.sites)
    add_msg_board('ERROR: No spin labeling site scan was performed yet.');
    delete(hObject);
    return;
end;

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];

set(handles.text_total_selected,'String','0');

handles.min_rotamers=10;

handles.num_pairs=20;
handles.basis_size=handles.num_pairs;

handles.max_dist=5;
handles.threshold=0.35;
handles.tolerance=0.1;
handles.constraint_sites=[];
handles.constraints=[];

set(handles.edit_min_rotamers,'String',sprintf('%i',handles.min_rotamers));
set(handles.edit_pairs,'String',sprintf('%i',handles.num_pairs));
set(handles.edit_basis_size,'String',sprintf('%i',handles.basis_size));
set(handles.edit_basis_size,'Enable','off');
set(handles.edit_max_dist,'String',sprintf('%4.1f',handles.max_dist));
set(handles.edit_threshold,'String',sprintf('%4.2f',handles.threshold));
set(handles.edit_constraint_tolerance,'String',sprintf('%4.2f',handles.tolerance));
set(handles.checkbox_matched_basis,'Value',1);

handles.scans=length(model.sites);
scan_list{1}='1';

if handles.scans>1,
    for k=2:handles.scans,
        scan_list{k}=sprintf('%i',k);
    end;
    set(handles.popupmenu_site_scan,'Enable','on');
else
    set(handles.popupmenu_site_scan,'Enable','off');
end;
set(handles.popupmenu_site_scan,'String',scan_list);

handles=setup_listboxes(handles);

handles.selection_list=zeros(0,3);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_sites wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_sites_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in popupmenu_site_scan.
function popupmenu_site_scan_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_site_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_site_scan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_site_scan

handles=setup_listboxes(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_site_scan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_site_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_selected.
function listbox_selected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_selected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selected

handles=pick_site_selected(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_unselect.
function pushbutton_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel=get(handles.listbox_selected,'Value');
list=get(handles.listbox_selected,'String');
selections=get(handles.listbox_selected,'UserData');
indy=selections(sel,:);
adr=list{sel};
num=length(list);
if sel<num,
    for k=sel:num-1,
        list{k}=list{k+1};
        selections(k,:)=selections(k+1,:);
    end;
end;
list=list(1:num-1);
selections=selections(1:num-1,:);
set(handles.listbox_selected,'Value',1);
set(handles.listbox_selected,'String',list);
set(handles.listbox_selected,'UserData',selections);
set(handles.text_selected,'String',sprintf('%i',num-1));

list=get(handles.listbox_unselected,'String');
unselections=get(handles.listbox_unselected,'UserData');
[num,n]=size(unselections);
unselections(num+1,:)=indy;
num=length(list);
list{num+1}=adr;
set(handles.listbox_unselected,'Value',num+1);
set(handles.listbox_unselected,'String',list);
set(handles.listbox_unselected,'UserData',unselections);
set(handles.text_unselected,'String',sprintf('%i',num+1));

handles=pick_site_selected(handles);
handles=pick_site_unselected(handles);

guidata(hObject,handles);

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel=get(handles.listbox_unselected,'Value');
list=get(handles.listbox_unselected,'String');
unselections=get(handles.listbox_unselected,'UserData');
indy=unselections(sel,:);
adr=list{sel};
num=length(list);
if sel<num,
    for k=sel:num-1,
        list{k}=list{k+1};
        unselections(k,:)=unselections(k+1,:);
    end;
end;
list=list(1:num-1);
unselections=unselections(1:num-1,:);
set(handles.listbox_unselected,'Value',1);
set(handles.listbox_unselected,'String',list);
set(handles.listbox_unselected,'UserData',unselections);
set(handles.text_unselected,'String',sprintf('%i',num-1));

list=get(handles.listbox_selected,'String');
selections=get(handles.listbox_selected,'UserData');
[num,n]=size(selections);
selections(num+1,:)=indy;
num=length(list);
list{num+1}=adr;
set(handles.listbox_selected,'Value',num+1);
set(handles.listbox_selected,'String',list);
set(handles.listbox_selected,'UserData',selections);
set(handles.text_selected,'String',sprintf('%i',num+1));

handles=pick_site_selected(handles);
handles=pick_site_unselected(handles);

guidata(hObject,handles);

% --- Executes on selection change in listbox_unselected.
function listbox_unselected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_unselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_unselected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_unselected

handles=pick_site_unselected(handles);

% --- Executes during object creation, after setting all properties.
function listbox_unselected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_unselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_rotamers_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_rotamers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_rotamers as text
%        str2double(get(hObject,'String')) returns contents of edit_min_rotamers as a double

[v,handles]=edit_update_MMM(handles,hObject,1,100,10,'%i',1);
handles.min_rotamers=v;
handles=setup_listboxes(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_min_rotamers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_rotamers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'select_sites.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selections=get(handles.listbox_selected,'UserData');
scan=get(handles.popupmenu_site_scan,'Value');

[m,n]=size(handles.selection_list);
[mm,n]=size(selections);

new_selections=zeros(size(selections));
poi=0;
for k=1:mm,
    new_item=true;
    if m>0,
        for kk=1:m,
            if handles.selection_list(kk,1)==scan && handles.selection_list(kk,2)==selections(k,1) && handles.selection_list(kk,3)==selections(k,2),
                new_item=false;
            end;
        end;
    end;
    if new_item,
        poi=poi+1;
        new_selections(poi,:)=selections(k,:);
    end;
end;
new_selection_list=zeros(m+poi,3);
if m>0,
    new_selection_list(1:m,:)=handles.selection_list;
end;
if poi>0,
    for k=m+1:m+poi,
        new_selection_list(k,1)=scan;
        new_selection_list(k,2:3)=new_selections(k-m,:);
    end;
    handles.selection_list=new_selection_list;
    set(handles.text_total_selected,'String',sprintf('%i',m+poi));
    guidata(hObject,handles);
end;

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

model.sites_selected=handles.selection_list;
delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button=questdlg('Do you really want to discard all selections?','All selections will be lost.','No','Yes','No');
if strcmp(button,'Yes'),
    delete(handles.figure1);
end;

function handles=setup_listboxes(handles,pick_sel,pick_unsel)

global model

if nargin<2,
    pick_sel=1;
end;
if nargin<3,
    pick_unsel=1;
end;

selected_list{1}='';
selection=zeros(1000,2);
selpoi=0;
unselected_list{1}='';
unselection=zeros(1000,2);
unselpoi=0;

scan=get(handles.popupmenu_site_scan,'Value');

sites=model.sites{scan};

n=length(sites);

for k=1:n,
    residues=sites(k).residue;
    m=length(residues);
    for kk=1:m,
        cindices=residues(kk).indices;
        address=mk_address(cindices,1);
        rotamers=residues(kk).rotamers;
        numrot=length(rotamers);
        if numrot>=handles.min_rotamers,
            selpoi=selpoi+1;
            selection(selpoi,1)=k;
            selection(selpoi,2)=kk;
            selected_list{selpoi}=address;
        else
            unselpoi=unselpoi+1;
            unselection(unselpoi,1)=k;
            unselection(unselpoi,2)=kk;
            unselected_list{unselpoi}=address;
        end;
    end;
end;

selection=selection(1:selpoi,:);
unselection=unselection(1:unselpoi,:);
set(handles.listbox_selected,'Value',pick_sel);
set(handles.listbox_selected,'String',selected_list);
set(handles.listbox_selected,'UserData',selection);
set(handles.listbox_unselected,'Value',pick_unsel);
set(handles.listbox_unselected,'String',unselected_list);
set(handles.listbox_unselected,'UserData',unselection);
set(handles.text_selected,'String',sprintf('%i',selpoi));
set(handles.text_unselected,'String',sprintf('%i',unselpoi));

handles=pick_site_selected(handles);
handles=pick_site_unselected(handles);

function handles=pick_site_selected(handles)

global model

partition_function_threshold=0.05;
loose_threshold=0.50;

% check if there is binding site and metal center information
if isfield(model,'keywords'),
    mc=tag2id('metal centers',model.keywords);
    bs=tag2id('binding sites',model.keywords);
else
    mc=[];
    bs=[];
end;

pick=get(handles.listbox_selected,'Value');
selection=get(handles.listbox_selected,'UserData');
[m,n]=size(selection);
if m<1,
    set(handles.text_picked_selected,'String','No sites available');
    set(handles.text_picked_selected_rotamers,'String','n.a.');
    set(handles.text_picked_selected_rmsd,'String','n.a.');
    set(handles.text_picked_selected_binding,'String','No specific binding reported.');
    set(handles.text_picked_selected_binding,'ForegroundColor','k');
    set(handles.text_selected_type,'String',' ');
    set(handles.text_selected_type,'ForegroundColor','k');
    return;
end;
sel=selection(pick,:);
scan=get(handles.popupmenu_site_scan,'Value');
sites=model.sites{scan};
residue=sites(sel(1)).residue(sel(2));
set(handles.text_picked_selected,'String',sprintf('Picked site information (%s)',residue.label));
numrot=length(residue.rotamers);
set(handles.text_picked_selected_rotamers,'String',sprintf('%i',numrot));
rmsd=NOpos_rmsd(residue.NOpos);
set(handles.text_picked_selected_rmsd,'String',sprintf('%4.2f',rmsd));
partition_function=residue.partition_function;

type=1;
if partition_function<partition_function_threshold,
    type=3;
end;
if type==1 && numrot<handles.min_rotamers,
    type=2;
end;
if type==1 && rmsd>loose_threshold,
    type=4;
end;

binding=false;
metal=false;
if ~isempty(bs) || ~isempty(mc), % test for binding site or metal center
    rindices=residue.indices;
    [msg,anno]=get_annotations(rindices);
    if ~isempty(bs) && ~isempty(anno),
        if ~isempty(find(anno.keywords==bs, 1)),
            binding=true;
        end;
    end;
    if ~isempty(mc),
        radr=mk_address(rindices);
        [msg,aindices]=get_object(radr,'children');
        [maa,naa]=size(aindices);
        for kaa=1:maa,
            [msg,anno]=get_annotations(aindices(kaa,:));
            if ~isempty(anno),
                if ~isempty(find(anno.keywords==mc, 1)),
                    aadr=mk_address(aindices(kaa,:));
                    paa=findstr(aadr,'.');
                    if ~isempty(paa) && paa<length(aadr),
                        metal=true;
                    end;
                end;
            end;
        end;
    end;
end;

if metal,
    set(handles.text_picked_selected_binding,'String','Metal binding residue.');
    set(handles.text_picked_selected_binding,'ForegroundColor','r');
elseif binding
    set(handles.text_picked_selected_binding,'String','Residue involved in binding site.');
    set(handles.text_picked_selected_binding,'ForegroundColor','m');
else
    set(handles.text_picked_selected_binding,'String','No specific binding reported.');
    set(handles.text_picked_selected_binding,'ForegroundColor','k');
end;

switch type
    case 1,
        msg='Well defined accessible site';
        col='g';
    case 2,
        msg='Restricted site, prediction may be poor.';
        col='y';
    case 3,
        msg='Tight site, labeling may fail.';
        col='r';
    case 4,
        msg='Loose site, may be imprecise.';
        col='y';
end;
set(handles.text_selected_type,'String',msg);
set(handles.text_selected_type,'ForegroundColor',col);

function handles=pick_site_unselected(handles)

global model

partition_function_threshold=0.05;
loose_threshold=0.50;

% check if there is binding site and metal center information
if isfield(model,'keywords'),
    mc=tag2id('metal centers',model.keywords);
    bs=tag2id('binding sites',model.keywords);
else
    mc=[];
    bs=[];
end;

pick=get(handles.listbox_unselected,'Value');
unselection=get(handles.listbox_unselected,'UserData');
[m,n]=size(unselection);
if m<1,
    set(handles.text_picked_unselected,'String','No sites available');
    set(handles.text_picked_unselected_rotamers,'String','n.a.');
    set(handles.text_picked_unselected_rmsd,'String','n.a.');
    set(handles.text_picked_unselected_binding,'String','No specific binding reported.');
    set(handles.text_picked_unselected_binding,'ForegroundColor','k');
    set(handles.text_unselected_type,'String',' ');
    set(handles.text_unselected_type,'ForegroundColor','k');
    return;
end;
sel=unselection(pick,:);
scan=get(handles.popupmenu_site_scan,'Value');
sites=model.sites{scan};
residue=sites(sel(1)).residue(sel(2));
set(handles.text_picked_unselected,'String',sprintf('Picked site information (%s)',residue.label));
numrot=length(residue.rotamers);
set(handles.text_picked_unselected_rotamers,'String',sprintf('%i',numrot));
rmsd=NOpos_rmsd(residue.NOpos);
set(handles.text_picked_unselected_rmsd,'String',sprintf('%4.2f',rmsd));
partition_function=residue.partition_function;

type=1;
if partition_function<partition_function_threshold,
    type=3;
end;
if type==1 && numrot<handles.min_rotamers,
    type=2;
end;
if type==1 && rmsd>loose_threshold,
    type=4;
end;

binding=false;
metal=false;
if ~isempty(bs) || ~isempty(mc), % test for binding site or metal center
    rindices=residue.indices;
    [msg,anno]=get_annotations(rindices);
    if ~isempty(bs) && ~isempty(anno),
        if ~isempty(find(anno.keywords==bs, 1)),
            binding=true;
        end;
    end;
    if ~isempty(mc),
        radr=mk_address(rindices);
        [msg,aindices]=get_object(radr,'children');
        [maa,naa]=size(aindices);
        for kaa=1:maa,
            [msg,anno]=get_annotations(aindices(kaa,:));
            if ~isempty(anno),
                if ~isempty(find(anno.keywords==mc, 1)),
                    aadr=mk_address(aindices(kaa,:));
                    paa=findstr(aadr,'.');
                    if ~isempty(paa) && paa<length(aadr),
                        metal=true;
                    end;
                end;
            end;
        end;
    end;
end;

if metal,
    set(handles.text_picked_unselected_binding,'String','Metal binding residue.');
    set(handles.text_picked_unselected_binding,'ForegroundColor','r');
elseif binding
    set(handles.text_picked_unselected_binding,'String','Residue involved in binding site.');
    set(handles.text_picked_unselected_binding,'ForegroundColor','m');
else
    set(handles.text_picked_unselected_binding,'String','No specific binding reported.');
    set(handles.text_picked_unselected_binding,'ForegroundColor','k');
end;

switch type
    case 1,
        msg='Well defined accessible site';
        col='g';
    case 2,
        msg='Restricted site, prediction may be poor.';
        col='y';
    case 3,
        msg='Tight site, labeling may fail.';
        col='r';
    case 4,
        msg='Loose site, may be imprecise.';
        col='y';
end;
set(handles.text_unselected_type,'String',msg);
set(handles.text_unselected_type,'ForegroundColor',col);

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


% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.selection_list=zeros(0,3);
set(handles.text_total_selected,'String','0');

% Update handles structure
guidata(hObject, handles);

function [shortlist,score,redundant,pairs,cutoff]=label_pair_score(maxn,site_list,CW_flag,corr_cut,max_dist,basis_size)
% function [shortlist,score,redundant]=label_pair_score(maxn,site_list,CW_flag,corr_cut,max_dist)
%
% determines up to maxn label pairs that are predicted to be most useful
% for characterizing state transitions via low-frequency modes of an
% elastic network model
%
% maxn      maximum number of residue pairs requested
% site_list list of label sites from which pairs can be created
% CW_flag   1: pair distances in CW EPR range are allowed, 0: they are not
%           allowed
% corr_cut  optional cutoff value for correlation between radical pairs
%           that are considered redundant
% max_dist  maximum pair distance allowed [nm]
% basis_size number of slow normal modes to be used, values smaller than
%            one are interpreted as full basis size (mode 7...N) 
% shortlist a short list of non-redundant residue pairs sorted by
%           descending score, 1st column: score, 2nd/3rd column residue
%           numbers, 4th column mean label-label distance, 5th column
%           standard deviation of label-label distance
% score     relative expected mean square distance changes for propagation
%           along low-frequency modes, sorted vector (descending) for
%           residue pairs
% redundant vector that shows whether a residue pair is (largely) redundant
%           with a previous one, if so the number of the previous (better)
%           pair is given, if not a zero is given
% pairs     indices into model.coarse(model.current_structure).Ca_coor for
%           the radical pairs
%
% the total number of returned pairs nn can be larger than the number maxn
% of requested non-redundant pairs
%
% for algorithm, see:
% W. Zheng, B. R. Brooks, Biophys. J. 2006, 90: 4327-4336
%
% G. Jeschke, 2010-2012


global model
global hMain

maxlow=10; % maximum number of low-frequency modes
corr_cutoff=0.3; % correlation cutoff for considering residue pair as redundant 
max_cutoff=0.75; % maximum correlation cutoff in adaptive cutoff selection
min_r_CW=1.0;
max_r_CW=1.8;
min_r_DEER=1.7;
max_r_DEER=8.0;

if maxlow<maxn,
    maxlow=maxn;
end;

if nargin>5,
    maxlow=basis_size;
end;

if nargin<3,
    CW_flag=1;
end;

if nargin>3 && ~isempty(corr_cut),
    corr_cutoff=corr_cut;
end;

if nargin>4,
    max_r_DEER=max_dist;
end;

% determine network indices of sites

[ms,n]=size(site_list);
correspondence=zeros(1,ms); % correspondence table between sites and network nodes

nindices=model.coarse(model.current_structure).indices;
[mn,n]=size(nindices);

add_msg_board(sprintf('Selecting %i pairs for a basis of %i modes.\n',maxn,maxlow));
add_msg_board('Indexing network');
for k=1:ms,
    cindices=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).indices;
    for kk=1:mn,
        diff=abs(cindices-nindices(kk,:));
        if diff==0,
            correspondence(k)=kk;
            break
        end;
    end;
end;

[m,n]=size(model.ANM(model.current_structure).u);
if maxlow<1,
    maxlow=m-6;
end;
if maxlow==m-6,
    corr_cutoff=0.1;
end;
network0=model.coarse(model.current_structure).Ca_coor;
network=network0;
scal=sqrt(model.ANM(model.current_structure).lambda(7));
for k=7:6+maxlow,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network+scal*mode'/sqrt(model.ANM(model.current_structure).lambda(k));
end;

% determine distances and their standard deviations for all pairs
add_msg_board('Analyzing pair distance distributions');
dmatrix=zeros(ms,ms);

diff=zeros((ms)*(ms-1)/2,1);
pairs=zeros((ms)*(ms-1)/2,2);
ind=0;

for k=1:ms-1,
    NOpos1=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).NOpos;
    for kk=k+1:ms,
        NOpos2=model.sites{site_list(kk,1)}(site_list(kk,2)).residue(site_list(kk,3)).NOpos;
        [rm,sr]=analyze_distribution(NOpos1,NOpos2);
        dmatrix(k,kk)=rm;
        dmatrix(kk,k)=sr;
        if correspondence(k)~=0 && correspondence(kk)~=0, % corresponding network nodes exist
            r0=norm(network0(correspondence(k),:)-network0(correspondence(kk),:));
            r=norm(network(correspondence(k),:)-network(correspondence(kk),:));
            allowed=false;
            if CW_flag,
                if rm>=min_r_CW && rm<=max_r_CW, 
                    allowed=true; 
                end;
            end;
            if rm>=min_r_DEER && rm<= max_r_DEER,
                allowed=true;
            end;
            if allowed,
                ind=ind+1;
                pairs(ind,1)=k;
                pairs(ind,2)=kk;
                % the following code tried to include predictions of
                % measurement precision, this turned out to be detrimental
    %             if CW_flag,
    %                 prec=min([prec_DEER(rm,sr),prec_CW(rm,sr)]);
    %             else
    %                 prec=prec_DEER(rm,sr);
    %             end;
                diff(ind)=(r-r0); %/prec;
            end;
        end;
    end;
end;

add_msg_board('Analyzing network motion');
diff=diff(1:ind);
pairs=pairs(1:ind,:);
[diff1,poi]=sort(diff.^2,'descend');
if maxn>ind,
    maxn=ind;
end;
n=0;

while n<maxn && corr_cutoff<max_cutoff,
    n=0;
    nn=0;
    score=zeros(1,ind);
    redundant=zeros(1,ind);
    shortlist=zeros(maxn,5);
    pointers=zeros(1,maxn);
    rijm=zeros(maxn,maxlow);
    rijm0=zeros(1,maxlow);
    slowmodes=zeros(maxlow,m/3,3);
    for k=7:6+maxlow,
        evec=model.ANM(model.current_structure).u(:,k);
        mode=reshape(evec,3,m/3);
        slowmodes(k-6,:,:)=mode';
    end;

    scal=model.ANM(model.current_structure).lambda(7);
    while n<maxn && nn <ind,
        nn=nn+1;
        score(nn)=diff1(nn);
        for k=7:6+maxlow,
            mode=reshape(slowmodes(k-6,:,:),m/3,3);
            i=correspondence(pairs(poi(nn),1));
            j=correspondence(pairs(poi(nn),2));
            r0=norm(network0(i,:)-network0(j,:));
            r=norm(network(i,:)+mode(i,:)-network(j,:)-mode(j,:));
            rijm0(k-6)=(r-r0)/sqrt(model.ANM(model.current_structure).lambda(k));
        end;
        rijm0=rijm0/norm(rijm0);
        if nn>1 && n>=1,
            maxcorr=0;
            for j=1:n,
                num=sum(rijm0.*rijm(j,:));
                denom=sum(rijm0.^2);
                corr=abs(num/denom);
                if corr>corr_cutoff && corr>maxcorr,
                    redundant(nn)=j;
                    maxcorr=corr;
                end;
            end;
        end;
        if redundant(nn)==0,
            n=n+1;
            shortlist(n,1)=score(nn);
            shortlist(n,2:3)=pairs(poi(nn),:);
            shortlist(n,4)=dmatrix(pairs(poi(nn),1),pairs(poi(nn),2));
            shortlist(n,5)=dmatrix(pairs(poi(nn),2),pairs(poi(nn),1));
            pointers(n)=poi(nn);
            rijm(n,:)=rijm0;
        end;
    end;
    shortlist=shortlist(1:n,:);
    score=score(1:nn);
    redundant=redundant(1:nn);
    pairs=pairs(1:nn,:);
    if n<maxn,
        add_msg_board(sprintf('Only %i pairs found so far, while %i pairs were requested.',n,maxn));
        corr_cutoff=corr_cutoff+0.05;
        if corr_cutoff<=max_cutoff,
            add_msg_board(sprintf('Increasing correlation cutoff to %5.2f.',corr_cutoff));
        else
            add_msg_board('Warning: Maximum correlation cutoff reached. Less pairs than requested.');
        end;
    end;
end;

cutoff=corr_cutoff;
add_msg_board(sprintf('Pair evaluation finished at correlation cutoff %5.2f\n.',corr_cutoff));

fname=pair_report(shortlist,site_list,corr_cutoff);
hMain.report_file=fname;
report_editor;

function shortlist=label_pair_score_thermal(maxn,site_list,CW_flag,max_dist,constraint_sites,tolerance)
% function [shortlist,score,redundant]=label_pair_score_thermal(maxn,site_list,CW_flag,max_dist,constraint_sites,tolerance)
%
% determines up to maxn label pairs that are predicted to be most useful
% for driving state transitions with energy equipartitioning to normal
% modes of an elastic network model
%
% maxn      maximum number of residue pairs requested
% site_list list of label sites from which pairs can be created
% CW_flag   1: pair distances in CW EPR range are allowed, 0: they are not
%           allowed
% max_dist  maximum pair distance allowed [nm]
% constraint_sites  list of existing constraint site pairs (indices to
%                   network points)
% tolerance         tolerance level for using existing constraints as
%                   compared to optimum constraints
%
% shortlist a short list of non-redundant residue pairs sorted by
%           descending score, 1st column: score, 2nd/3rd column residue
%           numbers, 4th column mean label-label distance, 5th column
%           standard deviation of label-label distance
%
% the total number of returned pairs nn can be larger than the number maxn
% of requested non-redundant pairs
%
% G. Jeschke, 2010-2012


global model
global hMain

min_r_CW=1.0;
max_r_CW=1.8;
min_r_DEER=1.7;
max_r_DEER=8.0;

if nargin<3,
    CW_flag=1;
end;

if nargin>3,
    max_r_DEER=max_dist;
end;

% determine network indices of sites

[ms,n]=size(site_list);
correspondence=zeros(1,ms); % correspondence table between sites and network nodes

nindices=model.coarse(model.current_structure).indices;
[mn,n]=size(nindices);

[m,n]=size(model.ANM(model.current_structure).u);
network0=model.coarse(model.current_structure).Ca_coor;
labels0=zeros(size(network0));

mask=zeros(1,m/3);
back_corr=zeros(1,m/3);
add_msg_board(sprintf('Selecting %i pairs with energy equipartitioning.\n',maxn));
add_msg_board('Indexing network');
for k=1:ms,
    cindices=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).indices;
    for kk=1:mn,
        diff=abs(cindices-nindices(kk,:));
        if diff==0,
            correspondence(k)=kk;
            back_corr(kk)=k;
            mask(kk)=1;
            break
        end;
    end;
end;

mask(1:2)=0;
mask(end-1:end)=0;
mmask=kron(mask,mask');



% determine distances and their standard deviations for all pairs
add_msg_board('Analyzing pair distance distributions');
dmatrix=zeros(ms,ms);
large_dmatrix=zeros(m/3,m/3);
large_smatrix=zeros(m/3,m/3);

for k=1:ms-1,
    NOpos1=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).NOpos;
    pop=NOpos1(:,4);
    NOmean=pop'*NOpos1(:,1:3);
    i1=correspondence(k);
    labels0(i1,:)=NOmean;
    for kk=k+1:ms,
        NOpos2=model.sites{site_list(kk,1)}(site_list(kk,2)).residue(site_list(kk,3)).NOpos;
        [rm,sr]=analyze_distribution(NOpos1,NOpos2);
        dmatrix(k,kk)=rm;
        dmatrix(kk,k)=sr;
        i2=correspondence(kk);
        large_dmatrix(i1,i2)=rm;
        large_dmatrix(i2,i1)=rm;
        large_smatrix(i1,i2)=sr;
        large_smatrix(i2,i1)=sr;
    end;
end;

if CW_flag,
    mmask=mmask.*(large_dmatrix>=min_r_CW);
else
    mmask=mmask.*(large_dmatrix>=min_r_DEER);
end;
mmask=mmask.*(large_dmatrix<= max_r_DEER);

add_msg_board('Determining best site pairs');
shortlist=zeros(maxn,5);
poi=0;
% dmat0=coor2dmat(network0); % distance matrix for initial structure
dmat0=coor2dmat(labels0); % label-to-label distance matrix for initial structure
[mc,~]=size(constraint_sites);
for k=7:6+maxn,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network0+mode'; % network change by move along this mode
    dmat=label_dmat(network0,network,labels0,mask); % corresponding label-to-label-distance matrix
    % dmat=coor2dmat(network); % corresponding distance matrix 
    dmat1=dmat-dmat0; % distance change matrix
    dmat=abs(dmat1).*mmask; % mask all pairs that should not be selected at all
    for k0=1:poi, % mask all pairs that were selected before and very close pairs
        k1=shortlist(k0,2);
        k2=shortlist(k0,3);
        dmat(k1-2:k1+2,k2-2:k2+2)=0;
        dmat(k2-2:k2+2,k1-2:k1+2)=0;
    end;
    [ma,k0]=max(dmat); % determine site pair that corresponds to 
    [m2,k2]=max(ma);   % maximum distance change for this mode
    k1=k0(k2);
    if mc>0, % check whether an existing constraint pair is of similar quality
        dc=zeros(1,mc);
        for kc=1:mc,
            if constraint_sites(kc,1)~=0 && constraint_sites(kc,2)~=0,
                dc(kc)=dmat(constraint_sites(kc,1),constraint_sites(kc,2));  
            end;
        end;
        [md,kc]=max(dc);
        if md>0 && md>(1-tolerance)*m2,
            k1=constraint_sites(kc,1);
            k2=constraint_sites(kc,2);
            add_msg_board(sprintf('Accepted pair from constraint list for mode %i.',k));
        end;
    end;
    poi=poi+1;
    shortlist(poi,1)=dmat1(k1,k2);
    shortlist(poi,2)=k1;
    shortlist(poi,3)=k2;
    shortlist(poi,4)=large_dmatrix(k1,k2);
    shortlist(poi,5)=large_smatrix(k1,k2);
end;

[msl,nsl]=size(shortlist);
for k=1:msl,
    shortlist(k,2)=back_corr(shortlist(k,2));
    shortlist(k,3)=back_corr(shortlist(k,3));
end;

fname=pair_report(shortlist,site_list);
hMain.report_file=fname;
report_editor;

function shortlist=label_pair_score_basis(maxn,site_list,CW_flag,max_dist,threshold,constraint_sites,tolerance)
% function [shortlist,score,redundant]=label_pair_score_thermal(maxn,site_list,CW_flag,max_dist,threshold,constraint_sites,tolerance)
%
% determines up to maxn label pairs that are predicted to be most useful
% for driving state transitions with energy equipartitioning to normal
% modes of an elastic network model
%
% maxn      maximum number of residue pairs requested
% site_list list of label sites from which pairs can be created
% CW_flag   1: pair distances in CW EPR range are allowed, 0: they are not
%           allowed
% max_dist  maximum pair distance allowed [nm]
% threshold only if the norm of the distance change vector for a site pair
%           is at least threshold times the maximum norm of a distance
%           change vector, the site pair is considered, defaults to 0.4
% constraint_sites  list of existing constraints (indices into network)
% tolerance         tolerance level for accepting existing constraints
%                   instead of an optimum constraint
%
% shortlist a short list of non-redundant residue pairs sorted by
%           descending score, 1st column: score, 2nd/3rd column residue
%           numbers, 4th column mean label-label distance, 5th column
%           standard deviation of label-label distance
%
% the total number of returned pairs nn can be larger than the number maxn
% of requested non-redundant pairs
%
% G. Jeschke, 2010-2012


global model
global hMain

min_r_CW=1.0;
max_r_CW=1.8;
min_r_DEER=1.7;
max_r_DEER=8.0;

if nargin<3,
    CW_flag=1;
end;

if nargin>3,
    max_r_DEER=max_dist;
end;

if nargin<5,
    threshold=0.4;
end;

% determine network indices of sites

[ms,n]=size(site_list);
correspondence=zeros(1,ms); % correspondence table between sites and network nodes

nindices=model.coarse(model.current_structure).indices;
[mn,n]=size(nindices);

[m,n]=size(model.ANM(model.current_structure).u);
network0=model.coarse(model.current_structure).Ca_coor;
labels0=zeros(size(network0));

mask=zeros(1,m/3);
back_corr=zeros(1,m/3);
add_msg_board(sprintf('Selecting %i pairs with energy equipartitioning.\n',maxn));
add_msg_board('Indexing network');
for k=1:ms,
    cindices=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).indices;
    for kk=1:mn,
        diff=abs(cindices-nindices(kk,:));
        if diff==0,
            correspondence(k)=kk;
            back_corr(kk)=k;
            mask(kk)=1;
            break
        end;
    end;
end;

mask(1:2)=0;
mask(end-1:end)=0;
mmask=kron(mask,mask');

% determine distances and their standard deviations for all pairs
add_msg_board('Analyzing pair distance distributions');
dmatrix=zeros(ms,ms);
large_dmatrix=zeros(m/3,m/3);
large_smatrix=zeros(m/3,m/3);

for k=1:ms-1,
    NOpos1=model.sites{site_list(k,1)}(site_list(k,2)).residue(site_list(k,3)).NOpos;
    pop=NOpos1(:,4);
    NOmean=pop'*NOpos1(:,1:3);
    i1=correspondence(k);
    labels0(i1,:)=NOmean;
    for kk=k+1:ms,
        NOpos2=model.sites{site_list(kk,1)}(site_list(kk,2)).residue(site_list(kk,3)).NOpos;
        [rm,sr]=analyze_distribution(NOpos1,NOpos2);
        dmatrix(k,kk)=rm;
        dmatrix(kk,k)=sr;
        i2=correspondence(kk);
        large_dmatrix(i1,i2)=rm;
        large_dmatrix(i2,i1)=rm;
        large_smatrix(i1,i2)=sr;
        large_smatrix(i2,i1)=sr;
    end;
end;

if CW_flag,
    mmask=mmask.*(large_dmatrix>=min_r_CW);
else
    mmask=mmask.*(large_dmatrix>=min_r_DEER);
end;
mmask=mmask.*(large_dmatrix<= max_r_DEER);

add_msg_board(sprintf('%i site pairs are considered.',(m/3)*(m/3-1)/2));
% make pair list
ind=0; % pair pointer
pairs=zeros((m/3)*(m/3-1)/2,2); % residue index pairs
rpairs=zeros((m/3)*(m/3-1)/2,1); % pair distances
pindices=zeros(size(mmask));
for i=3:m/3-3,
    for j=i+1:m/3-2,
        if mmask(i,j),
            ind=ind+1;
            pairs(ind,1)=i;
            pairs(ind,2)=j;
            rpairs(ind)=large_dmatrix(i,j);
            pindices(i,j)=ind;
        end;
    end;
end;
pairs=pairs(1:ind,:);
rpairs=rpairs(1:ind,:);
% Determine existing constraint pairs
[mc,~]=size(constraint_sites);
existing=zeros(1,mc);
epoi=0;
for k=1:mc,
    c1=constraint_sites(k,1);
    c2=constraint_sites(k,2);
    comp=repmat([c1,c2],ind,1);
    diff=sum(abs(comp-pairs),2);
    i1=find(diff==0);
    if isempty(i1),
    comp=repmat([c2,c1],ind,1);
        diff=sum(abs(comp-pairs),2);
        i1=find(diff==0);
    end;
    if ~isempty(i1),
        epoi=epoi+1;
        existing(epoi)=i1;
    end;
end;
existing=existing(1:epoi);
add_msg_board(sprintf('%i site pairs are within distance range.',ind));
if epoi>1,
    add_msg_board(sprintf('Of those, %i site pairs are existing constraints.',epoi));
end;

drvecs=zeros(ind,maxn); % matrix of distance changes

add_msg_board('Determining best site pairs');
% set up the matrix of distance changes for all site pairs with respect to
% all modes in the basis
dmat0=coor2dmat(labels0); % distance matrix for initial structure
%dmat0=coor2dmat(network0); % distance matrix for initial structure
for k=7:6+maxn,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network0+mode'; % network change by move along this mode
    dmat=label_dmat(network0,network,labels0,mask); % corresponding label-to-label-distance matrix
    % dmat=coor2dmat(network); % corresponding distance matrix 
    dmat1=dmat-dmat0; % distance change matrix
    for i=3:m/3-3,
        for j=i+1:m/3-2,
            if pindices(i,j),
                drvecs(pindices(i,j),k-6)=dmat1(i,j);
            end;
        end;
    end;
end;

% Determine site pairs that are above the threshold for norm of distance
% change vector
above_thresh=zeros(1,ind);
for k=1:ind,
    above_thresh(k)=sqrt(sum(drvecs(k,:).^2)); % norm of distance change vector
end;
above_thresh=above_thresh/max(above_thresh); % normalization
selected=find(above_thresh>=threshold);
thresh_ext=threshold;
while length(selected)>36000,
    thresh_ext=thresh_ext+0.05;
    selected=find(above_thresh>=thresh_ext);
end;
while length(selected)<maxn && thresh_ext>0.05,
    thresh_ext=thresh_ext-0.05;
    selected=find(above_thresh>=thresh_ext);
end;

add_msg_board(sprintf('%i site pairs are above distance change threshold of %5.2f.',length(selected),thresh_ext));

drvecsn=zeros(length(selected),maxn);
% Normalize distance change vectors and find existing constraint pairs
sel_existing=zeros(1,length(existing));
sepoi=0;
for ks=1:length(selected),
    drvecsn(ks,:)=drvecs(selected(ks),:)/norm(drvecs(selected(ks),:));
    for kc=1:length(existing),
        if existing(kc)==selected(ks),
            sepoi=sepoi+1;
            sel_existing(sepoi)=ks;
        end;
    end;
end;
sel_existing=sel_existing(1:sepoi);
add_msg_board(sprintf('Of these, %i site pairs are existing distance constraints',sepoi));

% Determine linear dependence of site pairs
mindep=1e6;
sellist=zeros(1,maxn);
selfom=sellist;
lindep=eye(length(selected));
for k1=1:length(selected)-1,
    for k2=k1+1:length(selected),
        lindep(k1,k2)=abs(sum(drvecsn(k1,:).*drvecsn(k2,:)));
        lindep(k2,k1)=lindep(k1,k2);
        % The first two site pairs are the one that are least dependent with
        % respect to the basis of modes
        if lindep(k1,k2)<mindep,
            mindep=lindep(k1,k2);
            sellist(1)=k1;
            sellist(2)=k2;
        end;
    end;
end;

fprintf(1,'Minimum linear dependence of first two site pairs: %5.3f\n',mindep);
selfom(2)=mindep;
poi=2;
% add additional site pairs that are least linearly dependent with respect
% to existing pairs
while poi<maxn,
    selmat=lindep(:,sellist(1:poi));
    selvec=max(selmat,[],2);
    for k=1:poi,
        selvec(sellist(k))=1;
    end;
    [mindep,sel]=min(selvec);
    mindep_ex=1e12;
    best=0;
    for kc=1:sepoi,
        if selvec(sel_existing(kc))<mindep_ex,
            mindep_ex=selvec(sel_existing(kc));
            best=sel_existing(kc);
        end;
    end;
    if best && mindep_ex*(1-tolerance)<=mindep,
        add_msg_board(sprintf('Site pair %i replaced by existing constraint with linear dependence ratio %5.3f.',poi+1,mindep_ex/mindep));
        mindep=mindep_ex;
        sel=best;
    end;
    poi=poi+1;
    sellist(poi)=sel;
    selfom(poi)=mindep;
    fprintf(1,'Maximum linear dependence of new site pair %i: %5.3f\n',poi,mindep);
end;

% Create short list from list of selections
shortlist=zeros(maxn,5);
for k=1:maxn,
    pair=selected(sellist(k));
    shortlist(k,1)=selfom(k);
    shortlist(k,2)=pairs(pair,1);
    shortlist(k,3)=pairs(pair,2);
    shortlist(k,4)=large_dmatrix(pairs(pair,1),pairs(pair,2));
    shortlist(k,5)=large_smatrix(pairs(pair,1),pairs(pair,2));
end;

[msl,nsl]=size(shortlist);
for k=1:msl,
    shortlist(k,2)=back_corr(shortlist(k,2));
    shortlist(k,3)=back_corr(shortlist(k,3));
end;

fname=pair_report(shortlist,site_list,threshold,1);
hMain.report_file=fname;
report_editor;


function prec=prec_DEER(r,sigr)
% predicts a precision measure for DEER measurements as a function of
% expected mean distance r and standard deviation (width of the distance
% distribution) sigr

cl_DEER=0.2;
rl_DEER=1.8;
dl_DEER=0.2;
cu_DEER=0.2;
ru_DEER=5;
du_DEER=1;
sc_sigr=0.2;
c0_DEER=0.05;

argu=(r-ru_DEER)/du_DEER;
argl=(rl_DEER-r)/dl_DEER;
prec=c0_DEER+cl_DEER*exp(argl)+cu_DEER*exp(argu)+sc_sigr*sigr;

function prec=prec_CW(r,sigr)
% predicts a precision measure for CW EPR distance measurements as a function of
% expected mean distance r and standard deviation (width of the distance
% distribution) sigr

cl_CW=0.4;
rl_CW=1;
dl_CW=0.15;
cu_CW=0.4;
ru_CW=2;
du_CW=0.1;
sc_sigr=0.2;
c0_CW=0.2;

argu=(r-ru_CW)/du_CW;
argl=(rl_CW-r)/dl_CW;
prec=c0_CW+cl_CW*exp(argl)+cu_CW*exp(argu)+sc_sigr*sigr;

% --- Executes on button press in pushbutton_pairs.
function pushbutton_pairs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

[m,n]=size(handles.selection_list);
max_pairs=m*(m-1)/2;
if max_pairs<handles.num_pairs,
    add_msg_board('ERROR. Not enough sites selected for pair selection.');
    add_msg_board('Use "Add sites" first.');
    return
end;
handles=make_ANM(hObject,handles);
if ~handles.ANM,
    add_msg_board('Pair selection failed since ANM could not be computed.');
    guidata(hObject,handles);
    return;
end;

constraints=handles.constraints;
if ~isfield(constraints,'DEER') || isempty(constraints.DEER),
    add_msg_board('Warning: No DEER constraints present. Constraint list ignored.');
    handles.constraint_sites=[];
    constraints.DEER=[];
end;

ms=length(constraints.DEER);
constraint_sites=zeros(ms,2); % correspondence table between constraint sites and network nodes

nindices=model.coarse(model.current_structure).indices;
[mn,~]=size(nindices);

add_msg_board('Indexing existing constraint sites');
for k=1:ms,
    cindices=resolve_address(constraints.DEER(k).adr1);
    for kk=1:mn,
        diff=abs(cindices-nindices(kk,:));
        if diff==0,
            constraint_sites(k,1)=kk;
            break
        end;
    end;
    cindices=resolve_address(constraints.DEER(k).adr2);
    for kk=1:mn,
        diff=abs(cindices-nindices(kk,:));
        if diff==0,
            constraint_sites(k,2)=kk;
            break
        end;
    end;
end;

CW_flag=get(handles.checkbox_CW,'Value');
algorithm=get(handles.popupmenu_algorithm,'Value');
if isfield(handles,'constraint_sites') && ~isempty(handles.constraint_sites) && algorithm==1,
    add_msg_board('Warning: Existing constraint list is ignored for Zheng-Brooks selection.');
end;

hf=gcf;
set(hf,'Pointer','watch');

switch algorithm,
    case 1
        [shortlist,score,redundant,pairs,cutoff]=label_pair_score(handles.num_pairs,handles.selection_list,CW_flag,[],handles.max_dist,handles.basis_size);
    case 2
        shortlist=label_pair_score_basis(handles.num_pairs,handles.selection_list,CW_flag,handles.max_dist,handles.threshold,constraint_sites,handles.tolerance);
        cutoff=handles.threshold;
    case 3
        shortlist=label_pair_score_thermal(handles.num_pairs,handles.selection_list,CW_flag,handles.max_dist,constraint_sites,handles.tolerance);
        cutoff=[];
    otherwise
        shortlist=label_pair_score_basis(handles.num_pairs,handles.selection_list,CW_flag,handles.max_dist,handles.threshold,constraint_sites,handles.tolerance);
        cutoff=handles.threshold;
end;
set(hf,'Pointer','arrow');
handles.shortlist=shortlist;
if algorithm==1,
    handles.score=score;
    handles.redundant=redundant;
    handles.pairs=pairs;
end;

button = questdlg('Do you want to generate constraints for a target structure?','Option for generating target constraint list','No');
if strcmpi(button,'Yes')
    answer = inputdlg('Please provide structure tag','Select target structure');
    stag = strtrim(char(answer));
    if ~strcmp(stag(1),'['), stag=['[' stag]; end;
    if ~strcmp(stag(end),']'), stag=[stag ']']; end;
    [indices,msg]=resolve_address(stag);
    if ~msg.error && length(indices)==1,
        set(hf,'Pointer','watch');
        fname=target_constraints(indices,handles.shortlist,handles.selection_list,handles.basis_size,handles.num_pairs,cutoff,algorithm);
        set(hf,'Pointer','arrow');
        if ~isempty(fname),
            button = questdlg('Do you want to open target constraint file?','Overwrite of template list by target constraint list','No');
            if strcmpi(button,'Yes'),
                hMain.report_file=fname;
                report_editor;
            end;
        end;
    else
        add_msg_board(sprintf('ERROR: Target constraint generation failed since target structure %s does not exist.',stag));
    end;
end;

guidata(hObject,handles);

function edit_pairs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pairs as text
%        str2double(get(hObject,'String')) returns contents of edit_pairs as a double

[v,handles]=edit_update_MMM(handles,hObject,4,100,10,'%i',1);
handles.num_pairs=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_pairs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=make_ANM(hObject,handles)

global model

handles.ANM=false;
min_lambda=1e3*eps; % minimum frequency of a non-degenerate mode

set(gcf,'Pointer','watch');

h = msgbox('Please be patient. This can take several minutes.','Site pairs are selected');

snum=model.current_structure;

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
    if isempty(Ca_coor) || length(masses)<2,
        set(gcf,'Pointer','arrow');
        add_msg_board('ERROR: Cannot create ANM, since current structure has less than two amino acid residues.');
        if ishandle(h),
            delete(h);
        end;        
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
end;


if ~isfield(model,'ANM') || length(model.ANM)<snum || isempty(model.ANM(snum).u),
    Hessian=setup_ANM_bonded(Ca_coor);
    contacts=[];
    [u,D]=eig(Hessian);
    [m,n]=size(Hessian);

    lambda=zeros(1,m);
    for k=1:m,
        lambda(k)=D(k,k);
    end;
    clear D
    model.ANM(snum).lambda=lambda;
    model.ANM(snum).u=u;
    model.ANM(snum).residues=m/3;
    model.ANM(snum).contacts=contacts;
    msf=zeros(1,m/3);
    suspicious=zeros(1,m);
    suspicious(1:6)=ones(1,6);
    fs=0;
    fs_flag=true;
    for k=7:m,
        if lambda(k)>min_lambda;
            mode=reshape(u(:,k),3,m/3);
            msf0=sum(mode.^2,1);
            [ma,poi]=max(abs(msf0));
            msf1=msf0;
            msf1(poi)=0;
            ma2=max(abs(msf1));
            if ma2>ma/3,
                msf=msf+msf0/lambda(k);
            else
                suspicious(k)=1;
                if fs_flag,
                    fs=k;
                    fs_flag=false;
                end;
            end;
        end;
    end;
    model.ANM(snum).msf=msf;
    model.ANM(snum).suspicious=suspicious;
end;


handles.ANM=true;

if ishandle(h),
    delete(h);
end;        

set(gcf,'Pointer','arrow');

function [rm,sr] = analyze_distribution(NOpos1,NOpos2)

pop1=NOpos1(:,4);
n1=length(pop1);
xyz1=NOpos1(:,1:3)/10; % divided by 10 for ? -> nm

pop2=NOpos2(:,4);
n2=length(pop2);
xyz2=NOpos2(:,1:3)/10; % divided by 10 for ? -> nm

xyz1sq = repmat(sum(xyz1.^2,2),1,n2);
xyz2sq = repmat(sum(xyz2.^2,2),1,n1).';
r = sqrt(abs(xyz1sq + xyz2sq - 2*xyz1*xyz2.'));
r = r(:);
pop = pop1*pop2.';
pop = pop(:);

rm = sum(r.*pop)/sum(pop);
dr = r-rm;
sr = sqrt(0.01+n1*n2*sum(pop.*dr.^2)/(sum(pop)*(n1*n2-1)));

% function [rm,sr]=analyze_distribution(NOpos1,NOpos2)
% 
% pop1=NOpos1(:,4);
% n1=length(pop1);
% xyz1=NOpos1(:,1:3)/10; % divided by 10 for ? -> nm
% pop2=NOpos2(:,4);
% n2=length(pop2);
% xyz2=NOpos2(:,1:3)/10; % divided by 10 for ? -> nm
% rij=zeros(1,n1*n2);
% popij=rij;
% poi=0;
% for k1=1:n1,
%     for k2=1:n2,
%       poi=poi+1;
%       popij(poi)=pop1(k1)*pop2(k2);
%       rij(poi)=norm(xyz1(k1,:)-xyz2(k2,:));
%     end;
% end;
% rm=sum(rij.*popij)/sum(popij);
% dr=rij-rm;
% sr=sqrt(0.01+n1*n2*sum(popij.*dr.^2)/(sum(popij)*(n1*n2-1)));


function fname=pair_report(shortlist,site_list,cutoff,thermal)

global general
global model

if nargin<4,
    thermal=false;
end;

snum=model.current_structure;
adr=mk_address(snum);

[m,n]=size(shortlist);
fname=[general.tmp_files 'suggested_pairs.dat'];
fid=fopen(fname,'w');
fprintf(fid,'%% MMM site pair suggestion for structure %s\n',adr);
if nargin>2 && ~thermal,
    fprintf(fid,'%% %i pairs were found at correlation cutoff %5.2f\n',m,cutoff);
elseif thermal
    fprintf(fid,'%% %i pairs were found for energy equipartitioning with threshold %5.2f.\n',m,cutoff);
else
    fprintf(fid,'%% %i pairs were found for energy equipartitioning.\n',m);
end;
fprintf(fid,'%% Site 1         Site 2            <r> (nm)  sigr (nm)  Score\n');
for k=1:m,
    ind1=site_list(shortlist(k,2),:);
    cind1=model.sites{ind1(1)}(ind1(2)).residue(ind1(3)).indices;
    [stag,ctag,modelnum,resnum]=mk_address_parts(cind1);
    adr1=sprintf('(%s)%i',ctag,resnum);
    % adr1=mk_address(cind1,1);
    while length(adr1)<15,
        adr1=[adr1 ' '];
    end;
    ind2=site_list(shortlist(k,3),:);
    cind2=model.sites{ind2(1)}(ind2(2)).residue(ind2(3)).indices;
    [stag,ctag,modelnum,resnum]=mk_address_parts(cind2);
    adr2=sprintf('(%s)%i',ctag,resnum);
    % adr2=mk_address(cind2,1);
    while length(adr2)<15,
        adr2=[adr2 ' '];
    end;
    fprintf(fid,'%15s  %15s %6.2f  %8.2f %%%12.4f\n',adr1,adr2,shortlist(k,4),shortlist(k,5),shortlist(k,1));
end;
fclose(fid);


% --- Executes on button press in checkbox_CW.
function checkbox_CW_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CW



function edit_max_dist_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_dist as text
%        str2double(get(hObject,'String')) returns contents of edit_max_dist as a double

[v,handles]=edit_update_MMM(handles,hObject,1.8,20,6,'%4.1f',0);
handles.max_dist=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_max_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_basis_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_basis_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_basis_size as text
%        str2double(get(hObject,'String')) returns contents of edit_basis_size as a double

[v,handles]=edit_update_MMM(handles,hObject,4,100,10,'%i',1);
handles.basis_size=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_basis_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_basis_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_matched_basis.
function checkbox_matched_basis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_matched_basis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_matched_basis

checked=get(hObject,'Value');
if checked,
    handles.basis_size=handles.num_pairs;
    set(handles.edit_basis_size,'String',sprintf('%i',handles.basis_size));
    set(handles.edit_basis_size,'Enable','off');
else
    set(handles.edit_basis_size,'Enable','on');
end;
guidata(hObject,handles);

function fname=target_constraints(snum,shortlist,site_list,basis_size,np,cutoff,algorithm)

global model

stag=mk_address(snum);
stag=stag(2:end-1);
add_msg_board(sprintf('Generating target constraints for structure [%s]',stag));

[m,n]=size(shortlist);
switch algorithm
    case 1
        fname=[stag sprintf('_%i_distances_ZB_%i.dat',m,round(100*cutoff))];
        alg_string='Zheng/Brooks';
    case 2
        fname=[stag sprintf('_%i_distances_TM_%i.dat',m,round(100*cutoff))];
        alg_string='thermal move';
    case 3
        fname=[stag sprintf('_%i_distances_MC.dat',m)];
        alg_string='maximum change';
    otherwise
        fname=[stag sprintf('_%i_distances_TM_%i.dat',m,round(100*cutoff))];
        alg_string='thremal move';
end;
fid=fopen(fname,'w');
fprintf(fid,'%% MMM simulated label-to-label distance constraints for target structure [%s]\n',stag);
fprintf(fid,'%% %s pair selection algorithm\n',alg_string);
switch algorithm
    case 1
        fprintf(fid,'%% Correlation cutoff: %5.2f\n',cutoff);
    case 2
        fprintf(fid,'%% Threshold: %5.2f\n',cutoff);
end;
fprintf(fid,'%% Requested constraints: %i; Generated constraints: %i\n',np,m);
fprintf(fid,'# BASIS %i\n',basis_size);
fprintf(fid,'# TARGET %s\n',stag);
fprintf(fid,'# DEER MTSL 175\n');
for k=1:m,
    ind1=site_list(shortlist(k,2),:);
    cind1=model.sites{ind1(1)}(ind1(2)).residue(ind1(3)).indices;
    [stag0,ctag1,modelnum,resnum1]=mk_address_parts(cind1);
    adr1=sprintf('[%s](%s)%i',stag,ctag1,resnum1);
    target_sites(k).adr1=adr1;
    % adr1=mk_address(cind1,1);
    while length(adr1)<15,
        adr1=[adr1 ' '];
    end;
    ind2=site_list(shortlist(k,3),:);
    cind2=model.sites{ind2(1)}(ind2(2)).residue(ind2(3)).indices;
    [stag0,ctag2,modelnum,resnum2]=mk_address_parts(cind2);
    adr2=sprintf('[%s](%s)%i',stag,ctag2,resnum2);
    target_sites(k).adr2=adr2;
    % adr2=mk_address(cind2,1);
end;
[DEER,cancelled]=label_target_structure(snum,target_sites);
if cancelled || isempty(DEER),
    add_msg_board('ERROR: Labeling of target structure failed. No constraint file was generated.');
    fname='';
    fclose(fid);
    return
end;
m=length(DEER);
for k=1:m,
    adr1=DEER(k).adr1;
    poi = strfind(adr1,']');
    poi=poi(end);
    adr1=adr1(poi+1:end);
%     poi = strfind(adr1,')');
%     resnum=str2double(adr1(poi+1:end));
%     adr1=sprintf('%s%i',adr1(1:poi),resnum+200);
    while length(adr1)<15,
        adr1=[adr1 ' '];
    end;
    adr2=DEER(k).adr2;
    poi = strfind(adr2,']');
    poi=poi(end);
    adr2=adr2(poi+1:end);
%     poi = strfind(adr2,')');
%     resnum=str2double(adr2(poi+1:end));
%     adr2=sprintf('%s%i',adr2(1:poi),resnum+200);
    while length(adr2)<15,
        adr2=[adr2 ' '];
    end;
    fprintf(fid,'%15s  %15s %6.2f  %8.2f\n',adr1,adr2,DEER(k).r,DEER(k).sigr);
end;
fprintf(fid,'# END\n');
fclose(fid);

function [DEER,cancelled]=label_target_structure(snum,target_sites)

global model
global hMain

cancelled=false;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

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
    for k=1:length(target_sites),
        adr1=target_sites(k).adr1;
        ind1=resolve_address(adr1);
        target_sites(k).indices1=ind1;
        if isempty(ind1),
            add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current target %s',mk_address(snum)));
            add_msg_board('Generation of distance constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,'MTSL'),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),'MTSL'),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}='MTSL';
                T_list(poi)=175;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.','MTSL',adr1));
            end;
        end;
        adr2=target_sites(k).adr2;
        ind2=resolve_address(adr2);
        target_sites(k).indices2=ind2;
        if isempty(ind2),
            add_msg_board(sprintf('ERROR: Constraint %i has second label at site %s',k,adr2));
            add_msg_board(sprintf('This site does not exist in current target %s',mk_address(snum)));
            add_msg_board('Generation of distance constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind2-lindices(l,:);
            if sum(abs(diff))==0  && strcmpi(labels(l).name,'MTSL'),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),'MTSL'),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                label_list{poi}='MTSL';
                T_list(poi)=175;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.','MTSL',adr2));
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(target_sites),
        adr1=target_sites(k).adr1;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),'MTSL'),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            label_list{poi}='MTSL';
            T_list(poi)=175;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
        adr2=target_sites.adr2;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),'MTSL'),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr2;
            label_list{poi}='MTSL';
            T_list(poi)=175;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr2));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

target_sites=add_label_information(model.sites,target_sites);

for k=1:length(target_sites),
    adr1=target_sites(k).adr1;
    ind1=resolve_address(adr1);
    adr2=target_sites(k).adr2;
    ind2=resolve_address(adr2);
    NOpos1=target_sites(k).NOpos1;
    NOpos2=target_sites(k).NOpos2;
    [rmean,sigr]=analyze_distribution(NOpos1,NOpos2);
    DEER(k).r=rmean;
    DEER(k).sigr=sigr;
    DEER(k).indices=[ind1;ind2];
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
end;


function target_sites=add_label_information(sites,target_sites)

global model

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            indices=sites{k0}(k1).residue(k).indices;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            for kk=1:length(target_sites),
                mismatch=sum(abs(indices-target_sites(kk).indices1));
                if ~mismatch,
                    target_sites(kk).NOpos1=NOpos;
                end;
                mismatch=sum(abs(indices-target_sites(kk).indices2));
                if ~mismatch,
                    target_sites(kk).NOpos2=NOpos;
                end;
            end;
        end;
    end;
end;

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

function dmat=label_dmat(network0,network,labels0,mask)

global model

[maxnum,m]=size(network);
labels1=zeros(size(labels0));
scarce=0;
cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network
% update of label coordinates
for k=1:maxnum,
    if mask(k), % this network node has a label attached
        xyz=labels0(k,:);
        local_template=zeros(5,3);
        local_template_0=zeros(5,3);
        % make a local template to fit rotation and translation
        poi=0;
        for kk=-2:2,
            if k+kk>0 && k+kk<=maxnum, % is addressed residue a network point?
                diff=cindices(k+kk,4)-cindices(k,4);
                if diff==kk, % is addressed residue part of a continuous segment?
                    poi=poi+1;
                    local_template(poi,:)=network(k+kk,:);
                    local_template_0(poi,:)=network0(k+kk,:);
                end;
            end;
        end;
        if poi>=3, % found sufficient number of points to determine local rotation and translation
            [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
            xyz=[xyz 1];
            xyz=transmat*xyz';
            xyz=xyz';
            labels1(k,:)=xyz(1:3);
        else
            labels1(k,:)=xyz+network(k,:)-network0(k,:);
            scarce=scarce+1;
        end;
    end;
end;

dmat=coor2dmat(labels1); % corresponding distance matrix 


% --- Executes on button press in checkbox_thermal.
function checkbox_thermal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_thermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_thermal



function edit_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_threshold as a double

[v,handles]=edit_update_MMM(handles,hObject,0.1,0.9,0.4,'%4.2f',0);
handles.threshold=v;
guidata(hObject,handles);

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


% --- Executes on selection change in popupmenu_algorithm.
function popupmenu_algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_algorithm


% --- Executes during object creation, after setting all properties.
function popupmenu_algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_constraints.
function pushbutton_constraints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global model

% Get constraint file name and read constraint file
my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load constraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Constraint loading cancelled by user');
    cd(my_path);
    guidata(hObject,handles);
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    constraints=rd_restraints(fullfile(pname,fname));
end;
cd(my_path);

if isfield(constraints,'DEER') && ~isempty(constraints.DEER),
    add_msg_board(sprintf('%i existing DEER constraints will be considered.',length(constraints.DEER)));
else
    add_msg_board('Warning: No DEER constraints found in file.');
end;
handles=rmfield(handles,'constraints');
handles.constraints=constraints;
% determine network indices of sites

guidata(hObject,handles);

function edit_constraint_tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_constraint_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_constraint_tolerance as text
%        str2double(get(hObject,'String')) returns contents of edit_constraint_tolerance as a double

[v,handles]=edit_update_MMM(handles,hObject,0,1,0.1,'%4.2f',0);
handles.tolerance=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_constraint_tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_constraint_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
