function varargout = monitor_distances(varargin)
% MONITOR_DISTANCES M-file for monitor_distances.fig
%      MONITOR_DISTANCES, by itself, creates a new MONITOR_DISTANCES or raises the existing
%      singleton*.
%
%      H = MONITOR_DISTANCES returns the handle to a new MONITOR_DISTANCES or the handle to
%      the existing singleton*.
%
%      MONITOR_DISTANCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONITOR_DISTANCES.M with the given input arguments.
%
%      MONITOR_DISTANCES('Property','Value',...) creates a new MONITOR_DISTANCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before monitor_distances_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to monitor_distances_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help monitor_distances

% Last Modified by GUIDE v2.5 13-Dec-2010 14:32:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @monitor_distances_OpeningFcn, ...
                   'gui_OutputFcn',  @monitor_distances_OutputFcn, ...
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


% --- Executes just before monitor_distances is made visible.
function monitor_distances_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to monitor_distances (see VARARGIN)

global graph_settings

% Choose default command line output for monitor_distances
handles.output = hObject;

handles.mode=1;
handles.sign=1;
handles.amplitude=1;
handles.pair_number=0;
handles.scaling=10; % 10 nm distance change at normal mode value 1

handles.computation.res=100;

handles.name_selection=10;
handles.named_color=graph_settings.colors(handles.name_selection,:);
set(handles.text_color,'ForegroundColor',handles.named_color);
set(handles.text_preselected,'ForegroundColor',handles.named_color);

set(handles.listbox_colors,'String',graph_settings.color_names);
set(handles.listbox_colors,'Value',handles.name_selection);

[handles,success]=mk_label_lists(handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

handles.experimental_distance=2.50;

handles.pairs(1).adr1='';
handles.pairs(1).indices1=[];
handles.pairs(1).adr2='';
handles.pairs(1).indices2=[];
handles.pairs(1).color='';
handles.pairs(1).rgb=[];
handles.pairs(1).r=[];

new_pairs{1}='<no pair selected>';
set(handles.listbox_pairs,'String',new_pairs);

if ~success,
   msgbox('Model must feature at least two spin labels or selected atoms','Distance monitoring impossible','error');
   figure1_CloseRequestFcn(hObject, eventdata, handles); 
   return
end;

handles=setup_preselection(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes monitor_distances wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = monitor_distances_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in listbox_site_1.
function listbox_site_1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_site_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_site_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_site_1

handles=setup_preselection(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_site_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_site_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_site_2.
function listbox_site_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_site_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_site_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_site_2

handles=setup_preselection(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_site_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_site_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_colors.
function listbox_colors_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_colors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_colors

handles=setup_preselection(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_colors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global graph_settings

selc=get(handles.listbox_colors,'Value');
colors=get(handles.listbox_colors,'String');
colorname=colors{selc};
rgb=graph_settings.colors(selc,:);

sel_1=get(handles.listbox_site_1,'Value');
addresses=get(handles.listbox_site_1,'String');
adr1=addresses{sel_1};

sel_2=get(handles.listbox_site_2,'Value');
addresses=get(handles.listbox_site_2,'String');
adr2=addresses{sel_2};

pair=sprintf('%s // %s ; %s',adr1,adr2,colorname);
pair_short_1=sprintf('%s // %s ; %s',adr1,adr2);
pair_short_2=sprintf('%s // %s ; %s',adr2,adr1);

all_pairs=get(handles.listbox_pairs,'String');

if handles.pair_number==0,
    handles.pairs(1).adr1=adr1;
    handles.pairs(1).indices1=resolve_address(adr1);
    handles.pairs(1).adr2=adr2;
    handles.pairs(1).indices2=resolve_address(adr2);
    handles.pairs(1).color=colorname;
    handles.pairs(1).rgb=rgb;
    if get(handles.checkbox_experimental,'Value'),
        handles.pairs(1).r=handles.experimental_distance;
    else
        handles.pair(1).r=[];
    end;
    handles.pair_number=1;
    all_pairs{1}=pair;
else
    found=false;
    for k=1:length(all_pairs),
        if sum(strfind(all_pairs{k},pair_short_1)) || sum(strfind(all_pairs{k},pair_short_2)),
            found=true;
            break;
        end;
    end;
    if found,
        add_msg_board('Warning: Same label pair cannot be selected twice!.');
    else
        poi=handles.pair_number+1;
        handles.pairs(poi).adr1=adr1;
        handles.pairs(poi).indices1=resolve_address(adr1);
        handles.pairs(poi).adr2=adr2;
        handles.pairs(poi).indices2=resolve_address(adr2);
        handles.pairs(poi).color=colorname;
        handles.pairs(poi).rgb=rgb;
        handles.pair_number=poi;
    if get(handles.checkbox_experimental,'Value'),
        handles.pairs(poi).r=handles.experimental_distance;
    else
        handles.pair(poi).r=[];
    end;
        all_pairs{poi}=pair;
    end;
end;
set(handles.listbox_pairs,'String',all_pairs);
handles=update(handles);

guidata(hObject,handles);



% --- Executes on button press in pushbutton_minus.
function pushbutton_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selc=get(handles.listbox_colors,'Value');
colors=get(handles.listbox_colors,'String');
colorname=colors{selc};

sel_1=get(handles.listbox_site_1,'Value');
addresses=get(handles.listbox_site_1,'String');
adr1=addresses{sel_1};

sel_2=get(handles.listbox_site_2,'Value');
addresses=get(handles.listbox_site_2,'String');
adr2=addresses{sel_2};

pair_short_1=sprintf('%s // %s ; %s',adr1,adr2);
pair_short_2=sprintf('%s // %s ; %s',adr2,adr1);

all_pairs=get(handles.listbox_pairs,'String');

if handles.pair_number==0,
    add_msg_board('Warning: Cannot unselect label pair from empty selection list.');
else
    found=false;
    for k=1:length(all_pairs),
        if sum(strfind(all_pairs{k},pair_short_1)) || sum(strfind(all_pairs{k},pair_short_2)),
            found=true;
            unsel=k;
            break;
        end;
    end;
    if ~found,
        add_msg_board('Warning: Label pair cannot be unselected, because it is not in selection list.');
        new_pairs=all_pairs;
    else
        poi=0;
        for k=1:length(all_pairs),
            if k~=unsel,
                poi=poi+1;
                handles.pairs(poi).adr1=handles.pairs(k).adr1;
                handles.pairs(poi).indices1=handles.pairs(k).indices1;
                handles.pairs(poi).adr2=handles.pairs(k).adr2;
                handles.pairs(poi).indices2=handles.pairs(k).indices2;
                handles.pairs(poi).color=handles.pairs(k).color;
                handles.pairs(poi).rgb=handles.pairs(k).rgb;
                handles.pairs(poi).r=handles.pairs(k).r;
                new_pairs{poi}=all_pairs{k};
            end;
        end;
        handles.pair_number=poi;
    end;
end;
if handles.pair_number==0,
    new_pairs{1}='<no pair selected>';
    handles.pairs(1).adr1='';
    handles.pairs(1).indices1=[];
    handles.pairs(1).adr2='';
    handles.pairs(1).indices2=[];
    handles.pairs(1).color='';
    handles.pairs(1).rgb=[];
    handles.pairs(1).r=[];
end;
set(handles.listbox_pairs,'Value',1);
set(handles.listbox_pairs,'String',new_pairs);
handles=update(handles);
guidata(hObject,handles);


% --- Executes on selection change in listbox_pairs.
function listbox_pairs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_pairs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_pairs


% --- Executes during object creation, after setting all properties.
function listbox_pairs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_compute.
function pushbutton_compute_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=compute(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

function [handles,success]=mk_label_lists(handles)

global model

success=0;

indices=resolve_address('*');

if ~isfield(model,'labels')
    labels=0;
    if isempty(indices),
        return
    end;
end;

sites=0;

if isfield(model,'labels') && ~isempty(model.labels),
    sites=length(model.labels);
    labels=sites;
    for k=1:length(model.labels),
        label_list{k}=model.labels(k).adr;
    end;
end;


% make list of selected atoms
if ~isempty(indices),
    [m,n]=size(indices);
    poi=0;
    handles.atoms=zeros(m,5);
    for ko=1:m, % loop over all objects
        idepth=length(find(indices(ko,:)>0)); % determine type of current object
        if idepth==5,
            poi=poi+1;
            cindices=indices(ko,1:idepth);
            address=mk_address(cindices,true);
            label_list{poi+sites}=address;
            handles.atoms(poi,:)=cindices;
            handles.atom_adr{poi}=address;
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
    set(handles.listbox_site_1,'String',' ');
    set(handles.listbox_site_2,'String',' ');
    set(handles.listbox_site_1,'Value',1);
    set(handles.listbox_site_2,'Value',1);
    set(handles.text_preselected,'String','<not enough labels and selected atoms>');
    return;
end;

set(handles.listbox_site_1,'String',label_list);
set(handles.listbox_site_1,'Value',1);
set(handles.listbox_site_2,'String',label_list);
set(handles.listbox_site_2,'Value',2);

function handles=setup_preselection(handles)

global graph_settings

selc=get(handles.listbox_colors,'Value');
colors=get(handles.listbox_colors,'String');
colorname=colors{selc};
rgb=graph_settings.colors(selc,:);
set(handles.text_preselected,'ForegroundColor',rgb);
set(handles.text_color,'ForegroundColor',rgb);

sel_1=get(handles.listbox_site_1,'Value');
addresses=get(handles.listbox_site_1,'String');
adr1=addresses{sel_1};

sel_2=get(handles.listbox_site_2,'Value');
addresses=get(handles.listbox_site_2,'String');
adr2=addresses{sel_2};

if sel_1==sel_2,
    set(handles.pushbutton_plus,'Enable','off');
    set(handles.pushbutton_minus,'Enable','off');
    set(handles.text_preselected,'String','<Warning: pair cannot consist of identical sites>');
else
    set(handles.pushbutton_plus,'Enable','on');
    set(handles.pushbutton_minus,'Enable','on');
    set(handles.text_preselected,'String',sprintf('%s // %s ; %s',adr1,adr2,colorname));
end;

function handles=update(handles)
% Display update

hf1=figure(1); clf;
set(gca,'FontSize',10);
hold on;
set(hf1,'Name','Distance distributions and mean distances');
if handles.pair_number<1,
    return
end;
hs1=subplot(2,1,1);
hold on;
title('Distance distributions');
xlabel('distance [nm]');
rmin=1;
rmax=6;
ymin=1e6;
ymax=0;
hdistr=zeros(1,handles.pair_number);
hmean=zeros(1,handles.pair_number);
rmean=zeros(handles.pair_number,21);
% Display distance distributions
for k=1:handles.pair_number,
    [rax,act_distr]=pair_distribution(handles,handles.pairs(k).adr1,handles.pairs(k).adr2);
    rmean(k,1)=sum(rax.*act_distr)/sum(act_distr);
    minpoi=find(act_distr>0.01*max(act_distr),1);
    maxpoi=find(act_distr>0.01*max(act_distr),1,'last');
    if rax(minpoi)<rmin, rmin=rax(minpoi); end;
    if rax(maxpoi)>rmax, rmax=rax(maxpoi); end;
    if min(act_distr)<ymin, ymin=min(act_distr); end;
    if max(act_distr)>ymax, ymax=max(act_distr); end;    
    hdistr(k)=plot(rax,act_distr,'Color',handles.pairs(k).rgb); % store handle
    if ~isempty(handles.pairs(k).r),
        plot([handles.pairs(k).r,handles.pairs(k).r],[0,1.1*max(act_distr)],':','Color',handles.pairs(k).rgb);
    end;
end;
ymin=0.9*ymin;
if ymin>0, ymin=0; end;
ymax=1.1*ymax;
axis([rmin,rmax,ymin,ymax]);

hs2=subplot(2,1,2);
hold on;
title('Mean distances');
ylabel('distance [nm]');
xlabel('phase angle [degree]');
x=0;
r_exp=zeros(1,handles.pair_number);
% Display mean distances
for k=1:handles.pair_number,
    y=rmean(k,1);
    % plot(x,y,'x','Color',handles.pairs(k).rgb);
    hmean(k)=plot(x,y,'Color',handles.pairs(k).rgb); % store handle
    if ~isempty(handles.pairs(k).r),
        plot([-5,185],[handles.pairs(k).r,handles.pairs(k).r],':','Color',handles.pairs(k).rgb);
        r_exp(k)=handles.pairs(k).r;
    end;
end;
axis([-5,185,0.9*min([min(rmean(:,1)),min(r_exp)]),1.1*max([max(rmean(:,1)),max(r_exp)])]);
% set(hmean(1),'XData',[0,18],'YData',[rmean(1,1),1.2*rmean(1,1)]);

handles.hdistr=hdistr;
handles.hmean=hmean;
handles.rmean=rmean;
handles.axes_mean=hs2;
handles.axes_distr=hs1;
handles.distr_max=ymax;

function [rax,distr]=pair_distribution(handles,adr1,adr2,network0,network1,residue1,residue2)
% Distance distribution function for a pair of spin labels, 
% if optional arguments network0 and network1 are given all NO positions
% are first transformed from the ones corresponding to coarse-grained 
% network0 to those corresponding to network1, 
% in this case residue1 and residue2 are indices of the labelled residues
% in the coarse-grained (C_alpha) networks, and update_labels is called

global model

sig=0.1; % Gaussian broadening of the distance distribution

if isfield(model,'labels'),
    for k=1:length(model.labels),
        if strcmp(model.labels(k).adr,adr1),
            NOpos1=model.labels(k).NOpos;
            if nargin>3,
                NOpos1=update_labels(NOpos1,network0,network1,residue1);
            end;
        end;
        if strcmp(model.labels(k).adr,adr2),
            NOpos2=model.labels(k).NOpos;
            if nargin>3,
                NOpos2=update_labels(NOpos2,network0,network1,residue2);
            end;
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
            if nargin>3,
                NOpos1=update_labels(NOpos1,network0,network1,residue1);
            end;
        end;
        if strcmp(handles.atom_adr{k},adr2),
            indices=handles.atoms(k,:);
            [msg0,xyz]=get_atom(indices,'xyz');
            [msg0,pop]=get_atom(indices,'populations');
            NOpos2=[xyz pop];
            if nargin>3,
                NOpos2=update_labels(NOpos2,network0,network1,residue2);
            end;
        end;
    end;
end;

[rax,distr]=get_distribution(NOpos1,NOpos2,sig);



function edit_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amplitude as text
%        str2double(get(hObject,'String')) returns contents of edit_amplitude as a double


[v,handles]=edit_update_MMM(handles,hObject,0.01,100,1.00,'%4.2f',0);
handles.amplitude=v;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amplitude (see GCBO)
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

entry=strcat(help_files,'monitor_distances.html');
webcall(entry,'-helpbrowser');

function edit_mode_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mode_number as text
%        str2double(get(hObject,'String')) returns contents of edit_mode_number as a double

global model

[v,handles]=edit_update_MMM(handles,hObject,1,3*model.ANM(model.current_structure).residues-6,handles.mode,'%4i',1);
handles.mode=v;
set(handles.figure1,'Name',sprintf('Monitoring distance changes in normal mode %i',v));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_mode_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_mode_minus.
function pushbutton_mode_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mode_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.mode>1,
    handles.mode=handles.mode-1;
end;
set(handles.edit_mode_number,'String',sprintf('%4i',handles.mode));
set(handles.figure1,'Name',sprintf('Monitoring distance changes in normal mode %i',handles.mode));

guidata(hObject,handles);

% --- Executes on button press in pushbutton_mode_plus.
function pushbutton_mode_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mode_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.mode<3*model.ANM(model.current_structure).residues-6,
    handles.mode=handles.mode+1;
end;
set(handles.edit_mode_number,'String',sprintf('%4i',handles.mode));
set(handles.figure1,'Name',sprintf('Monitoring distance changes in normal mode %i',handles.mode));

guidata(hObject,handles);

% --- Executes on button press in togglebutton_sign.
function togglebutton_sign_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_sign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_sign

if get(hObject,'Value'),
    set(hObject,'String','+');
    handles.sign=1;
else
    set(hObject,'String','-');
    handles.sign=-1;
end;
guidata(hObject,handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


function handles=compute(handles)

global model

Ca_coor=model.coarse(model.current_structure).Ca_coor;
rindices=model.coarse(model.current_structure).indices;
[m,n]=size(rindices);

% make index list of pairs into coarse-grained (C_alpha) model
pair_pointers=zeros(handles.pair_number,2);
for k=1:handles.pair_number,
    cindices1=handles.pairs(k).indices1;
    cindices2=handles.pairs(k).indices2;
    for kk=1:m,
        if ~sum(abs(cindices1(1:4)-rindices(kk,:))),
            pair_pointers(k,1)=kk;
            continue;
        end;
    end;
    for kk=1:m,
        if ~sum(abs(cindices2(1:4)-rindices(kk,:))),
            pair_pointers(k,2)=kk;
            continue;
        end;
    end;
end;

x0=Ca_coor(:,1);
y0=Ca_coor(:,2);
z0=Ca_coor(:,3);

m=model.ANM(model.current_structure).residues;
evec=model.ANM(model.current_structure).u(:,handles.mode+6);
mode=reshape(evec,3,m);
mode=mode';

msf=model.ANM(model.current_structure).msf;
Bfactors=model.coarse(model.current_structure).Bfactors;
if sum(Bfactors)==0,
    Bfactors=20*ones(size(BFactors));
end;
[msf,p]=lin_fit(msf,Bfactors);
Bfactor_scale=p(1);

freq_sq=model.ANM(model.current_structure).lambda(:,handles.mode+6);
scale=sqrt(Bfactor_scale/freq_sq);

amp=scale*handles.sign*handles.amplitude;


phases=180*linspace(0,1,handles.computation.res+1);
sin_trace=sin(pi*phases/180);
rmean=zeros(handles.pair_number,handles.computation.res+1);

axes(handles.axes_distr);
axis([0,6,-0.1*handles.distr_max,1.5*handles.distr_max]);

t=timer('TimerFcn',@(x,y)fprintf(1,''),'StartDelay',0.05);
for k=1:handles.computation.res+1,
    start(t);
    phax=phases(1:k);
    x=x0+mode(:,1)*sin_trace(k)*amp;
    y=y0+mode(:,2)*sin_trace(k)*amp;
    z=z0+mode(:,3)*sin_trace(k)*amp;
    for kk=1:handles.pair_number,
        i1=pair_pointers(kk,1);
        i2=pair_pointers(kk,2);
        adr1=handles.pairs(kk).adr1;
        adr2=handles.pairs(kk).adr2;
        [rax,act_distr]=pair_distribution(handles,adr1,adr2,Ca_coor,[x,y,z],i1,i2);
        set(handles.hdistr(kk),'XData',rax,'YData',act_distr);
        rmean(kk,k)=sum(rax.*act_distr)/sum(act_distr);
        xyz1=[x(i1),y(i1),z(i1)];
        xyz2=[x(i2),y(i2),z(i2)];
        % rmean(kk,k)=norm(xyz1-xyz2)/10;
        set(handles.hmean(kk),'XData',phax,'YData',rmean(kk,1:k));
    end;
    drawnow;
    wait(t);
end;
dmin=min(min(rmean));
dmax=max(max(rmean));

axes(handles.axes_mean);
axis tight;
v=axis;
axis([v(1),v(2),0.9*v(3),1.1*v(4)]);


function NOpos1=update_labels(NOpos0,network0,network1,residue)
% updates spin label or atom coordinates after a change of C_alpha
% coordinates of the network by applying local transformations
% NOpos0    n*4 coordinate set (xyz+population) before change
% network0  m*3 C_alpha coordinates before change
% network1  m*3 C_alpha coordinates after change
% residue   residue number (in network) corresponding to NOpos coordinates
% NOpos1    n*4 coordinate set after change

global model

cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network

local_template=zeros(5,3);
local_template_0=zeros(5,3);

[m,n]=size(network0);
xyz=NOpos0(:,1:3);
[n,n0]=size(xyz);
NOpos1=NOpos0;

% make a local template to fit rotation and translation
poi=0;
for kk=-2:2,
    if residue+kk>0 && residue+kk<=m, % is addressed residue a network point?
        diff=cindices(residue+kk,4)-cindices(residue,4);
        if diff==kk, % is addressed residue part of a continuous segment?
            poi=poi+1;
            local_template(poi,:)=network1(residue+kk,:);
            local_template_0(poi,:)=network0(residue+kk,:);
        end;
    end;
end;
poi=2;
if poi>=3, % found sufficient number of points to determine local rotation and translation
    [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
    xyz=[xyz ones(n,1)];
    xyz=transmat*xyz';
    xyz=xyz';
    NOpos1(:,1:3)=xyz(:,1:3);
else
    shift=network1(residue,:)-network0(residue,:);
    NOpos1(:,1:3)=xyz+repmat(shift,n,1);
end;

function [y1,p,rmsd,r]=lin_fit(x,y)

[y0,poi]=sort(y,'ascend');
x0=x(poi);
p=polyfit(x0,y0,1);
y1=polyval(p,x);
rmsd=y1-y;
rmsd=sqrt(sum(rmsd.*rmsd)/(length(rmsd)-1));
r0=corrcoef(x,y);
r=r0(1,2);



function edit_experimental_Callback(hObject, eventdata, handles)
% hObject    handle to edit_experimental (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_experimental as text
%        str2double(get(hObject,'String')) returns contents of edit_experimental as a double

[v,handles]=edit_update_MMM(handles,hObject,0.1,10,2.5,'%4.2f',0);
handles.experimental_distance=v;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_experimental_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_experimental (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_experimental.
function checkbox_experimental_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_experimental (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_experimental

if get(hObject,'Value'),
    set(handles.edit_experimental,'ForegroundColor',[0,0,0]);
else
    set(handles.edit_experimental,'ForegroundColor',[0.8,0.8,0.8]);
end;

guidata(hObject,handles);
