function varargout = assign_TM(varargin)
% ASSIGN_TM M-file for assign_TM.fig
%      ASSIGN_TM, by itself, creates a new ASSIGN_TM or raises the existing
%      singleton*.
%
%      H = ASSIGN_TM returns the handle to a new ASSIGN_TM or the handle to
%      the existing singleton*.
%
%      ASSIGN_TM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGN_TM.M with the given input arguments.
%
%      ASSIGN_TM('Property','Value',...) creates a new ASSIGN_TM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before assign_TM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to assign_TM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help assign_TM

% Last Modified by GUIDE v2.5 05-Jan-2010 12:31:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @assign_TM_OpeningFcn, ...
                   'gui_OutputFcn',  @assign_TM_OutputFcn, ...
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


% --- Executes just before assign_TM is made visible.
function assign_TM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to assign_TM (see VARARGIN)

% Choose default command line output for assign_TM
handles.output = hObject;

global model
% global MMM_icon
global hMain

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

[handles,success,TM_exist]=setup_secondary(handles);

if ~success,
    add_msg_board('ERROR: No secondary structure elements are defined.');
    delete(hObject);
end;

snum=model.current_structure;
adr=mk_address(snum);
set(hObject,'Name',sprintf('Assign TM helices or strands for structure %s',adr));

if TM_exist,
    set(handles.pushbutton_normal,'Enable','on');
else
    set(handles.pushbutton_normal,'Enable','off');
end;

handles.gobjects=[];

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes assign_TM wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = assign_TM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'assign_TM.html');
webcall(entry,'-helpbrowser');


% --- Executes on selection change in listbox_all_secondary.
function listbox_all_secondary_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_all_secondary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_all_secondary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_all_secondary


% --- Executes during object creation, after setting all properties.
function listbox_all_secondary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_all_secondary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_TM.
function listbox_TM_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_TM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_TM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_TM


% --- Executes during object creation, after setting all properties.
function listbox_TM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_TM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_pick.
function pushbutton_pick_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TM_indices=get(handles.listbox_TM,'UserData');
[nTM,n]=size(TM_indices);
TM_list=get(handles.listbox_TM,'String');
poi=length(TM_list);
if poi>0 && isempty(TM_list{poi}),
    poi=0;
end;
poi0=poi;

sec_list=get(handles.listbox_all_secondary,'String');
sec_indices=get(handles.listbox_all_secondary,'UserData');
sel=get(handles.listbox_all_secondary,'Value');

if ~isempty(sel),
    for k=1:length(sel),
        cindices=sec_indices(sel(k),:);
        found=0;
        if nTM>0,
            for kk=1:nTM,
                if sum(abs(cindices-TM_indices(kk,:)))==0,
                    add_msg_board(sprintf('Double selection of %s ignored.',TM_list{kk}));
                    found=1;
                    break;
                end;
            end;
        end;
        if ~found,
            poi=poi+1;
            TM_list{poi}=sec_list{sel(k)};
            if poi==1,
                TM_indices=sec_indices(sel(k),:);
            else
                TM_indices=[TM_indices;sec_indices(sel(k),:)];
            end;
            assign_it(sec_indices(sel(k),:),1);
        end;
    end;
end;

if poi>poi0,
    set(handles.listbox_TM,'String',TM_list);
    set(handles.listbox_TM,'UserData',TM_indices(1:poi,:));
    set(handles.listbox_TM,'Value',1);
end;

if poi>0,
    set(handles.pushbutton_normal,'Enable','on');
else
    set(handles.pushbutton_normal,'Enable','off');
end;

guidata(hObject,handles);

% --- Executes on button press in pushbutton_unpick.
function pushbutton_unpick_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_unpick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TM_indices=get(handles.listbox_TM,'UserData');
[nTM,n]=size(TM_indices);
TM_list=get(handles.listbox_TM,'String');
poi=length(TM_list);
if poi>0 && isempty(TM_list{poi}),
    poi=0;
end;
poi0=poi;

sel=get(handles.listbox_TM,'Value');

TM_indices0=TM_indices;
poi=0;
if ~isempty(sel),
    for k=1:poi0,
        for kk=1:length(sel),
            if k~=sel(kk),
                poi=poi+1;
                TM_indices(poi,:)=TM_indices0(k,:);
                new_TM_list{poi}=TM_list{k};
            else
                assign_it(TM_indices(k,:),0);
            end;
        end;
    end;
end;

if poi>0,
    set(handles.pushbutton_normal,'Enable','on');
    set(handles.listbox_TM,'UserData',TM_indices(1:poi,:));
    set(handles.listbox_TM,'String',new_TM_list);
    set(handles.listbox_TM,'Value',1);
else
    set(handles.pushbutton_normal,'Enable','off');
    set(handles.listbox_TM,'UserData',[]);
    set(handles.listbox_TM,'String',{''});
    set(handles.listbox_TM,'Value',1);
end;

guidata(hObject,handles);

% --- Executes on button press in pushbutton_normal.
function pushbutton_normal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if ~isempty(handles.gobjects),
    delete(handles.gobjects);
end;

handles.gobjects=[];

TM_indices=get(handles.listbox_TM,'UserData');
[nTM,n]=size(TM_indices);

axes(hMain.axes_model);

mvec=zeros(1,3);
mp0=zeros(1,3);
for k=1:nTM,
    [p0,vec]=get_director(TM_indices(k,:));
    if k>1,
        nvec2=vec/norm(vec);
        if sum(nvec1.*nvec2)<0,
            vec=-vec;
        end;
    else
        nvec1=vec/norm(vec);
    end;
    mp0=mp0+p0;
    mvec=mvec+vec;
    h=plot3([p0(1)-vec(1)/2,p0(1)+vec(1)/2],[p0(2)-vec(2)/2,p0(2)+vec(2)/2],[p0(3)-vec(3)/2,p0(3)+vec(3)/2],'r','LineWidth',1);
    handles.gobjects=[handles.gobjects h];
end;
p0=mp0/nTM;
vec=mvec/nTM;
nvec=vec/norm(vec);
add_msg_board(sprintf('Approximate membrane normal has a tilt of %5.1f? with respect to current z axis.',180*acos(nvec(3))/pi));
h=plot3([p0(1)-vec(1)/2,p0(1)+vec(1)/2],[p0(2)-vec(2)/2,p0(2)+vec(2)/2],[p0(3)-vec(3)/2,p0(3)+vec(3)/2],'r','LineWidth',2);
handles.gobjects=[handles.gobjects h];

guidata(hObject,handles);


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if ~isempty(handles.gobjects),
    delete(handles.gobjects);
end;

TM_indices=get(handles.listbox_TM,'UserData');
[nTM,n]=size(TM_indices);

p0=[];
vec=[];

if nTM>0,
    mvec=zeros(1,3);
    mp0=zeros(1,3);
    for k=1:nTM,
        [p0,vec]=get_director(TM_indices(k,:));
        if k>1,
            nvec2=vec/norm(vec);
            if sum(nvec1.*nvec2)<0,
                vec=-vec;
            end;
        else
            nvec1=vec/norm(vec);
        end;
        mp0=mp0+p0;
        mvec=mvec+vec;
    end;
    p0=mp0/nTM;
    vec=mvec/nTM;
end;

hMain.origin=p0;
hMain.z_axis=vec;

delete(handles.figure1);


% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if ~isempty(handles.gobjects),
    delete(handles.gobjects);
end;

hMain.origin=[];
hMain.z_axis=[];

delete(handles.figure1);

% --- Executes on button press in pushbutton_unpick_all.
function pushbutton_unpick_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_unpick_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TM_indices=get(handles.listbox_TM,'UserData');
[nTM,n]=size(TM_indices);
for k=1:nTM,
    assign_it(TM_indices(k,:),0);
end;

set(handles.pushbutton_normal,'Enable','off');
set(handles.listbox_TM,'UserData',[]);
set(handles.listbox_TM,'String',{''});
set(handles.listbox_TM,'Value',1);

guidata(hObject,handles);

function [handles,success,TM_exist]=setup_secondary(handles)

global model

success=0;
TM_exist=0;

snum=model.current_structure;
sadr=mk_address(snum);
chains=length(model.structures{snum});

poi=0;
TM_poi=0;

sec_indices=zeros(200,4);
TM_indices=zeros(200,4);

% list all helices
if chains>0,
    for cnum=1:chains,
        sadr=mk_address([snum cnum]);
        if isfield(model.structures{snum}(cnum),'helix_defs') && ~isempty(model.structures{snum}(cnum).helix_defs),
            for kk=1:length(model.structures{snum}(cnum).helix_defs),
                range=model.structures{snum}(cnum).helix_defs{kk}.range;
                poi=poi+1;
                sec_list{poi}=sprintf('%s<H.%s>; %i-%i',sadr,model.structures{snum}(cnum).helix_defs{kk}.name,range(1),range(2));
                sec_indices(poi,:)=[snum cnum kk 1];
                success=1;
                if isfield(model.structures{snum}(cnum).helix_defs{kk},'TM'),
                    if model.structures{snum}(cnum).helix_defs{kk}.TM,
                        TM_poi=TM_poi+1;
                        TM_list{TM_poi}=sprintf('%s<H.%s>; %i-%i',sadr,model.structures{snum}(cnum).helix_defs{kk}.name,range(1),range(2));
                        TM_indices(TM_poi,:)=[snum cnum kk 1];
                        TM_exist=1;
                    end;
                end;
            end;
        end;
    end;
end;
    
% list all strands
if chains>0,
    for cnum=1:chains,
        if isfield(model.structures{snum}(cnum),'sheet_defs') && ~isempty(model.structures{snum}(cnum).sheet_defs),
            for kk=1:length(model.structures{snum}(cnum).sheet_defs),
                range=model.structures{snum}(cnum).sheet_defs{kk}.range;
                poi=poi+1;
                sec_list{poi}=sprintf('%s<E.%s>; %i-%i',sadr,model.structures{snum}(cnum).sheet_defs{kk}.name,range(1),range(2));
                sec_indices(poi,:)=[snum cnum kk 2];
                success=1;
                if isfield(model.structures{snum}(cnum).sheet_defs{kk},'TM'),
                    if model.structures{snum}(cnum).sheet_defs{kk}.TM,
                        TM_poi=TM_poi+1;
                        TM_list{TM_poi}=sprintf('%s<E.%s>; %i-%i',sadr,model.structures{snum}(cnum).sheet_defs{kk}.name,range(1),range(2));
                        TM_indices(TM_poi,:)=[snum cnum kk 2];
                        TM_exist=1;
                    end;
                end;
            end;
        end;
    end;
end;

if success,
    set(handles.listbox_all_secondary,'String',sec_list);
    set(handles.listbox_all_secondary,'UserData',sec_indices(1:poi,:));
    set(handles.listbox_all_secondary,'Value',1);
end;

if TM_exist,
    set(handles.listbox_TM,'String',TM_list);
    set(handles.listbox_TM,'UserData',TM_indices(1:TM_poi,:));
    set(handles.listbox_TM,'Value',1);
else
    set(handles.listbox_TM,'String',{''});
    set(handles.listbox_TM,'UserData',[]);
    set(handles.listbox_TM,'Value',1);
end;

function assign_it(indices,flag)
% assign a helix or sheet as transmembrane or other
%
% indices   internal indices of helix or sheet
% flag      1 transmembrane, 0 other

global model
global graph_settings

adr=mk_address(indices(1:2));

switch indices(4)
    case 1
        range=model.structures{indices(1)}(indices(2)).helix_defs{indices(3)}.range;
        model.structures{indices(1)}(indices(2)).helix_defs{indices(3)}.TM=flag;
        if flag,
            col=graph_settings.TM_helix_color;
        else
            col=graph_settings.helix_color;
        end;
    case 2
        range=model.structures{indices(1)}(indices(2)).sheet_defs{indices(3)}.range;
        model.structures{indices(1)}(indices(2)).sheet_defs{indices(3)}.TM=flag;
        if flag,
            col=graph_settings.TM_sheet_color;
        else
            col=graph_settings.sheet_color;
        end;
end;

adr=sprintf('%s{:}%i-%i',adr,range(1),range(2));
indices=resolve_address(adr);
set_object(indices,'ribboncolor',{col});

function [p0,vec]=get_director(indices)

global model

p0=[];
vec=[];

adr=mk_address(indices(1:2));

switch indices(4)
    case 1
        range=model.structures{indices(1)}(indices(2)).helix_defs{indices(3)}.range;
    case 2
        range=model.structures{indices(1)}(indices(2)).sheet_defs{indices(3)}.range;
end;

adr=sprintf('%s{:}%i-%i',adr,range(1),range(2));
indices=resolve_address(adr);
[m,n]=size(indices);
if m==0,
    add_msg_board('Warning: Index mismatch for helix or strand director.');
else
    coor=zeros(m,3);
    for k=1:m,
        adr=mk_address(indices(k,:));
        newadr=sprintf('%s.CA',adr);
        [msg,xyz]=get_object(newadr,'coor');
        coor(k,:)=xyz;
    end;
    [p0,vec]=rmsd_line_3D(coor);
end;


% --- Executes on button press in pushbutton_hide_normal.
function pushbutton_hide_normal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_hide_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.gobjects),
    delete(handles.gobjects);
end;

handles.gobjects=[];

guidata(hObject,handles);
