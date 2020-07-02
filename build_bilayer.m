function varargout = build_bilayer(varargin)
% BUILD_BILAYER M-file for build_bilayer.fig
%      BUILD_BILAYER, by itself, creates a new BUILD_BILAYER or raises the existing
%      singleton*.
%
%      H = BUILD_BILAYER returns the handle to a new BUILD_BILAYER or the handle to
%      the existing singleton*.
%
%      BUILD_BILAYER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BUILD_BILAYER.M with the given input arguments.
%
%      BUILD_BILAYER('Property','Value',...) creates a new BUILD_BILAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before build_bilayer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to build_bilayer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help build_bilayer

% Last Modified by GUIDE v2.5 12-Jan-2018 12:39:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @build_bilayer_OpeningFcn, ...
                   'gui_OutputFcn',  @build_bilayer_OutputFcn, ...
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


% --- Executes just before build_bilayer is made visible.
function build_bilayer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to build_bilayer (see VARARGIN)

% Choose default command line output for build_bilayer
handles.output = hObject;

global model
% global MMM_icon
global hMain

handles.probe_radius=1.9; % based on L. Adamian, V. Nanda, W. F. De Grado, J. Liang
                  % Proteins, 2005, 59:496-509

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

snum=model.current_structure;

handles.width=60;
handles.layer=35;
handles.z=0;
handles.theta=0;
handles.phi=0;

if isfield(model.info{snum},'bilayer'),
    if ~isempty(model.info{snum}.bilayer.graphics),
        for k=1:3,
            if ishandle(model.info{snum}.bilayer.graphics(k)),
                delete(model.info{snum}.bilayer.graphics(k));
            end;
        end;
        model.info{snum}.bilayer.graphics=[];
    end;
    if ~isempty(model.info{snum}.bilayer.width)
        handles.width=model.info{snum}.bilayer.width;
    end;
    if ~isempty(model.info{snum}.bilayer.thickness)
        handles.layer=model.info{snum}.bilayer.thickness;
    end;
end;

set(handles.slider_midpoint,'Min',handles.z-10);
set(handles.slider_midpoint,'Max',handles.z+10);
set(handles.slider_midpoint,'Value',handles.z);
set(handles.slider_thickness,'Value',handles.layer);
set(handles.edit_midpoint,'String',sprintf('%4.1f',handles.z));
set(handles.edit_thickness,'String',sprintf('%4.1f',handles.layer));

set(handles.checkbox_outer,'Value',1);
set(handles.checkbox_center,'Value',1);

handles=make_planes(handles);

[handles.TM_indices,handles.propensities]=get_TM_residues;

[fname,access]=accessibility_report(handles.probe_radius,handles.TM_indices);

handles.accessibility=access(:,6);

handles=display_energy(handles);

set(handles.uibuttongroup_energy,'SelectionChangeFcn',@selcbk);

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes build_bilayer wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = build_bilayer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

snum=model.current_structure;

model.info{snum}.bilayer.graphics=handles.gobjects;
model.info{snum}.bilayer.width=handles.width;
model.info{snum}.bilayer.thickness=handles.layer;
model.info{snum}.bilayer.show_outer=get(handles.checkbox_outer,'Value');
model.info{snum}.bilayer.show_center=get(handles.checkbox_center,'Value');

if abs(handles.phi)>1e-3 || abs(handles.theta)>1e-3 || abs(handles.z)>1e-2,
    transmat=affine('Euler',[handles.phi,handles.theta,0]);
    transmat2=affine('translation',[0,0,-handles.z]);
    transform_structure(model.current_structure,{transmat,transmat2},true);
    hMain.graph_update=1;
end;

delete(handles.figure1);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for k=1:3,
    delete(handles.gobjects(k));
end;

delete(handles.figure1);

% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'build_bilayer.html');
webcall(entry,'-helpbrowser');


% --- Executes on button press in checkbox_outer.
function checkbox_outer_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_outer

hb=handles.gobjects(2);
ht=handles.gobjects(3);

if get(hObject,'Value'),
    set(hb,'Visible','on');
    set(ht,'Visible','on');
else
    set(hb,'Visible','off');
    set(ht,'Visible','off');
end;

guidata(hObject,handles);


% --- Executes on button press in checkbox_center.
function checkbox_center_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_center

hc=handles.gobjects(1);

if get(hObject,'Value'),
    set(hc,'Visible','on');
else
    set(hc,'Visible','off');
end;

guidata(hObject,handles);

% --- Executes on button press in pushbutton_expand.
function pushbutton_expand_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.width=sqrt(2)*handles.width;
for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);

guidata(hObject,handles),

% --- Executes on button press in pushbutton_shrink.
function pushbutton_shrink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shrink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.width=handles.width/sqrt(2);
for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);

guidata(hObject,handles),

function handles=make_planes(handles)

global hMain
global graph_settings

x=[-handles.width/2 handles.width/2 handles.width/2 -handles.width/2];
y=[-handles.width/2 -handles.width/2 handles.width/2 handles.width/2];
z=handles.z*ones(1,4);

transmat=affine('Euler',[0,-handles.theta,-handles.phi]);

for k=1:4,
    coor=[x(k) y(k) z(k)];
    coor=affine_trafo_point(coor,transmat);
    x(k)=coor(1);
    y(k)=coor(2);
    z(k)=coor(3);
end;

axes(hMain.axes_model);
hc=patch(x,y,z,graph_settings.bilayer_midplane_color);
hb=patch(x,y,z-handles.layer/2*ones(1,4),graph_settings.bilayer_outer_color);
ht=patch(x,y,z+handles.layer/2*ones(1,4),graph_settings.bilayer_outer_color);

set(hc,'FaceAlpha',0.25);
set(hb,'FaceAlpha',0.25);
set(ht,'FaceAlpha',0.25);
handles.gobjects=[hc hb ht];

if get(handles.checkbox_center,'Value'),
    set(hc,'Visible','on');
else
    set(hc,'Visible','off');
end;

if get(handles.checkbox_outer,'Value'),
    set(hb,'Visible','on');
    set(ht,'Visible','on');
else
    set(hb,'Visible','off');
    set(ht,'Visible','off');
end;


% --- Executes on slider movement.
function slider_midpoint_Callback(hObject, eventdata, handles)
% hObject    handle to slider_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.z=get(hObject,'Value');
set(handles.edit_midpoint,'String',sprintf('%4.1f',handles.z));

for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);
handles=display_energy(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_midpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_midpoint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_midpoint as text
%        str2double(get(hObject,'String')) returns contents of edit_midpoint as a double

[v,handles]=edit_update_MMM(handles,hObject,-100,100,0,'%4.1f',0);

handles.z=v;

set(handles.slider_midpoint,'Min',v-10);
set(handles.slider_midpoint,'Max',v+10);
set(handles.slider_midpoint,'Value',v);

for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);
handles=display_energy(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_midpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to slider_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.layer=get(hObject,'Value');
set(handles.edit_thickness,'String',sprintf('%4.1f',handles.layer));

for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);
handles=display_energy(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thickness as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness as a double

[v,handles]=edit_update_MMM(handles,hObject,30,50,40,'%4.1f',0);

handles.layer=v;

set(handles.slider_thickness,'Value',v);

for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);
handles=display_energy(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [TM_indices,propensities]=get_TM_residues

global model
global residue_defs

TM_exist=0;

propensities=[];

snum=model.current_structure;
sadr=mk_address(snum);
chains=length(model.structures{snum});

TM_poi=0;


sadr=mk_address(model.current_structure);
adr=sprintf('%s(:){1}',sadr);
cindices=resolve_address(adr);
[msg,TM_indices]=get_object(cindices,'children');
if iscell(TM_indices),
    TM_indices=vertcat(TM_indices{:});
end;

[TM_poi,n]=size(TM_indices);

propensities=zeros(TM_poi,2);

for k=1:TM_poi,
    myindices=TM_indices(k,:);
    [message,tag]=get_residue(myindices,'name');
    if strcmpi(tag,'R1A'),
        tag='ILE';
    end;
    resnum=tag2id(tag,upper(residue_defs.restags));
    if ~isempty(resnum),
        propensities(k,1)=residue_defs.residues(resnum).TMLIPH;
        propensities(k,2)=residue_defs.residues(resnum).TMLIPC;
        propensities(k,3)=residue_defs.residues(resnum).TMBH;
        propensities(k,4)=residue_defs.residues(resnum).TMBC;
    end;
end;

function handles=display_energy(handles)

propensities=[];
if get(handles.radiobutton_bundle,'Value'),
    propensities=handles.propensities(:,1:2);
elseif get(handles.radiobutton_barrel,'Value')
    propensities=handles.propensities(:,3:4);
end;

if isempty(propensities),
    set(handles.text_energy,'String','n.a.');
    add_msg_board('Warning: Energy function not defined.');
    return;
end;

[G,msg]=G_propensity(handles.TM_indices,handles.z,handles.layer,handles.accessibility,propensities);
if get(handles.checkbox_fit_normal,'Value'),
    [G,msg]=G_propensity_trafo(handles.TM_indices,handles.z,handles.layer,handles.accessibility,propensities,handles.theta,handles.phi);
end;
if ~isempty(G),
    set(handles.text_energy,'String',sprintf('%6.2f',G));
else
    set(handles.text_energy,'String','n.a.');
    add_msg_board(msg.text);
    drawnow;
end;


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propensities=[];
if get(handles.radiobutton_bundle,'Value'),
    propensities=handles.propensities(:,1:2);
elseif get(handles.radiobutton_barrel,'Value')
    propensities=handles.propensities(:,3:4);
end;

if isempty(propensities),
    set(handles.text_energy,'String','n.a.');
    add_msg_board('ERROR: Fit impossible. Energy function not defined.');
    return;
end;

set(gcf,'Pointer','watch');
add_msg_board('Performing grid search');
drawnow;
z_grid=linspace(-10,10,21);
d_grid=linspace(25,50,21);
if get(handles.checkbox_fit_thickness,'Value')==0,
    z_grid=linspace(-10,10,101);
    d_grid=handles.layer;
end;
minG=1e6;
opt_z=handles.z;
opt_d=handles.layer;
for k=1:length(z_grid),
    if mod(k,2)==0,
        add_msg_board(sprintf('%i%% completed, best energy is %4.1f kcal/mol',round(100*(k-1)/length(z_grid)),minG));
        add_msg_board(sprintf('optimum midpoint is %4.1f ?, optimum thickness %4.1f ?',opt_z,opt_d));
        drawnow;
    end;
    z=z_grid(k);
    for kk=1:length(d_grid),
        d=d_grid(kk);
        G=G_propensity(handles.TM_indices,z,d,handles.accessibility,propensities);
        if G<minG,
            minG=G;
            opt_z=z;
            opt_d=d;
        end;
    end;
end;
add_msg_board('--- Bilayer grid search for midpoint and thickness completed ---');
add_msg_board(sprintf('Preliminary midpoint z coordinate: %5.2f',opt_z));
add_msg_board(sprintf('Preliminary bilayer thickness    : %5.2f',opt_d));
drawnow
if get(handles.checkbox_fit_normal,'Value'),
    add_msg_board('Now doing grid search for membrane normal direction');
    th_grid=linspace(-pi/4,pi/4,21);
    phi_grid=linspace(-pi,pi,41);
    minG=1e6;
    opt_theta=handles.theta;
    opt_phi=handles.phi;
    for k=1:length(th_grid),
        add_msg_board(sprintf('%i%% completed, best energy is %4.1f kcal/mol',round(100*(k-1)/length(th_grid)),minG));
        add_msg_board(sprintf('optimum theta is %4.1f?, optimum phi %4.1f?',180*opt_theta/pi,180*opt_phi/pi));
        drawnow;
        theta=th_grid(k);
        for kk=1:length(phi_grid),
            phi=phi_grid(kk);
            G=G_propensity_trafo(handles.TM_indices,opt_z,opt_d,handles.accessibility,propensities,theta,phi);
            if G<minG,
                minG=G;
                opt_theta=theta;
                opt_phi=phi;
            end;
        end;
    end;
    add_msg_board('--- Bilayer grid search for membrane normal direction completed ---');
    add_msg_board(sprintf('Preliminary membrane normal theta angle: %4.1f?',180*opt_theta/pi));
    add_msg_board(sprintf('Preliminary membrane normal phi angle  : %4.1f?',180*opt_phi/pi));
    drawnow
end;
v0=[opt_z,opt_d];
if get(handles.checkbox_fit_normal,'Value'),
    if get(handles.checkbox_fit_thickness,'Value'),
        v0=[v0 opt_theta opt_phi];
        add_msg_board('Fitting midpoint, thickness, and membrane normal direction...');
    else
        v0=[opt_z opt_theta opt_phi];
        add_msg_board('Fitting midpoint and membrane normal direction...');
    end;
elseif get(handles.checkbox_fit_thickness,'Value'),
    add_msg_board('Fitting midpoint and thickness...');
else
    v0=opt_z;
    add_msg_board('Fitting midpoint...');    
end;
options=optimset('fminsearch');
options=optimset(options,'TolFun',0.05,'TolX',0.05);
drawnow;
[v1,Gmin]=fminsearch(@G_fit,v0,options,handles.TM_indices,handles.accessibility,propensities,handles);
add_msg_board('--- Bilayer fit completed ---');
if Gmin>0,
    add_msg_board('WARNING: Positive transfer energy indicates inappropriate propensities.');
    add_msg_board('         Please also check whether membrane normal could be wrong.');
end;
switch length(v1)
    case 1
        handles.z=v1;
        add_msg_board(sprintf('Optimum midpoint z coordinate: %5.2f',v1(1)));
    case 2
        handles.z=v1(1);
        handles.layer=v1(2);
        add_msg_board(sprintf('Optimum midpoint z coordinate: %5.2f',v1(1)));
        add_msg_board(sprintf('Optimum bilayer thickness    : %5.2f',v1(2)));
    case 3
        handles.z=v1(1);
        handles.theta=v1(2);
        handles.phi=v1(3);
        add_msg_board(sprintf('Optimum midpoint z coordinate: %5.2f',v1(1)));
        add_msg_board(sprintf('Optimum membrane normal angle theta: %4.1f?',180*v1(2)/pi));
        add_msg_board(sprintf('Optimum membrane normal angle phi  : %4.1f?',180*v1(3)/pi));
    case 4
        handles.z=v1(1);
        handles.layer=v1(2);
        handles.theta=v1(3);
        handles.phi=v1(4);
        add_msg_board(sprintf('Optimum midpoint z coordinate: %5.2f',v1(1)));
        add_msg_board(sprintf('Optimum bilayer thickness    : %5.2f',v1(2)));
        add_msg_board(sprintf('Optimum membrane normal angle theta: %4.1f?',180*v1(3)/pi));
        add_msg_board(sprintf('Optimum membrane normal angle phi  : %4.1f?',180*v1(4)/pi));
end;
set(gcf,'Pointer','arrow');

set(handles.slider_midpoint,'Min',handles.z-10);
set(handles.slider_midpoint,'Max',handles.z+10);
set(handles.slider_midpoint,'Value',handles.z);
set(handles.slider_thickness,'Value',handles.layer);
set(handles.edit_midpoint,'String',sprintf('%4.1f',handles.z));
set(handles.edit_thickness,'String',sprintf('%4.1f',handles.layer));
set(handles.text_theta,'String',sprintf('%4.1f',180*handles.theta/pi));
set(handles.text_phi,'String',sprintf('%4.1f',180*handles.phi/pi));

for k=1:3,
    delete(handles.gobjects(k));
end;
handles=make_planes(handles);
handles=display_energy(handles);

guidata(hObject,handles);

function G=G_fit(v,indices,accessibility,propensities,handles)

switch length(v)
    case 1
        G=G_propensity(indices,v,handles.layer,accessibility,propensities);
    case 2
        G=G_propensity(indices,v(1),v(2),accessibility,propensities);
    case 3
        G=G_propensity_trafo(indices,v(1),handles.layer,accessibility,propensities,v(2),v(3));
    case 4
        G=G_propensity_trafo(indices,v(1),v(2),accessibility,propensities,v(3),v(4));
end;


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_barrel.
function radiobutton_barrel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_barrel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=display_energy(handles);
guidata(hObject,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_bundle.
function radiobutton_bundle_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_bundle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=display_energy(handles);
guidata(hObject,handles);

function selcbk(source,eventdata)

handles=guidata(source);
handles=display_energy(handles);
guidata(source,handles);


% --- Executes on button press in checkbox_fit_normal.
function checkbox_fit_normal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fit_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fit_normal


% --- Executes on button press in checkbox_fit_thickness.
function checkbox_fit_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fit_thickness
