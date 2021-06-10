function varargout = fit_from_template(varargin)
% fit_from_template M-file for fit_from_template.fig
%      fit_from_template, by itself, creates a new fit_from_template or raises the existing
%      singleton*.
%
%      H = fit_from_template returns the handle to a new fit_from_template or the handle to
%      the existing singleton*.
%
%      fit_from_template('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in fit_from_template.M with the given input arguments.
%
%      fit_from_template('Property','Value',...) creates a new fit_from_template or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fit_from_template_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fit_from_template_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fit_from_template

% Last Modified by GUIDE v2.5 31-Aug-2015 18:04:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fit_from_template_OpeningFcn, ...
                   'gui_OutputFcn',  @fit_from_template_OutputFcn, ...
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


% --- Executes just before fit_from_template is made visible.
function fit_from_template_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fit_from_template (see VARARGIN)

% global MMM_icon
global hMain
global model
global ENM_param

handles.current_structure=model.current_structure;

min_lambda=1e3*eps; % minimum frequency of a non-degenerate mode
handles.darken=0.75; % darken chain colorscheme colors for better visibility and match with lighted ribbon model

% animation control parameters
handles.animation.fps=15;
handles.animation.amp=1;
handles.animation.res=handles.animation.fps;
handles.animation.cycles=5;
handles.animation.scaling=3;
handles.animation.rotation=1;

% fit parameters

handles.ensemble=1;
handles.uncertainty=0;
handles.exclude=true;

% Choose default command line output for fit_from_template
handles.output = hObject;

set(hMain.figure,'Pointer','watch');

drawnow;

h = msgbox('Please be patient. This can take several minutes.','Anisotropic network model is computed');

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];
hMain.fit_plot=true;
hMain.fit_axes=handles.axes_model;

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

snum=model.current_structure;
adr=mk_address(snum);
set(handles.figure1,'Name',sprintf('Fit transition from template structure %s',adr));

set(handles.text_ANM_type,'String',ENM_param.parametrization);
handles.stored_parametrization=ENM_param.parametrization;

mhandles=guidata(hMain.figure);
locked=get(mhandles.uitoggletool_lock,'State');

if strcmpi(locked,'on'),
    set(handles.checkbox_thermal,'Enable','off');
    set(handles.popupmenu_fit_mode,'Enable','off');
else
    set(handles.checkbox_thermal,'Enable','on');
    set(handles.popupmenu_fit_mode,'Enable','on');
end;

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


if ~isfield(model,'ANM') || length(model.ANM)<snum || isempty(model.ANM(snum).u),
    Hessian=setup_ANM_bonded(Ca_coor);
    % Hessian=setup_ANM_poly(Ca_coor);
    contacts=[];
    % [Hessian,contacts]=setup_ANM(Ca_coor,restypes,rindices); % residue-specific force constants, not recommended
    
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
    degenerate=sum(lambda<=min_lambda);
%    set(handles.text_info,'String',sprintf('%i modes are suspicious. First suspicious mode: %i',sum(suspicious)-degenerate,fs-degenerate));
else
    lambda=model.ANM(snum).lambda;
    m=model.ANM(snum).residues;
    msf=model.ANM(snum).msf;
end;

nu=sqrt(lambda(7:end));
inv_nu=ones(size(nu))./nu;
inv_nu=inv_nu/inv_nu(1);
handles.inv_nu=inv_nu;
[mi,handles.maxbas]=min(abs(inv_nu-ENM_param.tif));
set(handles.text_info,'String',sprintf('%i modes contribute with >=%4.1f%% weighting of slowest mode.',handles.maxbas,100*ENM_param.tif));

% Diagnostic plot on distribution of (inverse) mode frequencies
% figure(1); clf;
% modenum=1:length(inv_nu);
% plot(modenum,inv_nu,'r');
% hold on
% plot(modenum,ENM_param.tif*ones(size(modenum)),'r');
% axis([0,200,0,1.1]);

[msf,p]=lin_fit(msf,Bfactors);
handles.Bfactor_scale=p(1);

handles.degenerate=sum(lambda<=min_lambda);
handles.mode=1;

axes(handles.axes_plot);
plot(lambda,'k');
hold on;
ma=max(lambda(1+handles.degenerate:100+handles.degenerate));
plot([handles.mode+handles.degenerate,handles.mode+handles.degenerate],[-0.05*ma,1.05*ma],'b:');
axis([handles.degenerate+1,handles.degenerate+100,-0.075*ma,1.075*ma]);
xlabel('Mode number');
ylabel('Normalized frequency');
set(handles.text_auxiliary_msg,'String',sprintf('Norm. freq.: %6.4f',lambda(handles.mode+handles.degenerate)));

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

m=model.ANM(model.current_structure).residues;
evec=model.ANM(model.current_structure).u(:,handles.mode+handles.degenerate);
mode=reshape(evec,3,m);
handles.current_mode=mode;

if ishandle(h),
    delete(h);
end;

handles=set_animation_control(handles);

set(hMain.figure,'Pointer','arrow');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fit_from_template wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fit_from_template_OutputFcn(hObject, eventdata, handles) 
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

entry=strcat(help_files,'fit_from_template.html');
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

ENM_param.parametrization=handles.stored_parametrization;
ENM_param=set_ANM(ENM_param);

delete(handles.figure1);


function [sc,rmsd]=optimum_scale(y,x)
% Optimum scaling of data y_i, so that the root mean square deviation
% rmsd = sum((sc*y_i-x_i)^2) is minimized
% vectors x and y must have the same length
%
% y     vector of data to be scaled
% x     vector of data that should be fitted
% sc    scaling factor
% rmsd  minimized root mean square deviation
%
% G. Jeschke, 2010

sc=sum(x.*y)/sum(y.^2);
diff=sc*y-x;
rmsd=sqrt(sum(diff.^2)/(length(diff)-1));



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


% --- Executes on button press in pushbutton_animate.
function pushbutton_animate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_animate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
center=model.info{model.current_structure}.center;
cam_up=get(hMain.axes_model,'CameraUpVector');

axes(handles.axes_model);

M=animate(hObject,handles,true);

set(handles.text_info,'String','Playing animation...');
set(gcf,'Pointer','watch');

chains=length(handles.wire);

for c=1:chains,
    set(handles.wire(c),'visible','off');
end;
drawnow;

handles=query_radiobuttons(handles);

units=get(handles.axes_model,'Units');
set(handles.axes_model,'Units','pixels');
loc=get(handles.axes_model,'Position');
set(handles.axes_model,'Units',units);

if handles.animation.rotation==1,
    movie(gcf,M,handles.animation.cycles,handles.animation.fps,loc);
else
    movie(gcf,M,1,handles.animation.fps,loc);
end;

set(handles.axes_model,'CameraTarget',cam_tar);
set(handles.axes_model,'CameraPosition',cam_pos);
set(handles.axes_model,'CameraUpVector',cam_up);
camlookat(handles.axes_model);

for c=1:chains,
    set(handles.wire(c),'visible','on');
end;
set(handles.text_info,'String','Animation finished.');
set(gcf,'Pointer','arrow');
guidata(hObject,handles);

function F=animate(hObject,handles,record,full)

global model

set(handles.text_info,'String','Recording animation...');
set(gcf,'Pointer','watch');

if nargin<4,
    full=false;
end;

if ~record,
    F=[];
end;

axes(handles.axes_model);
if get(handles.checkbox_black,'Value'),
    bckg=get(handles.uipanel_animation,'BackgroundColor');
    frg=get(handles.uipanel_animation,'ForegroundColor');
    set(handles.uipanel_animation,'BackgroundColor','k');
    set(handles.uipanel_animation,'ForegroundColor','r');
end;

units=get(handles.axes_model,'Units');
set(handles.axes_model,'Units','pixels');
rect=get(handles.axes_model,'Position');
set(handles.axes_model,'Units',units);

handles=query_radiobuttons(handles);

cycles=handles.animation.cycles;
angle_step=360/(cycles*handles.animation.res+1);

Ca_coor=model.coarse(model.current_structure).Ca_coor;
chains=length(handles.wire);

if record && ~full && handles.animation.rotation==1,
    cycles=1;
end;

[mf,nf]=size(model.fit.dxmat);
frames=nf/3;

fr=0;
for k0=1:cycles,
    x=Ca_coor(:,1);
    y=Ca_coor(:,2);
    z=Ca_coor(:,3);
    for k=1:2*frames-1,
        fr=fr+1;
        if k>frames,
            kf=2*frames-k;
            sgn=-1;
        else
            kf=k;
            sgn=1;
        end;
        dx=sgn*model.fit.dxmat(:,3*kf-2:3*kf);
        switch handles.animation.rotation,
            case 2
                camorbit(handles.axes_model,0,angle_step);
            case 3
                camorbit(handles.axes_model,angle_step,0);
        end;
        for c=1:chains,
            ia=handles.chains(c,2);
            ie=handles.chains(c,3);
            x(ia:ie)=x(ia:ie)+dx(ia:ie,1);
            y(ia:ie)=y(ia:ie)+dx(ia:ie,2);
            z(ia:ie)=z(ia:ie)+dx(ia:ie,3);
            set(handles.wire(c),'XData',x(ia:ie),'YData',y(ia:ie),'ZData',z(ia:ie));
        end;
        drawnow;
        if record,
            F(fr)=getframe(gcf,rect);
        end;
    end;
end;

set(handles.axes_model,'DrawMode','normal');
if get(handles.checkbox_black,'Value'),
    set(handles.uipanel_animation,'BackgroundColor',bckg);
    set(handles.uipanel_animation,'ForegroundColor',frg);
end;
for c=1:chains,
    ia=handles.chains(c,2);
    ie=handles.chains(c,3);
    x=Ca_coor(ia:ie,1);
    y=Ca_coor(ia:ie,2);
    z=Ca_coor(ia:ie,3);
    set(handles.wire(c),'XData',x,'YData',y,'ZData',z);
end;
drawnow;
set(handles.text_info,'String','Animation recorded.');
set(gcf,'Pointer','arrow');
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
 

% --- Executes on button press in pushbutton_cycles_minus.
function pushbutton_cycles_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cycles_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.animation.cycles>1,
    handles.animation.cycles=handles.animation.cycles-1;
end;
set(handles.edit_cycles,'String',sprintf('%i',handles.animation.cycles));
guidata(hObject,handles);


% --- Executes on button press in pushbutton_cycles_plus.
function pushbutton_cycles_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cycles_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.animation.cycles<100,
    handles.animation.cycles=handles.animation.cycles+1;
end;
set(handles.edit_cycles,'String',sprintf('%i',handles.animation.cycles));
guidata(hObject,handles);



function edit_cycles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cycles as text
%        str2double(get(hObject,'String')) returns contents of edit_cycles as a double

[v,handles]=edit_update_MMM(handles,hObject,1,100,5,'%i',1);
handles.animation.cycles=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=set_animation_control(handles)

set(handles.edit_cycles,'String',sprintf('%i',handles.animation.cycles));
switch handles.animation.rotation,
    case 1
        set(handles.radiobutton_rot_none,'Value',1);
    case 2
        set(handles.radiobutton_rot_x,'Value',1);
    case 3
        set(handles.radiobutton_rot_z,'Value',1);
end;

function handles=query_radiobuttons(handles)

if get(handles.radiobutton_rot_none,'Value'),
    handles.animation.rotation=1;
end;
if get(handles.radiobutton_rot_x,'Value'),
    handles.animation.rotation=2;
end;
if get(handles.radiobutton_rot_z,'Value'),
    handles.animation.rotation=3;
end;


% --- Executes on button press in pushbutton_record.
function pushbutton_record_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain
global general

my_path=pwd;
cd(general.reports);

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
cam_up=get(hMain.axes_model,'CameraUpVector');

axes(handles.axes_model);

[filename, pathname] = uiputfile('*.avi', 'Save animation in AVI format');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Save as AVI cancelled by user');
else
    reset_user_paths(pathname);
    general.reports=pathname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Animation is saved as AVI file: %s',fname);
    add_msg_board(msg);
end

M=animate(hObject,handles,true,true);

snum=model.current_structure;
adr=mk_address(snum);
movie2avi(M,fname,'fps',handles.animation.fps,'videoname',sprintf('MMM animation of normal mode %i of %s',handles.mode+handles.degenerate,adr),'compression','Cinepak');

set(handles.axes_model,'CameraTarget',cam_tar);
set(handles.axes_model,'CameraPosition',cam_pos);
set(handles.axes_model,'CameraUpVector',cam_up);
camlookat(handles.axes_model);

set(handles.text_info,'String','Animation finished.');
set(gcf,'Pointer','arrow');
cd(my_path);
guidata(hObject,handles);

% --- Executes on button press in checkbox_initial_final.
function checkbox_initial_final_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_initial_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_initial_final

global model

if get(hObject,'Value'),
    col1=[0.6,0,0];
    col2=[0,0.6,0];
    Ca_coor=model.coarse(model.current_structure).Ca_coor;
    Ca_coor_end=model.fit.Ca_coor;
    axes(handles.axes_model);
    handles.initial_wire=line(Ca_coor(:,1),Ca_coor(:,2),Ca_coor(:,3),'color',col1,'LineWidth',1);
    handles.final_wire=line(Ca_coor_end(:,1),Ca_coor_end(:,2),Ca_coor_end(:,3),'color',col2,'LineWidth',1);
else
    set(handles.initial_wire,'Visible','off');
    set(handles.final_wire,'Visible','off');
end;
guidata(hObject,handles),

% --- Executes on button press in checkbox_black.
function checkbox_black_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_black (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_black


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% in generation of an ensemble, all direct C_alpha-C_alpha constraints are
% always used, only DEER constraints can be ommitted, a warning is issued
% if the number of DEER constraints is insufficient to generate the
% specified ensemble

global model
global ENM_param
global hMain
global general

requested_cpu=general.cpu;

thermal=get(handles.checkbox_thermal,'Value');

% if thermal,
%     store_cycles=ENM_param.cycles;
%     ENM_param.cycles=400;
% end;

correspondence=[];
target=[];
test_mode=false;
overlap_mode=false;
repack=false;

fit_mode=get(handles.popupmenu_fit_mode,'Value');
switch fit_mode,
    case 1
        ENM_param.diagonalize=false;
        ENM_param.reorientate=true;
        ENM_param.fix_local=2;
    case 2
        ENM_param.diagonalize=false;
        ENM_param.reorientate=false;
        ENM_param.fix_local=2;
    case 3
        ENM_param.diagonalize=true;
        ENM_param.reorientate=false;
        ENM_param.fix_local=2;
    case 4
        ENM_param.diagonalize=true;
        ENM_param.reorientate=false;
        ENM_param.fix_local=0;
    case 5
        ENM_param.diagonalize=true;
        ENM_param.reorientate=false;
        ENM_param.fix_local=0;
end;
if ~isempty(handles.target),
    [correspondence,target]=align_template_target(mk_address(model.current_structure),handles.target);
    if ~isempty(correspondence),
        test_mode=true;
    end;
end;
network0=model.coarse(model.current_structure).Ca_coor;

if test_mode && overlap_mode,
    [overlap,correlation,collectivity]=ANM_merit(network0(correspondence(1,:),:),target(correspondence(2,:),:));
    sumoverlap=sqrt(cumsum(overlap.^2));
    [maxoverlap,bestmode]=max(overlap);
    [mini,num80]=min(abs(sumoverlap-0.8));
    [mini,num90]=min(abs(sumoverlap-0.9));
    [mini,num95]=min(abs(sumoverlap-0.95));
    [mini,num99]=min(abs(sumoverlap-0.99));
    figure(13); clf;
    plot(overlap,'k');
    hold on;
    plot(sumoverlap,'r');
    axis([0,length(overlap)+1,-0.05,1.05]);
    disp(sprintf('Collectivity of structural change is %5.3f.',collectivity));
    disp(sprintf('Best mode %i has overlap %5.3f and correlation %5.3f with structural change.',bestmode,maxoverlap,correlation(bestmode)));
    disp(sprintf('%i modes cover 80%% of structural change.',num80));
    disp(sprintf('%i modes cover 90%% of structural change.',num90));
    disp(sprintf('%i modes cover 95%% of structural change.',num95));
    disp(sprintf('%i modes cover 99%% of structural change.',num99));
    disp(sprintf('10 modes cover %4.1f%% of structural change.',100*sumoverlap(10)));
    disp(sprintf('20 modes cover %4.1f%% of structural change.',100*sumoverlap(20)));
end;

if test_mode,
    rmsd0=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));
    fprintf(1,'Initial structure r.m.s.d. %4.2f ?\n',rmsd0);
else
    rmsd0=1e6;
end;

if ~isempty(handles.displacements) && (~isempty(handles.DEER) || ~isempty(handles.direct)),
    add_msg_board('Warning: Displacement constraints currently cannot be combined with other constraints.');
    add_msg_board('Only displacement constraints are used.');
end;

if handles.ensemble>1 && ~isempty(handles.displacements),
    add_msg_board('Warning: Ensemble generation not implemented for displacement restraints.');
    add_msg_board('Single structure will be generated.');
end;

% generate information for ensemble generation, if requested 
if handles.ensemble>1 && handles.exclude && isempty(handles.displacements),
    if isempty(handles.DEER) || length(handles.DEER)<2,
        add_msg_board('Warning: Ensemble can be generated only from DEER restraints.');
        add_msg_board('Only a single structure will be generated.');
        esize=1;
    else
        num_DEER=length(handles.DEER);
        DEER_diff=zeros(1,num_DEER); % vector of normalized distance restraint differences from template structure
        for k=1:num_DEER,
            r0=norm(handles.DEER(k).xyz2-handles.DEER(k).xyz1)/10;
            DEER_diff(k)=abs(handles.DEER(k).r-r0)/handles.DEER(k).sigr;
        end;
        [DEER_diff_sort,pointers]=sort(DEER_diff,'descend');
        % generate constraint index vectors for ensemble
        p=num_DEER;
        esize=1;
        constraint_indices=zeros(handles.ensemble,num_DEER);
        constraint_indices(1,:)=pointers; % first structure has all constraints
        maxomit=0;
        while p>1 && esize<handles.ensemble,
            p=p-1;
            missing=num_DEER-p; % number of missing constraints
            if missing>3,
                add_msg_board('Warning. No more than 3 DEER restraints can be omitted.');
                add_msg_board(sprintf('Ensemble size restricted to %i, whereas %i was requested.',esize,handles.ensemble));
                break;
            end;
            if p<2,
                add_msg_board('Warning. At least 2 DEER restraints must remain.');
                add_msg_board(sprintf('Ensemble size restricted to %i, whereas %i was requested.',esize,handles.ensemble));
                break;
            end;
            switch missing,
                case 1
                    maxomit=1;
                    pp=num_DEER;
                    while esize<handles.ensemble && pp>0,
                        esize=esize+1;
                        constraint_indices(esize,1:num_DEER-1)=[pointers(1:pp-1) pointers(pp+1:num_DEER)]; 
                        pp=pp-1;
                    end;
                case 2
                    maxomit=2;
                    citemp=zeros((num_DEER*(num_DEER-1))/2,num_DEER);
                    merit=zeros(1,(num_DEER*(num_DEER-1))/2);
                    poi=0;
                    for k=1:num_DEER-1,
                        for kk=k+1:num_DEER,
                            curr=[pointers(1:k-1) pointers(k+1:kk-1) pointers(kk+1:end)];
                            poi=poi+1;
                            citemp(poi,1:num_DEER-2)=curr;
                            merit(poi)=sum(DEER_diff(curr).^2);
                        end;
                    end;
                    [sortmerit,secpointers]=sort(merit,'descend');
                    newsets=handles.ensemble-esize;
                    if newsets<length(secpointers),
                        secpointers=secpointers(1:newsets);
                    end;
                    constraint_indices(esize+1:esize+length(secpointers),:)=citemp(secpointers,:);
                    esize=esize+length(secpointers);
                case 3
                    maxomit=3;
                    citemp=zeros(num_DEER*(num_DEER-1)*(num_DEER-2)/6,num_DEER);
                    merit=zeros(1,num_DEER*(num_DEER-1)*(num_DEER-2)/6);
                    poi=0;
                    for k=1:num_DEER-2,
                        for kk=k+1:num_DEER-1,
                            for kkk=kk+1:num_DEER,
                                curr=[pointers(1:k-1) pointers(k+1:kk-1) pointers(kk+1:kkk-1) pointers(kkk+1:end)];
                                poi=poi+1;
                                citemp(poi,1:num_DEER-3)=curr;
                                merit(poi)=sum(DEER_diff(curr).^2);
                            end;
                        end;
                    end;
                    [sortmerit,secpointers]=sort(merit,'descend');
                    newsets=handles.ensemble-esize;
                    if newsets<length(secpointers),
                        secpointers=secpointers(1:newsets);
                    end;
                    constraint_indices(esize+1:esize+length(secpointers),:)=citemp(secpointers,:);
                    esize=esize+length(secpointers);
            end;
        end;
        add_msg_board(sprintf('Maximum number of omitted restraints: %i. Ensemble size: %i',maxomit,esize));
        constraint_indices=constraint_indices(1:esize,:);
    end;
else
    constraint_indices=1:length(handles.DEER);
    if handles.exclude,
        esize=1;
    else
        esize=handles.ensemble;
        constraint_indices=repmat(constraint_indices,esize,1);
    end;
end;

direct=handles.direct;
stored_basis=ENM_param.fit_basis;
ENM_param.fit_basis=handles.basis;


all_rmsd=zeros(1,esize);
all_coverage=zeros(1,esize);
all_distance_rmsd=zeros(1,esize);
all_best_fit=zeros(1,esize);
all_best_drmsd=zeros(1,esize);
all_collectivity=zeros(1,esize);
fom_vec=zeros(1,esize);


parallel_flag=true;

exclude_mode=handles.exclude;
ivec0=1:length(handles.DEER);
DEER0=handles.DEER;
uncertainty=handles.uncertainty;
displacements=handles.displacements;
[mi,maxbas]=min(abs(handles.inv_nu-handles.tif));

my_current_structure=model.current_structure;

model_copy=model;

ENM_param_copy=ENM_param;

if parallel_flag, 
    hm=msgbox('Ensemble generation is running. This window will close on completion.','Parallel computation');
end;

options.maxbas=maxbas;
options.sat_cycle=handles.mmax;
options.overfit=get(handles.checkbox_overfit,'Value');
mb=maxbas;
if mb<handles.basis,
    mb=handles.basis;
end;

add_msg_board(sprintf('Active space extension from %i to %i modes in %i iterations.',handles.basis,mb,handles.mmax));

my_snum=model.current_structure;
model_numbers=zeros(1,esize);

set(hMain.figure,'Pointer','watch');

% ### change back to parfor ###
for k=1:esize, % loop for generating an ensemble
    c=clock;
    fprintf(1,'Started trial %i of %i trials at %i:%i:%i on %i.%i.%i.\n',k,esize,c(4),c(5),round(c(6)),c(3),c(2),c(1));
    etrial=k;
    if esize==1,
        etrial=0;
    end;
%     ivec=constraint_indices(k,:);
%     ivec=ivec(ivec>0);
    if exclude_mode && esize>1, % exclude some existing constraints for current model, if requested
        select=rand(1,length(ivec0));
        [ma,poi]=max(select);
        if poi==1,
            ivec=ivec0(2:end);
        elseif poi==length(ivec0),
            ivec=ivec0(1:end-1);
        else
            ivec=[ivec0(1:poi-1) ivec0(poi+1:end)];
        end;
        DEER=DEER0(ivec);
    else
        DEER=DEER0;
    end;
    if uncertainty>0, % add random deviations to distance constraints
        for kd=1:length(DEER),
            DEER(kd).r=DEER(kd).r+uncertainty*randn;
            % fprintf(fid,'Trial %k constraint: %8s%8s%5.2f\n',DEER(kd).adr1,DEER(kd).adr2,DEER(kd).r);
        end;
    end;
    add_msg_board(sprintf('Initializing fit for trial %i/%i.',k,esize));
    if isempty(displacements) && fit_mode~=5,
        if ~thermal,
            [fom,rmsd,network,dxmat,DEER,best,diagnostics] = fit_by_ANM(network0,DEER,direct,test_mode,correspondence,target,model_copy,ENM_param_copy);
        else
            [fom,rmsd,network,dxmat,DEER,best,diagnostics] = fit_by_ANM_transition(network0,DEER,direct,options,test_mode,correspondence,target,etrial,model_copy,ENM_param_copy);
        end;
    elseif ~isempty(displacements)
        [fom,rmsd,network,dxmat,best] = fit_displacements_by_ANM(network0,displacements,test_mode,correspondence,target);
        diagnostics=[];
    else
        [model_snum,fom,rmsd,network,dxmat,DEER,best,diagnostics] = fit_by_ANM_and_force_field(my_snum,network0,DEER,direct,options,test_mode,correspondence,target,etrial,model_copy,ENM_param_copy);
        model_numbers(k)=model_snum;
    end;
    if test_mode && isfield(diagnostics,'ys'),
        if ~isempty(diagnostics.ys),
            add_msg_board(sprintf('Coverage at last iteration: %5.3f',(rmsd0-diagnostics.ys(end))/rmsd0));
        end;
    end;
    fit_info{k}=diagnostics;
    DEERs{k}=DEER;
    networks{k}=network;
    dxmats{k}=dxmat;
    all_distance_rmsd(k)=rmsd;
    fom_vec(k)=fom;

    if test_mode,
        [rms,coor2]=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
        diff=network(correspondence(1,:),:)-network0(correspondence(1,:),:);
        [mn,nn]=size(diff);
        change=reshape(diff',1,3*mn);
        change=change';
        kappa_DR=mode_collectivity(change);
        all_rmsd(k)=rms;
        all_coverage(k)=(rmsd0-rms)/rmsd0;
        all_best_fit(k)=best.rmsd;
        all_best_drmsd(k)=best.drmsd;
        all_collectivity(k)=kappa_DR;
        % disp_local_change(target,network);
    end;

%     if repack,
%         newtag2=sprintf('tp%i',maxs+1);
%         curr_structures=curr_structures+1;
%         [snum1,infile]=repacked_copy(snum,newtag2,k);
%     end;
end;
snum=model.current_structure;
set(hMain.figure,'Pointer','arrow');

handles.fit_info=fit_info;
for k=1:esize,
    maxs=length(model.structures);
    newtag=sprintf('tr%i',maxs+1);
    if k==1,
        [snum,message]=copy_structure(model.current_structure,newtag,networks{k});
    else
        [snum1,message]=copy_structure(model.current_structure,newtag,networks{k},k,snum);
    end;
    if repack,
        newtag2=sprintf('tp%i',maxs+1);
        [snum1,infile]=repacked_copy(snum,newtag2,k);
    end;
end;

DEER=DEERs{end};
network=networks{end};
dxmat=dxmats{end};
fom=fom_vec(end);
rmsd=all_distance_rmsd(end);
handles.network_fit.Ca_coor=network;
handles.network_fit.fom=fom;
handles.network_fit.DEER=DEER;
handles=update_3D_model(handles,network,DEER);

ENM_param.fit_basis=stored_basis;
model.fit.Ca_coor=network;
model.fit.rmsd=rmsd;
model.fit.fom=fom;
model.fit.dxmat=dxmat;
fid=fopen(sprintf('%s_ensemble_20.log',handles.target),'w');
fprintf(fid,'Ensemble logfile for target %s\n\n',handles.target);

for k=1:esize,
    fprintf(fid,'\n*** Trial %i ***\n\n',k);
    DEER=DEERs{k};
    fprintf(fid,'Constraints\n');
    for kd=1:length(DEER),
        fprintf(fid,'%12s %12s %5.2f\n',DEER(kd).adr1,DEER(kd).adr2,DEER(kd).r); 
    end;
    add_msg_board(sprintf('Trial %i: New figure of merit : %6.4f',k,fom_vec(k)));
    add_msg_board(sprintf('Trial %i: New distance r.m.s.d: %6.4f',k,all_distance_rmsd(k)));
    if test_mode,
        add_msg_board(sprintf('Trial %i: Initial r.m.s.d. w.r.t. target structure: %5.2f ?',k,rmsd0));
        add_msg_board(sprintf('Trial %i: Final   r.m.s.d. w.r.t. target structure: %5.2f ?',k,all_rmsd(k)));
        add_msg_board(sprintf('Trial %i: Coverage of the structural change: %5.3f\n',k,all_coverage(k)));
        fprintf(fid,'\nFinal   r.m.s.d. w.r.t. target structure: %5.2f ?\n',all_rmsd(k));
        fprintf(fid,'Distance constraint r.m.s.d.     : %5.2f\n',all_distance_rmsd(k));
        fprintf(fid,'Coverage of the structural change: %5.3f\n',all_coverage(k));
    end;
end;

if test_mode,
    if esize>1,
        fprintf(fid,'\n\nMean final   r.m.s.d. w.r.t. target structure  : %5.2f ?\n',mean(all_rmsd));
        fprintf(fid,'Std. dev. of r.m.s.d. w.r.t. target structure  : %5.2f ?\n',std(all_rmsd));
        fprintf(fid,'Mean distance constraint r.m.s.d.              : %5.2f ?\n',mean(all_distance_rmsd));
        fprintf(fid,'Std. dev. of distance constraint r.m.s.d.      : %5.2f ?\n',std(all_distance_rmsd));
        fprintf(fid,'Mean coverage of the structural change         : %5.3f\n',mean(all_coverage));
        fprintf(fid,'Std. dev. of coverage of the structural change : %5.3f\n',std(all_coverage));
        fprintf(fid,'Worst        coverage of the structural change : %5.3f\n',min(all_coverage));
        fprintf(fid,'Best         coverage of the structural change : %5.3f\n',max(all_coverage));
        fprintf(fid,'Mean structure rmsd for best approach to target: %5.2f\n',mean(all_best_fit));
        fprintf(fid,'Mean distance rmsd for best approach to target : %5.2f\n',mean(all_best_drmsd));
        fprintf(fid,'Mean collectivity of the fitted change         : %5.3f\n',mean(all_collectivity));
        fprintf(fid,'Std. dev. of collectivity of the fitted change : %5.3f\n\n',std(all_collectivity));
        fprintf(fid,'Extracted fit    : Mean r.m.s.d. w.r.t. target : %5.2f ?\n',mean(all_rmsd));
        fprintf(fid,'                 : Max. r.m.s.d. w.r.t. target : %5.2f ?\n',max(all_rmsd));

        diff=all_coverage-mean(all_coverage);
        [~,poi]=min(abs(diff));
        fprintf(fid,'The most typical result was obtained in trial %i.\n',poi);
        fprintf(fid,'Deviation of this result from mean coverage is: %5.3f\n\n',diff(poi));
    end;
    fclose(fid);
end;

if test_mode,
    finame=sprintf('%s_ensemble_20_fit_info',handles.target);
    save(finame,'fit_info');
end;


% if repack
%     handles.current_structure=snum1;
% else
%     handles.current_structure=snum;
% end;

if parallel_flag,
    delete(hm);
    add_msg_board('Parallel ensemble computation is finished.');
end;

set(handles.pushbutton_save,'Enable','on');
guidata(hObject,handles);

% --- Executes on button press in pushbutton_fit.
function pushbutton_parametrize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% in generation of an ensemble, all direct C_alpha-C_alpha constraints are
% always used, only DEER constraints can be ommitted, a warning is issued
% if the number of DEER constraints is insufficient to generate the
% specified ensemble

global model
global ENM_param
global hMain
global general

bas_levels=[0.20,0.15,0.1,0.075,0.05];
nb=length(bas_levels);
bas_sizes=zeros(size(bas_levels));

for k=1:nb,
    [mi,bas_sizes(k)]=min(abs(handles.inv_nu-bas_levels(k)));
end;

sat_cycles=[15,20,25,30,35];
ns=length(sat_cycles);

bas_sizes=repmat(bas_sizes,ns,1);
sat_cycles=repmat(sat_cycles',1,nb);

ind1=1:nb;
ind2=1:ns;
ind1=repmat(ind1,1,ns);
ind2=repmat(ind2',1,nb);
ind2=reshape(ind2',1,ns*nb);

esize=ns*nb;

requested_cpu=general.cpu;

fprintf(1,'--- Parametrization grid ---\n');
for k=1:length(ind1), 
    fprintf(1,'Basis size: %i. Sat. cycle %i.\n',bas_sizes(ind1(k),ind2(k)),sat_cycles(ind1(k),ind2(k))); 
    options.maxbas=bas_sizes(ind1(k),ind2(k));
    options.sat_cyc=sat_cycles(ind1(k),ind2(k));
    option_vec{k}=options;
end;

thermal=true;

store_cycles=ENM_param.cycles;
ENM_param.cycles=400;

correspondence=[];
target=[];

ENM_param.diagonalize=true;
ENM_param.reorientate=false;
ENM_param.fix_local=0;

test_mode=false;
if ~isempty(handles.target),
    [correspondence,target]=align_template_target(mk_address(model.current_structure),handles.target);
    if ~isempty(correspondence),
        test_mode=true;
    end;
end;
network0=model.coarse(model.current_structure).Ca_coor;

if test_mode,
    rmsd0=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));
    fprintf(1,'Initial structure r.m.s.d. %4.2f ?\n',rmsd0);
end;

direct=handles.direct;
stored_basis=ENM_param.fit_basis;
ENM_param.fit_basis=handles.basis;

all_rmsd=zeros(1,esize);
all_coverage=zeros(1,esize);
all_distance_rmsd=zeros(1,esize);
all_best_fit=zeros(1,esize);
all_best_drmsd=zeros(1,esize);
all_collectivity=zeros(1,esize);
all_limit_fit=zeros(1,esize);
all_final_coverage=zeros(1,esize);
fom_vec=zeros(1,esize);

parallel_flag  = false;
if requested_cpu>1 && esize>1
    parallel_flag=true;
end;

exclude_mode=handles.exclude;
ivec0=1:length(handles.DEER);
DEER0=handles.DEER;
uncertainty=handles.uncertainty;
displacements=handles.displacements;
maxbas=handles.maxbas;

my_current_structure=model.current_structure;

model_copy=model;

ENM_param_copy=ENM_param;

if parallel_flag, 
    hm=msgbox('Ensemble generation is running. This window will close on completion.','Parallel computation');
end;

options.maxbas=maxbas;
DEER=handles.DEER;

parfor k=1:esize, % loop for generating an ensemble
    c=clock;
    fprintf(1,'Started trial %i of %i trials at %i:%i:%i on %i.%i.%i.\n',k,esize,c(4),c(5),round(c(6)),c(3),c(2),c(1));
    etrial=k;
    if esize==1,
        etrial=0;
    end;
    [fom,rmsd,network,dxmat,~,best,diagnostics] = fit_by_ANM_transition(network0,DEER,direct,option_vec{k},test_mode,correspondence,target,etrial,model_copy,ENM_param_copy);
    if isfield(diagnostics,'ys') && test_mode,
        all_final_coverage(k)=(rmsd0-diagnostics.ys(end))/rmsd0;
    end;
    fit_info{k}=diagnostics;
    DEERs{k}=DEER;
    networks{k}=network;
    dxmats{k}=dxmat;
    all_distance_rmsd(k)=rmsd;
    fom_vec(k)=fom;

    if test_mode,
        [rms,coor2]=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
        diff=network(correspondence(1,:),:)-network0(correspondence(1,:),:);
        [mn,nn]=size(diff);
        change=reshape(diff',1,3*mn);
        change=change';
        kappa_DR=mode_collectivity(change);
        all_rmsd(k)=rms;
        all_coverage(k)=(rmsd0-rms)/rmsd0;
        all_best_fit(k)=best.rmsd;
        all_best_drmsd(k)=best.drmsd;
        all_collectivity(k)=kappa_DR;
        if ~isempty(handles.DEER),
            all_limit_fit(k)=best.limit_rmsd;
        end;
    end;

end;
snum=model.current_structure;

handles.fit_info=fit_info;
for k=1:esize,
    maxs=length(model.structures);
    newtag=sprintf('tr%i',maxs+1);
    if k==1,
        [snum,message]=copy_structure(model.current_structure,newtag,networks{k});
    else
        [snum1,message]=copy_structure(model.current_structure,newtag,networks{k},k,snum);
    end;
end;

DEER=DEERs{end};
network=networks{end};
dxmat=dxmats{end};
fom=fom_vec(end);
rmsd=all_distance_rmsd(end);
handles.network_fit.Ca_coor=network;
handles.network_fit.fom=fom;
handles.network_fit.DEER=DEER;
handles=update_3D_model(handles,network,DEER);

ENM_param.fit_basis=stored_basis;
model.fit.Ca_coor=network;
model.fit.rmsd=rmsd;
model.fit.fom=fom;
model.fit.dxmat=dxmat;
fid=fopen(sprintf('%s_parametrization_50b20_edENM.log',handles.target),'w');
fprintf(fid,'Parametrization file for target %s\n\n',handles.target);

fid2=fopen(sprintf('%s_parametrization_50b20_edENM.dat',handles.target),'w');
fprintf(fid2,'%% (basis size, saturation cycle)\n%%');
for k=1:length(ind1),
    fprintf(fid2,' (%3i,%3i)',bas_sizes(ind1(k),ind2(k)),sat_cycles(ind1(k),ind2(k))); 
end;
fprintf(fid2,'\n');
for k=1:length(ind1),
    fprintf(fid2,'%10.3f',all_coverage(k)); 
end;
fprintf(fid2,'\n\n');
fclose(fid2);

fid2=fopen(sprintf('%s_parametrization_50b20_final_edENM.dat',handles.target),'w');
fprintf(fid2,'%% (basis size, saturation cycle)\n%%');
for k=1:length(ind1),
    fprintf(fid2,' (%3i,%3i)',bas_sizes(ind1(k),ind2(k)),sat_cycles(ind1(k),ind2(k))); 
end;
fprintf(fid2,'\n');
for k=1:length(ind1),
    fprintf(fid2,'%10.3f',all_final_coverage(k)); 
end;
fprintf(fid2,'\n\n');
fclose(fid2);

for k=1:esize,
    fprintf(fid,'\n*** Trial %i ***\n\n',k);
    fprintf(fid,'%% (basis size, saturation cycle) = (%3i,%3i)\n',bas_sizes(ind1(k),ind2(k)),sat_cycles(ind1(k),ind2(k)));
    add_msg_board(sprintf('Trial %i: New figure of merit : %6.4f',k,fom_vec(k)));
    add_msg_board(sprintf('Trial %i: New distance r.m.s.d: %6.4f',k,all_distance_rmsd(k)));
    if test_mode,
        add_msg_board(sprintf('Trial %i: Initial r.m.s.d. w.r.t. target structure: %5.2f ?',k,rmsd0));
        add_msg_board(sprintf('Trial %i: Final   r.m.s.d. w.r.t. target structure: %5.2f ?',k,all_rmsd(k)));
        add_msg_board(sprintf('Trial %i: Coverage of the structural change: %5.3f\n',k,all_coverage(k)));
        fprintf(fid,'\nFinal   r.m.s.d. w.r.t. target structure: %5.2f ?\n',all_rmsd(k));
        fprintf(fid,'Distance constraint r.m.s.d.     : %5.2f\n',all_distance_rmsd(k));
        fprintf(fid,'Coverage of the structural change: %5.3f\n',all_coverage(k));
    end;
end;

fclose(fid);

if thermal,
    ENM_param.cycles=store_cycles;
end;

if parallel_flag,
    delete(hm);
    add_msg_board('Parallel ensemble computation is finished.');
end;

set(handles.pushbutton_save,'Enable','on');
guidata(hObject,handles),


function [y1,p,rmsd,r]=lin_fit(x,y)

[y0,poi]=sort(y,'ascend');
x0=x(poi);
p=polyfit(x0,y0,1);
y1=polyval(p,x);
rmsd=y1-y;
rmsd=sqrt(sum(rmsd.*rmsd)/(length(rmsd)-1));
r0=corrcoef(x,y);
r=r0(1,2);


% --- Executes on button press in pushbutton_step.
function pushbutton_step_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

modes=8+handles.degenerate;
coeff=3;
network0=model.coarse(model.current_structure).Ca_coor;
DEER=handles.DEER;
fom0=handles.fom;
grad=zeros(1,100);
sgn=zeros(1,100);
for k=1:100,
    modes=k+handles.degenerate;
    coeff=2;
    [fom1,network,DEER] = propagate_ANM(coeff,modes,network0,DEER);
    [fom2,network,DEER] = propagate_ANM(-coeff,modes,network0,DEER);
    grad(k)=fom0-min([fom1,fom2]);
    if fom1<fom2,
        sgn(k)=coeff;
    else
        sgn(k)=-coeff;
    end;
end;
[ma,poi]=max(grad);
modes=poi+handles.degenerate;
coeff=sgn(poi);
[fom,network,DEER] = propagate_ANM(coeff,modes,network0,DEER);
figure(2); clf;
plot(grad,'k.');
disp(sprintf('New figure of merit: %6.4f',fom));

handles=update_3D_model(handles,network,DEER);

guidata(hObject,handles);

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global ENM_param
global general

snum=model.current_structure;
axes(handles.axes_model);

my_path=pwd;
cd(general.restraint_files);

[fname,pname,findex]=uigetfile('*.dat','Load restraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Restraint loading cancelled by user');
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    handles.basis=ENM_param.fit_basis;
    restraints=rd_restraints(fullfile(pname,fname));
    handles.ensemble=restraints.ensemble;
    handles.uncertainty=restraints.uncertainty;
    handles.exclude=restraints.exclude;
    handles.target=restraints.target;
    if isfield(restraints,'basis'),
        handles.basis=restraints.basis;
    end;
    if isfield(restraints,'tif'),
        handles.tif=restraints.tif;
    else
        handles.tif=ENM_param.tif;
    end;
    if isfield(restraints,'mmax'),
        handles.mmax=restraints.mmax;
    else
        handles.mmax=ENM_param.mmax;
    end;
    [DEER,cancelled]=process_DEER_restraints(restraints);
    if cancelled,
        add_msg_board('Processing of DEER constraints cancelled.');
        return
    end;
    handles.DEER=DEER;
    [direct,cancelled]=process_direct_restraints(restraints);
    if cancelled,
        add_msg_board('Processing of direct constraints cancelled.');
        return
    end;
    handles.direct=direct;
    [displacements,cancelled]=process_displacement_restraints(restraints);
    if cancelled,
        add_msg_board('Processing of displacement constraints cancelled.');
        return
    end;
    handles.displacements=displacements;
    set(handles.checkbox_restraints,'Enable','on');
    set(handles.checkbox_restraints,'Value',1);
    c=length(handles.wire);
    for k=1:c,
        set(handles.wire(k),'Color',[0.25,0.25,0.25]);
    end;
    handles.restraint_graphics=hgtransform;
    if ~isempty(DEER),
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
                ind1 = cindices(1,:);
                if handles.DEER(k).type(1) == 1
                    ipoi = 1;
                elseif handles.DEER(k).type(1) == 2
                    ind1(4) = round((cindices(1,4)+cindices(2,4))/2);
                    ipoi = 2;
                end
                diff=ind1-model.coarse(snum).indices(l,:);
                if sum(abs(diff))==0,
                    handles.DEER(k).res1=l;
                    f1=true;
                    x=[handles.DEER(k).xyz1(1) model.coarse(snum).Ca_coor(l,1)];
                    y=[handles.DEER(k).xyz1(2) model.coarse(snum).Ca_coor(l,2)];
                    z=[handles.DEER(k).xyz1(3) model.coarse(snum).Ca_coor(l,3)];
                    handles.DEER(k).rl1=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
                end;
                ind2 = cindices(ipoi+1,:);
                if length(handles.DEER(k).type) > 1 && handles.DEER(k).type(2) == 2
                    ind2(4) = round((cindices(ipoi+1,4)+cindices(ipoi+2,4))/2);
                end
                diff=ind2-model.coarse(snum).indices(l,:);
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
                add_msg_board(sprintf('Warning: Residue for first label of restraint %i not in network model.',k));
                min_approach = 1e6;
                closest = 0;
                atom_xyz = handles.DEER(k).xyz1;
                for l2=1:length(model.coarse(snum).indices)
                    dxyz = norm(atom_xyz - model.coarse(snum).Ca_coor(l2,:));
                    if dxyz < min_approach
                        min_approach = dxyz;
                        closest = l2;
                    end
                end
                handles.DEER(k).res1=closest;
                x=[handles.DEER(k).xyz1(1) model.coarse(snum).Ca_coor(closest,1)];
                y=[handles.DEER(k).xyz1(2) model.coarse(snum).Ca_coor(closest,2)];
                z=[handles.DEER(k).xyz1(3) model.coarse(snum).Ca_coor(closest,3)];
                handles.DEER(k).rl2=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
                add_msg_board(sprintf('Approximating by the closest Calpha atom at %5.1f %s',min_approach,197));
            end;
            if ~f2,
                add_msg_board(sprintf('Warning: Residue for second label of restraint %i not in network model.',k));
                min_approach = 1e6;
                closest = 0;
                atom_xyz = handles.DEER(k).xyz2;
                for l2=1:length(model.coarse(snum).indices)
                    dxyz = norm(atom_xyz - model.coarse(snum).Ca_coor(l2,:));
                    if dxyz < min_approach
                        min_approach = dxyz;
                        closest = l2;
                    end
                end
                handles.DEER(k).res2=closest;
                x=[handles.DEER(k).xyz2(1) model.coarse(snum).Ca_coor(closest,1)];
                y=[handles.DEER(k).xyz2(2) model.coarse(snum).Ca_coor(closest,2)];
                z=[handles.DEER(k).xyz2(3) model.coarse(snum).Ca_coor(closest,3)];
                handles.DEER(k).rl2=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
                add_msg_board(sprintf('Approximating by the closest Calpha atom at %5.1f %s',min_approach,197));
            end;
        end;
        dvec=dvec/sqrt(sum(sum(dvec.^2)));
    elseif ~isempty(handles.direct),
        dvec=zeros(length(handles.direct),3);
        network0=model.coarse(snum).Ca_coor;
        for k=1:length(handles.direct),
            xyz1=network0(handles.direct(k,1),:);
            xyz2=network0(handles.direct(k,2),:);
            line(xyz1(1),xyz1(2),xyz1(3),...
                'Parent',handles.restraint_graphics,'Marker','.','Color','b');
            line(xyz2(1),xyz2(2),xyz2(3),...
                'Parent',handles.restraint_graphics,'Marker','.','Color','b');
            x=[xyz1(1) xyz2(1)];
            y=[xyz1(2) xyz2(2)];
            z=[xyz1(3) xyz2(3)];
            r0=norm(xyz1-xyz2)/10;
            dvec(k,:)=(handles.direct(k,4)-r0)*(xyz1-xyz2)/(10*r0);
            det=abs(r0-handles.direct(k,4))/handles.direct(k,5);
            if det>2,
                col='r';
            elseif det>1,
                col=[255,190,0]/255; % yellow-orange
            else
                col='g';
            end;
            line(x,y,z,'Parent',handles.restraint_graphics,'Color',col,'LineWidth',1.5,'LineStyle',':');
        end;
    end;
end;
if ~isempty(DEER),
    overlaps=zeros(1,100);
    for k=1:100,
        evec=model.ANM(model.current_structure).u(:,k+6);
        m=length(evec)/3;
        mode=reshape(evec,3,m);
        mdvec=zeros(length(handles.DEER),3);
        for kk=1:length(handles.DEER),
            dchange=mode(:,handles.DEER(kk).res1)-mode(:,handles.DEER(kk).res2);
            mdvec(kk,:)=dchange';
        end;
        mdvec=mdvec/sqrt(sum(sum(mdvec.^2)));
        overlaps(k)=abs(sum(sum(dvec.*mdvec)));
    end;
    % figure(1); clf;
    % plot(overlaps,'k.');
    % figure(2); clf;
    % plot(cumsum(overlaps),'k');
    axes(handles.axes_plot);
    cla;
    set(handles.text_auxiliary,'String','Constraint matching');
    ma=0;
    rmsd=0;
    fom=0; % figure of merit
    for k=1:length(handles.DEER),
        r0=norm(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/10;
        rmsd=rmsd+(r0-handles.DEER(k).r)^2;
        det=abs(r0-handles.DEER(k).r)/handles.DEER(k).sigr;
        fom=fom+det^2;
        if det>2,
            col='r';
        elseif det>1,
            col=[255,190,0]/255; % yellow-orange
        else
            col='g';
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
end;
set(handles.pushbutton_fit,'Enable','on');
set(handles.pushbutton_parametrize,'Enable','on');
cd(my_path);
guidata(hObject,handles);

function [DEER,cancelled]=process_DEER_restraints(restraints)

global model
global hMain

cancelled=false;

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

snum=model.current_structure;

if isfield(restraints,'PDB'),
    if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
        button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
        if strcmp(button,'No'),
            cancelled=true;
            DEER=[];
            return
        end;
    end;
end;

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
type_list = zeros(1,200);
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
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('ERROR: Restraint %i has first label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            add_msg_board('Processing of DEER constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
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
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                labelname = complabel(restraints.DEER(k).label,1); 
                label_list{poi} = labelname;
                if length(restraints.DEER(k).type) == 1
                    type_list(poi,:) = restraints.DEER(k).type;
                else
                    type_list(poi,:) = restraints.DEER(k).type(1);
                end
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s (%i) at site %s will be generated.',restraints.DEER(k).label,type_list(poi),adr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        ind2=resolve_address(adr2);
        if isempty(ind2),
            add_msg_board(sprintf('ERROR: Restraint %i has second label at site %s',k,adr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            add_msg_board('Processing of DEER restraints cancelled');
            cancelled=true;
            DEER=[];
            return;
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
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                labelname = complabel(restraints.DEER(k).label,2); 
                label_list{poi} = labelname;
                T_list(poi)=restraints.DEER(k).T;
                if length(restraints.DEER(k).type) == 1
                    type_list(poi) = restraints.DEER(k).type;
                else
                    type_list(poi) = restraints.DEER(k).type(2);
                end
                add_msg_board(sprintf('Rotamers for label %s (%i) at site %s will be generated.',restraints.DEER(k).label,type_list(poi),adr2));
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),complabel(restraints.DEER(k).label,1))
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            labelname = complabel(restraints.DEER(k).label,1);
            label_list{poi} = labelname;
            if length(restraints.DEER(k).type) == 1
                type_list(poi) = restraints.DEER(k).type;
            else
                type_list(poi) = restraints.DEER(k).type(1);
            end
            T_list(poi)=restraints.DEER(k).T;
            add_msg_board(sprintf('Rotamers for label %s (%i) at site %s will be generated.',restraints.DEER(k).label,type_list(poi),adr1));
        end;
        adr2=restraints.DEER(k).adr2;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),complabel(restraints.DEER(k).label,2)),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr2;
            labelname = complabel(restraints.DEER(k).label,2);
            label_list{poi} = labelname;
            T_list(poi)=restraints.DEER(k).T;
            if length(restraints.DEER(k).type) == 1
                type_list(poi) = restraints.DEER(k).type;
            else
                type_list(poi) = restraints.DEER(k).type(2);
            end
            add_msg_board(sprintf('Rotamers for label %s (%i) at site %s will be generated.',restraints.DEER(k).label,type_list(poi),adr2));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' ')
        if type_list(k) == 1
            command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
            hMain.store_undo=false;
            cmd(hMain,command);
        elseif type_list(k) == 2
            compadr = to_do_list{k};
            spoi = strfind(compadr,'|');
            adr1 = compadr(1:spoi-1);
            adr2 = compadr(spoi+1:end);
            command=sprintf('bilabel %s %s %s %i',adr1,adr2,label_list{k});
            hMain.store_undo=false;
            cmd(hMain,command);
        elseif type_list(k) == 3
            indices = resolve_address(to_do_list{k});
            atom_adr = sprintf('%s.%s',to_do_list{k},label_list{k});
            [~,xyz] = get_object(atom_adr,'coor');
            j = length(model.sites);
            res = length(model.sites{j}(1).residue) + 1;
            model.sites{j}(1).residue(res).label = label_list{k};
            model.sites{j}(1).residue(res).T = 298;
            model.sites{j}(1).residue(res).indices = indices;
            model.sites{j}(1).residue(res).NOpos = [xyz 1];
            model.sites{j}(1).residue(res).rotamers(1).pop = 1;
            model.sites{j}(1).residue(res).rotamers(1).coor = xyz;
        end
    end;
end;

if isfield(model,'sites')
    labels=label_information(model.sites);
else
    labels = [];
end
if isfield(model,'bisites')
    bilabels = bilabel_information(model.bisites);
else
    bilabels = [];
end

for k=1:length(restraints.DEER),
    adr1=restraints.DEER(k).adr1;
    ind1=resolve_comp_address(adr1);
    [m1,~] = size(ind1);
    adr2=restraints.DEER(k).adr2;
    ind2=resolve_comp_address(adr2);
    [m2,~] = size(ind2);
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).indices=[ind1;ind2];
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
    DEER(k).type = restraints.DEER(k).type;
    f1=false;
    f2=false;
    for l=1:length(labels)
        if m1 == 1
            diff1=ind1-labels(l).indices;
            if sum(abs(diff1))==0,
                f1=true;
                DEER(k).xyz1=labels(l).xyz;
                DEER(k).rmsd1=labels(l).rmsd;
            end;
        end
        if m2 == 1
            diff2=ind2-labels(l).indices;
            if sum(abs(diff2))==0,
                f2=true;
                DEER(k).xyz2=labels(l).xyz;
                DEER(k).rmsd2=labels(l).rmsd;
            end;
        end
    end;
    for l=1:length(bilabels),
        if m1 == 2
            diff1a=ind1-bilabels(l).indices;
            diff1b=flipud(ind1)-bilabels(l).indices;
            if sum(sum(abs(diff1a)))==0 || sum(sum(abs(diff1b)))==0
                f1=true;
                DEER(k).xyz1=bilabels(l).xyz;
                DEER(k).rmsd1=bilabels(l).rmsd;
            end
        end
        if m2 == 2
            diff2a=ind2-bilabels(l).indices;
            diff2b=flipud(ind2)-bilabels(l).indices;
            if sum(sum(abs(diff2a)))==0 || sum(sum(abs(diff2b)))==0
                f2=true;
                DEER(k).xyz2=bilabels(l).xyz;
                DEER(k).rmsd2=bilabels(l).rmsd;
            end;
        end
    end;
    if ~f1 || ~f2,
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gjeschke@ethz.ch');
        cancelled=true;
        DEER=[];
        return;
    end;
end;


function [direct,cancelled]=process_direct_restraints(restraints)

global model

cancelled=false;
if ~isfield(restraints,'direct'),
    direct=[];
    return;
end;

snum=model.current_structure;

md=length(restraints.direct);
direct=zeros(md,5);

cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);
for k=1:md,
    adr1=restraints.direct(k).adr1;
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Specified residue %s does not exist in template structure.',adr1));
        cancelled=true;
        direct=[];
        return;
    end;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    direct(k,1)=net1;
    adr2=restraints.direct(k).adr2;
    ind2=resolve_address(adr2);
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    direct(k,2)=net2;    
    direct(k,4)=restraints.direct(k).r;
    direct(k,5)=restraints.direct(k).sigr;
    if net1==0 || net2==0,
        add_msg_board('ERROR: Constrained C_alpha-C_alpha distance does not exist in network model.');
        cancelled=true;
        direct=[];
        return;
    end;
end;

function [displacements,cancelled]=process_displacement_restraints(restraints)

global model

cancelled=false;
if ~isfield(restraints,'displace'),
    displacements=[];
    return;
end;

snum=model.current_structure;

md=length(restraints.displace);
displacements=zeros(md,5);

cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);
for k=1:md,
    adr=restraints.displace(k).adr;
    ind1=resolve_address(adr);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Specified residue %s does not exist in template structure.',adr));
        cancelled=true;
        displacements=[];
        return;
    end;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    displacements(k,1)=net1;
    displacements(k,2:4)=restraints.displace(k).vec;
    displacements(k,5)=restraints.displace(k).weight;
    if net1==0,
        add_msg_board('ERROR: Constrained displacement does not exist in network model.');
        cancelled=true;
        displacements=[];
        return;
    end;
end;
weights=displacements(:,5);
displacements(:,5)=weights/sum(weights);

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
            if ~isempty(id)
                labels(poi).name=label_defs.residues(id).short_name;
            else
                labels(poi).name = sites{k0}(k1).residue(k).label;
            end
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd = NOpos_rmsd(NOpos);
        end;
    end;
end;

function bilabels = bilabel_information(bisites)

global model
global ligand_libraries

poi=0;
for k0=1:length(bisites),
    for k1=1:length(bisites{k0}.sites),
        poi=poi+1;
        bilabels(poi).indices = bisites{k0}.sites{k1}.indices;
        for kl = 1:length(ligand_libraries)
            if strcmpi(bisites{k0}.sites{k1}.label,ligand_libraries(kl).label) || strcmpi(bisites{k0}.sites{k1}.label,ligand_libraries(kl).tc)
                bilabels(poi).name=ligand_libraries(kl).tc;
            end
        end
        bilabels(poi).T = 298;
        coor = model.bisites{k0}.sites{k1}.coor;
        x=sum(coor(:,1).*coor(:,4));
        y=sum(coor(:,2).*coor(:,4));
        z=sum(coor(:,3).*coor(:,4));
        bilabels(poi).xyz=[x y z];
        bilabels(poi).rmsd = model.bisites{k0}.sites{k1}.rmsd;
    end
end

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
if rmsd < 1e-2
    rmsd = 1e-2;
end

% --- Executes on button press in checkbox_restraints.
function checkbox_restraints_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_restraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_restraints

if get(hObject,'Value'),
    set(handles.restraint_graphics,'Visible','on');
else
    set(handles.restraint_graphics,'Visible','off');
end;
guidata(hObject,handles);

function handles = update_3D_model(handles,network,DEER)

global model 

snum=model.current_structure;

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
        x=[handles.DEER(k).xyz1(1) network(handles.DEER(k).res1,1)];
        y=[handles.DEER(k).xyz1(2) network(handles.DEER(k).res1,2)];
        z=[handles.DEER(k).xyz1(3) network(handles.DEER(k).res1,3)];
        set(handles.DEER(k).rl2,'XData',x,'YData',y,'ZData',z);
    end;
    if ~f2,
        add_msg_board(sprintf('Warning: Residue for second label of restraint %i not in network model.',k));
        x=[handles.DEER(k).xyz1(1) network(handles.DEER(k).res2,1)];
        y=[handles.DEER(k).xyz1(2) network(handles.DEER(k).res2,2)];
        z=[handles.DEER(k).xyz1(3) network(handles.DEER(k).res2,3)];
        set(handles.DEER(k).rl2,'XData',x,'YData',y,'ZData',z);
    end;
end;
axes(handles.axes_plot);
cla;
set(handles.text_auxiliary,'String','Constraint matching');
ma=0;
rmsd=0;
fom=0; % figure of merit
if ~isempty(DEER),
    for k=1:length(DEER),
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

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.reports);

[filename, pathname] = uiputfile('*.mat', 'Save fitted network model in internal format');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Save of fitted network model model cancelled by user');
else
    reset_user_paths(pathname);
    general.reports=pathname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Fitted network model model is saved as .mat file: %s',fname);
    add_msg_board(msg);
end

network_fit=handles.network_fit;
save(fname,'network_fit');
cd(my_path);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_compare.
function pushbutton_compare_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global target

test_mode=true;
overlap_mode=true;

load 1anf_coarse_grained
target=coarse.Ca_coor;
% target=target(1:end-2,:);

network0=model.coarse(model.current_structure).Ca_coor;
% network0=network0(3:end,:);

if test_mode && overlap_mode,
    [overlap,correlation,collectivity,rmsd,similarity]=ANM_merit(network0,target,20);
    sumoverlap=sqrt(cumsum(overlap.^2));
    [maxoverlap,bestmode]=max(overlap);
    [mini,num80]=min(abs(sumoverlap-0.8));
    [mini,num90]=min(abs(sumoverlap-0.9));
    [mini,num95]=min(abs(sumoverlap-0.95));
    [mini,num99]=min(abs(sumoverlap-0.99));
    figure(13); clf;
    plot(overlap,'k');
    hold on;
    plot(sumoverlap,'r');
    axis([0,length(overlap)+1,-0.05,1.05]);
    disp(sprintf('r.m.s.d. between the two structures is %5.3f ?',rmsd));
    disp(sprintf('Collectivity of structural change is %5.3f.',collectivity));
    disp(sprintf('Best mode %i has overlap %5.3f and correlation %5.3f with structural change.',bestmode,maxoverlap,correlation(bestmode)));
    disp(sprintf('%i modes cover 80%% of structural change.',num80));
    disp(sprintf('%i modes cover 90%% of structural change.',num90));
    disp(sprintf('%i modes cover 95%% of structural change.',num95));
    disp(sprintf('%i modes cover 99%% of structural change.',num99));
    disp(sprintf('10 modes cover %4.1f%% of structural change.',100*sumoverlap(10)));
    disp(sprintf('20 modes cover %4.1f%% of structural change.',100*sumoverlap(20)));
    for k=1:20,
         disp(sprintf('Mode %i of template is most similar to mode %i of target with dot product %6.4f',k,similarity(k,1),similarity(k,2)));
    end;
end;

return

slowmodes=20; % number of slowest modes to be used in fit
sigr=3; % standard deviation of distances, coresponds to expected fit quality

snum=model.current_structure;

[m,n]=size(model.coarse(snum).Ca_coor);

mask=zeros(m,1);

for k=1:m,
    ind=model.coarse(snum).indices(k,:);
    [stag,ctag,modelnum,resnum,icode]=mk_address_parts(ind);
    adr1=sprintf('[%s](%s){%i}%i%s','3FH6',ctag,modelnum,resnum,icode);
    [ind2,message]=resolve_address(adr1);
    if ~isempty(ind2) && message.error==0,
        mask(k)=1;
    end;
end;

[shortlist,score,redundant,pairs]=residue_pair_score(slowmodes,0.3,slowmodes,mask);

outname=strcat('2R6G_shortlist.dat');
wfile=fopen(outname,'w');
fprintf(wfile,'%20s%20s%12s\n','Residue 1','Residue 2','Distance');
[mm,nn]=size(shortlist);
for k=1:mm,
    ind1=model.coarse(snum).indices(shortlist(k,2),:);
    adr1=mk_address(ind1);
    ind2=model.coarse(snum).indices(shortlist(k,3),:);
    adr2=mk_address(ind2);
    fprintf(wfile,'%20s%20s%12.2f\n',adr1,adr2,shortlist(k,4));
end;
fclose(wfile);

outname=strcat('2R6G_to_3FH6_only_B_20_restraints.dat');
wfile=fopen(outname,'w');
fprintf(wfile,'%% Calpha-Calpha restraints for driving ATP-bound structure 2R6G chain F\n');
fprintf(wfile,'%% towards apo structure 3FH6 along the %i slowest modes\n',slowmodes);
fprintf(wfile,'# PDB	2R6G\n');
fprintf(wfile,'# basis %i\n',slowmodes);
fprintf(wfile,'# direct\n');
[mm,nn]=size(shortlist);
for k=1:mm,
    ind1=model.coarse(snum).indices(shortlist(k,2),:);
    [stag1,ctag1,modelnum1,resnum1,icode1]=mk_address_parts(ind1);
    adr1=sprintf('[%s](%s){%i}%i%s.CA','3FH6',ctag1,modelnum1,resnum1,icode1);
    [message,xyz1]=get_object(adr1,'coor');
    ind2=model.coarse(snum).indices(shortlist(k,3),:);
    [stag2,ctag2,modelnum2,resnum2,icode2]=mk_address_parts(ind2);
    adr2=sprintf('[%s](%s){%i}%i%s.CA','3FH6',ctag2,modelnum2,resnum2,icode2);
    [message,xyz2]=get_object(adr2,'coor');
    r=norm(xyz1-xyz2);
    fprintf(wfile,'(%s)%i%s    (%s)%i%s%6.2f%6.2f  %% %6.2f\n',ctag1,resnum1,icode1,ctag2,resnum2,icode2,r/10,sigr/10,shortlist(k,4)/10);
end;
fprintf(wfile,'# END\n');
fclose(wfile);

return

modes=handles.degenerate+1:handles.degenerate+10;

network0=model.coarse(model.current_structure).Ca_coor;


% load 1LST_coarse_grained
load 1LST_coarse_grained
% load coarse_grained_1L97A
% target=coarse.Ca_coor(1:162,:); % 1:162
target=coarse.Ca_coor;

% [shortlist,score,redundant,pairs]=residue_pair_score(10);
% [mm,nn]=size(shortlist);
% restraints=zeros(mm,4);
% for k=1:mm,
%     restraints(k,1:2)=shortlist(k,2:3);
%     restraints(k,3)=norm(network0(shortlist(k,2),:)-network0(shortlist(k,3),:));
%     restraints(k,4)=norm(target(shortlist(k,2),:)-target(shortlist(k,3),:));
% end;
% save 1A8E_to_1BP5_restraints restraints
% 
% return

disp_local_change(target,network0);


load 2LAO_to_1LST_restraints
% swap=restraints(:,3);
% restraints(:,3)=restraints(:,4);
% restraints(:,4)=swap;
% load 2LAO_to_1LST_restraints
% sigr=4;
% [m,n]=size(restraints(:,4));
% uncertain=sigr*rand(m,n);
% restraints(:,4)=restraints(:,4)+uncertain;

Hessian=setup_ANM_bonded(network0);
% Hessian=setup_ANM_poly(network0);
[model.ANM(model.current_structure).u,D]=eig(Hessian);
clear D;

% profile on,
% u=reorient_ANM(model.ANM(model.current_structure).u,network0,target);
% profile viewer;
% mauz=max(max(abs(u-model.ANM(model.current_structure).u))),
% return

rmsd0=rmsd_superimpose(target,network0);
% [rmsd,network] = transition_by_ANM(modes,network0,target);
[fom,rmsd,network,dxmat] = fit_by_ANM_1(network0,restraints);
model.fit.Ca_coor=network;
model.fit.rmsd=rmsd;
model.fit.fom=fom;
model.fit.dxmat=dxmat;
disp(sprintf('Final r.m.s.d. of distance restraints: %6.3f ?',rmsd));

[rms,coor2]=rmsd_superimpose(target,network);
disp(sprintf('Initial r.m.s.d. w.r.t. target structure: %6.3f ?',rmsd0));
disp(sprintf('Final   r.m.s.d. w.r.t. target structure: %6.3f ?',rms));

disp_local_change(target,network);

[snum,message]=copy_structure(model.current_structure,'+tes',network);

figure(3); clf;
plot3(target(:,1),target(:,2),target(:,3),'k','LineWidth',1);
hold on
axis equal
plot3(coor2(:,1),coor2(:,2),coor2(:,3),'r','LineWidth',1);

return

figure(1); clf;
plot3(coarse.Ca_coor(:,1),coarse.Ca_coor(:,2),coarse.Ca_coor(:,3),'k','LineWidth',1);
hold on
axis equal



load fit_1L97A_from_1L96
rms=rmsd_superimpose(target,network_fit.Ca_coor);
disp(sprintf('r.m.s.d. between fitted and target structure is %5.2f ?',rms));

[rms,coor,transmat]=rmsd_superimpose(target(1:50,:),network_fit.Ca_coor(1:50,:));
xyz=network_fit.Ca_coor;
[mm,nn]=size(xyz);
xyz=[xyz ones(mm,1)];
xyz=transmat*xyz';
xyz=xyz';

plot3(xyz(:,1),xyz(:,2),xyz(:,3),'g','LineWidth',1);

load coarse_grained_1L96
rms=rmsd_superimpose(target,coarse.Ca_coor);
disp(sprintf('r.m.s.d. between original and target structure is %5.2f ?',rms));
[rms,coor,transmat]=rmsd_superimpose(target(1:50,:),coarse.Ca_coor(1:50,:));
xyz=coarse.Ca_coor;
[mm,nn]=size(xyz);
xyz=[xyz ones(mm,1)];
xyz=transmat*xyz';
xyz=xyz';
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'r','LineWidth',1);


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
ENM_param.parametrization=handles.stored_parametrization;
ENM_param=set_ANM(ENM_param);
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

Ca_coor=model.coarse(snum).Ca_coor;
if isfield(handles,'network_fit'),
    Ca_coor=handles.network_fit.Ca_coor;
end;

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
    if isfield(handles,'network_fit'),
        DEER=handles.network_fit.DEER;
    end;
    for k=1:c,
        set(wires(k),'Color',[0.25,0.25,0.25]);
    end;
    for k=1:length(DEER),
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

if isfield(handles,'DEER'),
    DEER=handles.DEER;
    if isfield(handles,'network_fit'),
        DEER=handles.network_fit.DEER;
    end;
    figure; clf;
    hold on;
    set(gca,'FontSize',16);
    ma=0;
    rmsd=0;
    fom=0; % figure of merit
    for k=1:length(DEER),
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
        plot(k,r0,'o','Color',col,'MarkerFaceColor',col);
    end;
    axis([0,length(DEER)+1,0,1.05*ma]);
    xlabel('Constraint number');
    ylabel('Distance (nm)');
end;


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
function pushbutton_drive_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parametrize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global ENM_param
global hMain

thermal=get(handles.checkbox_thermal,'Value');

drive_mode=get(handles.popupmenu_drive_mode,'Value');
switch drive_mode,
    case {1,2,3}
        update_mode=0;
    case 4
        update_mode=1;
    case 5
        update_mode=2;
end;

correspondence=[];
target=[];
test_mode=true;
overlap_mode=false;
repack=false;

if ~isempty(handles.target),
    [correspondence,target]=align_template_target(mk_address(model.current_structure),handles.target);
    if ~isempty(correspondence),
        test_mode=true;
    end;
end;
network0=model.coarse(model.current_structure).Ca_coor;

if test_mode && overlap_mode,
    [overlap,correlation,collectivity]=ANM_merit(network0(correspondence(1,:),:),target(correspondence(2,:),:));
    sumoverlap=sqrt(cumsum(overlap.^2));
    [maxoverlap,bestmode]=max(overlap);
    [mini,num80]=min(abs(sumoverlap-0.8));
    [mini,num90]=min(abs(sumoverlap-0.9));
    [mini,num95]=min(abs(sumoverlap-0.95));
    [mini,num99]=min(abs(sumoverlap-0.99));
    figure(13); clf;
    plot(overlap,'k');
    hold on;
    plot(sumoverlap,'r');
    axis([0,length(overlap)+1,-0.05,1.05]);
    disp(sprintf('Collectivity of structural change is %5.3f.',collectivity));
    disp(sprintf('Best mode %i has overlap %5.3f and correlation %5.3f with structural change.',bestmode,maxoverlap,correlation(bestmode)));
    disp(sprintf('%i modes cover 80%% of structural change.',num80));
    disp(sprintf('%i modes cover 90%% of structural change.',num90));
    disp(sprintf('%i modes cover 95%% of structural change.',num95));
    disp(sprintf('%i modes cover 99%% of structural change.',num99));
    disp(sprintf('10 modes cover %4.1f%% of structural change.',100*sumoverlap(10)));
    disp(sprintf('20 modes cover %4.1f%% of structural change.',100*sumoverlap(20)));
end;

if test_mode,
    target_sel=target(correspondence(2,:),:);
    template_sel=network0(correspondence(1,:),:);
    [rmsd0,target_sel]=rmsd_superimpose(template_sel,target_sel);
    fprintf(1,'Initial structure r.m.s.d. %4.2f ?\n',rmsd0);
end;

stored_basis=ENM_param.fit_basis;
ENM_param.fit_basis=handles.basis;

switch drive_mode
    case 1
        do_ANM_analytics(template_sel,target_sel);
    case 2
        tic;
        [rmsd_0,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,0,thermal);
        tElapsed=toc;
        fprintf('Drive with original modes took %6.1f s\n',tElapsed);
        tic;
        [rmsd_1,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,1,thermal);
        tElapsed=toc;
        fprintf('Drive with reorientated modes took %6.1f s\n',tElapsed);
        tic;
        [rmsd_2,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,2,thermal);
        tElapsed=toc;
        fprintf('Drive with recomputed modes took %6.1f s\n',tElapsed);
        figure(2); clf;
        plot(rmsd_0,'r');
        hold on;
        plot(rmsd_1,'b');
        plot(rmsd_2,'k');
        rmsd=rmsd_2;
    case 3
        [rmsd,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,0,thermal);
    case 4
        [rmsd,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,1,thermal);
    case 5
        [rmsd,network_a] = drive_by_ANM(template_sel,target_sel,handles.basis,2,thermal);
end;

if drive_mode>1,
    network=network0;
    network(correspondence(1,:),:)=network_a;
    ENM_param.fit_basis=stored_basis;
    add_msg_board(sprintf('New distance r.m.s.d: %6.4f',rmsd(end)));

    handles.network_fit.Ca_coor=network;
    DEER=[];
    handles=update_3D_model(handles,network,DEER);

    model.fit.Ca_coor=network;
    model.fit.rmsd=rmsd;
    [rms,coor2]=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
    add_msg_board(sprintf('Initial r.m.s.d. w.r.t. target structure: %6.3f ?',rmsd0));
    add_msg_board(sprintf('Final   r.m.s.d. w.r.t. target structure: %6.3f ?',rms));
    maxs=length(model.structures);
    newtag=sprintf('tr%i',maxs+1);
    [snum,message]=copy_structure(model.current_structure,newtag,network);
    if repack,
        newtag2=sprintf('tp%i',maxs+1);
        [snum1,infile]=repacked_copy(snum,newtag2,k);
    end;

    if repack
        handles.current_structure=snum1;
    else
        handles.current_structure=snum;
    end;
end;

set(handles.pushbutton_save,'Enable','on');
guidata(hObject,handles),


% --- Executes on selection change in popupmenu_drive_mode.
function popupmenu_drive_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_drive_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_drive_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_drive_mode


% --- Executes during object creation, after setting all properties.
function popupmenu_drive_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_drive_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_fit_mode.
function popupmenu_fit_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_fit_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_fit_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_fit_mode


% --- Executes during object creation, after setting all properties.
function popupmenu_fit_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_fit_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_thermal.
function checkbox_thermal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_thermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_thermal

global ENM_param

if get(hObject,'Value'),
    ENM_param.parametrization=handles.stored_parametrization;
    ENM_param=set_ANM(ENM_param);
    set(handles.popupmenu_fit_mode,'Value',4); % recompute/no local constraints
    ENM_param.diagonalize=true;
    ENM_param.reorientate=false;
    ENM_param.fix_local=0;
else
    ENM_param.parametrization='Jeschke';
    ENM_param=set_ANM(ENM_param);
    set(handles.popupmenu_fit_mode,'Value',3); % recompute/local constraints i, i+1
    ENM_param.diagonalize=true;
    ENM_param.reorientate=false;
    ENM_param.fix_local=2;
end;

set(handles.text_ANM_type,'String',ENM_param.parametrization);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_diagnostics.
function pushbutton_diagnostics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_diagnostics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.fit_diagnostics=handles.fit_info;
inspect_ensemble_fit;


% --- Executes on button press in checkbox_overfit.
function checkbox_overfit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_overfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_overfit


% --- Executes on button press in pushbutton_save_coor_dx.
function pushbutton_save_coor_dx_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_coor_dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

Ca_coor = model.coarse(model.current_structure).Ca_coor;
dxmat = model.fit.dxmat;

[m,n] = size(Ca_coor);
[mt,nt] = size(dxmat);
if mt ~= m,
    add_msg_board('ERROR: Coordinate change matrix does not fit Calpha trace length. Aborting.'); 
    return
end;
frames = nt/n;
if frames - floor(frames) > 3*eps,
    add_msg_board('ERROR: Coordinate change matrix has wrong dimension. Aborting.'); 
    return
end;
[filename, pathname] = uiputfile('*.xyz', 'Save animation raw data in xyz format');
fname=fullfile(pathname, filename);
fid = fopen(fname,'wt');
fprintf(fid,'%i\n',m);
fprintf(fid,'Initial Calpha trace\n');
for k = 1:m,
    fprintf(fid,'C%16.6f%16.6f%16.6f\n',Ca_coor(k,:));
end;
for kf = 1:frames,
    fprintf(fid,'%i\n',m);
    fprintf(fid,'Trajectory frame %i\n',kf);
    dx = dxmat(:,3*kf-2:3*kf);
    Ca_coor = Ca_coor + dx;
    for k = 1:m,
        fprintf(fid,'C%16.6f%16.6f%16.6f\n',Ca_coor(k,:));
    end;
end;
fclose(fid);

function labelname = complabel(labelstr,num)

spoi = strfind(labelstr,'|');
if isempty(spoi)
    labelname = labelstr;
elseif num == 1
    labelname = labelstr(1:spoi-1);
else
    labelname = labelstr(spoi+1:end);
end

function indices = resolve_comp_address(compadr)

spoi = strfind(compadr,'|');
if isempty(spoi)
    indices = resolve_address(compadr);
else
    ind1 = resolve_address(compadr(1:spoi-1));
    ind2 = resolve_address(compadr(spoi+1:end));
    if ~isempty(ind1) && ~isempty(ind2) && length(ind1) == length(ind2)
        indices = [ind1;ind2];
    else
        indices = [];
    end
end