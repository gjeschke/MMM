function varargout = ANM(varargin)
% ANM M-file for ANM.fig
%      ANM, by itself, creates a new ANM or raises the existing
%      singleton*.
%
%      H = ANM returns the handle to a new ANM or the handle to
%      the existing singleton*.
%
%      ANM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANM.M with the given input arguments.
%
%      ANM('Property','Value',...) creates a new ANM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ANM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ANM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ANM

% Last Modified by GUIDE v2.5 06-Apr-2012 09:56:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ANM_OpeningFcn, ...
                   'gui_OutputFcn',  @ANM_OutputFcn, ...
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


% --- Executes just before ANM is made visible.
function ANM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ANM (see VARARGIN)

global MMM_icon
global hMain
global model

test_mode=false;

min_lambda=1e3*eps; % minimum frequency of a non-degenerate mode
handles.darken=0.75; % darken chain colorscheme colors for better visibility and match with lighted ribbon model

% animation control parameters
handles.animation.fps=15;
handles.animation.amp=1;
handles.animation.res=handles.animation.fps;
handles.animation.cycles=5;
handles.animation.scaling=3;
handles.animation.rotation=1;

handles.special_colors=false;

% Choose default command line output for ANM
handles.output = hObject;

set(hMain.figure,'Pointer','watch');

drawnow;

h = msgbox('Please be patient. This can take several minutes.','Anisotropic network model is computed');

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];
hMain.ANM_plot=true;
hMain.ANM_axes=handles.axes_model;

handles.crs1=[];
handles.crs2=[];

if ~exist('model','var') || ~isfield(model,'current_structure'),
    add_msg_board('Error: Cannot set up ANM, since no current structure exists.');
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
set(handles.figure1,'Name',sprintf('Anisotropic network model for structure %s',adr));

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
    if sum(Bfactors)==0,
        Bfactors=20*ones(size(Bfactors));
    end;
    if isempty(Ca_coor) || length(masses)<2,
        set(hMain.figure,'Pointer','arrow');
        add_msg_board('Error: Cannot set up ANM, since current structure has less than two amino acid residues.');
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
    tic;
    % [Hessian,contacts]=setup_ANM(Ca_coor);
    % Hessian=setup_ANM_poly(Ca_coor);
    Hessian=setup_ANM_bonded(Ca_coor);
    contacts=[];
    % [Hessian,contacts]=setup_ANM(Ca_coor,restypes,rindices); % residue-specific force constants, not recommended
    used=toc;
    add_msg_board(sprintf('Setting up the ANM took %6.2f s',used));
    tic;
    [u,D]=eig(Hessian);
    used=toc;
    add_msg_board(sprintf('Diagonalizing the ANM took %6.2f s',used));
    [m,n]=size(Hessian);

    model.ANM(snum).kappa=mode_collectivity(u,model.coarse(snum).masses);
    lambda=zeros(1,m);
    for k=1:m,
        lambda(k)=D(k,k);
    end;
    clear D
    model.ANM(snum).lambda=lambda;
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
    if test_mode,
        D=zeros(m,m);
        for k=degenerate+1:m,
            D=D+u(:,k)*u(:,k)'/lambda(k);
        end;
        corr=zeros(m/3,m/3);
        Bfacs=zeros(1,m/3);
        for i=1:m/3,
            basi=3*(i-1);
            for j=1:m/3,
                basj=3*(j-1);
                submat=D(basi+1:basi+3,basj+1:basj+3);
                corr(i,j)=trace(submat);
                if i==j, Bfacs(i)=corr(i,j); end;
            end;
        end;
        clear D
    end;
    model.ANM(snum).u=u;
    clear u
%     figure(11); clf;
%     corr=100*corr/max(max(corr));
%     image(corr);
%     set(gca,'YDir','normal');
%     axis equal
%     axis tight
%     greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
%     greyscale=flipud(greyscale');
%     colormap(greyscale.^8); 
%     figure(8); clf;
%     plot(Bfacs);
else
    lambda=model.ANM(snum).lambda;
    m=model.ANM(snum).residues;
    msf=model.ANM(snum).msf;
end;

nu=sqrt(lambda(7:end));
inv_nu=ones(size(nu))./nu;
inv_nu=inv_nu/inv_nu(1);

figure(7); clf;
set(gca,'FontSize',14);
plot(inv_nu,'k');
hold on
[~,poi1]=min(abs(inv_nu-exp(-1)));
plot([poi1,poi1],[0,exp(-1)],'r');
[~,poi2]=min(abs(inv_nu-0.2));
plot([poi2,poi2],[0,0.2],'b');
[~,poi3]=min(abs(inv_nu-0.1));
plot([poi3,poi3],[0,0.1],'g');
title(sprintf('%s: 1/e-point: %i, 0.2-point: %i, 0.1-point: %i',adr,poi1,poi2,poi3));
xlabel('Mode number');
ylabel('Normalized inverse mode frequency');

[msf,p]=lin_fit(msf,Bfactors);
handles.Bfactor_scale=p(1);

handles.degenerate=sum(lambda<=min_lambda);
handles.xa=1;
handles.ya=1;
handles.xe=m;
handles.ye=m;
handles.mode=1;
handles.crsr_right=[];
handles.crsr_left=[];
handles.range=[];

if handles.degenerate>6,
    msg='Network not fully connected.';
else
    msg='Network fully connected.';
end;
set(handles.text_info,'String',sprintf('%s %i modes are suspicious.',msg,sum(model.ANM(snum).suspicious)-handles.degenerate));

set(handles.text_collectivity,'String',sprintf('%6.4f',mean(model.ANM(snum).kappa)));

% figure(1); clf;
lambda0=lambda(1+handles.degenerate:end);
invl=ones(size(lambda0))./lambda0;
% plot(invl);
% set(gca,'FontSize',14);
% axis([0,100,0,1.05*max(invl)]);
coverage=cumsum(invl)/sum(invl);
handles.coverage=coverage;
axes(handles.axes_plot);
cla;
plot(coverage,'k');
hold on;
[mi,poi]=min(abs(coverage-0.8));
plot([poi,poi],[0,1],'r:');
plot([1,length(coverage)],[0.8,0.8],'r:');
plot([1,length(coverage)],[0.9,0.9],'b:');
plot([1,length(coverage)],[0.95,0.95],'g:');
axis([1,length(coverage),0,1.1]); 
xlabel('Mode number');
ylabel('Displacement amplitude coverage');
set(handles.text_selected_2,'String',sprintf('80%% covered by %i modes',poi));

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
    % col=handles.darken*color_grade(handles.chains(k,1),poi);
    col=color_grade(handles.chains(k,1),poi);
    if handles.special_colors,
        switch k
            case 3
                col=[1,0,0];
            case 4
                col=[0,128,0]/255;            
            case 5
                col=[128,238,128]/255;
        end;
    end;
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

handles.crsr_mode=0; % selection mode

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ANM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ANM_OutputFcn(hObject, eventdata, handles) 
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

entry=strcat(help_files,'ANM.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.ANM_plot=false;
delete(handles.figure1);


% --- Executes on selection change in popupmenu_plot.
function popupmenu_plot_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot

handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);

function handles=update_auxiliary_plot(hObject,handles)

global model

lambda=model.ANM(model.current_structure).lambda;
my_lambda=lambda(handles.mode+handles.degenerate);
axes(handles.axes_plot);
cla;

mode=handles.current_mode';
msf0=sqrt(sum(mode.^2,2));
[maxmsf,poi]=max(msf0);
vec0=mode(poi,:);
vec0=vec0/norm(vec0);
[m,n]=size(mode);
displace=zeros(1,m);
for k=1:m,
    displace(k)=sum(vec0.*mode(k,:));
end;
plot(displace,'k');
axis([0,length(displace)+1,1.1*min(displace),1.1*max(displace)]);
set(handles.text_selected_2,'String',sprintf('Norm. squared frq. %8.6f',1e3*my_lambda));
% plot(lambda/lambda(handles.degenerate+1),'k');
% hold on;
% ma=max(lambda(1+handles.degenerate:100+handles.degenerate))/lambda(handles.degenerate+1);
% plot([handles.mode+handles.degenerate,handles.mode+handles.degenerate],[-0.05*ma,1.05*ma],'b:');
% axis([handles.degenerate+1,handles.degenerate+100,-0.075*ma,1.075*ma]);
% set(handles.text_selected_2,'String',sprintf('Norm. square freq.: %6.4f',lambda(handles.mode+handles.degenerate)/lambda(handles.degenerate+1)));


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes_model_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton_restraints.
function pushbutton_restraints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_model);
handles.saturation=handles.saturation/2;
colormap(handles.greyscale.^handles.saturation);
guidata(hObject,handles);

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


% --- Executes on button press in pushbutton_detach_plot.
function pushbutton_detach_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

figure(2);
clf;

lambda=model.ANM(model.current_structure).lambda;
plot(lambda/lambda(handles.degenerate+1),'k');
hold on;
ma=max(lambda(1+handles.degenerate:100+handles.degenerate))/lambda(handles.degenerate+1);
plot([handles.mode+handles.degenerate,handles.mode+handles.degenerate],[-0.05*ma,1.05*ma],'b:');
axis([handles.degenerate+1,handles.degenerate+100,-0.075*ma,1.075*ma]);
set(handles.text_selected_2,'String',sprintf('Norm. square freq.: %6.4f',lambda(handles.mode+handles.degenerate)/lambda(handles.degenerate+1)));


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

% --- Executes on button press in pushbutton_minus.
function pushbutton_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.mode>1,
    handles.mode=handles.mode-1;
end;
set(handles.edit_mode_number,'String',sprintf('%5i',handles.mode));
set(handles.text_collectivity,'String',sprintf('%6.4f',model.ANM(model.current_structure).kappa(handles.mode+handles.degenerate)));
m=model.ANM(model.current_structure).residues;
evec=model.ANM(model.current_structure).u(:,handles.mode+handles.degenerate);
mode=reshape(evec,3,m);
handles.current_mode=mode;
handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.mode<model.ANM(model.current_structure).residues-1,
    handles.mode=handles.mode+1;
end;
set(handles.edit_mode_number,'String',sprintf('%5i',handles.mode));
set(handles.text_collectivity,'String',sprintf('%6.4f',model.ANM(model.current_structure).kappa(handles.mode+handles.degenerate)));
m=model.ANM(model.current_structure).residues;
evec=model.ANM(model.current_structure).u(:,handles.mode+handles.degenerate);
mode=reshape(evec,3,m);
handles.current_mode=mode;
handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);


function edit_mode_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mode_number as text
%        str2double(get(hObject,'String')) returns contents of edit_mode_number as a double

global model

[v,handles]=edit_update_MMM(handles,hObject,1,model.ANM(model.current_structure).residues-1,handles.mode,'%5i',1);
handles.mode=v;
set(handles.text_collectivity,'String',sprintf('%6.4f',model.ANM(model.current_structure).kappa(handles.mode+handles.degenerate)));
m=model.ANM(model.current_structure).residues;
evec=model.ANM(model.current_structure).u(:,handles.mode+handles.degenerate);
mode=reshape(evec,3,m);
handles.current_mode=mode;
handles=update_auxiliary_plot(hObject,handles);
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

if nargin<4,
    full=false;
end;

set(handles.text_info,'String','Recording animation...');
set(gcf,'Pointer','watch');

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

if get(handles.checkbox_quick_dirty,'Value'),
    set(handles.axes_model,'DrawMode','fast');
else
    set(handles.axes_model,'DrawMode','normal');
end;

cycles=handles.animation.cycles;
angle_step=360/(cycles*handles.animation.res+1);

Ca_coor=model.coarse(model.current_structure).Ca_coor;
chains=length(handles.wire);

x0=Ca_coor(:,1);
y0=Ca_coor(:,2);
z0=Ca_coor(:,3);
mode=handles.current_mode';
% figure(13); clf;
% plot(mode);
deg=handles.degenerate;
axes(handles.axes_model);
amp=handles.animation.amp;

freq_sq=model.ANM(model.current_structure).lambda(:,handles.mode+handles.degenerate);
freq0_sq=model.ANM(model.current_structure).lambda(:,1+handles.degenerate);
% scale=(freq0_sq/freq_sq);
scale=sqrt(handles.Bfactor_scale/freq_sq);
msf0=sqrt(sum(mode.^2,2));
switch handles.animation.scaling
    case 2
        amp=amp/max(msf0);
    case 3
        amp=amp*scale;
end;

sin_trace=sin(2*pi*linspace(0,1,handles.animation.res+1));
sin_trace=sin_trace(1:end-1);

if record && handles.animation.rotation==1 && ~full,
    cycles=1;
end;

fr=0;
for k0=1:cycles,
    for k=1:handles.animation.res,
        fr=fr+1;
        switch handles.animation.rotation,
            case 2
                camorbit(handles.axes_model,0,angle_step);
            case 3
                camorbit(handles.axes_model,angle_step,0);
        end;
        for c=1:chains,
            ia=handles.chains(c,2);
            ie=handles.chains(c,3);
            x=x0(ia:ie)+mode(ia:ie,1)*sin_trace(k)*amp;
            y=y0(ia:ie)+mode(ia:ie,2)*sin_trace(k)*amp;
            z=z0(ia:ie)+mode(ia:ie,3)*sin_trace(k)*amp;
            if sin_trace(k)==max(sin_trace),
                Cacoor_max(:,1)=x;
                Cacoor_max(:,2)=y;
                Cacoor_max(:,3)=z;
            end;
            if sin_trace(k)==min(sin_trace),
                Cacoor_min(:,1)=x;
                Cacoor_min(:,2)=y;
                Cacoor_min(:,3)=z;
            end;
            set(handles.wire(c),'XData',x,'YData',y,'ZData',z);
        end;
        drawnow;
        if record,
            F(fr)=getframe(gcf,rect);
        end;
    end;
end;

if get(handles.checkbox_copy_max_min,'Value'),
    copy_structure(model.current_structure,'maxamp',Cacoor_max);
    copy_structure(model.current_structure,'minamp',Cacoor_min);
end;


set(handles.axes_model,'DrawMode','normal');
if get(handles.checkbox_black,'Value'),
    set(handles.uipanel_animation,'BackgroundColor',bckg);
    set(handles.uipanel_animation,'ForegroundColor',frg);
end;
for c=1:chains,
    ia=handles.chains(c,2);
    ie=handles.chains(c,3);
    x=x0(ia:ie);
    y=y0(ia:ie);
    z=z0(ia:ie);
    set(handles.wire(c),'XData',x,'YData',y,'ZData',z); drawnow;
end;
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


% --- Executes on button press in pushbutton_steps_minus.
function pushbutton_steps_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_steps_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.animation.res>10,
    handles.animation.res=handles.animation.res-1;
end;
set(handles.edit_steps,'String',sprintf('%i',handles.animation.res));
guidata(hObject,handles);

% --- Executes on button press in pushbutton_steps_plus.
function pushbutton_steps_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_steps_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.animation.res<200,
    handles.animation.res=handles.animation.res+1;
end;
set(handles.edit_steps,'String',sprintf('%i',handles.animation.res));
guidata(hObject,handles);


function edit_steps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_steps as text
%        str2double(get(hObject,'String')) returns contents of edit_steps as a double

[v,handles]=edit_update_MMM(handles,hObject,10,200,36,'%i',1);
handles.animation.res=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to slider_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

v=get(hObject,'Value');
handles.animation.amp=v;
set(handles.text_amp,'String',sprintf('%5.2f',handles.animation.amp));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function handles=set_animation_control(handles)

set(handles.slider_amplitude,'Value',handles.animation.amp);
set(handles.text_amp,'String',sprintf('%5.1f',handles.animation.amp));
set(handles.edit_steps,'String',sprintf('%i',handles.animation.res));
set(handles.edit_cycles,'String',sprintf('%i',handles.animation.cycles));
switch handles.animation.scaling,
    case 1
        set(handles.radiobutton_scale_none,'Value',1);
    case 2
        set(handles.radiobutton_scale_max,'Value',1);
    case 3
        set(handles.radiobutton_scale_inv_freq,'Value',1);
end;
switch handles.animation.rotation,
    case 1
        set(handles.radiobutton_rot_none,'Value',1);
    case 2
        set(handles.radiobutton_rot_x,'Value',1);
    case 3
        set(handles.radiobutton_rot_z,'Value',1);
end;

function handles=query_radiobuttons(handles)

if get(handles.radiobutton_scale_none,'Value'),
    handles.animation.scaling=1;
end;
if get(handles.radiobutton_scale_max,'Value'),
    handles.animation.scaling=2;
end;
if get(handles.radiobutton_scale_inv_freq,'Value'),
    handles.animation.scaling=3;
end;
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

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
cam_up=get(hMain.axes_model,'CameraUpVector');

axes(handles.axes_model);

[filename, pathname] = uiputfile('*.avi', 'Save animation in AVI format');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Save as AVI cancelled by user');
else
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
guidata(hObject,handles);

% --- Executes on button press in checkbox_quick_dirty.
function checkbox_quick_dirty_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_quick_dirty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_quick_dirty


% --- Executes on button press in checkbox_black.
function checkbox_black_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_black (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_black


% --- Executes on button press in pushbutton_test.
function pushbutton_test_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

m=model.ANM(model.current_structure).residues;

hh=gcf;

bas_size = inputdlg({'Basis size (number of slow modes) for covariance matrix computation'},'Specify how many slow modes to include',1,{sprintf('%i',3*m-handles.degenerate)},'on');
bs=str2double(bas_size);
if isnan(bs)|| bs<=0 || bs>3*m-handles.degenerate,
    msgbox(sprintf('Basis size must be a number between 1 and %i',m-handles.degenerate),sprintf('Reverting to basis size %i',m-handles.degenerate),'warn');
    bs=3*m-handles.degenerate;
else
    bs=round(bs);
    if bs<1, bs=1; end;
    if bs>3*m-handles.degenerate, bs=3*m-handles.degenerate; end;
end;

button = questdlg('This might take very long. Answer "No" if you are not sure, "Cancel" to return immediately.','Perform full covariance matrix buildup analysis?','No');


set(hh,'Pointer','watch');
drawnow;
if strcmpi(button,'Yes'),
    snum=model.current_structure;
    corr=zeros(m,m);
    for k=handles.degenerate+1:m-1,
        evec=model.ANM(model.current_structure).u(:,k);
        mode1=reshape(evec,3,m);
        D0=zeros(m,m);
        for l=1:m,
            for ll=l+1:m,
                D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                D0(ll,l)=D0(l,ll);
            end;
        end;
        corr=corr+D0/model.ANM(model.current_structure).lambda(k);
    end;
    if isfield(model,'GNM') && length(model.GNM)>=snum,
        GNM=true;
        corriso=model.GNM(snum).fluctuations;
        [m0,n0]=size(corriso);
        for l=1:m0,
            corriso(l,l)=0;
        end;
        sc=sum(sum(corriso.*corr))/sum(sum(corr.*corr));
        isodiff=corriso-sc*corr;
        agree=1-sum(sum(isodiff.^2))/sum(sum(corriso.^2));
    else
        GNM=false;
    end;
    sqsum=sum(sum(corr.^2));
    err_left=1;
    bs=0;
    coverage=zeros(1,m);
    while err_left>0.01 && bs<m-handles.degenerate,
        bs=bs+1;
        tcorr=zeros(m,m);
        for k=handles.degenerate+1:handles.degenerate+bs-1,
            evec=model.ANM(model.current_structure).u(:,k);
            mode1=reshape(evec,3,m);
            D0=zeros(m,m);
            for l=1:m,
                for ll=l+1:m,
                    D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                    D0(ll,l)=D0(l,ll);
                end;
            end;
            tcorr=tcorr+D0/model.ANM(model.current_structure).lambda(k);
        end;
        sc=sum(sum(corr.*tcorr))/sum(sum(tcorr.*tcorr));
        if isnan(sc), sc=1; end;
        chisq=sum(sum((corr-sc*tcorr).^2));
        err_left=chisq/sqsum;
        coverage(bs)=1-err_left;
        add_msg_board(sprintf('%i modes cover %5.1f%% of covariance.',bs,100*(1-err_left)));
    end;
    coverage=coverage(1:bs);
    axes(handles.axes_plot);
    cla;
    plot(coverage,'k');
    hold on;
    [mi,poi]=min(abs(coverage-0.95));
    plot([poi,poi],[0,1],'r:');
    plot([1,length(coverage)],[0.8,0.8],'r:');
    plot([1,length(coverage)],[0.9,0.9],'b:');
    plot([1,length(coverage)],[0.95,0.95],'g:');
    axis([1,length(coverage),0,1.1]); 
    xlabel('Mode number');
    ylabel('Covariance coverage');
    set(handles.text_selected_2,'String',sprintf('95%% covered by %i modes',poi));
    add_msg_board(sprintf('95%% of covariance covered by %i modes',poi));
    [mi,poi]=min(abs(coverage-0.90));
    add_msg_board(sprintf('90%% of covariance covered by %i modes',poi));
    if bs>=10,
        add_msg_board(sprintf('10 modes cover %5.1f%% of covariance',100*coverage(10)));
    end;
    if bs>=20,
        add_msg_board(sprintf('20 modes cover %5.1f%% of covariance',100*coverage(20)));
    end;
    if GNM,
        add_msg_board(sprintf('%5.1f%% agreement of ANM to GNM covariance matrix.',100*agree));
    end;
else
    corr=zeros(m,m);
    for k=handles.degenerate+1:handles.degenerate+bs-1,
        evec=model.ANM(model.current_structure).u(:,k);
        mode1=reshape(evec,3,m);
        D0=zeros(m,m);
        for l=1:m,
            for ll=l+1:m,
                D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                D0(ll,l)=D0(l,ll);
            end;
        end;
        corr=corr+D0/model.ANM(model.current_structure).lambda(k);
    end;
end;
    ma=max(max(corr));
for k=1:m,
    corr(k,k)=ma;
end;
set(hh,'Pointer','arrow');
figure;
corr=255*corr/max(max(corr));
image(corr);
set(gca,'YDir','normal');
axis equal
axis tight
greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
greyscale=flipud(greyscale');
colormap(greyscale.^16); 

% 
% m=model.ANM(model.current_structure).residues;
% kappa=zeros(1,m);
% disp(handles.degenerate);
% for k=handles.degenerate+1:m,
%     mode=reshape(model.ANM(model.current_structure).u(:,k),3,m);
%     msf0=sum(mode.^2,1)/model.ANM(model.current_structure).lambda(k);
%     alpha=1/sum(msf0);
%     kappa(k)=exp(-sum(alpha*msf0.*log(alpha*msf0)))/m;
% end;
% 
% figure(4); clf;
% plot(kappa,'k.');

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


% --- Executes on button press in radiobutton_scale_none.
function radiobutton_scale_none_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_scale_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_scale_none

ma=get(hObject,'Value');
if ma,
    set(handles.radiobutton_scale_max,'Value',0);
    set(handles.radiobutton_scale_inv_freq,'Value',0);
end;
% --- Executes on button press in radiobutton_scale_max.
function radiobutton_scale_max_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_scale_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_scale_max
ma=get(hObject,'Value');
if ma,
    set(handles.radiobutton_scale_none,'Value',0);
    set(handles.radiobutton_scale_inv_freq,'Value',0);
end;


% --- Executes on button press in radiobutton_scale_inv_freq.
function radiobutton_scale_inv_freq_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_scale_inv_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_scale_inv_freq
ma=get(hObject,'Value');
if ma,
    set(handles.radiobutton_scale_none,'Value',0);
    set(handles.radiobutton_scale_max,'Value',0);
end;


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

threshold=0.4;

cflag=get(handles.checkbox_constraints,'Value');
tflag=get(handles.checkbox_thermal,'Value');
vflag=get(handles.checkbox_vector,'Value');

if cflag,
    sigr=3; % standard deviation of distances, corresponds to expected fit quality
    [m,n]=size(model.ANM(model.current_structure).u);
    if tflag,
        bas_size='';
        bs=[];
        C_cutoff=[];
    elseif vflag,
        bas_size='';
        bs=[];
        C_cutoff=-1;        
    else
        answer = inputdlg({'Number of normal modes in restricted basis'},'Specify basis size for Zheng/Brooks fitting',1,{sprintf('%i',m-6)},'on');
        bas_size = strtrim(char(answer));
        bs=str2double(bas_size);
        if isnan(bs)
            bs=20;
        elseif bs>m-6,
            bs=m-6;
        elseif bs<1,
            bs=20;
        end;
        if bs==m-6,
            C_cutoff=0.1;
        else
            C_cutoff=0.3;
        end;
    end;
    answer = inputdlg({'Number of direct constraints to be generated'},'Specify how many pair restraints to generate',1,{sprintf('%i',20)},'on');
    restraints = strtrim(char(answer));
    np=str2double(restraints);
    if isnan(np)
        np=20;
    end;
    answer = inputdlg('Please provide structure tag','Select target structure');
    stag = strtrim(char(answer));
    if ~strcmp(stag(1),'['), stag=['[' stag]; end;
    if ~strcmp(stag(end),']'), stag=[stag ']']; end;
    [indices,msg]=resolve_address(stag);
    if isempty(indices) || msg.error,
        add_msg_board('ERROR: Target structure does not exist.');
        return;
    end;
    
    if vflag,
        cutoff=-1;
        bs=-1;
        threshold=0.4;
        % shortlist=residue_pair_score_basis(np,[],stag,threshold);
        shortlist=residue_pair_score_vectors(np,[],stag);
    elseif tflag,
        cutoff=[];
        bs=[];
        shortlist=residue_pair_score_thermal(np,[],stag);
    else
        [shortlist,score,redundant,pairs,cutoff]=residue_pair_score(np,C_cutoff,bs,[],stag);
    end;
    if ~isempty(shortlist),
        target_constraints(shortlist,np,bs,cutoff,stag,threshold);
    end;
else
    [filename, pathname] = uiputfile('*.mat', 'Save coarse-graining information in internal format');
    if isequal(filename,0) || isequal(pathname,0)
        add_msg_board('Save of coarse-grained model cancelled by user');
    else
        fname=fullfile(pathname, filename);
        msg=sprintf('Coarse-grained model is saved as .mat file: %s',fname);
        add_msg_board(msg);
    end
    coarse=model.coarse(model.current_structure);
    ENM=model.ANM(model.current_structure);
    save(fname,'coarse','ENM');
end;

guidata(hObject,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global hMain

hMain.ANM_plot=false;
delete(hObject);


% --- Executes on button press in pushbutton_mode_cov.
function pushbutton_mode_cov_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mode_cov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

mnum=handles.mode+handles.degenerate;
m=model.ANM(model.current_structure).residues;

corr=zeros(m,m);

evec=model.ANM(model.current_structure).u(:,mnum);
mode1=reshape(evec,3,m);
corr=zeros(m,m);
for l=1:m,
    for ll=l+1:m,
        corr(l,ll)=sum(mode1(:,l).*mode1(:,ll));
        corr(ll,l)=corr(l,ll);
    end;
end;
ma=max(max(corr));
for k=1:m,
    corr(k,k)=ma;
end;

figure;
corr=255*corr/max(max(corr));
image(corr);
set(gca,'YDir','normal');
axis equal
axis tight
greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
greyscale=flipud(greyscale');
colormap(greyscale.^8); 
title(sprintf('Mode covariance for mode %i',mnum-handles.degenerate));

set(handles.text_collectivity,'String',sprintf('%6.4f',mean(model.ANM(model.current_structure).kappa)));
fprintf(1,'Coverage up to mode %i is %6.2f%%\n',mnum-handles.degenerate,100*handles.coverage(mnum-handles.degenerate));


% --- Executes on button press in pushbutton_monitor.
function pushbutton_monitor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

monitor_distances;

guidata(hObject,handles);


% --- Executes on button press in checkbox_copy_max_min.
function checkbox_copy_max_min_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_copy_max_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_copy_max_min


% --- Executes on button press in pushbutton_user.
function pushbutton_user_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

m=model.ANM(model.current_structure).residues;

hh=gcf;

bas_size = inputdlg({'Basis size (number of slow modes) for covariance matrix computation'},'Specify how many slow modes to include',1,{sprintf('%i',3*m-handles.degenerate)},'on');
bs=str2double(bas_size);
if isnan(bs)|| bs<=0 || bs>3*m-handles.degenerate,
    msgbox(sprintf('Basis size must be a number between 1 and %i',m-handles.degenerate),sprintf('Reverting to basis size %i',m-handles.degenerate),'warn');
    bs=3*m-handles.degenerate;
else
    bs=round(bs);
    if bs<1, bs=1; end;
    if bs>3*m-handles.degenerate, bs=3*m-handles.degenerate; end;
end;

button = questdlg('This might take very long. Answer "No" if you are not sure, "Cancel" to return immediately.','Perform full covariance matrix buildup analysis?','No');


set(hh,'Pointer','watch');
drawnow;
if strcmpi(button,'Yes'),
    snum=model.current_structure;
    corr=zeros(m,m);
    for k=handles.degenerate+1:m-1,
        evec=model.ANM(model.current_structure).u(:,k);
        mode1=reshape(evec,3,m);
        D0=zeros(m,m);
        for l=1:m,
            for ll=l+1:m,
                D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                D0(ll,l)=D0(l,ll);
            end;
        end;
        corr=corr+D0/model.ANM(model.current_structure).lambda(k);
    end;
    if isfield(model,'GNM') && length(model.GNM)>=snum,
        GNM=true;
        corriso=model.GNM(snum).fluctuations;
        [m0,n0]=size(corriso);
        for l=1:m0,
            corriso(l,l)=0;
        end;
        sc=sum(sum(corriso.*corr))/sum(sum(corr.*corr));
        isodiff=corriso-sc*corr;
        agree=1-sum(sum(isodiff.^2))/sum(sum(corriso.^2));
    else
        GNM=false;
    end;
    sqsum=sum(sum(corr.^2));
    err_left=1;
    bs=0;
    coverage=zeros(1,m);
    while err_left>0.01 && bs<m-handles.degenerate,
        bs=bs+1;
        tcorr=zeros(m,m);
        for k=handles.degenerate+1:handles.degenerate+bs-1,
            evec=model.ANM(model.current_structure).u(:,k);
            mode1=reshape(evec,3,m);
            D0=zeros(m,m);
            for l=1:m,
                for ll=l+1:m,
                    D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                    D0(ll,l)=D0(l,ll);
                end;
            end;
            tcorr=tcorr+D0/model.ANM(model.current_structure).lambda(k);
        end;
        sc=sum(sum(corr.*tcorr))/sum(sum(tcorr.*tcorr));
        if isnan(sc), sc=1; end;
        chisq=sum(sum((corr-sc*tcorr).^2));
        err_left=chisq/sqsum;
        coverage(bs)=1-err_left;
        add_msg_board(sprintf('%i modes cover %5.1f%% of covariance.',bs,100*(1-err_left)));
    end;
    coverage=coverage(1:bs);
    axes(handles.axes_plot);
    cla;
    plot(coverage,'k');
    hold on;
    [mi,poi]=min(abs(coverage-0.95));
    plot([poi,poi],[0,1],'r:');
    plot([1,length(coverage)],[0.8,0.8],'r:');
    plot([1,length(coverage)],[0.9,0.9],'b:');
    plot([1,length(coverage)],[0.95,0.95],'g:');
    axis([1,length(coverage),0,1.1]); 
    xlabel('Mode number');
    ylabel('Covariance coverage');
    set(handles.text_selected_2,'String',sprintf('95%% covered by %i modes',poi));
    add_msg_board(sprintf('95%% of covariance covered by %i modes',poi));
    [mi,poi]=min(abs(coverage-0.90));
    add_msg_board(sprintf('90%% of covariance covered by %i modes',poi));
    if bs>=10,
        add_msg_board(sprintf('10 modes cover %5.1f%% of covariance',100*coverage(10)));
    end;
    if bs>=20,
        add_msg_board(sprintf('20 modes cover %5.1f%% of covariance',100*coverage(20)));
    end;
    if GNM,
        add_msg_board(sprintf('%5.1f%% agreement of ANM to GNM covariance matrix.',100*agree));
    end;
else
    corr=zeros(m,m);
    for k=handles.degenerate+1:handles.degenerate+bs-1,
        evec=model.ANM(model.current_structure).u(:,k);
        mode1=reshape(evec,3,m);
        D0=zeros(m,m);
        for l=1:m,
            for ll=l+1:m,
                D0(l,ll)=sum(mode1(:,l).*mode1(:,ll));
                D0(ll,l)=D0(l,ll);
            end;
        end;
        corr=corr+D0/model.ANM(model.current_structure).lambda(k);
    end;
end;
    ma=max(max(corr));
for k=1:m,
    corr(k,k)=ma;
end;
set(hh,'Pointer','arrow');
figure;
corr=255*corr/max(max(corr));
image(corr);
set(gca,'YDir','normal');
axis equal
axis tight
greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
greyscale=flipud(greyscale');
colormap(greyscale.^16); 

% 
% m=model.ANM(model.current_structure).residues;
% kappa=zeros(1,m);
% disp(handles.degenerate);
% for k=handles.degenerate+1:m,
%     mode=reshape(model.ANM(model.current_structure).u(:,k),3,m);
%     msf0=sum(mode.^2,1)/model.ANM(model.current_structure).lambda(k);
%     alpha=1/sum(msf0);
%     kappa(k)=exp(-sum(alpha*msf0.*log(alpha*msf0)))/m;
% end;
% 
% figure(4); clf;
% plot(kappa,'k.');


% --- Executes on button press in checkbox_constraints.
function checkbox_constraints_Callback(~, eventdata, handles)
% hObject    handle to checkbox_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_constraints

function fname=target_constraints(shortlist,np,basis_size,cutoff,stag,threshold)

global model

stag=stag(2:end-1);
add_msg_board(sprintf('Generating target constraints for structure [%s]',stag));

vector_mode=false;
[m,n]=size(shortlist);
if ~isempty(basis_size),
    if basis_size>0,
        fname=[stag sprintf('_%i_direct_constraints_ZB_basis_%i.dat',m,basis_size)];
    else
        fname=[stag sprintf('_%i_direct_constraints_TM_%i.dat',m,round(100*threshold))];
        % fname=[stag sprintf('_%i_direct_constraints_GO.dat',m)];
        vector_mode=true;
    end;
else
    % fname=[stag sprintf('_%i_direct_constraints_equipartitioning.dat',m)];
    fname=[stag sprintf('_%i_direct_constraints_MC.dat',m)];
end;
fid=fopen(fname,'w');
fprintf(fid,'%% MMM simulated direct Calpha constraints for target structure [%s]\n',stag);
fprintf(fid,'%% Requested constraints: %i; Generated constraints: %i\n',np,m);
if ~isempty(cutoff) && cutoff>0,
    fprintf(fid,'%% Correlation cutoff: %5.2f\n',cutoff);
elseif isempty(basis_size),
    fprintf(fid,'%% Equipartitioning mode\n');
else
    fprintf(fid,'%% Vector mode\n');    
    % fprintf(fid,'%% Geometrically optimized mode\n');    
end;

if ~isempty(basis_size) && basis_size>0,
    fprintf(fid,'# BASIS %i\n',basis_size);
end;
fprintf(fid,'# TARGET %s\n',stag);
fprintf(fid,'# DIRECT\n');
rindices=model.coarse(model.current_structure).indices;
for k=1:m,
    cind1=rindices(shortlist(k,2),:);
    [stag0,ctag,modelnum,resnum]=mk_address_parts(cind1);
    adr1=sprintf('(%s)%i',ctag,resnum);
    while length(adr1)<15,
        adr1=[adr1 ' '];
    end;
    cind2=rindices(shortlist(k,3),:);
    [stag0,ctag,modelnum,resnum]=mk_address_parts(cind2);
    adr2=sprintf('(%s)%i',ctag,resnum);
    while length(adr2)<15,
        adr2=[adr2 ' '];
    end;
    if vector_mode,
        fprintf(fid,'%15s  %15s %6.2f  %8.2f  %% %5.3f\n',adr1,adr2,shortlist(k,4)/10,0.03,shortlist(k,1));
    else
        fprintf(fid,'%15s  %15s %6.2f  %8.2f\n',adr1,adr2,shortlist(k,4)/10,0.03);
    end;
end;
fprintf(fid,'# END\n');
fclose(fid);


% --- Executes on button press in checkbox_thermal.
function checkbox_thermal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_thermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_thermal

on_flag=get(hObject,'Value');
if on_flag,
    set(handles.checkbox_vector,'Value',0);
end;

guidata(hObject,handles);

% --- Executes on button press in checkbox_vector.
function checkbox_vector_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vector

on_flag=get(hObject,'Value');
if on_flag,
    set(handles.checkbox_thermal,'Value',0);
end;

guidata(hObject,handles);
